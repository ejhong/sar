"""R2: the depth axis is a property of the sub-aperture bank, not of the ground.

Two depth scales appear here and they differ by a factor of two. Because the steering
wavenumbers are uniformly spaced, |h(z)| is exactly periodic with 2 pi / dKz: past that the
output is the same output again. The derivative protocol's real two-component fit additionally
cannot tell z from -z, so its own unambiguous range is half of that. Figures and metrics below
report the smaller, more generous number, `nyquist_depth`, and name it as such.

With real orbit state vectors the virtual baseline of one ICEYE dwell is enormous: the
platform moves about 168 km during the 24.5 s aperture. The steering wavenumbers are then
spaced so widely that the depth axis repeats every few metres. Published depths of hundreds
of metres are only reachable by sweeping a small fraction of the aperture, or by choosing a
large "sound wavelength"; both are free parameters of the processing.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import PRODUCTS, BANKS, RESULTS, GIZA_PATCHES

# Only the three published-style designs belong in this comparison; the split-dwell banks
# defined alongside them are a cross-geometry tool, not a depth-axis proposal.
COMPARED = ('reference', 'paper', 'v17')
from sarsim import viz
from sarsim.dwell import DwellProduct
from sarsim.orbit import lla_to_ecef, perpendicular_baselines, steering_wavenumbers
from sarsim.tomo import nyquist_depth

TEST = 'r02_depth_axis'
CLAIMS = [(648.0, 'Khafre wells, 648 m'), (2000.0, 'underground city, 2 km')]


def geometry_at(product, lat, lon, height):
    acq = product.acq
    row, col = product.geolocate(lat, lon, height)
    eph = product.ephemeris()
    t_zd = float(np.clip(product.zero_doppler_time_s(row), *eph.span))
    target = lla_to_ecef(lat, lon, height)
    rate = float(product.doppler_rate_hz_s(col))
    return acq, eph, t_zd, target, rate, row, col


def nyquist_for(eph, t_zd, target, sweep_s, K, lam_s):
    lo, hi = eph.span
    times = np.clip(t_zd + np.linspace(-sweep_s / 2, sweep_s / 2, K), lo, hi)
    b, r, inc = perpendicular_baselines(eph, times, target, t_zd)
    kz = np.sort(steering_wavenumbers(b, r, inc, lam_s))
    return nyquist_depth(kz), r, inc, float(b.max() - b.min())


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    product = DwellProduct(PRODUCTS['giza'])
    _, _, lat, lon, height, _ = GIZA_PATCHES[0][:6]
    acq, eph, t_zd, target, rate, row, col = geometry_at(product, lat, lon, height)
    aperture = product.aperture_time_s(col)
    lam_s = 0.48
    figs = []

    # ---- (a) Nyquist depth against the swept fraction of the aperture
    fracs = np.geomspace(0.002, 1.0, 140)
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    rows_a = []
    for i, K in enumerate([20, 50, 100]):
        zn = np.array([nyquist_for(eph, t_zd, target, f * aperture, K, lam_s)[0] for f in fracs])
        ax.plot(fracs * 100, zn, color=viz.SERIES[i], lw=2, label=f'K = {K} sub-apertures')
        rows_a.append({'K': K, 'nyquist_full_aperture_m': float(zn[-1])})
    for depth, label in CLAIMS:
        ax.axhline(depth, color=viz.MUTED, lw=0.9, ls=':')
        ax.text(0.22, depth * 1.12, label, color=viz.INK2, fontsize=8.5)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('percentage of the 24.5 s aperture actually swept')
    ax.set_ylabel('unambiguous depth range of the real fit (m)')
    ax.set_title('How deep the axis can reach before it folds back', loc='left')
    ax.legend(loc='upper right')
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r02_nyquist.png'), 'caption':
                 '<b>The depth ceiling is bought by discarding aperture.</b> Using real ICEYE state vectors at Khafre, the unambiguous depth range of the steering basis is plotted against the fraction of the 24.5 s aperture the sub-aperture bank sweeps. The range plotted is the one the derivative protocol\'s real two-component fit can use; the output repeats exactly at twice that depth. Sweeping the whole aperture gives a few metres, because the platform moves 168 km and the wavenumbers are spaced very widely. Reaching the published depths requires sweeping about one percent of the aperture, which throws away the very time diversity the method is supposed to exploit.'})

    # ---- (b) the three bank designs on the real geometry
    rows_b = []
    for name in COMPARED:
        bank = BANKS[name]
        plan = bank.plan(acq, 4096, rate)
        zn, r, inc, span = nyquist_for(eph, t_zd, target, max(plan['span_s'], 1e-3), bank.K, lam_s)
        rows_b.append({'bank': name, 'K': bank.K, 'width_frac': bank.width_frac,
                       'shift_hz': bank.shift_hz, 'sub_aperture_s': plan['sub_aperture_s'],
                       'swept_s': plan['span_s'], 'swept_percent': 100 * plan['span_s'] / aperture,
                       'sampling_hz': plan['sampling_hz'], 'pair_lag_ms': 1e3 * plan['pair_lag_s'],
                       'overlap_within_pair': plan['overlap_within_pair'],
                       'baseline_span_m': span, 'nyquist_depth_m': float(zn)})

    # ---- (b2) the derivative protocol's shipped geometry against the real one
    V17 = {'K': 20, 'baseline_span_m': 1000.0, 'slant_range_m': 600000.0,
           'incidence_deg': 35.0, 'lam_s_m': 0.2405679137, 'z_max_m': 100.0}
    assumed = np.linspace(-V17['baseline_span_m'] / 2, V17['baseline_span_m'] / 2, V17['K'])
    kz_assumed = np.sort(steering_wavenumbers(assumed, V17['slant_range_m'], V17['incidence_deg'], V17['lam_s_m']))
    real_span = [r for r in rows_b if r['bank'] == 'v17'][0]['baseline_span_m']
    real_b = np.linspace(-real_span / 2, real_span / 2, V17['K'])
    kz_real = np.sort(steering_wavenumbers(real_b, V17['slant_range_m'], V17['incidence_deg'], V17['lam_s_m']))
    protocol = {'assumed_baseline_span_m': V17['baseline_span_m'],
                'assumed_unambiguous_depth_m': float(nyquist_depth(kz_assumed)),
                'real_baseline_span_m': float(real_span),
                'real_unambiguous_depth_m': float(nyquist_depth(kz_real)),
                'shipped_depth_grid_max_m': V17['z_max_m'],
                'baseline_underestimate_factor': float(real_span / V17['baseline_span_m'])}

    # ---- (c) lambda_s freedom
    lam_grid = np.geomspace(0.05, 20.0, 200)
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    for i, name in enumerate(COMPARED):
        bank = BANKS[name]
        plan = bank.plan(acq, 4096, rate)
        base = nyquist_for(eph, t_zd, target, max(plan['span_s'], 1e-3), bank.K, 1.0)[0]
        ax.plot(lam_grid, base * lam_grid, color=viz.SERIES[i], lw=2,
                label=f"{name}: sweeps {100 * plan['span_s'] / aperture:.1f}% of the aperture")
    for depth, label in CLAIMS:
        ax.axhline(depth, color=viz.MUTED, lw=0.9, ls=':')
        ax.text(0.055, depth * 1.15, label, color=viz.INK2, fontsize=8.5)
    ax.axvline(0.48, color=viz.SERIES[7], lw=1, ls='--')
    ax.text(0.5, 4, 'paper value 0.48 m', color=viz.SERIES[7], fontsize=8.5, rotation=90, va='bottom')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('declared "sound wavelength" $\\lambda_s$ (m)')
    ax.set_ylabel('unambiguous depth range (m)')
    ax.set_title('The same measurement, relabelled: depth scales linearly with a free constant', loc='left')
    ax.legend(loc='lower right', fontsize=8.5)
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r02_lambda.png'), 'caption':
                 '<b>Nothing in the data sets the scale.</b> The depth axis is proportional to the declared sound wavelength, which no measurement in the method constrains. For each of the three bank designs the reachable depth is drawn against that constant. Any published depth can be produced by choosing it, and the 2022 paper and the 2026 derivative protocol use values that differ by a factor of two.'})

    for row_ in rows_a:
        row_['exact_repetition_period_m'] = 2 * row_['nyquist_full_aperture_m']
    for row_ in rows_b:
        row_['exact_repetition_period_m'] = 2 * row_['nyquist_depth_m']
    protocol['real_exact_repetition_period_m'] = 2 * protocol['real_unambiguous_depth_m']
    m = {'aperture_s': aperture, 'row': row, 'col': col, 'lam_s_m': lam_s,
         'depth_convention': ('nyquist_depth is pi/dKz, the range the real two-component fit can use; '
                              'the complex output repeats exactly at twice that'),
         'nyquist_full_aperture': rows_a, 'banks': rows_b, 'derivative_protocol_v17': protocol,
         'percent_of_aperture_for_648m_K50': float(100 * 3.386 * 49 / 648 / aperture),
         'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 12, 'tag': 'artifact', 'eyebrow': '02 · Depth axis',
        'title': 'The reachable depth is set by the sub-aperture bank, not by the ground',
        'question': 'With real orbit geometry, how deep can this method see before its depth axis repeats?',
        'finding': (f"A few metres. During the 24.5 s aperture the platform moves 168 km, so a bank that uses the "
                    f"whole aperture spaces its steering wavenumbers very widely: at K = 50 the depth axis folds "
                    f"back at {rows_a[1]['nyquist_full_aperture_m']:.1f} m and repeats exactly at "
                    f"{2 * rows_a[1]['nyquist_full_aperture_m']:.1f} m. Depths of hundreds of metres are only "
                    f"reachable by sweeping about one percent of the aperture, or by enlarging the declared sound "
                    f"wavelength, which no measurement constrains. The same calculation applies to the public "
                    f"derivative protocol: its shipped example assumes a one-kilometre aperture span, where the real "
                    f"span for that bank design is {protocol['real_baseline_span_m']:,.0f} m, so its unambiguous "
                    f"depth is {protocol['real_unambiguous_depth_m']:.0f} m rather than "
                    f"{protocol['assumed_unambiguous_depth_m']:.0f} m, and its shipped depth grid runs past that "
                    f"limit."),
        'limitations': ('This is a statement about the ambiguity period of the steering basis under the published '
                        'geometry, computed from real state vectors. It does not by itself show what a given '
                        'published figure did, because the bank parameters used there were never disclosed. The '
                        'comparison with the derivative protocol uses the placeholder geometry in its shipped '
                        'example configuration, which its documentation asks the user to replace with real values; '
                        'the point is what happens when they are replaced.'),
        'method': ('Real ICEYE X33 state vectors, Khafre RPC position and WGS84 target, aperture-centre reference. '
                   'Baselines are the along-track platform offsets projected perpendicular to the line of sight, '
                   'exactly as the published steering formula requires. Nyquist depth is pi divided by the median '
                   'spacing of the resulting steering wavenumbers.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({'banks': rows_b, 'nyquist': rows_a}, indent=1))


if __name__ == '__main__':
    main()
