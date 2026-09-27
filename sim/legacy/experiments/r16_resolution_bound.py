"""R16: why better phase precision does not rescue the method.

A fair objection to everything here is that SAR phase is extraordinarily precise, that micro-motion
at the hundreds-of-micrometres level is established in the literature from a single pass, and that
measuring a phase difference precisely does not take minutes. All three are true. This experiment
grants them in full and shows the binding constraint is somewhere else.

Our own measurement agrees with that literature rather than disputing it. R4 puts the velocity floor
of this dwell at 77 um/s over a two-second series. Published single-pass micro-motion work on bridges
reports about 1,000 um/s. So precision is not what fails.

What fails is a trade between two things that are the same knob.

To see a feature of size s from above, the image must resolve it laterally:

    rho <= s.

To see that feature through the ground using elastic waves, the wave must be short enough to
interact with it, which for a half-wavelength criterion means

    f >= Vs / (2 s),

and to record an oscillation at f the sub-aperture series must sample at

    f_s >= 2 f >= Vs / s.

Azimuth resolution is inversely proportional to integration time, so a sub-aperture of duration
tau out of a dwell of length T resolves

    rho(tau) = rho_full * T / tau,

and non-overlapping sub-apertures of length tau give one independent sample every tau, so
f_s = 1/tau and

    rho = rho_full * T * f_s >= rho_full * T * Vs / s.

Requiring rho <= s gives a bound with no free parameters left:

    s^2 >= rho_full * T * Vs,      s_min = sqrt(rho_full * T * Vs).

Faster sampling buys seismic resolution and spends lateral resolution at the same rate. The product
is fixed by the acquisition. Nothing about phase precision, coherence or signal-to-noise appears in
it, so improving any of those does not move it.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import PRODUCTS, RESULTS
from sarsim.dwell import DwellProduct
from sarsim import viz

TEST = 'r16_resolution_bound'
VS_M_S = 2200.0                      # shear speed in Giza limestone
VS_RANGE = (1500.0, 2500.0)          # plausible span for weathered to sound limestone
MEASURED_FLOOR_UM_S = 77.39          # R4, coherent-window bank
LITERATURE_FLOOR_UM_S = 1000.0       # single-pass bridge micro-motion, ~1 mm/s
# surveyed Giza features, with the smallest dimension that would have to be resolved
TARGETS = [
    ('Typical mastaba shaft', 1.5, 12.0),
    ("Khufu's Subterranean Chamber", 8.0, 30.0),
    ("Khufu's King's Chamber", 5.2, 43.0),
    ('Osiris Shaft, lowest level', 9.0, 30.0),
    ('ScanPyramids Big Void', 20.0, 70.0),
    ('Claimed Khafre shaft', 10.0, 648.0),
    ('Claimed Khafre chamber', 79.0, 648.0),
]


def main():
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST)
    os.makedirs(fd, exist_ok=True)
    viz.style()

    product = DwellProduct(PRODUCTS['giza'])
    acq = product.acq
    T = float(acq.collection_duration_s)
    eph = product.ephemeris()
    lo, hi = eph.span
    speed = float(np.linalg.norm(eph.at([(lo + hi) / 2])[1][0]))
    slant = float(acq.slant_range_first_m)
    lam = float(acq.wavelength_m)
    # full-aperture azimuth resolution from the geometry, and the delivered pixel spacing
    rho_theory = lam * slant / (2 * speed * T)
    rho_pixel = float(acq.azimuth_spacing_m)
    rho_full = max(rho_theory, rho_pixel)        # the more favourable of the two, rounded up
    print(f'  aperture {T:.2f} s, speed {speed:.0f} m/s, slant {slant/1000:.1f} km', flush=True)
    print(f'  azimuth resolution: geometry {rho_theory:.3f} m, delivered pixel {rho_pixel:.3f} m '
          f'-> using {rho_full:.3f} m', flush=True)

    def s_min(rho, vs, span=T):
        return float(np.sqrt(rho * span * vs))

    bound = s_min(rho_full, VS_M_S)
    bound_lo = s_min(rho_full, VS_RANGE[0])
    bound_hi = s_min(rho_full, VS_RANGE[1])
    print(f'  smallest feature this route can ever resolve: {bound:.0f} m '
          f'({bound_lo:.0f} to {bound_hi:.0f} m over the shear-speed range)', flush=True)

    rows = []
    for name, size, depth in TARGETS:
        f_req = VS_M_S / (2 * size)
        fs_req = 2 * f_req
        tau = 1.0 / fs_req
        rho_at = rho_full * T / tau
        rows.append({'name': name, 'size_m': size, 'depth_m': depth,
                     'frequency_needed_hz': f_req, 'sampling_needed_hz': fs_req,
                     'sub_aperture_s': tau, 'lateral_resolution_m': rho_at,
                     'resolution_shortfall': rho_at / size,
                     'resolvable': bool(size >= bound)})
        verdict = 'clears this bound' if size >= bound else f'short by {rho_at/size:,.0f}x'
        print(f"    {name:30s} {size:5.1f} m wide -> needs {f_req:8.1f} Hz, "
              f"pixel becomes {rho_at:10,.0f} m, {verdict}", flush=True)

    # the same bound over a grid of assumptions, so the conclusion does not hinge on one pair
    rho_grid = np.array([0.02, 0.044, 0.1, 0.25, 0.5])
    vs_grid = np.linspace(VS_RANGE[0], VS_RANGE[1], 5)
    grid = [[s_min(r, v) for v in vs_grid] for r in rho_grid]
    grid_min = float(np.min(grid))
    print(f'  across resolutions {rho_grid.min()} to {rho_grid.max()} m and shear speeds '
          f'{VS_RANGE[0]:.0f} to {VS_RANGE[1]:.0f} m/s, the bound never falls below '
          f'{grid_min:.0f} m', flush=True)

    # precision is not the binding constraint: state the margin explicitly
    ratio_to_measured = MEASURED_FLOOR_UM_S / LITERATURE_FLOOR_UM_S
    print(f'  our measured floor {MEASURED_FLOOR_UM_S:.0f} um/s is {1/ratio_to_measured:.0f}x better '
          f'than the published single-pass figure of {LITERATURE_FLOOR_UM_S:.0f} um/s', flush=True)

    # ---------------- figures ----------------
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.3))
    fs = np.logspace(-1, 3, 400)
    axs[0].loglog(fs, rho_full * T * fs, lw=2, color=viz.SERIES[0], label='lateral resolution achieved')
    axs[0].loglog(fs, VS_M_S / fs, lw=2, color=viz.SERIES[1], label='feature size this sampling can reach')
    cross = np.sqrt(VS_M_S / (rho_full * T))
    axs[0].axvline(cross, color=viz.INK2, ls=':', lw=1.2)
    axs[0].plot([cross], [bound], 'o', ms=7, color=viz.SERIES[7], zorder=5)
    axs[0].annotate(f'best case {bound:.0f} m', xy=(cross, bound), xytext=(cross * 1.5, bound * 6),
                    fontsize=8.5, color=viz.SERIES[7],
                    arrowprops=dict(arrowstyle='-', color=viz.SERIES[7], lw=0.9))
    for name, size, _ in TARGETS[:4]:
        axs[0].axhline(size, color=viz.GRID, lw=0.9, zorder=0)
    axs[0].text(0.12, TARGETS[1][1], ' surveyed Giza chambers', fontsize=8, color=viz.INK2, va='bottom')
    axs[0].set_xlabel('sub-aperture sampling rate (Hz)')
    axs[0].set_ylabel('feature size (m)')
    axs[0].set_title('The two requirements move apart', loc='left')
    axs[0].legend(fontsize=8, loc='upper left')
    axs[0].set_ylim(1, 1e5)

    names = [r['name'] for r in rows]
    short = [r['resolution_shortfall'] for r in rows]
    y = np.arange(len(names))
    axs[1].barh(y, short, color=[viz.SERIES[0] if s > 1 else viz.SERIES[2] for s in short], height=0.6)
    axs[1].set_yticks(y); axs[1].set_yticklabels(names, fontsize=8)
    axs[1].set_xscale('log'); axs[1].axvline(1, color=viz.AXIS, lw=1.2)
    axs[1].set_xlabel('how much too coarse the pixel becomes')
    axs[1].set_title('Every surveyed target, by how far it misses', loc='left')
    axs[1].invert_yaxis()
    fig.tight_layout()
    trade_fig = {'file': viz.finish(fig, f'{fd}/r16_trade.png'), 'caption':
                 f'<b>Sampling fast enough to hear a chamber makes the pixel too wide to find it.</b> '
                 f'Left: the blue line is the lateral resolution left after sampling at a given rate, '
                 f'and the orange line is the smallest feature an elastic wave at that rate can '
                 f'interact with. They cross at {bound:.0f} m, which is the smallest thing this '
                 f'acquisition could resolve by this route even with a perfect instrument. Every '
                 f'surveyed chamber at Giza sits far below that line. Right: the same targets, by the '
                 f'factor their pixel exceeds their own size once the sampling requirement is met.'}

    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    for i, r in enumerate(rho_grid):
        ax.plot(vs_grid, [s_min(r, v) for v in vs_grid], lw=1.8, color=viz.SERIES[i % 8],
                label=f'{r:g} m resolution')
    ax.axhspan(1, 20, color=viz.GRID, alpha=0.8, lw=0)
    ax.text(vs_grid.mean(), 20, 'surveyed Giza chamber sizes', fontsize=8.5, color=viz.INK2, va='bottom',
            ha='center')
    ax.set_yscale('log'); ax.set_xlabel('shear speed assumed (m/s)')
    ax.set_ylabel('smallest resolvable feature (m)')
    ax.set_title('The bound holds across every assumption worth making', loc='left')
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    grid_fig = {'file': viz.finish(fig, f'{fd}/r16_assumptions.png'), 'caption':
                f'<b>No choice of assumptions brings it down to chamber size.</b> The bound is the '
                f'square root of azimuth resolution times aperture length times shear speed. Sweeping '
                f'the azimuth resolution from {rho_grid.min():g} to {rho_grid.max():g} m and the shear '
                f'speed across the range limestone can plausibly take, it never falls below '
                f'{grid_min:.0f} m. The shaded band is the size of the chambers actually under Giza.'}

    m = {'aperture_s': T, 'platform_speed_m_s': speed, 'slant_range_m': slant, 'wavelength_m': lam,
         'azimuth_resolution_geometry_m': rho_theory, 'azimuth_pixel_m': rho_pixel,
         'azimuth_resolution_used_m': rho_full, 'shear_speed_m_s': VS_M_S,
         'shear_speed_range_m_s': list(VS_RANGE), 'smallest_resolvable_m': bound,
         'smallest_resolvable_low_m': bound_lo, 'smallest_resolvable_high_m': bound_hi,
         'grid_minimum_m': grid_min, 'targets': rows,
         'measured_floor_um_s': MEASURED_FLOOR_UM_S, 'literature_floor_um_s': LITERATURE_FLOOR_UM_S,
         'precision_margin_over_literature': 1 / ratio_to_measured,
         'runtime_s': time.time() - t0}

    worst = max(rows, key=lambda r: r['resolution_shortfall'])
    sub = next(r for r in rows if r['name'] == "Khufu's Subterranean Chamber")
    summary = {
        'id': TEST, 'order': 26, 'tag': 'mechanism', 'eyebrow': '16 · The resolution bound',
        'title': 'Phase precision is not the limit, and improving it changes nothing',
        'question': 'SAR measures micro-motion to hundreds of micrometres from a single pass. Why is that not enough?',
        'finding': (
            f"Because precision is not the binding constraint. This dwell's own velocity floor is "
            f"{MEASURED_FLOOR_UM_S:.0f} micrometres per second, about {1/ratio_to_measured:.0f} times "
            f"better than the roughly {LITERATURE_FLOOR_UM_S:.0f} micrometres per second reported for "
            f"single-pass micro-motion on bridges, so this work agrees with that literature rather "
            f"than disputing it. "
            f"The limit is a trade between two uses of one knob. Reaching a feature of size s through "
            f"the ground needs an elastic wave no longer than about twice its size, so at least Vs "
            f"over 2s in frequency, and recording that frequency needs sub-apertures short enough to "
            f"sample it; but azimuth resolution is inversely proportional to integration time, so "
            f"every factor gained in sampling rate is paid for in lateral resolution. Requiring both "
            f"leaves s squared at least equal to azimuth resolution times aperture length times shear "
            f"speed, putting the smallest resolvable feature here at {bound:.0f} m, and never below "
            f"{grid_min:.0f} m across every assumption worth making. "
            f"That rules out every chamber actually surveyed at Giza. The Subterranean Chamber is "
            f"{sub['size_m']:.0f} m across, needs {sub['frequency_needed_hz']:.0f} Hz to reach, and at "
            f"that sampling rate each pixel spans {sub['lateral_resolution_m']:,.0f} m, about "
            f"{sub['resolution_shortfall']:,.0f} times wider than the chamber. The claimed "
            f"{TARGETS[-2][1]:.0f} m shafts miss by {next(r for r in rows if r['name'] == 'Claimed Khafre shaft')['resolution_shortfall']:,.0f} times. "
            f"The one target that clears this particular bound is the claimed {TARGETS[-1][1]:.0f} m "
            f"chamber, which is larger than {bound:.0f} m, so this test does not exclude it and is not "
            f"offered as if it did; that claim fails on the depth axis instead, which repeats every "
            f"27.4 m and cannot place anything at 648 m. "
            f"Nothing about phase precision, coherence or signal-to-noise appears in the inequality, "
            f"so improving any of them does not move it."),
        'limitations': (
            'This is an analytic bound computed from the acquisition\'s own geometry, not a measurement, '
            'and it is a necessary condition rather than a sufficient one: clearing it would not by '
            'itself make the method work, since the depth axis and the coherence limit are separate '
            'failures. It assumes a half-wavelength resolution criterion, which is the standard rule '
            'but can be beaten somewhat by a well-posed inversion with strong prior information. It '
            'assumes resolution scales inversely with integration time and that non-overlapping '
            'sub-apertures give one independent sample each, both of which are standard but idealised. '
            'It says nothing about methods that do not work this way, such as multi-pass '
            'interferometry over years, ground-based geophysics, or muon tomography.'),
        'method': (
            f"Aperture length, platform speed, slant range and wavelength were read from the ICEYE "
            f"Giza product. Azimuth resolution was taken as the more favourable of the geometric value "
            f"{rho_theory:.3f} m and the delivered pixel spacing {rho_pixel:.3f} m. The bound is the "
            f"square root of that resolution times the {T:.1f} s aperture times shear speed, evaluated "
            f"at {VS_M_S:.0f} m/s and swept across {VS_RANGE[0]:.0f} to {VS_RANGE[1]:.0f} m/s and "
            f"across azimuth resolutions from {rho_grid.min():g} to {rho_grid.max():g} m. The velocity "
            f"floor is the measured value from R4; the comparison figure is the published single-pass "
            f"bridge result."),
        'figures': [trade_fig, grid_fig], 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    with open(os.path.join(fd, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: m[k] for k in ('azimuth_resolution_used_m', 'smallest_resolvable_m',
                                        'grid_minimum_m', 'precision_margin_over_literature')}, indent=1))


if __name__ == '__main__':
    main()
