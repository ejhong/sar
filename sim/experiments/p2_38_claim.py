"""P2-38 · The claim computed: what one image could hold of the structure announced 1,220 m under Khafre.

    uv run python experiments/p2_38_claim.py

Until now the claimed structure was an estimate by scaling (P2-25: its imprint taken as 0.36 of the bench room's), and a
reader goes there first. This computes it, with the benchmark's layers (P2-36, BENCHMARK.md).

The structure: the deeper of the claims (bench-khafre-claim: a void 80 m on a side, its centre 1,220 m below Khafre's
base; its size is representative, the announcement gives depths, not sizes). At 15 times its size below the surface it
is a point source: its imprint on the ground's motion is Eshelby's void (an inclusion equivalent to the hole, for a
sphere of the same volume) seen through Okada's point source in a half-space, exact for a small sphere far from the
surface (the finite-size correction is of order (40 / 1,220)^2); a cube's moment differs from the sphere's by a factor
near one, carried as an assumed allowance of 1.5 in amplitude. Loaded by Giza's regional microseisms (0.1-0.3 Hz,
measured 67 km east), whose 15 km wavelength strains the ground to that depth almost as at the surface; each direction
carrying 1/24 of the power with independent phases, the expected squared motion taken (as P2-36's quiet case). A lorry's
surface wave does not reach it (at 69 Hz it decays within a few tens of metres); the regional trembling is the only
excitation that loads it.

Layers, as P2-36: any reader of the one image (the certificate line by line, each line's reference holding the
microseisms' envelope for every realisation and the receiver noise, fully developed speckle, averaged over the field's
realisations through a bound linear in each line's phase energy); the oracle told every scatterer's reflectivity
(the SNR per cell calibrated from the acquisition, sarsim.radiometry; Khafre's own measured brightness over its 320 m
footprint from the desktop's power map, the site's median ground elsewhere until a 2 km map of summed power arrives). The imprint spreads over kilometres, so the whole 5 km image is observed
and the analytic field covers it; beyond the image nothing is counted (the image is what a method reads).
The eight claimed shafts reach the pyramid's base and are not point sources; they are not computed here.
"""
import json
import time
from pathlib import Path

import numpy as np

from katabasis.runs import Run, load
from katabasis.seismic.analytic import moment_surface_displacement, void_moment

RID = 'p2_38_claim'
SITES = Path(__file__).resolve().parents[2] / 'sites'
GRID_HALF_M = 3600.0          # covers the 5 km image's corners
GRID_STEP_M = 10.0
CUBE_ALLOWANCE = 1.5          # amplitude: a cube's moment against the equal-volume sphere's (assumed)
ALLOWANCES = {'local ambient level over the regional': 10.0, 'site amplification': 3.0}   # assumed, as P2-29


def module(name, path):
    import importlib.util
    import sys
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def structure():
    site = json.loads((SITES / 'bench-khafre-claim' / 'site.json').read_text())
    f = next(x for x in site['features'] if x['id'] == 'deep-structure')
    return {'centre': f['shape']['centre'], 'size': f['shape']['size'], 'source': f['source']}


def rock():
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    r = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    val = lambda q: r[q]['value'] if isinstance(r[q], dict) else r[q]
    vp, vs, rho = val('vp'), val('vs'), val('rho')
    mu, lam = rho * vs ** 2, rho * vp ** 2 - 2 * rho * vs ** 2
    nu = lam / (2 * (lam + mu))
    E = mu * (3 * lam + 2 * mu) / (lam + mu)
    return lam, mu, nu, E


def imprint(g, st):
    """The imprint's slant-range increase per unit strain of a wave towards each of 24 azimuths, on a grid covering the
    image, and e = eps0^2 <M_phi^2> is formed by the caller."""
    lam, mu, nu, E = rock()
    depth = -st['centre'][2]
    vol = float(np.prod(st['size']))
    a = E / (1 - nu * nu)
    moments = [np.diag([a, nu * a, 0.0]), np.diag([nu * a, a, 0.0]), np.array([[0, 2 * mu, 0], [2 * mu, 0, 0], [0, 0, 0.0]])]
    xs = np.arange(-GRID_HALF_M, GRID_HALF_M + 1e-9, GRID_STEP_M)
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    los = np.asarray(g.los_enu)
    L = []
    for Mc in moments:
        u = moment_surface_displacement(void_moment(Mc, vol, nu), depth, X - st['centre'][0], Y - st['centre'][1], lam, mu)
        L.append(-np.tensordot(los, u, axes=(0, 0)))
    az = np.deg2rad(np.arange(24) * 7.5)
    maps = [np.sin(p) ** 2 * L[0] + np.cos(p) ** 2 * L[1] + np.sin(p) * np.cos(p) * L[2] for p in az]
    return xs, maps, {'depth_m': depth, 'volume_m3': vol, 'poisson': nu}


def main():
    m36 = module('p236', Path(__file__).with_name('p2_36_benchmark.py'))
    g = m36.geometry()
    frame = m36.image_frame(g)
    k0 = 4 * np.pi / g.lam
    p30 = load('p2_30_finite_certificates')
    micro_env = p30['background_included']['los_peak_envelope_m']
    # every realisation of the 24-wave field: its peak at most sqrt(24) times the envelope (P2-36)
    Qc_micro = frame['n_band_along'] * (np.sqrt(24) * k0 * micro_env) ** 2
    ql = m36.quiet_level()
    st = structure()
    params = {'structure': st, 'grid_half_m': GRID_HALF_M, 'grid_step_m': GRID_STEP_M, 'cube_allowance': CUBE_ALLOWANCE,
              'allowances': ALLOWANCES, 'allowance_status': 'assumed', 'quiet_level': ql,
              'model': "Eshelby's void (equal-volume sphere) through Okada's point source in a half-space",
              'snr_db_nominal': m36.FROZEN['snr_db_nominal'], 'snr_status': 'calibrated from the acquisition (sarsim.radiometry)',
              'benchmark_frozen_hash': m36.FROZEN_HASH}
    with Run(RID, 'The claim computed: the structure announced 1,220 m under Khafre', params) as run:
        t0 = time.time()
        xs, maps, info = imprint(g, st)
        eps0 = ql['strain_amplitude']
        e = eps0 ** 2 * np.mean(np.square(maps), axis=0)
        peak = float(eps0 * max(np.abs(M).max() for M in maps))
        rms_peak = float(np.sqrt(e.max()))
        half = float((e >= e.max() / 2).sum() * GRID_STEP_M ** 2)
        print(f"  imprint: peak {peak:.2e} m, footprint at half power {half:.3g} m^2 ({time.time() - t0:.0f} s)", flush=True)
        S = m36.line_integrals(g, frame, e, xs, xs, 1e9, None)
        Q = m36.phase_energy(g, S, 0.0)                       # expected over the field's realisations
        e_max = (eps0 / np.sqrt(len(maps)) * np.sum(np.abs(maps), axis=0)) ** 2   # every realisation
        Q_max = m36.phase_energy(g, m36.line_integrals(g, frame, e_max, xs, xs, 1e9, None), 0.0)
        Qc = np.full_like(Q, Qc_micro)
        ens = lambda amp, s2: m36.l1_ensemble(Q * amp ** 2, Q_max * amp ** 2, Qc, s2)
        integ = float(S.sum() * frame['dg'])
        rows = {}
        for label, amp in [('as computed', 1.0), (f'cube allowance x{CUBE_ALLOWANCE:g}', CUBE_ALLOWANCE),
                           ('and the assumed allowances (x10 level, x3 amplification)',
                            CUBE_ALLOWANCE * float(np.prod(list(ALLOWANCES.values()))))]:
            r = {'amplitude_factor': amp, 'l1': {}, 'oracle': {}}
            for snr in m36.FROZEN['snr_db_sweep']:
                c = ens(amp, 10 ** (-snr / 10))                # averaged over realisations, a linear bound
                o = m36.oracle(g, integ * amp ** 2, snr, 0.0)
                o['state'] = m36.state(o['tv_upper'])
                r['l1'][f'{snr:g}'] = c
                r['oracle'][f'{snr:g}'] = o
            s2 = 10 ** (-m36.FROZEN['snr_db_nominal'] / 10)
            r['l1_growth_to_target'] = m36.growth(lambda x: ens(amp * x, s2)['tv_upper'],
                                                  m36.AN['target_tv'])
            o30 = r['oracle'][f"{m36.FROZEN['snr_db_nominal']:g}"]
            r['oracle_growth_to_target'] = float(np.sqrt(m36.oracle_d2_target() / o30['e_delta2']))
            rows[label] = r
        # Khafre's measured brightness over its 320 m footprint (the desktop's power map), in place of the site's median
        # ground there: an interim step until a map of summed power over 2 km replaces the median ground elsewhere
        pm = np.load(SITES / 'acquisitions' / f"{m36.FROZEN['acquisition']}-khafre-power.npz")
        s0 = pm['calibration_factor'] * pm['power_dn2'] * np.sin(np.deg2rad(pm['incidence_deg']))[None, :]
        from sarsim.radiometry import snr_per_cell
        snr_info = snr_per_cell(m36.FROZEN['acquisition'])
        bright = float(np.mean(s0)) / 10 ** (snr_info['ground_sigma0_median_db'] / 10)
        X, Y = np.meshgrid(xs, xs, indexing='ij')
        frac = float(e[(np.abs(X) <= 160) & (np.abs(Y) <= 160)].sum() / e.sum())
        o_local = m36.oracle(g, integ * (1 + (bright - 1) * frac), m36.FROZEN['snr_db_nominal'], 0.0)
        khafre_map = {'footprint_mean_sigma0_db': float(10 * np.log10(np.mean(s0))),
                      'footprint_over_site_median': bright, 'share_of_imprint_in_footprint': frac,
                      'oracle_with_footprint_brightness': o_local,
                      'note': 'the pyramid footprint (320 m) at its measured mean sigma0, the rest at the site median; '
                              'a 2 km map of summed power would replace the rest'}
        bench = load('p2_36_benchmark')['quiet']['presence (room)']
        base = rows['as computed']
        n = f"{m36.FROZEN['snr_db_nominal']:g}"
        worst = rows['and the assumed allowances (x10 level, x3 amplification)']
        finding = (
            f"The structure announced 1,220 m under Khafre, computed rather than scaled: a void 80 m on a side (the size "
            f"representative) as Eshelby's void seen through Okada's point source, loaded by Giza's regional microseisms. Its "
            f"imprint on the ground's line-of-sight motion peaks at {peak:.1e} m and spreads over {half / 1e6:.2f} km^2 at half "
            f"power, so the whole image holds it. Any reader of the one image, bounded line by line with the microseisms in "
            f"each line's reference and averaged over the field's realisations: at most {base['l1'][n]['tv_upper']:.1e} above "
            f"its false-alarm rate ({base['l1'][n]['state']}; {base['l1'][n]['worst_realisation_tv_upper']:.1e} at the worst "
            f"realisation); the oracle told every scatterer's reflectivity, at the calibrated {float(n):.1f} dB: "
            f"{base['oracle'][n]['tv_upper']:.1e}. The motion would have to grow {base['l1_growth_to_target']:.1e} times for "
            f"one image and {base['oracle_growth_to_target']:.1e} times for the oracle; with a cube's moment allowed 1.5 times "
            f"the sphere's, the local level ten times the regional and site amplification three times, still "
            f"{worst['l1_growth_to_target']:.1e} and {worst['oracle_growth_to_target']:.1e}. For comparison the benchmark's "
            f"room under a 5 m roof: {bench['l1_ensemble']['tv_upper']:.1e} and {bench['oracle'][n]['tv_upper']:.1e}. Khafre's "
            f"measured footprint, {10 * np.log10(np.mean(s0)):.1f} dB on average against the site's "
            f"{snr_info['ground_sigma0_median_db']:.1f} dB, holds {100 * frac:.0f}% of the imprint and moves the oracle to "
            f"{o_local['tv_upper']:.1e}. A lorry's "
            f"wave does not reach that depth. The eight claimed shafts are not computed here; beyond Khafre's footprint the oracle "
            f"takes the site's median ground until a 2 km map of the image's summed power replaces it.")
        run.save({'khafre_power_map': khafre_map, 'snr': snr_info, 'imprint': {'peak_los_m': peak, 'rms_over_directions_peak_m': rms_peak, 'half_power_area_m2': half,
                              'integral_e_m4': integ, **info},
                  'rows': rows, 'microseism_common_phase_energy_per_line': Qc_micro, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
