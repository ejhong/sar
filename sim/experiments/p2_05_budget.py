"""P2-05 · The budget: the smallest imprint one dwell could detect, against the imprint there is.

    uv run python experiments/p2_05_budget.py

Everything here is measured in P2-01, P2-03 and P2-04; this experiment only combines it, generously to
the claim, and checks the combination by simulation once.

Detection threshold. A tracker that responds to a fraction r of the motion and whose readings wander by
sigma per look (P2-03, looks W = 1.96 s long) is allowed every one of the N = T / W independent looks the
dwell holds (P2-01), a 3 sigma detection, and the look's gain g at the imprint's frequency (P2-03):

    A_min = 3 sigma / (r g sqrt(N))

which treats every look as a full measurement of the amplitude, better than a real sinusoid allows. For
motion at 15 Hz (traffic) the looks must be shorter than 1 / 2f: sigma grows as 1 / W (a shorter look
resolves less) and N is the looks during a truck's passing (4 s).

Targets: ordinary ground by complex registration (the published protocols, R4), by magnitudes (the best
tracker found), a lone bright scatterer (simulated, P2-03), and a corner reflector 20 dB brighter, a
designed target no archive image contains, whose sigma scales as 1 / sqrt(signal-to-clutter).

Imprints: P2-04's line-of-sight velocity for each ambient case, on the one-chamber bench (a 6 m chamber 15 m
down, shallower and more favourable than anything claimed).

Check: a bright scatterer swaying at 0.2 Hz with the amplitude the formula calls a 3 sigma detection, over
many clutter realisations, recovered by the magnitude tracker from one image each.
"""
import numpy as np

from katabasis.runs import Run, load
from sarsim import synthesize
from sarsim.acquisition import DwellGeometry
from sarsim.looks import look_gain, look_masks, looks, velocity_series

import importlib.util
from pathlib import Path
_spec = importlib.util.spec_from_file_location('p203', Path(__file__).with_name('p2_03_what_a_dwell_sees.py'))
p203 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p203)

SIGMA_K = 3.0
TRUCK_PASS_S = 4.0
REFLECTOR_GAIN_DB = 20.0
CHECK_REALISATIONS = 24
CHECK_F = 0.2


def scr_of_point(g):
    """Signal-to-clutter of P2-03's bright scatterer in a 1.96 s look: its peak intensity over the clutter's mean."""
    shape = (1024, 48)
    img = synthesize(p203.ground(g, shape, np.random.default_rng(21), bright=p203.POINT_AMP), g, shape)
    clutter = synthesize(p203.ground(g, shape, np.random.default_rng(21)), g, shape)
    m = look_masks(g, shape[0], [0.0], p203.W)
    L, C = looks(img, m)[0], looks(clutter, m)[0]
    return float(np.abs(L).max() ** 2 / np.mean(np.abs(C) ** 2))


def check_threshold(g, A):
    """Fitted 0.2 Hz amplitude from one image's magnitude tracking, with and without a sway of amplitude A."""
    shape = (1024, 48)
    centres = np.arange(-11.0, 11.001, p203.W)            # independent looks across the dwell
    m = look_masks(g, shape[0], centres, p203.W)
    rows, cols = slice(shape[0] // 2 - 64, shape[0] // 2 + 64), slice(shape[1] // 2 - 4, shape[1] // 2 + 5)
    X = np.stack([np.cos(2 * np.pi * CHECK_F * centres), np.sin(2 * np.pi * CHECK_F * centres),
                  np.ones_like(centres)], 1)
    fit = lambda s: np.hypot(*np.linalg.lstsq(X, s, rcond=None)[0][:2])
    null, sig = [], []
    for k in range(CHECK_REALISATIONS):
        scene = lambda: p203.ground(g, shape, np.random.default_rng(100 + k), bright=p203.POINT_AMP)
        phase = np.random.default_rng(500 + k).uniform(0, 2 * np.pi)
        motion = lambda x, y, z, t: (A / (2 * np.pi * CHECK_F)) * np.sin(2 * np.pi * CHECK_F * t + phase)[:, None] \
            * (x == 0.0)[None, :]
        null.append(fit(velocity_series(looks(synthesize(scene(), g, shape), m), g, rows, cols, envelope=True)))
        sig.append(fit(velocity_series(looks(synthesize(scene(), g, shape, motion=motion), m), g, rows, cols,
                                       envelope=True)))
    null, sig = np.array(null), np.array(sig)
    thr = np.quantile(null, 0.95)
    return {'amplitude_m_s': A, 'looks': len(centres), 'null_fits_m_s': null, 'signal_fits_m_s': sig,
            'null_rms_m_s': float(np.sqrt(np.mean(null ** 2) / 2)),       # per quadrature: the fitted amplitude's sigma
            'null_95th_m_s': float(thr), 'detected_fraction': float((sig > thr).mean()),
            'signal_over_null_rms': float(np.mean(sig) / np.sqrt(np.mean(null ** 2)))}


def main():
    geo = load('p2_01_geometry')['sites']['giza-20250827']
    sees = load('p2_03_what_a_dwell_sees')
    imp = load('p2_04_chamber_imprint')
    g = DwellGeometry.from_record('giza-20250827')
    W = p203.W
    T = geo['aperture_s']
    N = T / W
    params = {'look_s': W, 'aperture_s': T, 'independent_looks': N, 'detection_sigma': SIGMA_K,
              'truck_pass_s': TRUCK_PASS_S, 'reflector_gain_db': REFLECTOR_GAIN_DB,
              'check': {'realisations': CHECK_REALISATIONS, 'f_hz': CHECK_F}}
    with Run('p2_05_budget', 'The smallest imprint one dwell could detect', params) as run:
        scr = scr_of_point(g)
        real, sim = sees['real']['block'], sees['motion']
        targets = [
            {'target': 'ordinary ground, complex registration (published protocols, R4)',
             'sigma_m_s': real['floor_complex_um_s'] * 1e-6, 'response': real['gain_complex'],
             'response_se': real['gain_complex_se'], 'source': 'P2-03, real Giza texture'},
            {'target': 'ordinary ground, magnitudes',
             'sigma_m_s': real['floor_envelope_um_s'] * 1e-6, 'response': real['gain_envelope'],
             'response_se': real['gain_envelope_se'], 'source': 'P2-03, real Giza texture'},
            {'target': 'a lone bright scatterer', 'sigma_m_s': sim['point']['floor_envelope_um_s'] * 1e-6,
             'response': sim['point']['gain_envelope'], 'response_se': sim['point']['gain_envelope_se'],
             'signal_to_clutter_db': 10 * np.log10(scr), 'source': 'P2-03, simulated on the real geometry'},
            {'target': f'a corner reflector {REFLECTOR_GAIN_DB:.0f} dB brighter (designed; in no archive image)',
             'sigma_m_s': sim['point']['floor_envelope_um_s'] * 1e-6 * 10 ** (-REFLECTOR_GAIN_DB / 20),
             'response': 1.0, 'response_se': 0.0, 'signal_to_clutter_db': 10 * np.log10(scr) + REFLECTOR_GAIN_DB,
             'source': 'the bright scatterer, sigma scaled as 1 / sqrt(signal-to-clutter)'},
        ]
        rows = []
        for case in imp['cases']:
            f = case['frequency_hz']
            for tg in targets:
                r = tg['response']
                if r - 2 * tg['response_se'] <= 0:
                    A_min, note = None, 'no response distinguishable from zero'
                elif f * W < 0.5:
                    A_min = SIGMA_K * tg['sigma_m_s'] / (r * look_gain(f, W) * np.sqrt(N))
                    note = f'{W} s looks, {N:.1f} of them'
                else:
                    w = 1 / (2 * f)
                    transient = 'truck' in case['case']
                    n = (TRUCK_PASS_S if transient else T) / w
                    A_min = SIGMA_K * tg['sigma_m_s'] * (W / w) / (r * look_gain(f, w) * np.sqrt(n))
                    note = (f'{w * 1e3:.0f} ms looks (resolving {g.look_resolution(w):.0f} m), {n:.0f} of them'
                            + (' during a truck\'s pass' if transient else ''))
                rows.append({'case': case['case'], 'target': tg['target'], 'imprint_m_s': case['imprint_los_velocity_m_s'],
                             'motion_m_s': case['vertical_velocity_m_s'], 'a_min_m_s': A_min,
                             'shortfall': None if A_min is None else A_min / case['imprint_los_velocity_m_s'],
                             'shake_needed_m_s': None if A_min is None else
                             case['vertical_velocity_m_s'] * A_min / case['imprint_los_velocity_m_s'],
                             'note': note})
        bright = targets[2]
        A_check = SIGMA_K * bright['sigma_m_s'] / (bright['response'] * look_gain(CHECK_F, W) * np.sqrt(N))
        check = check_threshold(g, A_check)
        check['predicted_sigma_m_s'] = A_check / SIGMA_K
        micro = [r for r in rows if r['case'] == imp['cases'][0]['case']]
        ground = micro[1]
        point, refl = micro[2], micro[3]
        finding = (
            f"One dwell falls short of the chamber by a factor of {ground['shortfall']:.0e} on ordinary ground and "
            f"{point['shortfall']:.0e} even on a bright scatterer directly above it. Under Giza's microseisms the "
            f"chamber's imprint is {ground['imprint_m_s'] * 1e6:.1e} um/s; the smallest a dwell could detect, "
            f"allowing every look, a 3-sigma threshold and the best tracker, is {ground['a_min_m_s'] * 1e6:.0f} um/s on "
            f"the plateau's own texture and {point['a_min_m_s'] * 1e6:.0f} um/s on a bright point "
            f"({refl['a_min_m_s'] * 1e6:.0f} um/s on a corner reflector no archive image contains); the published "
            f"protocols' complex registration responds to about 1% of the motion. For the plateau the ground would have to shake at "
            f"{ground['shake_needed_m_s']:.0f} m/s, beyond the strongest earthquakes recorded, and keep shaking for the "
            f"whole dwell. The threshold is checked by simulation: without motion the bright point's fitted sway "
            f"scatters by {check['null_rms_m_s'] * 1e6:.0f} um/s against the {check['predicted_sigma_m_s'] * 1e6:.0f} um/s "
            f"the formula uses, and a sway at the threshold is found in {check['detected_fraction'] * 100:.0f}% of "
            f"{CHECK_REALISATIONS} realisations at a 5% false-alarm rate.")
        run.save({'targets': targets, 'rows': rows, 'check': check, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
