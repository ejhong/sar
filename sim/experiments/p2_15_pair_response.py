"""P2-15 · What the six-second pairs register of a vibration, measured rather than assumed.

    uv run --with scikit-image python experiments/p2_15_pair_response.py

P2-13 quoted a boxcar average, sinc(f T), for how much of a motion each 6 s image of a pair keeps. That is the transfer
function of a rectangular average, not of a registration: P2-03's bright point, tracked by magnitudes through 1.96 s
looks, came back at 21% at 0.5 Hz where the boxcar gives 2%. Here the exact configuration is measured instead. The
gated reconstruction's own bank (its frequency_masks: 317 pairs of 32,330 Hz masks 404 Hz apart, leaping 101 Hz,
over the central 64,660 Hz) is cut from simulated images on the 2022 X13 geometry, and each pair is registered by its
own register_one (32 x 32 cores, upsampled 100 times, rounded to 0.01 px).

Scenes: one bright scatterer alone; the same scatterer, 30 dB above clutter; a 20 m block of clutter. Motion is put
into the image the way the lab's synthesizer and injections do it, as a phase history through the local Doppler-time
mapping (t = f / Ka): line-of-sight sinusoids of 2 and 20 mm/s at the ten frequencies the gate's modes stand for, and
at 0.1, 0.26, 0.5 and 3.66 Hz. This is an approximation to a vibrating scene (no range migration, occlusion or
changing reflectivity), shared by every motion test in this lab.

Scored, before any gate: the amplitude at the planted frequency of the pair shifts (moved minus still, azimuth and
range), against an instantaneous tracker between each pair's two moments and against the boxcar; and, as a positive
control that the motion is in the images, each image read against its motionless twin, which no real dwell has. After
the first gate: how many of the 268 windows pass for the moved and the still scene, and at which modes.
"""
import json

import numpy as np
import scipy.fft as sfft

import p2_13_gates as m
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.looks import inject_region_motion

RID = 'p2_15_pair_response'
SHAPE = (8192, 144)                      # 18.6 Hz bins: the 101 Hz leap is 5.4 bins; 144 columns hold the 128 px frame
SPEEDS = [2e-3, 2e-2]                    # m/s along the line of sight
EXTRA_HZ = [0.1, 0.26, 0.5, 3.66]
SCR_DB = 30.0                            # the bright point's peak over the clutter's mean intensity
BLOCK_M = 10.0                           # half-length of the moving block of clutter, metres of azimuth


def band_noise(g, shape, rng):
    """Complex Gaussian clutter confined to the image's processed band in both axes: fully developed speckle."""
    z = (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)
    Z = np.fft.fft2(z)
    fa = np.fft.fftfreq(shape[0], d=1.0 / g.prf)
    kr = np.fft.fftfreq(shape[1], d=g.dr)
    Z *= (np.abs(fa) <= g.doppler_band_hz / 2)[:, None] * (np.abs(kr) <= g.kr_band / 2)[None, :]
    out = np.fft.ifft2(Z)
    return (out / np.sqrt(np.mean(np.abs(out) ** 2))).astype(np.complex64)


def point(g, shape):
    """One scatterer at the image centre, band-limited in both axes, peak 1."""
    z = np.zeros(shape, np.complex64)
    z[shape[0] // 2, shape[1] // 2] = 1.0
    out = np.fft.ifft2(np.fft.fft2(z) * (np.abs(np.fft.fftfreq(shape[0], d=1.0 / g.prf)) <= g.doppler_band_hz / 2)[:, None]
                       * (np.abs(np.fft.fftfreq(shape[1], d=g.dr)) <= g.kr_band / 2)[None, :])
    return (out / np.abs(out).max()).astype(np.complex64)


def spectrum(img):
    return sfft.fftshift(sfft.fft(img, axis=0, workers=8), axes=0)


def band(S, mask):
    return sfft.ifft(sfft.ifftshift(S * mask[:, None], axes=0), axis=0, workers=8)


def register_series(img, masks, reg, gated, still=None):
    """Pair shifts [K, 2] (range, azimuth px) at the image centre; with `still`, each moved reference image is also
    read against the still one's (the paired reading, available only in simulation)."""
    r0, c0 = SHAPE[0] // 2, SHAPE[1] // 2
    S, S0 = spectrum(img), (spectrum(still) if still is not None else None)
    Y, P = [], []
    for k in range(masks.shape[0]):
        a, b = band(S, masks[k, 0]), band(S, masks[k, 1])
        Y.append(gated.register_one(a, b, r0, c0, reg)[0])
        if S0 is not None:
            P.append(gated.register_one(band(S0, masks[k, 0]), a, r0, c0, reg)[0])
    return np.asarray(Y), (np.asarray(P) if P else None)


def fit(series, t, f):
    """Amplitude and phase of each component at frequency f, with a constant: [2] amplitudes and the 2 x 2 map."""
    A = np.stack([np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t), np.ones_like(t)], 1)
    coef, *_ = np.linalg.lstsq(A, series, rcond=None)
    ellipse = coef[:2].T                                     # rows: range, azimuth; columns: cos, sin
    s = np.linalg.svd(ellipse, compute_uv=False)
    return np.hypot(coef[0], coef[1]), float(s[1] / max(s[0], 1e-15))


def main():
    g = DwellGeometry.from_record(m.ACQ)
    profile = json.loads(m.PROFILE.read_text())
    bank = profile['processing']['filter_bank']
    reg = profile['processing']['registration']
    gated = m.load_gated()
    masks, records = gated.frequency_masks(SHAPE[0], 1.0 / g.prf, profile)
    t_pair = np.array([r['pair_center_hz'] for r in records]) / g.Ka_signed          # each pair's moment, seconds
    t_ref = np.array([r['reference_center_hz'] for r in records]) / g.Ka_signed
    T = bank['mask_width_hz'] / g.Ka                                                   # one image's span
    dt = bank['bshift_hz'] / g.Ka                                                      # the pair's two moments apart
    window_s = 50 * bank['k_leap_hz'] / g.Ka
    mode_hz = [k / window_s for k in range(1, 11)]
    freqs = sorted(set([round(f, 4) for f in mode_hz + EXTRA_HZ]))
    az_m_per_m_s = 2 * g.V / (g.lam * g.Ka)
    params = {'acquisition': m.ACQ, 'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE),
              'shape': SHAPE, 'speeds_m_s': SPEEDS, 'frequencies_hz': freqs, 'scr_db': SCR_DB, 'block_half_m': BLOCK_M}
    with Run(RID, 'What the six-second pairs register of a vibration, measured rather than assumed', params) as run:
        rng = np.random.default_rng(15)
        clutter = band_noise(g, SHAPE, rng)
        pt = point(g, SHAPE)
        amp = np.sqrt(10 ** (SCR_DB / 10))
        rows = np.arange(SHAPE[0])
        block = (np.abs(rows - SHAPE[0] // 2) * g.dx < BLOCK_M).astype(float)
        scenes = {
            'point': (pt, None, np.ones(SHAPE[0])),
            'point_on_clutter': (amp * pt, clutter, np.ones(SHAPE[0])),
            'texture': (clutter, None, block),
        }
        cfg = {'gates': profile['gates']}
        out = {'geometry': {'image_span_s': T, 'pair_dt_s': dt, 'window_s': window_s, 'mode_hz': mode_hz,
                            'azimuth_m_per_m_s': az_m_per_m_s, 'pixel_m': [g.dx, g.dr]}, 'scenes': {}}
        for name, (moving, fixed, weight) in scenes.items():
            still = moving + (fixed if fixed is not None else 0)
            Y0, _ = register_series(still, masks, reg, gated)
            gate0, _ = gated.modal_gate(Y0[None].astype(np.float64), cfg)
            rows_out = []
            for v in SPEEDS:
                for f in freqs:
                    d = lambda t, f=f, v=v: (v / (2 * np.pi * f)) * np.sin(2 * np.pi * f * t)
                    moved = inject_region_motion(moving, g, weight, d) + (fixed if fixed is not None else 0)
                    Y1, P1 = register_series(moved, masks, reg, gated, still=still)
                    amp_pair, ratio = fit(Y1 - Y0, t_pair, f)
                    amp_paired, _ = fit(P1, t_ref, f)
                    inst = az_m_per_m_s * v * 2 * abs(np.sin(np.pi * f * dt)) / g.dx
                    box = inst * abs(np.sinc(f * T))
                    gate1, _ = gated.modal_gate(Y1[None].astype(np.float64), cfg)
                    won = gate1['window_mode'][0][gate1['window_pass'][0]]
                    rows_out.append({
                        'v_m_s': v, 'f_hz': f, 'mode_equivalent': f * window_s,
                        'pair_amp_px': {'range': float(amp_pair[0]), 'azimuth': float(amp_pair[1])},
                        'trace_width_ratio': ratio,
                        'instantaneous_px': float(inst), 'boxcar_px': float(box),
                        'gain_vs_instantaneous': float(amp_pair[1] / inst),
                        'paired_amp_px': float(amp_paired[1]),
                        'paired_boxcar_px': float(az_m_per_m_s * v * abs(np.sinc(f * T)) / g.dx),
                        'windows_passing': int(gate1['window_pass'].sum()),
                        'modes_passing': {int(k): int(n) for k, n in zip(*np.unique(won, return_counts=True))},
                    })
                    print(f"  {name} {v * 1e3:.0f} mm/s {f:6.3f} Hz: pair {amp_pair[1]:.4f} px (boxcar {box:.4f}, "
                          f"instantaneous {inst:.4f}); paired {amp_paired[1]:.4f} px; windows {int(gate1['window_pass'].sum())}",
                          flush=True)
            out['scenes'][name] = {'still_windows_passing': int(gate0['window_pass'].sum()),
                                   'still_exact_zero': float(np.mean(Y0 == 0)), 'rows': rows_out}

        # the finding, from the measurements
        def row(scene, v, f):
            return next(r for r in out['scenes'][scene]['rows'] if r['v_m_s'] == v and abs(r['f_hz'] - f) < 1e-3)
        pr = [row('point', 2e-3, f) for f in (0.26, 3.66)]
        m1 = row('point', 2e-2, round(mode_hz[0], 4))
        tx = [row('texture', 2e-3, f) for f in (0.26, 3.66)]
        ratios = [r['trace_width_ratio'] for s in out['scenes'].values() for r in s['rows'] if r['v_m_s'] == 2e-2]
        finding = (
            f"Through the reconstruction's own masks and registration, on simulated images of the 2022 geometry, a bright "
            f"point swaying at 2 mm/s along the line of sight moves the pair shifts by {pr[0]['pair_amp_px']['azimuth']:.4f} "
            f"px at 0.26 Hz and {pr[1]['pair_amp_px']['azimuth']:.4f} px at 3.66 Hz, where a boxcar average predicts "
            f"{pr[0]['boxcar_px']:.4f} and {pr[1]['boxcar_px']:.4f} px and a tracker reading each pair's two moments "
            f"instantaneously {pr[0]['instantaneous_px']:.4f} and {pr[1]['instantaneous_px']:.4f} px; at 20 mm/s and mode 1 "
            f"({mode_hz[0]:.2f} Hz) it registers {m1['gain_vs_instantaneous'] * 100:.1f}% of the instantaneous tracker's "
            f"shift. Read against its motionless twin, each image does carry the motion: {m1['paired_amp_px']:.3f} px "
            f"against a boxcar {m1['paired_boxcar_px']:.3f} px at mode 1. A 20 m block of texture at 2 mm/s moves the pair "
            f"shifts by {tx[0]['pair_amp_px']['azimuth']:.4f} and {tx[1]['pair_amp_px']['azimuth']:.4f} px. The registered "
            f"traces are {min(ratios):.2g} to {max(ratios):.2g} as wide as they are long.")
        run.save({**out, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
