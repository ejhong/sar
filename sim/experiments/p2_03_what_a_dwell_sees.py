"""P2-03 · What one dwell can see of the ground's motion.

    uv run python experiments/p2_03_what_a_dwell_sees.py

Every velocity series in the first investigation, and every trajectory in the published protocols, comes
from registering looks of one image against each other. This experiment asks what such a registration can
see, on the real Giza geometry (P2-01) and on the real Giza image.

A. Coherence between looks. Two looks cut from one image share their overlapping spectrum sample for
   sample, so their complex coherence is a function of the image's power spectrum alone: 1 - delta/W for a
   flat spectrum. The first investigation's R9 measured 0.98 s to half coherence and 1.76 s to a tenth with
   looks 1.96 s wide. The same numbers come out of simulated desert with no decorrelation of any kind, and
   out of the real crop's own power spectrum to six decimals: R9 measured the looks, not the ground.
B. Motion. A 20 m stretch of simulated ground, and separately one bright scatterer, move at 2 mm/s,
   0.1 Hz, on the real geometry. Sixteen looks 1.96 s long, sampled at 8 Hz (R4's bank), are registered
   against the centre look of the same image, as one dwell allows, by complex correlation (the published
   protocols, R4) and by correlation of magnitudes. Reading each look against the same look of the
   motionless scene, possible only in simulation, shows the motion the data hold.
C. The look's own filter. A look averages motion over its duration, gain sinc(f W). The bright point's
   motion is recovered at frequencies from 0.05 to 3 Hz with 1.96 s looks; R7 injected at 1, 2 and 3 Hz.
D. The real image. The same block motion is imposed on a real crop of open plateau west of the pyramids,
   and the same trackers run on genuine radar texture.
"""
import json
from pathlib import Path

import numpy as np

from katabasis.runs import Run, RESULTS
from sarsim import synthesize
from sarsim.acquisition import DwellGeometry
from sarsim.looks import (coherence, coherence_from_power, inject_region_motion, look_gain, look_masks, looks,
                          paired_velocity, velocity_series)
from sarsim.scene import Scatterers

RID = 'p2_03_what_a_dwell_sees'
PRODUCT = Path.home() / 'tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5'
DESERT_WEST = (29.97600, 31.11800, 80 + 15.5)      # legacy real_common.GIZA_PATCHES 'desert_west': lat, lon, height
R9 = RESULTS.parent / 'legacy' / 'results' / 'r09_coherence' / 'summary.json'
W = 1.96                                           # s: the first investigation's looks (8% of the band)
R4_CENTRES = 2.5 + np.linspace(-0.93, 0.93, 16)    # R4's sixteen looks at about 8 Hz, placed 2.5 s after the
                                                   # aperture centre, where the 0.1 Hz sway changes fastest
SHAPE = (2048, 128)
F_MOTION, V_MOTION = 0.1, 2e-3                     # Hz, m/s
BLOCK_HALF = 10.0                                  # m of azimuth
POINT_AMP = 300.0                                  # a bright scatterer's amplitude against clutter scatterers of 1
LOWPASS_F = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0]


def ground(g, shape, rng, bright=None):
    """Uniform random scatterers, one per full-resolution cell on average, and optionally one bright one at 0."""
    Lx, Lr = shape[0] * g.dx, shape[1] * g.dr
    n = int((Lx / g.resolution) * (Lr / (0.886 / g.kr_band)))
    x, r, a = rng.uniform(-Lx / 2, Lx / 2, n), rng.uniform(-Lr / 2, Lr / 2, n), np.ones(n)
    if bright is not None:
        x, r, a = np.r_[x, 0.0], np.r_[r, 0.0], np.r_[a, bright]
    N = len(x)
    z = np.zeros(N)
    return Scatterers(x=x, y=r / np.sin(g.theta), z=z, amp=a, phase=rng.uniform(0, 2 * np.pi, N), iso=np.ones(N),
                      flash=z, nu0=z, sig_nu=np.ones(N), vib_amp=z, vib_freq=z, vib_phase=z, label=np.zeros(N, int))


def sway(f, v, select):
    """Line-of-sight displacement v/(2 pi f) sin(2 pi f t) for scatterers where select(x) is true."""
    return lambda x, y, z, t: (v / (2 * np.pi * f)) * np.sin(2 * np.pi * f * t)[:, None] * select(x)[None, :]


def look_average(f, v, centres, width):
    """The velocity a look sees: v cos(2 pi f t) averaged over the look."""
    return v * np.cos(2 * np.pi * f * np.asarray(centres)) * look_gain(f, width)


def response(recovered, injected, noise):
    """Least-squares gain of recovered on injected, and its standard error from the motionless scene's own series
    (the same estimator on the same image without the motion: all it holds is noise)."""
    inj = np.asarray(injected, float)
    return float(np.dot(recovered, inj) / np.dot(inj, inj)), float(np.std(noise) / np.linalg.norm(inj))


def part_a(g):
    rng = np.random.default_rng(11)
    img = synthesize(ground(g, SHAPE, rng), g, SHAPE)
    deltas = np.arange(-7, 8) * 0.3714                     # R9's separations
    m = look_masks(g, SHAPE[0], np.r_[0.0, deltas], W)
    L = looks(img, m)
    sim = [float(coherence(L[0], L[k + 1])) for k in range(len(deltas))]
    r9 = json.loads(R9.read_text())['metrics']
    meas = [(row['delta_s'], row['coherence']) for row in r9['curves']['desert_west'] if abs(row['delta_s']) < 2.7]
    return {'delta_s': deltas, 'simulated': sim, 'overlap': np.clip(1 - np.abs(deltas) / W, 0, None),
            'r9_measured': {'delta_s': [d for d, _ in meas], 'coherence': [c for _, c in meas]},
            'r9_half_s': [p['half_coherence_s'] for p in r9['patches']],
            'r9_tenth_s': [p['tenth_coherence_s'] for p in r9['patches']], 'look_s': W}


def part_b(g):
    plain = lambda: ground(g, SHAPE, np.random.default_rng(12))                  # ordinary ground only
    lit = lambda: ground(g, SHAPE, np.random.default_rng(12), bright=POINT_AMP)  # the same, plus one bright scatterer
    m = look_masks(g, SHAPE[0], R4_CENTRES, W)
    L0 = looks(synthesize(plain(), g, SHAPE), m)
    LB = looks(synthesize(plain(), g, SHAPE, motion=sway(F_MOTION, V_MOTION,
                                                         lambda x: (np.abs(x) < BLOCK_HALF).astype(float))), m)
    P0 = looks(synthesize(lit(), g, SHAPE), m)
    LP = looks(synthesize(lit(), g, SHAPE, motion=sway(F_MOTION, V_MOTION, lambda x: (x == 0).astype(float))), m)
    n0 = SHAPE[0] // 2
    inj = look_average(F_MOTION, V_MOTION, R4_CENTRES, W)
    rel = inj - inj[len(inj) // 2]
    out = {'look_times_s': R4_CENTRES, 'injected_absolute': inj, 'injected_relative': rel}
    regions = {'ground': (slice(n0 - 180, n0 + 180), slice(0, SHAPE[1]), L0, LB),
               'point': (slice(n0 - 64, n0 + 64), slice(SHAPE[1] // 2 - 4, SHAPE[1] // 2 + 5), P0, LP)}
    for name, (rows, cols, S, L) in regions.items():
        out[name] = trackers(S, L, g, rows, cols, inj, rel)
    return out


def trackers(S, L, g, rows, cols, inj, rel):
    """The three readings of one region: complex and magnitude registration within one image (moved minus
    motionless, so only the motion's effect remains), and the paired reading against the motionless twin."""
    c0, e0 = velocity_series(S, g, rows, cols), velocity_series(S, g, rows, cols, envelope=True)
    cpx = velocity_series(L, g, rows, cols) - c0
    env = velocity_series(L, g, rows, cols, envelope=True) - e0
    par = paired_velocity(S, L, g, rows, cols)
    gc, sc = response(cpx, rel, c0)
    ge, se = response(env, rel, e0)
    gp, sp = response(par, inj, paired_velocity(S, S, g, rows, cols))
    return {'complex': cpx, 'envelope': env, 'paired': par, 'gain_complex': gc, 'gain_complex_se': sc,
            'gain_envelope': ge, 'gain_envelope_se': se, 'gain_paired': gp,
            'floor_complex_um_s': float(np.std(c0) * 1e6), 'floor_envelope_um_s': float(np.std(e0) * 1e6)}


def part_c(g):
    shape = (1024, 48)
    centres = np.arange(-10.0, 10.001, 0.125)      # 8 looks a second: every tested frequency is below their Nyquist
    m = look_masks(g, shape[0], centres, W)
    rows, cols = slice(shape[0] // 2 - 64, shape[0] // 2 + 64), slice(shape[1] // 2 - 4, shape[1] // 2 + 5)
    scene = lambda: ground(g, shape, np.random.default_rng(13), bright=POINT_AMP)
    L0 = looks(synthesize(scene(), g, shape), m)
    base = velocity_series(L0, g, rows, cols, envelope=True)
    rows_out = []
    for f in LOWPASS_F:
        L = looks(synthesize(scene(), g, shape, motion=sway(f, V_MOTION, lambda x: (x == 0).astype(float))), m)
        rec = velocity_series(L, g, rows, cols, envelope=True) - base
        # fit a cos + b sin + c at the known frequency; the series is relative to the centre look (a constant)
        A = np.stack([np.cos(2 * np.pi * f * centres), np.sin(2 * np.pi * f * centres), np.ones_like(centres)], 1)
        coef, *_ = np.linalg.lstsq(A, rec, rcond=None)
        gain = float(np.hypot(coef[0], coef[1]) / V_MOTION) * np.sign(coef[0])
        rows_out.append({'f_hz': f, 'recovered_gain': gain, 'sinc_gain': float(look_gain(f, W))})
    return {'look_s': W, 'look_centres_s': centres, 'rows': rows_out,
            'r7_frequencies_hz': [1.0, 2.0, 3.0], 'r7_gains': [float(look_gain(f, W)) for f in (1.0, 2.0, 3.0)]}


def part_d(g):
    from sarsim.dwell import DwellProduct
    p = DwellProduct(PRODUCT)
    row, col = p.geolocate(*DESERT_WEST)
    crop, (r0, c0) = p.crop(row, col, 4096, 1536)
    # A: R9's coherence on the real crop, measured and from its power spectrum alone
    deltas = np.array([0.0, 0.3714, 0.7427, 1.1141, 1.4855, 1.8568, 2.2282])
    m = look_masks(g, crop.shape[0], deltas, W)
    L = looks(crop, m)
    real_a = {'delta_s': deltas, 'measured': [float(coherence(L[0], L[k])) for k in range(len(deltas))],
              'from_power_spectrum': [float(coherence_from_power(crop, m[0], m[k])) for k in range(len(deltas))]}
    del L
    # D: block motion on real texture (256 range columns are ample)
    sub = crop[:, 640:896].copy()
    x = (np.arange(sub.shape[0]) - sub.shape[0] // 2) * g.dx
    edge = np.clip((BLOCK_HALF - np.abs(x)) / 1.5 + 0.5, 0, 1)
    w = edge ** 2 * (3 - 2 * edge)
    d = lambda t: (V_MOTION / (2 * np.pi * F_MOTION)) * np.sin(2 * np.pi * F_MOTION * t)
    moved = inject_region_motion(sub, g, w, d)
    mm = look_masks(g, sub.shape[0], R4_CENTRES, W)
    L0, L1 = looks(sub, mm), looks(moved, mm)
    n0 = sub.shape[0] // 2
    rows, cols = slice(n0 - 180, n0 + 180), slice(0, sub.shape[1])
    inj = look_average(F_MOTION, V_MOTION, R4_CENTRES, W)
    rel = inj - inj[len(inj) // 2]
    return {'crop_origin': [int(r0), int(c0)], 'crop_shape': list(crop.shape), 'coherence': real_a,
            'identity_max_error': float(np.max(np.abs(np.array(real_a['measured']) - real_a['from_power_spectrum']))),
            'block': {'injected_absolute': inj, 'injected_relative': rel,
                      **trackers(L0, L1, g, rows, cols, inj, rel)}}


def main():
    g = DwellGeometry.from_record('giza-20250827')
    params = {'geometry': 'giza-20250827', 'look_s': W, 'looks': len(R4_CENTRES), 'image_px': SHAPE,
              'motion': {'f_hz': F_MOTION, 'v_m_s': V_MOTION, 'block_half_m': BLOCK_HALF},
              'bright_amplitude': POINT_AMP, 'lowpass_hz': LOWPASS_F, 'real_crop': 'desert_west (legacy R1/R9 patch)'}
    with Run(RID, 'What one dwell can see of the ground\'s motion', params) as run:
        a, b, c = part_a(g), part_b(g), part_c(g)
        if PRODUCT.is_file():
            d = part_d(g)
            d['origin'] = 'computed from the product in this run'
        else:
            prev = json.loads((RESULTS / RID / 'summary.json').read_text())['real']
            d = {**prev, 'origin': 'carried from the last desktop run (the product is not on this machine)'}
        gr, pt = b['ground'], b['point']
        r7 = ', '.join(f"{abs(r['recovered_gain']) * 100:.0f}" for r in c['rows'] if r['f_hz'] in (1.0, 2.0, 3.0))
        finding = (
            f"One dwell is nearly blind to ordinary ground moving. Registering its looks against each other, as every "
            f"single-image method does, recovers {gr['gain_complex'] * 100:.0f} ± {gr['gain_complex_se'] * 100:.0f}% of a "
            f"20 m stretch of ground's motion by complex correlation and {gr['gain_envelope'] * 100:.0f} ± "
            f"{gr['gain_envelope_se'] * 100:.0f}% by magnitudes, on simulated desert on the real geometry, and "
            f"{d['block']['gain_complex'] * 100:.0f} ± {d['block']['gain_complex_se'] * 100:.0f}% and "
            f"{d['block']['gain_envelope'] * 100:.0f} ± {d['block']['gain_envelope_se'] * 100:.0f}% "
            f"on the real Giza image, though the motion is in the data ({gr['gain_paired'] * 100:.0f}% and "
            f"{d['block']['gain_paired'] * 100:.0f}% read against a motionless twin, which no real dwell has). "
            f"Looks cut from one image share their spectrum sample for sample, so their correlation locks at zero, and "
            f"motion common to a patch of speckle leaves the image's statistics unchanged; magnitudes catch only what is point-like in the ground. A lone bright point's "
            f"envelope follows its motion in full ({pt['gain_envelope'] * 100:.0f} ± {pt['gain_envelope_se'] * 100:.0f}%). The first investigation's coherence "
            f"curve is this identity: simulated desert with no decorrelation reproduces it, and the real crop's power "
            f"spectrum reproduces its own coherence to {d['identity_max_error']:.0e}. A look also averages motion over "
            f"its 1.96 s and passes little above 0.3 Hz: a bright point swaying at R7's 1, 2 and 3 Hz comes back at "
            f"{r7}% "
            f"of its speed.")
        run.save({'coherence': a, 'motion': b, 'lowpass': c, 'real': d, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
