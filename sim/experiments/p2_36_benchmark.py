"""P2-36 · One auditable comparison: what one image could tell of presence, location and shape, and what readers achieve.

    uv run python experiments/p2_36_benchmark.py

A frozen benchmark (BENCHMARK.md), after two independent reviews asked for one fully specified comparison before any
parameter ranges: the same source in every world, a valid reference for each world, the observation and its tails
declared, the scattering and noise model stated, bounds for presence, horizontal location and shape, and implemented
detectors' achieved performance beside the bounds, with what each is told recorded.

Worlds (identical ground, scatterers and source; only the cavity differs):
  W0  no cavity;  WA  a room 6 m on a side, roof 5 m down;  WB  an L-shaped tunnel 2.5 m square, roof 5 m down;
  WA6 the same room 6 m east;  WC  a room 10 m on a side, roof 5 m down (revision 2: the case earlier runs left open).
  Revision 3 changes only the analysis (the worlds' hash is unchanged): see BENCHMARK.md's erratum.
  Pairs: presence (W0-WA, W0-WB, W0-WC), horizontal location (WA-WA6), shape (WA-WB); absolute depth is not scored.
Excitations:
  quiet   Giza's regional microseisms (0.1-0.3 Hz, measured 67 km east), the imprint by the solver settling under a
          uniform strain (P2-28), waves from every direction; the oracle ensemble-averaged over the field (its expected
          squared motion, no peak allowance).
  strong  one vertical point force 15 m west of the site (a lorry), the same force in every world, solved by the
          elastic solver on a box recorded to 70 m; its amplitude scaled so the vertical velocity at the site without a
          cavity is the FTA's truck-over-a-bump level at 15 m; one declared frequency (the rule: the frequency in the
          FTA band where the presence (room) pair's oracle on the measured 70 m disc is largest, found and then frozen);
          every other frequency of the band, sampled every 2 Hz, reported beside it.
Observation: the whole image (the 2025 dwell's 5 km square). The solver's field is used within 70 m; beyond it an
  assumed envelope (a surface wave carried undiminished to the image's edge), a declared assumption, not a bound; the
  results are also given within 70 m and 39 m, and with the envelope larger and attenuated.
Scattering and noise: fully developed speckle (white circular Gaussian scatterers on the image's pixel grid, texture
  power one per cell) and white receiver noise at an SNR per cell that is conditional (sarsim.radiometry): the site's
  median sigma0 measured in the image over ICEYE's specified noise floor, the product carrying none; headline the best
  specified for Dwell Fine (documentation 6.0.8, -23.7 dB), the most favourable to a detector; every other specified
  value (documentation 6.0.0's -18 to -15 dB, 6.0.8's worst), ground at 0 dB and the measured lower bound swept. The Doppler-to-time relation (checked pulse by
  pulse, P2-25).
Layers:
  L1  unknown reflectivity (told the excitation and both worlds' motion, not the texture): Theorem C per along-track
      line (sarsim.finite.per_line_certificate), each line's reference holding its own common motion (the lorry's wave,
      the microseisms' worst-case envelope, world 0's cavity) and the receiver noise, KL summed over lines (independent
      with the range band widened, which can only add information); checked against the exact KL on lines. The quiet
      case averaged over the field's realisations through a bound linear in each line's phase energy (l1_ensemble).
  L2  the oracle told the reflectivity too: E Delta^2 = 8 sum SNR <sin^2(k du / 2)> (independent uniform scatterer
      phases; an upper bound on the average over scenes by Jensen), and the complete coherent echo difference for fixed
      scenes pulse by pulse (P2-29's echo model, complex128) against it.
Implemented detectors:
  D1  the oracle's test on fixed scenes: Gaussian in the receiver noise, its ROC exactly Phi(Phi^-1(alpha) + Delta)
      with Delta from the pulse-by-pulse echo difference.
  D2  the score test (known excitation and difference pattern, unknown texture): predicted from the exact Fisher
      information within 70 m; achieved on synthesised images (one noise law, threshold from calibration grounds, the
      evaluation grounds disjoint from them, every rate with its interval).
  D3  the published method (P2-07's pipeline) on the benchmark's strong-case images: P2-37.
Every allowance is labelled assumed or empirical, and the conclusion is shown under larger ones.
"""
import copy
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import ndtr, ndtri

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run, load, memo
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source
from sarsim import finite
from sarsim import information as inf
from sarsim.acquisition import DwellGeometry
from sarsim.radiometry import snr_per_cell

RID = 'p2_36_benchmark'
SITES = Path(__file__).resolve().parents[2] / 'sites'
_SNR = snr_per_cell('giza-20250827')
FROZEN = {
    'worlds': {
        'W0': [],
        'WA': [[[0.0, 0.0, -8.0], [6.0, 6.0, 6.0]]],
        'WB': [[[-5.0, 0.0, -6.25], [18.0, 2.5, 2.5]], [[2.75, 7.375, -6.25], [2.5, 17.25, 2.5]]],
        'WA6': [[[6.0, 0.0, -8.0], [6.0, 6.0, 6.0]]],
        'WC': [[[0.0, 0.0, -10.0], [10.0, 10.0, 10.0]]],     # v2: the room the earlier runs left open (P2-22, P2-26)
    },
    'pairs': {'presence (room)': ['W0', 'WA'], 'presence (L tunnel)': ['W0', 'WB'], 'location (6 m)': ['WA', 'WA6'],
              'shape (room or L tunnel)': ['WA', 'WB'], 'presence (10 m room)': ['W0', 'WC']},
    'strong_source': {'kind': 'vertical point force', 'east_m': -15.0, 'north_m': 0.0, 'level': 'FTA truck over a bump, 15 m'},
    'solver': {'extent': [[-100.0, 100.0], [-100.0, 100.0], [-100.0, 3.0]], 'h': 1.0, 'pml_m': 20.0, 'half_m': 70.0,
               'record_s': 0.5, 'f_top': 120.0, 'record_every': 4},
    'report_disc_m': 39.0,                    # reported beside the whole image, as the recorded 70 m disc is
    # the SNR per cell is conditional: the site's median sigma0 measured in the image over a specified noise floor
    # (sarsim.radiometry); the headline is the best specified for the product's mode (Dwell Fine, documentation 6.0.8),
    # the most favourable to a detector; every other specified value, and the measured lower bound, in the sweep
    'snr_db_headline': round(_SNR['headline_db'], 2),
    'snr_db_bright': round(_SNR['bright_db'], 2),        # ground at 0 dB (about the brightest natural ground, assumed)
    'snr_db_sweep': sorted({round(_SNR['measured_floor_lower_db'], 2), *[round(sc['snr_db'], 2) for sc in _SNR['scenarios']],
                            round(_SNR['headline_db'], 2), round(_SNR['bright_db'], 2), 30.0, 40.0}),   # 30, 40 dB: bright persistent ground
    'snr_status': _SNR['status'],
    'revision': 'v3, after the eighth review: the noise floor as conditional scenarios by documentation version and mode '
                '(headline Dwell Fine 6.0.8), the score test\'s evaluation grounds disjoint from its calibration grounds, '
                'analysis caches keyed on their inputs, the worst realisation\'s pass-mean over a finite window, the '
                'covariance floor falling back on the receiver noise',
    'acquisition': 'giza-20250827',
    'rock': 'limestone-mokattam',
    'analysis': {
        'f_step_hz': 1.0,                    # the strong case's frequency grid over the FTA band
        'recorded_disc_m': 70.0,             # the solver's field measured within it; the tail enveloped beyond
        'tail_ring_m': [55.0, 70.0],         # the ring whose largest r-weighted energy per azimuth bin sets the envelope
        'azimuth_bins': 36,
        'quiet_disc_m': 58.0,                # P2-28's receivers reach 60 m
        'quiet_ring_m': [46.0, 58.0],
        'quiet_energy_exponent': 2.8,        # |u|^2 ~ r^-2.8 beyond: the slowest decay measured on a 100 m domain (the
                                             # horizontal component, fourth review A3), empirical
        'directions_quiet': 24,
        'near_chance_tv': 0.05, 'target_tv': 0.9,
        'check_prf_hz': 400.0, 'check_realisations': 12, 'check_patch_m': [6.0, 1.5],
        'verify_line_px': 2048, 'verify_amplifications': [1.0, 10.0, 100.0, 1000.0],
        'grid_allowance_amplitude': 0.11,    # P2-26's largest change on halving the grid (surface incidence), empirical
        'tail_multipliers': [1.0, 2.0, 5.0], 'attenuation_q': [50.0, 20.0],
        # the score test's grounds: disjoint streams, asserted (v3: the v2 null evaluation grounds were calibration grounds)
        'mc_amplifications': [1.0, 800.0],
        'mc_seeds': {'noise_law': [35900, 35901], 'sign': [35950, 35951], 'calibration': [36000, 36400],
                     'null_evaluation': [37000, 37100], 'cavity_evaluation': [38000, 38100]},
    },
}
FROZEN_HASH = hashlib.sha1(json.dumps(FROZEN, sort_keys=True).encode()).hexdigest()[:12]
# what defines the images (the worlds, the lorry, the solver, the geometry, the band's grid): the same in every revision
# that changes only the analysis, so the solver's fields, and P2-37's images, carry over
WORLDS_HASH = hashlib.sha1(json.dumps({k: FROZEN[k] for k in ('worlds', 'strong_source', 'solver', 'acquisition', 'rock')}
                                      | {'f_step_hz': FROZEN['analysis']['f_step_hz']}, sort_keys=True).encode()).hexdigest()[:12]
CACHE = RESULTS / 'cache' / 'p2_36'


def module(name, path):
    """Another experiment's module, registered under `name` so numba's cache can find it again."""
    import sys
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ the strong case: one source, every world

def medium_for(boxes):
    site = load_site('bench-void')
    raw = copy.deepcopy(site.raw)
    f0 = raw['features'][0]
    feats = []
    for k, (c, s) in enumerate(boxes):
        f = copy.deepcopy(f0)
        f['id'] = f'part_{k}'
        f['shape'] = dict(f['shape'], centre=list(c), size=list(s))
        f['fill'] = 'air'
        feats.append(f)
    raw['features'] = feats
    sv = FROZEN['solver']
    g = Grid.covering(*[tuple(e) for e in sv['extent']], sv['h'])
    return g, Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))


def pulse(nt, dt, f_top):
    s = 1.0 / (2 * np.pi * f_top)
    t = np.arange(nt) * dt
    return np.exp(-0.5 * ((t - 6 * s) / s) ** 2).astype(np.float32)


def solve(world):
    """The surface's velocity (3 components) on the recorded grid for one world under the point force; kept on disk under
    a hash of everything that defines it."""
    sv = FROZEN['solver']
    key = {'what': 'p2_36 strong', 'boxes': FROZEN['worlds'][world], 'solver': sv, 'source': FROZEN['strong_source']}
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / (hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16] + '.npz')
    if path.is_file():
        z = np.load(path)
        return z['t'], z['traces'], z['xs'], z['ys'], json.loads(str(z['info']))
    g, med = medium_for(FROZEN['worlds'][world])
    h = sv['h']
    sim = Simulation(med, pml_width=int(round(sv['pml_m'] / h)), f0=30.0, pml_vp=3300.0)
    nt = int(np.ceil(sv['record_s'] / sim.dt))
    w = pulse(nt + 1, sim.dt, sv['f_top'])
    pick = lambda a: a[np.abs(a) <= sv['half_m']]
    xs, ys = pick(g.x), pick(g.y)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    top = np.argmax(med.solid, axis=2)
    ix, iy = np.searchsorted(g.x, xs), np.searchsorted(g.y, ys)
    ztop = g.origin[2]
    Z = ztop - top[np.ix_(ix, iy)] * h
    rec = Receivers(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1))
    src = FROZEN['strong_source']
    isx, isy = int(np.argmin(np.abs(g.x - src['east_m']))), int(np.argmin(np.abs(g.y - src['north_m'])))
    zs = ztop - top[isx, isy] * h
    sources = [Source((float(g.x[isx]), float(g.y[isy]), zs), w, 'force', (0.0, 0.0, -1.0))]
    t0 = time.time()
    res = sim.run(sources, rec, nt, record_every=sv['record_every'])
    info = {'runtime_s': round(time.time() - t0, 1), 'steps': nt, 'dt_s': sim.dt, 'cells': int(np.prod(g.shape))}
    np.savez(path, t=res.t, traces=res.traces.astype(np.float32), xs=xs, ys=ys, info=json.dumps(info))
    print(f"  solved {world}: {info['runtime_s']:.0f} s", flush=True)
    return res.t, res.traces.astype(np.float32), xs, ys, info


# ------------------------------------------------------------------ fields on the ground, per world

AN = FROZEN['analysis']


def geometry():
    return DwellGeometry.from_record(FROZEN['acquisition'])


def pass_seconds(g):
    return g.nu_band * g.V / g.Ka


def strong_fields(g):
    """Per world, the slant-range increase [m, complex amplitude] at each frequency of the FTA band on the recorded grid,
    under the one force scaled at each frequency so that the vertical velocity amplitude at the site without a cavity is
    sqrt(2) times the FTA's rms level for a truck over a bump at 15 m: that harmonic held for the whole pass."""
    path = CACHE / f'fields_{WORLDS_HASH}.npz'
    if path.is_file():
        z = np.load(path)
        return {'f': z['f'], 'xs': z['xs'], 'ys': z['ys'], 'scale': z['scale'], 'D': {w: z[w] for w in FROZEN['worlds']},
                'decay': json.loads(str(z['decay']))}
    cul = json.loads((SITES / 'ambient.json').read_text())['cultural']
    band, v_rms = cul['band_hz']['value'], cul['bus_or_truck_over_bump']['value']
    fs = np.arange(band[0], band[1] + 0.5 * AN['f_step_hz'], AN['f_step_hz'])
    los = np.asarray(g.los_enu)
    D, vz, decay, Wf = {}, {}, {}, None
    for w in FROZEN['worlds']:
        t, tr, xs, ys, info = solve(w)
        if Wf is None:
            n_src = info['steps'] + 1
            tw = np.arange(n_src) * info['dt_s']
            wl = pulse(n_src, info['dt_s'], FROZEN['solver']['f_top']).astype(np.float64)
            Wf = (wl @ np.exp(-2j * np.pi * np.outer(tw, fs))) * info['dt_s']
        E = np.exp(-2j * np.pi * np.outer(t, fs)) * (t[1] - t[0])
        V = (tr.reshape(-1, len(t)).astype(np.float64) @ E).reshape(tr.shape[0], 3, len(fs)) / Wf
        e = np.sum(tr.astype(np.float64) ** 2, axis=(0, 1))
        decay[w] = float(e[int(0.9 * len(t)):].sum() / e.sum())
        D[w] = (-np.einsum('rcf,c->fr', V, los) / (2j * np.pi * fs[:, None])).reshape(len(fs), len(xs), len(ys))
        vz[w] = V[:, 2, :]
    i0 = int(np.argmin(np.abs(xs))) * len(ys) + int(np.argmin(np.abs(ys)))
    scale = np.sqrt(2) * v_rms / np.abs(vz['W0'][i0])
    for w in D:
        D[w] = D[w] * scale[:, None, None]
    np.savez(path, f=fs, xs=xs, ys=ys, scale=scale, decay=json.dumps(decay), **D)
    return {'f': fs, 'xs': xs, 'ys': ys, 'scale': scale, 'D': D, 'decay': decay}


def quiet_level():
    """The regional microseisms as P2-04 and P2-25 take them: strain amplitude hv sqrt(2) V / c, the power split evenly
    over the directions with independent phases."""
    p204 = load('p2_04_chamber_imprint')
    row = p204['cases'][0]
    return {'case': row['case'], 'rms_velocity_m_s': row['vertical_velocity_m_s'], 'f_hz': row['frequency_hz'],
            'phase_speed_m_s': row['phase_speed_m_s'], 'hv': p204['hv'],
            'strain_amplitude': p204['hv'] * np.sqrt(2) * row['vertical_velocity_m_s'] / row['phase_speed_m_s']}


def quiet_energy(K, g, eps0):
    """e(x) = eps0^2 <M_phi(x)^2> over the directions, so that the expected pass-mean squared line-of-sight motion is
    e / 2: M_phi the imprint's slant-range increase per unit strain of a wave towards phi."""
    xy = K['xy']
    xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
    ix, iy = np.searchsorted(xs, xy[:, 0]), np.searchsorted(ys, xy[:, 1])
    los = np.asarray(g.los_enu)
    az = np.deg2rad(np.arange(AN['directions_quiet']) * 180.0 / AN['directions_quiet'])
    acc = np.zeros(len(xy))
    maps = []
    for p in az:
        m = -(K['Kxx'] * np.sin(p) ** 2 + K['Kyy'] * np.cos(p) ** 2 + K['Kxy'] * np.sin(p) * np.cos(p)) @ los
        acc += m * m
        M = np.zeros((len(xs), len(ys)))
        M[ix, iy] = m
        maps.append(M)
    e = np.zeros((len(xs), len(ys)))
    e[ix, iy] = eps0 ** 2 * acc / len(az)
    return xs, ys, e, maps


def kernel_difference(Ka, Kb):
    if Ka is None:
        return Kb
    if Kb is None:
        return {'xy': Ka['xy'], **{k: -Ka[k] for k in ('Kxx', 'Kyy', 'Kxy')}}
    return {'xy': Kb['xy'], **{k: Kb[k] - Ka[k] for k in ('Kxx', 'Kyy', 'Kxy')}}


# ------------------------------------------------------------------ envelopes beyond what was measured

def ring_envelope(e, xs, ys, ring, r_ref, exponent, centre=(0.0, 0.0), bins=None):
    """A(theta) such that e <= A(theta) (r_ref / r)^exponent beyond r_ref, from the largest r-weighted value in each
    azimuth bin of the ring (each bin also takes its neighbours' largest, so a lobe between bins is not missed)."""
    bins = bins or AN['azimuth_bins']
    X, Y = np.meshgrid(xs - centre[0], ys - centre[1], indexing='ij')
    R, T = np.hypot(X, Y), np.arctan2(Y, X)
    sel = (R >= ring[0]) & (R <= ring[1])
    b = ((T[sel] + np.pi) / (2 * np.pi) * bins).astype(int) % bins
    v = e[sel] * (R[sel] / r_ref) ** exponent
    A = np.zeros(bins)
    np.maximum.at(A, b, v)
    return np.maximum(A, np.maximum(np.roll(A, 1), np.roll(A, -1)))


def envelope_fn(A, r_ref, exponent, centre=(0.0, 0.0), atten=None):
    bins = len(A)

    def f(E, N):
        X, Y = E - centre[0], N - centre[1]
        R = np.maximum(np.hypot(X, Y), 1e-9)
        b = ((np.arctan2(Y, X) + np.pi) / (2 * np.pi) * bins).astype(int) % bins
        out = A[b] * (r_ref / R) ** exponent
        if atten is not None:
            out = out * np.exp(-atten * np.maximum(R - r_ref, 0.0))
        return out
    return f


def sigma_profile(e, xs, ys, step=5.0, r_max=70.0):
    """Energy per metre of radius on annuli: constant for an unattenuated surface wave, falling for body waves."""
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    R = np.hypot(X, Y)
    dA = float(xs[1] - xs[0]) * float(ys[1] - ys[0])
    edges = np.arange(10.0, r_max + 1e-9, step)
    prof = [float(e[(R >= a) & (R < a + step)].sum() * dA / step) for a in edges[:-1]]
    mid = edges[:-1] + step / 2
    sel = mid >= 40
    slope = float(np.polyfit(np.log(mid[sel]), np.log(prof)[sel], 1)[0]) if sel.sum() > 1 else None
    return {'r_mid_m': mid.tolist(), 'sigma_m3': prof, 'slope_40_70': slope,
            'outer_over_40m': float(prof[-1] / prof[int(np.argmin(np.abs(mid - 42.5)))])}


# ------------------------------------------------------------------ the image's lines along track

def image_frame(g):
    src = json.loads((SITES / 'acquisitions' / f"{FROZEN['acquisition']}.json").read_text())['source']
    n_along, n_lines = src['shape']
    dg = g.dr / np.sin(g.theta)
    rho = (np.arange(n_lines) - n_lines / 2 + 0.5) * dg
    in_a, _, _ = inf.band_masks(g, (n_along, 1))
    return {'rho': rho, 'dg': dg, 'half_along_m': n_along * g.dx / 2, 'n_along': int(n_along),
            'n_band_along': int(in_a.sum()), 'n_lines': int(n_lines)}


def line_integrals(g, frame, emap, xs, ys, r_in, env, centre=(0.0, 0.0)):
    """S[l] = integral along image line l of e (m^3): the measured map within r_in of `centre`, the envelope beyond, to
    the image's edge. Trapezoids: for lines within 200 m of the site, 0.25 m steps within 200 m along track and 2 m
    beyond; for lines farther out, where the envelope varies on the scale of its distance, 4 m steps."""
    a_hat, r_hat = np.asarray(g.along_track_en), np.asarray(g.ground_range_en)
    H = frame['half_along_m']

    def grid(parts):
        a = np.unique(np.concatenate(parts))
        w = np.zeros_like(a)
        d = np.diff(a)
        w[:-1] += d / 2
        w[1:] += d / 2
        return a, w
    coarse = np.arange(200.0, H, 2.0)
    near_grid = grid([-coarse, np.arange(-200.0, 200.0 + 1e-9, 0.25), coarse, [-H, H]])
    far_grid = grid([np.arange(-H, H, 4.0), [H]])
    interp = RegularGridInterpolator((xs, ys), emap, bounds_error=False, fill_value=0.0) if emap is not None else None
    rho_all = frame['rho']
    S = np.zeros(len(rho_all))
    for sel, (a, wts) in ((np.abs(rho_all) <= 200.0, near_grid), (np.abs(rho_all) > 200.0, far_grid)):
        idx = np.nonzero(sel)[0]
        for i0 in range(0, len(idx), 400):
            ii = idx[i0:i0 + 400]
            rho = rho_all[ii]
            E = a[None, :] * a_hat[0] + rho[:, None] * r_hat[0]
            N = a[None, :] * a_hat[1] + rho[:, None] * r_hat[1]
            R = np.hypot(E - centre[0], N - centre[1])
            inside = R <= r_in
            if env is None:
                val = np.zeros_like(E)
            else:
                val = env(E, N)
                val[inside] = 0.0
            if interp is not None and inside.any():
                val[inside] = interp(np.stack([E[inside], N[inside]], 1))
            S[ii] = val @ wts
    return S


def disc_integral(emap, xs, ys, r):
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    return float(emap[np.hypot(X, Y) <= r].sum() * (xs[1] - xs[0]) * (ys[1] - ys[0]))


# ------------------------------------------------------------------ L1: the per-line certificate

def phase_energy(g, S, osc):
    """Per line, sum over its pixels of (n_b / n) <4 sin^2(phi / 2)> <= band_frac / dx (k0^2 / 2)(1 + osc) S: phi = k0 x
    the slant-range change, e / 2 its pass-mean square times (1 + osc), osc bounding the pass-mean of the harmonic's
    oscillating part."""
    k0 = 4 * np.pi / g.lam
    return g.band_frac / g.dx * k0 ** 2 / 2 * (1 + osc) * np.asarray(S)


def l1_certificate(Q, Qc, sigma2, amplify=1.0):
    """Theorem C per line, KL summed over the lines (sarsim.finite.per_line_certificate): q^2 = Q / floor, floor =
    sigma^2 + (1 - sqrt(Qc))^2 (the line's own common motion in its reference, and its receiver noise), KL <= rho^2 /
    (2 (1 - rho)), rho = 2 q + q^2. `amplify` scales the differential motion's amplitude (Q by its square)."""
    Q = np.asarray(Q, float) * amplify ** 2
    c = finite.per_line_certificate(Q, Qc, sigma2)
    c['lines_touched'] = int(np.sum(Q > 1e-6 * Q.max())) if Q.max() > 0 else 0
    return c


def growth(tv_of, target):
    """The amplitude factor at which tv_of(factor) reaches target (bisection in log; monotone)."""
    lo, hi = 1e-6, 1e9
    if tv_of(hi) < target:
        return None
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if tv_of(mid) >= target:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1 + 1e-6:
            break
    return float(hi)


def state(tv):
    return 'near chance' if tv < AN['near_chance_tv'] else ('95% at 5% excluded' if tv < AN['target_tv'] else 'unresolved')


# ------------------------------------------------------------------ L2: the oracle, ensemble and fixed scenes

def oracle(g, integral_e, snr_db, osc):
    """E Delta^2 = 8 sum_j SNR_j <sin^2(k0 d_j / 2)> <= snr cells k0^2 (1 + osc) integral e dA (P2-29), snr the SNR per
    cell over the processed band's share; TV <= 2 Phi(sqrt(E Delta^2) / 2) - 1 (Jensen: an upper bound on the average
    over scenes of the exact conditional TV)."""
    k0 = 4 * np.pi / g.lam
    snr = 10 ** (snr_db / 10) / g.band_frac
    d2 = snr * inf.cells_per_m2(g) * k0 ** 2 * (1 + osc) * integral_e
    tv = float(2 * ndtr(np.sqrt(d2) / 2) - 1)
    return {'e_delta2': float(d2), 'tv_upper': tv, 'found_at_5pct_upper': min(1.0, 0.05 + tv)}


def oracle_d2_target():
    return float(2 * ndtri((1 + AN['target_tv']) / 2)) ** 2


def interp_complex(D, xs, ys, E, N):
    f = lambda Z: RegularGridInterpolator((xs, ys), Z, bounds_error=False, fill_value=0.0)(np.stack([E, N], -1))
    return f(D.real) + 1j * f(D.imag)


def fixed_scenes(g, D0, D1, xs, ys, f_hz):
    """The complete coherent echo difference for fixed scenes, pulse by pulse (P2-29's echo model, complex128): patches of
    scatterers at full density where the pair differs most, each world's whole motion (the lorry's wave and the cavity's
    scattering) on every scatterer, |mu_1 - mu_0|^2 exactly, against the ensemble formula sum_j E_j <4 sin^2(k0 dd_j / 2)>
    (independent uniform phases); with the common wave as simulated and ten times larger."""
    from sarsim.echo import EchoSetup, window
    m29 = module('p229', Path(__file__).with_name('p2_29_oracle.py'))
    st = EchoSetup.from_geometry(g, prf=AN['check_prf_hz'])
    t = st.times
    k0 = 4 * np.pi / g.lam
    Lu, Lr = AN['check_patch_m']
    n = 2 * int((Lu / g.resolution) * (Lr / st.rho))
    r_lo, dr_s, n_s = window(st, t, 3.0, Lu, Lr)
    half = int(np.ceil(4 * st.rho / dr_s))
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    dd = np.abs(D1 - D0) * (np.hypot(X, Y) <= 30.0)
    i, j = np.unravel_index(np.argmax(dd), dd.shape)
    c = (float(xs[i]), float(ys[j]))
    a_hat, r_hat = np.asarray(g.along_track_en), np.asarray(g.ground_range_en)
    rows = []
    for scale_common in (1.0, 10.0):
        for rep in range(AN['check_realisations']):
            rng = np.random.default_rng(900 + rep)
            u, r = rng.uniform(-Lu / 2, Lu / 2, n), rng.uniform(-Lr / 2, Lr / 2, n)
            amp = np.exp(2j * np.pi * rng.random(n)) * rng.rayleigh(1 / np.sqrt(2), n) * np.sqrt(2)
            psi = rng.uniform(0, 2 * np.pi)
            ga, gr = u * g.V / st.V_eff, r / np.sin(g.theta)
            E, N = c[0] + ga * a_hat[0] + gr * r_hat[0], c[1] + ga * a_hat[1] + gr * r_hat[1]
            z0, z1 = interp_complex(D0, xs, ys, E, N), interp_complex(D1, xs, ys, E, N)
            ph = np.exp(1j * (2 * np.pi * f_hz * t + psi))[:, None]
            d0 = scale_common * np.real(z0[None, :] * ph)
            d1 = d0 + np.real((z1 - z0)[None, :] * ph)
            ex, inc = m29._echo_pair(t, st.V_eff, st.R0, u, r, amp.real.copy(), amp.imag.copy(), d0, d1, st.lam, st.rho,
                                     r_lo, dr_s, n_s, half)
            Ej = m29._energy(t, st.V_eff, st.R0, u, r, st.rho, r_lo, dr_s, n_s, half) * np.abs(amp) ** 2
            formula = float(np.sum(Ej * np.mean(4 * np.sin(k0 * (d1 - d0) / 2) ** 2, axis=0)))
            rows.append({'common_scale': scale_common, 'realisation': rep, 'exact_over_formula': ex / formula,
                         'no_cross_terms_over_formula': inc / formula})
    summ = {}
    for sc in (1.0, 10.0):
        v = np.array([x['exact_over_formula'] for x in rows if x['common_scale'] == sc])
        summ[f'common_x{sc:g}'] = {'mean': float(v.mean()), 'se': float(v.std(ddof=1) / np.sqrt(len(v))),
                                   'spread': float(v.std(ddof=1)), 'min': float(v.min()), 'max': float(v.max())}
    return {'patch_centre_site_m': c, 'patch_m': [Lu, Lr], 'scatterers': n, 'pulses': len(t), 'prf_hz': st.prf,
            'f_hz': f_hz, 'rows': rows, 'summary': summ}


def verify_line(g, D0, D1, xs, ys, f_hz, sigma2, micro_phase):
    """The per-line certificate against the exact KL on one line of the model (cyclic, VERIFY px, through the pair's
    largest difference within 20 m), each world's whole motion and the regional microseisms' envelope as a uniform phase
    in both covariances, white receiver noise: KL(CN(0, C1) || CN(0, C0)) from the eigenvalues of the whitened change
    (finite_complex_kl), at the real level and with the difference amplified."""
    n = AN['verify_line_px']
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    dd = np.abs(D1 - D0) * (np.hypot(X, Y) <= 20.0)
    i, j = np.unravel_index(np.argmax(dd), dd.shape)
    a_hat = np.asarray(g.along_track_en)
    a = (np.arange(n) - n / 2) * g.dx
    E, N = xs[i] + a * a_hat[0], ys[j] + a * a_hat[1]
    z0, z1 = interp_complex(D0, xs, ys, E, N), interp_complex(D1, xs, ys, E, N)
    in_a, _, nu = inf.band_masks(g, (n, 1))
    tb = g.nu_to_time(nu)[in_a]
    fb = np.fft.fftfreq(n)[in_a]
    k0 = 4 * np.pi / g.lam
    Emat = np.exp(-2j * np.pi * np.outer(fb, np.arange(n))) / np.sqrt(n)
    ph = np.exp(1j * 2 * np.pi * f_hz * (tb[:, None] + a[None, :] / g.V))
    phi0 = k0 * np.real(z0[None, :] * ph) + micro_phase
    dphi = k0 * np.real((z1 - z0)[None, :] * ph)
    B0 = Emat * np.exp(-1j * phi0)
    Qc = float(np.sum(np.abs(Emat * np.expm1(-1j * phi0)) ** 2))
    C0 = B0 @ B0.conj().T + sigma2 * np.eye(len(fb))
    out = []
    for amp in AN['verify_amplifications']:
        Dm = B0 * np.expm1(-1j * amp * dphi)
        dC = finite.stable_covariance_difference(B0, Dm)
        ex = finite.finite_complex_kl(C0, dC)
        L = np.linalg.cholesky(C0)
        from scipy.linalg import solve_triangular
        q_true = float(np.linalg.norm(solve_triangular(L, Dm, lower=True)))
        Q = float(np.sum(np.abs(Dm) ** 2))
        cert = l1_certificate([Q], [Qc], sigma2)
        out.append({'amplification': amp, 'kl_exact': ex['kl'], 'kl_certificate': cert['kl_upper'],
                    'exact_over_certificate': ex['kl'] / cert['kl_upper'] if cert['kl_upper'] else None,
                    'q_true': q_true, 'q_bound': float(np.sqrt(Q / finite.line_floor(Qc, sigma2))), 'Q': Q})
    return {'line_px': n, 'line_m': n * g.dx, 'through_site_m': [float(xs[i]), float(ys[j])], 'band_bins': int(len(fb)),
            'common_phase_energy': Qc, 'floor': float(sigma2 + (1 - np.sqrt(Qc)) ** 2), 'rows': out}


def achieved_score_test(g, dD, xs, ys, f_hz, r_max, taper_m=10.0):
    """What the locally most powerful test told the differential pattern and its time (not the texture) reaches on the
    real image's band: its deflection d = sqrt F, F the exact Fisher information on the image's pixels (fisher_grid,
    noise-free, the reference taken as white: the common wave changes F only at second order), within r_max, tapered to
    zero over the last taper_m; TV ~ 2 Phi(d / 2) - 1 for a weak signal. Achieved by an implemented statistic
    (sarsim.information.score, checked against the exact test in P2-25), not a bound."""
    from scipy.ndimage import map_coordinates
    dg = g.dr / np.sin(g.theta)
    na, nr = int(2 * r_max / g.dx), int(2 * r_max / dg)
    A_ax, R_ax = (np.arange(na) - na / 2) * g.dx, (np.arange(nr) - nr / 2) * dg
    a_hat, r_hat = np.asarray(g.along_track_en), np.asarray(g.ground_range_en)
    Aa, Rr = np.meshgrid(A_ax, R_ax, indexing='ij')
    E, N = Aa * a_hat[0] + Rr * r_hat[0], Aa * a_hat[1] + Rr * r_hat[1]
    ci, cj = (E - xs[0]) / (xs[1] - xs[0]), (N - ys[0]) / (ys[1] - ys[0])
    K = map_coordinates(dD.real, [ci, cj], order=3, mode='constant') + 1j * map_coordinates(dD.imag, [ci, cj], order=3,
                                                                                           mode='constant')
    K *= np.clip((r_max - np.hypot(E, N)) / taper_m, 0, 1)
    A, B = inf.fisher_grid(K, g, f_hz)
    tv = lambda F: float(2 * ndtr(np.sqrt(F) / 2) - 1)
    return {'fisher_mean_over_phase': A, 'fisher_swing': B, 'deflection_mean': float(np.sqrt(A)),
            'tv_achieved_mean_phase': tv(A), 'tv_achieved_best_phase': tv(A + B), 'r_max_m': r_max, 'taper_m': taper_m}


def seed_streams(seeds):
    """The score test's grounds as named, half-open seed ranges [a, b); asserted pairwise disjoint, so no ground serves
    two roles (the v2 run drew its null evaluation grounds from inside its calibration range)."""
    rng = {k: range(a, b) for k, (a, b) in seeds.items()}
    names = list(rng)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            assert not set(rng[a]) & set(rng[b]), f'seed streams {a} and {b} overlap'
    return rng


def binomial_interval(k, n, level=0.95):
    """Clopper-Pearson interval for a rate k / n."""
    from scipy.stats import beta
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k))
    return [lo, hi]


def auc_se(auc, n1, n0):
    """Hanley and McNeil's standard error of an AUC from n1 positives and n0 negatives."""
    q1, q2 = auc / (2 - auc), 2 * auc * auc / (1 + auc)
    return float(np.sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc * auc) + (n0 - 1) * (q2 - auc * auc)) / (n1 * n0)))


def score_mc(g, D0, D1, xs, ys, f_hz, sigma2, amps=(1.0, 800.0), seeds=None, shape=(1024, 64)):
    """D2 implemented: the score statistic for the pair's difference, told its pattern and timing but not the speckle,
    on synthesised images (sarsim.synthesize, complex128) of a cyclic patch about the largest difference: scatterers on
    the pixel grid (white circular Gaussian), world 0 moving with its whole motion (the lorry's wave and its own cavity),
    world 1 the same plus the difference amplified `amp` times. One noise law for every image (white, its variance fixed
    once from a separate ground's motionless image at the headline SNR); the 5% threshold set on the calibration grounds
    of world 0; the evaluation on grounds disjoint from them and from each other (seed_streams asserts it), world 0's and
    world 1's, so the AUC compares independent images and the detection rate is read at a threshold not fitted to them;
    each achieved rate with its Clopper-Pearson interval and the AUC with its standard error. Against the weak-signal
    prediction from the exact Fisher information on the same patch (fisher_grid) and the per-line certificate there."""
    from sarsim import synthesize
    from sarsim.scene import Scatterers
    m31 = module('p231', Path(__file__).with_name('p2_31_shape_test.py'))
    Nx, Nr = shape
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    dd = np.abs(D1 - D0) * (np.hypot(X, Y) <= 20.0)
    i, j = np.unravel_index(np.argmax(dd), dd.shape)
    c = np.array([xs[i], ys[j]])
    a_hat, r_hat = np.asarray(g.along_track_en), np.asarray(g.ground_range_en)
    k0 = 4 * np.pi / g.lam
    xa = (np.arange(Nx) - Nx // 2) * g.dx
    yr = (np.arange(Nr) - Nr // 2) * g.dr / np.sin(g.theta)
    XA, YR = np.meshgrid(xa, yr, indexing='ij')
    site = lambda x, y: (c[0] + x * a_hat[0] + y * r_hat[0], c[1] + x * a_hat[1] + y * r_hat[1])
    Ep, Np = site(XA, YR)
    z0p, z1p = interp_complex(D0, xs, ys, Ep, Np), interp_complex(D1, xs, ys, Ep, Np)
    dK1 = z1p - z0p                                           # the unit difference; the statistic's template
    ph_x = np.exp(1j * 2 * np.pi * f_hz * XA / g.V)
    taus = [lambda t: np.cos(2 * np.pi * f_hz * t), lambda t: np.sin(2 * np.pi * f_hz * t)]
    n = Nx * Nr
    zero = np.zeros(n)

    def ground(sd):
        rng = np.random.default_rng(sd)
        sc = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
        scat = Scatterers(x=XA.ravel(), y=YR.ravel(), z=zero, amp=np.abs(sc), phase=np.angle(sc), iso=np.ones(n),
                          flash=zero, nu0=zero, sig_nu=np.ones(n), vib_amp=zero, vib_freq=zero, vib_phase=zero,
                          label=np.zeros(n, int))
        psi = float(rng.uniform(0, 2 * np.pi))
        noise = (rng.normal(size=(Nx, Nr)) + 1j * rng.normal(size=(Nx, Nr))) / np.sqrt(2)
        return scat, psi, noise

    # one noise law: its variance from a separate ground's motionless image
    S = seed_streams(seeds)
    scat_c, _, _ = ground(S['noise_law'][0])
    p_sig = float(np.mean(np.abs(synthesize(scat_c, g, shape, dtype=np.complex128)) ** 2))
    noise_sd = float(np.sqrt(sigma2 * p_sig))

    def statistic(sd, amp):
        scat, psi, noise = ground(sd)

        def motion(x, y, z, t):
            E, N = site(x, y)
            z0 = interp_complex(D0, xs, ys, E, N)
            zz = z0 + amp * (interp_complex(D1, xs, ys, E, N) - z0)
            return np.real(zz[None, :] * np.exp(1j * (2 * np.pi * f_hz * (t[:, None] + x[None, :] / g.V) + psi)))
        img = synthesize(scat, g, shape, motion=motion, dtype=np.complex128) + noise_sd * noise
        q = m31.score_fields(img, g, taus)
        ps = np.exp(1j * psi)
        K1, K2 = np.real(dK1 * ph_x * ps), -np.imag(dK1 * ph_x * ps)     # d = K1 cos(2 pi f t) + K2 sin(2 pi f t)
        return float(k0 * np.sum(K1 * q[0] + K2 * q[1]))
    cal = [statistic(sd, 0.0) for sd in S['calibration']]
    ev0 = [statistic(sd, 0.0) for sd in S['null_evaluation']]
    # the statistic's sign (the model's phase convention) from one separate ground, never from evaluation images
    sd_sign = S['sign'][0]
    sign = 1.0 if statistic(sd_sign, max(amps)) >= statistic(sd_sign, 0.0) else -1.0
    rows = []
    for amp in amps:
        t_cal, t0 = sign * np.array(cal), sign * np.array(ev0)
        t1 = sign * np.array([statistic(sd, amp) for sd in S['cavity_evaluation']])
        thr = float(np.quantile(t_cal, 0.95))
        auc = float(np.mean(t1[:, None] > t0[None, :]))
        A_f, _ = inf.fisher_grid(amp * dK1, g, f_hz, pad=(Nx, Nr))
        d = float(np.sqrt(A_f))
        osc = 1 / (2 * np.pi * f_hz * pass_seconds(g))
        Q = g.band_frac * k0 ** 2 / 2 * (1 + osc) * np.sum(np.abs(amp * dK1) ** 2, axis=0)
        Qc = g.band_frac * k0 ** 2 / 2 * (1 + osc) * np.sum(np.abs(z0p) ** 2, axis=0)
        cert = l1_certificate(Q, Qc, sigma2)
        from scipy.stats import fisher_exact, mannwhitneyu
        k1, k0, n1, n0 = int(np.sum(t1 > thr)), int(np.sum(t0 > thr)), len(t1), len(t0)
        rows.append({'amplification': amp, 'auc_achieved': auc, 'auc_predicted': float(ndtr(d / np.sqrt(2))),
                     'auc_se': auc_se(auc, n1, n0),
                     'auc_p_one_sided': float(mannwhitneyu(t1, t0, alternative='greater').pvalue),
                     'found_at_threshold_achieved': k1 / n1, 'found_count': [k1, n1],
                     'found_interval_95': binomial_interval(k1, n1),
                     'false_alarms_at_threshold_achieved': k0 / n0, 'false_alarm_count': [k0, n0],
                     'false_alarm_interval_95': binomial_interval(k0, n0),
                     'found_over_false_alarms_p_one_sided': float(fisher_exact([[k1, n1 - k1], [k0, n0 - k0]],
                                                                               alternative='greater')[1]),
                     'found_at_005_predicted': float(ndtr(ndtri(0.05) + d)), 'deflection_predicted': d,
                     'tv_predicted': float(2 * ndtr(d / 2) - 1), 'certificate_tv_upper_patch': cert['tv_upper'],
                     'sign': sign, 'mean_shift_over_sd': float((t1.mean() - t0.mean()) / t0.std(ddof=1))})
        print(f"    score MC x{amp:g}: AUC {auc:.3f} (predicted {rows[-1]['auc_predicted']:.3f}), found "
              f"{rows[-1]['found_at_threshold_achieved']:.2f} at a threshold giving {rows[-1]['false_alarms_at_threshold_achieved']:.2f} "
              f"false alarms (predicted {rows[-1]['found_at_005_predicted']:.2f} at 0.05)", flush=True)
    return {'patch_px': list(shape), 'patch_m': [Nx * g.dx, Nr * g.dr / np.sin(g.theta)], 'centre_site_m': c.tolist(),
            'calibration_grounds': len(S['calibration']), 'evaluation_grounds_each': len(S['null_evaluation']),
            'cavity_evaluation_grounds': len(S['cavity_evaluation']), 'seeds': seeds, 'noise_sd': noise_sd, 'rows': rows}

# ------------------------------------------------------------------ the quiet case: each world settling under strain

def host_moduli():
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    rock = next(m for m in rows_ if m['id'] == FROZEN['rock'])
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    vp, vs, rho = val('vp'), val('vs'), val('rho')
    return rho * vp ** 2 - 2 * rho * vs ** 2, rho * vs ** 2


def static_kernels(world):
    """The world's imprint on the surface per unit horizontal strain (east, north, up per eps_ee, eps_nn, eps_en), less
    the same ground without a cavity: P2-28's settling (P2-31's layouts), cached with P2-28's kernels."""
    m28 = module('p228', Path(__file__).with_name('p2_28_two_depths.py'))
    m31 = module('p231', Path(__file__).with_name('p2_31_shape_test.py'))
    boxes = [(tuple(c), tuple(s)) for c, s in FROZEN['worlds'][world]]
    if not boxes:
        return None
    K, info = m31.layout_kernels(m28, f'bench36_{world}', boxes, host_moduli())
    return K


# ------------------------------------------------------------------ one pair, one excitation: every layer

def lines_for(g, frame, emap, xs, ys, r_in, ring, exponent, osc, mult=1.0, atten=None):
    A = ring_envelope(emap, xs, ys, ring, r_in, exponent)
    env = envelope_fn(A * mult, r_in, exponent, atten=atten)
    S = line_integrals(g, frame, emap, xs, ys, r_in, env)
    S_disc = line_integrals(g, frame, emap, xs, ys, r_in, None)
    return {'Q': phase_energy(g, S, osc), 'Q_disc': phase_energy(g, S_disc, osc),
            'area_integral_m4': float(S.sum() * frame['dg']),
            'disc_integral_m4': disc_integral(emap, xs, ys, r_in),
            'obs_disc_integral_m4': disc_integral(emap, xs, ys, FROZEN['report_disc_m'])}


def snr_key(v):
    return f'{v:g}'


def layers(g, frame, Q, Qc, integ, osc, extra=None):
    """L1 at every SNR of the sweep (the noise only raises the floor), the oracle at every SNR, their margins."""
    out = {'l1': {}, 'oracle': {}}
    for snr in FROZEN['snr_db_sweep']:
        s2 = 10 ** (-snr / 10)
        c = l1_certificate(Q, Qc, s2)
        c['state'] = state(c['tv_upper'])
        out['l1'][snr_key(snr)] = c
        o = oracle(g, integ['area_integral_m4'], snr, osc)
        o['state'] = state(o['tv_upper'])
        out['oracle'][snr_key(snr)] = o
    nom = FROZEN['snr_db_headline']
    s2 = 10 ** (-nom / 10)
    out['l1_growth_to_target'] = growth(lambda a: l1_certificate(Q, Qc, s2, a)['tv_upper'], AN['target_tv'])
    on = out['oracle'][snr_key(nom)]
    out['oracle_growth_to_target'] = float(np.sqrt(oracle_d2_target() / on['e_delta2'])) if on['e_delta2'] > 0 else None
    # where the oracle's upper bound stops excluding 95% at 5%: not a demonstration that any detector gets there
    out['oracle_bound_stops_excluding_at_snr_db'] = (nom + 10 * np.log10(oracle_d2_target() / on['e_delta2'])
                                                     if on['e_delta2'] > 0 else None)
    out['l1_disc'] = l1_certificate(integ['Q_disc'], Qc, s2)
    out['oracle_disc'] = oracle(g, integ['disc_integral_m4'], nom, osc)
    out['oracle_report_disc'] = oracle(g, integ['obs_disc_integral_m4'], nom, osc)
    out['tail_share_of_integral'] = 1 - integ['disc_integral_m4'] / integ['area_integral_m4']
    return out


def l1_ensemble(Q_mean, Q_max, Qc_worst, sigma2):
    """The quiet case averaged over the excitation's realisations, through a bound linear in each line's phase energy.
    For every realisation with Q_l <= Q_max_l and each line's common motion within its worst-case envelope (floor
    lambda_l from line_floor), rho = q (2 + q) <= q (2 + q_max) and 1 - rho >= 1 - rho_max, so
    KL_l <= c_l Q_l / lambda_l with c_l = (2 + q_max)^2 / (2 (1 - rho_max)): linear in Q_l. Averaging, E KL <=
    sum_l c_l E[Q_l] / lambda_l, and E TV <= min(sqrt(E KL / 2), sqrt(1 - exp(-E KL))) since both are concave in KL.
    The worst realisation's own certificate (at Q_max) is returned beside it."""
    Qc_worst = np.broadcast_to(np.asarray(Qc_worst, float), np.shape(Q_mean))
    lam = np.array([finite.line_floor(qc, sigma2) or 0.0 for qc in np.ravel(Qc_worst)]).reshape(np.shape(Q_mean))
    if np.any(lam <= 0):
        return {'kl_upper': None, 'tv_upper': 1.0, 'informative': False}
    qm = np.sqrt(np.asarray(Q_max) / lam)
    rho_m = qm * (2 + qm)
    if np.any(rho_m >= 1):
        return {'kl_upper': None, 'tv_upper': 1.0, 'informative': False}
    c = (2 + qm) ** 2 / (2 * (1 - rho_m))
    kl = float(np.sum(c * np.asarray(Q_mean) / lam))
    tv = min(1.0, np.sqrt(kl / 2), np.sqrt(-np.expm1(-kl)))
    worst = finite.per_line_certificate(Q_max, Qc_worst, sigma2)
    return {'kl_upper': kl, 'tv_upper': float(tv), 'c_max': float(c.max()), 'floor_min': float(lam.min()),
            'worst_realisation_tv_upper': worst['tv_upper'], 'state': state(tv)}


_S0_CACHE = {}


def strong_pair(g, frame, F, fi, w0, w1, Qc_micro, quiet_Q=None, mult=1.0, atten=None, grid=1.0):
    """The strong case for one pair at one frequency: the differential map and the world-0 map, per line."""
    f = float(F['f'][fi])
    osc = 1 / (2 * np.pi * f * pass_seconds(g))
    xs, ys = F['xs'], F['ys']
    D0, D1 = F['D'][w0][fi], F['D'][w1][fi]
    R = AN['recorded_disc_m']
    e = np.abs(D1 - D0) ** 2 * grid ** 2
    diff = lines_for(g, frame, e, xs, ys, R, AN['tail_ring_m'], 1.0, osc, mult, atten)
    key = (w0, fi, mult, atten, grid)
    if key not in _S0_CACHE:                   # world 0's own motion, shared by every pair with the same world 0
        e0 = np.abs(D0) ** 2 * grid ** 2
        src = (FROZEN['strong_source']['east_m'], FROZEN['strong_source']['north_m'])
        A0 = ring_envelope(e0, xs, ys, (40.0, 55.0), 40.0, 1.0, centre=src)
        _S0_CACHE[key] = line_integrals(g, frame, e0, xs, ys, R, envelope_fn(A0 * mult, 40.0, 1.0, centre=src, atten=atten))
    S0 = _S0_CACHE[key]
    Qc = (np.sqrt(phase_energy(g, S0, osc)) + np.sqrt(Qc_micro)) ** 2
    Q = diff['Q'] if quiet_Q is None else (np.sqrt(diff['Q']) + np.sqrt(quiet_Q)) ** 2
    return Q, Qc, diff, osc


def main():
    g = geometry()
    frame = image_frame(g)
    k0 = 4 * np.pi / g.lam
    p30 = load('p2_30_finite_certificates')
    micro_env = p30['background_included']['los_peak_envelope_m']
    micro_phase = k0 * micro_env
    # every realisation of the 24-wave field: its peak at most sqrt(24) times the envelope of one equal-power wave set
    micro_worst = np.sqrt(AN['directions_quiet']) * micro_phase
    Qc_micro = frame['n_band_along'] * micro_worst ** 2
    ql = quiet_level()
    params = {'frozen': FROZEN, 'frozen_hash': FROZEN_HASH, 'image_frame': {k: v for k, v in frame.items() if k != 'rho'},
              'pass_s': pass_seconds(g), 'microseism_los_peak_envelope_m': micro_env,
              'microseism_envelope_source': 'P2-30 background_included (three regional harmonics), times sqrt(24) for '
              'every realisation of the 24 directions', 'quiet_level': ql, 'snr': _SNR}
    with Run(RID, 'One auditable comparison: presence, location and shape, bounds beside achieved', params) as run:
        F = strong_fields(g)
        fs = F['f']
        print(f"  fields: record decays to {max(F['decay'].values()):.1e} of its energy in the last 10%", flush=True)
        # the frequency rule: the presence (room) pair's oracle on the measured disc, largest over the band
        disc = [disc_integral(np.abs(F['D']['WA'][i] - F['D']['W0'][i]) ** 2, F['xs'], F['ys'], AN['recorded_disc_m'])
                for i in range(len(fs))]
        fi = int(np.argmax(disc))
        f_star = float(fs[fi])
        print(f"  frozen frequency: {f_star:.0f} Hz", flush=True)
        # the quiet case
        kern = {w: memo(RID, f'kern_{w}', lambda w=w: static_kernels(w), __file__, version='kern-v1',
                        inputs={'boxes': FROZEN['worlds'][w], 'rock': FROZEN['rock']})
                for w in FROZEN['worlds']}
        eps0 = ql['strain_amplitude']
        osc_q = 1 / (2 * np.pi * ql['f_hz'] * pass_seconds(g))
        quiet, strong, curves = {}, {}, {}
        for pname, (w0, w1) in FROZEN['pairs'].items():
            t0 = time.time()
            dK = kernel_difference(kern[w0], kern[w1])
            xq, yq, eq, maps = quiet_energy(dK, g, eps0)
            integ_q = lines_for(g, frame, eq, xq, yq, AN['quiet_disc_m'], AN['quiet_ring_m'],
                                AN['quiet_energy_exponent'], 0.0)
            Qq = integ_q['Q']                   # expected over the field's realisations: the pass-mean of cos^2 is 1/2
            # every realisation: |dd| <= sum over directions of eps0 / sqrt(n) |M_phi|, and over a finite pass the mean of
            # cos^2(w t + phi) is at most (1 + 1 / (w T)) / 2 whatever the phase (osc_q; the exact |sin w T| / (w T) is smaller)
            e_max = (eps0 / np.sqrt(len(maps)) * np.sum(np.abs(maps), axis=0)) ** 2
            Qq_max = lines_for(g, frame, e_max, xq, yq, AN['quiet_disc_m'], AN['quiet_ring_m'],
                               AN['quiet_energy_exponent'], osc_q)['Q']
            Qc_q = np.full_like(Qq, Qc_micro)
            if kern[w0] is not None:           # world 0's own imprint in its reference too, at its worst realisation
                _, _, _, maps0 = quiet_energy(kern[w0], g, eps0)
                e0max = (eps0 / np.sqrt(len(maps0)) * np.sum(np.abs(maps0), axis=0)) ** 2
                i0q = lines_for(g, frame, e0max, xq, yq, AN['quiet_disc_m'], AN['quiet_ring_m'],
                                AN['quiet_energy_exponent'], osc_q)
                Qc_q = (np.sqrt(Qc_q) + np.sqrt(i0q['Q'])) ** 2
            lay = layers(g, frame, Qq, Qc_q, integ_q, 0.0)
            lay['l1_ensemble'] = l1_ensemble(Qq, Qq_max, Qc_q, 10 ** (-FROZEN['snr_db_headline'] / 10))
            ach = memo(RID, f'achieved_quiet_{pname}', lambda: _achieved_quiet(g, xq, yq, maps, eps0, ql['f_hz']),
                       __file__, version='achq-v1', inputs={'maps': maps, 'eps0': eps0, 'f': ql['f_hz'],
                                                             'disc': AN['quiet_disc_m'], 'acq': FROZEN['acquisition']})
            quiet[pname] = {'pair': [w0, w1], 'integrals': {k: v for k, v in integ_q.items() if k not in ('Q', 'Q_disc')},
                            **lay, 'achieved_score_test': ach,
                            'sigma_profile': sigma_profile(eq, xq, yq, step=4.0, r_max=AN['quiet_disc_m'])}
            # the strong case at the frozen frequency, and every frequency of the band for the curve
            Q, Qc, integ, osc = strong_pair(g, frame, F, fi, w0, w1, Qc_micro, quiet_Q=Qq_max)
            lay = layers(g, frame, Q, Qc, integ, osc)
            D0, D1 = F['D'][w0][fi], F['D'][w1][fi]
            sens = {}
            for label, kw in [('grid allowance (x1.11 amplitude, empirical)', {'grid': 1 + AN['grid_allowance_amplitude']})] + \
                    [(f'tail envelope x{m:g}', {'mult': m}) for m in AN['tail_multipliers'] if m != 1.0] + \
                    [(f'attenuation Q = {qa:g} (illustrative)', {'atten': 2 * np.pi * f_star / (qa * 1690.8)})
                     for qa in AN['attenuation_q']]:
                Qs, Qcs, ints, _ = strong_pair(g, frame, F, fi, w0, w1, Qc_micro, quiet_Q=Qq_max, **kw)
                c = l1_certificate(Qs, Qcs, 10 ** (-FROZEN['snr_db_headline'] / 10))
                o = oracle(g, ints['area_integral_m4'], FROZEN['snr_db_headline'], osc)
                sens[label] = {'l1_tv_upper': c['tv_upper'], 'l1_state': state(c['tv_upper']),
                               'oracle_tv_upper': o['tv_upper'], 'oracle_state': state(o['tv_upper'])}
            fields_in = {'worlds': WORLDS_HASH, 'pair': [w0, w1], 'f': f_star}
            ach = memo(RID, f'achieved_strong_{pname}', lambda: achieved_score_test(g, D1 - D0, F['xs'], F['ys'], f_star,
                                                                                  AN['recorded_disc_m']),
                       __file__, version='achs-v1', inputs={**fields_in, 'disc': AN['recorded_disc_m']})
            s2_head = 10 ** (-FROZEN['snr_db_headline'] / 10)
            ver = memo(RID, f'verify_{pname}', lambda: verify_line(g, D0, D1, F['xs'], F['ys'], f_star, s2_head, micro_phase),
                       __file__, version='ver-v1', inputs={**fields_in, 'sigma2': s2_head, 'micro_phase': micro_phase,
                                                            'px': AN['verify_line_px'], 'amps': AN['verify_amplifications']})
            fx = memo(RID, f'fixed_{pname}', lambda: fixed_scenes(g, D0, D1, F['xs'], F['ys'], f_star), __file__,
                      version='fix-v1', inputs={**fields_in, 'prf': AN['check_prf_hz'], 'n': AN['check_realisations'],
                                                'patch': AN['check_patch_m']})
            strong[pname] = {'pair': [w0, w1], 'f_hz': f_star, 'integrals': {k: v for k, v in integ.items() if k not in ('Q', 'Q_disc')},
                             **lay, 'sensitivity': sens, 'achieved_score_test': ach, 'verify_line': ver,
                             'fixed_scenes': fx,
                             'sigma_profile': sigma_profile(np.abs(D1 - D0) ** 2, F['xs'], F['ys'])}
            curves[pname] = memo(RID, f'curve_{pname}', lambda: _curve(g, frame, F, w0, w1, Qc_micro, Qq_max), __file__,
                                 version='curve-v2', inputs={**fields_in, 'Qc_micro': Qc_micro, 'Qq_max': Qq_max,
                                                             'snr': [FROZEN['snr_db_headline'], FROZEN['snr_db_bright']],
                                                             'tail': [AN['tail_ring_m'], AN['recorded_disc_m']]})
            nk = snr_key(FROZEN['snr_db_headline'])
            print(f"  {pname}: quiet L1 {quiet[pname]['l1_ensemble']['tv_upper']:.2e} (averaged), oracle "
                  f"{quiet[pname]['oracle'][nk]['tv_upper']:.2e}; strong at {f_star:.0f} Hz L1 "
                  f"{strong[pname]['l1'][nk]['tv_upper']:.3f} (within 70 m {strong[pname]['l1_disc']['tv_upper']:.3f}), "
                  f"oracle {strong[pname]['oracle'][nk]['tv_upper']:.3f} ({time.time() - t0:.0f} s)", flush=True)
        w0, w1 = FROZEN['pairs']['presence (room)']
        s2_head = 10 ** (-FROZEN['snr_db_headline'] / 10)
        mc = memo(RID, 'score_mc', lambda: score_mc(g, F['D'][w0][fi], F['D'][w1][fi], F['xs'], F['ys'], f_star,
                                                    s2_head, amps=AN['mc_amplifications'], seeds=AN['mc_seeds']),
                  __file__, version='mc-v4', inputs={'worlds': WORLDS_HASH, 'f': f_star, 'sigma2': s2_head,
                                                     'amps': AN['mc_amplifications'], 'seeds': AN['mc_seeds']})
        strong['presence (room)']['score_monte_carlo'] = mc
        run.save({'snr_db_headline': FROZEN['snr_db_headline'], 'snr_db_bright': FROZEN['snr_db_bright'],
                  'snr_key_headline': snr_key(FROZEN['snr_db_headline']), 'snr_key_bright': snr_key(FROZEN['snr_db_bright']),
                  'f_star_hz': f_star, 'frequency_rule': 'the presence (room) pair\'s oracle on the measured 70 m disc, '
                  'largest over the FTA band', 'disc_integral_by_f': {'f_hz': fs.tolist(), 'presence_room_m4': disc},
                  'field_decay_last_10pct': F['decay'], 'scale_per_unit_force': {'f_hz': fs.tolist(),
                                                                                  'scale': F['scale'].tolist()},
                  'microseism_common_phase_energy_per_line': Qc_micro, 'quiet': quiet, 'strong': strong,
                  'curves': curves, 'finding': finding(quiet, strong, curves, f_star)})
        res = {'f_star_hz': f_star, 'frequency_rule': 'the presence (room) pair\'s oracle on the measured disc, largest',
               'quiet': quiet, 'strong': strong, 'curves': curves}
        (RESULTS / RID / 'report.md').write_text(report(res))
        print(finding(quiet, strong, curves, f_star))


def _achieved_quiet(g, xq, yq, maps, eps0, f_hz):
    A = []
    for M in maps:
        a = achieved_score_test(g, eps0 * M.astype(complex), xq, yq, f_hz, AN['quiet_disc_m'], taper_m=8.0)
        A.append(a['fisher_mean_over_phase'])
    Fm = float(np.mean(A))
    return {'fisher_expected': Fm, 'deflection': float(np.sqrt(Fm)), 'tv_achieved': float(2 * ndtr(np.sqrt(Fm) / 2) - 1),
            'directions': len(maps), 'r_max_m': AN['quiet_disc_m']}


def _curve(g, frame, F, w0, w1, Qc_micro, Qq):
    rows = []
    for i in range(0, len(F['f']), 2):
        Q, Qc, integ, osc = strong_pair(g, frame, F, i, w0, w1, Qc_micro, quiet_Q=Qq)
        c = l1_certificate(Q, Qc, 10 ** (-FROZEN['snr_db_headline'] / 10))
        o = oracle(g, integ['area_integral_m4'], FROZEN['snr_db_headline'], osc)
        od = oracle(g, integ['disc_integral_m4'], FROZEN['snr_db_headline'], osc)
        og = oracle(g, integ['area_integral_m4'], FROZEN['snr_db_bright'], osc)
        rows.append({'f_hz': float(F['f'][i]), 'l1_tv_upper': c['tv_upper'], 'oracle_tv_upper': o['tv_upper'],
                     'oracle_bright_tv_upper': og['tv_upper'],
                     'oracle_disc_tv_upper': od['tv_upper'], 'tail_share': 1 - integ['disc_integral_m4'] /
                     integ['area_integral_m4']})
    return rows


def _fmt(x):
    if x is None:
        return 'none'
    return f'{x:.2g}' if (abs(x) < 1e-3 or abs(x) >= 1e4) else (f'{x:.3f}' if abs(x) < 1 else f'{x:.3g}')


def finding(quiet, strong, curves, f_star):
    n = snr_key(FROZEN['snr_db_headline'])
    gen = FROZEN['snr_db_bright']
    names = list(strong)
    s = strong
    pct = lambda v: f"{100 * v:.1f}%"
    q_ens = max(r['l1_ensemble']['tv_upper'] for r in quiet.values())
    q_worst = max(r['l1_ensemble']['worst_realisation_tv_upper'] for r in quiet.values())
    q_or = max(r['oracle'][n]['tv_upper'] for r in quiet.values())
    tail = min(s[k]['tail_share_of_integral'] for k in names), max(s[k]['tail_share_of_integral'] for k in names)
    l1 = {k: s[k]['l1'][n]['tv_upper'] for k in names}
    l1_disc = {k: s[k]['l1_disc']['tv_upper'] for k in names}
    band = {k: max(r['l1_tv_upper'] for r in curves[k]) for k in names}
    q50 = {k: v['l1_tv_upper'] for k in names for lab, v in s[k]['sensitivity'].items() if 'Q = 50' in lab}
    x5 = {k: v['l1_tv_upper'] for k in names for lab, v in s[k]['sensitivity'].items() if 'x5' in lab}
    orn = {k: s[k]['oracle'][n]['tv_upper'] for k in names}
    org = {k: s[k]['oracle'][snr_key(gen)]['tv_upper'] for k in names}
    sweep = FROZEN['snr_db_sweep']
    or_old = {k: s[k]['oracle'][snr_key(12.66)]['tv_upper'] for k in names} if snr_key(12.66) in s[names[0]]['oracle'] else None
    or40 = {k: s[k]['oracle'][snr_key(max(sweep))]['tv_upper'] for k in names}
    stop = {k: s[k]['oracle_bound_stops_excluding_at_snr_db'] for k in names}
    gr = min(s[k]['l1_growth_to_target'] or np.inf for k in names)
    rng_ = lambda d: f"{min(d.values()):.3f} to {max(d.values()):.3f}"
    fx = s[names[0]]['fixed_scenes']['summary']
    ver = [x['exact_over_certificate'] for k in names for x in s[k]['verify_line']['rows'][:1]]
    mcr = s[names[0]].get('score_monte_carlo', {})
    mc = mcr.get('rows')
    wc = s.get('presence (10 m room)')
    single_ok = max(band.values()) < AN['target_tv'] and max(x5.values()) < AN['target_tv'] and q_worst < AN['target_tv']
    ci = lambda x: f"{100 * x[0]:.0f} to {100 * x[1]:.0f}%"
    mc_txt = ''
    if mc:
        hi, lo = mc[-1], mc[0]
        mc_txt = (f" Achieved: the score test, implemented on synthesised images with one noise law, its 5% threshold set "
                  f"on {mcr['calibration_grounds']} calibration grounds and read on {mcr['evaluation_grounds_each']} + "
                  f"{mcr['cavity_evaluation_grounds']} evaluation grounds disjoint from them: with the difference amplified "
                  f"{hi['amplification']:g} times it finds {pct(hi['found_at_threshold_achieved'])} (95% interval "
                  f"{ci(hi['found_interval_95'])}) at {pct(hi['false_alarms_at_threshold_achieved'])} false alarms "
                  f"({ci(hi['false_alarm_interval_95'])}), AUC {hi['auc_achieved']:.2f} +- {hi['auc_se']:.2f} (predicted "
                  f"{hi['auc_predicted']:.2f}); at the real level {pct(lo['found_at_threshold_achieved'])} found at "
                  f"{pct(lo['false_alarms_at_threshold_achieved'])} false alarms (one-sided p {lo['found_over_false_alarms_p_one_sided']:.2f} "
                  f"that it finds more than it falsely alarms), AUC {lo['auc_achieved']:.2f} +- {lo['auc_se']:.2f} "
                  f"(p {lo['auc_p_one_sided']:.2f}): at chance.")
    return (
        f"One auditable comparison (BENCHMARK.md, revision 3): five worlds identical but for the cavity (none; a room 6 m on "
        f"a side under a 5 m roof; an L-shaped tunnel at that depth; the room 6 m east; and the 10 m room under a 5 m roof "
        f"that earlier runs left open), the same lorry 15 m west at the FTA's truck-over-a-bump level held at {f_star:.0f} Hz "
        f"all pass, the same scatterers and noise law, one geology (uniform limestone, no attenuation in the solver), fully "
        f"developed speckle. The whole image is observed: the solver's field within 70 m, and beyond it an assumed envelope, "
        f"a surface wave carried undiminished to the image's edge, which supplies {100 * tail[0]:.0f} to {100 * tail[1]:.0f}% "
        f"of the strong case's signal. The signal-to-noise ratio is conditional: the site's backscatter measured in the "
        f"image over ICEYE's specified noise floor, the products carrying none; the headline takes the best specified for "
        f"the product's mode (Dwell Fine, documentation 6.0.8), {FROZEN['snr_db_headline']:.1f} dB, the most favourable to "
        f"a detector, and every other specified value is in the sweep. For any method reading the one image, told the "
        f"excitation and both worlds' whole motion but not the speckle (the certificate line by line, each line's reference "
        f"holding the lorry's wave, the microseisms' worst-case envelope and the noise), detection may exceed false alarms "
        f"by at most {rng_(l1_disc)} within 70 m; with the envelope, {rng_(l1)}"
        + (f" (the 10 m room {wc['l1'][n]['tv_upper']:.3f})" if wc else '')
        + f"; with the scattered wave attenuated at Q = 50, a realistic rock, {rng_(q50)}; with the envelope five times "
        f"larger, {rng_(x5)}; the largest over the sampled band (every 2 Hz) {rng_(band)}. Found at 5% false alarms is at most "
        f"{pct(0.05 + max(band.values()))} for any of these questions; the lorry's differential motion would have to grow "
        f"{gr:.0f} times to reach 95% at 5%. The oracle told every scatterer's reflectivity as well: at "
        f"{FROZEN['snr_db_headline']:.1f} dB {rng_(orn)}"
        + (f" ({rng_(or_old)} at 12.7 dB, documentation 6.0.0's best)" if or_old else '')
        + f", at {gen:.1f} dB (ground at 0 dB) {rng_(org)}; its upper bound stops excluding 95% at 5% only at "
        f"{min(stop.values()):.0f} to {max(stop.values()):.0f} dB, which no specified ICEYE noise floor reaches over natural "
        f"ground. The one open corner of the strong case is bright persistent ground: at {max(sweep):.0f} dB, which stands "
        f"for a building, a corner reflector or the pyramid's edge, the oracle's bound is {rng_(or40)} for every strong pair "
        f"and the speckle layer does not apply, so a shallow room beside a lorry under a bright persistent target is not "
        f"excluded by any layer here; that is the real-acquisition case (diffuse ground with bright scatterers) still to "
        f"compute. Under Giza's regional microseisms, averaged over the field's realisations through a bound linear in each "
        f"line's phase energy, one image is at most {q_ens:.1e} above chance ({q_worst:.1e} at the worst realisation, its "
        f"pass-mean bounded over the finite window) and the oracle {q_or:.1e}. Checks: the exact divergence on lines of the "
        f"model is {min(ver):.2f} to {max(ver):.2f} of the certificate; pulse by pulse on fixed scenes the coherent echo "
        f"difference is {fx['common_x1']['mean']:.3f} +- {fx['common_x1']['se']:.3f} of the ensemble formula, the same "
        f"with the lorry's wave ten times larger."
        + mc_txt
        + (" So in this benchmark one image does not reach 95% at 5% for presence, a 6 m location or shape, quiet or beside a "
           "lorry, under every tabulated allowance, and at every specified noise floor over natural ground neither does the "
           "oracle's bound. " if single_ok else " Some allowances leave cases open; see the report. ")
        + "The scope is the benchmark's: five layouts, one geology, the declared excitation, fully developed speckle "
        "(bright points are covered by the oracle only), the assumed envelope beyond 70 m, a specified rather than measured "
        "noise floor; a shallow room is not by itself an upper bound on every larger or deeper structure.")


def report(res):
    """The benchmark's tables, generated from the summary."""
    n = snr_key(FROZEN['snr_db_headline'])
    L = [f"# P2-36 benchmark (revision 3): generated report", '',
         f"Frozen parameters hash `{FROZEN_HASH}` (worlds `{WORLDS_HASH}`, unchanged since revision 2); frozen frequency "
         f"{res['f_star_hz']:.0f} Hz ({res['frequency_rule']}). SNR conditional: the site's median sigma0 measured in the "
         f"image over ICEYE's specified noise floor; headline {FROZEN['snr_db_headline']:.2f} dB (Dwell Fine's best, "
         f"documentation 6.0.8), {FROZEN['snr_db_bright']:.1f} dB for ground at 0 dB, every other specified value and the "
         f"measured lower bound in the sweep. States: near chance (TV < 0.05), 95% at 5% excluded (TV < 0.9), "
         f"unresolved. Each entry bounds detection rate minus false-alarm rate; 'predicted' is a weak-signal calculation, "
         f"'achieved' a simulated detector.", '']
    for exc in ('quiet', 'strong'):
        L += [f"## {exc}", '', '| pair | L1 whole image | state | L1 within 70 m | growth to 0.9 | oracle | state | '
              'oracle bound stops excluding 95/5 at | oracle within 70 m | predicted (score test, within 70 m) | tail share |'
              + (' L1 averaged over realisations (linear bound) | worst realisation |' if exc == 'quiet' else ''),
              '|---|---|---|---|---|---|---|---|---|---|---|' + ('---|---|' if exc == 'quiet' else '')]
        for k, r in res[exc].items():
            ach = r['achieved_score_test']
            a_ = ach.get('tv_achieved_best_phase', ach.get('tv_achieved'))
            row = (f"| {k} | {_fmt(r['l1'][n]['tv_upper'])} | {r['l1'][n]['state']} | {_fmt(r['l1_disc']['tv_upper'])} | "
                   f"{_fmt(r['l1_growth_to_target'])} | {_fmt(r['oracle'][n]['tv_upper'])} | {r['oracle'][n]['state']} | "
                   f"{_fmt(r['oracle_bound_stops_excluding_at_snr_db'])} dB | {_fmt(r['oracle_disc']['tv_upper'])} | "
                   f"{_fmt(a_)} | {100 * r['tail_share_of_integral']:.0f}% |")
            if exc == 'quiet':
                e = r['l1_ensemble']
                row += f" {_fmt(e['tv_upper'])} | {_fmt(e['worst_realisation_tv_upper'])} |"
            L.append(row)
        L += ['', "SNR sweep (L1 / oracle):", '', '| pair | ' + ' | '.join(f'{d:g} dB' for d in FROZEN['snr_db_sweep'])
              + ' |', '|---|' + '---|' * len(FROZEN['snr_db_sweep'])]
        for k, r in res[exc].items():
            L.append(f"| {k} | " + ' | '.join(f"{_fmt(r['l1'][snr_key(d)]['tv_upper'])} / {_fmt(r['oracle'][snr_key(d)]['tv_upper'])}"
                                             for d in FROZEN['snr_db_sweep']) + ' |')
        L.append('')
    L += ['## strong: allowances and sensitivity (the envelope is a declared assumption; these are robustness tests, '
          'not its verification)', '', '| pair | case | L1 | oracle |', '|---|---|---|---|']
    for k, r in res['strong'].items():
        for lab, v in r['sensitivity'].items():
            L.append(f"| {k} | {lab} | {_fmt(v['l1_tv_upper'])} ({v['l1_state']}) | {_fmt(v['oracle_tv_upper'])} "
                     f"({v['oracle_state']}) |")
    L += ['', '## strong: checks', '']
    for k, r in res['strong'].items():
        fx, ver, sp = r['fixed_scenes']['summary'], r['verify_line'], r['sigma_profile']
        L.append(f"- {k}: fixed scenes exact/formula {fx['common_x1']['mean']:.3f} ± {fx['common_x1']['se']:.3f} "
                 f"(common wave ×10: {fx['common_x10']['mean']:.3f} ± {fx['common_x10']['se']:.3f}); one line, exact KL / "
                 f"certificate " + ', '.join(f"×{x['amplification']:g}: {_fmt(x['exact_over_certificate'])}" for x in ver['rows'])
                 + f"; line floor {ver['floor']:.4f}; σ(r) slope over 40–70 m {sp['slope_40_70']:.2f}, σ(67.5)/σ(42.5) "
                 f"{sp['outer_over_40m']:.2f}")
    mc = res['strong'].get('presence (room)', {}).get('score_monte_carlo')
    if mc:
        L += ['', f"## strong: the score test achieved on synthesised images (presence, room; a cyclic patch "
              f"{mc['patch_m'][0]:.0f} × {mc['patch_m'][1]:.0f} m; threshold from {mc['calibration_grounds']} calibration "
              f"grounds, evaluation on {mc['evaluation_grounds_each']} + {mc['cavity_evaluation_grounds']} grounds disjoint "
              f"from them, seed streams {mc['seeds']})", '',
              '| amplification | AUC achieved ± se (p) | AUC predicted | found, achieved [95%] | false alarms, achieved [95%] | '
              'found > false alarms, p | found at 5%, predicted | patch certificate |', '|---|---|---|---|---|---|---|---|']
        for x in mc['rows']:
            iv = lambda v: f"[{v[0]:.2f}, {v[1]:.2f}]"
            L.append(f"| ×{x['amplification']:g} | {x['auc_achieved']:.3f} ± {x['auc_se']:.3f} ({x['auc_p_one_sided']:.2g}) | "
                     f"{x['auc_predicted']:.3f} | {x['found_count'][0]}/{x['found_count'][1]} {iv(x['found_interval_95'])} | "
                     f"{x['false_alarm_count'][0]}/{x['false_alarm_count'][1]} {iv(x['false_alarm_interval_95'])} | "
                     f"{x['found_over_false_alarms_p_one_sided']:.2g} | {x['found_at_005_predicted']:.2f} | "
                     f"{_fmt(x['certificate_tv_upper_patch'])} |")
    L += ['', '## strong: over the sampled FTA band (every 2 Hz; not a proved maximum between samples)', '',
          '| f (Hz) | ' + ' | '.join(res['curves']) + ' |', '|---|' + '---|' * len(res['curves'])]
    names = list(res['curves'])
    for i, row in enumerate(res['curves'][names[0]]):
        L.append(f"| {row['f_hz']:.0f} | " + ' | '.join(
            f"{_fmt(res['curves'][k][i]['l1_tv_upper'])} / {_fmt(res['curves'][k][i]['oracle_tv_upper'])}" for k in names) + ' |')
    L += ['', 'Each cell of the band table: L1 / oracle at the headline SNR.', '']
    return '\n'.join(L)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'solve':
        for w in FROZEN['worlds']:
            solve(w)
    if len(sys.argv) > 1 and sys.argv[1] == 'static':
        for w in FROZEN['worlds']:
            static_kernels(w)
    if len(sys.argv) == 1:
        main()
