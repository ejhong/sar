"""P2-36 · One auditable comparison: what one image could tell of presence, location and shape, and what readers achieve.

    uv run python experiments/p2_36_benchmark.py

A frozen benchmark (BENCHMARK.md), after two independent reviews asked for one fully specified comparison before any
parameter ranges: the same source in every world, a valid reference for each world, the observation and its tails
declared, the scattering and noise model stated, bounds for presence, horizontal location and shape, and implemented
detectors' achieved performance beside the bounds, with what each is told recorded.

Worlds (identical ground, scatterers and source; only the cavity differs):
  W0  no cavity;  WA  a room 6 m on a side, roof 5 m down;  WB  an L-shaped tunnel 2.5 m square, roof 5 m down;
  WA6 the same room 6 m east.
  Pairs: presence (W0-WA, W0-WB), horizontal location (WA-WA6), shape (WA-WB); absolute depth is not scored.
Excitations:
  quiet   Giza's regional microseisms (0.1-0.3 Hz, measured 67 km east), the imprint by the solver settling under a
          uniform strain (P2-28), waves from every direction; the oracle ensemble-averaged over the field (its expected
          squared motion, no peak allowance).
  strong  one vertical point force 15 m west of the site (a lorry), the same force in every world, solved by the
          elastic solver on a box recorded to 70 m; its amplitude scaled so the vertical velocity at the site without a
          cavity is the FTA's truck-over-a-bump level at 15 m; one declared frequency (the rule: the frequency in the
          FTA band where the presence pair's ensemble oracle is largest, found and then frozen).
Observation: the image's scatterers within a disk of radius R_OBS about the site (declared); the scatterers outside
  add a tail, bounded separately (information outside is at most added: the two sets are independent in the model).
Scattering and noise: fully developed speckle (white circular Gaussian scatterers on the image's pixel grid, texture
  power one per cell) and white receiver noise at SNR per cell (30 dB nominal, assumed: sigma0 / NESZ for the brightest
  natural ground over ICEYE Dwell's best specified NESZ, -26.7 dB; the acquisition's own calibration and noise records
  are not in this repository), swept 20 to 60 dB. The Doppler-to-time relation (checked pulse by pulse, P2-25).
Layers:
  L1  unknown reflectivity (told the excitation and both worlds' motion, not the texture): the exact KL per
      along-track line, each world's own covariance (its common motion and receiver noise in it), summed over lines
      (independent with the range band widened, which can only add information); no floor, no Fisher approximation.
  L2  the oracle told the reflectivity too: the complete coherent echo difference for fixed scenes (the synthesizer in
      complex128), Delta^2 = 2 beta |dz|^2 / sigma_z^2, and its ensemble average checked against
      E Delta^2 = 8 sum SNR <sin^2(k du / 2)>, which holds under independent uniform scatterer phases.
Implemented detectors:
  D1  the oracle's likelihood-ratio test on noisy images (Monte Carlo), against 2 Phi(Delta / 2) - 1.
  D2  the Neyman-Pearson test of L1 (known excitation, unknown texture): its exact ROC from the generalised eigenvalues.
  D3  the published method (P2-07's pipeline) read blind, as an ordinary mapper, on noisy images of the strong case.
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

RID = 'p2_36_benchmark'
SITES = Path(__file__).resolve().parents[2] / 'sites'
FROZEN = {
    'worlds': {
        'W0': [],
        'WA': [[[0.0, 0.0, -8.0], [6.0, 6.0, 6.0]]],
        'WB': [[[-5.0, 0.0, -6.25], [18.0, 2.5, 2.5]], [[2.75, 7.375, -6.25], [2.5, 17.25, 2.5]]],
        'WA6': [[[6.0, 0.0, -8.0], [6.0, 6.0, 6.0]]],
    },
    'pairs': {'presence (room)': ['W0', 'WA'], 'presence (L tunnel)': ['W0', 'WB'], 'location (6 m)': ['WA', 'WA6'],
              'shape (room or L tunnel)': ['WA', 'WB']},
    'strong_source': {'kind': 'vertical point force', 'east_m': -15.0, 'north_m': 0.0, 'level': 'FTA truck over a bump, 15 m'},
    'solver': {'extent': [[-100.0, 100.0], [-100.0, 100.0], [-100.0, 3.0]], 'h': 1.0, 'pml_m': 20.0, 'half_m': 70.0,
               'record_s': 0.5, 'f_top': 120.0, 'record_every': 4},
    'r_obs_m': 39.0,
    'snr_db_nominal': 30.0,
    'snr_db_sweep': [20.0, 30.0, 40.0, 50.0, 60.0],
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
    },
}
FROZEN_HASH = hashlib.sha1(json.dumps(FROZEN, sort_keys=True).encode()).hexdigest()[:12]
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
    path = CACHE / f'fields_{FROZEN_HASH}.npz'
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


def score_mc(g, D0, D1, xs, ys, f_hz, sigma2, amps=(1.0, 800.0), seeds=40, shape=(1024, 64)):
    """D2 implemented: the score statistic for the pair's difference, told its pattern and timing but not the speckle,
    applied to synthesised images (sarsim.synthesize, complex128) of a cyclic patch about the largest difference:
    scatterers on the pixel grid (white circular Gaussian), world 0 moving with its whole motion (the lorry's wave and
    its own cavity), world 1 the same plus the difference amplified `amp` times, the same scatterers in both, white
    receiver noise. Its empirical AUC and detection at 5% false alarms against the weak-signal prediction from the exact
    Fisher information on the same patch (fisher_grid), and against the per-line certificate on the same patch."""
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
    ph_x = np.exp(1j * 2 * np.pi * f_hz * XA / g.V)
    rows = []
    for amp in amps:
        dK = amp * (z1p - z0p)
        A_f, B_f = inf.fisher_grid(dK, g, f_hz, pad=(Nx, Nr))
        taus = [lambda t: np.cos(2 * np.pi * f_hz * t), lambda t: np.sin(2 * np.pi * f_hz * t)]
        T = {0: [], 1: []}
        for sd in range(seeds + 1):
            rng = np.random.default_rng(3600 + sd)
            sc = (rng.normal(size=Nx * Nr) + 1j * rng.normal(size=Nx * Nr)) / np.sqrt(2)
            n = Nx * Nr
            zero = np.zeros(n)
            scat = Scatterers(x=XA.ravel(), y=YR.ravel(), z=zero, amp=np.abs(sc), phase=np.angle(sc), iso=np.ones(n),
                              flash=zero, nu0=zero, sig_nu=np.ones(n), vib_amp=zero, vib_freq=zero, vib_phase=zero,
                              label=np.zeros(n, int))
            psi = rng.uniform(0, 2 * np.pi)
            noise = (rng.normal(size=(Nx, Nr)) + 1j * rng.normal(size=(Nx, Nr))) / np.sqrt(2)
            for w in (0, 1):
                def motion(x, y, z, t, w=w):
                    E, N = site(x, y)
                    zz = interp_complex(D0, xs, ys, E, N) + w * amp * (interp_complex(D1, xs, ys, E, N) -
                                                                      interp_complex(D0, xs, ys, E, N))
                    return np.real(zz[None, :] * np.exp(1j * (2 * np.pi * f_hz * (t[:, None] + x[None, :] / g.V) + psi)))
                img = synthesize(scat, g, shape, motion=motion, dtype=np.complex128)
                img = img + np.sqrt(sigma2 * np.mean(np.abs(img) ** 2)) * noise
                ps = np.exp(1j * psi)
                q = m31.score_fields(img, g, taus)
                Kp1, Kp2 = np.real(dK * ph_x * ps), -np.imag(dK * ph_x * ps)     # d = K1 cos(2 pi f t) + K2 sin
                T[w].append(float(k0 * np.sum(Kp1 * q[0] + Kp2 * q[1])))
        # the first seed fixes the statistic's sign (the model's phase convention), the rest are scored
        sign = 1.0 if T[1][0] >= T[0][0] else -1.0
        t0, t1 = sign * np.array(T[0][1:]), sign * np.array(T[1][1:])
        auc = float(np.mean(t1[:, None] > t0[None, :]))
        thr = np.quantile(t0, 0.95)
        tpr = float(np.mean(t1 > thr))
        d = float(np.sqrt(A_f))
        osc = 1 / (2 * np.pi * f_hz * pass_seconds(g))
        Q = g.band_frac * k0 ** 2 / 2 * (1 + osc) * np.sum(np.abs(dK) ** 2, axis=0)        # per range line
        Qc = g.band_frac * k0 ** 2 / 2 * (1 + osc) * np.sum(np.abs(z0p) ** 2, axis=0)
        cert = l1_certificate(Q, Qc, sigma2)
        rows.append({'amplification': amp, 'auc_empirical': auc, 'auc_predicted': float(ndtr(d / np.sqrt(2))),
                     'tpr_at_005_empirical': tpr, 'tpr_at_005_predicted': float(ndtr(ndtri(0.05) + d)),
                     'deflection_predicted': d, 'tv_predicted': float(2 * ndtr(d / 2) - 1),
                     'certificate_tv_upper_patch': cert['tv_upper'], 'sign': sign,
                     'mean_shift_over_sd': float((t1.mean() - t0.mean()) / t0.std(ddof=1))})
        print(f"    score MC x{amp:g}: AUC {auc:.3f} (predicted {rows[-1]['auc_predicted']:.3f}), shift/sd "
              f"{rows[-1]['mean_shift_over_sd']:.2f} (predicted d {d:.3f}), patch certificate {cert['tv_upper']:.3g}",
              flush=True)
    return {'patch_px': list(shape), 'patch_m': [Nx * g.dx, Nr * g.dr / np.sin(g.theta)], 'centre_site_m': c.tolist(),
            'seeds': seeds, 'rows': rows}


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
            'obs_disc_integral_m4': disc_integral(emap, xs, ys, FROZEN['r_obs_m'])}


def layers(g, frame, Q, Qc, integ, osc, extra=None):
    """L1 at every SNR of the sweep (the noise only raises the floor), the oracle at every SNR, their margins."""
    out = {'l1': {}, 'oracle': {}}
    for snr in FROZEN['snr_db_sweep']:
        s2 = 10 ** (-snr / 10)
        c = l1_certificate(Q, Qc, s2)
        c['state'] = state(c['tv_upper'])
        out['l1'][f'{snr:g}'] = c
        o = oracle(g, integ['area_integral_m4'], snr, osc)
        o['state'] = state(o['tv_upper'])
        out['oracle'][f'{snr:g}'] = o
    s2 = 10 ** (-FROZEN['snr_db_nominal'] / 10)
    out['l1_growth_to_target'] = growth(lambda a: l1_certificate(Q, Qc, s2, a)['tv_upper'], AN['target_tv'])
    o30 = out['oracle'][f"{FROZEN['snr_db_nominal']:g}"]
    out['oracle_growth_to_target'] = float(np.sqrt(oracle_d2_target() / o30['e_delta2'])) if o30['e_delta2'] > 0 else None
    out['oracle_snr_db_for_target'] = (FROZEN['snr_db_nominal'] + 10 * np.log10(oracle_d2_target() / o30['e_delta2'])
                                       if o30['e_delta2'] > 0 else None)
    out['l1_disc'] = l1_certificate(integ['Q_disc'], Qc, s2)
    out['oracle_disc'] = oracle(g, integ['disc_integral_m4'], FROZEN['snr_db_nominal'], osc)
    out['oracle_obs_disc'] = oracle(g, integ['obs_disc_integral_m4'], FROZEN['snr_db_nominal'], osc)
    out['tail_share_of_integral'] = 1 - integ['disc_integral_m4'] / integ['area_integral_m4']
    return out


def strong_pair(g, frame, F, fi, w0, w1, Qc_micro, quiet_Q=None, mult=1.0, atten=None, grid=1.0):
    """The strong case for one pair at one frequency: the differential map and the world-0 map, per line."""
    f = float(F['f'][fi])
    osc = 1 / (2 * np.pi * f * pass_seconds(g))
    xs, ys = F['xs'], F['ys']
    D0, D1 = F['D'][w0][fi], F['D'][w1][fi]
    R = AN['recorded_disc_m']
    e = np.abs(D1 - D0) ** 2 * grid ** 2
    diff = lines_for(g, frame, e, xs, ys, R, AN['tail_ring_m'], 1.0, osc, mult, atten)
    e0 = np.abs(D0) ** 2 * grid ** 2
    src = (FROZEN['strong_source']['east_m'], FROZEN['strong_source']['north_m'])
    A0 = ring_envelope(e0, xs, ys, (40.0, 55.0), 40.0, 1.0, centre=src)
    S0 = line_integrals(g, frame, e0, xs, ys, R, envelope_fn(A0 * mult, 40.0, 1.0, centre=src, atten=atten))
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
    Qc_micro = frame['n_band_along'] * micro_phase ** 2
    ql = quiet_level()
    params = {'frozen': FROZEN, 'frozen_hash': FROZEN_HASH, 'image_frame': {k: v for k, v in frame.items() if k != 'rho'},
              'pass_s': pass_seconds(g), 'microseism_los_peak_envelope_m': micro_env,
              'microseism_envelope_source': 'P2-30 background_included (three regional harmonics)',
              'quiet_level': ql, 'snr_status': 'assumed: sigma0 about 0 dB over ICEYE Dwell best NESZ -26.7 dB; the '
              'acquisition\'s calibration and noise records are not in this repository'}
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
        kern = {w: memo(RID, f'kern_{w}', lambda w=w: static_kernels(w), __file__, version='kern-v1')
                for w in FROZEN['worlds']}
        eps0 = ql['strain_amplitude']
        quiet, strong, curves = {}, {}, {}
        for pname, (w0, w1) in FROZEN['pairs'].items():
            t0 = time.time()
            dK = kernel_difference(kern[w0], kern[w1])
            xq, yq, eq, maps = quiet_energy(dK, g, eps0)
            integ_q = lines_for(g, frame, eq, xq, yq, AN['quiet_disc_m'], AN['quiet_ring_m'],
                                AN['quiet_energy_exponent'], 0.0)
            Qq = integ_q['Q']
            Qc_q = np.full_like(Qq, Qc_micro)
            if kern[w0] is not None:           # world 0's own imprint in its reference too
                _, _, e0q, _ = quiet_energy(kern[w0], g, eps0)
                i0q = lines_for(g, frame, e0q, xq, yq, AN['quiet_disc_m'], AN['quiet_ring_m'],
                                AN['quiet_energy_exponent'], 0.0)
                Qc_q = (np.sqrt(Qc_q) + np.sqrt(i0q['Q'])) ** 2
            lay = layers(g, frame, Qq, Qc_q, integ_q, 0.0)
            ach = memo(RID, f'achieved_quiet_{pname}', lambda: _achieved_quiet(g, xq, yq, maps, eps0, ql['f_hz']),
                       __file__, version='achq-v1')
            quiet[pname] = {'pair': [w0, w1], 'integrals': {k: v for k, v in integ_q.items() if k not in ('Q', 'Q_disc')},
                            **lay, 'achieved_score_test': ach,
                            'sigma_profile': sigma_profile(eq, xq, yq, step=4.0, r_max=AN['quiet_disc_m'])}
            # the strong case at the frozen frequency, and every frequency of the band for the curve
            Q, Qc, integ, osc = strong_pair(g, frame, F, fi, w0, w1, Qc_micro, quiet_Q=Qq)
            lay = layers(g, frame, Q, Qc, integ, osc)
            D0, D1 = F['D'][w0][fi], F['D'][w1][fi]
            sens = {}
            for label, kw in [('grid allowance (x1.11 amplitude, empirical)', {'grid': 1 + AN['grid_allowance_amplitude']})] + \
                    [(f'tail envelope x{m:g}', {'mult': m}) for m in AN['tail_multipliers'] if m != 1.0] + \
                    [(f'attenuation Q = {qa:g} (illustrative)', {'atten': 2 * np.pi * f_star / (qa * 1690.8)})
                     for qa in AN['attenuation_q']]:
                Qs, Qcs, ints, _ = strong_pair(g, frame, F, fi, w0, w1, Qc_micro, quiet_Q=Qq, **kw)
                c = l1_certificate(Qs, Qcs, 10 ** (-FROZEN['snr_db_nominal'] / 10))
                o = oracle(g, ints['area_integral_m4'], FROZEN['snr_db_nominal'], osc)
                sens[label] = {'l1_tv_upper': c['tv_upper'], 'l1_state': state(c['tv_upper']),
                               'oracle_tv_upper': o['tv_upper'], 'oracle_state': state(o['tv_upper'])}
            ach = memo(RID, f'achieved_strong_{pname}', lambda: achieved_score_test(g, D1 - D0, F['xs'], F['ys'], f_star,
                                                                                  AN['recorded_disc_m']),
                       __file__, version='achs-v1')
            ver = memo(RID, f'verify_{pname}', lambda: verify_line(g, D0, D1, F['xs'], F['ys'], f_star,
                                                                  10 ** (-FROZEN['snr_db_nominal'] / 10), micro_phase),
                       __file__, version='ver-v1')
            fx = memo(RID, f'fixed_{pname}', lambda: fixed_scenes(g, D0, D1, F['xs'], F['ys'], f_star), __file__,
                      version='fix-v1')
            strong[pname] = {'pair': [w0, w1], 'f_hz': f_star, 'integrals': {k: v for k, v in integ.items() if k not in ('Q', 'Q_disc')},
                             **lay, 'sensitivity': sens, 'achieved_score_test': ach, 'verify_line': ver,
                             'fixed_scenes': fx,
                             'sigma_profile': sigma_profile(np.abs(D1 - D0) ** 2, F['xs'], F['ys'])}
            curves[pname] = memo(RID, f'curve_{pname}', lambda: _curve(g, frame, F, w0, w1, Qc_micro, Qq), __file__,
                                 version='curve-v1')
            print(f"  {pname}: quiet L1 {quiet[pname]['l1']['30']['tv_upper']:.2e}, oracle "
                  f"{quiet[pname]['oracle']['30']['tv_upper']:.2e}; strong at {f_star:.0f} Hz L1 "
                  f"{strong[pname]['l1']['30']['tv_upper']:.3f}, oracle {strong[pname]['oracle']['30']['tv_upper']:.3f} "
                  f"({time.time() - t0:.0f} s)", flush=True)
        w0, w1 = FROZEN['pairs']['presence (room)']
        mc = memo(RID, 'score_mc', lambda: score_mc(g, F['D'][w0][fi], F['D'][w1][fi], F['xs'], F['ys'], f_star,
                                                    10 ** (-FROZEN['snr_db_nominal'] / 10), amps=(1.0, 800.0), seeds=24),
                  __file__, version='mc-v1')
        strong['presence (room)']['score_monte_carlo'] = mc
        run.save({'f_star_hz': f_star, 'frequency_rule': 'the presence (room) pair\'s oracle on the measured 70 m disc, '
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
        c = l1_certificate(Q, Qc, 10 ** (-FROZEN['snr_db_nominal'] / 10))
        o = oracle(g, integ['area_integral_m4'], FROZEN['snr_db_nominal'], osc)
        od = oracle(g, integ['disc_integral_m4'], FROZEN['snr_db_nominal'], osc)
        rows.append({'f_hz': float(F['f'][i]), 'l1_tv_upper': c['tv_upper'], 'oracle_tv_upper': o['tv_upper'],
                     'oracle_disc_tv_upper': od['tv_upper'], 'tail_share': 1 - integ['disc_integral_m4'] /
                     integ['area_integral_m4']})
    return rows


def _fmt(x):
    if x is None:
        return 'none'
    return f'{x:.2g}' if (abs(x) < 1e-3 or abs(x) >= 1e4) else (f'{x:.3f}' if abs(x) < 1 else f'{x:.3g}')


def finding(quiet, strong, curves, f_star):
    n = f"{FROZEN['snr_db_nominal']:g}"
    q_l1 = max(r['l1'][n]['tv_upper'] for r in quiet.values())
    q_or = max(r['oracle'][n]['tv_upper'] for r in quiet.values())
    q_gr = min(r['l1_growth_to_target'] or np.inf for r in quiet.values())
    q_go = min(r['oracle_growth_to_target'] or np.inf for r in quiet.values())
    s = strong
    names = list(s)
    l1 = {k: s[k]['l1'][n] for k in names}
    orc = {k: s[k]['oracle'][n] for k in names}
    worst_l1_band = {k: max(r['l1_tv_upper'] for r in curves[k]) for k in names}
    worst_or_band = {k: max(r['oracle_tv_upper'] for r in curves[k]) for k in names}
    gr = min(s[k]['l1_growth_to_target'] or np.inf for k in names)
    fx = s[names[0]]['fixed_scenes']['summary']
    ver = s[names[0]]['verify_line']['rows'][0]
    ach = {k: s[k]['achieved_score_test']['tv_achieved_best_phase'] for k in names}
    tail = max(s[k]['tail_share_of_integral'] for k in names)
    sens_worst = max(v['l1_tv_upper'] for k in names for v in s[k]['sensitivity'].values())
    snr_need = {k: s[k]['oracle_snr_db_for_target'] for k in names}
    listing = lambda d, f=_fmt: '; '.join(f"{k} {f(v)}" for k, v in d.items())
    single_ok = q_l1 < AN['target_tv'] and max(worst_l1_band.values()) < AN['target_tv'] and sens_worst < AN['target_tv']
    oracle_ok = max(orc[k]['tv_upper'] for k in names) < AN['target_tv']
    conclusion = (("So in this benchmark one image under fully developed speckle does not reach 95% at 5% for presence, a "
                   "6 m location or shape, quiet or favourable, whatever reads it; only a detector told the reflectivity "
                   + ("could, and only at an SNR above the assumed one." if oracle_ok else "reaches it at the assumed SNR."))
                  if single_ok else
                  "So in this benchmark the favourable case is not excluded for a single image under every allowance; "
                  "the pairs and allowances where it is not are listed in the report.")
    states = lambda d: '; '.join(f"{k} {v['tv_upper']:.3f} ({v['state']})" for k, v in d.items())
    return (
        f"One auditable comparison (BENCHMARK.md): four worlds identical but for the cavity (none; a room 6 m on a side, roof "
        f"5 m down; an L-shaped tunnel at the same depth; the room 6 m east), the same source, scatterers and noise law in "
        f"each, the whole image observed, the tail beyond the measured 70 m an assumed envelope. Quiet ground (Giza's "
        f"regional microseisms, expected over the field): every pair near chance on every layer, any single-image reader "
        f"at most {_fmt(q_l1)} above its false-alarm rate, the oracle told every scatterer's reflectivity at most {_fmt(q_or)} "
        f"at {n} dB (assumed); the differential motion would have to grow {_fmt(q_gr)} times (single image) and "
        f"{_fmt(q_go)} times (oracle) to reach 95% found at 5%. Favourable strong shaking (a truck over a bump 15 m away, "
        f"its FTA level held at {f_star:.0f} Hz for the whole pass, the scattered wave carried unattenuated to the image's "
        f"edge, {100 * tail:.0f}% of it beyond the measured disc): for any method reading the one image, told the "
        f"excitation and both worlds' motion but not the speckle, the per-line certificate, each line's reference holding "
        f"the lorry's wave, gives {states(l1)}; the worst over the band {listing(worst_l1_band)}; under the tabulated "
        f"allowances at most {_fmt(sens_worst)}; the motion would have to grow {_fmt(gr)} times. The oracle told every "
        f"reflectivity, a privileged detector (a perfect reference image), gives at {n} dB {states(orc)} (worst over the "
        f"band {listing(worst_or_band)}) and reaches 0.9 at {listing(snr_need, lambda v: f'{v:.0f} dB')}. Pulse by pulse on "
        f"fixed scenes the coherent echo difference is {fx['common_x1']['mean']:.3f} +- {fx['common_x1']['se']:.3f} of the "
        f"ensemble formula ({fx['common_x10']['mean']:.3f} +- {fx['common_x10']['se']:.3f} with the lorry's wave ten times "
        f"larger); on one line of the model the exact KL is {_fmt(ver['exact_over_certificate'])} of the certificate. The "
        f"score test told the pattern (achieved, within 70 m) reaches {listing(ach)}"
        + (lambda mc: (f"; implemented on synthesised images it reaches the predicted deflection (mean shift "
                       f"{mc[-1]['mean_shift_over_sd']:.2f} spreads against {mc[-1]['deflection_predicted']:.2f} predicted, "
                       f"with the difference amplified {mc[-1]['amplification']:g} times; AUC {mc[0]['auc_empirical']:.2f} at the "
                       f"real level)" if mc else ''))(s[names[0]].get('score_monte_carlo', {}).get('rows'))
        + ". " + conclusion
        + " It rests on the elastic model without attenuation, the tail envelope, the speckle model (bright points are "
        "covered by the oracle only) "
        f"and the assumed SNR."
    )


def report(res):
    """The benchmark's tables, generated from the summary."""
    n = f"{FROZEN['snr_db_nominal']:g}"
    L = [f"# P2-36 benchmark: generated report", '', f"Frozen parameters hash `{FROZEN_HASH}`; frozen frequency "
         f"{res['f_star_hz']:.0f} Hz ({res['frequency_rule']}). SNR {n} dB assumed. States: near chance (TV < 0.05), "
         f"95% at 5% excluded (TV < 0.9), unresolved.", '']
    for exc in ('quiet', 'strong'):
        L += [f"## {exc}", '', '| pair | L1 single image, TV | state | growth to 0.9 | L1 within 70 m | oracle, TV | state | '
              'SNR for 0.9 | oracle within 70 m | oracle within 39 m | achieved (score test) | tail share |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for k, r in res[exc].items():
            ach = r['achieved_score_test']
            a = ach.get('tv_achieved_best_phase', ach.get('tv_achieved'))
            L.append(f"| {k} | {_fmt(r['l1'][n]['tv_upper'])} | {r['l1'][n]['state']} | {_fmt(r['l1_growth_to_target'])} | "
                     f"{_fmt(r['l1_disc']['tv_upper'])} | {_fmt(r['oracle'][n]['tv_upper'])} | {r['oracle'][n]['state']} | "
                     f"{_fmt(r['oracle_snr_db_for_target'])} dB | {_fmt(r['oracle_disc']['tv_upper'])} | "
                     f"{_fmt(r['oracle_obs_disc']['tv_upper'])} | {_fmt(a)} | {100 * r['tail_share_of_integral']:.0f}% |")
        L += ['', f"SNR sweep (L1 and oracle, TV):", '', '| pair | ' + ' | '.join(f'{d:g} dB' for d in FROZEN['snr_db_sweep'])
              + ' |', '|---|' + '---|' * len(FROZEN['snr_db_sweep'])]
        for k, r in res[exc].items():
            L.append(f"| {k} | " + ' | '.join(f"{_fmt(r['l1'][f'{d:g}']['tv_upper'])} / {_fmt(r['oracle'][f'{d:g}']['tv_upper'])}"
                                             for d in FROZEN['snr_db_sweep']) + ' |')
        L.append('')
    L += ['## strong: allowances', '', '| pair | allowance | L1 | oracle |', '|---|---|---|---|']
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
        L += ['', f"## strong: the score test on synthesised images (presence, room; {mc['seeds']} grounds, a cyclic patch "
              f"{mc['patch_m'][0]:.0f} × {mc['patch_m'][1]:.0f} m)", '', '| amplification | AUC empirical | AUC predicted | '
              'shift / sd | predicted deflection | found at 5%, empirical | predicted | patch certificate |',
              '|---|---|---|---|---|---|---|---|']
        for x in mc['rows']:
            L.append(f"| ×{x['amplification']:g} | {x['auc_empirical']:.3f} | {x['auc_predicted']:.3f} | "
                     f"{x['mean_shift_over_sd']:.2f} | {x['deflection_predicted']:.3f} | {x['tpr_at_005_empirical']:.2f} | "
                     f"{x['tpr_at_005_predicted']:.2f} | {_fmt(x['certificate_tv_upper_patch'])} |")
    L += ['', '## strong: over the FTA band (every 2 Hz)', '', '| f (Hz) | ' + ' | '.join(res['curves']) + ' |',
          '|---|' + '---|' * len(res['curves'])]
    names = list(res['curves'])
    for i, row in enumerate(res['curves'][names[0]]):
        L.append(f"| {row['f_hz']:.0f} | " + ' | '.join(
            f"{_fmt(res['curves'][k][i]['l1_tv_upper'])} / {_fmt(res['curves'][k][i]['oracle_tv_upper'])}" for k in names) + ' |')
    L += ['', 'Each cell of the band table: L1 / oracle at the nominal SNR.', '']
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
