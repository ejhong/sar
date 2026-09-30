"""P2-29 · The oracle: what an ideal detector could learn, told more than any method knows.

    uv run python experiments/p2_29_oracle.py

P2-25 bounds every method under a model of ordinary ground (fully developed speckle, no reference image). This bounds
an oracle that needs no model of the ground's texture at all. It is told the exact reflectivity of every scatterer, the
background motion, the realisation of the shaking and both physical models (with the room and without); it lacks only
which of the two is true, and the receiver noise. Given all that, the raw echoes under the two hypotheses are Gaussian
with the same covariance (the noise) and means that differ by the echo of the room-specific motion alone, and the best
test's performance is exact, with no Pinsker step:

    detection rate - false-alarm rate <= TV = 2 Phi(Delta / 2) - 1,    Delta^2 = 2 |mu_1 - mu_0|^2 / sigma^2

(white circular receiver noise, variance sigma^2 per raw sample). Any method, told less, does no better: conditional on
what the oracle knows, the method is a test on the same echoes, and averaging over what it does not know cannot raise its
TV above the average of the oracle's; the oracle's TV is concave in Delta^2, so E TV <= 2 Phi(sqrt(E Delta^2) / 2) - 1.

E Delta^2 needs only that the scatterers' phases are independent and uniform (any amplitudes, any texture, bright points
included), for then the cross terms vanish:

    E Delta^2 = 2 sum_j (E_j / sigma^2) < |exp(-i k0 d1_j(t)) - exp(-i k0 d0_j(t))|^2 > = 8 sum_j SNR_j < sin^2(k0 delta_j(t) / 2) >

with E_j scatterer j's echo energy, the average over its pulses, and delta = d1 - d0 the room-specific motion: the
background motion, however large, cancels exactly. No Doppler-to-time relation enters; the echoes are pulse by pulse.
Over distributed ground sum_j SNR_j = SNR x cells, SNR the image's signal-to-noise ratio per resolution cell (sigma0 / NESZ),
raised by the processed band's share since the raw echoes hold the whole illuminated band.

1. The formula against pulse-by-pulse raw echoes (sarsim.echo's geometry) on small scenes: random scatterers on the real
   dwell, a room-specific motion of 0.8 rad and a background motion of the same size, and ten times it; Delta^2 exactly,
   realisation by realisation, against the formula.
2. The chamber: P2-25's cases (the static imprint under the ambient levels; P2-26's dynamic imprint under the urban
   background and a truck over a bump 15 m away, the scattered wave carried unattenuated across the scene, and within 32 m
   alone), with SNR 30 dB, above the brightest natural ground (sigma0 about 0 dB, assumed) over the best Dwell NESZ
   (-26.7 dB, ICEYE's product specification); and a corner reflector 50 dB over the ground on the imprint's peak.
3. Margins: how many times the room-specific motion would have to grow for the oracle to reach 95% found at 5% false
   alarms (Delta = 3.29; the phase saturating only raises it), and what survives stated, assumed allowances (the local level
   ten times the regional, site amplification three times, SNR 10 dB higher). A case whose margin does not survive is open.
4. Depth: the bench room 30 m down against rooms 15 m down, small and spread so that their marks on the ground match it
   (a non-negative combination of one small room's kernel, the lab's solver settling under strain as in P2-28), and the
   dilute superposition checked by solving two small rooms together. What is left of the imprint to tell the depths apart
   bounds, for any method, how much more signal depth needs than detection.
"""
import copy
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from numba import njit, prange
from scipy.optimize import nnls
from scipy.special import ndtr

from katabasis.runs import RESULTS, Run, load
from sarsim import information as inf
from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup, window

RID = 'p2_29_oracle'
SITES = Path(__file__).resolve().parents[2] / 'sites'
SNR_DB = 30.0
NESZ_BEST_DB = -26.7                  # ICEYE product specification, Dwell (Table 2-11), best of its range
REFLECTOR_DB = 50.0                   # over the ground, as P2-23
D_TARGET = 3.29                       # Delta for 95% found at 5% false alarms (2 Phi(D / 2) - 1 = 0.9)
ALLOWANCES = {'local ambient level over the regional': 10.0, 'site amplification': 3.0,
              'SNR 10 dB higher': float(np.sqrt(10.0))}
CHECK_PRF = 400.0
CHECK_REALISATIONS = 12
SMALL = {'centre': [0, 0, -15], 'size': [2, 2, 2], 'fill': 'air'}
SHIFT_CELLS = 3                       # basis rooms 6 m apart (the receivers are 2 m apart)
SHIFT_MAX = 6                         # to 36 m from the axis
INNER_M = 40.0


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tv_exact(delta):
    return float(2 * ndtr(np.asarray(delta, float) / 2) - 1)


# ------------------------------------------------------------------ 1. the formula against raw echoes

@njit(parallel=True, fastmath=False, cache=True)
def _echo_pair(t, V_eff, R0, u, r, a_re, a_im, d0, d1, lam, rho, r_lo, dr_s, n_s, half):
    """Range-compressed echoes of the same scatterers moved by d0 [pulse, scatterer] and by d1, their difference's
    energy, and the same without cross terms (each scatterer's own difference energy)."""
    n_p = t.shape[0]
    k4 = 4.0 * np.pi / lam
    diff_e = np.zeros(n_p)
    incoh = np.zeros(n_p)
    for n in prange(n_p):
        buf = np.zeros(n_s, np.complex128)
        x_p = V_eff * t[n]
        acc = 0.0
        for k in range(u.shape[0]):
            dx = x_p - u[k]
            rr = R0 + r[k]
            R = np.sqrt(rr * rr + dx * dx)
            c0 = complex(a_re[k], a_im[k]) * np.exp(-1j * k4 * (R + d0[n, k]))
            c1 = complex(a_re[k], a_im[k]) * np.exp(-1j * k4 * (R + d1[n, k]))
            dc = c1 - c0
            pos = (R - r_lo[n]) / dr_s
            i0 = int(np.floor(pos))
            for i in range(i0 - half, i0 + half + 1):
                if i < 0 or i >= n_s:
                    continue
                x = (i - pos) * dr_s / rho
                w = 1.0 if abs(x) < 1e-12 else np.sin(np.pi * x) / (np.pi * x)
                buf[i] += dc * w
                acc += (dc.real * dc.real + dc.imag * dc.imag) * w * w
        diff_e[n] = np.sum(buf.real * buf.real + buf.imag * buf.imag)
        incoh[n] = acc
    return diff_e.sum(), incoh.sum()


@njit(parallel=True, fastmath=False, cache=True)
def _energy(t, V_eff, R0, u, r, rho, r_lo, dr_s, n_s, half):
    """Each scatterer's echo energy per unit amplitude, summed over the pulses."""
    out = np.zeros(u.shape[0])
    for k in prange(u.shape[0]):
        acc = 0.0
        for n in range(t.shape[0]):
            dx = V_eff * t[n] - u[k]
            rr = R0 + r[k]
            R = np.sqrt(rr * rr + dx * dx)
            pos = (R - r_lo[n]) / dr_s
            i0 = int(np.floor(pos))
            for i in range(i0 - half, i0 + half + 1):
                if i < 0 or i >= n_s:
                    continue
                x = (i - pos) * dr_s / rho
                w = 1.0 if abs(x) < 1e-12 else np.sin(np.pi * x) / (np.pi * x)
                acc += w * w
        out[k] = acc
    return out


def check_raw(g):
    """Delta^2 sigma^2 / 2 = |mu_1 - mu_0|^2 from the echoes, against sum_j E_j <4 sin^2(k0 delta_j / 2)>."""
    st = EchoSetup.from_geometry(g, prf=CHECK_PRF)
    t = st.times
    k0 = 4 * np.pi / g.lam
    Lu, Lr = 3.0, 1.5
    n = 2 * int((Lu / g.resolution) * (Lr / st.rho))
    r_lo, dr_s, n_s = window(st, t, 3.0, Lu, Lr)
    half = int(np.ceil(4 * st.rho / dr_s))
    L, f_room, f_bg = 0.5, 5.0, 0.2
    D = 0.8 / k0
    rows = []
    for bg_scale in (1.0, 10.0):
        for rep in range(CHECK_REALISATIONS):
            rng = np.random.default_rng(700 + rep)
            u, r = rng.uniform(-Lu / 2, Lu / 2, n), rng.uniform(-Lr / 2, Lr / 2, n)
            a = np.exp(2j * np.pi * rng.random(n)) * rng.rayleigh(1 / np.sqrt(2), n) * np.sqrt(2)
            tt = np.broadcast_to(t[:, None], (len(t), n))                  # each pulse's absolute time
            bg = bg_scale * D * (np.sin(2 * np.pi * f_bg * tt + 0.7) + 0.3 * u[None, :] * np.cos(2 * np.pi * 0.37 * tt))
            room = D * np.exp(-(u ** 2 + r ** 2) / (2 * L * L))[None, :] * np.sin(2 * np.pi * f_room * t[:, None] + 0.3)
            ex, inc = _echo_pair(t, st.V_eff, st.R0, u, r, a.real.copy(), a.imag.copy(), bg, bg + room, st.lam, st.rho,
                                 r_lo, dr_s, n_s, half)
            E = _energy(t, st.V_eff, st.R0, u, r, st.rho, r_lo, dr_s, n_s, half) * np.abs(a) ** 2
            formula = float(np.sum(E * np.mean(4 * np.sin(k0 * room / 2) ** 2, axis=0)))
            rows.append({'background_scale': bg_scale, 'realisation': rep, 'exact_over_formula': ex / formula,
                         'no_cross_terms_over_formula': inc / formula})
        sel = [x['exact_over_formula'] for x in rows if x['background_scale'] == bg_scale]
        print(f"  raw echoes, background x{bg_scale:g}: exact / formula {np.mean(sel):.4f} +- "
              f"{np.std(sel) / np.sqrt(len(sel)):.4f}", flush=True)
    summ = {}
    for bg_scale in (1.0, 10.0):
        sel = np.array([x['exact_over_formula'] for x in rows if x['background_scale'] == bg_scale])
        summ[f'background_x{bg_scale:g}'] = {'mean': float(sel.mean()), 'se': float(sel.std() / np.sqrt(len(sel))),
                                             'spread': float(sel.std())}
    return {'scatterers': n, 'pulses': len(t), 'prf_hz': CHECK_PRF, 'patch_m': [Lu, Lr], 'room_width_m': L,
            'room_phase_rad': 0.8, 'room_f_hz': f_room, 'background_phase_rad': 0.8, 'rows': rows, 'summary': summ}


# ------------------------------------------------------------------ 2. the chamber; 3. margins

def near_bounds(maps, name, g):
    """4 N <Phi^2> per unit incident velocity amplitude, per frequency and direction, for P2-26's map within 39 m
    (tapered as P2-25 tapers it): 2 cells k0^2 integral |K|^2 dA, K = H / (2 pi i f)."""
    xs, ys, fs, H = maps[f'{name}_x'], maps[f'{name}_y'], maps[f'{name}_f'], maps[f'{name}_H']
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    tp = np.clip((39.0 - np.hypot(X, Y)) / 7.0, 0, 1)
    dA = float(xs[1] - xs[0]) ** 2
    k0 = 4 * np.pi / g.lam
    K2 = (np.abs(H) ** 2 * (tp ** 2)[None, :, :, None]).sum(axis=(1, 2)) * dA / (2 * np.pi * fs[None, :]) ** 2
    return fs, 2 * inf.cells_per_m2(g) * k0 ** 2 * K2.T            # [f, direction]


def chamber(g, p225, maps, amb):
    snr = 10 ** (SNR_DB / 10) / g.band_frac
    k0 = 4 * np.pi / g.lam
    rows = []
    add = lambda label, d2, kind, note='': rows.append({'case': label, 'kind': kind, 'delta': float(np.sqrt(d2)),
                                                        'tv_exact': tv_exact(np.sqrt(d2)),
                                                        'tv_pinsker': float(min(np.sqrt(d2) / 2, 1.0)),
                                                        'growth_to_target': float(D_TARGET / np.sqrt(d2)), 'note': note})
    for r in p225['static']:
        add(f"{r['case']} (static imprint)", 0.5 * snr * r['fisher_bound_any_frequency'], 'ambient')
    cul = amb['cultural']
    v_urban, band = cul['urban_background']['value'], cul['band_hz']['value']
    v_truck = cul['bus_or_truck_over_bump']['value']
    for d in p225['dynamic']:
        fs, near = near_bounds(maps, d['case'], g)
        far = np.asarray(d['far_bound_per_v2'])
        nd = near.shape[1]
        sel = (fs >= band[0]) & (fs <= band[1])
        df = float(np.mean(np.diff(fs)))
        amp2 = 2 * v_urban ** 2 / (band[1] - band[0]) * df / nd
        tot = (np.sqrt(near) + np.sqrt(far)[:, None]) ** 2
        add(f"urban background, all pass ({d['case']})", 0.5 * snr * amp2 * tot[sel].sum(), 'cultural',
            'scattered wave unattenuated across the scene')
        add(f"urban background, all pass, within 32 m ({d['case']})", 0.5 * snr * amp2 * near[sel].sum(), 'cultural')
        k, j = np.unravel_index(np.argmax(np.where(sel[:, None], tot, 0)), tot.shape)
        add(f"truck over a bump 15 m away, all pass, at {fs[k]:.0f} Hz ({d['case']})", 0.5 * snr * 2 * v_truck ** 2 * tot[k, j],
            'cultural', 'scattered wave unattenuated across the scene')
        kn, jn = np.unravel_index(np.argmax(np.where(sel[:, None], near, 0)), near.shape)
        add(f"truck over a bump 15 m away, all pass, within 32 m, at {fs[kn]:.0f} Hz ({d['case']})",
            0.5 * snr * 2 * v_truck ** 2 * near[kn, jn], 'cultural')
    # a corner reflector on the imprint's peak: its own SNR, <4 sin^2(k0 d / 2)> <= (k0 d)^2 / 2 for d = D cos(.)
    for p in p225['points']:
        if p['scr_db'] < 40:
            continue
        snr_pt = 10 ** ((SNR_DB + REFLECTOR_DB) / 10) / g.band_frac
        add(f"corner reflector on the peak, {p['shaking']}", 2 * snr_pt * (k0 * p['peak_displacement_m']) ** 2 / 2,
            'helped')
    allowance = float(np.prod(list(ALLOWANCES.values())))
    for r in rows:
        r['growth_after_allowances'] = r['growth_to_target'] / allowance
        r['open'] = bool(r['growth_after_allowances'] <= 1.0)
    return rows, snr, allowance


# ------------------------------------------------------------------ 4. depth

def los_maps(K, los, n):
    """The three strain components' line-of-sight maps on the receiver grid [3, n, n] (x east first, y north)."""
    xy = K['xy']
    xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
    ix, iy = np.searchsorted(xs, xy[:, 0]), np.searchsorted(ys, xy[:, 1])
    out = np.zeros((3, len(xs), len(ys)))
    for c, name in enumerate(('Kxx', 'Kyy', 'Kxy')):
        out[c, ix, iy] = -(K[name] @ los)
    return out, xs, ys


def shifted(M, di, dj):
    out = np.zeros_like(M)
    n0, n1 = M.shape[1:]
    si = slice(max(di, 0), n0 + min(di, 0))
    sj = slice(max(dj, 0), n1 + min(dj, 0))
    ti = slice(max(-di, 0), n0 + min(-di, 0))
    tj = slice(max(-dj, 0), n1 + min(-dj, 0))
    out[:, si, sj] = M[:, ti, tj]
    return out


def direction_weights():
    """The mean over wave azimuths of |Kxx s^2 + Kyy c^2 + Kxy s c|^2 as a quadratic form in the three maps: its
    Cholesky factor, so that least squares in the transformed maps is least squares averaged over directions."""
    Q = np.array([[3 / 8, 1 / 8, 0.0], [1 / 8, 3 / 8, 0.0], [0.0, 0.0, 1 / 8]])
    return np.linalg.cholesky(Q).T


def depth(g, host, m28):
    los = np.asarray(g.los_enu)
    K30, _ = m28.room_kernels('hollow_30', m28.ROOMS['hollow_30'], host)
    Ks, info_s = m28.room_kernels('small_15', SMALL, host)
    T, xs, ys = los_maps(K30, los, None)
    S, _, _ = los_maps(Ks, los, None)
    W = direction_weights()
    inner = (np.abs(xs)[:, None] <= INNER_M) & (np.abs(ys)[None, :] <= INNER_M)
    vec = lambda M: np.einsum('ij,jab->iab', W, M)[:, inner].ravel()
    shifts = [(i, j) for i in range(-SHIFT_MAX, SHIFT_MAX + 1) for j in range(-SHIFT_MAX, SHIFT_MAX + 1)]
    B = np.stack([vec(shifted(S, i * SHIFT_CELLS, j * SHIFT_CELLS)) for i, j in shifts], 1)
    b = vec(T)
    w, res = nnls(B, b, maxiter=50 * B.shape[1])
    delta = float(res / np.linalg.norm(b))
    # the single best small room at the axis, and the 15 m bench room alone, for comparison
    K15, _ = m28.room_kernels('hollow_15', m28.ROOMS['hollow_15'], host)
    T15, _, _ = los_maps(K15, los, None)
    b15 = vec(T15)
    best_single = float(np.linalg.norm(b - b15 * (b15 @ b) / (b15 @ b15)) / np.linalg.norm(b))
    # the dilute superposition: two small rooms 6 m apart solved together against the sum of the two alone
    P, _, _ = los_maps(pair_kernels(m28, host), los, None)
    two = S + shifted(S, SHIFT_CELLS, 0)
    interaction = float(np.linalg.norm(vec(P) - vec(two)) / np.linalg.norm(vec(P)))
    vol_small = float(np.prod(SMALL['size']))
    used = w > 1e-6 * w.max()
    return {'residual_fraction': delta, 'signal_factor_depth_over_detection': 1 / delta,
            'best_single_bench_room_residual': best_single, 'rooms_used': int(used.sum()),
            'largest_weight_small_rooms': float(w.max()), 'total_volume_m3': float(w.sum() * vol_small),
            'volume_30m_room_m3': float(np.prod(m28.ROOMS['hollow_30']['size'])),
            'weights': [{'east_m': float(i * SHIFT_CELLS * 2.0), 'north_m': float(j * SHIFT_CELLS * 2.0), 'small_rooms': float(x)}
                        for (i, j), x in zip(shifts, w) if x > 1e-6 * w.max()],
            'pair_interaction_fraction': interaction, 'small_room_runs': info_s, 'inner_m': INNER_M,
            'spacing_m': SHIFT_CELLS * 2.0}


def pair_kernels(m28, host):
    """Two small rooms 6 m apart, at the axis and 6 m east, solved together (the dilute superposition's check)."""
    def media_pair():
        from katabasis.compose import Grid, load_site, voxelise
        from katabasis.compose.site import parse_site
        from katabasis.seismic.elastic3d import Medium
        site = load_site('bench-void')
        raw = copy.deepcopy(site.raw)
        f0 = raw['features'][0]
        feats = []
        for k, x in enumerate((0.0, 2.0 * SHIFT_CELLS)):
            f = copy.deepcopy(f0)
            f['id'] = f'small_{k}'
            f['shape'] = dict(f['shape'], centre=[float(x), 0.0, SMALL['centre'][2]], size=SMALL['size'])
            f['fill'] = 'air'
            feats.append(f)
        raw['features'] = feats
        g = Grid.covering((-m28.HALF_DOMAIN, m28.HALF_DOMAIN), (-m28.HALF_DOMAIN, m28.HALF_DOMAIN), (m28.BOTTOM, m28.TOP), m28.H)
        return g, Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))
    from katabasis.seismic.arrays import surface_grid
    g, full = media_pair()
    rec = surface_grid(g, full.solid, m28.HALF_REC, m28.SPACING)
    K1, _ = m28.cached('pair_15', lambda: m28.settle(full, rec, host, 'pair_15'))
    K0, _ = m28.cached('without', lambda: m28.settle(m28.media(None)[1], rec, host, 'without'))
    return {'xy': K1['xy'], **{k: K1[k] - K0[k] for k in ('Kxx', 'Kyy', 'Kxy')}}


def main():
    g = DwellGeometry.from_record('giza-20250827')
    amb = json.loads((SITES / 'ambient.json').read_text())
    p225 = load('p2_25_information_bound')
    maps = np.load(RESULTS / 'p2_26_imprint_spectrum' / 'maps.npz')
    m28 = module('p228', Path(__file__).with_name('p2_28_two_depths.py'))
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    rock = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    host = (val('rho') * val('vp') ** 2 - 2 * val('rho') * val('vs') ** 2, val('rho') * val('vs') ** 2)
    params = {'acquisition': 'giza-20250827', 'snr_db': SNR_DB, 'nesz_best_db': NESZ_BEST_DB,
              'nesz_source': 'ICEYE product specification, imaging modes, Table 2-11 (Dwell): -26.7 to -15.6 dB',
              'sigma0_bright_natural_db': 0.0, 'sigma0_status': 'assumed', 'reflector_db': REFLECTOR_DB,
              'target_delta': D_TARGET, 'allowances': ALLOWANCES, 'allowance_status': 'assumed',
              'check': {'prf_hz': CHECK_PRF, 'realisations': CHECK_REALISATIONS},
              'depth': {'small_room': SMALL, 'spacing_m': SHIFT_CELLS * 2.0, 'reach_m': SHIFT_MAX * SHIFT_CELLS * 2.0}}
    with Run(RID, 'The oracle: what an ideal detector could learn, told more than any method knows', params) as run:
        t0 = time.time()
        raw = check_raw(g)
        print(f"  raw-echo check: {time.time() - t0:.0f} s", flush=True)
        rows, snr, allowance = chamber(g, p225, maps, amb)
        for r in rows:
            print(f"  {r['case']}: Delta {r['delta']:.3g}, TV {r['tv_exact']:.3g}, grow x{r['growth_to_target']:.3g}"
                  f"{' (open)' if r['open'] else ''}", flush=True)
        t0 = time.time()
        dep = depth(g, host, m28)
        print(f"  depth: residual {dep['residual_fraction']:.3g} ({time.time() - t0:.0f} s)", flush=True)
        amb_rows = [r for r in rows if r['kind'] == 'ambient']
        cult = [r for r in rows if r['kind'] == 'cultural']
        opened = [r['case'] for r in rows if r['open']]
        worst = max(cult, key=lambda r: r['tv_exact'])
        sb = raw['summary']
        finding = (
            f"An oracle told the exact reflectivity, the background motion, the shaking's realisation and both physical models, "
            f"lacking only which is true and the receiver noise, is bounded exactly: its detection rate exceeds its false-alarm "
            f"rate by 2 Phi(Delta / 2) - 1, and no method told less does better, whatever the ground's texture. Pulse by "
            f"pulse on the real dwell the formula for Delta holds ({sb['background_x1']['mean']:.3f} +- "
            f"{sb['background_x1']['se']:.3f} of it; {sb['background_x10']['mean']:.3f} +- {sb['background_x10']['se']:.3f} "
            f"with the background ten times larger, which cancels). At {SNR_DB:.0f} dB per resolution cell, under the ambient "
            f"levels the bench room gives the oracle at most {max(r['tv_exact'] for r in amb_rows):.1e}: its imprint would "
            f"have to grow {min(r['growth_to_target'] for r in amb_rows):.1e} times, {min(r['growth_after_allowances'] for r in amb_rows):.1e} "
            f"after allowing the local level ten times the regional, site amplification three times and 10 dB more SNR. "
            f"The open regime is strong nearby shaking known exactly, the scattered wave carried unattenuated across the "
            f"scene: {worst['case']} gives the oracle {worst['tv_exact']:.2f}"
            + (f", enough within the model for 95% found at 5% false alarms (it would need only {worst['growth_to_target']:.2f} "
               f"of its motion)" if worst['growth_to_target'] < 1 else
               f", {worst['growth_to_target']:.1f} times short of 95% found at 5% false alarms")
            + f"; {len(opened)} cases are open after the allowances. "
            + (f"Depth is carried by the mark's shape in this model: rooms 15 m down, small and spread, leave "
               f"{100 * dep['residual_fraction']:.0f}% of the mark of a room 30 m down unmatched (the bench room 15 m down, "
               f"{100 * dep['best_single_bench_room_residual']:.0f}%), so telling those two depths apart takes only "
               f"{dep['signal_factor_depth_over_detection']:.1f} times the signal that detecting either takes; the "
               f"superposition of small rooms holds to {100 * dep['pair_interaction_fraction']:.1f}%."
               if dep['residual_fraction'] > 0.3 else
               f"Depth is another matter: rooms 15 m down, small and spread, match the mark of a room 30 m down to within "
               f"{100 * dep['residual_fraction']:.1f}% (the bench room itself only to "
               f"{100 * dep['best_single_bench_room_residual']:.0f}%), so telling the two depths apart takes "
               f"{dep['signal_factor_depth_over_detection']:.0f} times the signal that detecting either takes, for any method; "
               f"the superposition of small rooms holds to {100 * dep['pair_interaction_fraction']:.1f}%."))
        run.save({'raw_echo_check': raw, 'chamber': rows, 'snr_per_cell_raw': snr, 'allowance_product': allowance,
                  'open_cases': opened, 'depth': dep, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
