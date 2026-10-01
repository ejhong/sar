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
import copy
import json
import time
from pathlib import Path

import numpy as np

from katabasis.runs import RESULTS, Run, load, memo
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


# ------------------------------------------------------------------ the shafts: one open shaft by the static solver

SHAFT = {'radius_m': 5.0, 'length_claimed_m': 640.0,
         'positions': [[x, y] for y in (-45.0, 45.0) for x in (-67.5, -22.5, 22.5, 67.5)]}
SHAFT_SOLVER = {'h': 2.0, 'half_m': 120.0, 'bottom_m': -260.0, 'top_m': 8.0, 'rec_half_m': 90.0, 'rec_spacing_m': 2.0,
                'duration_s': 1.4, 'damping': 30.0, 'truncations_m': [100.0, 200.0]}
CACHE = RESULTS / 'cache' / 'p2_38'


def shaft_medium(length):
    """Flat ground of the bench's limestone (the pyramid omitted, its load and its scattering) with one open shaft of the
    claim's radius from the surface to `length` (None: no shaft), on the solver's grid."""
    from katabasis.compose import Grid, load_site, voxelise
    from katabasis.compose.site import parse_site
    from katabasis.seismic.elastic3d import Medium
    site = load_site('bench-void')
    raw = copy.deepcopy(site.raw)
    f0 = raw['features'][0]
    raw['features'] = []
    if length is not None:
        f = copy.deepcopy(f0)
        f['id'] = 'shaft'
        f['shape'] = {'type': 'cylinder', 'centre': [0.0, 0.0, -length / 2 + 2.0], 'radius': SHAFT['radius_m'],
                      'height': length + 4.0}                     # open to the air above the ground
        f['fill'] = 'air'
        raw['features'] = [f]
    sv = SHAFT_SOLVER
    g = Grid.covering((-sv['half_m'], sv['half_m']), (-sv['half_m'], sv['half_m']), (sv['bottom_m'], sv['top_m']), sv['h'])
    return g, Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))


def shaft_settle(length, comp):
    """The surface's displacement (east, north, up) per unit horizontal strain `comp` ('xx' or 'xy', the tensor
    component) when the strained ground with the shaft comes to rest (P2-28's loading, relax_strain); kept on disk."""
    import hashlib
    key = json.dumps({'solver': SHAFT_SOLVER, 'radius': SHAFT['radius_m'], 'length': length, 'comp': comp}, sort_keys=True)
    path = CACHE / f"shaft_{hashlib.sha1(key.encode()).hexdigest()[:16]}.npz"
    if path.is_file():
        z = np.load(path)
        return z['xy'], z['u'], json.loads(str(z['info']))
    from katabasis.seismic.arrays import surface_grid
    m28 = module('p228', Path(__file__).with_name('p2_28_two_depths.py'))
    m36 = module('p236', Path(__file__).with_name('p2_36_benchmark.py'))
    g, med = shaft_medium(length)
    rec = surface_grid(g, med.solid, SHAFT_SOLVER['rec_half_m'], SHAFT_SOLVER['rec_spacing_m'])
    m28.DURATION = SHAFT_SOLVER['duration_s']
    t0 = time.time()
    u, drift = m28.relax_strain(med, {comp: m28.STRAIN}, rec, m36.host_moduli(), damping=SHAFT_SOLVER['damping'])
    info = {'runtime_s': round(time.time() - t0, 1), 'last_step_change': drift, 'cells': int(np.prod(g.shape))}
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(path, xy=rec[:, :2], u=u / m28.STRAIN, info=json.dumps(info))
    print(f"  shaft {length} m, eps_{comp}: {info['runtime_s']:.0f} s (last step changed {drift:.1e})", flush=True)
    return rec[:, :2], u / m28.STRAIN, info


def shaft_solves():
    out = {}
    for comp in ('xx', 'xy'):
        xy, u0, i0 = shaft_settle(None, comp)
        for L in SHAFT_SOLVER['truncations_m']:
            if comp == 'xy' and L != max(SHAFT_SOLVER['truncations_m']):
                continue
            _, u1, i1 = shaft_settle(L, comp)
            out[(comp, L)] = {'xy': xy, 'K': u1 - u0, 'info': {'shaft': i1, 'without': i0}}
    return out


# ------------------------------------------------------------------ fields at the microseisms' frequency

NEAR_HALF_M = 200.0          # the near grid round the shafts, 2 m
NEAR_STEP_M = 2.0
R_CAL = (60.0, 84.0)         # the ring where the solver's near field calibrates the column model (its receivers reach 90 m)
COLUMN_DZ = 5.0
TAPER_M = 10.0               # the column model's smoothing; read only beyond R_CAL
Q_DYNAMIC = (100.0, 1000.0)  # the rock's damping in the dynamic field: realistic, and nearly none
PASSES = ('giza-20250827', 'giza-20220715')


def props():
    lam, mu, nu, E = rock()
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    r = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    val = lambda q: r[q]['value'] if isinstance(r[q], dict) else r[q]
    return {'vp': val('vp'), 'vs': val('vs'), 'rho': val('rho'), 'lam': lam, 'mu': mu, 'nu': nu}


def unit_strains(lam, mu):
    """The three unit horizontal strains (eps_xx, eps_yy, the tensor eps_xy) in plane stress, as P2-04 and P2-28 load."""
    zz = -lam / (lam + 2 * mu)
    return [np.diag([1.0, 0.0, zz]), np.diag([0.0, 1.0, zz]), np.array([[0, 1.0, 0], [1.0, 0, 0], [0, 0, 0]])]


def profiles(kind, omega, Q, ks, rs, length=None):
    """Radial profiles [3 strains][orders, 3, r] of the deep void (a sphere of the structure's volume at its depth) or of
    one shaft (a column of Eshelby cylinder moments every COLUMN_DZ to `length`), per unit strain."""
    from katabasis.seismic.analytic import void_moment_from_strain
    from katabasis.seismic.halfspace import radial_coefficients, radial_profiles
    pr = props()
    cns = []
    for eps in unit_strains(pr['lam'], pr['mu']):
        if kind == 'void':
            st = structure()
            M = void_moment_from_strain(eps, float(np.prod(st['size'])), pr['lam'], pr['mu'], 'sphere')
            cns.append(radial_coefficients(M, -st['centre'][2], omega, pr['vp'], pr['vs'], pr['rho'], Q, ks))
        else:
            m = void_moment_from_strain(eps, np.pi * SHAFT['radius_m'] ** 2, pr['lam'], pr['mu'], 'cylinder')
            z = np.arange(COLUMN_DZ / 2, length, COLUMN_DZ)
            cns.append(radial_coefficients(np.repeat((m * COLUMN_DZ)[None], len(z), 0), z, omega, pr['vp'], pr['vs'],
                                           pr['rho'], Q, ks, taper_m=TAPER_M))
    return radial_profiles(cns, ks, rs)


def radial_grid(r_max):
    return np.unique(np.concatenate([np.arange(0.0, 400.0, 0.5), np.arange(400.0, r_max + 5.0, 5.0)]))


def solver_near(sol):
    """Interpolators for one shaft's static near field per unit strain (xx and xy from the solver, yy by the shaft's
    symmetry under a quarter turn), each giving [3, ...] at (x, y) relative to the shaft."""
    from scipy.interpolate import RegularGridInterpolator
    L = max(SHAFT_SOLVER['truncations_m'])
    out = []
    for comp in ('xx', 'xy'):
        d = sol[(comp, L)]
        xy, K = d['xy'], d['K']
        xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
        A = np.zeros((len(xs), len(ys), 3))
        A[np.searchsorted(xs, xy[:, 0]), np.searchsorted(ys, xy[:, 1])] = K
        out.append(RegularGridInterpolator((xs, ys), A, bounds_error=False, fill_value=0.0))
    fxx, fxy = out

    def f_xx(x, y):
        return np.moveaxis(fxx(np.stack([x, y], -1)), -1, 0)

    def f_yy(x, y):                       # K_yy(p) = R90 K_xx(R90^T p): R90^T (x, y) = (y, -x); R90 (e, n) = (-n, e)
        k = f_xx(y, -x)
        return np.stack([-k[1], k[0], k[2]])

    def f_xy(x, y):
        return np.moveaxis(fxy(np.stack([x, y], -1)), -1, 0)
    return [f_xx, f_yy, f_xy]


def shaft_field(X, Y, F, rs, near, c_cal):
    """[3 strains][3, ...] at (X, Y): the eight shafts, each the solver's near field within R_CAL[1] and the column model
    times c_cal beyond, summed in phase (they are 135 m apart under a wave 8 km long)."""
    from katabasis.seismic.halfspace import field_from_profiles
    out = [np.zeros((3,) + X.shape, complex) for _ in range(3)]
    for (x0, y0) in SHAFT['positions']:
        dx, dy = X - x0, Y - y0
        inner = np.hypot(dx, dy) < R_CAL[1]
        for c in range(3):
            u = c_cal * field_from_profiles(F[c], rs, dx, dy)
            if inner.any():
                u[:, inner] = near[c](dx[inner], dy[inner])
            out[c] += u
    return out


def void_field(X, Y, F, rs):
    from katabasis.seismic.halfspace import field_from_profiles
    st = structure()
    return [field_from_profiles(F[c], rs, X - st['centre'][0], Y - st['centre'][1]) for c in range(3)]


def los_maps(u3, los):
    """The 24 directions' line-of-sight maps (slant-range increase per unit strain amplitude), complex."""
    L = [-np.tensordot(los, u, axes=(0, 0)) for u in u3]
    az = np.deg2rad(np.arange(24) * 7.5)
    return [np.sin(p) ** 2 * L[0] + np.cos(p) ** 2 * L[1] + np.sin(p) * np.cos(p) * L[2] for p in az]


def frame_for(g, name):
    import json as _j
    from sarsim import information as inf
    src = _j.loads((SITES / 'acquisitions' / f'{name}.json').read_text())['source']
    n_along, n_lines = src['shape']
    dg = g.dr / np.sin(g.theta)
    rho = (np.arange(n_lines) - n_lines / 2 + 0.5) * dg
    in_a, _, _ = inf.band_masks(g, (n_along, 1))
    return {'rho': rho, 'dg': dg, 'half_along_m': n_along * g.dx / 2, 'n_along': int(n_along),
            'n_band_along': int(in_a.sum()), 'n_lines': int(n_lines)}


def brightness(name, X, Y, ground_db):
    """sigma0 over the site median at (X, Y) from the pass's 2 km map of summed power (5 m cells in image geometry,
    placed through its tie points); NaN outside the map or the image."""
    from scipy.interpolate import RegularGridInterpolator
    z = np.load(SITES / 'acquisitions' / f'{name}-khafre-power-2km.npz')
    tx, ty = z['tie_x_m'][0], z['tie_y_m'][:, 0]
    pts = np.stack([Y.ravel(), X.ravel()], 1)
    li = RegularGridInterpolator((ty, tx), z['tie_line'], bounds_error=False, fill_value=np.nan)(pts)
    sa = RegularGridInterpolator((ty, tx), z['tie_sample'], bounds_error=False, fill_value=np.nan)(pts)
    P = z['power_sum_dn2']
    i = np.floor((li - z['line0']) / z['lines_per_cell'])
    j = np.floor((sa - z['sample0']) / z['samples_per_cell'])
    ok = np.isfinite(i) & np.isfinite(j) & (i >= 0) & (i < P.shape[0]) & (j >= 0) & (j < P.shape[1])
    s0 = np.full(pts.shape[0], np.nan)
    ii, jj = i[ok].astype(int), j[ok].astype(int)
    s0[ok] = (z['calibration_factor'] * P[ii, jj] / (z['lines_per_cell'] * z['samples_per_cell'])
              * np.sin(np.deg2rad(z['incidence_deg'][jj])))
    cells = z['calibration_factor'] * P / (z['lines_per_cell'] * z['samples_per_cell']) * np.sin(np.deg2rad(z['incidence_deg']))[None, :]
    ref = 10 ** (ground_db / 10)
    return (s0 / ref).reshape(X.shape), {'map_max_over_median': float(cells.max() / ref),
                                         'map_p999_over_median': float(np.quantile(cells, 0.999) / ref),
                                         'map_median_db': float(10 * np.log10(np.median(cells))),
                                         'cells': int(cells.size)}


def main():
    from scipy.interpolate import RegularGridInterpolator
    from sarsim.acquisition import DwellGeometry
    from sarsim.radiometry import snr_per_cell
    from katabasis.seismic.analytic import rayleigh_speed
    from katabasis.seismic.halfspace import k_grid
    m36 = module('p236', Path(__file__).with_name('p2_36_benchmark.py'))
    p30 = load('p2_30_finite_certificates')
    micro_env = p30['background_included']['los_peak_envelope_m']   # sqrt(1 + hv^2) x vertical peak: any line of sight
    ql = m36.quiet_level()
    eps0 = ql['strain_amplitude']
    omega = 2 * np.pi * ql['f_hz']
    pr = props()
    vr = rayleigh_speed(pr['vp'], pr['vs'])
    st = structure()
    depth = -st['centre'][2]
    params = {'structure': st, 'shafts': SHAFT, 'shaft_solver': SHAFT_SOLVER, 'grid_half_m': GRID_HALF_M,
              'grid_step_m': GRID_STEP_M, 'near_half_m': NEAR_HALF_M, 'near_step_m': NEAR_STEP_M, 'r_cal_m': R_CAL,
              'column_dz_m': COLUMN_DZ, 'taper_m': TAPER_M, 'q_dynamic': Q_DYNAMIC, 'f_hz': ql['f_hz'],
              'rock': pr, 'cube_allowance': CUBE_ALLOWANCE, 'allowances': ALLOWANCES, 'allowance_status': 'assumed',
              'quiet_level': ql, 'passes': PASSES,
              'model': "the deep structure as Eshelby's sphere of its volume, the shafts as columns of Eshelby cylinder "
                       "moments calibrated by the static solver near the surface, both through a damped half-space's "
                       "point-source response at the microseisms' frequency (katabasis.seismic.halfspace)",
              'benchmark_frozen_hash': m36.FROZEN_HASH}
    with Run(RID, 'The claim computed: the structures announced under Khafre', params) as run:
        t0 = time.time()
        rs = radial_grid(np.hypot(GRID_HALF_M, GRID_HALF_M))
        xf = np.arange(-GRID_HALF_M, GRID_HALF_M + 1e-9, GRID_STEP_M)
        Xf, Yf = np.meshgrid(xf, xf, indexing='ij')
        xn = np.arange(-NEAR_HALF_M, NEAR_HALF_M + 1e-9, NEAR_STEP_M)
        Xn, Yn = np.meshgrid(xn, xn, indexing='ij')
        # the deep structure: static (Okada) against dynamic, each damping
        # one grid for every damping: the union of each one's refinement about the pole and the branch points
        ks_v = np.unique(np.concatenate([k_grid(omega, pr['vp'], pr['vs'], vr, depth, rs[-1], Q, n_base=6000)
                                         for Q in Q_DYNAMIC]))
        Fv = {Q: memo(RID, f'void_{Q:g}', lambda Q=Q: profiles('void', omega, Q, ks_v, rs), __file__, version='void-v1',
                      inputs={'st': st, 'omega': omega, 'Q': Q, 'ks': ks_v, 'rs': rs, 'rock': pr})
              for Q in Q_DYNAMIC}
        print(f"  deep structure profiles ({time.time() - t0:.0f} s)", flush=True)
        # the shafts: the solver near the surface, the column model static (to the solver's truncation) and dynamic (to
        # the claimed length)
        sol = shaft_solves()
        near = solver_near(sol)
        ks_s = np.unique(np.concatenate([k_grid(omega, pr['vp'], pr['vs'], vr, TAPER_M, rs[-1], Q, k_top=5.0, n_base=10000)
                                         for Q in Q_DYNAMIC]))
        L_sol = max(SHAFT_SOLVER['truncations_m'])
        om_s = 0.01 * pr['vs'] / SHAFT['length_claimed_m']          # quasi-static: k_s L = 0.01
        Fs_static = memo(RID, 'shaft_static', lambda: profiles('shaft', om_s, 100.0, ks_s, rs, L_sol), __file__,
                         version='shaft-v1', inputs={'omega': om_s, 'ks': ks_s, 'rs': rs, 'L': L_sol, 'shaft': SHAFT,
                                                     'dz': COLUMN_DZ, 'taper': TAPER_M, 'rock': pr})
        Fs_static_full = memo(RID, 'shaft_static_full', lambda: profiles('shaft', om_s, 100.0, ks_s, rs,
                                                                          SHAFT['length_claimed_m']), __file__,
                              version='shaft-v1', inputs={'omega': om_s, 'ks': ks_s, 'rs': rs, 'shaft': SHAFT,
                                                          'dz': COLUMN_DZ, 'taper': TAPER_M, 'rock': pr})
        Fs = {Q: memo(RID, f'shaft_{Q:g}', lambda Q=Q: profiles('shaft', omega, Q, ks_s, rs, SHAFT['length_claimed_m']),
                      __file__, version='shaft-v1', inputs={'omega': omega, 'Q': Q, 'ks': ks_s, 'rs': rs, 'shaft': SHAFT,
                                                            'dz': COLUMN_DZ, 'taper': TAPER_M, 'rock': pr})
              for Q in Q_DYNAMIC}
        print(f"  shaft profiles ({time.time() - t0:.0f} s)", flush=True)
        # calibration: the solver's near field against the column model (same truncation, static) on the ring
        from katabasis.seismic.halfspace import field_from_profiles
        ang = np.deg2rad(np.arange(0, 360, 5.0))
        rr = np.linspace(R_CAL[0], R_CAL[1], 7)
        RA, AA = np.meshgrid(rr, ang, indexing='ij')
        cx, cy = RA * np.cos(AA), RA * np.sin(AA)
        ratios = {}
        for c, name in ((0, 'xx'), (2, 'xy')):
            us = near[c](cx, cy)
            uc = field_from_profiles(Fs_static[c], rs, cx, cy).real
            ratios[name] = float(np.sqrt(np.sum(us ** 2) / np.sum(uc ** 2)))
        c_cal = max(1.0, max(ratios.values()))
        d100, d200 = sol[('xx', 100.0)]['K'], sol[('xx', 200.0)]['K']
        rec = sol[('xx', 200.0)]['xy']
        ring = (np.hypot(rec[:, 0], rec[:, 1]) >= R_CAL[0]) & (np.hypot(rec[:, 0], rec[:, 1]) <= R_CAL[1])
        truncation = {'peak_100_over_200': float(np.abs(d100).max() / np.abs(d200).max()),
                      'ring_rms_100_over_200': float(np.sqrt(np.sum(d100[ring] ** 2) / np.sum(d200[ring] ** 2))),
                      'peak_m_per_strain_200': float(np.abs(d200).max())}
        print(f"  calibration: solver over column on the ring {ratios}; truncation 100/200 m {truncation}", flush=True)
        # the maps per unit strain, both grids, every model
        def maps_for(F_void, F_shaft, dyn=True):
            fv = [void_field(X, Y, F_void, rs) for X, Y in ((Xn, Yn), (Xf, Yf))]
            fs = [shaft_field(X, Y, F_shaft, rs, near, c_cal) for X, Y in ((Xn, Yn), (Xf, Yf))]
            return fv, fs
        fields = {Q: maps_for(Fv[Q], Fs[Q]) for Q in Q_DYNAMIC}
        print(f"  fields on the grids ({time.time() - t0:.0f} s)", flush=True)
        passes = {}
        for name in PASSES:
            g = DwellGeometry.from_record(name)
            fr = frame_for(g, name)
            los = np.asarray(g.los_enu)
            k0 = 4 * np.pi / g.lam
            osc_q = 1 / (2 * np.pi * ql['f_hz'] * m36.pass_seconds(g))
            snr = snr_per_cell(name)
            head = round(snr['headline_db'], 2)
            sweep = sorted({round(snr['measured_floor_lower_db'], 2), *[round(x['snr_db'], 2) for x in snr['scenarios']],
                            round(head, 2), round(snr['bright_db'], 2)})
            Qc_micro = fr['n_band_along'] * (np.sqrt(24) * k0 * micro_env) ** 2
            wn, binfo = brightness(name, Xn, Yn, snr['ground_sigma0_median_db'])
            wf, _ = brightness(name, Xf, Yf, snr['ground_sigma0_median_db'])
            w_rem = binfo['map_max_over_median']             # beyond the map, no brighter than its brightest cell
            wn_med, wf_med = np.where(np.isfinite(wn), wn, 1.0), np.where(np.isfinite(wf), wf, 1.0)   # or at the median
            wn, wf = np.where(np.isfinite(wn), wn, w_rem), np.where(np.isfinite(wf), wf, w_rem)

            def energies(Q, part):
                fv, fs = fields[Q]
                out = []
                for gi in range(2):
                    mv, ms = los_maps(fv[gi], los), los_maps(fs[gi], los)
                    if part == 'deep structure':
                        a = [np.abs(m) for m in mv]
                    elif part == 'shafts':
                        a = [np.abs(m) for m in ms]
                    else:                                    # both, their amplitudes added (an upper bound on the sum)
                        a = [np.abs(m1) + np.abs(m2) for m1, m2 in zip(mv, ms)]
                    e = eps0 ** 2 * np.mean(np.square(a), axis=0)
                    e_max = (eps0 / np.sqrt(len(a)) * np.sum(a, axis=0)) ** 2
                    out.append((e, e_max, float(eps0 * max(m.max() for m in a))))
                return out

            def lines(en, ef, w=None):
                wn_, wf_ = (1.0, 1.0) if w is None else w
                interp = RegularGridInterpolator((xf, xf), ef * wf_, bounds_error=False, fill_value=0.0)
                env = lambda E, N: interp(np.stack([E, N], -1))
                return m36.line_integrals(g, fr, en * wn_, xn, xn, NEAR_HALF_M - 10.0, env)

            res = {}
            for part in ('deep structure', 'shafts', 'claim (both)'):
                best = None
                for Q in Q_DYNAMIC:
                    (en, en_max, pk_n), (ef, ef_max, pk_f) = energies(Q, part)
                    S = lines(en, ef)
                    integ = float(S.sum() * fr['dg'])
                    if best is None or integ > best['integ']:
                        best = {'Q': Q, 'S': S, 'integ': integ, 'S_max': lines(en_max, ef_max),
                                'S_w': lines(en, ef, (wn, wf)), 'S_w_med': lines(en, ef, (wn_med, wf_med)),
                                'peak': max(pk_n, pk_f), 'rms_peak': float(np.sqrt(max(en.max(), ef.max()))),
                                'half_area': float((ef >= ef.max() / 2).sum() * GRID_STEP_M ** 2)}
                    best.setdefault('integ_by_q', {})[f'{Q:g}'] = integ
                integ_w = float(best['S_w'].sum() * fr['dg'])
                integ_w_med = float(best['S_w_med'].sum() * fr['dg'])
                Qm = m36.phase_energy(g, best['S'], 0.0)
                Qx = m36.phase_energy(g, best['S_max'], osc_q)
                # the share of the weighted integral from beyond the map (taken at its brightest cell)
                rows = {}
                allow = [('as computed', 1.0, 1.0, 1.0)]
                if part != 'shafts':
                    allow.append((f'cube allowance x{CUBE_ALLOWANCE:g} on the deep structure', CUBE_ALLOWANCE, 1.0, 1.0))
                amb = float(np.prod(list(ALLOWANCES.values())))
                allow.append(('and the ambient allowances (x10 local level, x3 site amplification), background too',
                              CUBE_ALLOWANCE if part != 'shafts' else 1.0, amb, amb))
                for label, cube, a_sig, a_bg in allow:
                    amp = a_sig * (cube if part == 'deep structure' else 1.0)
                    if part == 'claim (both)' and cube != 1.0:
                        # the cube scales only the deep structure's share: bound (sqrt(e_s) + cube sqrt(e_v))^2 by cube^2 e
                        amp = a_sig * cube
                    Qc = np.full_like(Qm, Qc_micro * a_bg ** 2)
                    r = {'signal_amplitude_factor': amp, 'background_amplitude_factor': a_bg, 'l1': {}, 'oracle': {},
                         'oracle_power_weighted': {}}
                    for sdb in sweep:
                        s2 = 10 ** (-sdb / 10)
                        r['l1'][f'{sdb:g}'] = m36.l1_ensemble(Qm * amp ** 2, Qx * amp ** 2, Qc, s2)
                        o = m36.oracle(g, best['integ'] * amp ** 2, sdb, 0.0)
                        o['state'] = m36.state(o['tv_upper'])
                        r['oracle'][f'{sdb:g}'] = o
                        ow = m36.oracle(g, integ_w * amp ** 2, sdb, 0.0)
                        ow['state'] = m36.state(ow['tv_upper'])
                        r['oracle_power_weighted'][f'{sdb:g}'] = ow
                    s2 = 10 ** (-head / 10)
                    r['l1_growth_to_target'] = m36.growth(lambda x: m36.l1_ensemble(Qm * (amp * x) ** 2, Qx * (amp * x) ** 2,
                                                                                    Qc, s2)['tv_upper'], m36.AN['target_tv'])
                    ow = r['oracle_power_weighted'][f'{head:g}']
                    r['oracle_growth_to_target'] = float(np.sqrt(m36.oracle_d2_target() / ow['e_delta2']))
                    rows[label] = r
                res[part] = {'q_used': best['Q'], 'integral_e_m4': best['integ'], 'integral_by_q_m4': best['integ_by_q'],
                             'integral_power_weighted_m4': integ_w, 'integral_power_weighted_beyond_at_median_m4': integ_w_med,
                             'oracle_power_weighted_beyond_at_median': m36.oracle(g, integ_w_med, head, 0.0),
                             'peak_los_m': best['peak'],
                             'rms_over_directions_peak_m': best['rms_peak'],
                             'half_power_area_m2': best['half_area'], 'rows': rows}
            # static against dynamic for the deep structure (Okada's field, as the first version)
            xs_s, maps_s, _ = imprint(g, st)
            e_s = eps0 ** 2 * np.mean(np.square(maps_s), axis=0)
            integ_static = float(m36.line_integrals(g, fr, e_s, xs_s, xs_s, 1e9, None).sum() * fr['dg'])
            passes[name] = {'snr': snr, 'snr_sweep_db': sweep, 'headline_db': head, 'brightness': {**binfo,
                            'beyond_map_taken_as': 'the brightest cell of the map'},
                            'microseism_common_phase_energy_per_line': Qc_micro, 'osc_window': osc_q,
                            'deep_structure_static_integral_m4': integ_static,
                            'deep_structure_dynamic_over_static': res['deep structure']['integral_e_m4'] / integ_static,
                            'parts': res}
            print(f"  {name}: deep structure dynamic/static integral "
                  f"{passes[name]['deep_structure_dynamic_over_static']:.2f}; claim L1 "
                  f"{res['claim (both)']['rows']['as computed']['l1'][f'{head:g}']['tv_upper']:.2e}, oracle (weighted) "
                  f"{res['claim (both)']['rows']['as computed']['oracle_power_weighted'][f'{head:g}']['tv_upper']:.2e} "
                  f"({time.time() - t0:.0f} s)", flush=True)
        bench = load('p2_36_benchmark')['quiet']['presence (room)']
        out = {'passes': passes, 'calibration': {'solver_over_column_on_ring': ratios, 'c_cal': c_cal,
                                                 'truncation': truncation},
               'shaft_solver_info': {f'{k[0]}_{k[1]:g}': v['info'] for k, v in sol.items()},
               'bench_room_quiet': {'l1_ensemble': bench['l1_ensemble'], 'oracle': bench['oracle']}}
        out['finding'] = finding(out)
        run.save(out)
        print(out['finding'])


def finding(out):
    P = out['passes']
    p25 = P['giza-20250827']
    h = f"{p25['headline_db']:g}"
    c = p25['parts']['claim (both)']['rows']
    base, worst = c['as computed'], list(c.values())[-1]
    dv = p25['parts']['deep structure']
    sh = p25['parts']['shafts']
    p22 = P.get('giza-20220715')
    t = out['calibration']['truncation']
    f22 = ''
    if p22:
        h22 = f"{p22['headline_db']:g}"
        b22 = p22['parts']['claim (both)']['rows']['as computed']
        f22 = (f" The 2022 pass (Dwell, its own {p22['headline_db']:.1f} dB): one image at most "
               f"{b22['l1'][h22]['tv_upper']:.1e}, the oracle {b22['oracle_power_weighted'][h22]['tv_upper']:.1e}.")
    return (
        f"The structures announced under Khafre, computed: the eight shafts (10 m across, layout representative, open at "
        f"the surface with the pyramid omitted, which is generous) and the deep structure (a void 80 m on a side at 1,220 m, "
        f"its size representative), loaded by Giza's regional microseisms at 0.2 Hz, strained at most as at the surface. "
        f"The shafts: the lab's static solver for one shaft near the surface (truncated at 100 and 200 m, the near field "
        f"changing by {abs(1 - t['ring_rms_100_over_200']) * 100:.1f}% between them), a column of Eshelby cylinder moments "
        f"to their full 640 m beyond, scaled by {out['calibration']['c_cal']:.2f} where the two meet; the deep structure as "
        f"Eshelby's sphere of its volume; both through a damped half-space's response at 0.2 Hz, so the field radiated to "
        f"the image's edge is counted (the deep structure's dynamic integral {p25['deep_structure_dynamic_over_static']:.2f} "
        f"times its static one). Imprint on the line of sight: the shafts peak at {sh['peak_los_m']:.1e} m, the deep "
        f"structure at {dv['peak_los_m']:.1e} m. On the 2025 pass (Dwell Fine, {p25['headline_db']:.1f} dB, the best "
        f"specified noise floor), any reader of the one image, bounded line by line with the microseisms in each line's "
        f"reference and averaged over the field's realisations: at most {base['l1'][h]['tv_upper']:.1e} above its "
        f"false-alarm rate ({base['l1'][h]['worst_realisation_tv_upper']:.1e} at the worst realisation); the oracle told "
        f"every scatterer's reflectivity, weighted by the image's own measured brightness over 2 km "
        f"(beyond, its brightest cell): {base['oracle_power_weighted'][h]['tv_upper']:.1e}. The motion would have to grow "
        f"{base['l1_growth_to_target']:.1e} times for one image and {base['oracle_growth_to_target']:.1e} for the oracle; "
        f"with a cube's moment 1.5 times the sphere's (the solver's cube measured 1.22 to 1.36 times the sphere for the "
        f"bench room), the local level ten times the regional and site amplification three times, applied to the "
        f"background as well, still {worst['l1_growth_to_target']:.1e} and {worst['oracle_growth_to_target']:.1e}."
        + f22 +
        f" These allowances are sensitivities, not enclosures of missing physics. Not computed: the pyramid's own load and "
        f"scattering, other excitations (a lorry's body waves, local sources), a layered geology; the speckle layer's "
        f"model is not validated for the real acquisition by the brightness map, which enters only the oracle.")


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'shafts':
        shaft_solves()
    else:
        main()
