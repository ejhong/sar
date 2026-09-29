"""P2-28 · Two depths, hollow and solid: does the published method's picture follow the room?

    uv run python experiments/p2_28_two_depths.py

A depth imager must put a room 30 m down deeper than a room 15 m down, and must tell a hollow room from a solid block of
the same size. P2-06 shows by algebra that the published focusing puts motion at a depth set by its frequency; P2-07 runs
the chain whole over the bench's room. Here the chain runs whole over three rooms under the same shaking:

- the bench's 6 m room, 15 m down (hollow);
- the same room 30 m down (hollow);
- a 6 m block of granite, 15 m down (solid, stiffer than the limestone round it).

Each room's imprint on the ground's motion is found as P2-04 finds it, by letting the lab's elastic solver settle, here
under a uniform horizontal strain rather than a stress, so that a solid inclusion is loaded as a passing wave loads it
(stress lambda(x) tr(eps) + 2 mu(x) eps in every cell; for a hollow room the same as P2-04's loading, which is checked).
The chain is P2-07's: simulated open desert, Giza's measured microseisms from every direction, the real dwell, and the
published method, with a lambda_s of 1.45 m so that its depth axis passes below both rooms (the published 0.48 m ends at
13.6 m; another value only relabels the same picture, P2-06).

The imprints are boosted, each to the same largest phase over the room (0.5 and 2 rad; P2-07's output changes visibly
from about 0.2 rad), and at one common gain, which keeps the rooms' own ratio. A null: the 15 m room's change to the
image at 2 rad with each pixel's phase at random.

Stated before the run: if the method images depth, the change the 30 m room makes will gather near 30 m and the 15 m
room's near 15 m, and the granite block's change will differ from the hollow room's more than the null's does; if the
depth is the shaking's frequency, the changes will spread alike over the whole axis, whatever the room.
"""
import copy
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from katabasis.ambient.field import microseisms, rayleigh_hv
from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run
from katabasis.seismic.arrays import surface_grid
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation
from sarsim import synthesize
from sarsim.acquisition import DwellGeometry

RID = 'p2_28_two_depths'
SITES = Path(__file__).resolve().parents[2] / 'sites'
ROOMS = {
    'hollow_15': {'centre': [0, 0, -15], 'size': [6, 6, 6], 'fill': 'air'},
    'hollow_30': {'centre': [0, 0, -30], 'size': [6, 6, 6], 'fill': 'air'},
    'granite_15': {'centre': [0, 0, -15], 'size': [6, 6, 6], 'fill': 'granite'},
}
H = 1.0
HALF_DOMAIN = 90.0
BOTTOM = -90.0
TOP = 4.0
HALF_REC = 60.0
SPACING = 2.0
TAPER = (50.0, 59.0)
STRAIN = 1e-4                       # the problem is linear; any small value will do
DURATION = 0.6
LAM_S = 1.45
PHASES = (0.5, 2.0)
NEAR_M = 15.0
CACHE = RESULTS / 'cache' / 'p2_28_kernels'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def media(room):
    site = load_site('bench-void')
    raw = copy.deepcopy(site.raw)
    if room is None:
        raw['features'] = []
    else:
        raw['features'][0]['shape'] = dict(raw['features'][0]['shape'], centre=room['centre'], size=room['size'])
        raw['features'][0]['fill'] = room['fill']
    g = Grid.covering((-HALF_DOMAIN, HALF_DOMAIN), (-HALF_DOMAIN, HALF_DOMAIN), (BOTTOM, TOP), H)
    return g, Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))


def relax_strain(medium, eps, receivers, host, damping=30.0, pml_width=12):
    """Displacement (n, 3; east, north, up) at `receivers` when the ground, strained uniformly by the horizontal `eps`
    {'xx', 'yy', 'xy'} (plane stress in the host rock: eps_zz = -lambda (eps_xx + eps_yy) / (lambda + 2 mu) of the host),
    comes to rest: every cell starts at its own stress lambda(x) tr(eps) + 2 mu(x) eps, so a hollow room (lambda = mu = 0)
    is loaded as P2-04 loads it and a solid inclusion is loaded as a passing wave loads it."""
    sim = Simulation(medium, pml_width=pml_width, f0=5.0)
    sxx, syy, szz, sxy, sxz, syz = sim.s
    lam_h, mu_h = host
    exx, eyy, exy = eps.get('xx', 0.0), eps.get('yy', 0.0), eps.get('xy', 0.0)
    ezz = -lam_h * (exx + eyy) / (lam_h + 2 * mu_h)
    lam, l2m = medium.lam.astype(np.float64), medium.lam2mu.astype(np.float64)
    tr = exx + eyy + ezz
    sxx[...] = (lam * tr + (l2m - lam) * exx).astype(np.float32)
    syy[...] = (lam * tr + (l2m - lam) * eyy).astype(np.float32)
    szz[...] = (lam * tr + (l2m - lam) * ezz).astype(np.float32)
    sxy[...] = (2 * medium.mu_xy.astype(np.float64) * exy).astype(np.float32)
    nt = int(math.ceil(DURATION / sim.dt))
    decay = np.float32(math.exp(-damping * sim.dt))

    def damp(it, s):
        for v in s.v:
            v *= decay
    res = sim.run([], Receivers(np.asarray(receivers, float)), nt, on_step=damp)
    u = np.cumsum(res.traces.astype(np.float64), axis=2) * sim.dt
    return u[:, :, -1], float(np.abs(u[:, :, -1] - u[:, :, -2]).max() / max(np.abs(u[:, :, -1]).max(), 1e-30))


def cached(name, compute):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f'{name}.npz'
    if path.is_file():
        K = np.load(path)
        return {k: K[k] for k in K.files if k != 'info'}, json.loads(str(K['info']))
    out, info = compute()
    np.savez_compressed(path, **out, info=json.dumps(info))
    return out, info


def settle(medium, rec, host, label):
    """The surface's displacement per unit strain for eps_ee, eps_nn and eps_en (the tensor component, as P2-04's Kxy)."""
    out, info = {'xy': rec[:, :2]}, {}
    for comp, key in (('xx', 'Kxx'), ('yy', 'Kyy'), ('xy', 'Kxy')):
        t0 = time.time()
        u, drift = relax_strain(medium, {comp: STRAIN}, rec, host)
        out[key] = u / STRAIN
        info[comp] = {'runtime_s': round(time.time() - t0, 1), 'last_step_change': drift}
        print(f"  {label} eps_{comp}: {info[comp]['runtime_s']:.0f} s", flush=True)
    return out, info


def room_kernels(name, room, host):
    """The room's imprint alone: the surface's displacement per unit strain with the room, less the same without it (the
    uniformly strained ground at rest, but for the absorbing edges); cached."""
    g, full = media(room)
    rec = surface_grid(g, full.solid, HALF_REC, SPACING)
    K1, i1 = cached(name, lambda: settle(full, rec, host, name))
    K0, i0 = cached('without', lambda: settle(media(None)[1], rec, host, 'without'))
    return {'xy': K1['xy'], **{k: K1[k] - K0[k] for k in ('Kxx', 'Kyy', 'Kxy')}}, {'room': i1, 'without': i0}


def interpolators(K, los):
    xy = K['xy']
    xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
    ix, iy = np.searchsorted(xs, xy[:, 0]), np.searchsorted(ys, xy[:, 1])
    out = []
    for name in ('Kxx', 'Kyy', 'Kxy'):
        M = np.zeros((len(xs), len(ys)))
        M[ix, iy] = K[name] @ los
        out.append(RegularGridInterpolator((xs, ys), M, bounds_error=False, fill_value=0.0))
    return out


def los_rms_map(K, los):
    """The imprint's line-of-sight displacement per unit strain, rms over the waves' directions (P2-04's map)."""
    az = np.deg2rad(np.arange(0, 180, 7.5))
    maps = [-(K['Kxx'] * np.sin(p) ** 2 + K['Kyy'] * np.cos(p) ** 2 + K['Kxy'] * np.sin(p) * np.cos(p)) @ los for p in az]
    return np.sqrt(np.mean(np.square(maps), axis=0))


def main():
    m07 = module('p207', Path(__file__).with_name('p2_07_whole_chain.py'))
    m07.TAPER = TAPER
    g = DwellGeometry.from_record('giza-20250827')
    amb = json.loads((SITES / 'ambient.json').read_text())
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    val = lambda m, q: m[q]['value'] if isinstance(m[q], dict) else m[q]
    rock = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    vp, vs, rho = val(rock, 'vp'), val(rock, 'vs'), val(rock, 'rho')
    host = (rho * vp ** 2 - 2 * rho * vs ** 2, rho * vs ** 2)
    los = np.asarray(g.los_enu)
    level = amb['regional']['microseism_vertical_0.1_0.3_hz']['value']
    speed = amb['regional']['microseism_phase_speed']['value']
    params = {'rooms': ROOMS, 'grid_m': H, 'domain_half_m': HALF_DOMAIN, 'domain_bottom_m': BOTTOM,
              'receivers': f'{2 * HALF_REC:.0f} m square at {SPACING:.0f} m', 'taper_m': TAPER, 'lam_s_m': LAM_S,
              'phases_rad': PHASES, 'geometry': 'giza-20250827', 'image_px': m07.SHAPE,
              'microseisms_m_s': level, 'band_hz': [0.1, 0.3], 'phase_speed_m_s': speed}
    with Run(RID, 'Two depths, hollow and solid: does the published picture follow the room?', params) as run:
        kern, check = {}, {}
        for name, room in ROOMS.items():
            K, info = room_kernels(name, room, host)
            rms = los_rms_map(K, los)
            kern[name] = {'K': K, 'rms': rms, 'info': info}
        # the hollow 15 m room under strain against P2-04's published map (loaded by stress on a smaller domain)
        p204 = json.loads((RESULTS / 'p2_04_chamber_imprint' / 'summary.json').read_text())
        pub = np.asarray(p204['imprint_map_los_per_strain_m'])
        xy = kern['hollow_15']['K']['xy']
        sel = (np.abs(xy[:, 0]) <= 40) & (np.abs(xy[:, 1]) <= 40)
        mine = kern['hollow_15']['rms'][sel]
        check['hollow_15_vs_p2_04'] = {'peak_ratio': float(mine.max() / pub.max()),
                                       'note': "P2-04's domain is 120 m wide and 64 m deep, this one 180 m and 94 m"}
        rooms = {}
        for name, k in kern.items():
            r = k['rms']
            X, Y = xy[:, 0], xy[:, 1]
            rooms[name] = {'peak_los_per_strain_m': float(r.max()),
                           'half_max_area_m2': float(np.sum(r >= r.max() / 2) * SPACING ** 2),
                           'half_max_radius_m': float(np.hypot(X, Y)[r >= r.max() / 2].max()),
                           'centre_uz_per_strain_m': float(k['K']['Kxx'][np.argmin(np.hypot(X, Y)), 2]
                                                           + k['K']['Kyy'][np.argmin(np.hypot(X, Y)), 2]),
                           'kernel_runs': k['info']}
        # the imprints themselves at the surface, signed (a wave along the look direction, as P2-07's map): how alike the
        # rooms' marks on the ground are before any method reads them
        phi = np.arctan2(g.ground_range_en[0], g.ground_range_en[1])
        sgn = {n: -(k['K']['Kxx'] * np.sin(phi) ** 2 + k['K']['Kyy'] * np.cos(phi) ** 2
                    + k['K']['Kxy'] * np.sin(phi) * np.cos(phi)) @ los for n, k in kern.items()}
        surface_corr = {'hollow_15_vs_hollow_30': float(np.corrcoef(sgn['hollow_15'], sgn['hollow_30'])[0, 1]),
                        'hollow_15_vs_granite_15': float(np.corrcoef(sgn['hollow_15'], sgn['granite_15'])[0, 1])}
        print(f"  rooms: {json.dumps({n: round(v['peak_los_per_strain_m'], 4) for n, v in rooms.items()})}; "
              f"surface maps {json.dumps({k: round(v, 3) for k, v in surface_corr.items()})}", flush=True)

        # the chain, as P2-07 runs it
        rng = np.random.default_rng(71)
        f = np.linspace(0.1, 0.3, 21)
        field = microseisms(f, np.full(len(f), level ** 2 / 0.2), np.random.default_rng(72), per_bin=4, speed=speed,
                            hv=rayleigh_hv(vp, vs))
        scat = m07.scene(g, rng)
        SH = m07.SHAPE
        rows = np.arange(48, SH[0] - 48, m07.STRIDE[0])
        cols = np.arange(32, SH[1] - 32, m07.STRIDE[1])
        RR, CC = np.meshgrid(rows, cols, indexing='ij')
        x_rows = (rows - SH[0] // 2) * g.dx
        foot = np.hypot(*np.meshgrid(x_rows, (cols - SH[1] // 2) * g.dr / np.sin(g.theta), indexing='ij')).ravel() < NEAR_M
        m07.LAM_S = LAM_S
        t_fine = np.linspace(-12.0, 12.0, 2401)
        interp = {n: interpolators(k['K'], los) for n, k in kern.items()}
        unit = {n: m07.Shaking(g, field, interp[n], 1.0).imprint_phase(t_fine) for n in ROOMS}
        t0 = time.time()
        img_B = synthesize(scat, g, SH, motion=m07.Shaking(g, field, interp['hollow_15'], 0.0), dtype=np.complex128)
        T_B, z = m07.run_method(img_B, g, RR.ravel(), CC.ravel())
        print(f"  without: {time.time() - t0:.0f} s; depth axis to {z[-1]:.1f} m", flush=True)
        cases, images, T = [], {}, {}
        for ph in PHASES:
            for n in ROOMS:
                cases.append((f'{n}_phase_{ph:g}', n, ph / unit[n]))
        common = PHASES[-1] / unit['hollow_15']
        for n in ('hollow_30', 'granite_15'):
            cases.append((f'{n}_common_gain', n, common))
        for label, n, G in cases:
            t0 = time.time()
            img = synthesize(scat, g, SH, motion=m07.Shaking(g, field, interp[n], G), dtype=np.complex128)
            images[label] = img
            T[label], _ = m07.run_method(img, g, RR.ravel(), CC.ravel())
            print(f"  {label}: gain {G:.3g}, {time.time() - t0:.0f} s", flush=True)
        d = images[f'hollow_15_phase_{PHASES[-1]:g}'] - img_B
        img_null = img_B + np.abs(d) * np.exp(2j * np.pi * np.random.default_rng(2028).random(d.shape))
        T['null'], _ = m07.run_method(img_null, g, RR.ravel(), CC.ravel())

        # where the change goes, over depth and over the ground
        in_band = lambda lo, hi: (z >= lo) & (z <= hi)
        results = {}
        for label in T:
            ch = T[label] - T_B
            a = np.abs(ch)
            prof = a[foot].mean(axis=0) / T_B.mean()
            inner = prof[3:-3]
            results[label] = {
                'over_room': float(a[foot].mean() / T_B[foot].mean()), 'elsewhere': float(a[~foot].mean() / T_B[~foot].mean()),
                'change_profile': prof, 'peak_depth_m': float(z[int(np.argmax(prof))]),
                'spread': [float(inner.min() / inner.mean()), float(inner.max() / inner.mean())],
                'share_12_18_m': float(prof[in_band(12, 18)].sum() / prof.sum()),
                'share_27_33_m': float(prof[in_band(27, 33)].sum() / prof.sum()),
                'axis_share_12_18_m': float(in_band(12, 18).mean()), 'axis_share_27_33_m': float(in_band(27, 33).mean())}
        corr = lambda a, b: float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
        pairs = {}
        for ph in PHASES:
            s = f'_phase_{ph:g}'
            c15, c30, cg = (T[f'{n}{s}'] - T_B for n in ROOMS)
            pairs[f'{ph:g}'] = {'profiles_15_vs_30': corr(results[f'hollow_15{s}']['change_profile'],
                                                          results[f'hollow_30{s}']['change_profile']),
                                'profiles_hollow_vs_granite': corr(results[f'hollow_15{s}']['change_profile'],
                                                                   results[f'granite_15{s}']['change_profile']),
                                'maps_15_vs_30': corr(c15[foot], c30[foot]),
                                'maps_hollow_vs_granite': corr(c15[foot], cg[foot])}
        pairs['null'] = {'maps_hollow_vs_null': corr((T[f'hollow_15_phase_{PHASES[-1]:g}'] - T_B)[foot], (T['null'] - T_B)[foot]),
                         'profiles_hollow_vs_null': corr(results[f'hollow_15_phase_{PHASES[-1]:g}']['change_profile'],
                                                         results['null']['change_profile'])}
        pairs['surface_imprints'] = surface_corr
        np.savez_compressed(Path(run.dir) / 'tomograms.npz', without=T_B.astype(np.float32), z=z, rows=rows, cols=cols,
                            **{k: v.astype(np.float32) for k, v in T.items()})

        hi = f'_phase_{PHASES[-1]:g}'
        top = pairs[f'{PHASES[-1]:g}']
        r15, r30, rg = (results[f'{n}{hi}'] for n in ROOMS)
        follows = r30['share_27_33_m'] > 2 * r30['axis_share_27_33_m'] and r15['share_12_18_m'] > 2 * r15['axis_share_12_18_m']
        finding = (
            ("The published method's picture does not follow the room's depth. " if not follows else
             "The published method's picture follows the room's depth. ")
            + f"Over simulated desert shaking at Giza's measured "
            f"microseism level, on the real dwell, with a depth axis from 0 to {z[-1]:.0f} m, the change each room's imprint "
            f"makes to the picture (each boosted to {PHASES[-1]:g} rad of phase over the room) spreads over the whole axis: "
            f"the 15 m room puts {100 * r15['share_12_18_m']:.0f}% of it between 12 and 18 m and the 30 m room "
            f"{100 * r30['share_27_33_m']:.0f}% between 27 and 33 m, where {100 * r15['axis_share_12_18_m']:.0f}% and "
            f"{100 * r30['axis_share_27_33_m']:.0f}% of the axis lie; their largest changes sit at {r15['peak_depth_m']:.1f} m "
            f"and {r30['peak_depth_m']:.1f} m, and their depth profiles correlate at {top['profiles_15_vs_30']:.2f}. "
            f"A granite block where the hollow room was changes the picture over the room with a map that correlates with the "
            f"hollow room's at {top['maps_hollow_vs_granite']:.2f} (the null, random phases at the same pixels: "
            f"{pairs['null']['maps_hollow_vs_null']:.2f}), as their marks on the ground correlate at "
            f"{surface_corr['hollow_15_vs_granite_15']:.2f} before any method reads them: the picture carries the surface's "
            f"motion, sign and all, spread over every depth (depth profiles correlating at "
            f"{top['profiles_hollow_vs_granite']:.2f}). "
            + ('Stated before the run, a depth imager would gather each room\'s change at its depth; this does. '
               if follows else 'Stated before the run, a depth imager would gather each room\'s change at its depth; '
               'neither room\'s does.')
            + f" The imprints came from the lab's solver settling under a uniform strain; for the hollow 15 m room its "
            f"peak is {check['hollow_15_vs_p2_04']['peak_ratio']:.3f} of P2-04's (a larger domain).")
        run.save({'rooms': rooms, 'check': check, 'surface_imprint_correlation': surface_corr, 'unit_gain_phase_rad': unit, 'depth_axis_m': z,
                  'results': results, 'pairs': pairs, 'follows_depth': bool(follows), 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
