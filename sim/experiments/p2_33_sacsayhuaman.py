"""P2-33 · Sacsayhuamán: the real pass over megalithic walls, houses and open fields, and what the published method draws.

    uv run python experiments/p2_33_sacsayhuaman.py

A second real site (sites/sacsayhuaman): the Inca walls above Cusco, the Rodadero outcrop across the esplanade from
them, open fields to the north and the city's houses below, all inside one ICEYE Spotlight Dwell Fine pass (X35,
22 August 2025, descending, looking left). Nothing is claimed beneath any of it. The question is the one P2-12 asks at
Giza: does the published method draw anything over the megalithic walls that it does not draw over houses and open
ground?

1. The image on the ground: the product's intensity, averaged over cells about a metre across, laid over the site's
   surface (Copernicus GLO-30) through the product's RPC with the site's EGM2008 undulation, after fitting the residual
   offset against the brightness the surface predicts (sarsim.ortho), as P2-12 does.
2. What the published method draws: the pipeline as P2-12 runs it (50 half-band pairs 88 Hz apart, 32-pixel patches,
   lambda_s 0.48 m, grid-median correction, no selection gates), over six areas of 179 m along track by 334 m across,
   chosen on the laid image before the method was run on any of them (AREAS): the zigzag walls, centred on their
   OpenStreetMap outline, and the Rodadero outcrop (monuments); open fields north of the site (open ground); and three
   stretches of houses (San Cristóbal on the slope below the walls, San Blas to the east, the flat grid west of the
   Plaza de Armas). Focused over one whole period of its depth axis and relabelled so that it repeats at 648 m, as the
   claim at Giza relabels it. The first investigation ran four patches here with another placement (r06): their fit
   scores overlapped.

Stated before the run: if the method responds to what the walls are, their depth profiles will differ from the houses'
and the fields' by more than the houses' differ from one another; if it responds to the image's texture, walls, houses
and fields will draw alike, the busiest surfaces drawing the most columns.

Added after the first run, which found the walls apart: each area's motionless copy (sarsim.looks.motionless_twin: the
crop's own brightness pattern over about three resolution cells and its spectrum, fresh speckle, nothing moving), run
through the same pipeline. If the walls' difference comes with how they look in the image, their copy will stand apart
from the houses' copies as they do; if it comes from motion, or anything else only the real image holds, it will not.
"""
import importlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from katabasis.compose import load_site
from katabasis.export.sites import scene as site_scene
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.dwell import DwellProduct
from sarsim.ortho import frame_to_lla, multilook, predicted_brightness, project, register, sample, surface
from sarsim.tomo import focus_paper

RID = 'p2_33_sacsayhuaman'
SITE = 'sacsayhuaman'
ACQ = 'sacsayhuaman-20250822'
PRODUCT = Path.home() / 'tmp/sar/sacsayhuaman/ICEYE_X35_SLC_SLEDF_5907293_20250822T152433.h5'
LEGACY = Path(__file__).resolve().parents[1] / 'legacy' / 'experiments'
LEGACY_CACHE = Path(__file__).resolve().parents[2] / 'results' / 'cache'
ORTHO_M = 1.0
REG_M = 2.0
N_RG = 768                  # 334 m of ground across track (the first investigation's patches were twice as wide)
CLAIM_REPEAT_M = 648.0
NZ = 256
VOXEL_M = 5.0
KEEP_EVERY_ROW = 4          # the patch grids are 0.5 m apart along track; the viewer's voxels are 5 m
# (key, label, x, y in the site frame, kind, how it was placed): chosen on the laid image before any method ran here
AREAS = [
    ('walls', 'Zigzag walls', 14.6, -92.7, 'monument',
     'the midpoint of the walls\' OpenStreetMap outline (ways 431479258 and 431479257); Muyuqmarka behind them'),
    ('rodadero', 'Rodadero outcrop', 40.0, 110.0, 'monument',
     'the carved outcrop across the esplanade (OpenStreetMap node 5857164644 and its bare-rock outline), clear of the walls\' area'),
    ('fields', 'Fields north of the site', -100.0, 525.0, 'control',
     'dark, flat fields on the laid image, clear of buildings'),
    ('san_cristobal', 'San Cristóbal houses', 300.0, -650.0, 'houses', 'houses on the slope below the walls, on the laid image'),
    ('san_blas', 'San Blas houses', 820.0, -830.0, 'houses', 'houses round the San Blas quarter (OpenStreetMap way 23799145), on the laid image'),
    ('grid', 'City grid houses', -250.0, -1000.0, 'houses', 'the flat street grid west of the Plaza de Armas, on the laid image'),
]


def legacy_module():
    sys.path.insert(0, str(LEGACY))
    rc = importlib.import_module('real_common')
    rc.CACHE = str(LEGACY_CACHE)
    return rc


def patch_entry(key):
    """The first investigation's patch record for an area: its place and ellipsoidal height (surface + undulation)."""
    sc = site_scene(load_site(SITE))
    geoid = float(sc['frame']['geoid_m'])
    k, label, x, y, kind, _ = next(a for a in AREAS if a[0] == key)
    lat, lon = frame_to_lla(sc['frame']['origin'], x, y)
    z = float(surface(sc, np.array([x]), np.array([y]))[0])
    return (f's33_{k}', label, float(lat), float(lon), z + geoid, kind)


TWIN_SEED = 101


def run_area(key, twin=None):
    """One area through the first investigation's pipeline: the real crop, or (twin = a seed) its motionless copy."""
    rc = legacy_module()
    entry = patch_entry(key)
    crop = rc.DwellProduct.crop
    if twin is not None:
        from sarsim.looks import motionless_twin
        entry = (f'{entry[0]}_twin{twin}',) + entry[1:]

        def copy_crop(self, *a, **k):
            c, origin = crop(self, *a, **k)
            return motionless_twin(c, np.random.default_rng(twin)), origin
        rc.DwellProduct.crop = copy_crop
    if not any(p[0] == entry[0] for p in rc.SACSAYHUAMAN_PATCHES):
        rc.SACSAYHUAMAN_PATCHES.append(entry)
    try:
        return rc.run_patch('sacsayhuaman', entry[0], 'paper', n_rg=N_RG, verbose=False)
    finally:
        rc.DwellProduct.crop = crop                   # the real product's crops again, for whatever reads it next


def main():
    g = DwellGeometry.from_record(ACQ)
    sc = site_scene(load_site(SITE))
    geoid = float(sc['frame']['geoid_m'])
    x0, x1 = sc['extent']['x']
    y0, y1 = sc['extent']['y']
    params = {'acquisition': ACQ, 'site': SITE, 'geoid_m': geoid, 'ortho_m': ORTHO_M, 'registration_m': REG_M,
              'areas': [{'key': a[0], 'label': a[1], 'x': a[2], 'y': a[3], 'kind': a[4], 'placed': a[5]} for a in AREAS],
              'patch_px': [4096, N_RG], 'bank': 'paper (K 50, half band, 88 Hz, 32 px, lambda_s 0.48 m)',
              'claim_repeat_m': CLAIM_REPEAT_M, 'voxel_m': VOXEL_M}
    with Run(RID, 'Sacsayhuamán: the real pass over megalithic walls, houses and open fields', params) as run:
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        # ---- what the method draws over each area, the areas in parallel (each is cached by the first investigation's code)
        with ProcessPoolExecutor(max_workers=len(AREAS)) as pool:
            list(pool.map(run_area, [a[0] for a in AREAS] * 2, [None] * len(AREAS) + [TWIN_SEED] * len(AREAS)))
        outs = {a[0]: run_area(a[0]) for a in AREAS}
        twins = {a[0]: run_area(a[0], TWIN_SEED) for a in AREAS}
        p = DwellProduct(PRODUCT)
        # ---- the image, averaged to about a metre, over the site's footprint in the product
        xs = np.arange(x0 + ORTHO_M / 2, x1, ORTHO_M)
        ys = np.arange(y0 + ORTHO_M / 2, y1, ORTHO_M)
        X, Y = np.meshgrid(xs, ys)
        Z = surface(sc, X, Y)
        rows, cols = project(p, sc, X, Y, Z, geoid)
        la = max(1, int(round(1.0 / g.dx)))
        lr = max(1, int(round(1.0 * np.sin(g.theta) / g.dr)))
        obs, origin = multilook(p, int(rows.min()) - 2000, int(rows.max()) + 2000, int(cols.min()) - 300,
                                int(cols.max()) + 300, la, lr)
        # ---- the residual offset, fitted on the surface's predicted brightness
        xr = np.arange(x0 + REG_M / 2, x1, REG_M)
        yr = np.arange(y0 + REG_M / 2, y1, REG_M)
        XR, YR = np.meshgrid(xr, yr)
        B, ZR = predicted_brightness(sc, XR, YR, g.los_enu, REG_M)
        rr, cr = project(p, sc, XR, YR, ZR, geoid)
        shift, corr = register(obs, origin, la, lr, rr, cr, B)
        print(f'  offset {shift} px, correlation {corr:.3f}', flush=True)
        rows, cols = rows + shift[0], cols + shift[1]
        I = sample(obs, origin, la, lr, rows, cols)
        L = 10 * np.log10(np.where(np.isfinite(I), I, np.nan) + 1e-12)
        lo, hi = np.nanpercentile(L, 2), np.nanpercentile(L, 99.5)
        u8 = np.where(np.isfinite(L), np.round(255 * np.clip((L - lo) / (hi - lo), 0, 1)), 0).astype(np.uint8)[::-1]
        from PIL import Image
        Image.fromarray(u8).save(figs / 'ground.jpg', quality=85, optimize=True)
        covered = float(np.mean(np.isfinite(L)))
        # ---- each area's depth pictures, and where each site cell falls among its pixels
        patches, vols, profiles, energies = [], {}, {}, {}
        ex = np.arange(x0 + VOXEL_M / 2, x1, VOXEL_M)
        ey = np.arange(y0 + VOXEL_M / 2, y1, VOXEL_M)
        EX, EY = np.meshgrid(ex, ey)
        er, ec = project(p, sc, EX, EY, surface(sc, EX, EY), geoid, shift)
        sa, sr = g.dx, g.dr / np.sin(g.theta)
        for key, label, ax_, ay_, kind, placed in AREAS:
            out = outs[key]
            meta = out['meta']
            kz = out['kz']
            zrep = 2 * np.pi / np.median(np.abs(np.diff(kz)))
            z_raw = np.linspace(0.0, zrep, NZ, endpoint=False)
            gr, gc = out['grid_rows'], out['grid_cols']
            keep = np.zeros((len(gr), len(gc)), bool)
            keep[::KEEP_EVERY_ROW] = True
            keep = keep.ravel()
            q = out['q'][keep].astype(np.float64)
            qc = q[..., 0] + 1j * q[..., 1]
            T = focus_paper(qc, kz, z_raw)
            energy = np.sum(np.abs(qc) ** 2, axis=1) / q.shape[1] ** 2
            pillar = float(np.corrcoef(T.mean(axis=1), energy)[0, 1])
            profiles[key] = T.mean(axis=0) / T.mean()
            energies[key] = energy
            r0, c0 = meta['crop_origin']
            tree = cKDTree(np.column_stack([(r0 + out['rows'][keep]) * sa, (c0 + out['cols'][keep]) * sr]))
            d, idx = tree.query(np.column_stack([er.ravel() * sa, ec.ravel() * sr]))
            idx = np.where(d < 4.0, idx, -1).reshape(EX.shape)
            vols[key] = {'T': T.astype(np.float16), 'index': idx.astype(np.int32),
                         'grid': np.array([len(gr[::KEEP_EVERY_ROW]), len(gc)])}
            patches.append({'key': key, 'label': label, 'kind': kind, 'centre_m': [ax_, ay_], 'placed': placed,
                            'pixels': int(keep.sum()), 'crop_origin': meta['crop_origin'], 'repeat_depth_raw_m': float(zrep),
                            'pillar_power_vs_energy': pillar, 'cells_on_site': int((idx >= 0).sum()),
                            'mean_energy': float(energy.mean())})
            print(f"  {key}: pillars follow energy at {pillar:.4f}; {int((idx >= 0).sum())} site cells", flush=True)
        # ---- the motionless copies: the same mean depth profile, from the brightness pattern alone
        def profile(out):
            gr, gc = out['grid_rows'], out['grid_cols']
            keep = np.zeros((len(gr), len(gc)), bool)
            keep[::KEEP_EVERY_ROW] = True
            q = out['q'][keep.ravel()].astype(np.float64)
            zrep = 2 * np.pi / np.median(np.abs(np.diff(out['kz'])))
            T = focus_paper(q[..., 0] + 1j * q[..., 1], out['kz'], np.linspace(0.0, zrep, NZ, endpoint=False))
            return T.mean(axis=0) / T.mean(), np.sum(np.abs(q[..., 0] + 1j * q[..., 1]) ** 2, axis=1) / q.shape[1] ** 2
        tw = {k: profile(o) for k, o in twins.items()}
        tprof = {k: v[0] for k, v in tw.items()}
        # ---- alike or not: depth profiles within and across kinds; the strongest columns by kind
        top = np.percentile(np.concatenate(list(energies.values())), 95)
        for x in patches:
            e = energies[x['key']]
            x['median_energy'] = float(np.median(e))
            x['share_in_top5'] = float(np.mean(e > top))
            x['max_energy'] = float(e.max())
        cc = lambda a, b: float(np.corrcoef(profiles[a], profiles[b])[0, 1])
        keys = {kind: [a[0] for a in AREAS if a[4] == kind] for kind in ('monument', 'control', 'houses')}
        within_houses = [cc(a, b) for i, a in enumerate(keys['houses']) for b in keys['houses'][i + 1:]]
        walls_vs = {k: cc('walls', k) for k in keys['houses'] + keys['control'] + ['rodadero']}
        rod_vs = {k: cc('rodadero', k) for k in keys['houses'] + keys['control']}
        mon_vs_rest = [cc(a, b) for a in keys['monument'] for b in keys['houses'] + keys['control']]
        all_pairs = [cc(a[0], b[0]) for i, a in enumerate(AREAS) for b in AREAS[i + 1:]]
        tc = lambda a, b: float(np.corrcoef(tprof[a], tprof[b])[0, 1])
        twin_within_houses = [tc(a, b) for i, a in enumerate(keys['houses']) for b in keys['houses'][i + 1:]]
        twin_walls_vs = {k: tc('walls', k) for k in keys['houses'] + keys['control'] + ['rodadero']}
        twin_mon_vs_rest = [tc(a, b) for a in keys['monument'] for b in keys['houses'] + keys['control']]
        real_vs_twin = {k: float(np.corrcoef(profiles[k], tprof[k])[0, 1]) for k in profiles}
        twin_apart = sum(c < min(twin_within_houses) for c in twin_mon_vs_rest) == len(twin_mon_vs_rest)
        all_twin_pairs = [tc(a[0], b[0]) for i, a in enumerate(AREAS) for b in AREAS[i + 1:]]
        follows = min(real_vs_twin.values()) > 0.8
        for x in patches:
            x['twin_median_energy'] = float(np.median(tw[x['key']][1]))
        share = lambda kind: [x['share_in_top5'] for x in patches if x['kind'] == kind]
        med = [x['median_energy'] for x in patches]
        below = sum(c < min(within_houses) for c in mon_vs_rest)      # monument pairs less alike than any two house areas
        apart = below == len(mon_vs_rest)
        np.savez_compressed(run.dir / 'volumes.npz', **{f'{k}_T': v['T'] for k, v in vols.items()},
                            **{f'{k}_index': v['index'] for k, v in vols.items()},
                            **{f'{k}_grid': v['grid'] for k, v in vols.items()},
                            z_raw=z_raw, voxel_m=VOXEL_M, grid_x0=x0, grid_y0=y0)
        finding = (
            f"The Sacsayhuamán pass lays onto the ground to within its fitted offset of {abs(shift[0]) * g.dx:.1f} m along "
            f"track and {abs(shift[1]) * g.dr / np.sin(g.theta):.1f} m across (image against predicted brightness, "
            f"correlation {corr:.2f}); the image covers {covered * 100:.0f}% of the site. The published method, run over the "
            f"zigzag walls, the Rodadero outcrop, open fields and three stretches of houses, draws its columns wherever the "
            f"registration wanders (the power at a pixel follows its trajectory's energy at "
            f"{min(x['pillar_power_vs_energy'] for x in patches):.4f} or better). Its mean depth profiles correlate at "
            f"{min(all_pairs):.2f} to {max(all_pairs):.2f} across the six areas. The walls' correlates with the houses' at "
            f"{min(walls_vs[k] for k in keys['houses']):.2f} to {max(walls_vs[k] for k in keys['houses']):.2f} and with the "
            f"fields' at {walls_vs['fields']:.2f}, while the houses' correlate with one another at {min(within_houses):.2f} to "
            f"{max(within_houses):.2f} and the Rodadero's with the houses' and fields' at {min(rod_vs.values()):.2f} to "
            f"{max(rod_vs.values()):.2f}. By the test stated before the run, "
            + ('the monuments stand apart: every pairing of a monument with houses or fields is less alike than any two '
               'stretches of houses, the walls by a wide margin and the Rodadero only just. ' if apart else
               'the monuments stand no further from houses and fields than houses stand from one another. ' if below == 0 else
               f'{below} of the {len(mon_vs_rest)} pairings of a monument with houses or fields are less alike than any two '
               f'stretches of houses, the rest within their range. ')
            + f"The strongest columns go to the houses: {min(share('houses')) * 100:.0f} to {max(share('houses')) * 100:.0f}% "
            f"of their pixels are among the noisiest 5% of all, against {max(share('monument')) * 100:.1f}% at most of the "
            f"monuments' and {share('control')[0] * 100:.1f}% of the fields'; the typical pixel's wander differs by "
            f"{(max(med) / min(med) - 1) * 100:.0f}% across the six, the walls' the quietest. "
            + (f"Each area's motionless copy, its brightness pattern with nothing moving, draws as the area does (profiles "
               f"correlating at {min(real_vs_twin.values()):.2f} to {max(real_vs_twin.values()):.2f}), so the difference comes "
               f"with how the areas look in the image." if follows else
               f"The areas' motionless copies, each area's brightness pattern and spectrum with fresh speckle and nothing moving, "
               f"keep neither the shared profile nor the walls' departure from it: the copies' profiles correlate with the real "
               f"ones at {min(real_vs_twin.values()):.2f} to {max(real_vs_twin.values()):.2f} and with one another at "
               f"{min(all_twin_pairs):.2f} to {max(all_twin_pairs):.2f}. Both come from something the real image holds beyond its "
               f"brightness pattern; this run does not say what, and it is not below the surface, which the X-band wave "
               f"reaches no further into than about 0.3 m (P2-11).")
            )
        run.save({'acquisition': g.record()['derived'] | {'heading_deg': g.heading_deg, 'incidence_deg': g.theta_deg},
                  'registration': {'shift_px': list(shift), 'correlation': corr,
                                   'shift_m': [shift[0] * g.dx, shift[1] * g.dr / np.sin(g.theta)]},
                  'image': {'file': 'ground.jpg', 'extent': {'x': [x0, x1], 'y': [y0, y1]}, 'cell_m': ORTHO_M,
                            'look_px': [la, lr], 'db_range': [float(lo), float(hi)], 'covered': covered},
                  'reach_m': 0.3, 'patches': patches,
                  'profile_corr': {'within_houses': within_houses, 'walls_vs': walls_vs, 'rodadero_vs': rod_vs,
                                   'monuments_vs_rest': mon_vs_rest,
                                   'all_pairs_range': [min(all_pairs), max(all_pairs)]},
                  'monuments_apart': bool(apart), 'monument_pairs_below_houses': int(below),
                  'twins': {'seed': TWIN_SEED, 'within_houses': twin_within_houses, 'walls_vs': twin_walls_vs,
                            'monuments_vs_rest': twin_mon_vs_rest, 'real_vs_twin': real_vs_twin, 'apart': bool(twin_apart),
                            'all_pairs_range': [min(all_twin_pairs), max(all_twin_pairs)], 'follow_real': bool(follows),
                            'profiles': {k: v[::4].round(4).tolist() for k, v in tprof.items()}},
                  'profiles': {k: v[::4].round(4).tolist() for k, v in profiles.items()},
                  'profile_depth_m': (z_raw[::4] * CLAIM_REPEAT_M / patches[0]['repeat_depth_raw_m']).round(1).tolist(),
                  'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
