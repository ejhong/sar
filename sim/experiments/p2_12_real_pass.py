"""P2-12 · The real pass over Giza, laid on the ground.

    uv run python experiments/p2_12_real_pass.py

What a viewer given a radar product can honestly show: where the satellite was, the image it made resampled onto the
ground, how far into the ground its wave reaches, and what the published method draws from it.

1. The pass: the acquisition record read from the product (P2-01).
2. The image on the ground: the product's intensity, averaged over cells about a metre across, resampled onto the Giza
   site's frame over its surface (terrain and pyramids), after fitting the residual geolocation offset by correlating the
   image with the brightness the surface predicts (sarsim.ortho). Layover stays in: a pyramid's top lands on the ground
   in front of it, as the radar records it.
3. The reach: X band gets about 0.3 m into the driest sand and less into rock (P2-11).
4. What the method draws: the published pipeline (50 half-band pairs 88 Hz apart, 32-pixel patches, lambda_s 0.48 m, as
   the first investigation ran it on this image) over the three pyramids and two stretches of open plateau, focused over
   one whole period of its depth axis and relabelled as the claim relabels it, so the axis repeats at 648 m (P2-08).
   Tested: is the picture over the monuments any different from the picture over empty desert?
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

from katabasis.compose import load_site
from katabasis.export.sites import scene as site_scene
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.dwell import DwellProduct
from sarsim.ortho import multilook, predicted_brightness, project, register, sample, surface
from sarsim.tomo import focus_paper

RID = 'p2_12_real_pass'
ACQ = 'giza-20250827'
PRODUCT = Path.home() / 'tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5'
LEGACY = Path(__file__).resolve().parents[1] / 'legacy' / 'experiments'
LEGACY_CACHE = Path(__file__).resolve().parents[2] / 'results' / 'cache'    # where the first investigation cached its runs
GEOID_M = 15.5              # ellipsoidal minus orthometric height at Giza (the first investigation's value); refitted below
ORTHO_M = 1.0
REG_M = 2.0
PATCHES = [('khafre', 'monument'), ('khufu', 'monument'), ('menkaure', 'monument'),
           ('desert_west', 'control'), ('desert_south', 'control')]
CLAIM_REPEAT_M = 648.0
NZ = 256
VOXEL_M = 5.0
KEEP_EVERY_ROW = 4          # the patch grids are 0.5 m apart in azimuth; the viewer's voxels are 5 m


def legacy_patch(key):
    sys.path.insert(0, str(LEGACY))
    rc = importlib.import_module('real_common')
    rc.CACHE = str(LEGACY_CACHE)
    return rc.run_patch('giza', key, 'paper')


def main():
    g = DwellGeometry.from_record(ACQ)
    sc = site_scene(load_site('giza'))
    x0, x1 = sc['extent']['x']
    y0, y1 = sc['extent']['y']
    params = {'acquisition': ACQ, 'site': 'giza', 'geoid_m': GEOID_M, 'ortho_m': ORTHO_M, 'registration_m': REG_M,
              'patches': PATCHES, 'bank': 'paper (K 50, half band, 88 Hz, 32 px, lambda_s 0.48 m)',
              'claim_repeat_m': CLAIM_REPEAT_M, 'voxel_m': VOXEL_M}
    with Run(RID, 'The real pass over Giza, laid on the ground', params) as run:
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        p = DwellProduct(PRODUCT)
        # ---- the image, averaged to about a metre, over the site's footprint in the product
        xs = np.arange(x0 + ORTHO_M / 2, x1, ORTHO_M)
        ys = np.arange(y0 + ORTHO_M / 2, y1, ORTHO_M)
        X, Y = np.meshgrid(xs, ys)
        Z = surface(sc, X, Y)
        rows, cols = project(p, sc, X, Y, Z, GEOID_M)
        la = max(1, int(round(1.0 / g.dx)))
        lr = max(1, int(round(1.0 * np.sin(g.theta) / g.dr)))
        pad = 2000
        obs, origin = multilook(p, int(rows.min()) - pad, int(rows.max()) + pad, int(cols.min()) - 300,
                                int(cols.max()) + 300, la, lr)
        print(f'  multilooked {obs.shape} cells of {la} x {lr} px', flush=True)
        # ---- the residual offset, fitted on the surface's predicted brightness
        xr = np.arange(x0 + REG_M / 2, x1, REG_M)
        yr = np.arange(y0 + REG_M / 2, y1, REG_M)
        XR, YR = np.meshgrid(xr, yr)
        B, ZR = predicted_brightness(sc, XR, YR, g.los_enu, REG_M)
        rr, cr = project(p, sc, XR, YR, ZR, GEOID_M)
        shift, corr = register(obs, origin, la, lr, rr, cr, B)
        print(f'  offset {shift} px, correlation {corr:.3f}', flush=True)
        # ---- the image on the ground
        rows, cols = rows + shift[0], cols + shift[1]
        I = sample(obs, origin, la, lr, rows, cols)
        L = 10 * np.log10(np.where(np.isfinite(I), I, np.nan) + 1e-12)
        lo, hi = np.nanpercentile(L, 2), np.nanpercentile(L, 99.5)
        u8 = np.round(255 * np.clip((L - lo) / (hi - lo), 0, 1))
        u8 = np.where(np.isfinite(L), u8, 0).astype(np.uint8)[::-1]          # north up
        from PIL import Image
        Image.fromarray(u8).save(figs / 'ground.jpg', quality=85, optimize=True)
        covered = float(np.mean(np.isfinite(L)))
        # ---- what the method draws, patch by patch
        patches, vols = [], {}
        profiles = {}
        for key, kind in PATCHES:
            out = legacy_patch(key)
            meta = out['meta']
            kz = out['kz']
            zrep = 2 * np.pi / np.median(np.abs(np.diff(kz)))
            z_raw = np.linspace(0.0, zrep, NZ, endpoint=False)
            gr, gc = out['grid_rows'], out['grid_cols']
            keep = np.zeros((len(gr), len(gc)), bool)
            keep[::KEEP_EVERY_ROW] = True
            keep = keep.ravel()
            q = out['q'][keep].astype(np.float64)
            T = focus_paper(q[..., 0] + 1j * q[..., 1], kz, z_raw)
            energy = np.sum(np.abs(q[..., 0] + 1j * q[..., 1]) ** 2, axis=1) / q.shape[1] ** 2
            pillar = float(np.corrcoef(T.mean(axis=1), energy)[0, 1])
            prof = T.mean(axis=0) / T.mean()
            profiles[key] = prof
            # where each site cell falls among the patch's pixels
            r0, c0 = meta['crop_origin']
            prow = r0 + out['rows'][keep]
            pcol = c0 + out['cols'][keep]
            ex = np.arange(x0 + VOXEL_M / 2, x1, VOXEL_M)
            ey = np.arange(y0 + VOXEL_M / 2, y1, VOXEL_M)
            EX, EY = np.meshgrid(ex, ey)
            er, ec = project(p, sc, EX, EY, surface(sc, EX, EY), GEOID_M, shift)
            # nearest patch pixel in image space (rows and columns scaled to metres on the ground)
            from scipy.spatial import cKDTree
            sa, sr = g.dx, g.dr / np.sin(g.theta)
            tree = cKDTree(np.column_stack([prow * sa, pcol * sr]))
            d, idx = tree.query(np.column_stack([er.ravel() * sa, ec.ravel() * sr]))
            idx = np.where(d < 4.0, idx, -1).reshape(EX.shape)
            vols[key] = {'T': T.astype(np.float16), 'index': idx.astype(np.int32)}
            patches.append({'key': key, 'label': meta['label'], 'kind': kind, 'pixels': int(keep.sum()),
                            'crop_origin': meta['crop_origin'], 'repeat_depth_raw_m': float(zrep),
                            'pillar_power_vs_energy': pillar, 'cells_on_site': int((idx >= 0).sum()),
                            'mean_energy': float(energy.mean())})
            print(f"  {key}: pillars follow energy at {pillar:.4f}; {int((idx >= 0).sum())} site cells", flush=True)
        mon = np.mean([profiles[k] for k, kind in PATCHES if kind == 'monument'], axis=0)
        con = np.mean([profiles[k] for k, kind in PATCHES if kind == 'control'], axis=0)
        same = float(np.corrcoef(mon, con)[0, 1])
        spread = [float(np.corrcoef(profiles[a], profiles[b])[0, 1]) for i, (a, _) in enumerate(PATCHES)
                  for b, _ in PATCHES[i + 1:]]
        np.savez_compressed(run.dir / 'volumes.npz', **{f'{k}_T': v['T'] for k, v in vols.items()},
                            **{f'{k}_index': v['index'] for k, v in vols.items()},
                            z_raw=z_raw, voxel_m=VOXEL_M, grid_x0=x0, grid_y0=y0)
        finding = (
            f"The real pass lays onto the ground to within its fitted offset of {abs(shift[0]) * g.dx:.1f} m along track "
            f"and {abs(shift[1]) * g.dr / np.sin(g.theta):.1f} m across (image against predicted brightness, correlation "
            f"{corr:.2f}); the image covers {covered * 100:.0f}% of the site. Its wave reaches at most about 0.3 m into dry "
            f"sand and less into rock. The published method, run on this image, draws the same kind of picture over the "
            f"three pyramids as over two stretches of empty plateau: the mean depth profiles over the monuments and over "
            f"the controls correlate at {same:.3f} (patch against patch {min(spread):.3f} to {max(spread):.3f}), and on every "
            f"patch the power at a pixel follows that pixel's trajectory energy ({min(x['pillar_power_vs_energy'] for x in patches):.4f} "
            f"or better), so its columns mark where the registration wanders, not what lies below.")
        run.save({'acquisition': g.record()['derived'] | {'heading_deg': g.heading_deg, 'incidence_deg': g.theta_deg},
                  'registration': {'shift_px': list(shift), 'correlation': corr,
                                   'shift_m': [shift[0] * g.dx, shift[1] * g.dr / np.sin(g.theta)]},
                  'image': {'file': 'ground.jpg', 'extent': {'x': [x0, x1], 'y': [y0, y1]}, 'cell_m': ORTHO_M,
                            'look_px': [la, lr], 'db_range': [float(lo), float(hi)], 'covered': covered},
                  'reach_m': 0.3, 'patches': patches, 'monument_vs_control_profile_corr': same,
                  'patch_profile_corr_range': [min(spread), max(spread)],
                  'profiles': {k: v[::4].round(4).tolist() for k, v in profiles.items()},
                  'profile_depth_m': (z_raw[::4] * CLAIM_REPEAT_M / patches[0]['repeat_depth_raw_m']).round(1).tolist(),
                  'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
