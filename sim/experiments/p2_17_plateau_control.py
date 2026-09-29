"""P2-17 · The same raster over open plateau: an in-image control for P2-16.

    uv run --with scikit-image python experiments/p2_17_plateau_control.py

P2-16 laid 95 east-west lines of 101 positions across the Great Pyramid and stacked the gated reconstruction's curtains
into a volume, beside a motionless copy of the image. A synthetic copy is one kind of null; the other is the real image
itself where no monument stands. Here the identical raster (300 m by 282 m, lines and positions 3 m apart) is laid on
open plateau south-west of Menkaure, centred at (-470, -600) m in the site frame, chosen from the image as open sand
crossed by a track, clear of the pyramids and their queens' pyramids, before the method was run there. The lines follow
the terrain; everything else is P2-16's: the same unchanged pipeline, profile, supports (one position and four), RPC
placement with the EGM96 undulation, and the same stacking into a 3 m volume. Nothing is claimed under this ground,
though open plateau is not surveyed ground.
"""
import json
from pathlib import Path

import numpy as np

import p2_13_gates as m
import p2_16_khufu_volume as k
from katabasis.compose import load_site
from katabasis.export.sites import scene as site_scene
from katabasis.runs import Run
from sarsim.dwell import DwellProduct
from sarsim.ortho import frame_to_lla, project, surface

RID = 'p2_17_plateau_control'
CENTRE = (-470.0, -600.0)


def main():
    sc = site_scene(load_site('giza'))
    p = DwellProduct(m.PRODUCT)
    cx, cy = CENTRE
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=k.SUPPORTS)
    ys = np.arange(-k.HALF_NS_M, k.HALF_NS_M + 1e-6, k.STEP_M)
    xs_trace = np.linspace(k.HALF_EW_M, -k.HALF_EW_M, k.TRACE_POINTS)
    xs = xs_trace[np.linspace(0, k.TRACE_POINTS - 1, 101).round().astype(int)]
    params = {'acquisition': m.ACQ, 'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE),
              'centre_m': list(CENTRE), 'lines': len(ys), 'positions_per_line': 101, 'step_m': k.STEP_M,
              'geoid_m': m.GEOID_M, 'supports': k.SUPPORTS}
    with Run(RID, 'The same raster over open plateau: an in-image control for P2-16', params) as run:
        o = sc['frame']['origin']
        geom = {'title': 'The P2-16 raster laid on open plateau south-west of Menkaure',
                'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
                'tracks': {}}
        heights, z_surf = [], np.zeros((len(ys), 101))
        for j, y in enumerate(ys):
            X = cx + xs_trace
            Y = np.full_like(X, cy + y)
            Z = surface(sc, X, Y)
            r, c = project(p, sc, X, Y, Z, m.GEOID_M)
            heights.append(float(np.mean(Z)) + m.GEOID_M)
            z_surf[j] = surface(sc, cx + xs, np.full_like(xs, cy + y))
            lat, lon = frame_to_lla(o, X[[0, -1]], Y[[0, -1]])
            geom['tracks'][f'L{j:03d}'] = {
                'native_trace_col_row': np.stack([np.round(c), np.round(r)], 1).astype(int).tolist(),
                'geographic_endpoint_lat_lon': [[float(lat[0]), float(lon[0])], [float(lat[1]), float(lon[1])]]}
        geom['projection'] = {'height_m': float(np.mean(heights))}
        gdir = Path(__file__).resolve().parents[1] / 'data' / RID
        gdir.mkdir(parents=True, exist_ok=True)
        (gdir / 'geometry.json').write_text(json.dumps(geom))
        prof = dict(profile)
        prof['input'] = dict(profile['input'], hdf5_path=str(m.PRODUCT), geometry_path=str(gdir / 'geometry.json'))
        (gdir / 'profile.json').write_text(json.dumps(prof, indent=1))

        gated = m.load_gated()
        gated.render_tomograms = lambda *a, **kw: None
        gated.run_full(prof, gdir / 'profile.json', gdir / 'real')
        audit = np.load(gdir / 'real' / 'biondi_v1_8_audit.npz')
        z_m = audit['z_m']
        focus = {P: audit[f'focus_p{P}'] for P in k.SUPPORTS}
        supported = {P: int(np.isfinite(focus[P]).any(axis=2).sum()) for P in k.SUPPORTS}
        z_top = float(z_surf.max())
        z_bottom = float(z_surf.min() - z_m[-1])
        zc = z_top - k.STEP_M * np.arange(int(np.ceil((z_top - z_bottom) / k.STEP_M)) + 1)
        vols = {f'real_p{P}': k.curtains_to_volume(focus[P], z_m, z_surf, zc) for P in k.SUPPORTS}
        np.savez_compressed(Path(run.dir) / 'volumes.npz', **{kk: v.astype(np.float16) for kk, v in vols.items()},
                            x=cx + xs[::-1], y=cy + ys, z=zc, z_surface=z_surf[:, ::-1].T, step=k.STEP_M)
        stats = {kk: {'voxels_with_score': int(np.isfinite(v).sum()),
                      'share_above_0_5': float(np.mean(v[np.isfinite(v)] >= 0.5)) if np.isfinite(v).any() else None}
                 for kk, v in vols.items()}
        khufu = json.loads((Path(__file__).resolve().parents[1] / 'results' / 'p2_16_khufu_volume' / 'summary.json').read_text())
        kr = khufu['results']['real']['supported_positions']
        finding = (
            f"Laid on open plateau south-west of Menkaure, the P2-16 raster passes {supported[1]:,} positions at a support of "
            f"one and {supported[4]} at four in the real 2022 image, against {int(kr['1']):,} and {int(kr['4'])} across the "
            f"Great Pyramid; its volume at support one fills {stats['real_p1']['voxels_with_score']:,} voxels.")
        run.save({'supported_positions': supported, 'stats': stats, 'khufu_supported_positions': kr,
                  'grid': {'x0': float(cx + xs[-1]), 'y0': float(cy + ys[0]), 'z_top': z_top, 'step_m': k.STEP_M,
                           'shape': list(vols['real_p1'].shape)}, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
