"""P2-16 · The gated reconstruction over the Great Pyramid, in three dimensions.

    uv run --with scikit-image python experiments/p2_16_khufu_volume.py

The reconstruction draws tomograms along lines: at every position of a line it fits a score at every nominal depth, so
one line is a curtain hanging below the surface. Its own eleven lines on Khufu cross the east face from east to west,
5 m apart. Here the same kind of line is laid across the whole pyramid, from 150 m east of its centre to 150 m west,
every 3 m from 141 m south to 141 m north: 95 lines of 101 positions, 3 m apart, so that the curtains tile a grid.
Each line is placed through the product's RPC over the site's surface (terrain and pyramid) with the EGM96
undulation, as the reconstruction's corrected geometry places its own lines. The 2022 image is run through the
unchanged pipeline (source and profile hashed) and, beside it, a motionless twin of the same crop (its brightness and
spectrum, fresh speckle, nothing moving and nothing inside). Pixel support is reported at one position, the loosest,
and at four, the profile's own. Only the pipeline's summary figure is skipped (95 panels); every array is kept.

The curtains are stacked into a volume on a 3 m grid: each voxel holds the pipeline's conditional adjusted R^2 at that
nominal depth below the surface point above it, unsmoothed, empty where nothing passed. They are what the pipeline
computes, nominal-depth fit scores, which its own manifest says are not detection probabilities, physical motion,
cavity identity, penetration or absolute depth. As a check a believer would make, the scores inside the surveyed
chambers (King's, Queen's, the Grand Gallery, the Subterranean Chamber, the passages, the ScanPyramids voids) are
compared with the scores at the same depths elsewhere, in the real image and in its twin.
"""
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter1d

import p2_13_gates as m
from katabasis.compose import load_site
from katabasis.compose.shapes import contains
from katabasis.export.sites import scene as site_scene
from katabasis.runs import Run
from sarsim.dwell import DwellProduct
from sarsim.looks import motionless_twin
from sarsim.ortho import frame_to_lla, project, surface

RID = 'p2_16_khufu_volume'
STEP_M = 3.0
HALF_EW_M = 150.0                        # each line runs from 150 m east of the centre to 150 m west
HALF_NS_M = 141.0                        # lines every 3 m from 141 m south to 141 m north
TRACE_POINTS = 10001
SUPPORTS = [1, 4]
CASES = ['real', 'twin101']


def box_mask(feature, X, Y, Z):
    """Voxels whose centres fall inside a surveyed box grown by one voxel, by the site's own inside test."""
    sh = dict(feature['shape'], size=[s + STEP_M for s in feature['shape']['size']])
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    return contains(sh, pts).reshape(X.shape)


def curtains_to_volume(F, z_m, z_surf, zc):
    """Each position's depth profile hung below its surface point, sampled at the voxel centres after taking the largest
    score within half a voxel either way (so thin peaks are not lost); NaN where the pipeline gave no score."""
    dz = float(z_m[1] - z_m[0])
    half = int(round(STEP_M / 2 / dz))
    filled = np.where(np.isfinite(F), F, -np.inf)
    peak = maximum_filter1d(filled, size=2 * half + 1, axis=2, mode='nearest')
    peak = np.where(np.isfinite(peak), peak, np.nan)
    depth = z_surf[:, :, None] - zc[None, None, :]                  # [lines, positions, voxel layers]
    idx = np.round((depth - z_m[0]) / dz).astype(int)
    valid = (idx >= 0) & (idx < len(z_m))
    vals = np.take_along_axis(peak, np.clip(idx, 0, len(z_m) - 1), axis=2)
    vals = np.where(valid, vals, np.nan)
    return vals[:, ::-1, :].transpose(1, 0, 2).astype(np.float32)   # (x west to east, y south to north, z down)


def main():
    sc = site_scene(load_site('giza'))
    p = DwellProduct(m.PRODUCT)
    kh = next(s for s in sc['structures'] if s['id'] == 'khufu')['shape']
    cx, cy = kh['centre'][0], kh['centre'][1]
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=SUPPORTS)
    ys = np.arange(-HALF_NS_M, HALF_NS_M + 1e-6, STEP_M)
    xs_trace = np.linspace(HALF_EW_M, -HALF_EW_M, TRACE_POINTS)          # P0 in the east, as the reconstruction's lines
    xs = xs_trace[np.linspace(0, TRACE_POINTS - 1, 101).round().astype(int)]
    params = {'acquisition': m.ACQ, 'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE),
              'khufu_centre_m': [cx, cy], 'lines': len(ys), 'positions_per_line': 101, 'step_m': STEP_M,
              'geoid_m': m.GEOID_M, 'supports': SUPPORTS, 'twin_seed': 101}
    with Run(RID, 'The gated reconstruction over the Great Pyramid, in three dimensions', params) as run:
        o = sc['frame']['origin']
        geom = {'title': "Lines across Khufu, east to west, every 3 m (the lab's raster, placed as the corrected geometry is)",
                'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
                'tracks': {}}
        heights = []
        z_surf = np.zeros((len(ys), 101))
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
        prof = dict(profile)
        prof['input'] = dict(profile['input'], hdf5_path=str(m.PRODUCT), geometry_path=str(gdir / 'geometry.json'))
        # the pipeline's own outputs are reused only when its inputs are byte for byte what this run would give it
        text_geom, text_prof = json.dumps(geom), json.dumps(prof, indent=1)
        same_inputs = ((gdir / 'geometry.json').exists() and (gdir / 'geometry.json').read_text() == text_geom
                       and (gdir / 'profile.json').exists() and (gdir / 'profile.json').read_text() == text_prof)
        (gdir / 'geometry.json').write_text(text_geom)
        (gdir / 'profile.json').write_text(text_prof)

        gated = m.load_gated()
        original = gated.load_source_and_preflight
        gated.render_tomograms = lambda *a, **k: None                   # the 95-panel summary figure only
        results, focus, reused = {}, {}, {}
        for case in CASES:
            def wrapped(config, config_path, case=case):
                out = list(original(config, config_path))
                if case.startswith('twin'):
                    out[7] = motionless_twin(out[7], np.random.default_rng(int(case[4:])))
                return tuple(out)
            gated.load_source_and_preflight = wrapped
            saved = gdir / case / 'biondi_v1_8_audit.npz'
            if same_inputs and saved.exists():
                made = json.loads((gdir / case / 'manifest.json').read_text())
                reused[case] = {'audit_sha256': m.sha256(saved), 'pipeline_run_created': made.get('created_utc'),
                                'geometry_sha256': made.get('geometry_sha256'), 'config_sha256': made.get('config_sha256')}
            else:
                gated.run_full(prof, gdir / 'profile.json', gdir / case)
            audit = np.load(saved)
            z_m = audit['z_m']
            focus[case] = {P: audit[f'focus_p{P}'] for P in SUPPORTS}              # [lines, 101, depths]
            results[case] = {
                'supported_positions': {P: int(np.isfinite(focus[case][P]).any(axis=2).sum()) for P in SUPPORTS},
                'first_gate_windows': int(sum(audit[f'{n.lower()}_window_pass'].sum() for n in geom['tracks'])),
                'exact_zero': float(np.mean(audit['Y'] == 0))}
            print(f"  {case}: supported positions " + ', '.join(f"P{P} {results[case]['supported_positions'][P]}" for P in SUPPORTS),
                  flush=True)
        gated.load_source_and_preflight = original

        # the volume: each column hangs below its surface point; voxels on a 3 m grid, the score the largest in the voxel
        z_top = float(z_surf.max())
        z_bottom = float(z_surf.min() - z_m[-1])
        zc = z_top - STEP_M * np.arange(int(np.ceil((z_top - z_bottom) / STEP_M)) + 1)
        vols, stats = {}, {}
        for case in CASES:
            for P in SUPPORTS:
                vols[f'{case}_p{P}'] = curtains_to_volume(focus[case][P], z_m, z_surf, zc)
        X, Y, Z = np.meshgrid(cx + xs[::-1], cy + ys, zc, indexing='ij')
        feats = [f for f in sc['features'] if f['id'].startswith('khufu')]
        inside = np.zeros(X.shape, bool)
        for f in feats:
            inside |= box_mask(f, X, Y, Z)
        for key, V in vols.items():
            ok = np.isfinite(V)
            layer = np.repeat(np.arange(V.shape[2])[None, None, :], V.shape[0] * V.shape[1], 0).reshape(V.shape)
            # the same depths elsewhere: every voxel in the layers the chambers occupy, outside them
            layers = np.unique(layer[inside])
            same = ok & ~inside & np.isin(layer, layers)
            stats[key] = {'voxels_with_score': int(ok.sum()),
                          'mean_score_in_chambers': float(np.nanmean(V[inside & ok])) if (inside & ok).any() else None,
                          'chamber_voxels_with_score': int((inside & ok).sum()), 'chamber_voxels': int(inside.sum()),
                          'mean_score_same_depths_elsewhere': float(np.nanmean(V[same])) if same.any() else None,
                          'share_above_0_5': float(np.mean(V[ok] >= 0.5)) if ok.any() else None}
        np.savez_compressed(Path(run.dir) / 'volumes.npz', **{k: v.astype(np.float16) for k, v in vols.items()},
                            x=cx + xs[::-1], y=cy + ys, z=zc, z_surface=z_surf[:, ::-1].T, step=STEP_M)
        r, t = results['real'], results['twin101']
        s1r, s1t = stats['real_p1'], stats['twin101_p1']
        fmt = lambda v: 'none' if v is None else f'{v:.3f}'
        finding = (
            f"Across the Great Pyramid, 95 east-west lines of 101 positions 3 m apart in the 2022 image, the unchanged gated "
            f"pipeline supports {r['supported_positions'][1]:,} positions at a support of one and {r['supported_positions'][4]} "
            f"at four; its motionless twin supports {t['supported_positions'][1]:,} and {t['supported_positions'][4]}. Stacked "
            f"into a volume, the real image's scores at support one fill {s1r['voxels_with_score']:,} voxels and the twin's "
            f"{s1t['voxels_with_score']:,}. Inside the surveyed chambers and passages the real image scores {fmt(s1r['mean_score_in_chambers'])} "
            f"on average against {fmt(s1r['mean_score_same_depths_elsewhere'])} at the same depths elsewhere (the twin: "
            f"{fmt(s1t['mean_score_in_chambers'])} and {fmt(s1t['mean_score_same_depths_elsewhere'])}), with "
            f"{s1r['chamber_voxels_with_score']} of {s1r['chamber_voxels']} chamber voxels holding any score.")
        run.save({'results': results, 'stats': stats, 'reused_pipeline_outputs': reused, 'grid': {'x0': float(cx + xs[-1]), 'y0': float(cy + ys[0]),
                  'z_top': z_top, 'step_m': STEP_M, 'shape': list(vols['real_p1'].shape)},
                  'chambers': [f['id'] for f in feats], 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
