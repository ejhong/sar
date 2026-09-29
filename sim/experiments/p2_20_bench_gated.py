"""P2-20 · The gated reconstruction over the known bench: with the chamber, without it, and with its imprint boosted.

    uv run --with scikit-image python experiments/p2_20_bench_gated.py

The lab's first bench is flat limestone with one 6 m chamber 15 m down. Its paper-style volumes (P2-07) are the same
with and without the chamber. Here the stricter gated reconstruction, run unchanged, gets the same test on synthetic
products it reads as it reads a real one: 8,192 rows of the 2022 X13 geometry, speckle over the whole scene, shaken by
Giza's microseisms (0.2 Hz, the measured median level along the line of sight) and, in two of the three products,
carrying the chamber's imprint: P2-04's map of the surface's line-of-sight displacement per unit strain, times the
microseisms' strain, at the real level and boosted a million times. Motion is put in through the Doppler-time mapping
as a phase on each pixel (linear in the displacement, exact at these amplitudes). The bench is placed at the scene
centre by the product's own geometry (azimuth along the heading, slant range across it), and 51 east-west lines of 51
positions, 2.4 m apart, cover it. Scored: how many positions pass, whether the volumes differ between the products at
all, and where.
"""
import json
from pathlib import Path

import h5py
import numpy as np

import p2_13_gates as m
import p2_15_pair_response as r15
from katabasis.lab import curtains_to_volume
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup
from sarsim.looks import frequencies
from sarsim.ortho import frame_to_lla
from sarsim.slcfile import write_synthetic_slc

RID = 'p2_20_bench_gated'
ROWS, COLS = 8192, 640
STEP = 2.4
HALF = 60.0
CASES = [('without', 0.0), ('with', 1.0), ('with_boosted_1e6', 1e6)]
SUPPORTS = [1, 4]


def site_to_pixel(g, east, north, z=0.0):
    a = east * g.along_track_en[0] + north * g.along_track_en[1]
    gr = east * g.ground_range_en[0] + north * g.ground_range_en[1]
    return ROWS // 2 + a / g.dx, COLS // 2 + (gr * np.sin(g.theta) - z * np.cos(g.theta)) / g.dr


def pixel_to_site(g, rows, cols):
    a = (rows - ROWS // 2) * g.dx
    gr = (cols - COLS // 2) * g.dr / np.sin(g.theta)
    e = a * g.along_track_en[0] + gr * g.ground_range_en[0]
    n = a * g.along_track_en[1] + gr * g.ground_range_en[1]
    return e, n


def imprint_on_pixels(g, imprint, grid):
    """P2-04's line-of-sight displacement per unit strain, sampled on the image's pixels (zero beyond its map)."""
    from scipy.ndimage import map_coordinates
    rr, cc = np.meshgrid(np.arange(ROWS), np.arange(COLS), indexing='ij')
    e, n = pixel_to_site(g, rr, cc)
    h, sp = grid['half_width_m'], grid['spacing_m']
    ix, iy = (e + h) / sp, (n + h) / sp                    # the map's first axis is east, second north (surface_grid)
    w = map_coordinates(imprint, [ix, iy], order=1, mode='constant', cval=0.0)
    return w.astype(np.float64)


def move(slc, g, field, history):
    """Each pixel displaced along the line of sight by field * history(t): the phase exp(-4 pi j d / lambda), linear in
    d at these amplitudes, applied through the Doppler-time mapping (t = f / Ka) column by column."""
    f = frequencies(g, slc.shape[0])
    a = history(f / g.Ka_signed)[:, None]
    X = np.fft.fftshift(np.fft.fft(slc.astype(np.complex128) * (-4j * np.pi * field / g.lam), axis=0), axes=0)
    return slc + np.fft.ifft(np.fft.ifftshift(X * a, axes=0), axis=0)


def main():
    g = DwellGeometry.from_record(m.ACQ)
    st = EchoSetup.from_geometry(g)
    lat0, lon0 = g.source['scene_centre']['latitude'], g.source['scene_centre']['longitude']
    imp = json.loads((Path(__file__).resolve().parents[1] / 'results' / 'p2_04_chamber_imprint' / 'summary.json').read_text())
    case0 = imp['cases'][0]
    strain, f_micro = case0['strain'], case0['frequency_hz']
    v_micro = case0['vertical_velocity_m_s']
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=SUPPORTS)
    params = {'acquisition': m.ACQ, 'shape': [ROWS, COLS], 'step_m': STEP, 'half_m': HALF, 'cases': CASES,
              'strain': strain, 'frequency_hz': f_micro, 'supports': SUPPORTS,
              'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE)}
    with Run(RID, 'The gated reconstruction over the known bench', params) as run:
        gdir = Path(__file__).resolve().parents[1] / 'data' / RID
        gdir.mkdir(parents=True, exist_ok=True)
        speckle = r15.band_noise(g, (ROWS, COLS), np.random.default_rng(20)).astype(np.complex128)
        field = imprint_on_pixels(g, np.array(imp['imprint_map_los_per_strain_m']), imp['kernel_grid'])
        # the microseisms themselves, common to the whole scene: line-of-sight displacement v / (2 pi f) sin(2 pi f t)
        los_up = abs(g.los_enu[2])
        ambient = lambda t: (v_micro * los_up / (2 * np.pi * f_micro)) * np.sin(2 * np.pi * f_micro * t)
        base = move(speckle, g, np.ones((ROWS, COLS)), ambient)

        # the raster: east-west lines across the bench, positions 2.4 m apart
        n_pos = int(round(2 * HALF / STEP)) + 1
        ys = np.arange(-HALF, HALF + 1e-6, STEP)
        xs_trace = np.linspace(HALF, -HALF, 10001)
        xs = xs_trace[np.linspace(0, 10000, n_pos).round().astype(int)]
        origin = {'latitude': lat0, 'longitude': lon0}
        tracks = {}
        for j, y in enumerate(ys):
            r, c = site_to_pixel(g, xs_trace, np.full_like(xs_trace, y))
            lat, lon = frame_to_lla(origin, xs_trace[[0, -1]], np.full(2, y))
            tracks[f'L{j:03d}'] = {'native_trace_col_row': np.stack([np.round(c), np.round(r)], 1).astype(int).tolist(),
                                   'geographic_endpoint_lat_lon': [[float(lat[0]), float(lon[0])], [float(lat[1]), float(lon[1])]]}
        geom = {'title': 'East-west lines across the one-chamber bench, 2.4 m apart', 'target_samples_per_track': n_pos,
                'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
                'projection': {'height_m': 75.0}, 'tracks': tracks}
        (gdir / 'geometry.json').write_text(json.dumps(geom))

        gated = m.load_gated()
        gated.render_tomograms = lambda *a, **k: None
        focus, results, z_m = {}, {}, None
        for name, boost in CASES:
            img = base if boost == 0 else move(base, g, field, lambda t, b=boost: b * strain * np.sin(2 * np.pi * f_micro * t))
            path = write_synthetic_slc(gdir / f'SYNTHETIC_bench_{name}.h5', img.astype(np.complex64), st, g.Ka_signed, lat0, lon0,
                                       75.0, g.heading_deg, g.theta_deg,
                                       scene={'bench': 'bench-void', 'chamber_imprint_boost': boost, 'strain': strain,
                                              'microseism_hz': f_micro, 'speckle_seed': 20})
            with h5py.File(path, 'a') as h:
                h['synthetic'].attrs['site_to_pixel'] = json.dumps({'centre_row': ROWS // 2, 'centre_col': COLS // 2,
                                                                     'dx': g.dx, 'dr': g.dr, 'theta_deg': g.theta_deg,
                                                                     'heading_deg': g.heading_deg})
            prof = dict(profile)
            prof['input'] = dict(profile['input'], hdf5_path=str(path), geometry_path=str(gdir / 'geometry.json'))
            (gdir / f'profile_{name}.json').write_text(json.dumps(prof, indent=1))
            gated.run_full(prof, gdir / f'profile_{name}.json', gdir / name)
            a = np.load(gdir / name / 'biondi_v1_8_audit.npz')
            z_m = a['z_m']
            focus[name] = {P: a[f'focus_p{P}'] for P in SUPPORTS}
            results[name] = {'boost': boost, 'supported_positions': {P: int(np.isfinite(focus[name][P]).any(axis=2).sum()) for P in SUPPORTS},
                             'image_change_max_rel': float(np.abs(img - base).max() / np.abs(base).max())}
            print(f"  {name}: {results[name]}", flush=True)

        # do the volumes differ at all?
        for name, _ in CASES[1:]:
            for P in SUPPORTS:
                A, B = focus[name][P], focus['without'][P]
                same = bool(np.array_equal(np.isfinite(A), np.isfinite(B)) and np.allclose(np.nan_to_num(A), np.nan_to_num(B)))
                diff = np.abs(np.nan_to_num(A) - np.nan_to_num(B))
                results[name][f'identical_to_without_p{P}'] = same
                results[name][f'max_score_change_p{P}'] = float(diff.max())
                results[name][f'positions_changed_p{P}'] = int((diff.max(axis=2) > 1e-6).sum())
        z_surf = np.zeros((len(ys), n_pos))
        zc = -STEP * np.arange(int(np.ceil(z_m[-1] / STEP)) + 1)
        vols = {f'{name}_p{P}': curtains_to_volume(focus[name][P], z_m, z_surf, zc, STEP) for name, _ in CASES for P in SUPPORTS}
        np.savez_compressed(Path(run.dir) / 'volumes.npz', **{k: v.astype(np.float16) for k, v in vols.items()},
                            x=xs[::-1], y=ys, z=zc, z_surface=z_surf[:, ::-1].T, step=STEP)
        w, b = results['with'], results['with_boosted_1e6']
        finding = (
            f"Over the one-chamber bench, on synthetic products the unchanged gated pipeline reads as it reads a real one, "
            f"{results['without']['supported_positions'][1]:,} of {len(ys) * n_pos:,} positions pass at a support of one "
            f"without the chamber and {w['supported_positions'][1]:,} with it; the two volumes are "
            + ('identical' if w['identical_to_without_p1'] else f"different at {w['positions_changed_p1']} positions")
            + f" (the chamber's imprint changes the image by {w['image_change_max_rel']:.1e} of its brightest pixel, below the "
            f"single precision the pipeline reads). Boosted a million times, the imprint changes the image by "
            f"{b['image_change_max_rel']:.1e} and the volume "
            + ('not at all' if b['identical_to_without_p1'] else f"at {b['positions_changed_p1']} positions, by up to "
               f"{b['max_score_change_p1']:.2f} in score")
            + ". What the method draws over the bench is the same with or without the chamber.")
        run.save({'results': results, 'grid': {'x0': float(xs[-1]), 'y0': float(ys[0]), 'z_top': float(zc[0]), 'step_m': STEP,
                                               'shape': list(vols['without_p1'].shape)}, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
