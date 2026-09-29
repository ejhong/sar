"""P2-20 · The gated reconstruction over the known bench: with the chamber, without it, and with its imprint boosted.

    uv run --with scikit-image python experiments/p2_20_bench_gated.py

The lab's first bench is flat limestone with one 6 m chamber 15 m down. Its paper-style volumes (P2-07) are the same
with and without the chamber. Here the stricter gated reconstruction, run unchanged, gets the same test on synthetic
products it reads as it reads a real one: 8,192 rows of the 2022 X13 geometry, speckle over the whole scene, shaken as
Giza shakes both below the band the reconstruction reads and inside it: the measured microseisms (0.2 Hz) and the
measured 1-3 Hz level placed on the gates' second mode (2.14 Hz), each along the line of sight. The chamber's imprint is
P2-04's map of the surface's line-of-sight displacement per unit strain, tapered to zero over its outer 10 m so its
edge adds no step, times each band's strain, driven in phase with the ground's velocity: at the real level, and boosted
ten thousand and a hundred million times. A null product carries a random perturbation of the same size at the same
pixels as the real-level imprint. Motion enters through the Doppler-time mapping as a phase on each pixel, to second
order in the displacement (at the largest boost the phase stays under a quarter of a radian). The bench sits at the
scene centre by the product's own geometry, and 51 east-west lines of 51 positions, 2.4 m apart, cover it.

Stated before the run: if the reconstruction responds to the chamber, the positions whose fit scores change will gather
over it (within 12 m of its axis, where 3% of the positions lie) as the boost grows, with their largest changes at the
chamber's depth (12 to 18 m); if it does not, they will be as few and as scattered as the null's.
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
CASES = [('without', 0.0), ('with', 1.0), ('null', None), ('boost_1e4', 1e4), ('boost_1e8', 1e8)]
SUPPORTS = [1, 4]
TAPER_M = 10.0
NEAR_M = 12.0                         # "over the chamber": within twice its width of its axis
CHAMBER_DEPTH_M = (12.0, 18.0)
MODE = 2                              # the gates' mode the 1-3 Hz shaking is placed on


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


def tukey(n, spacing, width):
    """1 inside, falling to 0 at both ends over `width` metres as a raised cosine."""
    x = np.arange(n) * spacing
    edge = np.minimum(x, x[-1] - x)
    return np.where(edge < width, 0.5 * (1 - np.cos(np.pi * np.clip(edge / width, 0, 1))), 1.0)


def imprint_on_pixels(g, imprint, grid):
    """P2-04's line-of-sight displacement per unit strain, tapered at its edge, sampled on the image's pixels."""
    from scipy.ndimage import map_coordinates
    n, h, sp = grid['n'], grid['half_width_m'], grid['spacing_m']
    w = tukey(n, sp, TAPER_M)
    tapered = imprint * w[:, None] * w[None, :]
    rr, cc = np.meshgrid(np.arange(ROWS), np.arange(COLS), indexing='ij')
    e, nn = pixel_to_site(g, rr, cc)
    ix, iy = (e + h) / sp, (nn + h) / sp                   # the map's first axis is east, second north (surface_grid)
    return map_coordinates(tapered, [ix, iy], order=1, mode='constant', cval=0.0).astype(np.float64)


def move(slc, g, field, history):
    """Each pixel displaced along the line of sight by field * history(t): the phase exp(-4 pi j d / lambda) to second
    order in d, applied through the Doppler-time mapping (t = f / Ka) column by column."""
    f = frequencies(g, slc.shape[0])
    h = history(f / g.Ka_signed)[:, None]
    k = -4j * np.pi * field / g.lam
    spec = lambda x: np.fft.fftshift(np.fft.fft(x, axis=0), axes=0)
    back = lambda X: np.fft.ifft(np.fft.ifftshift(X, axes=0), axis=0)
    s = slc.astype(np.complex128)
    return s + back(spec(s * k) * h) + back(spec(s * k * k / 2) * h * h)


def main():
    g = DwellGeometry.from_record(m.ACQ)
    st = EchoSetup.from_geometry(g)
    lat0, lon0 = g.source['scene_centre']['latitude'], g.source['scene_centre']['longitude']
    imp = json.loads((Path(__file__).resolve().parents[1] / 'results' / 'p2_04_chamber_imprint' / 'summary.json').read_text())
    micro = next(c for c in imp['cases'] if c['case'].startswith('Giza microseisms'))
    band = next(c for c in imp['cases'] if c['case'].startswith('Giza 1-3 Hz'))
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=SUPPORTS)
    bank, gates = profile['processing']['filter_bank'], profile['gates']
    window_s = gates['window_pairs'] * bank['k_leap_hz'] / abs(g.Ka_signed)
    f_mode = MODE / window_s
    # each band: frequency, the strain its measured motion carries (hv V / c, P2-04), and its vertical velocity
    bands = [(micro['frequency_hz'], micro['strain'], micro['vertical_velocity_m_s']),
             (f_mode, band['strain'], band['vertical_velocity_m_s'])]
    params = {'acquisition': m.ACQ, 'shape': [ROWS, COLS], 'step_m': STEP, 'half_m': HALF, 'cases': CASES,
              'bands': [{'frequency_hz': f, 'strain': e, 'vertical_velocity_m_s': v} for f, e, v in bands],
              'band_sources': [micro['source'], band['source'] + f'; placed on the gates\' mode {MODE}'],
              'taper_m': TAPER_M, 'near_m': NEAR_M, 'supports': SUPPORTS,
              'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE)}
    with Run(RID, 'The gated reconstruction over the known bench', params) as run:
        gdir = Path(__file__).resolve().parents[1] / 'data' / RID
        gdir.mkdir(parents=True, exist_ok=True)
        speckle = r15.band_noise(g, (ROWS, COLS), np.random.default_rng(20)).astype(np.complex128)
        field = imprint_on_pixels(g, np.array(imp['imprint_map_los_per_strain_m']), imp['kernel_grid'])
        # the ground's own motion, common to the scene: line-of-sight displacement from each band's vertical velocity
        cos_t = abs(g.los_enu[2])
        ambient = lambda t: sum(v * cos_t / (2 * np.pi * f) * np.sin(2 * np.pi * f * t) for f, _, v in bands)
        # the chamber's imprint follows the strain, which moves with the ground's velocity
        strain = lambda t: sum(e * np.cos(2 * np.pi * f * t) for f, e, _ in bands)
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
        focus, results, images, z_m = {}, {}, {}, None
        for name, boost in CASES:
            if boost == 0:
                img = base
            elif boost is None:
                # the null: the real-level imprint's change, each pixel's phase replaced at random
                d = images['with'] - base
                img = base + np.abs(d) * np.exp(2j * np.pi * np.random.default_rng(2020).random(d.shape))
            else:
                img = move(base, g, field, lambda t, b=boost: b * strain(t))
            img32 = img.astype(np.complex64)
            images[name] = img
            path = write_synthetic_slc(gdir / f'SYNTHETIC_bench_{name}.h5', img32, st, g.Ka_signed, lat0, lon0,
                                       75.0, g.heading_deg, g.theta_deg,
                                       scene={'bench': 'bench-void', 'chamber_imprint_boost': boost, 'bands': params['bands'],
                                              'null': boost is None, 'speckle_seed': 20})
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
            b32 = images['without'].astype(np.complex64)
            results[name] = {'boost': boost, 'supported_positions': {P: int(np.isfinite(focus[name][P]).any(axis=2).sum()) for P in SUPPORTS},
                             'values_changed': int(((img32.real != b32.real) | (img32.imag != b32.imag)).sum()),
                             'image_change_max_rel': float(np.abs(img - base).max() / np.abs(base).max())}
            print(f"  {name}: {results[name]}", flush=True)

        # where the fit scores changed, against the product without the chamber
        X, Y = np.meshgrid(xs, ys)                                   # [line, position]: east, north
        near = np.hypot(X, Y) <= NEAR_M
        share = float(near.mean())
        in_band = (z_m >= CHAMBER_DEPTH_M[0]) & (z_m <= CHAMBER_DEPTH_M[1])
        for name, _ in CASES[1:]:
            for P in SUPPORTS:
                A, B = focus[name][P], focus['without'][P]
                diff = np.abs(np.nan_to_num(A) - np.nan_to_num(B))
                changed = diff.max(axis=2) > 1e-6
                flipped = np.isfinite(A).any(axis=2) != np.isfinite(B).any(axis=2)
                deepest = z_m[np.argmax(diff, axis=2)]
                results[name][f'p{P}'] = {
                    'identical': bool(not changed.any()),
                    'positions_changed': int(changed.sum()), 'positions_changed_near': int((changed & near).sum()),
                    'passes_flipped': int(flipped.sum()), 'passes_flipped_near': int((flipped & near).sum()),
                    'max_score_change': float(diff.max()),
                    'near_changes_at_chamber_depth': int((changed & near & in_band[np.argmax(diff, axis=2)]).sum()),
                    'depths_of_largest_change_near_m': [float(z) for z in deepest[changed & near]]}
        z_surf = np.zeros((len(ys), n_pos))
        zc = -STEP * np.arange(int(np.ceil(z_m[-1] / STEP)) + 1)
        vols = {f'{name}_p{P}': curtains_to_volume(focus[name][P], z_m, z_surf, zc, STEP) for name, _ in CASES for P in SUPPORTS}
        np.savez_compressed(Path(run.dir) / 'volumes.npz', **{k: v.astype(np.float16) for k, v in vols.items()},
                            x=xs[::-1], y=ys, z=zc, z_surface=z_surf[:, ::-1].T, step=STEP)

        # the finding, from the numbers
        w, nl, b4, b8 = (results[k]['p1'] for k in ('with', 'null', 'boost_1e4', 'boost_1e8'))
        n_all = len(ys) * n_pos
        v8 = 1e8 * float(np.abs(field).max()) * bands[1][1] * 2 * np.pi * f_mode
        gathers = lambda r: r['positions_changed_near'] >= 5 and r['positions_changed_near'] > 3 * share * r['positions_changed']
        where = lambda r: f"{r['positions_changed']:,} ({r['positions_changed_near']} within {NEAR_M:.0f} m of the chamber)"
        finding = (
            f"Over the one-chamber bench, on synthetic products the unchanged gated pipeline reads as it reads a real one, "
            f"shaken as Giza shakes ({bands[0][0]:.1f} Hz, and the measured 1-3 Hz level on its mode {MODE}, {f_mode:.2f} Hz), "
            f"{results['without']['supported_positions'][1]:,} of {n_all:,} positions pass at a support of one without the "
            f"chamber and {results['with']['supported_positions'][1]:,} with it. The chamber's real imprint changes "
            f"{results['with']['values_changed']:,} of the {2 * ROWS * COLS / 1e6:.1f} million values the pipeline reads and the "
            f"fit scores at {where(w)} of the positions, where {share * 100:.0f}% of positions lie so near; a random "
            f"perturbation of the same size at the same pixels changes {where(nl)}. Boosted ten thousand times the imprint "
            f"changes {where(b4)}; a hundred million times, {v8 * 1e3:.1f} mm/s of line-of-sight motion over the chamber at "
            f"{f_mode:.2f} Hz, {where(b8)}"
            + (f", {b8['near_changes_at_chamber_depth']} of them largest at the chamber's depth" if b8['positions_changed_near'] else '')
            + '. '
            + ('The changes gather over the chamber only at the largest boost.' if gathers(b8) and not gathers(b4) and not gathers(w)
               else 'The changes gather over the chamber from a boost of ten thousand.' if gathers(b4) and not gathers(w)
               else 'The changes gather over the chamber at its real level.' if gathers(w)
               else 'At no boost do the changes gather over the chamber: what the method draws over the bench does not depend on it.'))
        run.save({'results': results, 'near_share': share, 'mode_hz': f_mode, 'imprint_los_velocity_boost_1e8_m_s': v8,
                  'grid': {'x0': float(xs[-1]), 'y0': float(ys[0]), 'z_top': float(zc[0]), 'step_m': STEP,
                           'shape': list(vols['without_p1'].shape)}, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
