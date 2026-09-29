"""P2-24 · The gated reconstruction on the one shaking: the field test's products, with the chamber and without.

    uv run --with scikit-image python experiments/p2_24_gated_on_the_shaking.py

P2-23 finds the force at which corner reflectors laid across the bench's chamber let the favourable detector tell the
chamber from none. Here the gated reconstruction, run unchanged, reads products of the same scene: speckle, the real
image's natural bright points, and the row of reflectors, all moving with P2-22's vibrator at 10.69 Hz (the method's
own tenth mode). Two forces: P2-23's boundary for the reflectors, where the data demonstrably hold the chamber, and a
thousandth of it. At each, a product with the chamber, one without, and a null: the product without the chamber,
changed by random noise the size of the chamber's own change. 51 east-west lines of 51 positions, 2.4 m apart, as in
P2-20.

The motion enters exactly: each pixel's line-of-sight displacement Re(D exp(i w t)) is expanded in Bessel terms
(Jacobi-Anger), each separable in place and time, and applied through the Doppler-time mapping. Within 10 m of the
vibrator the displacement is held at its value 10 m out, where a real ground would not stay linear at these forces.

Stated before the run: if the reconstruction reads the chamber where the data hold it, the positions whose scores change
between the products with and without the chamber will gather over it (within 12 m of its axis) beyond the null's; if
not, what changes will be as scattered as the null's.
"""
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import jv

import p2_13_gates as m
import p2_20_bench_gated as b20
import p2_23_every_reader as e23
from katabasis.lab import curtains_to_volume
from katabasis.runs import RESULTS, Run
from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup
from sarsim.ortho import frame_to_lla
from sarsim.slcfile import write_synthetic_slc

RID = 'p2_24_gated_on_the_shaking'
ROWS, COLS = b20.ROWS, b20.COLS
STEP, HALF = b20.STEP, b20.HALF
SUPPORTS = [1, 4]
NEAR_M = 12.0
SOURCE_CLEAR_M = 10.0
FORCE_FRACTIONS = [1.0, 1e-3]


def los_displacement(g, name, force, f):
    """The line-of-sight complex displacement amplitude D at every pixel (site position from the product's geometry)."""
    rr, cc = np.meshgrid(np.arange(ROWS), np.arange(COLS), indexing='ij')
    e, n = b20.pixel_to_site(g, rr, cc)
    d = np.load(RESULTS / 'p2_22_one_shaking' / 'fields.npz')
    x, y, src = d['x'], d['y'], d['source']
    los = np.asarray(g.los_enu, float)
    V = sum(los[k] * d[f'U_{name}'][..., k] for k in range(3))
    re = RegularGridInterpolator((x, y), V.real, bounds_error=False, fill_value=0.0)
    im = RegularGridInterpolator((x, y), V.imag, bounds_error=False, fill_value=0.0)
    # held at its value 10 m from the vibrator inside that radius
    r = np.hypot(e - src[0], n - src[1])
    scale = np.where(r < SOURCE_CLEAR_M, SOURCE_CLEAR_M / np.maximum(r, 1e-6), 1.0)
    ex, ny = src[0] + (e - src[0]) * scale, src[1] + (n - src[1]) * scale
    pts = np.stack([ex.ravel(), ny.ravel()], 1)
    v = (re(pts) + 1j * im(pts)).reshape(ROWS, COLS)
    return force * v / (2j * np.pi * f)


def vibrate_exact(slc, g, D, f, nmax=None):
    """exp(-i k Re(D e^{i w t})) = sum_n (-i)^n J_n(k|D|) e^{i n (w t + arg D)}: each term a field in place times a history
    in time, applied through the Doppler-time mapping t = f_D / Ka."""
    k = 4 * np.pi / g.lam
    a = k * np.abs(D)
    nmax = int(np.ceil(a.max() + 4 * a.max() ** (1 / 3) + 4)) if nmax is None else nmax
    fd = np.fft.fftfreq(ROWS, d=1.0 / g.prf)
    t = fd / g.Ka_signed
    phi = np.angle(D)
    out = np.zeros_like(slc, dtype=np.complex128)
    for n in range(-nmax, nmax + 1):
        field = ((-1j) ** n) * jv(n, a) * np.exp(1j * n * phi)
        hist = np.exp(1j * n * 2 * np.pi * f * t)
        out += np.fft.ifft(np.fft.fft(slc * field, axis=0) * hist[:, None], axis=0)
    return out, nmax


def main():
    p23 = json.loads((RESULTS / 'p2_23_every_reader' / 'summary.json').read_text())
    f = float(p23['manifest']['params']['frequency_hz'])
    F90 = float(p23['satellite']['reflectors']['bench']['force_n'])
    g = DwellGeometry.from_record(m.ACQ)
    st = EchoSetup.from_geometry(g)
    lat0, lon0 = g.source['scene_centre']['latitude'], g.source['scene_centre']['longitude']
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=SUPPORTS)
    forces = [F90 * s for s in FORCE_FRACTIONS]
    params = {'acquisition': m.ACQ, 'shape': [ROWS, COLS], 'frequency_hz': f, 'forces_n': forces, 'supports': SUPPORTS,
              'source_clear_m': SOURCE_CLEAR_M, 'p2_23_commit': p23['manifest']['commit'], 'step_m': STEP, 'half_m': HALF,
              'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE)}
    with Run(RID, 'The gated reconstruction on the one shaking', params) as run:
        gdir = Path(__file__).resolve().parents[1] / 'data' / RID
        gdir.mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(24)
        mask = e23.band_mask(g)
        base = e23.speckle(g, rng, mask)
        # the scene: the reflectors across the chamber, and natural bright points drawn from the real image's population
        s_off = np.arange(-e23.REFLECTOR_HALF_M, e23.REFLECTOR_HALF_M + 1e-9, e23.REFLECTOR_STEP_M)
        cols_r = (COLS // 2 + np.round(s_off * np.sin(g.theta) / g.dr)).astype(int)
        for c in cols_r:
            base = base + e23.point_image(mask, ROWS // 2, int(c), e23.REFLECTOR_SCR_DB)
        dens = p23['satellite']['natural']['per_m2'][f'{min(e23.CAL_SCR_DB):g}']
        scr_pool = p23['satellite']['natural']['scr_db']
        n_nat = rng.poisson(dens * (2 * 60.0) ** 2)
        for _ in range(n_nat):
            xy = rng.uniform(-60.0, 60.0, 2)
            r, c = b20.site_to_pixel(g, xy[0], xy[1])
            scr = float(rng.choice(scr_pool)) if scr_pool else 35.0
            base = base + e23.point_image(mask, int(round(r)) % ROWS, int(round(c)) % COLS, scr)

        geom = b20_geometry(g, lat0, lon0)
        (gdir / 'geometry.json').write_text(json.dumps(geom))
        gated = m.load_gated()
        gated.render_tomograms = lambda *a, **k: None
        focus, results, z_m = {}, {}, None
        for F in forces:
            tag = f'F{F:.0e}'.replace('+', '')
            D0 = los_displacement(g, 'none', F, f)
            D1 = los_displacement(g, 'bench', F, f)
            img0, n0 = vibrate_exact(base, g, D0, f)
            img1, n1 = vibrate_exact(base, g, D1, f)
            dchange = img1 - img0
            null = img0 + np.abs(dchange) * np.exp(2j * np.pi * rng.random(dchange.shape))
            for case, img in (('without', img0), ('with', img1), ('null', null)):
                name = f'{tag}_{case}'
                path = write_synthetic_slc(gdir / f'SYNTHETIC_{name}.h5', img.astype(np.complex64), st, g.Ka_signed, lat0, lon0,
                                           75.0, g.heading_deg, g.theta_deg,
                                           scene={'bench': 'bench-void', 'vibrator_force_n': F, 'frequency_hz': f,
                                                  'chamber': case == 'with', 'null': case == 'null', 'reflectors': len(cols_r),
                                                  'natural_points': int(n_nat)})
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
                results[name] = {'force_n': F, 'bessel_terms': int(max(n0, n1)),
                                 'supported_positions': {P: int(np.isfinite(focus[name][P]).any(axis=2).sum()) for P in SUPPORTS}}
                print(f"  {name}: {results[name]}", flush=True)

        # where the scores changed, against the product without the chamber at the same force
        n_pos = int(round(2 * HALF / STEP)) + 1
        ys = np.arange(-HALF, HALF + 1e-6, STEP)
        xs = np.linspace(HALF, -HALF, 10001)[np.linspace(0, 10000, n_pos).round().astype(int)]
        X, Y = np.meshgrid(xs, ys)
        near = np.hypot(X, Y) <= NEAR_M
        share = float(near.mean())
        for F in forces:
            tag = f'F{F:.0e}'.replace('+', '')
            for case in ('with', 'null'):
                A, B = focus[f'{tag}_{case}'][1], focus[f'{tag}_without'][1]
                diff = np.abs(np.nan_to_num(A) - np.nan_to_num(B))
                changed = diff.max(axis=2) > 1e-6
                results[f'{tag}_{case}']['p1'] = {'positions_changed': int(changed.sum()),
                                                  'positions_changed_near': int((changed & near).sum()),
                                                  'max_score_change': float(diff.max())}
        z_surf = np.zeros((len(ys), n_pos))
        zc = -STEP * np.arange(int(np.ceil(z_m[-1] / STEP)) + 1)
        vols = {f'{k}_p{P}': curtains_to_volume(v[P], z_m, z_surf, zc, STEP) for k, v in focus.items() for P in SUPPORTS}
        np.savez_compressed(Path(run.dir) / 'volumes.npz', **{k: v.astype(np.float16) for k, v in vols.items()},
                            x=xs[::-1], y=ys, z=zc, z_surface=z_surf[:, ::-1].T, step=STEP)
        gathers = lambda r: r['positions_changed_near'] >= 5 and r['positions_changed_near'] > 3 * share * r['positions_changed']
        tag = f'F{forces[0]:.0e}'.replace('+', '')
        w, nl = results[f'{tag}_with']['p1'], results[f'{tag}_null']['p1']
        finding = (
            f"At the force where corner reflectors let the favourable detector find the bench's chamber ({forces[0]:.1e} N), "
            f"the gated reconstruction, unchanged, passes {results[f'{tag}_without']['supported_positions'][1]} positions at a "
            f"support of one without the chamber and {results[f'{tag}_with']['supported_positions'][1]} with it; the chamber "
            f"changes the scores at {w['positions_changed']} positions ({w['positions_changed_near']} within {NEAR_M:.0f} m of "
            f"it, where {share * 100:.0f}% of positions lie), random noise of the same size at {nl['positions_changed']} "
            f"({nl['positions_changed_near']}). "
            + ('Its changes gather over the chamber.' if gathers(w) and not gathers(nl) else
               'Its changes do not gather over the chamber: where the data hold the chamber, the reconstruction does not read it.'))
        run.save({'results': results, 'near_share': share, 'finding': finding,
                  'grid': {'x0': float(xs[-1]), 'y0': float(ys[0]), 'z_top': float(zc[0]), 'step_m': STEP,
                           'shape': list(vols[next(iter(vols))].shape)}})
        print(finding)


def b20_geometry(g, lat0, lon0):
    """P2-20's raster: 51 east-west lines of 51 positions across the bench, placed by the product's own geometry."""
    n_pos = int(round(2 * HALF / STEP)) + 1
    ys = np.arange(-HALF, HALF + 1e-6, STEP)
    xs_trace = np.linspace(HALF, -HALF, 10001)
    origin = {'latitude': lat0, 'longitude': lon0}
    tracks = {}
    for j, y in enumerate(ys):
        r, c = b20.site_to_pixel(g, xs_trace, np.full_like(xs_trace, y))
        lat, lon = frame_to_lla(origin, xs_trace[[0, -1]], np.full(2, y))
        tracks[f'L{j:03d}'] = {'native_trace_col_row': np.stack([np.round(c), np.round(r)], 1).astype(int).tolist(),
                               'geographic_endpoint_lat_lon': [[float(lat[0]), float(lon[0])], [float(lat[1]), float(lon[1])]]}
    return {'title': 'East-west lines across the one-chamber bench, 2.4 m apart', 'target_samples_per_track': n_pos,
            'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
            'projection': {'height_m': 75.0}, 'tracks': tracks}


if __name__ == '__main__':
    main()
