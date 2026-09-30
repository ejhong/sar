"""The lab's processing command: run a named method over an area of a product and add the result to the viewer.

    uv run --with scikit-image python -m katabasis.lab gated --site giza --centre 0,0 --name khafre
    uv run --with scikit-image python -m katabasis.lab gated --site giza --centre 324.26,349.19 --twin 101 --name khufu

Methods:

  gated   the stricter gated reconstruction, its source imported as published and hashed, its profile unchanged but for
          the input paths and the supports reported (1 and 4 by default). Lines run east to west across the area (or
          north to south, --lines ns), `step` metres apart, each placed through the product's RPC over the site's
          surface with the site's geoid undulation (its frame's geoid_m; Giza's 15.5 m by default), shifted by any
          residual offset given (--shift, in lines and samples); the curtains each line yields are stacked into a volume on a `step`-metre
          grid: the method's fit score at each nominal depth below the surface point above, unsmoothed, empty where
          nothing passed. The reconstruction compares each position with its neighbours along its line, so the way the
          lines are laid is one of the things a real feature should survive.

Every run writes sim/results/lab_<name>/ (summary.json with its manifest, volumes.npz), which `katabasis.export` adds
to the site's viewer. With --twin SEED the same raster is also run on a motionless twin of the crop.

The product defaults to the 2022 X13 dwell on the desktop. --product takes any product in the ICEYE SLC layout: a real
one is placed through its RPC, a synthetic one (sarsim.slcfile) through the flat mapping its file records, its scene
centre at the site's origin. The gated profile's filter bank was frozen for the 2022 product and is used as frozen; on a
product with another Doppler rate or band it may not suit, and the command does not adapt it.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter1d

from .compose import load_site
from .export.sites import scene as site_scene
from .runs import Run

ROOT = Path(__file__).resolve().parents[2]
GATED = Path.home() / 'tmp/sar/biondi_v18/pkg/Biondi_Protocol_v1.8'
MODULE = GATED / 'biondi_tomography_v1_8.py'
PROFILE = GATED / 'profiles/iceye_x13_khafre_w50.json'
PRODUCT = Path.home() / 'tmp/sar/giza2/ICEYE_X13_SLC_SLED_868226_20220715T235744.h5'
GEOID_M = 15.5
TRACE_POINTS = 10001


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_gated():
    spec = importlib.util.spec_from_file_location('gated', MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


WGS84_A, WGS84_E2 = 6378137.0, 6.69437999014e-3


def synthetic_mapping(product) -> dict | None:
    """A synthetic product's flat site-to-pixel mapping (sarsim.slcfile) with its scene centre's place, or None."""
    import h5py
    with h5py.File(product, 'r') as h:
        if 'synthetic' not in h or 'site_to_pixel' not in h['synthetic'].attrs:
            return None
        m = json.loads(h['synthetic'].attrs['site_to_pixel'])
        cc = h['coord_center'][...]
        lat = np.deg2rad(float(cc[2]))
        x, y, _ = h['synthetic']['scene_centre_ecef'][...]
        n = WGS84_A / np.sqrt(1 - WGS84_E2 * np.sin(lat) ** 2)
        m['origin'] = {'latitude': float(cc[2]), 'longitude': float(cc[3])}
        m['height_m'] = float(np.hypot(x, y) / np.cos(lat) - n)
        return m


def flat_project(m, x, y, z):
    """Row and column of site points in a synthetic product: azimuth along the heading, slant range to its right."""
    h, th = np.deg2rad(m['heading_deg']), np.deg2rad(m['theta_deg'])
    a = x * np.sin(h) + y * np.cos(h)
    gr = x * np.cos(h) - y * np.sin(h)
    return m['centre_row'] + a / m['dx'], m['centre_col'] + (gr * np.sin(th) - z * np.cos(th)) / m['dr']


def site_geoid(sc) -> float:
    """The site's undulation from orthometric to ellipsoidal height: its frame's geoid_m, Giza's by default."""
    return float(sc['frame'].get('geoid_m', GEOID_M))


def raster(product, sc, centre, half_ew, half_ns, step, positions=101, lines='ew', shift=(0.0, 0.0)):
    """Lines across a square of the site, east to west (or north to south), placed through the product's RPC over the
    site's surface, or through a synthetic product's own flat mapping. Returns the geometry, the positions along a line,
    each line's offset across, and the surface height at every [line, position]."""
    from sarsim.ortho import frame_to_lla, project, surface
    syn = synthetic_mapping(product)
    if syn is None:
        from sarsim.dwell import DwellProduct
        dp = DwellProduct(product)
    cx, cy = centre
    ew = lines == 'ew'
    offs = np.arange(-half_ns, half_ns + 1e-6, step) if ew else np.arange(-half_ew, half_ew + 1e-6, step)
    trace = np.linspace(half_ew, -half_ew, TRACE_POINTS) if ew else np.linspace(half_ns, -half_ns, TRACE_POINTS)
    sel = np.linspace(0, TRACE_POINTS - 1, positions).round().astype(int)
    pos = trace[sel]
    o = syn['origin'] if syn else sc['frame']['origin']
    geom = {'title': f"{'East-west' if ew else 'North-south'} lines {step:g} m apart across ({cx:g}, {cy:g})",
            'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
            'tracks': {}}
    heights, z_surf = [], np.zeros((len(offs), positions))
    for j, off in enumerate(offs):
        X = cx + trace if ew else np.full_like(trace, cx + off)
        Y = np.full_like(trace, cy + off) if ew else cy + trace
        Z = surface(sc, X, Y)
        if syn:
            r, c = flat_project(syn, X, Y, Z)
        else:
            r, c = project(dp, sc, X, Y, Z, site_geoid(sc), shift)
        heights.append(syn['height_m'] if syn else float(np.mean(Z)) + site_geoid(sc))
        z_surf[j] = Z[sel]
        lat, lon = frame_to_lla(o, X[[0, -1]], Y[[0, -1]])
        geom['tracks'][f'L{j:03d}'] = {
            'native_trace_col_row': np.stack([np.round(c), np.round(r)], 1).astype(int).tolist(),
            'geographic_endpoint_lat_lon': [[float(lat[0]), float(lon[0])], [float(lat[1]), float(lon[1])]]}
    geom['projection'] = {'height_m': float(np.mean(heights))}
    return geom, pos, offs, z_surf


def curtains_to_volume(F, z_m, z_surf, zc, step, lines='ew'):
    """Each position's depth profile hung below its surface point, sampled at the voxel centres after taking the largest
    score within half a voxel either way; NaN where the method gave no score. F is [line, position, depth], positions
    running east to west (or north to south); returns (x west to east, y south to north, z down)."""
    dz = float(z_m[1] - z_m[0])
    half = int(round(step / 2 / dz))
    filled = np.where(np.isfinite(F), F, -np.inf)
    peak = maximum_filter1d(filled, size=2 * half + 1, axis=2, mode='nearest')
    peak = np.where(np.isfinite(peak), peak, np.nan)
    depth = z_surf[:, :, None] - zc[None, None, :]
    idx = np.round((depth - z_m[0]) / dz).astype(int)
    valid = (idx >= 0) & (idx < len(z_m))
    vals = np.where(valid, np.take_along_axis(peak, np.clip(idx, 0, len(z_m) - 1), axis=2), np.nan)[:, ::-1, :]
    return (vals.transpose(1, 0, 2) if lines == 'ew' else vals).astype(np.float32)


def run_gated(name, site='giza', centre=(0.0, 0.0), half_ew=150.0, half_ns=141.0, step=3.0, supports=(1, 4),
              twin=None, product=PRODUCT, title=None, lines='ew', area=None, shift=(0.0, 0.0)):
    from sarsim.looks import motionless_twin
    if not re.fullmatch(r'[a-z0-9_]+', name):
        raise SystemExit('a run name is lower-case letters, digits and underscores')
    if lines not in ('ew', 'ns'):
        raise SystemExit('lines run east-west (ew) or north-south (ns)')
    along = half_ew if lines == 'ew' else half_ns
    if abs(2 * along / 100 - step) > 1e-9:
        flag = '--half-ew' if lines == 'ew' else '--half-ns'
        raise SystemExit(f'positions along a line are 2 x {flag[2:]} / 100 apart; for a {step:g} m grid use {flag} {50 * step:g}')
    sc = site_scene(load_site(site))
    profile = json.loads(PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=list(supports))
    rid = f'lab_{name}'
    params = {'method': 'gated', 'site': site, 'centre_m': list(centre), 'half_ew_m': half_ew, 'half_ns_m': half_ns,
              'step_m': step, 'supports': list(supports), 'twin_seed': twin, 'product': Path(product).name, 'lines': lines,
              'gated_module_sha256': sha256(MODULE), 'gated_profile_sha256': sha256(PROFILE), 'geoid_m': site_geoid(sc)}
    if area:
        params['area'] = area
    if any(shift):
        params['shift_px'] = [float(v) for v in shift]
    with Run(rid, title or f'The gated reconstruction across ({centre[0]:g}, {centre[1]:g}) on {site}', params) as run:
        geom, pos, offs, z_surf = raster(product, sc, centre, half_ew, half_ns, step, lines=lines, shift=shift)
        gdir = ROOT / 'sim' / 'data' / rid
        gdir.mkdir(parents=True, exist_ok=True)
        (gdir / 'geometry.json').write_text(json.dumps(geom))
        prof = dict(profile)
        prof['input'] = dict(profile['input'], hdf5_path=str(product), geometry_path=str(gdir / 'geometry.json'))
        (gdir / 'profile.json').write_text(json.dumps(prof, indent=1))
        gated = load_gated()
        original = gated.load_source_and_preflight
        gated.render_tomograms = lambda *a, **k: None
        cases = ['real'] + ([f'twin{twin}'] if twin is not None else [])
        vols, results = {}, {}
        for case in cases:
            def wrapped(config, config_path, case=case):
                out = list(original(config, config_path))
                if case.startswith('twin'):
                    out[7] = motionless_twin(out[7], np.random.default_rng(int(case[4:])))
                return tuple(out)
            gated.load_source_and_preflight = wrapped
            gated.run_full(prof, gdir / 'profile.json', gdir / case)
            audit = np.load(gdir / case / 'biondi_v1_8_audit.npz')
            z_m = audit['z_m']
            z_top = float(z_surf.max())
            zc = z_top - step * np.arange(int(np.ceil((z_top - (z_surf.min() - z_m[-1])) / step)) + 1)
            results[case] = {'supported_positions': {P: int(np.isfinite(audit[f'focus_p{P}']).any(axis=2).sum()) for P in supports}}
            for P in supports:
                vols[f'{case}_p{P}'] = curtains_to_volume(audit[f'focus_p{P}'], z_m, z_surf, zc, step, lines)
            print(f"  {case}: supported positions " + ', '.join(f"P{P} {results[case]['supported_positions'][P]}" for P in supports),
                  flush=True)
        gated.load_source_and_preflight = original
        ew = lines == 'ew'
        x = centre[0] + (pos[::-1] if ew else offs)
        y = centre[1] + (offs if ew else pos[::-1])
        np.savez_compressed(Path(run.dir) / 'volumes.npz', **{k: v.astype(np.float16) for k, v in vols.items()},
                            x=x, y=y, z=zc, z_surface=(z_surf[:, ::-1].T if ew else z_surf[:, ::-1]), step=step)
        n = len(offs) * len(pos)
        finding = (f"The gated reconstruction, unchanged, along {len(offs)} {'east-west' if ew else 'north-south'} lines {step:g} m apart across "
                   f"({centre[0]:g}, {centre[1]:g}) on {site}: " + '; '.join(
                       f"{case}: " + ', '.join(f"{results[case]['supported_positions'][P]:,} of {n:,} positions pass at support {P}"
                                               for P in supports) for case in cases) + '.')
        run.save({'method': 'gated', 'results': results, 'cases': cases,
                  'grid': {'x0': float(x[0]), 'y0': float(y[0]), 'z_top': float(zc[0]),
                           'step_m': step, 'shape': list(vols['real_p1' if 1 in supports else f'real_p{supports[0]}'].shape)},
                  'finding': finding})
        print(finding)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='katabasis.lab', description='Run a method over an area of a product and add it to the lab.')
    sub = ap.add_subparsers(dest='method', required=True)
    g = sub.add_parser('gated', help='the stricter gated reconstruction, run unchanged')
    g.add_argument('--name', required=True, help='run name: lab_<name> in sim/results')
    g.add_argument('--site', default='giza')
    g.add_argument('--centre', default='0,0', help='x,y in the site frame, metres')
    g.add_argument('--half-ew', type=float, default=150.0)
    g.add_argument('--half-ns', type=float, default=141.0)
    g.add_argument('--step', type=float, default=3.0)
    g.add_argument('--supports', default='1,4')
    g.add_argument('--twin', type=int, default=None, help='also run a motionless twin with this seed')
    g.add_argument('--product', default=str(PRODUCT))
    g.add_argument('--title', default=None)
    g.add_argument('--lines', default='ew', choices=['ew', 'ns'], help='lines east-west (default) or north-south')
    g.add_argument('--area', default=None, help='the place the lab names this run by (a name the other methods use)')
    g.add_argument('--shift', default='0,0', help='residual offset in lines,samples added to the RPC placement')
    a = ap.parse_args(argv)
    cx, cy = (float(v) for v in a.centre.split(','))
    run_gated(a.name, a.site, (cx, cy), a.half_ew, a.half_ns, a.step, tuple(int(s) for s in a.supports.split(',')),
              a.twin, Path(a.product), a.title, a.lines, a.area, tuple(float(v) for v in a.shift.split(',')))


if __name__ == '__main__':
    main()
