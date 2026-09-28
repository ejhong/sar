"""The satellite, for the viewer: one radar pass over a site, and what the published method draws from it.

bench-void (P2-07, P2-09): the satellite's track and beam from the real Giza acquisition record, the simulated image draped
on the ground, a grid of the image's virtual sensors with the true motion under each and what each reports, and the
published method's depth volume with the chamber and without it, resampled onto the site grid.

bench-khafre-claim (P2-08): the published method's volume from the real Khafre image, its depth relabelled so that the
steering's repeat falls at 648 m, placed where the claimed shafts stand.

Writes sites/<site>/radar.json (merged into scene.json by export_sites), the draped image, and radar volumes appended to
volumes.json with status 'radar' (drawn in cinnabar).
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

from ..runs import RESULTS, load
from .sites import DATA


def _geometry():
    from sarsim.acquisition import DwellGeometry
    return DwellGeometry.from_record('giza-20250827')


def _log_u8(T: np.ndarray) -> tuple[np.ndarray, list[float]]:
    """The published-style scale: log of the power over its median, 2nd to 99.5th percentile, as 0..255."""
    A = np.log10(np.maximum(T, 1e-30) / np.median(T))
    lo, hi = np.percentile(A, 2), np.percentile(A, 99.5)
    return np.round(255 * np.clip((A - lo) / (hi - lo), 0, 1)).astype(np.uint8), [float(lo), float(hi)]


def _resample(T3: np.ndarray, xr: np.ndarray, yr: np.ndarray, z: np.ndarray, g, half: float, h: float,
              zmax: float) -> tuple[np.ndarray, list[float], list[int]]:
    """A radar-grid volume T3[azimuth, range, depth] onto the site grid (east, north, depth), h metres a voxel."""
    ax = np.arange(-half + h / 2, half, h)
    zz = np.arange(h / 2, zmax, h)
    E, N, Z = np.meshgrid(ax, ax, zz, indexing='ij')
    a, r = g.site_to_radar(E, N)
    f = RegularGridInterpolator((xr, yr, z), T3, bounds_error=False, fill_value=np.nan)
    V = f(np.stack([a, r, Z], axis=-1))
    inside = np.isfinite(V)
    V[~inside] = np.nanmin(V) if inside.any() else 0.0
    return V, [float(ax[0]), float(ax[0]), 0.0], [len(ax), len(ax), len(zz)]


def _write_volume(d: Path, vid: str, V: np.ndarray, u8: np.ndarray, origin: list[float], h: float, info: dict) -> dict:
    nx, ny, nz = V.shape
    (d / f'vol-{vid}.u8').write_bytes(np.ascontiguousarray(u8.transpose(2, 1, 0)).tobytes())
    ztop = -h / 2                      # centre of the top layer, just below the ground at z = 0
    return {'id': vid, 'status': 'radar', 'file': f'vol-{vid}.u8', 'shape': [nx, ny, nz],
            'origin': [origin[0], origin[1], ztop], 'spacing': h, **info}


def _merge_volumes(d: Path, new: list[dict]):
    path = d / 'volumes.json'
    vols = json.loads(path.read_text()) if path.exists() else []
    vols = [v for v in vols if v.get('status') != 'radar'] + new
    path.write_text(json.dumps(vols, indent=1))


def export_bench(out: Path = DATA) -> dict | None:
    rid7, rid9 = 'p2_07_whole_chain', 'p2_09_virtual_sensors'
    vpath, spath = RESULTS / rid7 / 'volumes.npz', RESULTS / rid9 / 'sensors.npz'
    if not (vpath.exists() and spath.exists()):
        return None
    g = _geometry()
    s7, s9 = load(rid7), load(rid9)
    d = out / 'sites' / 'bench-void'
    d.mkdir(parents=True, exist_ok=True)
    # the method's volumes, with and without the chamber
    v = np.load(vpath)
    nr, nc = len(v['rows']), len(v['cols'])
    shape = s7['manifest']['params']['image_px']
    xr = (v['rows'] - shape[0] // 2) * float(v['dx'])
    yr = (v['cols'] - shape[1] // 2) * float(v['dr']) / np.sin(float(v['theta']))
    z = v['z']
    run = f"{rid7} · {s7['manifest']['date']} · {s7['manifest']['commit']}"
    vols = []
    both = {k: gaussian_filter(v[k].reshape(nr, nc, -1), (0.6, 0.6, 1.0)) for k in ('A', 'B')}
    ref = _log_u8(both['B'])[1]
    for key, label, cap in (('A', 'Satellite · the published method, with the chamber',
                             'its depth picture of this ground: pillars and streaks everywhere, the chamber nowhere'),
                            ('B', 'Satellite · the published method, without the chamber',
                             'the same ground with no chamber: the same picture, to one part in a hundred million')):
        V, origin, _ = _resample(both[key], xr, yr, z, g, 50.0, 1.0, float(z[-1]))
        A = np.log10(np.maximum(V, 1e-30) / np.median(both['B']))
        u8 = np.round(255 * np.clip((A - ref[0]) / (ref[1] - ref[0]), 0, 1)).astype(np.uint8)
        vols.append(_write_volume(d, f'radar-{key}', V, u8, origin, 1.0, {
            'label': label, 'method': 'One ICEYE dwell, sub-aperture pairs registered and focused as published (lambda_s 0.48 m)',
            'quantity': 'focused power, log scale', 'units': 'relative', 'range': ref, 'run': run, 'caption': cap}))
    _merge_volumes(d, vols)
    # the draped image: north up (rows reversed), rotated by the heading in the viewer
    src = RESULTS / rid7 / 'figs' / 'simulated_ground.png'
    from PIL import Image
    im = Image.open(src).convert('RGB').transpose(Image.FLIP_TOP_BOTTOM)
    im.save(d / 'radar-ground.jpg', quality=88, optimize=True)
    # the sensors, thinned to a readable grid
    sv = np.load(spath)
    gr, gc = (int(x) for x in sv['grid'])
    keep = np.zeros((gr, gc), bool)
    keep[::4, ::2] = True
    keep = keep.ravel()
    to_um = lambda a: np.round(a[keep] * 1e6, 1).tolist()
    radar = {
        'kind': 'bench',
        'acquisition': {'name': 'giza-20250827', 'heading_deg': g.heading_deg, 'incidence_deg': g.theta_deg,
                        'los_enu': list(g.los_enu), 'aperture_s': g.aperture_time, 'track_km': g.source['baseline_span_m'] / 1e3,
                        'slant_range_km': g.R0 / 1e3, 'along_track_en': g.along_track_en.tolist(),
                        'ground_range_en': g.ground_range_en.tolist()},
        'image': {'file': 'radar-ground.jpg', 'centre': [0.0, 0.0], 'width_m': s7['geometry']['ground_range_m'],
                  'height_m': s7['geometry']['azimuth_m'], 'rotation_deg': float(np.rad2deg(np.arctan2(g.ground_range_en[1], g.ground_range_en[0])))},
        'sensors': {'east': np.round(sv['east'][keep], 2).tolist(), 'north': np.round(sv['north'][keep], 2).tolist(),
                    'looks_s': np.round(sv['looks_s'], 3).tolist(), 'truth_um_s': to_um(sv['truth']),
                    'complex_um_s': to_um(sv['complex']), 'magnitude_um_s': to_um(sv['magnitude']),
                    'test_wave': s9['manifest']['params']['test_wave'],
                    'gains': s9['gains'], 'run': f"{rid9} · {s9['manifest']['date']} · {s9['manifest']['commit']}"},
        'volumes': {'with': 'radar-A', 'without': 'radar-B'},
    }
    (d / 'radar.json').write_text(json.dumps(radar, separators=(',', ':')))
    return {'site': 'bench-void', 'volumes': [x['id'] for x in vols], 'sensors': int(keep.sum())}


def export_claim(out: Path = DATA) -> dict | None:
    rid = 'p2_08_pillars_and_boxes'
    vpath = RESULTS / rid / 'khafre_volume.npz'
    if not vpath.exists():
        return None
    g = _geometry()
    s8 = load(rid)
    k = s8['sets']['khafre']
    v = np.load(vpath)
    T = v['T']
    nr, nc, nz = T.shape
    pa, pr = k['pixel_m']
    xr = (np.arange(nr) - nr / 2) * pa
    yr = (np.arange(nc) - nc / 2) * pr
    z = v['z_shown']
    T = gaussian_filter(T, (1.0, 1.0, 1.5))
    V, origin, shape = _resample(T, xr, yr, z, g, 300.0, 8.0, 1400.0)
    u8, rng = _log_u8(T)
    A = np.log10(np.maximum(V, 1e-30) / np.median(T))
    u8 = np.round(255 * np.clip((A - rng[0]) / (rng[1] - rng[0]), 0, 1)).astype(np.uint8)
    d = out / 'sites' / 'bench-khafre-claim'
    d.mkdir(parents=True, exist_ok=True)
    vol = _write_volume(d, 'radar-khafre', V, u8, origin, 8.0, {
        'label': 'Satellite · the published method on the real Khafre image',
        'method': 'The first investigation’s trajectories on the ICEYE Giza dwell, focused past the repeat depth and relabelled',
        'quantity': 'focused power, log scale', 'units': 'relative', 'range': rng,
        'run': f"{rid} · {s8['manifest']['date']} · {s8['manifest']['commit']}",
        'caption': f"pillars where the surface reading is noisiest; blocks at {k['repeat_depth_m_shown']:,.0f} m and "
                   f"{2 * k['repeat_depth_m_shown']:,.0f} m, where the depth axis repeats"})
    _merge_volumes(d, [vol])
    return {'site': 'bench-khafre-claim', 'volumes': [vol['id']]}


def export_radar(out: Path = DATA) -> list[dict]:
    return [r for r in (export_bench(out), export_claim(out)) if r]
