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
import re
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


def _write_volume(d: Path, vid: str, V: np.ndarray, u8: np.ndarray, origin: list[float], h: float, info: dict,
                  ztop: float | None = None) -> dict:
    nx, ny, nz = V.shape
    (d / f'vol-{vid}.u8').write_bytes(np.ascontiguousarray(u8.transpose(2, 1, 0)).tobytes())
    ztop = -h / 2 if ztop is None else ztop      # centre of the top layer (just below the ground at z = 0 on a bench)
    return {'id': vid, 'status': 'radar', 'file': f'vol-{vid}.u8', 'shape': [nx, ny, nz],
            'origin': [origin[0], origin[1], ztop], 'spacing': h, **info}


def _acquisition(g) -> dict:
    return {'name': g.name, 'heading_deg': g.heading_deg, 'incidence_deg': g.theta_deg, 'los_enu': list(g.los_enu),
            'aperture_s': g.aperture_time, 'track_km': g.source['baseline_span_m'] / 1e3, 'slant_range_km': g.R0 / 1e3,
            'along_track_en': g.along_track_en.tolist(), 'ground_range_en': g.ground_range_en.tolist(),
            'satellite': g.source['satellite'], 'date': g.source['collection_start'][:10]}


def _merge_volumes(d: Path, new: list[dict], drop=None):
    """Replace a site's radar volumes with `new`; with `drop`, only the existing volumes it selects."""
    path = d / 'volumes.json'
    vols = json.loads(path.read_text()) if path.exists() else []
    drop = drop or (lambda v: v.get('status') == 'radar')
    vols = [v for v in vols if not drop(v)] + new
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
        'acquisition': _acquisition(g),
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


def export_real(out: Path = DATA, site: str = 'giza') -> dict | None:
    """The real Giza pass (P2-12): the image on the ground, and the published method's volume at each patch inside the
    site, on one brightness scale so that monuments and open plateau compare fairly. Depth is measured from the rock
    surface and relabelled as the claim relabels it (the axis repeats at 648 m), down to the site's floor: Giza stops at
    its block's floor; Giza's depths runs the axis past its repeats, as the published pictures do."""
    rid = 'p2_12_real_pass'
    vpath = RESULTS / rid / 'volumes.npz'
    if not vpath.exists():
        return None
    from ..compose import load_site
    from .sites import scene as site_scene
    from sarsim.acquisition import DwellGeometry
    from sarsim.ortho import grid_height
    s = load(rid)
    sc = site_scene(load_site(site))
    g = DwellGeometry.from_record(s['manifest']['params']['acquisition'])
    v = np.load(vpath)
    h = float(v['voxel_m'])
    gx0, gy0 = float(v['grid_x0']), float(v['grid_y0'])
    nz_axis = len(v['z_raw'])
    per_index = float(s['manifest']['params']['claim_repeat_m']) / nz_axis
    shown = [p for p in s['patches'] if p['cells_on_site'] > 0]
    # smoothed as the claim's volume is (export_claim): a pixel or so across, a step and a half in depth, the depth axis
    # wrapping round, since it covers one whole period
    def smoothed(key):
        T = v[f'{key}_T'].astype(np.float32)
        nr, nc = (int(n) for n in v[f'{key}_grid'])
        return gaussian_filter(T.reshape(nr, nc, -1), (1.0, 1.0, 1.5), mode=('nearest', 'nearest', 'wrap')).reshape(nr * nc, -1)
    Ts = {p['key']: smoothed(p['key']) for p in shown}
    allT = np.concatenate([Ts[p['key']].ravel() for p in shown])
    med = float(np.median(allT))
    A = np.log10(np.maximum(allT, 1e-30) / med)
    lo, hi = float(np.percentile(A, 2)), float(np.percentile(A, 99.5))
    d = out / 'sites' / site
    d.mkdir(parents=True, exist_ok=True)
    zbot = sc['extent']['z'][0]
    run = f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}"
    vols = []
    for p in shown:
        key = p['key']
        idx = v[f'{key}_index']
        T = Ts[key]
        iy, ix = np.nonzero(idx >= 0)
        sub = idx[iy.min():iy.max() + 1, ix.min():ix.max() + 1]
        xs = gx0 + h / 2 + h * np.arange(ix.min(), ix.max() + 1)
        ys = gy0 + h / 2 + h * np.arange(iy.min(), iy.max() + 1)
        X, Y = np.meshgrid(xs, ys)
        ground = grid_height(sc['terrain'], X, Y) if sc['terrain']['kind'] == 'grid' else np.full(X.shape, sc['terrain']['z'])
        ztop = float(np.ceil(ground.max() / h) * h)
        nz = int((ztop - zbot) / h)
        zc = ztop - h / 2 - h * np.arange(nz)
        depth = ground[..., None] - zc[None, None, :]
        k = np.round(depth / per_index).astype(int) % nz_axis
        ok = (depth >= 0) & (sub[..., None] >= 0)
        vals = np.where(ok, T[np.maximum(sub, 0)[..., None], k], np.nan)
        a = np.log10(np.maximum(vals, 1e-30) / med)
        u8 = np.where(ok, np.round(255 * np.clip((a - lo) / (hi - lo), 0, 1)), 0).astype(np.uint8)
        u8 = u8.transpose(1, 0, 2)                                    # (x, y, z)
        control = p['kind'] == 'control'
        vols.append(_write_volume(d, f'radar-{key}', u8, u8, [float(xs[0]), float(ys[0])], h, {
            'label': f"Satellite · the published method {'over' if control else 'at'} {p['label']}{' (control)' if control else ''}",
            'method': 'The real 2025 ICEYE pass through the paper-style pipeline (50 half-band pairs, 32 px patches, grid-median '
                      'correction, no selection gates, lambda_s 0.48 m): focused power on a log scale, depth relabelled as the '
                      'claim does (repeat at 648 m), smoothed for display; not the gated reconstruction',
            'quantity': 'focused power, log scale', 'units': 'relative', 'range': [lo, hi], 'run': run,
            'caption': ('the same kind of picture over empty plateau' if control else
                        'columns where the registration wanders, as over empty plateau')}, ztop=ztop - h / 2))
    _merge_volumes(d, vols)
    src = RESULTS / rid / 'figs' / 'ground.jpg'
    shutil.copyfile(src, d / 'radar-ground.jpg')
    radar = {
        'kind': 'real',
        'acquisition': _acquisition(g),
        'image': {'kind': 'ortho', 'file': 'radar-ground.jpg', 'extent': s['image']['extent']},
        'reach_m': s['reach_m'],
        'volumes': [x['id'] for x in vols],
        # deep enough to see the whole depth axis: the image, draped on the ground, would hide what hangs beneath it
        **({'image_hidden': True} if sc['extent']['z'][1] - sc['extent']['z'][0] > 1000 else {}),
        'stats': {'monument_vs_control_profile_corr': s['monument_vs_control_profile_corr'],
                  'patch_profile_corr_range': s['patch_profile_corr_range'],
                  'pillar_power_vs_energy_min': min(x['pillar_power_vs_energy'] for x in s['patches'])},
        'run': run,
    }
    (d / 'radar.json').write_text(json.dumps(radar, separators=(',', ':')))
    return {'site': site, 'volumes': [x['id'] for x in vols]}


def export_khufu_gated(out: Path = DATA) -> dict | None:
    """The gated reconstruction across the Great Pyramid (P2-16): its fit scores hung below the surface on a 3 m grid, for
    the real 2022 image and a motionless copy, at supports one and four. Unsmoothed, empty where nothing passed, cut at
    the site's floor; drawn in gold beside the paper-style pipeline's cinnabar."""
    rid = 'p2_16_khufu_volume'
    vpath = RESULTS / rid / 'volumes.npz'
    rj = out / 'sites' / 'giza' / 'radar.json'
    if not (vpath.exists() and rj.exists()):
        return None
    from ..compose import load_site
    from .sites import scene as site_scene
    s = load(rid)
    sc = site_scene(load_site('giza'))
    v = np.load(vpath)
    h = float(v['step'])
    x, y, z = v['x'], v['y'], v['z']
    keep = z >= sc['extent']['z'][0] + h / 2
    d = out / 'sites' / 'giza'
    run = f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}"
    vols, entries = [], []
    for case, name in (('real', 'the real image'), ('twin101', 'a motionless copy')):
        for P in s['manifest']['params']['supports']:
            V = v[f'{case}_p{P}'].astype(np.float32)[:, :, keep]
            u8 = np.where(np.isfinite(V), 1 + np.round(254 * np.clip(np.nan_to_num(V, nan=0.0), 0, 1)), 0).astype(np.uint8)
            short = 'real' if case == 'real' else 'twin'
            vid = f'radar-khufu-gated-{short}-p{P}'
            vols.append(_write_volume(d, vid, V, u8, [float(x[0]), float(y[0])], h, {
                'label': f'Satellite · the gated reconstruction at Khufu, {name}, support {P}',
                'method': ("The 2022 ICEYE pass through the gated reconstruction's own code, unchanged (317 pairs; ellipse, "
                           f"Track10 and support-{P} gates), along 95 lines across the Great Pyramid 3 m apart"),
                'quantity': 'conditional adjusted R2 at each nominal depth below the surface', 'units': '0 to 1',
                'range': [0.0, 1.0], 'run': run, 'tint': 'gated',
                'caption': 'fit scores, unsmoothed; empty where nothing passed'}, ztop=float(z[keep][0])))
            entries.append({'id': vid, 'case': short, 'support': int(P)})
    # the same raster over open plateau (P2-17), where no monument stands: the in-image control
    cpath = RESULTS / 'p2_17_plateau_control' / 'volumes.npz'
    if cpath.exists():
        c = load('p2_17_plateau_control')
        w = np.load(cpath)
        cx_, cy_ = c['manifest']['params']['centre_m']
        ckeep = w['z'] >= sc['extent']['z'][0] + h / 2
        crun = f"p2_17_plateau_control · {c['manifest']['date']} · {c['manifest']['commit']}"
        for P in c['manifest']['params']['supports']:
            V = w[f'real_p{P}'].astype(np.float32)[:, :, ckeep]
            u8 = np.where(np.isfinite(V), 1 + np.round(254 * np.clip(np.nan_to_num(V, nan=0.0), 0, 1)), 0).astype(np.uint8)
            vid = f'radar-khufu-gated-plateau-p{P}'
            vols.append(_write_volume(d, vid, V, u8, [float(w['x'][0]), float(w['y'][0])], h, {
                'label': f'Satellite · the gated reconstruction over open plateau, support {P}',
                'method': "The same raster as at Khufu, laid on open plateau south-west of Menkaure, through the same unchanged code",
                'quantity': 'conditional adjusted R2 at each nominal depth below the surface', 'units': '0 to 1',
                'range': [0.0, 1.0], 'run': crun, 'tint': 'gated',
                'caption': 'the in-image control: no monument stands here'}, ztop=float(w['z'][ckeep][0])))
            entries.append({'id': vid, 'case': 'plateau', 'support': int(P),
                            'focus': [cx_, cy_, float(np.median(w['z_surface']))]})
    _merge_volumes(d, vols, drop=lambda q: q['id'].startswith('radar-khufu-gated'))
    radar = json.loads(rj.read_text())
    kh = next(t for t in sc['structures'] if t['id'] == 'khufu')['shape']
    st = s['stats']
    radar['gated'] = {
        'title': 'What the gated reconstruction computes at Khufu',
        'date': '2022',
        'volumes': entries,
        'focus': [kh['centre'][0], kh['centre'][1], kh['centre'][2] + 0.2 * kh['height']],
        'radius_m': 620.0,
        'chambers': {'real': [st['real_p1']['mean_score_in_chambers'], st['real_p1']['mean_score_same_depths_elsewhere']],
                     'twin': [st['twin101_p1']['mean_score_in_chambers'], st['twin101_p1']['mean_score_same_depths_elsewhere']]},
        'note': ("The stricter reconstruction's own code, unchanged, run on the 2022 pass along 95 lines across the pyramid, "
                 "3 m apart. Each voxel is its fit score at that nominal depth below the surface point above it: what the "
                 "pipeline computes, not a detection and not a validated depth. Beside it, the same calculation on a "
                 "motionless copy of the image. The satellite drawn above is the 2025 pass."),
        'run': run,
    }
    rj.write_text(json.dumps(radar, separators=(',', ':')))
    return {'site': 'giza', 'volumes': [e['id'] for e in entries]}


def export_bench_gated(out: Path = DATA) -> dict | None:
    """The gated reconstruction over the one-chamber bench (P2-20): its fit scores on synthetic products without the
    chamber, with it, with its imprint boosted, and with a random perturbation of the imprint's size, hung below the
    flat ground on a 2.4 m grid and cut at the bench's floor. Drawn with the satellite's other volumes."""
    rid = 'p2_20_bench_gated'
    vpath = RESULTS / rid / 'volumes.npz'
    rj = out / 'sites' / 'bench-void' / 'radar.json'
    if not (vpath.exists() and rj.exists()):
        return None
    from ..compose import load_site
    from .sites import scene as site_scene
    s = load(rid)
    sc = site_scene(load_site('bench-void'))
    v = np.load(vpath)
    h = float(v['step'])
    keep = v['z'] >= sc['extent']['z'][0] + h / 2
    d = out / 'sites' / 'bench-void'
    run = f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}"
    vols, entries = [], []
    for case, boost in s['manifest']['params']['cases']:
        kind = 'null' if boost is None else 'without' if boost == 0 else 'with' if boost == 1 else 'boosted'
        name = ("a random perturbation of the imprint's size" if kind == 'null' else 'without the chamber' if kind == 'without'
                else 'with the chamber' if kind == 'with' else f"with the chamber's imprint boosted {boost:,.0f} times")
        for P in s['manifest']['params']['supports']:
            V = v[f'{case}_p{P}'].astype(np.float32)[:, :, keep]
            u8 = np.where(np.isfinite(V), 1 + np.round(254 * np.clip(np.nan_to_num(V, nan=0.0), 0, 1)), 0).astype(np.uint8)
            vid = f'radar-gated-{case.replace("_", "-")}-p{P}'
            vols.append(_write_volume(d, vid, V, u8, [float(v['x'][0]), float(v['y'][0])], h, {
                'label': f'Satellite · the gated reconstruction over the bench, {name}, support {P}',
                'method': ("Synthetic products in the 2022 layout, shaken as Giza shakes, through the gated reconstruction's "
                           "own code, unchanged, along 51 lines 2.4 m apart"),
                'quantity': 'conditional adjusted R2 at each nominal depth below the surface', 'units': '0 to 1',
                'range': [0.0, 1.0], 'run': run, 'tint': 'gated',
                'caption': 'fit scores, unsmoothed; empty where nothing passed'}, ztop=float(v['z'][keep][0])))
            entries.append({'id': vid, 'case': kind, 'support': int(P), 'shaking': 'ambient',
                            **({'boost': float(boost)} if kind == 'boosted' else {})})
    # the one shaking (P2-24): the vibrator at the force where reflectors let the favourable detector find the chamber
    rid24 = 'p2_24_gated_on_the_shaking'
    v24path = RESULTS / rid24 / 'volumes.npz'
    if v24path.exists():
        s24 = load(rid24)
        v24 = np.load(v24path)
        keep24 = v24['z'] >= sc['extent']['z'][0] + h / 2
        run24 = f"{rid24} · {s24['manifest']['date']} · {s24['manifest']['commit']}"
        F = s24['manifest']['params']['forces_n'][0]
        tag = f'F{F:.0e}'.replace('+', '')
        for case, kind in (('without', 'without'), ('with', 'with'), ('null', 'null')):
            for P in s24['manifest']['params']['supports']:
                V = v24[f'{tag}_{case}_p{P}'].astype(np.float32)[:, :, keep24]
                u8 = np.where(np.isfinite(V), 1 + np.round(254 * np.clip(np.nan_to_num(V, nan=0.0), 0, 1)), 0).astype(np.uint8)
                vid = f'radar-gated-vibrator-{case}-p{P}'
                vols.append(_write_volume(d, vid, V, u8, [float(v24['x'][0]), float(v24['y'][0])], h, {
                    'label': f'Satellite · the gated reconstruction over the bench, shaken by a vibrator, {kind}, support {P}',
                    'method': ("Synthetic products of the bench shaken by P2-22's vibrator at the reflectors' boundary force, with the "
                               "reflectors and the real image's bright points, through the gated reconstruction's own code, unchanged"),
                    'quantity': 'conditional adjusted R2 at each nominal depth below the surface', 'units': '0 to 1',
                    'range': [0.0, 1.0], 'run': run24, 'tint': 'gated',
                    'caption': 'fit scores, unsmoothed; empty where nothing passed'}, ztop=float(v24['z'][keep24][0])))
                entries.append({'id': vid, 'case': kind, 'support': int(P), 'shaking': 'vibrator', 'note': s24['finding']})
    _merge_volumes(d, vols, drop=lambda q: q['id'].startswith('radar-gated-'))
    radar = json.loads(rj.read_text())
    radar['gated'] = {
        'title': 'What the gated reconstruction computes over the bench',
        'date': 'synthetic',
        'volumes': entries,
        'focus': [0.0, 0.0, -15.0],
        'radius_m': 200.0,
        'note': s['finding'],
        'run': run,
    }
    rj.write_text(json.dumps(radar, separators=(',', ':')))
    return {'site': 'bench-void', 'volumes': [e['id'] for e in entries]}


def export_lab(out: Path = DATA) -> list[dict]:
    """Every lab run (katabasis.lab, sim/results/lab_*): its volumes into its site's viewer, drawn in gold as the gated
    reconstruction's fit scores, and a list of runs in the site's radar.json for the satellite panel."""
    from ..compose import load_site
    from .sites import scene as site_scene
    runs = sorted(d for d in RESULTS.glob('lab_*') if (d / 'volumes.npz').exists() and (d / 'summary.json').exists())
    by_site: dict[str, list] = {}
    for d in runs:
        s = load(d.name)
        prm = s['manifest']['params']
        site = prm['site']
        sd = out / 'sites' / site
        rj = sd / 'radar.json'
        if not rj.exists():
            continue
        sc = site_scene(load_site(site))
        v = np.load(d / 'volumes.npz')
        h = float(v['step'])
        keep = v['z'] >= sc['extent']['z'][0] + h / 2
        run = f"{d.name} · {s['manifest']['date']} · {s['manifest']['commit']}"
        name = d.name[len('lab_'):]
        vols, entries = [], []
        for case in s['cases']:
            short = 'real' if case == 'real' else 'twin'
            for P in prm['supports']:
                V = v[f'{case}_p{P}'].astype(np.float32)[:, :, keep]
                u8 = np.where(np.isfinite(V), 1 + np.round(254 * np.clip(np.nan_to_num(V, nan=0.0), 0, 1)), 0).astype(np.uint8)
                vid = f'radar-lab-{name}-{short}-p{P}'
                vols.append(_write_volume(sd, vid, V, u8, [float(v['x'][0]), float(v['y'][0])], h, {
                    'label': f"Satellite · lab run {name}, {'the real image' if short == 'real' else 'a motionless copy'}, support {P}",
                    'method': s['finding'], 'quantity': 'conditional adjusted R2 at each nominal depth below the surface',
                    'units': '0 to 1', 'range': [0.0, 1.0], 'run': run, 'tint': 'gated',
                    'caption': 'fit scores, unsmoothed; empty where nothing passed'}, ztop=float(v['z'][keep][0])))
                entries.append({'id': vid, 'case': short, 'support': int(P)})
        _merge_volumes(sd, vols, drop=lambda q, name=name: q['id'].startswith(f'radar-lab-{name}-'))
        cx, cy = prm['centre_m']
        stamp = re.search(r'_(\d{4})\d{4}T', prm.get('product', ''))
        by_site.setdefault(site, []).append({
            'name': name, 'title': s['manifest']['title'], 'volumes': entries,
            'pass': stamp.group(1) if stamp else 'synthetic', 'lines': prm.get('lines', 'ew'),
            'focus': [cx, cy, float(np.median(v['z_surface']))], 'radius_m': 620.0, 'note': s['finding'], 'run': run})
    for site, labs in by_site.items():
        rj = out / 'sites' / site / 'radar.json'
        radar = json.loads(rj.read_text())
        radar['lab'] = labs
        rj.write_text(json.dumps(radar, separators=(',', ':')))
    return [{'site': site, 'volumes': [e['id'] for lab in labs for e in lab['volumes']]} for site, labs in by_site.items()]


def export_survives(out: Path = DATA) -> dict | None:
    """What survives a change of pass or of lines (P2-21): across Khafre, the lesser of two runs' fit scores wherever both
    scored, for the two layouts on the 2022 pass and for the two passes, added to the lab's runs as choices of pass and
    lines ('both')."""
    rid = 'p2_21_what_survives'
    vpath = RESULTS / rid / 'volumes.npz'
    rj = out / 'sites' / 'giza' / 'radar.json'
    if not (vpath.exists() and rj.exists()):
        return None
    from ..compose import load_site
    from .sites import scene as site_scene
    s = load(rid)
    sc = site_scene(load_site('giza'))
    v = np.load(vpath)
    h = float(v['step'])
    d = out / 'sites' / 'giza'
    run = f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}"
    radar = json.loads(rj.read_text())
    labs = [x for x in radar.get('lab', []) if not x['name'].startswith('agree_')]
    vols = []
    for pair, name, title, pas, lines in (
            ('ew22_ns22', 'agree_lines', 'Where the east-west and north-south layouts agree, across Khafre, from the 2022 image', '2022', 'both'),
            ('ew22_ew25', 'agree_pass', 'Where the 2022 and 2025 passes agree, across Khafre, east-west lines', 'both', 'ew')):
        V = v[f'{pair}_agree'].astype(np.float32)
        z = v[f'{pair}_z']
        keep = z >= sc['extent']['z'][0] + h / 2
        V = V[:, :, keep]
        u8 = np.where(np.isfinite(V), 1 + np.round(254 * np.clip(np.nan_to_num(V, nan=0.0), 0, 1)), 0).astype(np.uint8)
        vid = f"radar-agree-{pair.replace('_', '-')}"
        vols.append(_write_volume(d, vid, V, u8, [float(v[f'{pair}_x'][0]), float(v[f'{pair}_y'][0])], h, {
            'label': f'Satellite · {title}',
            'method': 'The lesser of two runs\' fit scores wherever both scored (P2-21), support 1',
            'quantity': 'conditional adjusted R2 at each nominal depth below the surface, where both runs scored',
            'units': '0 to 1', 'range': [0.0, 1.0], 'run': run, 'tint': 'gated',
            'caption': 'the lesser of two fit scores; empty unless both runs scored'}, ztop=float(z[keep][0])))
        labs.append({'name': name, 'title': title, 'pass': pas, 'lines': lines,
                     'volumes': [{'id': vid, 'case': 'real', 'support': 1}],
                     'focus': [0.0, 0.0, float(np.median(z[keep][:1]))], 'radius_m': 620.0, 'note': s['finding'], 'run': run})
    _merge_volumes(d, vols, drop=lambda q: q['id'].startswith('radar-agree-'))
    radar['lab'] = labs
    rj.write_text(json.dumps(radar, separators=(',', ':')))
    return {'site': 'giza', 'volumes': [x['id'] for x in vols]}


def export_radar(out: Path = DATA) -> list[dict]:
    return [r for r in (export_bench(out), export_bench_gated(out), export_claim(out), export_real(out), export_khufu_gated(out))
            if r] + export_lab(out) + [r for r in (export_survives(out), export_real(out, 'giza-deep')) if r]
