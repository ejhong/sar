"""Recovered volumes and survey layouts for the viewer.

Each published run that makes tomograms registers them with its site:
web/public/data/sites/<site>/volumes.json lists them, vol-<id>.u8 holds
each grid (x fastest, then y, then depth), surveys.json the stations.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..runs import RESULTS
from .sites import DATA

METHOD = {
    'crosshole': ('Travel-time tomography, crosshole', 'First arrivals between four boreholes, inverted for P-wave speed'),
    'surface': ('Travel-time tomography, surface', 'First arrivals from a hammer grid to a geophone grid'),
}


def export_p1_02(out: Path = DATA) -> list[dict]:
    rid = 'p1_02_traveltime'
    summ = RESULTS / rid / 'summary.json'
    if not summ.exists():
        return []
    s = json.loads(summ.read_text())
    arrays = np.load(RESULTS / rid / 'volumes.npz')
    site = s['manifest']['params']['site']
    d = out / 'sites' / site
    d.mkdir(parents=True, exist_ok=True)
    vols, surveys = [], []
    for v in s['volumes']:
        name = v['name']
        vid = f'tt-{name}'
        data = arrays[name]
        (d / f'vol-{vid}.u8').write_bytes(np.ascontiguousarray(data.transpose(2, 1, 0)).tobytes())
        sc = s['surveys'][name]['score']
        label, method = METHOD[name]
        vols.append({'id': vid, 'label': label, 'method': method, 'survey': name, 'status': 'tomogram',
                     'quantity': v['quantity'], 'units': v['units'], 'file': f'vol-{vid}.u8',
                     'shape': v['shape'], 'origin': v['origin'], 'spacing': v['spacing'], 'range': v['range'],
                     'run': f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}",
                     'caption': f"chamber recovered at {100 * sc['chamber_mean_rel']:+.1f}% · brightest = 12% slower than the rock or more"})
        sv = s['surveys'][name]
        surveys.append({'id': name, 'label': name, 'description': sv['description'],
                        'sources': sv['stations']['sources'], 'receivers': sv['stations']['receivers'],
                        'boreholes': sv.get('boreholes', [])})
    (d / 'volumes.json').write_text(json.dumps(vols, indent=1))
    (d / 'surveys.json').write_text(json.dumps(surveys, separators=(',', ':')))
    return vols


def _signed_u8(frames: np.ndarray, vmax: float) -> np.ndarray:
    q = np.sign(frames) * np.sqrt(np.clip(np.abs(frames) / vmax, 0, 1))
    return np.round(128 + 127 * q).astype(np.uint8)


def export_p1_03(out: Path = DATA) -> list[dict]:
    rid = 'p1_03_scattered'
    summ = RESULTS / rid / 'summary.json'
    if not summ.exists():
        return []
    s = json.loads(summ.read_text())
    arrays = np.load(RESULTS / rid / 'volumes.npz')
    site = s['manifest']['params']['site']
    d = out / 'sites' / site
    vols = json.loads((d / 'volumes.json').read_text()) if (d / 'volumes.json').exists() else []
    vols = [v for v in vols if v['id'] != 'rtm-surface']
    q = arrays['rtm']
    (d / 'vol-rtm-surface.u8').write_bytes(np.ascontiguousarray(q.transpose(2, 1, 0)).tobytes())
    v = s['volume']
    sc = s['score']
    vols.append({'id': 'rtm-surface', 'label': 'Echo imaging, surface', 'method': 'Reverse-time migration of the scattered wave',
                 'survey': 'surface', 'status': 'tomogram', 'quantity': v['quantity'], 'units': v['units'],
                 'file': 'vol-rtm-surface.u8', 'shape': v['shape'], 'origin': v['origin'], 'spacing': v['spacing'],
                 'range': [0, 1], 'run': f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}",
                 'caption': f"best case, rock known exactly · chamber {sc['contrast']:.1f}× brighter than elsewhere"})
    (d / 'volumes.json').write_text(json.dumps(vols, indent=1))
    mv = s['movie']
    movie, echo = arrays['movie'], arrays['echo']
    vmax = float(np.percentile(np.abs(movie), 99.8))
    emax = float(np.percentile(np.abs(echo), 99.8))
    waves = []
    for wid, frames, m, label, caption, kind in (
            ('motion', movie, vmax, 'A hammer blow', "the ground's vertical motion from the central shot", 'motion'),
            ('echo', echo, emax, "The chamber's echo alone", f'with minus without the chamber, amplified {vmax / emax:.0f}×', 'echo')):
        (d / f'wave-{wid}.u8').write_bytes(np.ascontiguousarray(_signed_u8(frames, m).transpose(0, 2, 1)).tobytes())
        go = mv['grid_origin']
        waves.append({'id': wid, 'label': label, 'caption': caption, 'file': f'wave-{wid}.u8',
                      'shape': [int(frames.shape[0]), int(frames.shape[1]), int(frames.shape[2])],
                      'origin': [go[0], go[1]], 'spacing': mv['spacing'], 'z': 0.0, 'dt_ms': mv['dt_ms'], 'kind': kind,
                      'shot': mv['shot']})
    (d / 'wavefields.json').write_text(json.dumps(waves, indent=1))
    return vols


def export_p2_22(out: Path = DATA) -> list[dict]:
    """The one shaking (P2-22): the vibrator's steady vertical motion over the bench with its chamber, and the chamber's
    imprint on it (with minus without), three cycles each, added to the site's waves."""
    rid = 'p2_22_one_shaking'
    fpath = RESULTS / rid / 'fields.npz'
    if not fpath.exists():
        return []
    s = json.loads((RESULTS / rid / 'summary.json').read_text())
    a = np.load(fpath)
    x, y, f = a['x'], a['y'], float(a['f_hz'])
    up0, up1 = a['U_none'][..., 2], a['U_bench'][..., 2]
    n = 72                                                      # three cycles, 24 frames each
    ph = np.exp(2j * np.pi * 3 * np.arange(n) / n)
    movie = np.real(up1[None] * ph[:, None, None])
    imprint = np.real((up1 - up0)[None] * ph[:, None, None])
    vmax = float(np.percentile(np.abs(movie), 99.0))              # the vibrator's own spot would swamp the scale
    emax = float(np.percentile(np.abs(imprint), 99.8))
    d = out / 'sites' / s['manifest']['params']['site']
    path = d / 'wavefields.json'
    waves = [w for w in (json.loads(path.read_text()) if path.exists() else []) if not w['id'].startswith('vibrator')]
    src = [float(v) for v in a['source']]
    for wid, frames, m, label, caption, kind in (
            ('vibrator', movie, vmax, 'A vibrator, steady',
             f'{f:.2f} Hz, 30 m from the chamber: the one shaking geophones and the satellite both read', 'motion'),
            ('vibrator-imprint', imprint, emax, "The chamber's imprint on it",
             f'with minus without the chamber, amplified {vmax / emax:.0f}×', 'echo')):
        (d / f'wave-{wid}.u8').write_bytes(np.ascontiguousarray(_signed_u8(frames, m).transpose(0, 2, 1)).tobytes())
        waves.append({'id': wid, 'label': label, 'caption': caption, 'file': f'wave-{wid}.u8',
                      'shape': [int(frames.shape[0]), int(frames.shape[1]), int(frames.shape[2])],
                      'origin': [float(x[0]), float(y[0])], 'spacing': float(x[1] - x[0]), 'z': 0.0,
                      'dt_ms': 1e3 / f / 24, 'kind': kind, 'shot': [src[0], src[1], src[2]]})
    path.write_text(json.dumps(waves, indent=1))
    return [{'id': w['id']} for w in waves if w['id'].startswith('vibrator')]


def export_p1_04(out: Path = DATA) -> list[dict]:
    rid = 'p1_04_ambient'
    summ = RESULTS / rid / 'summary.json'
    if not summ.exists():
        return []
    s = json.loads(summ.read_text())
    arrays = np.load(RESULTS / rid / 'volumes.npz')
    site = s['manifest']['params']['site']
    d = out / 'sites' / site
    vols = json.loads((d / 'volumes.json').read_text()) if (d / 'volumes.json').exists() else []
    vols = [v for v in vols if v['id'] != 'ambient-surface']
    q = arrays['ambient']
    (d / 'vol-ambient-surface.u8').write_bytes(np.ascontiguousarray(q.transpose(2, 1, 0)).tobytes())
    v = s['volume']
    sc = s['score']
    vols.append({'id': 'ambient-surface', 'label': 'Ambient noise, passive', 'method': 'Noise correlations as virtual sources, echoes migrated',
                 'survey': 'ambient', 'status': 'tomogram', 'quantity': v['quantity'], 'units': v['units'],
                 'file': 'vol-ambient-surface.u8', 'shape': v['shape'], 'origin': v['origin'], 'spacing': v['spacing'],
                 'range': [0, 1], 'run': f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}",
                 'caption': f"best case, rock known, endless record · chamber {sc['image_contrast_ensemble']:.1f}× brighter than elsewhere"})
    (d / 'volumes.json').write_text(json.dumps(vols, indent=1))
    svp = d / 'surveys.json'
    surveys = [x for x in (json.loads(svp.read_text()) if svp.exists() else []) if x['id'] != 'ambient']
    st = s['stations']
    surveys.append({'id': 'ambient', 'label': 'ambient', 'description': 'noise sources scattered over the ground, the surface geophone grid listening',
                    'sources': st['noise_sources'], 'receivers': st['receivers']})
    svp.write_text(json.dumps(surveys, separators=(',', ':')))
    return vols


def export_p1_05(out: Path = DATA) -> list[dict]:
    rid = 'p1_05_fwi'
    summ = RESULTS / rid / 'summary.json'
    if not summ.exists():
        return []
    s = json.loads(summ.read_text())
    arrays = np.load(RESULTS / rid / 'volumes.npz')
    site = s['manifest']['params']['site']
    d = out / 'sites' / site
    vols = json.loads((d / 'volumes.json').read_text()) if (d / 'volumes.json').exists() else []
    vols = [v for v in vols if not v['id'].startswith('fwi-')]
    v = s['volume']
    for name, label, quantity in (('vs', 'Full-waveform inversion, S speed', 'relative S-wave speed'),
                                  ('vp', 'Full-waveform inversion, P speed', 'relative P-wave speed')):
        vid = f'fwi-{name}'
        (d / f'vol-{vid}.u8').write_bytes(np.ascontiguousarray(arrays[name].transpose(2, 1, 0)).tobytes())
        sc = s['score'][name]
        vols.append({'id': vid, 'label': label, 'method': "Whole records fitted from a practitioner's start, 60 Hz then 120 Hz",
                     'survey': 'surface', 'status': 'tomogram', 'quantity': quantity, 'units': v['units'],
                     'file': f'vol-{vid}.u8', 'shape': v['shape'], 'origin': v['origin'], 'spacing': v['spacing'],
                     'range': v['range'], 'run': f"{rid} · {s['manifest']['date']} · {s['manifest']['commit']}",
                     'caption': f"chamber recovered at {100 * sc['chamber_mean_rel']:+.1f}% · brightest = 12% slower than the rock or more"})
    (d / 'volumes.json').write_text(json.dumps(vols, indent=1))
    return vols
