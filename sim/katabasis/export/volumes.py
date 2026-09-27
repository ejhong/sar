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
