"""Site scenes for the underworld viewer: web/public/data/sites/."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..compose import list_sites, load_site
from ..compose.materials import load_sources
from ..compose.site import Site

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'web' / 'public' / 'data'
SURFACE_N = 33           # samples per side for stratum surfaces


def _grid_block(x0: float, y0: float, dx: float, z: np.ndarray, digits: int = 2) -> dict:
    return {'x0': round(float(x0), 3), 'y0': round(float(y0), 3), 'dx': round(float(dx), 4),
            'nx': int(z.shape[0]), 'ny': int(z.shape[1]),
            'z': [round(float(v), digits) for v in z.ravel(order='F')]}   # x fastest


def scene(site: Site) -> dict:
    e = site.extent
    xs = np.linspace(e['x'][0], e['x'][1], SURFACE_N)
    ys = np.linspace(e['y'][0], e['y'][1], SURFACE_N)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    dx = xs[1] - xs[0]
    t = site.terrain
    if t.kind == 'flat':
        terrain = {'kind': 'flat', 'z': t.flat_z}
    else:
        terrain = {'kind': 'grid', **_grid_block(t.x[0], t.y[0], t.x[1] - t.x[0], t.z),
                   'source': (t.meta or {}).get('source', 'generated relief')}
    surface = t.height(X, Y)
    strata = []
    for s in site.strata:
        top = np.minimum(s.top_height(X, Y, t), surface)
        strata.append({'name': s.name, 'material': s.material,
                       'top': 'surface' if s.top == 'surface' else _grid_block(xs[0], ys[0], dx, top),
                       'heterogeneous': bool(s.heterogeneity)})
    mats = {m: site.materials[m].summary() for m in site.used_materials()}
    return {
        'id': site.id, 'name': site.name, 'kind': site.kind, 'summary': site.summary,
        'frame': site.frame, 'extent': e, 'terrain': terrain,
        'cover': site.cover, 'strata': strata, 'water_table': site.water_table,
        'materials': mats,
        'structures': site.structures,
        'features': site.features,
        'landmarks': site.raw.get('landmarks', []),
        'notes': site.notes, 'sources': site.sources,
        'volumes': [],
    }


def export_sites(out: Path = DATA) -> list[dict]:
    index = []
    for sid in list_sites():
        site = load_site(sid)
        sc = scene(site)
        d = out / 'sites' / sid
        d.mkdir(parents=True, exist_ok=True)
        vols = d / 'volumes.json'
        if vols.exists():                 # published runs register their volumes here
            sc['volumes'] = json.loads(vols.read_text())
        svs = d / 'surveys.json'
        if svs.exists():
            sc['surveys'] = json.loads(svs.read_text())
        wfs = d / 'wavefields.json'
        if wfs.exists():
            sc['wavefields'] = json.loads(wfs.read_text())
        rad = d / 'radar.json'
        if rad.exists():                  # the satellite's pass over the site (export/radar.py)
            sc['radar'] = json.loads(rad.read_text())
        (d / 'scene.json').write_text(json.dumps(sc, separators=(',', ':'), ensure_ascii=False))
        counts = {}
        for f in site.features:
            counts[f['status']] = counts.get(f['status'], 0) + 1
        radar = sum(v.get('status') == 'radar' for v in sc['volumes'])
        index.append({'id': sid, 'name': site.name, 'kind': site.kind, 'summary': site.summary,
                      'extent': site.extent, 'features': counts, 'volumes': len(sc['volumes']),
                      'instruments': {'geophones': len(sc['volumes']) - radar + len(sc.get('wavefields', [])),
                                      'satellite': radar}})
    order = {'test': 0, 'real': 1}
    index.sort(key=lambda s: (order[s['kind']], ['bench-void', 'bench-shafts', 'bench-khafre-claim'].index(s['id'])
                              if s['id'].startswith('bench') else 9, s['id']))
    (out / 'sites').mkdir(parents=True, exist_ok=True)
    (out / 'sites' / 'index.json').write_text(json.dumps({'sites': index}, indent=1, ensure_ascii=False))
    (out / 'materials.json').write_text(json.dumps(
        {'sources': load_sources(),
         'materials': {m.id: m.summary() | {'vp_quote': m.vp.quote, 'vs_note': m.vs.note}
                       for m in load_site(list_sites()[0]).materials.values()}}, indent=1, ensure_ascii=False))
    lib = Path(__file__).resolve().parents[3] / 'sites' / 'ambient.json'
    if lib.exists():                                  # the ambient library, verbatim
        (out / 'ambient.json').write_text(lib.read_text())
    return index
