"""Site files (sites/<id>/site.json): loading and validation."""
from __future__ import annotations

import copy
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import shapes
from .materials import Material, load_materials
from .terrain import Terrain, build_terrain

ROOT = Path(__file__).resolve().parents[3]
SITES_DIR = ROOT / 'sites'

KINDS = ('test', 'real')
FEATURE_KINDS = ('chamber', 'passage', 'shaft', 'void', 'column', 'tomb', 'corridor')
FEATURE_STATUS = ('truth', 'surveyed', 'claimed', 'representative')
FILLS = ('air', 'water', 'rubble', 'dry-sand-giza', 'wet-sand')


@dataclass
class Stratum:
    material: str
    top: object                  # float, 'surface', or {'plane': {...}}
    name: str = ''
    heterogeneity: dict | None = None

    def top_height(self, x: np.ndarray, y: np.ndarray, terrain: Terrain) -> np.ndarray:
        if self.top == 'surface':
            return terrain.height(x, y)
        if isinstance(self.top, (int, float)):
            return np.full(np.broadcast(x, y).shape, float(self.top))
        p = self.top['plane']
        az = math.radians(p.get('dip_azimuth_deg', 0.0))
        slope = math.tan(math.radians(p.get('dip_deg', 0.0)))
        return p['z0'] - slope * (np.asarray(x) * math.sin(az) + np.asarray(y) * math.cos(az))


@dataclass
class Site:
    id: str
    name: str
    kind: str
    summary: str
    extent: dict
    terrain: Terrain
    cover: list[dict]
    strata: list[Stratum]
    water_table: float | None
    structures: list[dict]
    features: list[dict]
    materials: dict[str, Material]
    frame: dict = field(default_factory=dict)
    sources: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    raw: dict = field(default_factory=dict, repr=False)
    directory: Path | None = None

    def used_materials(self) -> list[str]:
        ids = [c['material'] for c in self.cover] + [s.material for s in self.strata]
        ids += [s['material'] for s in self.structures] + [f['fill'] for f in self.features]
        if self.water_table is not None:
            ids += [self.materials[m].saturated_as for m in ids if self.materials[m].saturated_as]
        return list(dict.fromkeys(['air', *ids]))

    def feature(self, fid: str) -> dict:
        for f in self.features:
            if f['id'] == fid:
                return f
        raise KeyError(fid)


def _check(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def parse_site(raw: dict, directory: Path | None = None, materials: dict | None = None) -> Site:
    materials = materials or load_materials()
    sid = raw['id']
    _check(raw.get('kind') in KINDS, f'{sid}: kind must be one of {KINDS}')
    ext = raw['extent']
    for a in 'xyz':
        _check(len(ext[a]) == 2 and ext[a][0] < ext[a][1], f'{sid}: extent.{a} must be [lo, hi]')
    terrain = build_terrain(raw.get('terrain', {'type': 'flat'}), ext, directory or SITES_DIR / sid)

    cover = raw.get('cover', [])
    for c in cover:
        _check(c['material'] in materials, f'{sid}: unknown cover material {c["material"]!r}')
        _check(c['thickness'] > 0, f'{sid}: cover thickness must be positive')

    strata = [Stratum(s['material'], s['top'], s.get('name', ''), s.get('heterogeneity')) for s in raw['strata']]
    _check(len(strata) > 0, f'{sid}: at least one stratum')
    _check(strata[0].top == 'surface', f'{sid}: the first stratum starts at the surface')
    for s in strata:
        _check(s.material in materials, f'{sid}: unknown stratum material {s.material!r}')
        _check(not materials[s.material].is_fluid, f'{sid}: a stratum cannot be a fluid ({s.material})')

    ids = set()
    for group, allowed in (('structures', None), ('features', FEATURE_KINDS)):
        for f in raw.get(group, []):
            _check(f['id'] not in ids, f'{sid}: duplicate id {f["id"]!r}')
            ids.add(f['id'])
            shapes.validate(f['shape'])
            if group == 'features':
                _check(f.get('kind') in allowed, f'{sid}/{f["id"]}: kind must be one of {allowed}')
                _check(f.get('status') in FEATURE_STATUS, f'{sid}/{f["id"]}: status must be one of {FEATURE_STATUS}')
                f.setdefault('fill', 'air')
                _check(f['fill'] in materials, f'{sid}/{f["id"]}: unknown fill {f["fill"]!r}')
                _check(bool(f.get('source')), f'{sid}/{f["id"]}: every feature needs a source')
            else:
                _check(f['material'] in materials, f'{sid}/{f["id"]}: unknown material')

    wt = raw.get('water_table')
    return Site(sid, raw['name'], raw['kind'], raw.get('summary', ''), ext, terrain, cover, strata,
                None if wt is None else float(wt['z']), raw.get('structures', []), raw.get('features', []),
                materials, raw.get('frame', {}), raw.get('sources', {}), raw.get('notes', []), raw, directory)


def resolve(raw: dict, sites_dir: Path = SITES_DIR) -> tuple[dict, Path | None]:
    """A composition that extends another: the base's, with this one's keys over it and its `features_add`,
    `sources_add` and `notes_add` added to the base's. Returns the merged composition and the directory its terrain
    and other files come from (the base's), or None when it extends nothing."""
    base_id = raw.get('extends')
    if not base_id:
        return raw, None
    base, base_dir = resolve(json.loads((sites_dir / base_id / 'site.json').read_text()), sites_dir)
    out = copy.deepcopy(base)
    for k, v in raw.items():
        if k == 'extends':
            continue
        if k == 'features_add':
            out['features'] = out.get('features', []) + v
        elif k == 'sources_add':
            out['sources'] = {**out.get('sources', {}), **v}
        elif k == 'notes_add':
            out['notes'] = out.get('notes', []) + v
        else:
            out[k] = v
    return out, base_dir or sites_dir / base_id


def load_site(sid: str, sites_dir: Path = SITES_DIR) -> Site:
    d = sites_dir / sid
    raw, base_dir = resolve(json.loads((d / 'site.json').read_text()), sites_dir)
    return parse_site(raw, base_dir or d)


def list_sites(sites_dir: Path = SITES_DIR) -> list[str]:
    return sorted(p.parent.name for p in sites_dir.glob('*/site.json'))
