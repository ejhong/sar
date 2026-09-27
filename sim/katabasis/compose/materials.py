"""The materials library (sites/materials.json)."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MATERIALS_FILE = ROOT / 'sites' / 'materials.json'

STATUSES = {'sourced', 'derived', 'assumed', 'definition'}


@dataclass(frozen=True)
class Property:
    value: float
    range: tuple[float, float] | None
    status: str
    source: str | None = None
    quote: str | None = None
    note: str | None = None

    @classmethod
    def parse(cls, d: dict, sources: dict) -> 'Property':
        status = d.get('status', 'sourced' if 'source' in d else None)
        if status not in STATUSES:
            raise ValueError(f'property {d} has no source and no status')
        if status == 'sourced' and not d.get('quote'):
            raise ValueError(f'sourced property {d} needs its quoted line')
        src = d.get('source')
        if src is not None and src not in sources:
            raise ValueError(f'unknown source {src!r}')
        rng = tuple(d['range']) if 'range' in d else None
        if rng is not None and not (rng[0] <= d['value'] <= rng[1]):
            raise ValueError(f'value {d["value"]} outside its range {rng}')
        return cls(float(d['value']), rng, status, src, d.get('quote'), d.get('note'))


@dataclass(frozen=True)
class Material:
    id: str
    name: str
    family: str
    colour: str
    vp: Property
    vs: Property
    rho: Property
    saturated_as: str | None = None
    extra: dict = field(default_factory=dict, compare=False)

    @property
    def is_fluid(self) -> bool:
        return self.vs.value == 0.0

    def summary(self) -> dict:
        def p(x: Property):
            return {'value': x.value, 'range': list(x.range) if x.range else None, 'status': x.status,
                    'source': x.source}
        return {'id': self.id, 'name': self.name, 'family': self.family, 'colour': self.colour,
                'vp': p(self.vp), 'vs': p(self.vs), 'rho': p(self.rho)}


@lru_cache(maxsize=4)
def _load(path: str) -> tuple[dict, dict]:
    raw = json.loads(Path(path).read_text())
    sources = raw['sources']
    out = {}
    for mid, d in raw['materials'].items():
        out[mid] = Material(mid, d['name'], d['family'], d['colour'],
                            Property.parse(d['vp'], sources), Property.parse(d['vs'], sources),
                            Property.parse(d['rho'], sources), d.get('saturated_as'))
    for m in out.values():
        if m.saturated_as and m.saturated_as not in out:
            raise ValueError(f'{m.id}: saturated_as {m.saturated_as!r} is not a material')
        if not m.is_fluid and m.vs.value >= m.vp.value / 2 ** 0.5:
            raise ValueError(f'{m.id}: Vs must be below Vp/sqrt(2) for a positive Lame lambda')
    return out, sources


def load_materials(path: Path | str = MATERIALS_FILE) -> dict[str, Material]:
    return _load(str(path))[0]


def load_sources(path: Path | str = MATERIALS_FILE) -> dict:
    return _load(str(path))[1]
