"""Feature generators for the test benches whose features are many and regular.

    uv run python -m katabasis.compose.benches

rewrites the `features` of each site that names a generator, so every site
file stays self-contained and reviewable.
"""
from __future__ import annotations

import json

import numpy as np

from .site import SITES_DIR

SHAFT_SOURCE = ('Representative: dimensions after Idu (G 7102 C) and Hetepheres I (G 7000 X) in '
                'sim/legacy/catalog/giza/surveys.json; positions and depths chosen to span the 5-30 m band.')


def shaft_field() -> list[dict]:
    rng = np.random.default_rng(7)
    depths = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
    rng.shuffle(depths)
    feats = []
    k = 0
    for row, y in enumerate((-30.0, -6.0, 18.0)):
        for col, x in enumerate((-36.0, -12.0, 12.0, 36.0)):
            d = depths[k]
            jx, jy = rng.uniform(-4, 4, 2)
            cx, cy = round(x + jx, 1), round(y + jy, 1)
            side = 1.72 if d < 15 else 1.3
            feats.append({'id': f'shaft-{k + 1:02d}', 'name': f'Shaft {k + 1} ({d:.0f} m)', 'kind': 'shaft',
                          'status': 'representative', 'fill': 'air',
                          'shape': {'type': 'box', 'centre': [cx, cy, round(-d / 2 + 0.5, 2)],
                                    'size': [side, side, d + 1.0]},
                          'source': SHAFT_SOURCE, 'group': 'shafts'})
            yaw = float(rng.choice([0.0, 90.0]))
            feats.append({'id': f'chamber-{k + 1:02d}', 'name': f'Chamber {k + 1}', 'kind': 'tomb',
                          'status': 'representative', 'fill': 'air',
                          'shape': {'type': 'box', 'centre': [cx, round(cy + 2.4, 1), round(-d + 1.0, 2)],
                                    'size': [2.72, 5.22, 1.95], 'yaw_deg': yaw},
                          'source': SHAFT_SOURCE, 'group': 'shafts'})
            k += 1
    return feats


def khafre_claim() -> list[dict]:
    feats = []
    src = 'Claimed (egyptianstreets2025: "eight vertical structures extending over 2,100 feet"); layout and 10 m diameter representative.'
    for i, (x, y) in enumerate([(x, y) for y in (-45.0, 45.0) for x in (-67.5, -22.5, 22.5, 67.5)]):
        feats.append({'id': f'well-{i + 1}', 'name': f'Claimed shaft {i + 1}', 'kind': 'shaft', 'status': 'claimed',
                      'fill': 'air', 'shape': {'type': 'cylinder', 'centre': [x, y, -320.0], 'radius': 5.0, 'height': 640.0},
                      'source': src, 'group': 'claimed shafts'})
    feats.append({'id': 'deep-structure', 'name': 'Claimed deep structure', 'kind': 'void', 'status': 'claimed',
                  'fill': 'air', 'shape': {'type': 'box', 'centre': [0.0, 0.0, -1220.0], 'size': [80.0, 80.0, 80.0]},
                  'source': 'Claimed (egyptianstreets2025: "additional unidentified structures potentially lying at depths of 4,000 feet"); size representative.',
                  'group': 'claimed deep'})
    return feats


GENERATORS = {'katabasis.compose.benches.shaft_field': shaft_field,
              'katabasis.compose.benches.khafre_claim': khafre_claim}


def main() -> None:
    for path in sorted(SITES_DIR.glob('*/site.json')):
        raw = json.loads(path.read_text())
        gen = raw.get('generator')
        if gen:
            raw['features'] = GENERATORS[gen]()
            path.write_text(json.dumps(raw, indent=1, ensure_ascii=False) + '\n')
            print(f'{path.parent.name}: {len(raw["features"])} features')


if __name__ == '__main__':
    main()
