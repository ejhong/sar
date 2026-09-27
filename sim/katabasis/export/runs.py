"""Published runs: sim/results/<id>/summary.json -> web/public/data/runs/<id>.json."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from ..runs import RESULTS
from .sites import DATA

PUBLISHED = ['p1_01_forward_validation', 'p1_02_traveltime', 'p1_03_scattered', 'p1_04_ambient', 'm1_01_ambient_levels']


def export_runs(out: Path = DATA) -> list[str]:
    d = out / 'runs'
    d.mkdir(parents=True, exist_ok=True)
    done = []
    for rid in PUBLISHED:
        src = RESULTS / rid / 'summary.json'
        if src.exists():
            shutil.copyfile(src, d / f'{rid}.json')
            done.append(rid)
    (d / 'index.json').write_text(json.dumps(
        {rid: json.loads((d / f'{rid}.json').read_text())['manifest'] for rid in done}, indent=1))
    return done
