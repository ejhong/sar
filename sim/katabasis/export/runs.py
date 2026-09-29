"""Published runs: sim/results/<id>/summary.json -> web/public/data/runs/<id>.json."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from ..runs import RESULTS
from .sites import DATA

PUBLISHED = ['p1_01_forward_validation', 'p1_02_traveltime', 'p1_03_scattered', 'p1_04_ambient', 'm1_01_ambient_levels', 'p2_00_ambient_giza', 'p1_05_fwi',
             'p2_01_geometry', 'p2_03_what_a_dwell_sees', 'p2_04_chamber_imprint', 'p2_05_budget', 'p2_06_depth_is_frequency',
             'p2_07_whole_chain', 'p2_08_pillars_and_boxes', 'p2_09_virtual_sensors', 'p2_11_every_way_in', 'p2_12_real_pass',
             'p2_13_gates', 'p2_14_own_lines', 'p2_15_pair_response',
             'p2_16_khufu_volume']


def export_runs(out: Path = DATA) -> list[str]:
    d = out / 'runs'
    d.mkdir(parents=True, exist_ok=True)
    done = []
    for rid in PUBLISHED:
        src = RESULTS / rid / 'summary.json'
        if src.exists():
            shutil.copyfile(src, d / f'{rid}.json')
            figs = RESULTS / rid / 'figs'
            if figs.is_dir():                      # images a run renders for the site (P2-07)
                (d / rid).mkdir(exist_ok=True)
                for f in figs.glob('*.png'):
                    if f.stat().st_size > 250_000:     # photographic images (radar texture) travel as JPEG
                        from PIL import Image
                        Image.open(f).convert('RGB').save(d / rid / f'{f.stem}.jpg', quality=90, optimize=True)
                    else:
                        shutil.copyfile(f, d / rid / f.name)
            done.append(rid)
    (d / 'index.json').write_text(json.dumps(
        {rid: json.loads((d / f'{rid}.json').read_text())['manifest'] for rid in done}, indent=1))
    return done
