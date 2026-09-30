"""Published runs: sim/results/<id>/summary.json -> web/public/data/runs/<id>.json, and an index of every run's record
with the command that repeats it: an experiment's from its docstring, a lab run's from its parameters."""
from __future__ import annotations

import ast
import json
import shutil
from pathlib import Path

from ..runs import RESULTS
from .sites import DATA

PUBLISHED = ['p1_01_forward_validation', 'p1_02_traveltime', 'p1_03_scattered', 'p1_04_ambient', 'm1_01_ambient_levels', 'p2_00_ambient_giza', 'p1_05_fwi',
             'p2_01_geometry', 'p2_03_what_a_dwell_sees', 'p2_04_chamber_imprint', 'p2_05_budget', 'p2_06_depth_is_frequency',
             'p2_07_whole_chain', 'p2_08_pillars_and_boxes', 'p2_09_virtual_sensors', 'p2_11_every_way_in', 'p2_12_real_pass',
             'p2_13_gates', 'p2_14_own_lines', 'p2_15_pair_response',
             'p2_16_khufu_volume', 'p2_17_plateau_control',
             'p2_18_echo_validation', 'p2_19_known_vibration', 'p2_20_bench_gated', 'p2_21_what_survives',
             'p2_22_one_shaking', 'p2_23_every_reader', 'p2_24_gated_on_the_shaking', 'p2_33_sacsayhuaman',
             'p2_34_gated_depth_scale', 'p2_35_rungs_and_rooms']

EXPERIMENTS = Path(__file__).resolve().parents[2] / 'experiments'


def command(rid: str, manifest: dict) -> str | None:
    """How to make the run again, from sim/."""
    prm = manifest.get('params', {})
    if rid.startswith('lab_') and prm.get('method') == 'gated':
        cx, cy = prm['centre_m']
        parts = ['uv run --with scikit-image python -m katabasis.lab gated', f"--name {rid[len('lab_'):]}", f"--site {prm['site']}",
                 f'--centre {cx:g},{cy:g}', f"--half-ew {prm['half_ew_m']:g}", f"--half-ns {prm['half_ns_m']:g}",
                 f"--step {prm['step_m']:g}", '--supports ' + ','.join(str(x) for x in prm['supports'])]
        if prm.get('twin_seed') is not None:
            parts.append(f"--twin {prm['twin_seed']}")
        if prm.get('lines', 'ew') != 'ew':
            parts.append(f"--lines {prm['lines']}")
        if prm.get('product') and prm['product'] != 'ICEYE_X13_SLC_SLED_868226_20220715T235744.h5':
            parts.append(f"--product <path to {prm['product']}>")
        return ' '.join(parts)
    src = EXPERIMENTS / f'{rid}.py'
    if not src.exists():
        return None
    doc = ast.get_docstring(ast.parse(src.read_text())) or ''
    return next((line.strip() for line in doc.splitlines() if line.strip().startswith('uv run ')), None)


def export_runs(out: Path = DATA) -> list[str]:
    d = out / 'runs'
    d.mkdir(parents=True, exist_ok=True)
    done = []
    labs = sorted(p.name for p in RESULTS.glob('lab_*') if (p / 'summary.json').exists())
    for rid in PUBLISHED + labs:
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
    index = {}
    for rid in done:
        rec = json.loads((d / f'{rid}.json').read_text())
        man = rec['manifest']
        index[rid] = dict(man, command=command(rid, man), **({'reused': True} if rec.get('reused_pipeline_outputs') else {}))
    (d / 'index.json').write_text(json.dumps(index, indent=1))
    return done
