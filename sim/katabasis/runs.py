"""Runs: every published result names the run that made it.

    with Run('p1_01_forward_validation', 'The solver against exact solutions', params) as run:
        ...
        run.save(summary)

writes sim/results/<id>/summary.json with a manifest (date, commit, runtime,
parameters, environment). Results are generated, never edited by hand.
"""
from __future__ import annotations

import json
import math
import os
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'sim' / 'results'


def _commit() -> str:
    try:
        sha = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(['git', 'diff', '--quiet', 'HEAD', '--', 'sim/katabasis', 'sim/sarsim', 'sim/experiments'],
                                cwd=ROOT) != 0
        return sha + ('+' if dirty else '')
    except Exception:
        return 'unknown'


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating,)):
        return _clean(float(o))
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, float):
        if not math.isfinite(o):                   # a vacuous bound (infinite) or an undefined score travels as null
            return None
        return round(o, 6) if abs(o) >= 1e-3 or o == 0 else float(f'{o:.4g}')
    return o


class Run:
    def __init__(self, rid: str, title: str, params: dict | None = None):
        self.id = rid
        self.title = title
        self.params = params or {}
        self.dir = RESULTS / rid
        self.t0 = 0.0

    def __enter__(self) -> 'Run':
        self.dir.mkdir(parents=True, exist_ok=True)
        self.t0 = time.time()
        self.started = datetime.now(timezone.utc)
        return self

    def __exit__(self, *exc):
        return False

    def manifest(self) -> dict:
        import numba
        return {'id': self.id, 'title': self.title,
                'date': self.started.strftime('%Y-%m-%d'), 'started_utc': self.started.isoformat(timespec='seconds'),
                'runtime_s': round(time.time() - self.t0, 1), 'commit': _commit(),
                'params': self.params,
                'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                                'numba': numba.__version__, 'threads': numba.get_num_threads(),
                                'cpus': os.cpu_count()}}

    def save(self, summary: dict, name: str = 'summary.json') -> Path:
        out = {'manifest': self.manifest(), **summary}
        p = self.dir / name
        p.write_text(json.dumps(_clean(out), indent=1, ensure_ascii=False, allow_nan=False))
        print(f'wrote {p.relative_to(ROOT)}')
        return p


def load(rid: str) -> dict:
    return json.loads((RESULTS / rid / 'summary.json').read_text())


def memo(rid: str, name: str, compute, source: str | Path, version: str | None = None):
    """A long step of a run kept on disk, so a run cut short resumes where it stopped: the value of compute() under
    sim/results/cache/<rid>/, keyed by `name` and by the experiment's source file, so any change to the code recomputes.
    With `version`, keyed by it instead: the caller declares that the step is unchanged while the version is, so a run
    can be extended without recomputing what it already has. Only intermediate values are kept; the summary is always
    written by the run itself."""
    import hashlib
    import pickle
    head = version.encode() if version is not None else Path(source).read_bytes()
    key = hashlib.sha1(head + name.encode()).hexdigest()[:16]
    path = RESULTS / 'cache' / rid / f'{key}.pkl'
    if path.is_file():
        with open(path, 'rb') as fh:
            return pickle.load(fh)
    value = compute()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    with open(tmp, 'wb') as fh:
        pickle.dump(value, fh, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)
    return value
