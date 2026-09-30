"""P2-34 · The gated reconstruction's depth scale, pass by pass: why its pictures stand at two levels on the 2025 passes.

    uv run python experiments/p2_34_gated_depth_scale.py

The gated reconstruction fits each position's shifts, window by window, to one of ten modes and draws the fit at the
depth where the steering phase turns that many times across a window of 50 pairs: one turn spans 2 pi / (50 dKz), dKz
the spacing of the pass's steering wavenumbers. The scale repeats every 2 pi / dKz, and a fit cannot tell a phase ramp
from its reverse, so each mode's depth has a mirror that far down, less the depth itself. Both numbers come from the
pass's geometry through dKz, not from the ground. Read here from every gated run on disk, from its own audit (its
steering wavenumbers) and its own volume (where each position's best score sits below the surface point above it):

- the metres one turn spans and the depth at which the scale repeats;
- how many positions that pass put their best score within ten turns of the surface (the modes themselves) and how many
  within ten turns of the repeat (their mirrors): the two levels the pictures stand at when the repeat is shallower than
  the 300 m the method draws.
"""
import json
from pathlib import Path

import numpy as np

from katabasis.runs import RESULTS, Run, load

RID = 'p2_34_gated_depth_scale'
DATA = Path(__file__).resolve().parents[1] / 'data'
RUNS = ['p2_16_khufu_volume', 'p2_17_plateau_control', 'lab_khafre', 'lab_khafre_ns', 'lab_khafre_2025',
        'lab_sacsay_walls', 'lab_sacsay_grid', 'lab_sacsay_fields']
WINDOW = 50                  # pairs per window in the frozen profile
MODES = 10


def scale(rid):
    """One turn and the repeat, from the run's steering wavenumbers (every line's, the median spacing)."""
    a = np.load(DATA / rid / 'real' / 'biondi_v1_8_audit.npz', allow_pickle=True)
    kz = np.sort(a['kz_rad_per_m'], axis=1)
    dk = float(np.median(np.diff(kz, axis=1)))
    return 2 * np.pi / (WINDOW * dk), 2 * np.pi / dk


def levels(rid, turn, repeat):
    """Where each passing position's best score sits below the surface, at a support of one."""
    v = np.load(RESULTS / rid / 'volumes.npz')
    V = v['real_p1'].astype(np.float32)
    ok = np.isfinite(V)
    passed = ok.any(axis=2)
    k = np.argmax(np.where(ok, V, -np.inf), axis=2)
    best = (v['z_surface'].astype(np.float32) - v['z'][k])[passed]
    reach = MODES * turn
    shallow = int(np.sum(best <= reach))
    mirror = int(np.sum(np.abs(best - repeat) <= reach))
    return {'positions': int(passed.sum()), 'shallow': shallow, 'mirror': mirror, 'elsewhere': int(passed.sum()) - shallow - mirror,
            'drawn_to_m': float(np.nanmax(v['z_surface']) - v['z'][-1]),
            'best_depth_counts_per_6m': np.histogram(best, bins=np.arange(0, 306, 6))[0].tolist()}


def main():
    have = [r for r in RUNS if (DATA / r / 'real' / 'biondi_v1_8_audit.npz').exists() and (RESULTS / r / 'volumes.npz').exists()]
    params = {'runs': have, 'window_pairs': WINDOW, 'modes': MODES,
              'run_commits': {r: load(r)['manifest']['commit'] for r in have}}
    with Run(RID, 'The gated reconstruction\'s depth scale, pass by pass', params) as run:
        out = {}
        for rid in have:
            turn, repeat = scale(rid)
            product = Path(json.loads((DATA / rid / 'profile.json').read_text())['input']['hdf5_path']).name   # what it read
            out[rid] = {'product': product, 'turn_m': turn, 'repeat_m': repeat, **levels(rid, turn, repeat)}
            r = out[rid]
            print(f"  {rid}: one turn {turn:.2f} m, repeat {repeat:.1f} m; of {r['positions']} positions {r['shallow']} near the "
                  f"surface, {r['mirror']} at the mirror", flush=True)
        by = lambda f: sorted({round(x['repeat_m']) for x in out.values() if f(x['product'])})
        p22 = [x for x in out.values() if '2022' in x['product']]
        p25 = [x for x in out.values() if '2025' in x['product']]
        finding = (
            f"The gated reconstruction's depth comes from the pass, not the ground. On the 2022 pass one turn of its fit "
            f"spans {min(x['turn_m'] for x in p22):.2f} m and the scale repeats every {min(x['repeat_m'] for x in p22):.0f} m, "
            f"past the 300 m it draws, so each position's fit stands near the surface "
            f"({sum(x['shallow'] for x in p22):,} of {sum(x['positions'] for x in p22):,} positions within ten turns of it). "
            f"On the 2025 passes, over Giza and over Sacsayhuamán, one turn spans {min(x['turn_m'] for x in p25):.2f} m and the "
            f"scale repeats every {min(x['repeat_m'] for x in p25):.0f} to {max(x['repeat_m'] for x in p25):.0f} m, so the same "
            f"fits stand twice, near the surface ({sum(x['shallow'] for x in p25):,} of {sum(x['positions'] for x in p25):,} "
            f"positions) and at their mirror near {min(x['repeat_m'] for x in p25):.0f} m down "
            f"({sum(x['mirror'] for x in p25):,}): one choice drawn at two levels.") if p22 and p25 else 'Runs missing.'
        run.save({'runs': out, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
