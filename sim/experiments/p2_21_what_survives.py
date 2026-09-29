"""P2-21 · What survives a change of pass or of lines: the gated reconstruction across Khafre, three ways.

    uv run python experiments/p2_21_what_survives.py

The lab's command has run the gated reconstruction across the same square over Khafre three ways: from the 2022 pass
with its lines laid east-west (lab_khafre) and north-south (lab_khafre_ns), and from the 2025 pass, east-west
(lab_khafre_2025), each beside a motionless copy of its image. Something under Khafre would not move with the way the
lines are laid or with the pass that looked at it; what the processing makes of an image's texture would.

Stated before the run: if the reconstruction images something under Khafre, the positions it passes, and the depths it
scores there, will agree between the three runs more than chance allows, and more than each run agrees with its own
motionless copy; if its picture comes from each image and each layout, the overlap will sit near chance. Two runs on
the same image share its texture, so their agreement beyond chance is expected on that account alone; two passes do
not, which makes the pass comparison the sharper test.

Chance, for positions: a run passing a of the n positions both runs cover, and another b, share a b / n by chance.
For depths, at positions both pass: the share of such positions whose best depths lie within half the fit's cycle of
each other, against the share when one run's depths are paired with the other's at random (1,000 shuffles). Written:
the agreement volumes (the lesser of two scores wherever both runs scored) for the viewer, and the numbers.
"""
import json
from pathlib import Path

import numpy as np

from katabasis.runs import RESULTS, Run, load

RID = 'p2_21_what_survives'
RUNS = {'ew22': 'lab_khafre', 'ns22': 'lab_khafre_ns', 'ew25': 'lab_khafre_2025'}
PAIRS = [('ew22', 'ns22', 'lines'), ('ew22', 'ew25', 'pass'), ('ns22', 'ew25', 'both')]
CYCLE_M = 6.17                      # the depth one cycle of the gates' first mode spans (P2-13): half of it is "the same depth"
SHUFFLES = 1000


def load_run(rid):
    v = np.load(RESULTS / rid / 'volumes.npz')
    return {k: v[k].astype(np.float32) for k in v.files}


def common(a, b):
    """Index maps of two runs' volumes onto the positions both cover (the same 3 m lattice) and onto heights paired
    within half a voxel (each run hangs its grid from its own highest surface point)."""
    xs = np.intersect1d(np.round(a['x'], 3), np.round(b['x'], 3))
    ys = np.intersect1d(np.round(a['y'], 3), np.round(b['y'], 3))
    near = np.abs(a['z'][:, None] - b['z'][None, :]).argmin(axis=1)
    ok = np.abs(a['z'] - b['z'][near]) <= 0.5 * float(a['step'])
    ka, kb = np.flatnonzero(ok), near[ok]
    ix = lambda r, vals, key: np.searchsorted(np.round(r[key], 3), vals)
    return xs, ys, a['z'][ka], [(ix(a, xs, 'x'), ix(a, ys, 'y'), ka), (ix(b, xs, 'x'), ix(b, ys, 'y'), kb)]


def crop(vol, idx):
    i, j, k = idx
    return vol[np.ix_(i, j, k)]


def best_depth(V, zs, z_surface):
    """At each position that passed, the depth below its surface of its highest score; NaN where nothing passed."""
    filled = np.where(np.isfinite(V), V, -np.inf)
    k = np.argmax(filled, axis=2)
    passed = np.isfinite(V).any(axis=2)
    return np.where(passed, z_surface - zs[k], np.nan)


def compare(A, B, zs, zsurf, rng):
    pa, pb = np.isfinite(A).any(axis=2), np.isfinite(B).any(axis=2)
    n = pa.size
    both = pa & pb
    expected = pa.sum() * pb.sum() / n
    da, db = best_depth(A, zs, zsurf), best_depth(B, zs, zsurf)
    same = np.abs(da[both] - db[both]) <= CYCLE_M / 2
    null = [np.mean(np.abs(da[both] - rng.permutation(db[both])) <= CYCLE_M / 2) for _ in range(SHUFFLES)] if both.sum() > 1 else [np.nan]
    return {'positions': int(n), 'passing': [int(pa.sum()), int(pb.sum())], 'both': int(both.sum()),
            'chance_both': float(expected), 'ratio': float(both.sum() / expected) if expected else None,
            'same_depth_share': float(same.mean()) if same.size else None,
            'same_depth_chance': float(np.nanmean(null)), 'same_depth_chance_p95': float(np.nanpercentile(null, 95))}


def main():
    have = {k: (RESULTS / r / 'volumes.npz').exists() for k, r in RUNS.items()}
    if not all(have.values()):
        raise SystemExit(f'waiting for the lab runs: {have}')
    records = {k: load(r)['manifest'] for k, r in RUNS.items()}
    params = {'runs': RUNS, 'run_commits': {k: m['commit'] for k, m in records.items()}, 'pairs': PAIRS,
              'cycle_m': CYCLE_M, 'shuffles': SHUFFLES, 'support': 1}
    with Run(RID, 'What survives a change of pass or of lines', params) as run:
        rng = np.random.default_rng(21)
        vols = {k: load_run(r) for k, r in RUNS.items()}
        results, agreement, grids = {}, {}, {}
        for a, b, what in PAIRS:
            xs, ys, zs, (ia, ib) = common(vols[a], vols[b])
            zsurf = crop(vols[a]['z_surface'][..., None], (ia[0], ia[1], [0]))[..., 0]
            res = {}
            for case in ('real', 'twin101'):
                A, B = crop(vols[a][f'{case}_p1'], ia), crop(vols[b][f'{case}_p1'], ib)
                res[case] = compare(A, B, zs, zsurf, rng)
                if case == 'real':
                    agreement[f'{a}_{b}'] = np.where(np.isfinite(A) & np.isfinite(B), np.minimum(A, B), np.nan)
                    grids[f'{a}_{b}'] = {'x': xs, 'y': ys, 'z': zs}
            # each run against its own motionless copy, on the same positions
            for k, idx in ((a, ia), (b, ib)):
                res[f'{k}_vs_copy'] = compare(crop(vols[k]['real_p1'], idx), crop(vols[k]['twin101_p1'], idx), zs, zsurf, rng)
            results[f'{a}_{b}'] = dict(res, what=what)
            print(f"  {a} vs {b} ({what}): {res['real']}", flush=True)

        np.savez_compressed(Path(run.dir) / 'volumes.npz',
                            **{f'{k}_agree': v.astype(np.float16) for k, v in agreement.items()},
                            **{f'{k}_{ax}': g[ax] for k, g in grids.items() for ax in ('x', 'y', 'z')}, step=3.0)

        L, P = results['ew22_ns22']['real'], results['ew22_ew25']['real']
        C = results['ew22_ew25']['ew22_vs_copy']                  # an image and its motionless copy: nothing under it
        above = lambda r: r['ratio'] is not None and r['ratio'] > 1.5 and r['same_depth_share'] is not None \
            and r['same_depth_share'] > r['same_depth_chance_p95']
        beyond_copy = lambda r: r['ratio'] > C['ratio'] and r['same_depth_share'] > C['same_depth_share']
        pct = lambda x: f"{x * 100:.0f}%"
        finding = (
            f"Across Khafre, the gated reconstruction's lines laid north-south instead of east-west pass {L['passing'][1]:,} "
            f"positions against {L['passing'][0]:,}; {L['both']:,} pass both ways, {L['ratio']:.2f} times what chance gives, "
            f"and at {pct(L['same_depth_share'])} of those the best depths agree within half a cycle, against "
            f"{pct(L['same_depth_chance'])} by chance. The 2025 pass passes {P['passing'][1]:,}; {P['both']:,} of its positions "
            f"pass in 2022 too, {P['ratio']:.2f} times chance, their best depths agreeing at {pct(P['same_depth_share'])} "
            f"against {pct(P['same_depth_chance'])}. For scale, the 2022 image and its own motionless copy, which keeps only its "
            f"smoothed brightness and has nothing under it, share {C['ratio']:.2f} times chance and agree in depth at "
            f"{pct(C['same_depth_share'])}. "
            + ('Positions and depths survive a change of pass beyond chance and beyond what an image shares with its copy.'
               if above(P) and beyond_copy(P) else
               'Across the two passes they agree beyond chance, but no more than an image agrees with its own motionless copy.'
               if above(P) else
               'What passes, and at what depth, does not survive a change of pass beyond chance.'))
        run.save({'results': results, 'finding': finding,
                  'grids': {k: {ax: [float(g[ax][0]), float(g[ax][-1]), int(len(g[ax]))] for ax in ('x', 'y', 'z')} for k, g in grids.items()}})
        print(finding)


if __name__ == '__main__':
    main()
