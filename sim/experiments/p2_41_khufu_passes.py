"""P2-41 · What survives a change of pass across Khufu: P2-21's test, on the Great Pyramid.

    uv run python experiments/p2_41_khufu_passes.py

P2-16 ran the gated reconstruction across Khufu from the 2022 pass, along 95 east-west lines 3 m apart, beside a
motionless copy of the image; the lab's command has run the same lines from the 2025 pass (lab_khufu_2025), beside a
motionless copy of its own. Something inside or under Khufu would not move with the pass that looked at it; what the
processing makes of each image's texture would. Measured as P2-21 measures Khafre, with its own functions: the positions
both passes pass against chance (a b / n), the share of those whose best depths lie within half the fit's cycle of each
other against shuffled pairings, and each pass against its own motionless copy, which shares its texture and has nothing
under it. Written for the lab: where the passes agree, the lesser of their two fit scores wherever both scored, at each
support.
"""
import importlib.util
from pathlib import Path

import numpy as np

from katabasis.runs import Run, load

RID = 'p2_41_khufu_passes'
RUNS = {'2022': 'p2_16_khufu_volume', '2025': 'lab_khufu_2025'}
HERE = Path(__file__).resolve().parent


def p2_21():
    """P2-21's own functions (loading, the common lattice, the comparison), so both tests are the same test."""
    spec = importlib.util.spec_from_file_location('p2_21', HERE / 'p2_21_what_survives.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    t = p2_21()
    params = {'runs': RUNS, 'run_commits': {k: load(r)['manifest']['commit'] for k, r in RUNS.items()},
              'cycle_m': t.CYCLE_M, 'shuffles': t.SHUFFLES, 'test': 'p2_21_what_survives.compare'}
    with Run(RID, 'What survives a change of pass across Khufu', params) as run:
        rng = np.random.default_rng(41)
        a, b = t.load_run(RUNS['2022']), t.load_run(RUNS['2025'])
        xs, ys, zs, (ia, ib) = t.common(a, b)
        zsurf = t.crop(a['z_surface'][..., None], (ia[0], ia[1], [0]))[..., 0]
        supports = sorted(int(k.split('_p')[1]) for k in a if k.startswith('real_p') and f'real_p{k.split("_p")[1]}' in b)
        res, agree = {}, {}
        for P in supports:
            A, B = t.crop(a[f'real_p{P}'], ia), t.crop(b[f'real_p{P}'], ib)
            res[f'p{P}'] = {
                'passes': t.compare(A, B, zs, zsurf, rng),
                'copies': t.compare(t.crop(a[f'twin101_p{P}'], ia), t.crop(b[f'twin101_p{P}'], ib), zs, zsurf, rng),
                '2022_vs_copy': t.compare(A, t.crop(a[f'twin101_p{P}'], ia), zs, zsurf, rng),
                '2025_vs_copy': t.compare(B, t.crop(b[f'twin101_p{P}'], ib), zs, zsurf, rng)}
            agree[f'agree_p{P}'] = np.where(np.isfinite(A) & np.isfinite(B), np.minimum(A, B), np.nan).astype(np.float16)
            print(f"  support {P}: {res[f'p{P}']['passes']}", flush=True)
        np.savez_compressed(run.dir / 'volumes.npz', **agree, x=xs, y=ys, z=zs, z_surface=zsurf, step=float(a['step']))
        r1 = res['p1']
        P, C, S = r1['passes'], r1['copies'], r1['2022_vs_copy']
        pct = lambda v: f'{v * 100:.0f}%'
        finding = (
            f"Across Khufu, of {P['positions']:,} positions on the same lines the 2022 pass passes {P['passing'][0]:,} at a support of one "
            f"and the 2025 pass {P['passing'][1]:,}; {P['both']:,} pass in both, {P['ratio']:.2f} times chance, and at "
            f"{pct(P['same_depth_share'])} of those the best depths agree within half a cycle, against {pct(P['same_depth_chance'])} when "
            f"paired at random. The two motionless copies, with nothing moving and nothing inside, share positions at "
            f"{C['ratio']:.2f} times chance; the 2022 image and its own copy at {S['ratio']:.2f} times. "
            + ('Across the two passes they agree beyond chance, but no more than the images\' texture does on its own.'
               if P['ratio'] <= max(C['ratio'], S['ratio']) else
               'Across the two passes they agree more than the images\' texture does on its own.'))
        run.save({'supports': supports, 'results': res, 'agree_positions': {P: int(np.isfinite(agree[f'agree_p{P}'].astype(np.float32)).any(axis=2).sum()) for P in supports},
                  'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
