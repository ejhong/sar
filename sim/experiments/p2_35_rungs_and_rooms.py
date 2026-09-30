"""P2-35 · What would give the published pictures their rungs and rooms: a search over choices the paper does not state.

    uv run python experiments/p2_35_rungs_and_rooms.py

The published Khafre sections show wide columns crossed by clean horizontal rungs, with bright blocks low down. The
paper-style pipeline as the 2022 paper states it finds the same three ingredients in the real Khafre image (P2-08), but
there its bands stay beads that do not line up across a column's width, and its block at the repeat depth is faint
(about 1.4 times the mean). Keeping the motion the whole scene shares, instead of removing it, changes neither (checked
while designing this run). Here choices the paper does not state are varied on the same real crop (the first
investigation's Khafre patch of the 2025 pass, as P2-08 and P2-12 use it), and each picture is measured:

  patch   the registration patch, 32 px as the first investigation ran it, or 64 px, so neighbouring pixels share more of
          what they register;
  pool    each pixel's shifts averaged with its neighbours' before focusing, over 1, 3 x 3 or 7 x 7 targets (a coherent
          average: the shifts, not the pictures);
  bank    the paper's 50 half-band pairs 88 Hz apart, or 20 pairs 404 Hz apart (the derivative protocol v1.7's spacing).

Measured on the published-style section (the power averaged over a slab 31 lines thick along track, across the patch's
range, over just more than one period of the depth axis, relabelled so it repeats at 648 m):

  rungs   how alike the banding is across the width of the strongest column: the correlation, between neighbouring range
          positions, of their depth profiles once each profile's slow trend is removed, over the column's seven middle
          positions (24 m);
  room    the power at the repeat depth against the same column's mean away from the surface and the repeat.

Stated before the run: a setting makes rungs if neighbouring positions' banding correlates at 0.5 or more across the
strongest column, and makes a room if the repeat layer stands at least twice the column's mean. This searches for how the
published look could arise; it is not a test of the claim, and nothing in it depends on what lies below.
"""
import importlib
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d, uniform_filter

from katabasis.runs import Run
from sarsim.tomo import focus_paper

RID = 'p2_35_rungs_and_rooms'
LEGACY = Path(__file__).resolve().parents[1] / 'legacy' / 'experiments'
LEGACY_CACHE = Path(__file__).resolve().parents[2] / 'results' / 'cache'
PATCHES = [32, 64]
POOLS = [1, 3, 7]
BANKS = ['paper', 'v17']
SLAB = 15                    # lines either side of the middle along track: a slab 31 lines thick, as P2-08 renders it
NZ = 324
SHOW = 1.08                  # the section runs just past one repeat, as P2-08's does
COLUMN_HALF = 3              # the strongest column's middle seven range positions
TREND_M = 60.0               # the slow trend removed from each depth profile before comparing bands (relabelled metres)
RUNGS_AT, ROOM_AT = 0.5, 2.0


def trajectories(args):
    """The first investigation's patch run, with the registration patch and the bank set (cached under results/cache)."""
    patch, bank = args
    sys.path.insert(0, str(LEGACY))
    rc = importlib.import_module('real_common')
    rc.CACHE = str(LEGACY_CACHE)
    rc.PATCH_PX = patch
    return rc.run_patch('giza', 'khafre', bank, verbose=False)


def section(out, pool):
    """The published-style section and the depth axis (relabelled metres), with each pixel's shifts pooled first."""
    kz = out['kz']
    order = np.argsort(kz)
    kz = kz[order]
    nr, nc = len(out['grid_rows']), len(out['grid_cols'])
    q = out['q'][:, order, :].astype(np.float64).reshape(nr, nc, len(kz), 2)
    if pool > 1:
        q = uniform_filter(q, size=(pool, pool, 1, 1), mode='nearest')
    zrep = 2 * np.pi / np.median(np.diff(kz))
    z = np.linspace(0.0, SHOW * zrep, NZ)
    mid = nr // 2
    qc = q[mid - SLAB:mid + SLAB + 1, ..., 0] + 1j * q[mid - SLAB:mid + SLAB + 1, ..., 1]
    T = focus_paper(qc.reshape(-1, len(kz)), kz, z).reshape(2 * SLAB + 1, nc, NZ)
    return T.mean(axis=0), z * 648.0 / zrep, zrep


def measure(S, zm):
    """Rungs and room for the strongest column of a section S [range position, depth]."""
    dz = zm[1] - zm[0]
    body = (zm > 0.1 * 648) & (zm < 0.9 * 648)
    c = int(np.argmax(S[:, body].mean(axis=1)))
    c = min(max(c, COLUMN_HALF), S.shape[0] - COLUMN_HALF - 1)
    cols = np.arange(c - COLUMN_HALF, c + COLUMN_HALF + 1)
    L = np.log10(np.maximum(S[cols], 1e-30))
    B = L - gaussian_filter1d(L, TREND_M / dz, axis=1, mode='nearest')
    rung = float(np.mean([np.corrcoef(B[i, body], B[i + 1, body])[0, 1] for i in range(len(cols) - 1)]))
    i_rep = int(np.argmin(np.abs(zm - 648.0)))
    room = float(np.mean(S[cols, i_rep] / S[cols][:, body].mean(axis=1)))
    return {'column_m': float(c), 'rungs': rung, 'room': room, 'makes_rungs': rung >= RUNGS_AT, 'makes_room': room >= ROOM_AT}


def main():
    params = {'crop': 'the first investigation\'s Khafre patch (legacy real_common, giza/khafre, the 2025 pass)',
              'patches_px': PATCHES, 'pools': POOLS, 'banks': BANKS, 'slab_lines': 2 * SLAB + 1,
              'rungs_threshold': RUNGS_AT, 'room_threshold': ROOM_AT, 'trend_m': TREND_M, 'column_positions': 2 * COLUMN_HALF + 1}
    with Run(RID, 'What would give the published pictures their rungs and rooms', params) as run:
        combos = [(p, 'paper') for p in PATCHES] + [(32, 'v17')]
        with ProcessPoolExecutor(max_workers=len(combos)) as pool:
            outs = dict(zip(combos, pool.map(trajectories, combos)))
        rows, sections = [], {}
        for (patch, bank), out in outs.items():
            for pool_n in POOLS:
                S, zm, zrep = section(out, pool_n)
                m = measure(S, zm)
                label = f"{'the paper' if bank == 'paper' else '20 pairs'}, {patch} px, pooled {pool_n} x {pool_n}"
                rows.append({'patch_px': patch, 'bank': bank, 'pool': pool_n, 'label': label, 'repeat_raw_m': float(zrep), **m})
                sections[label] = (S, zm)
                print(f"  {label}: rungs {m['rungs']:.2f}, room {m['room']:.2f}", flush=True)
        # the sections, drawn as the published pictures are: log scale, smoothed, rainbow palette
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        n = len(sections)
        fig, axs = plt.subplots(1, n, figsize=(2.1 * n, 4.6))
        for ax, (label, (S, zm)) in zip(axs, sections.items()):
            A = gaussian_filter(np.log10(np.maximum(S.T, 1e-30)), (2.0, 1.0))
            lo, hi = np.percentile(A, 2), np.percentile(A, 99.5)
            ax.imshow(np.clip((A - lo) / (hi - lo), 0, 1), cmap='jet', aspect='auto', extent=(0, S.shape[0] * 3.46, zm[-1], 0))
            ax.axhline(648, color='w', lw=0.5, ls=':')
            ax.set_title(label.replace(', pooled', '\npooled'), fontsize=7)
            ax.set_xticks([])
            ax.tick_params(labelsize=6)
        plt.tight_layout()
        plt.savefig(figs / 'sections.png', dpi=110)
        base = next(r for r in rows if r['patch_px'] == 32 and r['bank'] == 'paper' and r['pool'] == 1)
        rungs = [r for r in rows if r['makes_rungs']]
        rooms = [r for r in rows if r['makes_room']]
        best_r = max(rows, key=lambda r: r['rungs'])
        best_m = max(rows, key=lambda r: r['room'])
        finding = (
            f"As the paper states it (50 pairs, 32 px patches, no pooling), the strongest column's banding correlates at "
            f"{base['rungs']:.2f} between neighbouring range positions and its repeat layer stands {base['room']:.1f} times its "
            f"mean. " + (f"Rungs (0.5 or more) come with {', '.join(r['label'] for r in rungs)}; "
                         if rungs else f"No setting tried lines the bands up into rungs (the best, {best_r['label']}, reaches "
                                       f"{best_r['rungs']:.2f}); ")
            + (f"a room (twice the mean or more) with {', '.join(r['label'] for r in rooms)}."
               if rooms else f"none makes a room at the repeat (the brightest, {best_m['label']}, {best_m['room']:.1f} times "
                             f"the mean).")
            + (" So the published look can come from choices the paper does not state." if rungs and rooms else
               " So part of the published look can come from choices the paper does not state; the rest needs others, or "
               "rendering beyond these." if rungs or rooms else
               " The published look needs choices the paper does not state beyond these, or rendering."))
        run.save({'rows': rows, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
