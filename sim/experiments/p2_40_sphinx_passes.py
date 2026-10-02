"""P2-40 · Whether the Sphinx volume's blobs stay put from one pass to the next.

    uv run python experiments/p2_40_sphinx_passes.py

The lab's command ran the gated reconstruction across the Sphinx with its author's Sphinx settings (Track 7, supports
1 to 3, 0 to 50 m) along 151 north-south lines 0.6 m apart, from the 2022 pass (lab_sphinx) and, on the same lines, from
the 2025 pass (lab_sphinx_2025), each beside a motionless copy of its image. The 2022 volume shows two dense blobs that
look like chambers. Stated before the 2025 run finished: a void stays where it is from pass to pass; a reading of
the image's texture reshuffles; rock standing above the surface the lines were laid on appears shifted toward the
satellite, by its height times the cotangent of the incidence, so between these passes (incidence 35.0 and 20.8 degrees,
both looking east-north-east) by 1.2 times its height: about 24 m for the statue's 20 m.

The blobs are taken from the 2022 real image as its two largest clusters of positions passing at a support of one
(neighbours touching, diagonals included). Measured here, in each of the four maps (each pass and its copy): the share
of each blob's positions passing, against the share elsewhere; at each blob, the depth of each passing position's best
score, in metres and in turns of that pass's fit (one turn spans 2 pi / (50 dKz), from each run's own steering
wavenumbers, as P2-34 reads them); for the blob that does not stay put, the share passing in 2025 where it would stand
if it were rock standing 0 to 30 m above the lines, carried toward the satellite by 1.2 times that height; and over the
whole square, the positions passing in both passes, against chance (a b / n) and against the same maps shifted by more
than their clumps (wrapping), with the two copies' agreement as the texture's own. The blobs' depths were added after
the four maps had been seen; a void would also stand at one depth.
"""
import hashlib
import json
import re
import subprocess
import urllib.request
import zipfile
from io import BytesIO
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from katabasis.runs import RESULTS, Run, load

RID = 'p2_40_sphinx_passes'
DATA = Path(__file__).resolve().parents[1] / 'data'
RUNS = {'2022': 'lab_sphinx', '2025': 'lab_sphinx_2025'}
WINDOW = 50                      # pairs per window in the frozen profile (P2-34)
HEIGHTS = np.arange(0.0, 30.1, 2.5)
NULL_SHIFTS, MIN_SHIFT = 400, 25 # shifted maps, each moved at least 25 positions (15 m) one way or the other
REPO = 'https://raw.githubusercontent.com/BiondiProtocol/Replication-and-Verification-Biondi-Protocol/c0cb293'
KMZ = 'ICEYE_X13_Great_Sphinx_Ditch_EastWoodAvoidance_GE.kmz'
ORIGIN = (29.976025, 31.130794)  # the Giza site frame's origin (katabasis.compose.build_giza)


def kmz_lines():
    """The statue's outline (the author's enclosure references) and his revised line, in the site frame."""
    from katabasis.compose.geo import to_local
    data = subprocess.run(['curl', '-sSfL', f'{REPO}/{urllib.request.quote(KMZ)}'], check=True, capture_output=True).stdout
    kml = zipfile.ZipFile(BytesIO(data)).read('doc.kml').decode()
    out = {}
    for m in re.finditer(r'<Placemark>(.*?)</Placemark>', kml, re.S):
        name = re.search(r'<name>(.*?)</name>', m.group(1)).group(1)
        c = re.search(r'<coordinates>(.*?)</coordinates>', m.group(1), re.S)
        if c and len(c.group(1).split()) > 2:
            p = np.array([[float(v) for v in x.split(',')] for x in c.group(1).split()])
            out[name] = to_local(p[:, 1], p[:, 0], *ORIGIN)[:, :2]
    return out, hashlib.sha256(data).hexdigest()


def acquisition(product):
    import h5py
    with h5py.File(Path.home() / 'tmp/sar' / product, 'r') as f:
        return {'incidence_deg': float(f['incidence_center'][()]), 'heading_deg': float(f['heading'][()]),
                'look_side': f['look_side'][()].decode()}


def main():
    sums = {k: load(r) for k, r in RUNS.items()}
    products = {'2022': 'giza2/ICEYE_X13_SLC_SLED_868226_20220715T235744.h5',
                '2025': 'giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5'}
    for k, s in sums.items():
        assert Path(products[k]).name == s['manifest']['params']['product'], (k, s['manifest']['params']['product'])
    params = {'runs': RUNS, 'run_commits': {k: s['manifest']['commit'] for k, s in sums.items()}, 'window_pairs': WINDOW,
              'heights_m': HEIGHTS.tolist(), 'null_shifts': NULL_SHIFTS, 'min_shift_positions': MIN_SHIFT, 'kmz': f'{REPO}/{KMZ}'}
    with Run(RID, "Whether the Sphinx volume's blobs stay put from one pass to the next", params) as run:
        V = {k: np.load(RESULTS / r / 'volumes.npz') for k, r in RUNS.items()}
        x, y, step = V['2022']['x'], V['2022']['y'], float(V['2022']['step'])
        assert np.allclose(V['2025']['x'], x) and np.allclose(V['2025']['y'], y)
        X, Y = np.meshgrid(x, y, indexing='ij')
        maps = {}
        for k in RUNS:
            for case, label in (('real', k), ('twin101', f'{k} copy')):
                maps[label] = np.isfinite(V[k][f'{case}_p1']).any(axis=2)
        n = maps['2022'].size
        acq = {k: acquisition(products[k]) for k in RUNS}
        turn = {}
        for k, r in RUNS.items():
            kz = np.sort(np.load(DATA / r / 'real' / 'biondi_v1_8_audit.npz', allow_pickle=True)['kz_rad_per_m'], axis=1)
            turn[k] = float(2 * np.pi / (WINDOW * np.median(np.diff(kz, axis=1))))

        outline, kmz_sha = kmz_lines()
        ref = np.vstack([outline[k] for k in outline if k.startswith(('North reference', 'South reference', 'West reference'))])
        west_x, east_x = float(ref[:, 0].min()), float(ref[:, 0].max())
        south = outline[next(k for k in outline if k.startswith('South reference'))]

        def where(c):
            """Where a point stands along the statue's outline (the author's references), west to east."""
            f = (c[0] - west_x) / (east_x - west_x)
            if f < 0.25:
                return 'west end'
            if f > 0.75:
                return 'east end'
            return 'south side' if abs(c[1] - np.interp(c[0], south[:, 0], south[:, 1])) < 5 else 'middle'

        # the blobs: the 2022 real image's two largest clusters
        lab, nl = ndi.label(maps['2022'], structure=np.ones((3, 3)))
        size = ndi.sum(maps['2022'], lab, range(1, nl + 1))
        order = np.argsort(size)[::-1][:2]
        blobs = []
        for i in order:
            m = lab == i + 1
            blobs.append({'mask': m, 'positions': int(m.sum()), 'centre_m': [float(X[m].mean()), float(Y[m].mean())],
                          'x_m': [float(X[m].min()), float(X[m].max())], 'y_m': [float(Y[m].min()), float(Y[m].max())],
                          'where': where((float(X[m].mean()), float(Y[m].mean())))})
        rest = ~(blobs[0]['mask'] | blobs[1]['mask'])
        for b in blobs:
            b['passing'] = {k: float(m[b['mask']].mean()) for k, m in maps.items()}
            b['depths'] = {}
            for k in RUNS:
                F = V[k]['real_p1'].astype(np.float32)
                ok = np.isfinite(F)
                best = (V[k]['z_surface'].astype(np.float32) - V[k]['z'][np.argmax(np.where(ok, F, -np.inf), axis=2)])
                d = best[b['mask'] & ok.any(axis=2)]
                q = np.percentile(d, [25, 50, 75]) if d.size else [np.nan] * 3
                b['depths'][k] = {'quartiles_m': [float(v) for v in q], 'median_turns': float(q[1] / turn[k])}
        elsewhere = {k: float(m[rest].mean()) for k, m in maps.items()}

        # rock above the lines: where each blob would stand in 2025, carried toward the satellite by 1.2 times its height
        cot = {k: 1 / np.tan(np.deg2rad(a['incidence_deg'])) for k, a in acq.items()}
        toward = np.deg2rad(np.mean([a['heading_deg'] for a in acq.values()]) + (90 if acq['2022']['look_side'] == 'right' else -90) + 180)
        u = np.array([np.sin(toward), np.cos(toward)])                     # east, north: from the ground toward the satellite
        for b in blobs:
            ii, jj = np.nonzero(b['mask'])
            b['if_rock'] = []
            for h in HEIGHTS:
                d = h * (cot['2025'] - cot['2022']) * u / step
                i2, j2 = np.round(ii + d[0]).astype(int), np.round(jj + d[1]).astype(int)
                inside = (i2 >= 0) & (i2 < len(x)) & (j2 >= 0) & (j2 < len(y))
                b['if_rock'].append({'height_m': float(h), 'shift_m': float(h * (cot['2025'] - cot['2022'])),
                                     'inside_share': float(inside.mean()),
                                     'passing_2025': float(maps['2025'][i2[inside], j2[inside]].mean()) if inside.any() else None})

        # the whole square: positions passing in both passes, against chance and against shifted maps
        rng = np.random.default_rng(0)

        def both(p, q):
            a, b_ = maps[p], maps[q]
            null = []
            while len(null) < NULL_SHIFTS:
                dx, dy = rng.integers(-len(x) // 2, len(x) // 2 + 1), rng.integers(-len(y) // 2, len(y) // 2 + 1)
                if abs(dx) < MIN_SHIFT and abs(dy) < MIN_SHIFT:
                    continue
                null.append(int((a & np.roll(np.roll(b_, dx, 0), dy, 1)).sum()))
            k = int((a & b_).sum())
            return {'both': k, 'chance': float(a.sum() * b_.sum() / n), 'ratio': float(k / (a.sum() * b_.sum() / n)),
                    'shifted_mean': float(np.mean(null)), 'shifted_p95': float(np.percentile(null, 95)),
                    'share_shifted_at_or_above': float(np.mean(np.array(null) >= k))}
        agreement = {f'{p} & {q}': both(p, q) for p, q in (('2022', '2025'), ('2022 copy', '2025 copy'), ('2022', '2025 copy'), ('2022 copy', '2025'))}

        # the figure: the four maps, the blobs outlined, the statue's outline and the author's line from his KMZ
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        fig, axs = plt.subplots(2, 2, figsize=(9.4, 6.3), sharex=True, sharey=True)
        ext = (x[0] - step / 2, x[-1] + step / 2, y[0] - step / 2, y[-1] + step / 2)
        for ax, label in zip(axs.ravel(), ['2022', '2022 copy', '2025', '2025 copy']):
            m = maps[label]
            ax.imshow(m.T, origin='lower', extent=ext, cmap='Greys', vmin=0, vmax=1.5, interpolation='nearest')
            for name, p in outline.items():
                if name.startswith(('North reference', 'South reference', 'West reference', 'Revised East: NE')):
                    ax.plot(p[:, 0], p[:, 1], color='#3a7bd0' if name.startswith('Revised') else '#9a9a9a', lw=0.9)
            for b, c in zip(blobs, ('#d0603a', '#e8a33a')):
                ax.contour(X, Y, b['mask'].astype(float), levels=[0.5], colors=[c], linewidths=1.0)
            title = {'2022': 'the 2022 pass', '2022 copy': 'its motionless copy', '2025': 'the 2025 pass, same lines',
                     '2025 copy': 'its motionless copy'}[label]
            ax.set_title(f'{title}: {m.sum():,} of {n:,} pass', fontsize=8.5)
            ax.set_aspect('equal')
            ax.tick_params(labelsize=7)
        for ax in axs[1]:
            ax.set_xlabel('metres east (site frame)', fontsize=7.5)
        for ax in axs[:, 0]:
            ax.set_ylabel('metres north', fontsize=7.5)
        plt.tight_layout()
        plt.savefig(figs / 'passes.png', dpi=130)

        big, small = blobs
        moved = max(big['if_rock'], key=lambda r: r['passing_2025'] or 0)
        name = lambda b: b['where']
        a22 = agreement['2022 & 2025']['ratio']
        acp = agreement['2022 copy & 2025 copy']['ratio']
        big_gone = big['passing']['2025'] < 2 * elsewhere['2025'] and moved['passing_2025'] < 2 * elsewhere['2025']
        copies_gather = all(small['passing'][c] >= 2 * elsewhere[c] for c in ('2022 copy', '2025 copy'))
        depth_moves = abs(small['depths']['2022']['quartiles_m'][1] - small['depths']['2025']['quartiles_m'][1]) > 1.0
        finding = (
            f"The 2022 volume's largest blob ({big['positions']} positions, at the statue's {name(big)}) does not come back on the 2025 "
            f"pass: {big['passing']['2025'] * 100:.0f}% of its positions pass there, where {elsewhere['2025'] * 100:.0f}% pass elsewhere; "
            f"carried toward the satellite as rock standing up to {HEIGHTS[-1]:.0f} m above the lines would be, it reaches at most "
            f"{moved['passing_2025'] * 100:.0f}% (at {moved['height_m']:.1f} m). The second ({small['positions']} positions, at the "
            f"statue's {name(small)}) stays put: {small['passing']['2025'] * 100:.0f}% of its positions pass in 2025, against "
            f"{elsewhere['2025'] * 100:.0f}% elsewhere, and the motionless copies, with nothing moving, gather there too "
            f"({small['passing']['2022 copy'] * 100:.0f}% and {small['passing']['2025 copy'] * 100:.0f}%, against "
            f"{elsewhere['2022 copy'] * 100:.0f}% and {elsewhere['2025 copy'] * 100:.0f}%). Its depth does not hold: its best scores "
            f"stand at a median {small['depths']['2022']['quartiles_m'][1]:.1f} m in 2022 and {small['depths']['2025']['quartiles_m'][1]:.1f} m "
            f"in 2025 ({small['depths']['2022']['median_turns']:.2f} and {small['depths']['2025']['median_turns']:.2f} turns of fits that "
            f"span {turn['2022']:.2f} and {turn['2025']:.2f} m a turn), where a void would stand at one depth. Across the square, {agreement['2022 & 2025']['both']:,} positions "
            f"pass in both passes against {agreement['2022 & 2025']['chance']:,.0f} by chance "
            f"({a22:.2f} times; shifted maps {agreement['2022 & 2025']['shifted_mean']:,.0f}), and the two motionless copies agree "
            f"{acp:.2f} times chance: " + ("the passes agree no more than the images' texture does on its own. " if a22 <= acp else
                                           "the passes agree more than the copies do. ")
            + ("Neither blob behaves as a void would." if big_gone and copies_gather and depth_moves else
               "Whether either blob behaves as a void would is not settled by these numbers."))
        run.save({'acquisitions': acq, 'turn_m': turn, 'elsewhere': elsewhere,
                  'blobs': [{k: v for k, v in b.items() if k != 'mask'} for b in blobs], 'agreement': agreement,
                  'kmz_sha256': kmz_sha, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
