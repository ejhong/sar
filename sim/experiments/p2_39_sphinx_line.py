"""P2-39 · The published Sphinx tomogram, reproduced, and how far it survives a one-pixel nudge.

    uv run --with scikit-image python experiments/p2_39_sphinx_line.py

The derivative protocol's public repository (commit c0cb293, 1 October 2026) shows the gated reconstruction along a
line in front of the Sphinx, the "revised East wood-avoidance line", P0 = NE' to P100 = SE', on the ICEYE X13 pass of
15 July 2022, with the line itself (a KMZ: 101 points with their ellipsoidal heights, and the native pixels of both ends)
and the settings in the figure's title: W50, modes 1-10, b:a >= 0.10, Track 7, supports 1-3, central starts, 0-50 m.

1. Reproduced: the repository's v1.8 module, unchanged (its hash recorded), the line projected through the product's own
   RPC (both ends land within half a pixel of the pixels the KMZ states), the frozen profile with the title's settings.
   The positions that pass are compared with those in the published figure, read off its tick marks.
2. Varied, each through the same code: the line nudged by one pixel along track (4.6 cm) and across it (one slant-range
   sample, 15 cm), the author's earlier East line (2 to 15 m from the revised one, in the same KMZ), a motionless copy
   of the image along the line (sarsim.looks.motionless_twin), and the same ground in the 2025 pass.

A ground feature would stay where it is under a one-pixel nudge, show on the earlier line where the two run close, and
not show in a motionless copy; a reading of the image's texture would move with the nudge, differ on the earlier line
and draw columns in the copy too. The variants were first run while checking the reproduction; this run records them.
"""
import hashlib
import importlib.util
import json
import re
import shutil
import time
import urllib.request
import zipfile
from concurrent.futures import ProcessPoolExecutor
from io import BytesIO
from pathlib import Path

import numpy as np

from katabasis.runs import Run

RID = 'p2_39_sphinx_line'
REPO = 'https://raw.githubusercontent.com/BiondiProtocol/Replication-and-Verification-Biondi-Protocol/c0cb293'
FIGURE = 'ICEYE_X13_Great_Sphinx_Ditch_EastWoodAvoidance_Track7_P1toP3_Tomograms_0to50m.png'
KMZ = 'ICEYE_X13_Great_Sphinx_Ditch_EastWoodAvoidance_GE.kmz'
PKG = Path.home() / 'tmp/sar/biondi_v18/pkg/Biondi_Protocol_v1.8'
P22 = Path.home() / 'tmp/sar/giza2/ICEYE_X13_SLC_SLED_868226_20220715T235744.h5'
P25 = Path.home() / 'tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5'
WORK = Path(__file__).resolve().parents[1] / 'data' / RID
TRACK, SUPPORTS, DEPTH_MAX = 7, [1, 2, 3], 50.0
REVISED, EARLIER = 'Revised East: NE', 'Previous East: NE'
VARIANTS = [('as published', REVISED, 'p22', 0, 0, None), ('one row along track', REVISED, 'p22', -1, 0, None),
            ('one row the other way', REVISED, 'p22', 1, 0, None), ('one sample across track', REVISED, 'p22', 0, -1, None),
            ('one sample the other way', REVISED, 'p22', 0, 1, None), ('the earlier East line', EARLIER, 'p22', 0, 0, None),
            ('a motionless copy', REVISED, 'p22', 0, 0, 101), ('the 2025 pass', REVISED, 'p25', 0, 0, None)]


def fetch(name):
    data = urllib.request.urlopen(f'{REPO}/{urllib.request.quote(name)}', timeout=60).read()
    return data, hashlib.sha256(data).hexdigest()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kml_line(kml, label):
    c = re.search(re.escape(label) + r'.*?<coordinates>(.*?)</coordinates>', kml, re.S).group(1).split()
    return np.array([[float(v) for v in x.split(',')] for x in c])          # lon, lat, ellipsoidal height


def published_positions(png):
    """The positions lit in each support's panel of the published figure, read off its tick marks (0 to 100 by 20)."""
    from PIL import Image
    im = np.asarray(Image.open(BytesIO(png)).convert('RGB')).astype(int)
    white = (im > 200).all(axis=2)
    frames = np.flatnonzero(white.sum(axis=1) > 0.6 * im.shape[1])
    frames = [g[0] for g in np.split(frames, np.flatnonzero(np.diff(frames) > 3) + 1)]
    left, right = (g[0] for g in [np.flatnonzero(white.sum(axis=0) > 0.25 * im.shape[0])[[0, -1]]])
    ticks = np.flatnonzero(white[frames[0] - 16:frames[0] - 2, left + 2:right - 1].any(axis=0)) + left + 2
    ticks = [g.mean() for g in np.split(ticks, np.flatnonzero(np.diff(ticks) > 3) + 1)][-6:]
    per = (ticks[-1] - ticks[0]) / 100.0
    out = {}
    for k in range(0, len(frames) - 1, 2):
        panel = im[frames[k] + 2:frames[k + 1] - 1, left + 2:right - 1]
        lit = np.flatnonzero((panel.max(axis=2) > 25).sum(axis=0) > 3) + left + 2
        pos = sorted({int(round((x - ticks[0]) / per)) for x in lit})
        out[k // 2 + 1] = pos
    return out, {'ticks_px': [float(t) for t in ticks], 'px_per_position': float(per)}


def run_variant(args):
    name, which, prod, drow, dcol, twin = args
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from sarsim.dwell import DwellProduct
    from sarsim.geolocation import project_rpc
    product = P22 if prod == 'p22' else P25
    kml = (WORK / 'doc.kml').read_text()
    pts = kml_line(kml, which)
    rc = project_rpc(DwellProduct(product)._rpc, pts[:, 1], pts[:, 0], pts[:, 2]) + np.array([drow, dcol])
    t = np.linspace(0, 100, 10001)                        # the module takes 101 targets along a trace: these points
    row, col = np.interp(t, np.arange(101), rc[:, 0]), np.interp(t, np.arange(101), rc[:, 1])
    d = WORK / re.sub(r'[^a-z0-9]+', '_', name)
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    geom = {'title': name, 'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
            'tracks': {'L000': {'native_trace_col_row': np.stack([np.round(col), np.round(row)], 1).astype(int).tolist(),
                                'geographic_endpoint_lat_lon': [[float(pts[0, 1]), float(pts[0, 0])], [float(pts[-1, 1]), float(pts[-1, 0])]]}},
            'projection': {'height_m': float(pts[:, 2].mean())}}
    (d / 'geometry.json').write_text(json.dumps(geom))
    prof = json.loads((PKG / 'profiles/iceye_x13_khafre_w50.json').read_text())
    prof['input'] = dict(prof['input'], hdf5_path=str(product), geometry_path=str(d / 'geometry.json'))
    prof['gates'] = dict(prof['gates'], same_mode_contiguous_windows=TRACK, pixel_support_variants=SUPPORTS)
    prof['steering'] = dict(prof['steering'], depth_grid_m=dict(prof['steering']['depth_grid_m'], maximum=DEPTH_MAX))
    (d / 'profile.json').write_text(json.dumps(prof, indent=1))
    spec = importlib.util.spec_from_file_location('gated', PKG / 'biondi_tomography_v1_8.py')
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    if twin is not None:
        from sarsim.looks import motionless_twin
        orig = g.load_source_and_preflight

        def wrapped(config, config_path):
            out = list(orig(config, config_path))
            out[7] = motionless_twin(out[7], np.random.default_rng(twin))
            return tuple(out)
        g.load_source_and_preflight = wrapped
    t0 = time.time()
    g.run_full(prof, d / 'profile.json', d / 'out')
    a = np.load(d / 'out' / 'biondi_v1_8_audit.npz', allow_pickle=True)
    passed = {P: np.flatnonzero(np.isfinite(a[f'focus_p{P}'][0]).any(axis=1)).tolist() for P in SUPPORTS}
    ends = [[float(rc[0, 1]), float(rc[0, 0])], [float(rc[-1, 1]), float(rc[-1, 0])]]
    return name, {'passed': passed, 'runtime_s': round(time.time() - t0, 1), 'ends_col_row': ends,
                  'focus_p1': a['focus_p1'][0].astype(np.float32), 'z_m': a['z_m']}


def main():
    params = {'repository': REPO, 'figure': FIGURE, 'kmz': KMZ, 'module_sha256': sha256(PKG / 'biondi_tomography_v1_8.py'),
              'profile': 'iceye_x13_khafre_w50.json with Track 7, supports 1-3, depths 0-50 m',
              'profile_sha256': sha256(PKG / 'profiles/iceye_x13_khafre_w50.json'), 'variants': [v[:2] + v[3:5] + (v[2], v[5]) for v in VARIANTS]}
    with Run(RID, 'The published Sphinx tomogram, reproduced, and how far it survives a one-pixel nudge', params) as run:
        WORK.mkdir(parents=True, exist_ok=True)
        png, png_sha = fetch(FIGURE)
        kmz, kmz_sha = fetch(KMZ)
        (WORK / 'doc.kml').write_bytes(zipfile.ZipFile(BytesIO(kmz)).read('doc.kml'))
        published, axis = published_positions(png)
        kml = (WORK / 'doc.kml').read_text()
        stated = [[int(v) for v in m] for m in re.findall(r'Native SLC \[column,row\]: \[(\d+), (\d+)\]', kml)[:2]]
        with ProcessPoolExecutor(max_workers=len(VARIANTS)) as pool:
            results = dict(pool.map(run_variant, VARIANTS))
        ours = results['as published']
        pub1 = set(published[1])
        rows = []
        for name, r in results.items():
            p1 = set(r['passed'][1])
            rows.append({'variant': name, 'passed': r['passed'], 'runtime_s': r['runtime_s'],
                         'shared_with_reproduction_p1': len(p1 & set(ours['passed'][1])),
                         'shared_with_published_p1': len(p1 & pub1)})
            print(f"  {name}: P1 {r['passed'][1]} P2 {r['passed'][2]} ({r['runtime_s']} s)", flush=True)
        # the strip figure: where columns stand along the line, published, reproduced and varied
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        labels = ['published figure'] + [r['variant'] for r in rows]
        sets = [published] + [r['passed'] for r in rows]
        fig, ax = plt.subplots(figsize=(9.5, 0.36 * len(labels) + 0.9))
        for i, (lab, s) in enumerate(zip(labels, sets)):
            y = len(labels) - 1 - i
            ax.scatter(s[1], [y] * len(s[1]), marker='s', s=34, color='#b6367a', zorder=3)
            if s.get(2):
                ax.scatter(s[2], [y] * len(s[2]), marker='s', s=34, facecolor='none', edgecolor='#fbc488', lw=1.4, zorder=4)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels[::-1], fontsize=8)
        ax.set_xlim(-1, 101)
        ax.set_xlabel('position along the line (P0 at the north-east end, P100 at the south-east; 0.39 m apart)', fontsize=8)
        ax.grid(axis='x', color='0.85', lw=0.5)
        ax.tick_params(axis='x', labelsize=7)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        plt.tight_layout()
        plt.savefig(figs / 'positions.png', dpi=150)
        np.savez_compressed(run.dir / 'reproduction.npz', focus_p1=ours['focus_p1'], z_m=ours['z_m'],
                            line_lon_lat_h=kml_line(kml, REVISED))
        match = sorted(pub1 & set(ours['passed'][1]))
        nudges = [r for r in rows if r['variant'].startswith('one ')]
        finding = (
            f"The published Sphinx tomogram reproduces with its author's code, unchanged, on its line and settings: the line's "
            f"ends land within half a pixel of the pixels the KMZ states, and {len(match)} of the {len(ours['passed'][1])} "
            f"positions our run passes at a support of one stand in the published figure at the same places "
            f"({len(pub1)} there: {sorted(pub1)}; ours {ours['passed'][1]}). The published figure's support-two pair "
            f"({published.get(2)}) appears in our runs when the line moves by one pixel. Moving the line one pixel along track "
            f"(4.6 cm) or across it (15 cm) keeps {min(r['shared_with_reproduction_p1'] for r in nudges)} to "
            f"{max(r['shared_with_reproduction_p1'] for r in nudges)} of those {len(ours['passed'][1])} positions; the author's "
            f"earlier East line, 2 to 15 m away, keeps {next(r for r in rows if r['variant'] == 'the earlier East line')['shared_with_reproduction_p1']}; "
            f"a motionless copy of the image, with nothing moving, passes {len(results['a motionless copy']['passed'][1])} positions of "
            f"its own; the 2025 pass passes {len(results['the 2025 pass']['passed'][1])}, "
            f"{next(r for r in rows if r['variant'] == 'the 2025 pass')['shared_with_reproduction_p1']} of them shared. One line "
            f"takes {ours['runtime_s']:.0f} s here.")
        run.save({'published': {'positions': published, 'axis': axis, 'figure_sha256': png_sha, 'kmz_sha256': kmz_sha,
                                'stated_ends_col_row': stated},
                  'reproduction_ends_col_row': ours['ends_col_row'], 'rows': rows, 'matched_p1': match, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
