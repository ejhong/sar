"""P2-37 · The benchmark read: the published method and a detector told the lorry, on images of the benchmark's worlds.

    uv run python experiments/p2_37_benchmark_read.py [worker k n]

P2-36 bounds what one image of the benchmark's worlds could hold (BENCHMARK.md); this reads such images. The same four
worlds (no cavity; a room 6 m on a side under a 5 m roof; an L-shaped tunnel at that depth; the room 6 m east), the same
lorry 15 m west at the FTA's truck-over-a-bump level held at P2-36's frozen frequency for the whole pass, each world's
whole line-of-sight motion from P2-36's solver fields (tapered to zero between 60 and 70 m, where the solver's record
ends), Giza's regional microseisms on top, synthesised on the real dwell (P2-07's scene and chain) in complex128.

Grounds. Independent: each its own speckle, its own microseism realisation and its own lorry timing, shared by every
world on that ground (paired comparisons within a ground; the grounds are the evaluation units).
Levels. The real one; and the cavity's difference (each world's motion less the no-cavity world's) amplified 200 and
1,000 times, labelled diagnostics: the first where P2-36's single-image bound within the recorded disc stops excluding
95% at 5%, the second a positive control where the told detector must succeed if the experiment can reveal a signal.
Readers.
1. The published method (P2-07's pipeline): its focused power summed over depth, read blind as a plan map; and its
   change against the same ground's no-cavity image (paired), a diagnostic it never has.
2. A detector told the lorry: the score statistic for each cavity world's difference pattern and timing (told the
   excitation and where to look, not the speckle), sarsim.information's score as P2-36 implements it.
Scored, stated before the run: presence (each image's largest normalised response, cavity images against the controls,
AUC); shape and location together, naming which of the three cavity worlds (chance one in three); location alone, the
room or the room 6 m east (chance one in two); shape alone, the room or the L tunnel (chance one in two); for the
published maps, the centroid of the top 5% against the true footprint's centroid, with its chance from the controls.
p-values are null tests with the true labels permuted within each ground. Amplified levels are diagnostics, not
physical predictions.
"""
import importlib.util
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.special import ndtr

from katabasis.ambient.field import microseisms, rayleigh_hv
from katabasis.runs import Run, load, memo
from sarsim import information as inf
from sarsim import synthesize

RID = 'p2_37_benchmark_read'
SITES = Path(__file__).resolve().parents[2] / 'sites'
GROUNDS = tuple(range(3701, 3709))
LEVELS = {'real': 1.0, 'x200': 200.0, 'x1000': 1000.0}
CAVITIES = ('WA', 'WB', 'WA6')
TAPER = (60.0, 70.0)
TOP_SHARE = 0.05
FIELD_SEED = 2000
READINGS = 'read-v1'      # the version of the scenes and readers; bump it if either changes


def module(name, path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


HERE = Path(__file__).parent
m07 = module('p207', HERE / 'p2_07_whole_chain.py')
m31 = module('p231', HERE / 'p2_31_shape_test.py')
m36 = module('p236', HERE / 'p2_36_benchmark.py')


class Lorry:
    """Each scatterer's slant-range increase: the world's line-of-sight field D_0 + amp (D_W - D_0) as a harmonic at f,
    each row at its own zero-Doppler time, tapered to zero where the solver's record ends; the microseisms added."""

    def __init__(self, g, xs, ys, Z, f_hz, psi, micro):
        self.g, self.f, self.psi, self.micro = g, f_hz, psi, micro
        self.re = RegularGridInterpolator((xs, ys), Z.real, bounds_error=False, fill_value=0.0)
        self.im = RegularGridInterpolator((xs, ys), Z.imag, bounds_error=False, fill_value=0.0)
        self.a, self.r = np.asarray(g.along_track_en), np.asarray(g.ground_range_en)

    def field(self, E, N):
        p = np.stack([E, N], -1)
        tp = np.clip((TAPER[1] - np.hypot(E, N)) / (TAPER[1] - TAPER[0]), 0, 1)
        return (self.re(p) + 1j * self.im(p)) * tp

    def __call__(self, x, y, z, t):
        E, N = x * self.a[0] + y * self.r[0], x * self.a[1] + y * self.r[1]
        Zs = self.field(E, N)
        ph = np.exp(1j * (2 * np.pi * self.f * (t[:, None] + x[None, :] / self.g.V) + self.psi))
        d = np.real(Zs[None, :] * ph)
        return d + (self.micro(x, y, z, t) if self.micro is not None else 0.0)


def setup():
    g = m36.geometry()
    F = m36.strong_fields(g)
    p36 = load('p2_36_benchmark')
    f_star = float(p36['f_star_hz'])
    fi = int(np.argmin(np.abs(F['f'] - f_star)))
    amb = json.loads((SITES / 'ambient.json').read_text())
    vp, vs = 3300.0, 1830.0
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    rock = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    vp, vs = val('vp'), val('vs')
    level = amb['regional']['microseism_vertical_0.1_0.3_hz']['value']
    speed = amb['regional']['microseism_phase_speed']['value']
    fq = np.linspace(0.1, 0.3, 21)

    def micro_for(seed):
        fld = microseisms(fq, np.full(len(fq), level ** 2 / 0.2), np.random.default_rng(FIELD_SEED + seed), per_bin=4,
                          speed=speed, hv=rayleigh_hv(vp, vs))
        zero = [lambda p: np.zeros(len(p))] * 3
        return m07.Shaking(g, fld, zero, 0.0)
    SH = m07.SHAPE
    xa = (np.arange(SH[0]) - SH[0] // 2) * g.dx
    xg = (np.arange(SH[1]) - SH[1] // 2) * g.dr / np.sin(g.theta)
    XA, XG = np.meshgrid(xa, xg, indexing='ij')
    PE, PN = m31.to_site(g, XA, XG)
    rows = np.arange(48, SH[0] - 48, m07.STRIDE[0])
    cols = np.arange(32, SH[1] - 32, m07.STRIDE[1])
    RR, CC = np.meshgrid(rows, cols, indexing='ij')
    return {'g': g, 'F': F, 'fi': fi, 'f': f_star, 'micro_for': micro_for, 'XA': XA, 'PE': PE, 'PN': PN, 'RR': RR,
            'CC': CC, 'GE': PE[RR, CC], 'GN': PN[RR, CC]}


def templates(S, n, psi):
    """The told detector's two quadrature templates for world n's difference, on the image's pixels (the same taper)."""
    g, F, fi = S['g'], S['F'], S['fi']
    D0, Dn = F['D']['W0'][fi], F['D'][n][fi]
    lor = Lorry(g, F['xs'], F['ys'], Dn - D0, S['f'], 0.0, None)
    dK = lor.field(S['PE'], S['PN'])
    ph = np.exp(1j * (2 * np.pi * S['f'] * S['XA'] / g.V + psi))
    k0 = 4 * np.pi / g.lam
    return k0 * np.real(dK * ph), -k0 * np.imag(dK * ph), dK


def reading(S, seed, key):
    g, F, fi = S['g'], S['F'], S['fi']
    rng = np.random.default_rng(seed)
    scat = m07.scene(g, rng)
    psi = float(np.random.default_rng(seed + 50000).uniform(0, 2 * np.pi))
    SH = m07.SHAPE
    if key == 'motionless':
        img = synthesize(scat, g, SH, dtype=np.complex128)
    else:
        world, lv = key.split('|')
        D0 = F['D']['W0'][fi]
        Z = D0 if world == 'W0' else D0 + LEVELS[lv] * (F['D'][world][fi] - D0)
        img = synthesize(scat, g, SH, motion=Lorry(g, F['xs'], F['ys'], Z, S['f'], psi, S['micro_for'](seed)),
                         dtype=np.complex128)
    T, z = m07.run_method(img, g, S['RR'].ravel(), S['CC'].ravel())
    taus = [lambda t: np.cos(2 * np.pi * S['f'] * t), lambda t: np.sin(2 * np.pi * S['f'] * t)]
    q = m31.score_fields(img, g, taus)
    stats = {}
    for n in CAVITIES:
        K1, K2, _ = templates(S, n, psi)
        stats[n] = float(np.sum(K1 * q[0] + K2 * q[1]))
    return {'seed': seed, 'key': key, 'stats': stats, 'blind': T.sum(axis=1).reshape(S['RR'].shape).astype(np.float32),
            'T': T.astype(np.float32), 'z': z}


def keys():
    return ['W0|real', 'motionless'] + [f'{n}|{lv}' for lv in LEVELS for n in CAVITIES]


def worker(k, n):
    S = setup()
    for i, seed in enumerate(GROUNDS):
        if i % n != k:
            continue
        for key in keys():
            t0 = time.time()
            memo(RID, f'reading-{seed}-{key}', lambda: reading(S, seed, key), __file__, version=READINGS)
            print(f"  ground {seed} {key}: {time.time() - t0:.0f} s", flush=True)


def footprints(S):
    boxes = m36.FROZEN['worlds']
    return {n: m31.plan_mask([(tuple(c), tuple(s)) for c, s in boxes[n]], S['GE'], S['GN']) for n in CAVITIES}


def perm_p_named(rows, reader, names):
    """P(at least the observed number named correctly) with the true labels permuted within each ground (every
    arrangement when few, else 20,000 random ones)."""
    by = {}
    for x in rows:
        if x[reader].get('named') is not None:
            by.setdefault(x['seed'], []).append((x['world'], x[reader]['named']))
    if not by:
        return None
    obs = sum(t == m for v in by.values() for t, m in v)
    groups = list(by.values())
    sizes = [len(v) for v in groups]
    n_arr = int(np.prod([math.factorial(s) for s in sizes]))
    rng = np.random.default_rng(37)
    if n_arr <= 50000:
        per = [[sum(t == m for t, m in zip(p, [m for _, m in v])) for p in itertools.permutations([t for t, _ in v])]
               for v in groups]
        tot = [sum(c) for c in itertools.product(*per)]
        return float(np.mean([t >= obs for t in tot]))
    hits = 0
    for _ in range(20000):
        s = 0
        for v in groups:
            truths = [t for t, _ in v]
            p = rng.permutation(len(truths))
            s += sum(truths[p[j]] == v[j][1] for j in range(len(v)))
        hits += s >= obs
    return float(hits / 20000)


def score(S, recs):
    feet = footprints(S)
    GE, GN = S['GE'], S['GN']
    ctrl = [r for r in recs if r['key'] in ('W0|real', 'motionless')]
    sd = {n: float(np.std([r['stats'][n] for r in ctrl], ddof=1)) for n in CAVITIES}
    pub_norm = lambda M: (M - np.median(M)) / (np.std(M) + 1e-30)
    t_none = {r['seed']: r['T'] for r in recs if r['key'] == 'W0|real'}

    def maps(r):
        out = {'published_blind': r['blind']}
        if r['seed'] in t_none and r['key'] != 'W0|real':
            out['published_paired'] = np.abs(r['T'] - t_none[r['seed']]).sum(axis=1).reshape(S['RR'].shape)
        return out

    def read_pub(M, choices):
        Mn = pub_norm(gaussian_filter(M, 1.0))
        ce, cn, _ = m31.top_centroid(Mn, GE, GN, TOP_SHARE)
        cors = {n: float(np.nan_to_num(np.corrcoef(Mn.ravel(), gaussian_filter(feet[n].astype(float), 1.0).ravel())[0, 1],
                                       nan=-1.0)) for n in choices}
        return Mn, (ce, cn), cors

    rows = []
    for r in recs:
        world, lv = (r['key'].split('|') if '|' in r['key'] else ('motionless', 'control'))
        if world == 'W0':
            lv = 'control'
        z = {n: r['stats'][n] / sd[n] for n in CAVITIES}
        row = {'seed': r['seed'], 'world': world, 'level': lv,
               'told': {'presence': max(z.values()), 'z': z, 'named': max(z, key=z.get) if world in CAVITIES else None}}
        if world in ('WA', 'WA6'):
            row['told']['location_named'] = 'WA' if z['WA'] >= z['WA6'] else 'WA6'
        if world in ('WA', 'WB'):
            row['told']['shape_named'] = 'WA' if z['WA'] >= z['WB'] else 'WB'
        for reader, M in maps(r).items():
            Mn, (ce, cn), cors = read_pub(M, CAVITIES)
            x = {'presence': float(Mn.max()), 'correlation': cors, 'named': max(cors, key=cors.get) if world in CAVITIES else None}
            if world in CAVITIES:
                x['centroid_error_m'] = float(np.hypot(ce - GE[feet[world]].mean(), cn - GN[feet[world]].mean()))
            if world in ('WA', 'WA6'):
                x['location_named'] = 'WA' if cors['WA'] >= cors['WA6'] else 'WA6'
            if world in ('WA', 'WB'):
                x['shape_named'] = 'WA' if cors['WA'] >= cors['WB'] else 'WB'
            row[reader] = x
        rows.append(row)

    def binary_p(rs, reader, field, pair):
        """Exact: labels swapped within each ground's pair."""
        by = {}
        for x in rs:
            if x['world'] in pair and field in x.get(reader, {}):
                by.setdefault(x['seed'], []).append(x[reader][field] == x['world'])
        obs = sum(sum(v) for v in by.values())
        tot = [sum(c) for c in itertools.product(*[[sum(v), len(v) - sum(v)] for v in by.values()])]
        return float(np.mean([t >= obs for t in tot])) if tot else None, obs, sum(len(v) for v in by.values())

    summary = {}
    ctrl_rows = [x for x in rows if x['level'] == 'control']
    chance_cent = []
    for r in ctrl:
        _, (ce, cn), _ = read_pub(r['blind'], CAVITIES)
        chance_cent += [float(np.hypot(ce - GE[feet[n]].mean(), cn - GN[feet[n]].mean())) for n in CAVITIES]
    for lv in LEVELS:
        cav = [x for x in rows if x['level'] == lv]
        s = {}
        for reader in ('told', 'published_blind', 'published_paired'):
            pc = [x[reader]['presence'] for x in cav if reader in x]
            pn = [x[reader]['presence'] for x in ctrl_rows if reader in x]
            auc = float(np.mean([[1.0 if a > b else 0.5 if a == b else 0.0 for b in pn] for a in pc])) if pc and pn else None
            named = [x[reader]['named'] == x['world'] for x in cav if reader in x and x[reader].get('named')]
            loc_p, loc_k, loc_n = binary_p(cav, reader, 'location_named', ('WA', 'WA6'))
            shp_p, shp_k, shp_n = binary_p(cav, reader, 'shape_named', ('WA', 'WB'))
            s[reader] = {'presence_auc_vs_controls': auc, 'named_of_three': float(np.mean(named)) if named else None,
                         'named_of_three_p': perm_p_named(cav, reader, CAVITIES) if named else None,
                         'location_correct': f'{loc_k} of {loc_n}', 'location_p': loc_p,
                         'shape_correct': f'{shp_k} of {shp_n}', 'shape_p': shp_p, 'images': len(pc)}
            if reader != 'told':
                ce = [x[reader]['centroid_error_m'] for x in cav if reader in x and 'centroid_error_m' in x[reader]]
                s[reader]['median_centroid_error_m'] = float(np.median(ce)) if ce else None
            else:
                s[reader]['median_presence_z'] = float(np.median(pc)) if pc else None
        summary[lv] = s
    summary['control_centroid_error_m'] = float(np.median(chance_cent))
    return rows, summary, sd


def predicted(S):
    """The told detector's weak-signal deflection per level on this image (exact Fisher information on the image's pixels,
    fisher_grid), for presence and for the location and shape pairs."""
    g = S['g']
    out = {}
    dK = {n: templates(S, n, 0.0)[2] for n in CAVITIES}
    pairs = {'presence (room)': dK['WA'], 'presence (L tunnel)': dK['WB'], 'location (6 m)': dK['WA6'] - dK['WA'],
             'shape (room or L tunnel)': dK['WB'] - dK['WA']}
    for name, K in pairs.items():
        A, B = inf.fisher_grid(K, g, S['f'])
        out[name] = {lv: {'deflection': float(a * np.sqrt(A)), 'tv': float(2 * ndtr(a * np.sqrt(A) / 2) - 1)}
                     for lv, a in LEVELS.items()}
    return out


def main():
    S = setup()
    params = {'benchmark_frozen_hash': m36.FROZEN_HASH, 'benchmark_worlds_hash': m36.WORLDS_HASH,
              'images_note': 'the images depend only on the worlds hash (worlds, lorry, solver, geometry), unchanged since '
                             'revision 2, whose fields they were read from', 'f_hz': S['f'], 'grounds': GROUNDS, 'levels': LEVELS,
              'taper_m': TAPER, 'top_share': TOP_SHARE, 'image_px': m07.SHAPE, 'geometry': m36.FROZEN['acquisition'],
              'microseisms': 'regional level, each ground its own realisation', 'amplified_levels_are_diagnostic': True}
    with Run(RID, 'The benchmark read: the published method and a detector told the lorry', params) as run:
        recs = []
        for seed in GROUNDS:
            for key in keys():
                t0 = time.time()
                recs.append(memo(RID, f'reading-{seed}-{key}', lambda: reading(S, seed, key), __file__,
                                 version=READINGS))
        rows, summary, sd = score(S, recs)
        pred = memo(RID, 'predicted', lambda: predicted(S), __file__, version='pred-v1')
        for lv in LEVELS:
            t, b, p = summary[lv]['told'], summary[lv]['published_blind'], summary[lv]['published_paired']
            print(f"  {lv}: told presence AUC {t['presence_auc_vs_controls']:.2f}, names {t['named_of_three']:.2f} "
                  f"(p {t['named_of_three_p']:.3g}), location {t['location_correct']}, shape {t['shape_correct']}; "
                  f"published blind AUC {b['presence_auc_vs_controls']:.2f}, names {b['named_of_three']:.2f} "
                  f"(p {b['named_of_three_p']:.3g}); paired names {p['named_of_three']:.2f}", flush=True)
        pos, real = summary['x1000'], summary['real']
        told_ok = pos['told']['named_of_three'] is not None and pos['told']['named_of_three'] >= 2 / 3
        fmt_p = lambda v: 'p < 0.001' if v is not None and v < 0.001 else f'p = {v:.2g}'
        finding = (
            f"The benchmark's worlds read as images (P2-36's solver fields, a lorry 15 m west at the FTA's truck-over-a-bump "
            f"level held at {S['f']:.0f} Hz, Giza's microseisms, {len(GROUNDS)} independent grounds). At the real level the "
            f"detector told the lorry separates cavities from controls with AUC {real['told']['presence_auc_vs_controls']:.2f} "
            f"and names the world in {100 * real['told']['named_of_three']:.0f}% (chance 33%); the published map read blind, "
            f"AUC {real['published_blind']['presence_auc_vs_controls']:.2f} and {100 * real['published_blind']['named_of_three']:.0f}%. "
            f"With the cavity's difference amplified 1,000 times (a diagnostic), the told detector names "
            f"{100 * pos['told']['named_of_three']:.0f}% ({fmt_p(pos['told']['named_of_three_p'])}, labels permuted within "
            f"grounds), places the room against the room 6 m east in {pos['told']['location_correct']} and tells room from "
            f"L tunnel in {pos['told']['shape_correct']}; the published map read blind names "
            f"{100 * pos['published_blind']['named_of_three']:.0f}% ({fmt_p(pos['published_blind']['named_of_three_p'])}), "
            f"separates cavities from controls with AUC {pos['published_blind']['presence_auc_vs_controls']:.2f} and centres "
            f"its top places {pos['published_blind']['median_centroid_error_m']:.1f} m from the footprint (controls "
            f"{summary['control_centroid_error_m']:.1f} m); its change against the no-cavity image, a diagnostic it never "
            f"has, names {100 * pos['published_paired']['named_of_three']:.0f}% (its presence statistic is not scored: the "
            f"synthetic twin has no receiver noise, so any difference shows, and it falls as the difference grows). Where "
            f"p is 1 the reader named the same world for every image of a ground, its answer unchanged by the cavity. "
            + ("The told detector recovers the imposed worlds at the positive control, so the experiment can reveal a real "
               "signal; what the published method does not read there it does not read in the image." if told_ok else
               "The told detector does not recover the worlds even at the positive control; the experiment is examined "
               "before any reader is judged.")
            + " Amplified levels are diagnostics, not physical predictions.")
        run.save({'summary': summary, 'rows': rows, 'null_sd_told': sd, 'predicted_told': pred, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'worker':
        worker(int(sys.argv[2]), int(sys.argv[3]))
    else:
        main()
