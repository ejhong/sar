"""P2-31 · The shape test: does one image, read by the published method or by a reference detector, recover where a
cavity is and what shape it has, even if its depth is wrong?

    uv run python experiments/p2_31_shape_test.py

The organising question has widened: does one complex SAR acquisition hold reproducible information about what is
underground, where it is, or its shape, even with depth uncertain? A method that drew a branching tunnel at one of two
possible depths would be a major success. So five claims are scored apart: presence, horizontal location, shape, relative
vertical geometry and absolute depth.

Scenes, generated physically and never through any inversion's depth-steering formula. The one-chamber bench's
limestone with, in turn: no cavity; an isolated room (6 m, centre 12 m down); an L-shaped tunnel (2.5 m square, arms of
18 m and 17 m, 10 m down); a branching (T) tunnel (a 28 m main run and a 15 m branch, 10 m down). Each cavity's imprint
on the ground's motion is the lab's elastic solver settling under a uniform horizontal strain (P2-28's loading, checked
against P2-04's). The unrelated ground is matched within each ground: three independent grounds, each with its own
speckle realisation and its own microseism realisation (Giza's measured level), shared by every layout and control on
that ground; the grounds are the evaluation units and the layouts on one ground are paired comparisons. The images come
from the lab's synthesizer on the real dwell (P2-07's chain) in complex128 (a first version synthesised in single
precision, whose rounding lies 20 to 40 times above the real-level imprints: its real-level row measured rounding).

Levels. The real one (the imprint as it is); and three amplified diagnostics, each cavity's imprint boosted to the same
largest phase, 0.3, 2 and 20 rad, labelled as such: physically shaped, not physically sized. At microseism frequencies
only the imprint's gradient survives in the image (P2-25), so 2 rad (millimetres of motion) sits near what any reader
needs (P2-27); 20 rad is the positive control. The information bound (P2-25,
P2-30) says the real level holds almost nothing for any reader; the diagnostics ask whether the experiment and each
reader can see a real signal where one is there.

Readers.
1. The published method (P2-07's pipeline, the paper's lambda_s): its focused power summed over depth, read blind from the
   single image as a plan map; and, as a diagnostic the method never has, its change against the same ground without the
   cavity (paired).
2. A reference detector told the shaking: the locally most powerful statistic for each candidate layout (the score, as
   P2-25 checks it), q_c(x) = Im z conj(z_eps_c) per strain component correlated with each layout's imprint templates;
   presence and shape by the largest normalised score among the three layouts at the true placement (it is told where
   to look), location by where each layout's statistic peaks when every placement is scanned.
   Normalised by the spread of the same maps over the images with no cavity and the motionless copies.
Controls: no cavity, and a motionless copy of the ground, for every realisation, read identically.

Scoring, stated before the run. Presence: each image's largest normalised response, cavity scenes against the controls.
Horizontal location: the distance from the reader's strongest place (the centroid of its top 5% for the published maps)
to the nearest point of the true plan footprint; added after the first run, because the tunnels' footprints reach
within a few metres of almost any central place (chance 3.4 m), the distance from that centroid to the footprint's
centroid, with its chance from the controls. Shape: the layout whose plan footprint the reader's map correlates with
best (chance one in three), the controls scored the same way. Extra structures: connected regions of the top 5% that
touch no true footprint within 5 m, counted as false positives. Depth, for the published method: one predeclared family of
global transformations (lambda_s of 0.24 or 0.48 m, and the axis's repeat and mirror branches), the same for every image;
the change's peak depth scored after the best member of that family, the controls given the same freedom. Relative
vertical geometry is not testable with these single-level layouts and is not claimed. If the reference detector cannot
recover the imposed layouts at the strongest diagnostic, the forward model or the acquisition's sensitivity is examined
before any reader is blamed.
"""
import copy
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import binary_dilation, gaussian_filter, label
from scipy.signal import fftconvolve

from katabasis.ambient.field import microseisms, rayleigh_hv
from katabasis.runs import Run, memo
from sarsim import synthesize
from sarsim import information as inf
from sarsim.acquisition import DwellGeometry

RID = 'p2_31_shape_test'
SITES = Path(__file__).resolve().parents[2] / 'sites'
BOXES = {
    'room': [((0.0, 0.0, -12.0), (6.0, 6.0, 6.0))],
    'L_tunnel': [((-5.0, 0.0, -10.0), (18.0, 2.5, 2.5)), ((2.75, 7.375, -10.0), (2.5, 17.25, 2.5))],
    'branching_tunnel': [((0.0, 0.0, -10.0), (28.0, 2.5, 2.5)), ((0.0, -8.625, -10.0), (2.5, 14.75, 2.5))],
}
SHAPES = list(BOXES)
SEEDS = (71, 171, 271)
LEVELS = {'real': None, 'diagnostic_0.3_rad': 0.3, 'diagnostic_2_rad': 2.0, 'diagnostic_20_rad': 20.0}
READINGS = 'readings-v2'      # the version of the scenes and readers; bump it if either changes (v2: complex128, a field per ground)
FIELD_SEED = 1000             # each ground's microseism realisation: default_rng(FIELD_SEED + seed)
TAPER = (45.0, 54.0)
TOP_SHARE = 0.05
EXTRA_M = 5.0
LAM_S_FAMILY = (0.24, 0.48)
DEPTH_FAMILY_SIZE = 2 * len(LAM_S_FAMILY) * 3        # both branches, three repeats, each lambda_s


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def media_boxes(m28, boxes):
    from katabasis.compose import Grid, load_site, voxelise
    from katabasis.compose.site import parse_site
    from katabasis.seismic.elastic3d import Medium
    site = load_site('bench-void')
    raw = copy.deepcopy(site.raw)
    f0 = raw['features'][0]
    feats = []
    for k, (c, s) in enumerate(boxes):
        f = copy.deepcopy(f0)
        f['id'] = f'part_{k}'
        f['shape'] = dict(f['shape'], centre=list(c), size=list(s))
        f['fill'] = 'air'
        feats.append(f)
    raw['features'] = feats
    g = Grid.covering((-m28.HALF_DOMAIN, m28.HALF_DOMAIN), (-m28.HALF_DOMAIN, m28.HALF_DOMAIN), (m28.BOTTOM, m28.TOP), m28.H)
    return g, Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))


def layout_kernels(m28, name, boxes, host):
    from katabasis.seismic.arrays import surface_grid
    g, full = media_boxes(m28, boxes)
    rec = surface_grid(g, full.solid, m28.HALF_REC, m28.SPACING)
    K1, i1 = m28.cached(f'shape_{name}', lambda: m28.settle(full, rec, host, f'shape_{name}'))
    K0, _ = m28.cached('without', lambda: m28.settle(m28.media(None)[1], rec, host, 'without'))
    return {'xy': K1['xy'], **{k: K1[k] - K0[k] for k in ('Kxx', 'Kyy', 'Kxy')}}, i1


def interpolators(K, los):
    xy = K['xy']
    xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
    ix, iy = np.searchsorted(xs, xy[:, 0]), np.searchsorted(ys, xy[:, 1])
    out = []
    for name in ('Kxx', 'Kyy', 'Kxy'):
        M = np.zeros((len(xs), len(ys)))
        M[ix, iy] = K[name] @ los
        out.append(RegularGridInterpolator((xs, ys), M, bounds_error=False, fill_value=0.0))
    return out


def plan_mask(boxes, east, north):
    inside = np.zeros(east.shape, bool)
    for (c, s) in boxes:
        inside |= (np.abs(east - c[0]) <= s[0] / 2) & (np.abs(north - c[1]) <= s[1] / 2)
    return inside


def to_site(g, x, y):
    a, r = g.along_track_en, g.ground_range_en
    return x * a[0] + y * r[0], x * a[1] + y * r[1]


def score_fields(img, g, taus):
    """q_c(x) = (2 / s2) Im z(x) conj(z_c(x)), z_c the image re-formed with each Doppler bin weighted by strain
    component c's history at that bin's time: the per-pixel field whose inner product with a layout's phase template is
    that layout's locally most powerful statistic (sarsim.information.score)."""
    Y = np.fft.fft2(img)
    in_a, in_r, nu = inf.band_masks(g, img.shape)
    band = in_a[:, None] & in_r[None, :]
    s2 = float(np.mean(np.abs(Y[band]) ** 2)) / img.size
    z = np.fft.ifft2(Y * band)
    t = g.nu_to_time(nu)
    out = []
    for tau in taus:
        zc = np.fft.ifft2(np.asarray(tau(t))[:, None] * Y * band)
        out.append(2.0 / s2 * np.imag(z * np.conj(zc)))
    return out


def corr_map(q, K):
    """The statistic for the template K placed at every pixel: sum_c (q_c correlated with K_c)."""
    return sum(fftconvolve(qc, Kc[::-1, ::-1], mode='same') for qc, Kc in zip(q, K))


def top_centroid(M, E, N, share=TOP_SHARE):
    thr = np.quantile(M, 1 - share)
    sel = M >= thr
    return float(E[sel].mean()), float(N[sel].mean()), sel


def distance_to(mask, E, N, e, n):
    if not mask.any():
        return float('nan')
    return float(np.min(np.hypot(E[mask] - e, N[mask] - n)))


def extras(sel, truth, step_m):
    grow = binary_dilation(truth, iterations=max(1, int(round(EXTRA_M / step_m)))) if truth.any() else truth
    lab, n = label(sel)
    return int(sum(1 for k in range(1, n + 1) if not (grow & (lab == k)).any()))


def main():
    m07 = module('p207', Path(__file__).with_name('p2_07_whole_chain.py'))
    m28 = module('p228', Path(__file__).with_name('p2_28_two_depths.py'))
    m07.TAPER = TAPER
    g = DwellGeometry.from_record('giza-20250827')
    los = np.asarray(g.los_enu)
    k0 = 4 * np.pi / g.lam
    amb = json.loads((SITES / 'ambient.json').read_text())
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    rock = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    vp, vs, rho = val('vp'), val('vs'), val('rho')
    host = (rho * vp ** 2 - 2 * rho * vs ** 2, rho * vs ** 2)
    level = amb['regional']['microseism_vertical_0.1_0.3_hz']['value']
    speed = amb['regional']['microseism_phase_speed']['value']
    params = {'layouts': {k: [{'centre': c, 'size': s} for c, s in v] for k, v in BOXES.items()}, 'seeds': SEEDS,
              'levels': LEVELS, 'taper_m': TAPER, 'top_share': TOP_SHARE, 'extra_structure_margin_m': EXTRA_M,
              'depth_family': {'lam_s_m': LAM_S_FAMILY, 'branches': 'the axis repeat and mirror'},
              'geometry': 'giza-20250827', 'image_px': m07.SHAPE, 'microseisms_m_s': level,
              'imprint': 'solver settling under uniform strain (P2-28)', 'amplified_levels_are_diagnostic': True}
    with Run(RID, 'The shape test: location and shape, with depth allowed to be wrong', params) as run:
        t0 = time.time()
        kern = {n: layout_kernels(m28, n, b, host) for n, b in BOXES.items()}
        print(f"  kernels: {time.time() - t0:.0f} s", flush=True)
        interp = {n: interpolators(k[0], los) for n, k in kern.items()}
        f = np.linspace(0.1, 0.3, 21)
        # each ground its own microseism realisation, so the grounds are independent scenes (the evaluation units)
        fields = {}
        def field_for(seed):
            if seed not in fields:
                fields[seed] = microseisms(f, np.full(len(f), level ** 2 / 0.2), np.random.default_rng(FIELD_SEED + seed),
                                           per_bin=4, speed=speed, hv=rayleigh_hv(vp, vs))
            return fields[seed]
        def taus_for(seed):
            fld = field_for(seed)
            eps = lambda tt, c: (lambda G: [G[0, 0], G[1, 1], 0.5 * (G[0, 1] + G[1, 0])][c])(fld.centre_series(tt)['g'])
            return [lambda tt, c=c: eps(tt, c) for c in range(3)]
        SH = m07.SHAPE
        # the image's pixels in the site frame, and the pipeline's positions
        xa = (np.arange(SH[0]) - SH[0] // 2) * g.dx
        xg = (np.arange(SH[1]) - SH[1] // 2) * g.dr / np.sin(g.theta)
        XA, XG = np.meshgrid(xa, xg, indexing='ij')
        PE, PN = to_site(g, XA, XG)
        rows = np.arange(48, SH[0] - 48, m07.STRIDE[0])
        cols = np.arange(32, SH[1] - 32, m07.STRIDE[1])
        RR, CC = np.meshgrid(rows, cols, indexing='ij')
        GE, GN = PE[RR, CC], PN[RR, CC]
        step_pub = float(np.hypot(*(np.array(to_site(g, m07.STRIDE[0] * g.dx, 0.0)))))
        truth_pix = {n: plan_mask(b, PE, PN) for n, b in BOXES.items()}
        truth_pub = {n: plan_mask(b, GE, GN) for n, b in BOXES.items()}
        # each layout's phase templates per unit strain on the pixels (the reference detector's), and its unit phase
        tmpl, unit = {}, {}
        t_fine = np.linspace(-12.0, 12.0, 2401)
        for n in SHAPES:
            L = m07.Shaking(g, field_for(SEEDS[0]), interp[n], 1.0).lk(PE.ravel(), PN.ravel())   # [pixels, 3] LOS per unit strain
            tmpl[n] = [(-k0 * L[:, c]).reshape(SH) for c in range(3)]
            for sd in SEEDS:
                unit[(sd, n)] = m07.Shaking(g, field_for(sd), interp[n], 1.0).imprint_phase(t_fine)
        # a compact template for location: the room's, cut to 24 m about its centre
        room_t = [T[SH[0] // 2 - int(24 / g.dx): SH[0] // 2 + int(24 / g.dx), SH[1] // 2 - int(24 * np.sin(g.theta) / g.dr): SH[1] // 2 + int(24 * np.sin(g.theta) / g.dr)]
                  for T in tmpl['room']]

        def reading(seed, key):
            """One image, made and read: kept on disk so a run cut short resumes (katabasis.runs.memo)."""
            scat = m07.scene(g, np.random.default_rng(seed))
            field = field_for(seed)
            # complex128: the real-level imprints (1e-9 rad) lie far below single precision's rounding
            if key == 'none':
                img = synthesize(scat, g, SH, motion=m07.Shaking(g, field, interp['room'], 0.0), dtype=np.complex128)
            elif key == 'motionless':
                img = synthesize(scat, g, SH, dtype=np.complex128)
            else:
                n, lv = key.split('|')
                G = 1.0 if LEVELS[lv] is None else LEVELS[lv] / unit[(seed, n)]
                img = synthesize(scat, g, SH, motion=m07.Shaking(g, field, interp[n], G), dtype=np.complex128)
            T, z = m07.run_method(img, g, RR.ravel(), CC.ravel())
            q = score_fields(img, g, taus_for(seed))
            # each layout's statistic at every placement; its value at the true placement (the scene's centre), and
            # the spread of its values over placements (used from the controls as the null)
            stats, spread, peak_at = {}, {}, {}
            for n in SHAPES:
                M = corr_map(q, tmpl[n])
                stats[n] = float(M[SH[0] // 2, SH[1] // 2])
                inner = M[SH[0] // 8: -SH[0] // 8, SH[1] // 8: -SH[1] // 8]
                spread[n] = float(np.std(inner))
                ia, ir = np.unravel_index(np.argmax(inner), inner.shape)
                peak_at[n] = (float(PE[ia + SH[0] // 8, ir + SH[1] // 8]), float(PN[ia + SH[0] // 8, ir + SH[1] // 8]))
            loc = corr_map(q, room_t).astype(np.float32)
            return {'seed': seed, 'key': key, 'scene': key.split('|')[0],
                    'level': key.split('|')[1] if '|' in key else 'control', 'stats': stats, 'spread': spread,
                    'peak_at': peak_at, 'loc': loc, 'blind': T.sum(axis=1).reshape(RR.shape), 'T': T, 'z': z}

        records = []
        for seed in SEEDS:
            keys = ['none', 'motionless'] + [f'{n}|{lv}' for n in SHAPES for lv in LEVELS]
            T_none = None
            for key in keys:
                t1 = time.time()
                r = memo(RID, f'reading-{seed}-{key}', lambda: reading(seed, key), __file__, version=READINGS)
                records.append(r)
                if key == 'none':
                    T_none = r['T']
                print(f"  seed {seed} {key}: {time.time() - t1:.0f} s", flush=True)
            for r in records:
                if r['seed'] == seed:
                    r['paired'] = np.abs(r['T'] - T_none).sum(axis=1).reshape(RR.shape)
                    r['paired_T'] = r['T'] - T_none

        # normalisation from the controls: spread of each layout's statistic over shifted placements in control images
        ctrl = [r for r in records if r['level'] == 'control']
        sd_loc = float(np.std(np.concatenate([r['loc'][200:-200, 20:-20].ravel() for r in ctrl])))
        sd_stat = {n: float(np.sqrt(np.mean([r['spread'][n] ** 2 for r in ctrl]))) for n in SHAPES}
        pub_norm = lambda M: (M - np.median(M)) / (np.std(M) + 1e-30)
        def depth_after_allowance(r, name):
            zt = BOXES[name][0][0][2] * -1.0
            foot = truth_pub[name].ravel()
            prof = np.abs(r['paired_T'])[foot].mean(axis=0) if foot.any() else np.zeros_like(r['z'])
            zpk = float(r['z'][int(np.argmax(prof))])
            Z = float(r['z'][-1])
            cands = []
            for lam in LAM_S_FAMILY:
                a = lam / 0.48
                for k in range(3):
                    cands += [a * (2 * k * Z + zpk), a * (2 * (k + 1) * Z - zpk)]
            return {'peak_axis_m': zpk, 'true_m': zt, 'best_error_after_allowance_m': float(min(abs(c - zt) for c in cands)),
                    'family_size': len(cands)}

        results = {}
        for lv in list(LEVELS) + ['control']:
            rs = [r for r in records if r['level'] == lv]
            out = []
            for r in rs:
                truth_s = r['scene'] if r['scene'] in SHAPES else None
                zst = {n: r['stats'][n] / sd_stat[n] for n in SHAPES}
                ref_shape = max(zst, key=zst.get)
                li = np.unravel_index(np.argmax(r['loc']), r['loc'].shape)
                # location: where the chosen layout's statistic peaks, against where the layout truly lies (the origin);
                # for the controls, the same distance is chance
                pe, pn = r['peak_at'][ref_shape]
                row = {'seed': r['seed'], 'scene': r['scene'], 'reference': {
                    'presence': max(zst.values()), 'z_by_layout': zst, 'shape': ref_shape,
                    'shape_correct': (ref_shape == truth_s) if truth_s else None,
                    'location_error_m': float(np.hypot(pe, pn)),
                    'room_template_peak_over_null': float(r['loc'][li] / sd_loc)}}
                for reader, M in (('published_blind', r['blind']), ('published_paired', r.get('paired'))):
                    if M is None:
                        continue
                    Mn = pub_norm(gaussian_filter(M, 1.0))
                    ce, cn, sel = top_centroid(Mn, GE, GN)
                    cors = {n: float(np.nan_to_num(np.corrcoef(Mn.ravel(), gaussian_filter(truth_pub[n].astype(float), 1.0).ravel())[0, 1], nan=-1.0))
                            for n in SHAPES}
                    best = max(cors, key=cors.get)
                    row[reader] = {'presence': float(Mn.max()), 'correlation_by_layout': cors, 'shape': best,
                                   'shape_correct': (best == truth_s) if truth_s else None,
                                   'location_error_m': distance_to(truth_pub[truth_s], GE, GN, ce, cn) if truth_s else None,
                                   'centroid_error_m': float(np.hypot(ce - GE[truth_pub[truth_s]].mean(), cn - GN[truth_pub[truth_s]].mean())) if truth_s else None,
                                   'extra_structures': extras(sel, truth_pub[truth_s] if truth_s else np.zeros_like(sel), step_pub)}
                # depth, for the published change (paired), under the predeclared family
                if r.get('paired_T') is not None and truth_s:
                    row['published_depth'] = depth_after_allowance(r, truth_s)
                elif r['key'] == 'motionless':
                    # the controls given the same freedom: the motionless copy's change over each layout's footprint
                    row['published_depth_control'] = [depth_after_allowance(r, n)['best_error_after_allowance_m'] for n in SHAPES]
                out.append(row)
            results[lv] = out

        def scene_permutation_p(rows_lv, reader):
            """P(at least the observed number named correctly) when the true layouts are permuted within each scene, over
            every arrangement: the scenes (seeds) are the independent units."""
            import itertools
            by_scene = {}
            for x in rows_lv:
                if reader in x and x[reader].get('shape') is not None:
                    by_scene.setdefault(x['seed'], []).append((x['scene'], x[reader]['shape']))
            if not by_scene:
                return None
            observed = sum(t == n for v in by_scene.values() for t, n in v)
            counts = []
            for v in by_scene.values():
                truths, named = [t for t, _ in v], [n for _, n in v]
                counts.append([sum(t == n for t, n in zip(perm, named)) for perm in itertools.permutations(truths)])
            totals = [sum(c) for c in itertools.product(*counts)]
            return float(np.mean([t >= observed for t in totals]))

        def summarise(rows_lv, reader, field_):
            v = [x[reader][field_] for x in rows_lv if reader in x and x[reader].get(field_) is not None]
            return v
        summary = {}
        ctrl_rows = results['control']
        for lv in LEVELS:
            s = {}
            cav = [x for x in results[lv] if x['scene'] in SHAPES]
            for reader in ('reference', 'published_blind', 'published_paired'):
                pres_c = [x[reader]['presence'] for x in cav if reader in x]
                pres_n = [x[reader]['presence'] for x in ctrl_rows if reader in x]
                auc = float(np.mean([[1.0 if a > b else 0.5 if a == b else 0.0 for b in pres_n] for a in pres_c])) if pres_c and pres_n else None
                sc = summarise(cav, reader, 'shape_correct')
                le = summarise(cav, reader, 'location_error_m')
                s[reader] = {'presence_auc_vs_controls': auc, 'shape_accuracy': float(np.mean(sc)) if sc else None,
                             'shape_chance': 1 / len(SHAPES), 'median_location_error_m': float(np.median(le)) if le else None,
                             # a null significance test with the scenes as units (the layouts in a scene share its
                             # ground): the true labels permuted within each scene, every arrangement; not a confidence
                             # interval for predictive performance
                             'shape_p_permutation': scene_permutation_p(cav, reader),
                             'images': len(sc), 'scenes': len({x['seed'] for x in cav})}
                if reader == 'reference':
                    # presence and shape are read at the true placement (told where to look); location by scanning
                    s[reader]['median_presence_z_at_true_place'] = float(np.median(pres_c)) if pres_c else None
                else:
                    ex = summarise(cav, reader, 'extra_structures')
                    s[reader]['mean_extra_structures'] = float(np.mean(ex)) if ex else None
                    cen = summarise(cav, reader, 'centroid_error_m')
                    s[reader]['median_centroid_error_m'] = float(np.median(cen)) if cen else None
            dep = [x['published_depth']['best_error_after_allowance_m'] for x in cav if 'published_depth' in x]
            s['published_depth_error_after_allowance_m'] = float(np.median(dep)) if dep else None
            summary[lv] = s
        # chance: the controls' top places against each layout's footprint
        chance_loc = []
        for r in ctrl:
            ce, cn, _ = top_centroid(pub_norm(gaussian_filter(r['blind'], 1.0)), GE, GN)
            chance_loc += [distance_to(truth_pub[n], GE, GN, ce, cn) for n in SHAPES]
        summary['control_location_error_to_layouts_m'] = float(np.median(chance_loc))
        summary['control_depth_error_after_allowance_m'] = float(np.median([e for x in ctrl_rows for e in x.get('published_depth_control', [])]))
        # chance for the centroid: the controls' top places against each layout's footprint centroid; for the paired map
        # the motionless copy's change (the no-cavity image's change against itself is identically zero)
        def chance_centroid(maps_):
            out_ = []
            for M in maps_:
                ce, cn, _ = top_centroid(pub_norm(gaussian_filter(M, 1.0)), GE, GN)
                out_ += [float(np.hypot(ce - GE[truth_pub[n]].mean(), cn - GN[truth_pub[n]].mean())) for n in SHAPES]
            return float(np.median(out_))
        summary['control_centroid_error_blind_m'] = chance_centroid([r['blind'] for r in ctrl])
        summary['control_centroid_error_paired_m'] = chance_centroid([r['paired'] for r in ctrl if r['key'] == 'motionless'])
        summary['control_reference_peak_distance_m'] = float(np.median([x['reference']['location_error_m'] for x in ctrl_rows]))
        np.savez_compressed(Path(run.dir) / 'maps.npz',
                            **{f"{r['seed']}_{r['key'].replace('|', '_')}_blind": r['blind'].astype(np.float32) for r in records},
                            **{f"{r['seed']}_{r['key'].replace('|', '_')}_paired": r['paired'].astype(np.float32) for r in records if 'paired' in r},
                            east=GE, north=GN, **{f'truth_{n}': truth_pub[n] for n in SHAPES})
        pos = summary['diagnostic_20_rad']
        mid = summary['diagnostic_2_rad']
        real = summary['real']
        ok_ref = (pos['reference']['shape_accuracy'] or 0) >= 2 / 3
        rp, rm = pos['reference'], mid['reference']
        bp, pp = pos['published_blind'], pos['published_paired']
        mp = mid['published_paired']
        finding = (
            f"Scenes generated through the forward model (the lab's solver for each cavity's imprint, the synthesizer for "
            f"the image), with the unrelated ground and shaking matched. At the positive control (each imprint boosted to "
            f"20 rad) the reference detector told the shaking and read where the cavity lies names the layout in "
            f"{100 * rp['shape_accuracy']:.0f}% of images (chance {100 / len(SHAPES):.0f}%); scanning every placement, its "
            f"strongest response lies a median {rp['median_location_error_m']:.1f} m from the layout's centre (the controls' "
            f"{summary['control_reference_peak_distance_m']:.0f} m). At 2 rad, near what any reader needs at these "
            f"frequencies (P2-27), it still names {100 * rm['shape_accuracy']:.0f}% when told where to look, but its "
            f"response there stands a median {rm['median_presence_z_at_true_place']:.1f} spreads above the placements' "
            f"scatter and a scan of the whole image peaks {rm['median_location_error_m']:.0f} m away: not knowing where "
            f"costs more than knowing what. The published method's plan map, read blind, names the layout in "
            f"{100 * bp['shape_accuracy']:.0f}% at 20 rad (p = {bp['shape_p_permutation']:.2f}, labels permuted within the three "
            f"scenes, which are the independent units; a null test, not a confidence interval) "
            f"and centres its top places {bp['median_centroid_error_m']:.1f} m "
            f"from the layout's centre (the controls' {summary['control_centroid_error_blind_m']:.1f} m); its change "
            f"against the same ground without the cavity, a diagnostic it never has, {100 * pp['shape_accuracy']:.0f}% and "
            f"{pp['median_centroid_error_m']:.1f} m (the motionless copy's change, "
            f"{summary['control_centroid_error_paired_m']:.1f} m); at 2 rad that change names {100 * mp['shape_accuracy']:.0f}% "
            f"(permutation p = {mp['shape_p_permutation']:.2f}, before allowing for the dozen comparisons) and centres "
            f"{mp['median_centroid_error_m']:.1f} m away: the picture changes near where the image changes (P2-28). The blind "
            f"map's largest value separates cavities from controls with AUC {bp['presence_auc_vs_controls']:.2f} (0.5 is chance) "
            f"at 20 rad. At the real level the reference detector names "
            f"{100 * real['reference']['shape_accuracy']:.0f}% and the published blind map "
            f"{100 * real['published_blind']['shape_accuracy']:.0f}%, chance. After the predeclared depth family the "
            f"published change's peak lies {pos['published_depth_error_after_allowance_m']:.1f} m from the true depth at 20 rad "
            f"and the motionless copy's {summary['control_depth_error_after_allowance_m']:.1f} m: with {DEPTH_FAMILY_SIZE} "
            f"candidate depths, closeness alone says little. "
            + ("The reference detector recovers the imposed layouts at the positive control, so the experiment can reveal "
               "a real signal." if ok_ref else "The reference detector does not recover the imposed layouts even at the positive control: "
               "the forward model or the acquisition's sensitivity must be examined before any reader is judged.")
            + " Amplified levels are diagnostics, not physical predictions.")
        run.save({'summary': summary, 'results': results, 'unit_gain_phase_rad': {f'{k[0]}|{k[1]}': v for k, v in unit.items()}, 'kernel_runs': {n: k[1] for n, k in kern.items()},
                  'normalisation': {'location_sd': sd_loc, 'layout_sd': sd_stat}, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
