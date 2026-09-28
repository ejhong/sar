"""P2-13 · Stricter gates: motionless images pass them; a planted vibration is not recovered.

    uv run --with scikit-image python experiments/p2_13_gates.py

Later reconstructions of the method add selection gates between the shifts and the depth focus: keep only windows of
50 pairs whose two shift components trace an ellipse at one of modes 1 to 10, with a minimum fit and shape; require ten
successive windows to agree on the mode; keep the family inside the central pairs; and require four or five neighbouring
pixels to pass at the same window. The depth focus then runs only where everything passed. A public implementation of
this gated pipeline, with a frozen profile for the ICEYE X13 pass of 15 July 2022 (317 pairs, masks 32,330 Hz wide and
404 Hz apart, registration to 0.01 px, lambda_s 0.2406 m, depths 0 to 300 m), is run here unchanged: its source is
imported as published and its hash recorded; only the image it reads is replaced.

Lines: four sides round Khafre's faces at 70.5 m above the base (the profile's height), placed through the product's RPC
over the site's surface with the residual offset fitted against predicted brightness (sarsim.ortho, within 80 m).

Five cases, the same code:
  real     the product's crop, as it is;
  twin     a motionless twin (sarsim.looks.motionless_twin): the crop's brightness pattern and spectrum, fresh speckle,
           nothing moving and nothing inside;
  planted  the real crop with two stretches of the faces vibrating along the line of sight at 2 mm/s (some 40,000 times
           Giza's microseisms): 0.26 Hz on the north face and 3.66 Hz on the south face, frequencies proposed for a 648 m
           column of air and for the Grand Gallery's length; put into the image itself (sarsim.looks.inject_region_motion);
  desert   the real crop, with the same loops laid on open ground south of Khafre (P2-12's in-block control, chosen from
           the image before any method was run there), where nothing is claimed.

Scored: the shifts (how many sit at exactly zero, their spread), how many windows pass each gate, where features survive,
each feature's winning mode against the depth of its best score (depth = mode x the depth one cycle spans in a window, or
its mirror in the repeat), and how far the planted vibrations move the shifts against what a tracker following them fully
would report. From the geometry: the frequencies the gate's modes stand for, how much of a motion at them a pair's image
keeps, and the shape a vibrating point's shifts trace (its azimuth shift follows the line-of-sight velocity, its range shift
the displacement) against the shape the gate keeps.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from katabasis.compose import load_site
from katabasis.export.sites import scene as site_scene
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.dwell import DwellProduct
from sarsim.looks import inject_region_motion, motionless_twin
from sarsim.ortho import frame_to_lla, multilook, predicted_brightness, project, register, surface

RID = 'p2_13_gates'
ACQ = 'giza-20220715'
PRODUCT = Path.home() / 'tmp/sar/giza2/ICEYE_X13_SLC_SLED_868226_20220715T235744.h5'
GATED = Path.home() / 'tmp/sar/biondi_v18/pkg/Biondi_Protocol_v1.8'
MODULE = GATED / 'biondi_tomography_v1_8.py'
PROFILE = GATED / 'profiles/iceye_x13_khafre_w50.json'
GEOID_M = 15.5
FIT_WINDOW_M = 80.0
FACE_HEIGHT_M = 70.5                 # the profile's own height
HEIGHTS = [20.0, 35.0, 50.0, 60.0, 70.5, 80.0, 95.0]
V_AMP = 2e-3
PLANTS = [('North', 0.26), ('South', 3.66)]
MIDDLE = (30, 71)
MARGIN = 40
TWIN_SEEDS = [101, 102]
DESERT = (-100.0, -640.0)            # P2-12's in-block control: open ground south of Khafre, (x, y) in the site frame
SIDES = ['North', 'East', 'South', 'West']
CASES = ['real'] + [f'twin{seed}' for seed in TWIN_SEEDS] + ['planted', 'desert']


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_gated():
    spec = importlib.util.spec_from_file_location('gated', MODULE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fit_offset(p, g, sc):
    x = np.arange(-600, 600, 2.0)
    X, Y = np.meshgrid(x, x)
    B, Z = predicted_brightness(sc, X, Y, g.los_enu, 2.0)
    rr, cc = project(p, sc, X, Y, Z, GEOID_M)
    la = int(round(1.0 / g.dx))
    lr = max(1, int(round(np.sin(g.theta) / g.dr)))
    obs, origin = multilook(p, int(rr.min()) - 2000, int(rr.max()) + 2000, int(cc.min()) - 300, int(cc.max()) + 300, la, lr)
    window = (int(FIT_WINDOW_M / g.dx), int(FIT_WINDOW_M * np.sin(g.theta) / g.dr))
    return register(obs, origin, la, lr, rr, cc, B, max_shift_px=window)


def face_geometry(p, sc, shift, height, at=None):
    """Khafre's four faces at a height above its base; with `at`, the same square laid on the ground centred there."""
    kh = next(s for s in sc['structures'] if s['id'] == 'khafre')['shape']
    cx, cy, zb = kh['centre']
    half = kh['base'] / 2
    h = half - height / (kh['height'] / half)
    if at is not None:
        cx, cy = at
    corners = {'NW': (cx - h, cy + h), 'NE': (cx + h, cy + h), 'SE': (cx + h, cy - h), 'SW': (cx - h, cy - h)}

    def px(pt):
        x, y = np.array([pt[0]]), np.array([pt[1]])
        z = np.array([zb + height]) if at is None else surface(sc, x, y)
        r, c = project(p, sc, x, y, z, GEOID_M, shift)
        return [float(c[0]), float(r[0])]

    def ll(pt):
        a, b = frame_to_lla(sc['frame']['origin'], pt[0], pt[1])
        return [float(a), float(b)]
    title = (f"Khafre's four faces at {height:g} m above the base" if at is None
             else f"The loop of Khafre's faces at {height:g} m, laid on open ground at ({cx:g}, {cy:g}) m")
    geom = {'title': title,
            'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
            'projection': {'height_m': FACE_HEIGHT_M}, 'sides': {}}
    for name, (a, b) in zip(SIDES, [('NW', 'NE'), ('NE', 'SE'), ('SE', 'SW'), ('SW', 'NW')]):
        geom['sides'][name] = {'P0_label': a, 'P_last_label': b, 'geographic_endpoint_lat_lon': [ll(corners[a]), ll(corners[b])],
                               'native_endpoints_col_row': [px(corners[a]), px(corners[b])]}
    return geom, float(h)


def plant(img, tracks, row0, col0, g):
    out = img.copy()
    regions = []
    for side, f in PLANTS:
        tgt = tracks[side]['target'][MIDDLE[0]:MIDDLE[1]]
        r0, r1 = int(tgt[:, 0].min()) - MARGIN - row0, int(tgt[:, 0].max()) + MARGIN - row0
        c0, c1 = int(tgt[:, 1].min()) - MARGIN - col0, int(tgt[:, 1].max()) + MARGIN - col0
        w = np.zeros(out.shape[0])
        w[r0:r1] = 1.0
        amp = V_AMP / (2 * np.pi * f)
        out[:, c0:c1] = inject_region_motion(out[:, c0:c1], g, w, lambda t, a=amp, fr=f: a * np.sin(2 * np.pi * fr * t))
        regions.append({'side': side, 'f_hz': f, 'displacement_amp_m': amp})
    return out, regions


def ideal_pair_response_px(g, f, mask_hz, bshift_hz):
    """What a tracker following a point fully would see between a pair's two images: each image is displaced by
    2 V_g / (lambda |Ka|) times the velocity averaged over its span (mask / |Ka| seconds, gain sinc(f T)); the two spans
    are bshift / |Ka| apart. Amplitude in native azimuth pixels."""
    T = mask_hz / abs(g.Ka_signed)
    dt = bshift_hz / abs(g.Ka_signed)
    scale = 2 * g.V / (g.lam * abs(g.Ka_signed))
    return float(scale * V_AMP * abs(np.sinc(f * T)) * 2 * abs(np.sin(np.pi * f * dt)) / g.dx)


def runs_of(feats):
    """Contiguous runs of passing positions on one side: what the method would call one feature."""
    runs = []
    for f in sorted(feats, key=lambda f: f['position']):
        if runs and f['position'] == runs[-1][-1]['position'] + 1:
            runs[-1].append(f)
        else:
            runs.append([f])
    return runs


def analyse(audit, kz_per_side):
    a = audit
    out = {'sides': {}}
    for s, side in enumerate(SIDES):
        key = side.lower()
        Y = a['Y'][s]
        F = a['focus_p4'][s]
        pos = np.flatnonzero(np.isfinite(F).any(axis=1))
        dk = float(np.median(np.abs(np.diff(kz_per_side[s]))))
        zres = 2 * np.pi / (50 * dk)
        zrep = 2 * np.pi / dk
        feats = []
        for p in pos:
            i = int(np.nanargmax(F[p]))
            st = int(a['winning_start_p4'][s][p, i])
            m = int(a[f'{key}_window_mode'][p, st]) if st >= 0 else -1
            feats.append({'position': int(p), 'mode': m, 'best_depth_m': float(a['z_m'][i]), 'score': float(F[p, i]),
                          'mode_depth_m': m * zres, 'mirror_depth_m': zrep - m * zres})
        out['sides'][side] = {'exact_zero': float(np.mean(Y == 0)), 'within_one_step': float(np.mean(np.abs(Y) <= 0.0100001)),
                              'rms_px': float(np.sqrt(np.nanmean(Y ** 2))),
                              'window_passes': int(a[f'{key}_window_pass'].sum()),
                              'track10_families': int(a[f'{key}_family_pass'].sum()),
                              'p4_cells': int(a['p4_window_gate'][s].sum()), 'features': feats,
                              'runs': [[f['position'] for f in r] for r in runs_of(feats)]}
        out['depth_per_cycle_m'] = zres
        out['repeat_m'] = zrep
    out['feature_positions'] = int(sum(len(v['features']) for v in out['sides'].values()))
    out['feature_runs'] = int(sum(len(v['runs']) for v in out['sides'].values()))
    return out


def plural(n, one, many=None):
    return f"{n} {one if n == 1 else (many or one + 's')}"


def main():
    g = DwellGeometry.from_record(ACQ)
    sc = site_scene(load_site('giza'))
    p = DwellProduct(PRODUCT)
    profile = json.loads(PROFILE.read_text())
    bank = profile['processing']['filter_bank']
    gates = profile['gates']
    params = {'acquisition': ACQ, 'gated_module_sha256': sha256(MODULE), 'gated_profile_sha256': sha256(PROFILE),
              'profile': profile['profile_name'], 'face_height_m': FACE_HEIGHT_M, 'geoid_m': GEOID_M,
              'fit_window_m': FIT_WINDOW_M, 'planted': PLANTS, 'planted_v_m_s': V_AMP, 'twin_seeds': TWIN_SEEDS,
              'heights_m': HEIGHTS, 'desert_centre_m': DESERT,
              'bank': {k: bank[k] for k in ('support_hz', 'mask_width_hz', 'bshift_hz', 'k_leap_hz', 'pair_count')}}
    with Run(RID, 'Stricter gates: motionless images pass them; a planted vibration is not recovered', params) as run:
        shift, corr = fit_offset(p, g, sc)
        gdir = Path(__file__).resolve().parents[1] / 'data' / RID       # the gated pipeline's own outputs: large, not in git
        if gdir.exists():
            import shutil
            shutil.rmtree(gdir)
        gdir.mkdir(parents=True)
        gated = load_gated()
        original = gated.load_source_and_preflight
        sets, plants = [], []
        for height in HEIGHTS:
            hdir = gdir / f'h{height:g}'
            hdir.mkdir()
            paths = {}
            for where, at in (('faces', None), ('desert', DESERT)):
                geom, half_at_h = face_geometry(p, sc, shift, height, at)
                (hdir / f'{where}.geometry.json').write_text(json.dumps(geom, indent=1))
                prof = dict(profile)
                prof['input'] = dict(profile['input'], hdf5_path=str(PRODUCT),
                                     geometry_path=str(hdir / f'{where}.geometry.json'))
                paths[where] = hdir / f'{where}.profile.json'
                paths[where].write_text(json.dumps(prof, indent=1))
            res, regions_here = {}, []
            for case in CASES:
                def wrapped(config, config_path, case=case):
                    out = list(original(config, config_path))
                    if case.startswith('twin'):
                        out[7] = motionless_twin(out[7], np.random.default_rng(int(case[4:])))
                    elif case == 'planted':
                        out[7], regions = plant(out[7], out[4], out[5][0], out[6][0], g)
                        regions_here.extend(regions)
                    return tuple(out)
                gated.load_source_and_preflight = wrapped
                prof_path = paths['desert' if case == 'desert' else 'faces']
                gated.run_full(json.loads(prof_path.read_text()), prof_path, hdir / case)
                audit = np.load(hdir / case / 'biondi_v1_8_audit.npz')
                res[case] = analyse(audit, audit['kz_rad_per_m'])
                print(f"  {height:g} m, {case}: feature positions {res[case]['feature_positions']}", flush=True)
            real = np.load(hdir / 'real' / 'biondi_v1_8_audit.npz')
            pl = np.load(hdir / 'planted' / 'biondi_v1_8_audit.npz')
            for reg in regions_here:
                sd = SIDES.index(reg['side'])
                d = pl['Y'][sd][MIDDLE[0]:MIDDLE[1]] - real['Y'][sd][MIDDLE[0]:MIDDLE[1]]
                ideal = ideal_pair_response_px(g, reg['f_hz'], bank['mask_width_hz'], bank['bshift_hz'])
                plants.append({**reg, 'height_m': height, 'changed_fraction': float(np.mean(d != 0)),
                               'rms_change_px': float(np.sqrt(np.mean(d ** 2))),
                               'rms_change_range_azimuth_px': [float(np.sqrt(np.mean(d[..., c] ** 2))) for c in (0, 1)],
                               'full_tracker_rms_px': ideal / np.sqrt(2)})
            sets.append({'height_m': height, 'face_half_width_m': half_at_h, 'cases': res})
        gated.load_source_and_preflight = original
        zgrid = np.load(gdir / f'h{HEIGHTS[0]:g}' / 'real' / 'biondi_v1_8_audit.npz')['z_m']

        total = lambda case, key: sum(st['cases'][case][key] for st in sets)
        stage = lambda case, key: sum(v[key] for st in sets for v in st['cases'][case]['sides'].values())
        twins = [f'twin{seed}' for seed in TWIN_SEEDS]
        counts = {'real': total('real', 'feature_positions'), 'twins': [total(t, 'feature_positions') for t in twins],
                  'planted': total('planted', 'feature_positions'), 'desert': total('desert', 'feature_positions'),
                  'positions_per_image': len(HEIGHTS) * 4 * 101}
        runs = {'real': total('real', 'feature_runs'), 'twins': [total(t, 'feature_runs') for t in twins],
                'planted': total('planted', 'feature_runs'), 'desert': total('desert', 'feature_runs')}
        first_gate = {c: stage(c, 'window_passes') for c in CASES}
        track10 = {c: stage(c, 'track10_families') for c in CASES}

        # what planting changed: features in the planted image that the real image does not have at the same place
        new_runs, kept = [], 0
        for st in sets:
            for sd, v in st['cases']['planted']['sides'].items():
                have = {q for r in st['cases']['real']['sides'][sd]['runs'] for q in r}
                for r in runs_of(v['features']):
                    pos = [f['position'] for f in r]
                    if have & set(pos):
                        kept += 1
                        continue
                    modes = sorted({f['mode'] for f in r})
                    new_runs.append({'height_m': st['height_m'], 'side': sd, 'positions': pos, 'modes': modes,
                                     'planted_hz': dict(PLANTS).get(sd)})

        # depth against the winning mode: m x the depth one cycle spans in a window, or its mirror in the repeat
        feats = [dict(f, side=sd, case=c, height_m=st['height_m']) for st in sets for c, r in st['cases'].items()
                 for sd, v in r['sides'].items() for f in v['features']]
        zc = sets[0]['cases']['real']['depth_per_cycle_m']
        zrep = sets[0]['cases']['real']['repeat_m']
        near = np.array([min(abs(f['best_depth_m'] - f['mode_depth_m']), abs(f['best_depth_m'] - f['mirror_depth_m']))
                         for f in feats])
        within = int(np.sum(near <= zc / 2))
        chance = float(np.mean([np.mean((np.abs(zgrid - f['mode_depth_m']) <= zc / 2)
                                        | (np.abs(zgrid - f['mirror_depth_m']) <= zc / 2)) for f in feats])) if feats else None

        # from the geometry: the gate's modes as frequencies, what a pair's image keeps of them, the shape a vibration traces
        span_s = bank['mask_width_hz'] / abs(g.Ka_signed)
        pair_s = bank['k_leap_hz'] / abs(g.Ka_signed)
        window_s = gates['window_pairs'] * pair_s
        mode_hz = [m / window_s for m in range(gates['mode_range'][0], gates['mode_range'][1] + 1)]
        az_m_per_m_s = 2 * g.V / (g.lam * abs(g.Ka_signed))        # azimuth shift per unit line-of-sight velocity
        gate = {'window_s': window_s, 'pair_s': pair_s, 'mode_hz': mode_hz,
                'look_gain': [float(abs(np.sinc(f * span_s))) for f in mode_hz],
                'trace_width_ratio': [float(g.dx / (2 * np.pi * f * az_m_per_m_s * g.dr)) for f in mode_hz],
                'axis_ratio_min': gates['ellipse_axis_ratio_min'],
                'minor_min_px': gates['ellipse_minor_semiaxis_min_native_pixel'],
                'azimuth_m_per_m_s': az_m_per_m_s, 'pixel_m': [g.dx, g.dr]}

        zero = float(np.mean([v['exact_zero'] for st in sets for v in st['cases']['real']['sides'].values()]))
        step = float(np.mean([v['within_one_step'] for st in sets for v in st['cases']['real']['sides'].values()]))
        shared = 1 - bank['bshift_hz'] / bank['mask_width_hz']
        by_f = {f: [x for x in plants if x['f_hz'] == f] for _, f in PLANTS}
        rms = {f: float(np.mean([x['rms_change_px'] for x in v])) for f, v in by_f.items()}
        full = {f: float(np.mean([x['full_tracker_rms_px'] for x in v])) for f, v in by_f.items()}
        n_new = len(new_runs)
        long_over_wide = 1 / max(gate['trace_width_ratio'])
        unit = 10 ** (int(np.floor(np.log10(long_over_wide))) - 1)
        thin = np.floor(long_over_wide / unit) * unit              # rounded down: "at least"
        new_text = '; '.join(f"at mode {', '.join(map(str, r['modes']))} "
                             f"({', '.join(f'{m / window_s:.1f}' for m in r['modes'])} Hz) on the "
                             f"{r['side'].lower()} face" + (f", planted at {r['planted_hz']} Hz" if r['planted_hz'] else '')
                             for r in new_runs)
        finding = (
            f"A gated version of the method, run unchanged on {len(HEIGHTS)} loops round Khafre's faces in the 2022 image "
            f"({counts['positions_per_image']:,} positions), passes {counts['real']} positions "
            f"({plural(runs['real'], 'feature')}) on the real image, "
            f"{' and '.join(str(n) for n in counts['twins'])} "
            f"({' and '.join(str(n) for n in runs['twins'])}) on two motionless twins of it, with nothing moving and nothing "
            f"inside, and {counts['desert']} ({runs['desert']}) on the same loops laid on open desert in the same image. "
            f"With the faces made to vibrate at 2 mm/s along the line of sight, some 40,000 times Giza's microseisms, it "
            f"passes {counts['planted']} ({runs['planted']}): "
            + (f"the real image's feature again" if kept == runs['real'] == 1
               else f"{kept} of the real image's {plural(runs['real'], 'feature')} again")
            + (f", and {plural(n_new, 'new one')}, {new_text}" if n_new else ', and none new') + ". "
            f"Each image of its pairs spans {span_s:.1f} s of the pass and the two share {shared * 100:.2f}% of their "
            f"spectrum: {zero * 100:.0f}% of its shifts are exactly zero and {step * 100:.0f}% within one 0.01 px step, and "
            f"the planted vibrations move them by {rms[PLANTS[0][1]]:.3f} px at {PLANTS[0][1]} Hz and "
            f"{rms[PLANTS[1][1]]:.3f} px at {PLANTS[1][1]} Hz, where a tracker following them fully would see "
            f"{full[PLANTS[0][1]]:.3f} and {full[PLANTS[1][1]]:.3f} px. Its modes stand for {mode_hz[0]:.2f} to "
            f"{mode_hz[-1]:.1f} Hz, where a {span_s:.1f} s image keeps at most {max(gate['look_gain']) * 100:.1f}% of a "
            f"motion; a vibrating point's shifts trace an ellipse at least "
            f"{thin:,.0f} times as long as it is wide, and the gate keeps only "
            f"those at most {1 / gate['axis_ratio_min']:g} times."
            + (f" Of the {len(feats)} passing positions, {within} have their best depth within half a cycle of their mode "
               f"times {zc:.2f} m or its mirror ({zrep:.1f} m less that), against {chance * 100:.0f}% by chance."
               if feats else ""))
        run.save({'registration': {'shift_px': list(shift), 'correlation': corr}, 'heights_m': HEIGHTS, 'sets': sets,
                  'features': feats, 'planted': plants, 'counts': counts, 'runs': runs, 'first_gate_windows': first_gate,
                  'track10_families': track10, 'planted_new': new_runs, 'planted_kept': kept,
                  'depth_per_cycle_m': zc, 'repeat_m': zrep, 'depth_within_half_cycle': within,
                  'depth_chance_fraction': chance, 'gate': gate, 'pair_span_s': span_s, 'pair_shared_fraction': shared,
                  'real_exact_zero': zero, 'real_within_one_step': step,
                  'planted_rms_px': rms, 'full_tracker_rms_px': full, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
