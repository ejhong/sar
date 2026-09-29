"""P2-14 · The reconstruction's own lines round Khafre's base, with the same controls.

    uv run --with scikit-image python experiments/p2_14_own_lines.py

P2-13 ran the gated pipeline along loops the lab drew on Khafre's faces. The profile it ships with names its own frozen
geometry: Khafre's four base corners (Petrie's, placed on a map by Dash) projected through the product's RPC at
H = 70.5 m, the base's height above sea level, given to the RPC directly. A second version of that file adds the EGM96
undulation (h = H + N, 85.9 m on the ellipsoid, which the RPC expects), moving the lines about 50 px in range. Both files
are run here as they were sent, their hashes recorded, next to the profile where it looks for them; the corrected
file's projection.height_m is a text note, which the pipeline's steering reads as a number, so it is set to the
corners' ellipsoidal height. A third geometry lays the corrected square on open ground south of Khafre, P2-12's
in-block control, through the same RPC over the site's surface.

Cases, the same code: the real crop; two motionless twins (seeds 101 and 102); the real crop with 2 mm/s line-of-sight
vibrations planted on the north side at 0.26 Hz and the south side at 3.66 Hz (positions 30 to 70, as in P2-13); and,
on the desert square, the real crop. The profile is unchanged except that pixel support is reported for one to five
contiguous positions: its own four and five, and the looser settings besides.

Scored: positions passing every gate at each support, the first gate's windows, each passing position's best depth
against its winning mode times the depth one cycle spans (or its mirror), what planting adds and removes, and, from
the pipeline's own steering record, the Doppler rate it uses and how far apart the four sides' zero-Doppler moments
fall.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

import p2_13_gates as m
from katabasis.compose import load_site
from katabasis.export.sites import scene as site_scene
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.dwell import DwellProduct
from sarsim.geolocation import project_rpc
from sarsim.looks import motionless_twin
from sarsim.ortho import M_PER_DEG_LAT, frame_to_lla, project, surface

RID = 'p2_14_own_lines'
DATA = m.GATED / 'data'
SEA = DATA / 'ICEYE_X13_Khafre_Petrie_Dash_FourSides_H70p5m_FrozenGeometry.json'
WGS = DATA / 'ICEYE_X13_Khafre_Petrie_Dash_FourSides_EGM96_WGS84_CorrectedGeometry.json'
SUPPORTS = [1, 2, 3, 4, 5]
DESERT = m.DESERT
SIDES = m.SIDES
PLAN = {'sea_level': ['real', 'twin101', 'twin102', 'planted'], 'wgs84': ['real', 'twin101', 'twin102', 'planted'],
        'desert': ['real']}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lla_to_frame(origin, lat, lon):
    lat0, lon0 = origin['latitude'], origin['longitude']
    return ((np.asarray(lon) - lon0) * M_PER_DEG_LAT * np.cos(np.deg2rad(lat0)),
            (np.asarray(lat) - lat0) * M_PER_DEG_LAT)


def desert_geometry(p, sc, wgs):
    """The corrected square moved from Khafre's base centre to open ground, heights from the site's surface."""
    o = sc['frame']['origin']
    geom = {'title': f"The reconstruction's Khafre base square laid on open ground at {DESERT} m", 'sides': {}}
    hs = []
    for name, side in wgs['sides'].items():
        ends_ll, ends_px = [], []
        for label in (side['P0_corner'], side['P_last_corner']):
            x, y = lla_to_frame(o, *wgs['corner_lat_lon'][label])
            x, y = np.array([x + DESERT[0]]), np.array([y + DESERT[1]])
            z = surface(sc, x, y)
            hs.append(float(z[0]) + m.GEOID_M)
            r, c = project(p, sc, x, y, z, m.GEOID_M)
            lat, lon = frame_to_lla(o, x, y)
            ends_ll.append([float(lat[0]), float(lon[0])])
            ends_px.append([float(c[0]), float(r[0])])
        geom['sides'][name] = {'geographic_endpoint_lat_lon': ends_ll, 'native_endpoints_col_row': ends_px}
    geom['projection'] = {'height_m': float(np.mean(hs))}
    return geom


def passing(audit, P):
    return {SIDES[s]: [int(q) for q in np.flatnonzero(np.isfinite(audit[f'focus_p{P}'][s]).any(axis=1))]
            for s in range(len(SIDES))}


def features(audit, P):
    """Each passing position's winning mode and best depth, against mode x the depth per cycle or its mirror."""
    z = audit['z_m']
    out = []
    for s, side in enumerate(SIDES):
        dk = float(np.median(np.abs(np.diff(audit['kz_rad_per_m'][s]))))
        zc, zrep = 2 * np.pi / (50 * dk), 2 * np.pi / dk
        F, ws = audit[f'focus_p{P}'][s], audit[f'winning_start_p{P}'][s]
        modes = audit[f'{side.lower()}_window_mode']
        for q in np.flatnonzero(np.isfinite(F).any(axis=1)):
            j = int(np.nanargmax(F[q]))
            st = int(ws[q, j])
            if st < 0:
                continue
            mo = int(modes[q, st])
            near = min(abs(z[j] - mo * zc), abs(z[j] - (zrep - mo * zc)))
            chance = float(np.mean((np.abs(z - mo * zc) <= zc / 2) | (np.abs(z - (zrep - mo * zc)) <= zc / 2)))
            out.append({'side': side, 'position': int(q), 'mode': mo, 'best_depth_m': float(z[j]),
                        'at_mode_depth': bool(near <= zc / 2), 'chance': chance})
    return out


def main():
    g = DwellGeometry.from_record(m.ACQ)
    sc = site_scene(load_site('giza'))
    p = DwellProduct(m.PRODUCT)
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=SUPPORTS)
    sea, wgs = json.loads(SEA.read_text()), json.loads(WGS.read_text())
    h_wgs = float(np.mean(list(wgs['corner_ellipsoid_height_WGS84_m'].values())))
    params = {'acquisition': m.ACQ, 'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE),
              'geometry_sea_level_sha256': sha256(SEA), 'geometry_wgs84_sha256': sha256(WGS),
              'wgs84_projection_height_m': h_wgs, 'supports': SUPPORTS, 'planted': m.PLANTS, 'planted_v_m_s': m.V_AMP,
              'planted_positions': list(m.MIDDLE), 'twin_seeds': m.TWIN_SEEDS, 'desert_centre_m': DESERT}
    with Run(RID, "The reconstruction's own lines round Khafre's base, with the same controls", params) as run:
        # the lab's RPC projection of the corrected corners against the file's own
        rpc_check = {}
        for label, (lat, lon) in wgs['corner_lat_lon'].items():
            rc = project_rpc(p._rpc, np.array([lat]), np.array([lon]),
                             np.array([wgs['corner_ellipsoid_height_WGS84_m'][label]]))
            side = next(s for s in wgs['sides'].values() if label in (s['P0_corner'], s['P_last_corner']))
            k = 0 if side['P0_corner'] == label else 1
            rpc_check[label] = [float(rc[0, 0] - side['raw_rpc_endpoint_row_col'][k][0]),
                                float(rc[0, 1] - side['raw_rpc_endpoint_row_col'][k][1])]
        geoms = {'sea_level': sea,
                 'wgs84': dict(wgs, projection=dict(wgs['projection'], height_m=h_wgs)),
                 'desert': desert_geometry(p, sc, wgs)}

        gdir = Path(__file__).resolve().parents[1] / 'data' / RID
        if gdir.exists():
            import shutil
            shutil.rmtree(gdir)
        gdir.mkdir(parents=True)
        gated = m.load_gated()
        original = gated.load_source_and_preflight
        res, steering = {}, {}
        for where, cases in PLAN.items():
            d = gdir / where
            d.mkdir()
            (d / 'geometry.json').write_text(json.dumps(geoms[where], indent=1))
            prof = dict(profile)
            prof['input'] = dict(profile['input'], hdf5_path=str(m.PRODUCT), geometry_path=str(d / 'geometry.json'))
            (d / 'profile.json').write_text(json.dumps(prof, indent=1))
            res[where] = {}
            for case in cases:
                def wrapped(config, config_path, case=case):
                    out = list(original(config, config_path))
                    if case.startswith('twin'):
                        out[7] = motionless_twin(out[7], np.random.default_rng(int(case[4:])))
                    elif case == 'planted':
                        out[7], _ = m.plant(out[7], out[4], out[5][0], out[6][0], g)
                    return tuple(out)
                gated.load_source_and_preflight = wrapped
                gated.run_full(prof, d / 'profile.json', d / case)
                audit = np.load(d / case / 'biondi_v1_8_audit.npz')
                a = m.analyse(audit, audit['kz_rad_per_m'])
                res[where][case] = {
                    'passing': {P: passing(audit, P) for P in SUPPORTS},
                    'positions': {P: sum(len(v) for v in passing(audit, P).values()) for P in SUPPORTS},
                    'features': {P: features(audit, P) for P in SUPPORTS},
                    'first_gate_windows': int(sum(v['window_passes'] for v in a['sides'].values())),
                    'track10_families': int(sum(v['track10_families'] for v in a['sides'].values())),
                    'exact_zero': float(np.mean([v['exact_zero'] for v in a['sides'].values()]))}
                if case == 'real':
                    st = json.loads((d / case / 'manifest.json').read_text())['steering_audit']
                    times = [np.datetime64(x['reference_time_utc'].replace('+00:00', '')) for x in st]
                    steering[where] = {'doppler_rate_hz_s': [x['doppler_rate_hz_per_s'] for x in st],
                                       'zero_doppler_spread_s': float((max(times) - min(times)) / np.timedelta64(1, 'us') / 1e6)}
                print(f"  {where}, {case}: " + ', '.join(f"P{P} {res[where][case]['positions'][P]}" for P in SUPPORTS), flush=True)
        gated.load_source_and_preflight = original

        # what planting adds and removes, against the real image on the same lines
        stretch = set(range(*m.MIDDLE))
        plant_hz = dict(m.PLANTS)
        rate = float(np.mean(np.abs(steering['wgs84']['doppler_rate_hz_s'])))
        window_s = 50 * profile['processing']['filter_bank']['k_leap_hz'] / rate     # seconds of the pass per window
        changes = {}
        for where in ('sea_level', 'wgs84'):
            changes[where] = {}
            for P in SUPPORTS:
                r, q = res[where]['real']['passing'][P], res[where]['planted']['passing'][P]
                modes = {(f['side'], f['position']): f['mode'] for f in res[where]['planted']['features'][P]}
                new = [{'side': sd, 'position': x, 'mode': modes.get((sd, x)), 'in_stretch': x in stretch,
                        'planted_hz': plant_hz.get(sd)} for sd in SIDES for x in sorted(set(q[sd]) - set(r[sd]))]
                lost = sum(len(set(r[sd]) - set(q[sd])) for sd in SIDES)
                changes[where][P] = {'new': new, 'lost': lost}
        # a planted frequency comes back only as the mode nearest to it, in its own stretch
        near_mode = {sd: f * window_s for sd, f in plant_hz.items()}
        returned = [dict(c, where=w, support=P) for w in changes for P in SUPPORTS for c in changes[w][P]['new']
                    if c['in_stretch'] and c['planted_hz'] and c['mode'] is not None
                    and abs(c['mode'] - near_mode[c['side']]) < 0.5]

        pooled = {P: [f for w in res for c in res[w] for f in res[w][c]['features'][P]] for P in SUPPORTS}
        depth = {P: {'n': len(v), 'at_mode_depth': int(sum(f['at_mode_depth'] for f in v)),
                     'chance': float(np.mean([f['chance'] for f in v])) if v else None} for P, v in pooled.items()}
        pos = {w: {c: res[w][c]['positions'] for c in res[w]} for w in res}
        span = profile['processing']['filter_bank']['mask_width_hz'] / rate

        W, S = pos['wgs84'], pos['sea_level']
        own = [P for P in SUPPORTS if P >= 4]
        nothing = all(W[c][P] == 0 for c in W for P in own) and all(pos['desert']['real'][P] == 0 for P in own)
        fmt = lambda d, P: f"{d['real'][P]}, twins {d['twin101'][P]} and {d['twin102'][P]}, planted {d['planted'][P]}"
        n_new_n = sum(1 for w in changes for P in SUPPORTS for c in changes[w][P]['new'] if c['in_stretch'] and c['side'] == 'North')
        finding = (
            "On the reconstruction's own lines round Khafre's base (Petrie's corners at 70.5 m above sea level, as sent "
            "and corrected to the ellipsoid), its unchanged pipeline "
            + (f"passes nothing at its own support of {own[0]} or {own[-1]} positions on the corrected lines, in any "
               f"image: real, two motionless twins, planted vibrations, open desert"
               if nothing else f"passes at support {own[0]}: " + fmt(W, own[0]))
            + f"; on the uncorrected lines, real {S['real'][4]}, twins {S['twin101'][4]} and {S['twin102'][4]}, planted "
            f"{S['planted'][4]}. With support loosened to one position, the corrected lines pass {fmt(W, 1)}, and open desert "
            f"in the same image {pos['desert']['real'][1]}; at two, {fmt(W, 2)}, desert {pos['desert']['real'][2]}. The "
            f"planted vibrations do not return at their frequencies: the 0.26 Hz is below every mode, and the {n_new_n} "
            f"positions that start passing in its stretch do so at modes {sorted({c['mode'] for w in changes for P in SUPPORTS for c in changes[w][P]['new'] if c['in_stretch'] and c['side'] == 'North'})}; "
            f"{len(returned)} position{'s' if len(returned) != 1 else ''} start{'' if len(returned) != 1 else 's'} passing "
            f"in the 3.66 Hz stretch at its nearest mode, and over all supports "
            f"{sum(len(changes[w][P]['new']) for w in changes for P in SUPPORTS)} positions start passing and "
            f"{sum(changes[w][P]['lost'] for w in changes for P in SUPPORTS)} stop. Of the {depth[1]['n']} positions "
            f"passing at support one, {depth[1]['at_mode_depth']} have their best depth within half a cycle of their mode "
            f"times the depth per cycle or its mirror, against {depth[1]['chance'] * 100:.0f}% by chance. Its own steering "
            f"uses a Doppler rate of {rate:,.0f} Hz/s at Khafre, so each {profile['processing']['filter_bank']['mask_width_hz']:,.0f} "
            f"Hz image holds {span:.2f} s of every pixel's pass, and puts the four sides' zero-Doppler moments within "
            f"{steering['wgs84']['zero_doppler_spread_s']:.3f} s of each other.")
        run.save({'rpc_check_row_col_px': rpc_check, 'steering': steering, 'doppler_rate_hz_s': rate,
                  'image_span_s': span, 'positions': pos,
                  'first_gate_windows': {w: {c: res[w][c]['first_gate_windows'] for c in res[w]} for w in res},
                  'exact_zero': {w: {c: res[w][c]['exact_zero'] for c in res[w]} for w in res},
                  'planted_changes': changes, 'planted_returned': returned, 'depth_vs_mode': depth,
                  'passing': {w: {c: res[w][c]['passing'] for c in res[w]} for w in res}, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
