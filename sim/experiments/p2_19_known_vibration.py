"""P2-19 · A known vibration through the unchanged pipeline: synthetic products, a corner reflector, a tracker beside it.

    uv run --with scikit-image python experiments/p2_19_known_vibration.py

A measurement control for the gated reconstruction, asked for by an independent review: before a method's silence
about a chamber means anything, it has to be shown responding to motion that is known to be there. Four synthetic
products are written in the ICEYE SLC layout (sarsim.slcfile), 8,192 rows by 384 columns on the 2022 geometry, and read
by the reconstruction's own loader, unchanged. Each holds static speckle (band-limited complex Gaussian, the product's
spectral support) and one corner reflector 50 dB above it at full resolution, simulated pulse by pulse (sarsim.echo)
and added to the speckle, which is exact since imaging is linear: still; vibrating along the line of sight at 5 mm/s at
1.07 Hz (the gate's mode 1) and at 3.21 Hz (mode 3); and at 20 mm/s at 0.26 Hz. Five lines of 101 positions run
across range through and beside the reflector, 4 rows apart.

Beside the pipeline, on the same files, a standard reading: the reflector's magnitude tracked through 50 ms looks,
which P2-18 showed recovers a 5 mm/s vibration. Scored: what the tracker recovers; the pipeline's shifts at the
positions whose patches hold the reflector, at the planted frequency, against the still product; and what its gates
pass there and elsewhere.
"""
import json
from pathlib import Path

import numpy as np

import p2_13_gates as m
import p2_15_pair_response as r15
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup, PointScene, simulate
from sarsim.looks import look_masks, looks, velocity_series
from sarsim.slcfile import write_synthetic_slc

RID = 'p2_19_known_vibration'
ROWS, COLS = 8192, 384
SCR_DB = 50.0
WIN = (512, 32)                     # half-size of the window the reflector is simulated in (rows, columns)
CASES = [('still', 0.0, 0.0), ('mode1_5mm', 5e-3, 1.0691), ('mode3_5mm', 5e-3, 3.2074), ('slow_20mm', 2e-2, 0.26)]
TRACK_ROWS = [-8, -4, 0, 4, 8]      # rows of the five lines relative to the reflector
SUPPORTS = [1, 2, 4, 5]


def reflector_image(st, amp, v, f):
    """The reflector alone, pulse by pulse, in a window about the image centre, pasted into a full-size image."""
    sc = PointScene.still([0.0], [0.0], [amp])
    if v:
        sc = sc.moving(np.array([True]), ar=v / (2 * np.pi * f), f=f)
    win = simulate(st, sc, 2 * WIN[0], 2 * WIN[1])
    full = np.zeros((ROWS, COLS), np.complex64)
    r0, c0 = ROWS // 2 - WIN[0], COLS // 2 - WIN[1]
    full[r0:r0 + 2 * WIN[0], c0:c0 + 2 * WIN[1]] = win
    return full


def geometry_json(lat, lon, height):
    tracks = {}
    for j, dr in enumerate(TRACK_ROWS):
        row = ROWS // 2 + dr
        tracks[f'L{j}'] = {'native_endpoints_col_row': [[COLS // 2 - 100, row], [COLS // 2 + 100, row]],
                           'geographic_endpoint_lat_lon': [[lat, lon - 1e-4], [lat, lon + 1e-4]]}
    return {'title': 'Five lines across range through a synthetic corner reflector', 'projection': {'height_m': height},
            'coordinate_order': 'Geographic endpoints are [latitude, longitude]. Native pixels are [column, row].',
            'tracks': tracks}


def main():
    g = DwellGeometry.from_record(m.ACQ)
    st = EchoSetup.from_geometry(g)
    lat, lon, height = g.source['scene_centre']['latitude'], g.source['scene_centre']['longitude'], 75.0
    profile = json.loads(m.PROFILE.read_text())
    profile['gates'] = dict(profile['gates'], pixel_support_variants=SUPPORTS)
    params = {'acquisition': m.ACQ, 'shape': [ROWS, COLS], 'scr_db': SCR_DB, 'cases': CASES, 'track_rows': TRACK_ROWS,
              'supports': SUPPORTS, 'gated_module_sha256': m.sha256(m.MODULE), 'gated_profile_sha256': m.sha256(m.PROFILE)}
    with Run(RID, 'A known vibration through the unchanged pipeline', params) as run:
        gdir = Path(__file__).resolve().parents[1] / 'data' / RID
        gdir.mkdir(parents=True, exist_ok=True)
        (gdir / 'geometry.json').write_text(json.dumps(geometry_json(lat, lon, height), indent=1))
        clutter = r15.band_noise(g, (ROWS, COLS), np.random.default_rng(19))
        amp = float(np.sqrt(10 ** (SCR_DB / 10)))
        gated = m.load_gated()
        gated.render_tomograms = lambda *a, **k: None
        W_look = 0.05
        centres = np.arange(-2.0, 2.0 + 1e-9, 0.025)
        sub_rows = 4096                                            # the tracker reads a 4,096-row window about the reflector
        masks_look = look_masks(g, sub_rows, centres, W_look)
        rs = slice(ROWS // 2 - sub_rows // 2, ROWS // 2 + sub_rows // 2)
        cs = slice(COLS // 2 - 8, COLS // 2 + 8)
        rr = slice(sub_rows // 2 - 400, sub_rows // 2 + 400)
        cc = slice(8 - 2, 8 + 3)
        results = {}
        audits = {}
        for name, v, f in CASES:
            img = clutter + reflector_image(st, amp, v, f)
            path = write_synthetic_slc(gdir / f'SYNTHETIC_{name}.h5', img, st, g.Ka_signed, lat, lon, height, g.heading_deg,
                                       g.theta_deg, scene={'speckle_seed': 19, 'reflector_scr_db': SCR_DB,
                                                           'vibration_los_m_s': v, 'vibration_hz': f})
            prof = dict(profile)
            prof['input'] = dict(profile['input'], hdf5_path=str(path), geometry_path=str(gdir / 'geometry.json'))
            (gdir / f'profile_{name}.json').write_text(json.dumps(prof, indent=1))
            gated.run_full(prof, gdir / f'profile_{name}.json', gdir / name)
            audits[name] = np.load(gdir / name / 'biondi_v1_8_audit.npz')
            # the standard reading on the same file
            tracker = None
            if v:
                vel = velocity_series(looks(img[rs, cs], masks_look), g, rr, cc, envelope=True)
                A = np.stack([np.cos(2 * np.pi * f * centres), np.sin(2 * np.pi * f * centres), np.ones_like(centres)], 1)
                coef, *_ = np.linalg.lstsq(A, vel, rcond=None)
                se = float(np.std(vel - A @ coef) * np.sqrt(2.0 / len(centres)))
                tracker = {'recovered_m_s': float(np.hypot(coef[0], coef[1])),
                           'truth_look_average_m_s': float(v * abs(np.sinc(f * W_look))), 'se_m_s': se}
            results[name] = {'v_m_s': v, 'f_hz': f, 'tracker': tracker}
            print(f"  {name}: tracker {tracker}", flush=True)

        # the pipeline at the reflector: positions whose 32 x 32 cores hold it, on the line through it
        still = audits['still']
        n_pos = still['Y'].shape[1]
        cols_of = np.round(np.linspace(COLS // 2 - 100, COLS // 2 + 100, n_pos)).astype(int)
        near = np.flatnonzero(np.abs(cols_of - COLS // 2) <= 8)
        far = np.flatnonzero(np.abs(cols_of - COLS // 2) > 40)
        profile_bank = profile['processing']['filter_bank']
        k = np.arange(profile_bank['pair_count'])
        centres_hz = -0.5 * profile_bank['support_hz'] + 0.5 * profile_bank['mask_width_hz'] + k * profile_bank['k_leap_hz'] \
            + 0.5 * profile_bank['bshift_hz']
        t_pair = centres_hz / g.Ka_signed
        centre_track = TRACK_ROWS.index(0)
        for name, v, f in CASES:
            a = audits[name]
            rec = results[name]
            if v:
                dY = a['Y'][centre_track][near] - still['Y'][centre_track][near]      # [positions, pairs, 2]
                amps = [r15.fit(dY[i], t_pair, f)[0] for i in range(len(near))]
                rec['pipeline_shift_amp_px'] = {'range': float(np.max([x[0] for x in amps])),
                                                'azimuth': float(np.max([x[1] for x in amps]))}
            wp = a[f'l{centre_track}_window_pass']
            modes = a[f'l{centre_track}_window_mode']
            rec['first_gate_windows_at_reflector'] = int(wp[near].sum())
            rec['first_gate_windows_per_position_far'] = float(wp[far].sum() / max(len(far), 1))
            won = modes[near][wp[near]]
            rec['modes_passing_at_reflector'] = {int(x): int(c) for x, c in zip(*np.unique(won, return_counts=True))}
            rec['passing_positions'] = {P: {f'L{j}': [int(q) for q in np.flatnonzero(np.isfinite(a[f'focus_p{P}'][j]).any(axis=1))]
                                            for j in range(len(TRACK_ROWS))} for P in SUPPORTS}
            rec['passes_at_reflector'] = {P: int(sum(np.isin(np.flatnonzero(np.isfinite(a[f'focus_p{P}'][j]).any(axis=1)), near).sum()
                                                  for j in range(len(TRACK_ROWS)))) for P in SUPPORTS}
            print(f"  {name}: first-gate windows at the reflector {rec['first_gate_windows_at_reflector']}, "
                  f"passes there {rec['passes_at_reflector']}", flush=True)
        moving = [(n, results[n]) for n, v, f in CASES if v]
        finding = (
            "Four synthetic products in the ICEYE layout, read by the reconstruction's own loader, hold static speckle and a "
            f"corner reflector {SCR_DB:.0f} dB above it, simulated pulse by pulse. A standard reading of the same files, the "
            "reflector's magnitude tracked through 50 ms looks, recovers "
            + '; '.join(f"{r['tracker']['recovered_m_s'] * 1e3:.2f} ± {r['tracker']['se_m_s'] * 1e3:.2f} mm/s of a look-averaged "
                        f"{r['tracker']['truth_look_average_m_s'] * 1e3:.2f} mm/s at {r['f_hz']:.2f} Hz" for _, r in moving)
            + ". Through the unchanged pipeline, the shifts at the positions holding the reflector move at the planted frequency by "
            + '; '.join(f"{r['pipeline_shift_amp_px']['azimuth']:.4f} px ({r['f_hz']:.2f} Hz)" for _, r in moving)
            + f", against a rounding step of 0.01 px; its gates pass {results['still']['passes_at_reflector'][1]} positions at "
            f"support one at the still reflector and "
            + ', '.join(f"{r['passes_at_reflector'][1]} ({r['f_hz']:.2f} Hz)" for _, r in moving)
            + " at the vibrating one, and at the profile's own support of four "
            + ', '.join(str(results[n]['passes_at_reflector'][4]) for n, _, _ in CASES) + " in the four products.")
        run.save({'results': results, 'near_positions': near.tolist(), 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
