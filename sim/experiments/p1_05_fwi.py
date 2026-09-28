"""P1-05 · Full-waveform inversion from a practitioner's start.

    uv run python experiments/p1_05_fwi.py [--reuse]

The records of P1-02's surface survey (25 hammer points, 81 three-component
geophones, 120 Hz, noise 30 dB below each record's peak, the same noisy
records P1-03 migrated) are fitted whole, by adjoint-state gradients of P
and S speed. Unlike P1-03, the rock is not given: the start is what a
practitioner would have, uniform rock at the P speed the surface travel
times gave (P1-02) and an S speed from the textbook ratio Vp/Vs = √3,
which is 3% fast here. Density is held at its assumed value.

One stage, the full band: at 120 Hz the start is only about 16 degrees of
phase out at the farthest geophones, so no low-frequency stage is needed to
avoid cycle skips, and below 60 Hz this source carries too little energy
above the noise (a 60 Hz stage, tried, fitted noise and the near field).
Traces within 4 m of their source leave the misfit (the near field, which
no practitioner fits), and the iterations stop when the misfit reaches the
noise level, estimated from each record's quiet first 2 ms. The data are
computed with the same solver and grid that made them (the "inverse
crime"); the noise and the wrong start soften it, and it is stated wherever
the result is shown.
"""
import json
import sys
import time

import numpy as np
from scipy.ndimage import map_coordinates

from katabasis.compose import Grid, load_site, voxelise
from katabasis.runs import RESULTS, Run
from katabasis.seismic.arrays import snap_to_ground
from katabasis.seismic.fwi import Problem, gaussian_sigma, lbfgs, lowpass, noise_misfit
from katabasis.seismic.picking import add_noise

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from p1_02_traveltime import CACHE as CACHE02, F0, SITE, SNR_DB, build, surveys  # noqa: E402

# the whole survey footprint, from the first layer of rock down: the start is wrong everywhere the waves go,
# so everywhere they go may change (the absorbing layers begin 20 m out on the fine grid)
BOX = ((-19.0, 19.0), (-19.0, 19.0), (-30.0, -0.2))
STAGES = [
    {'name': 'full', 'h': 0.5, 'lowpass_hz': None, 'iterations': 10, 'extent': 28.0, 'pml': 16},
]
QUIET_S = 0.002             # the start of every record, before any arrival: its noise level
CACHE = RESULTS / 'cache' / 'p1_05'


def index_box(g: Grid):
    x0, y0, ztop = g.origin
    h = g.spacing
    return (slice(int(round((BOX[0][0] - x0) / h)), int(round((BOX[0][1] - x0) / h)) + 1),
            slice(int(round((BOX[1][0] - y0) / h)), int(round((BOX[1][1] - y0) / h)) + 1),
            slice(int(round((ztop - BOX[2][1]) / h)), int(round((ztop - BOX[2][0]) / h)) + 1))


def start_model(site, g, vp0, vs0):
    """Uniform rock at the practitioner's speeds, with the site's air and assumed density."""
    raw = json.loads(json.dumps(site.raw))
    raw['features'] = []
    from katabasis.compose.site import parse_site
    m = voxelise(parse_site(raw, site.directory), g)
    solid = ~m.air
    return (np.where(solid, vp0, 0.0), np.where(solid, vs0, 0.0), np.where(solid, m.rho, 0.0), m.air)


def score(g: Grid, box, vp, vs, site) -> dict:
    xs, ys, zs = g.x[box[0]], g.y[box[1]], g.z[box[2]]
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    ch = site.feature('chamber')['shape']
    c, half = np.array(ch['centre']), np.array(ch['size']) / 2
    inside = (np.abs(X - c[0]) <= half[0]) & (np.abs(Y - c[1]) <= half[1]) & (np.abs(Z - c[2]) <= half[2])
    ring = ~inside & (np.abs(X - c[0]) <= 12) & (np.abs(Y - c[1]) <= 12) & (np.abs(Z - c[2]) <= 12)
    out = {}
    for name, v, true in (('vp', vp, 3300.0), ('vs', vs, 1830.0)):
        bg = float(np.median(v))
        rel = v / bg - 1
        k = np.unravel_index(np.argmin(rel), rel.shape)
        out[name] = {'background_m_s': bg, 'background_true_m_s': true,
                     'chamber_mean_rel': float(rel[inside].mean()), 'chamber_min_rel': float(rel[inside].min()),
                     'surroundings_std_rel': float(rel[ring].std()),
                     'contrast_to_noise': float(-rel[inside].mean() / max(rel[ring].std(), 1e-6)),
                     'slowest_cell': [float(X[k]), float(Y[k]), float(Z[k])],
                     'slowest_offset_from_chamber_m': float(np.linalg.norm(np.array([X[k], Y[k], Z[k]]) - c))}
    return out


def main():
    smoke = '--smoke' in sys.argv                     # two shots, one iteration per stage: checks the plumbing
    site, g_f, _, med = build()
    sv = surveys(g_f, med)[1]
    z = np.load(CACHE02 / 'surface.npz')
    full, dt_f, w_f = z['traces'], float(z['dt']), z['wavelet']
    ns, nr = full.shape[:2]
    obs_f = add_noise(full.reshape(ns * nr, 3, -1), SNR_DB, np.random.default_rng(3)).reshape(full.shape).astype(np.float64)
    tt = json.loads((RESULTS / 'p1_02_traveltime' / 'summary.json').read_text())
    vp0 = float(tt['surveys']['surface']['score']['background_m_s'])
    vs0 = vp0 / np.sqrt(3.0)
    params = {'site': SITE, 'survey': sv.description, 'f0_hz': F0, 'snr_db': SNR_DB,
              'start': {'vp_m_s': vp0, 'vs_m_s': vs0, 'from': 'P1-02 surface travel times; Vs = Vp / sqrt(3)'},
              'stages': STAGES, 'box': {'x': BOX[0], 'y': BOX[1], 'z': BOX[2]},
              'optimiser': 'L-BFGS (memory 5), illumination preconditioner, 1-cell smoothing of the gradient, muted within 1.5-3 m of stations',
              'misfit': 'L2 on the velocity records, traces within 4 m of their source excluded; stop at the noise level',
              'inverse_crime': 'the full-band stage uses the grid and solver that made the records'}
    with Run('p1_05_fwi' + ('_smoke' if smoke else ''), "Full-waveform inversion from a practitioner's start", params) as run:
        prev = None
        stages_out = []
        for st in STAGES:
            h = st['h']
            e = st['extent']
            g = Grid.covering((-e, e), (-e, e), (-42.0, 2.0 + h), h)       # at least three air cells at any spacing
            vp, vs, rho, air = start_model(site, g, vp0, vs0)
            dec = int(round(h / g_f.spacing))
            sigma = gaussian_sigma(st['lowpass_hz']) if st['lowpass_hz'] else 0.0
            # records low-passed before decimation (no aliased noise); synthetics are filtered inside the misfit
            obs = (lowpass(obs_f, sigma, dt_f) if sigma else obs_f)[..., ::dec]
            w = w_f.astype(np.float64)[::dec]
            box = index_box(g)
            assert not air[box].any() and air[box[0], box[1], box[2].start - 1].all(), 'the box must start at the first rock layer'
            src = snap_to_ground(g, ~air, sv.sources[:, :2])
            rec = snap_to_ground(g, ~air, sv.receivers[:, :2])
            if smoke:
                src, obs = src[:2], obs[:2]
            prob = Problem(g, rho, air, vp, vs, box, dt_f * dec, src, rec, obs, np.ascontiguousarray(w), F0, pml=st['pml'],
                           vp_bounds=(0.35, 1.12), vs_bounds=(0.35, 1.12),
                           sigma_s=sigma, obs_filtered=True)
            start = None
            if prev is not None:                       # carry the coarse result onto this grid
                pg, pdvp, pdvs = prev
                xs, ys, zs = g.x[box[0]], g.y[box[1]], g.z[box[2]]
                X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
                bxs, bys, bzs = pg.x[index_box(pg)[0]], pg.y[index_box(pg)[1]], pg.z[index_box(pg)[2]]
                idx = np.stack([(X - bxs[0]) / pg.spacing, (Y - bys[0]) / pg.spacing, (bzs[0] - Z) / pg.spacing])
                start = (map_coordinates(pdvp, idx, order=1, mode='nearest'), map_coordinates(pdvs, idx, order=1, mode='nearest'))
            t0 = time.time()
            print(f"  stage {st['name']}: grid {g.shape} ({np.prod(g.shape) / 1e6:.2f} M cells), box {prob.vp0[box].shape}", flush=True)
            target = noise_misfit(prob.obs, prob.weights, prob.dt, QUIET_S) if not st['lowpass_hz'] else None
            res = lbfgs(prob, 1 if smoke else st['iterations'], start=start, label=st['name'], target=target)
            vpi, vsi = prob.speeds(res['dvp'], res['dvs'])
            sc = score(g, box, vpi[box], vsi[box], site)
            print(f"  stage {st['name']}: misfit {res['history'][-1]['misfit']:.3f}, Vp chamber {100 * sc['vp']['chamber_mean_rel']:+.1f}% "
                  f"(CNR {sc['vp']['contrast_to_noise']:.1f}), Vs chamber {100 * sc['vs']['chamber_mean_rel']:+.1f}% "
                  f"(CNR {sc['vs']['contrast_to_noise']:.1f}) ({time.time() - t0:.0f} s)", flush=True)
            stages_out.append({'name': st['name'], 'spacing_m': h, 'lowpass_hz': st['lowpass_hz'], 'history': res['history'], 'score': sc,
                               'noise_misfit': None if target is None else target / res['misfit_start'],
                               'traces_kept': int(prob.weights.sum()), 'traces': int(prob.weights.size)})
            CACHE.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(CACHE / f"{st['name']}.npz", dvp=res['dvp'], dvs=res['dvs'])
            prev = (g, res['dvp'], res['dvs'])

        # the final model: slices and volumes
        xs, ys, zs = g.x[box[0]], g.y[box[1]], g.z[box[2]]
        c = site.feature('chamber')['shape']['centre']
        j = int(np.argmin(np.abs(ys - c[1])))
        k = int(np.argmin(np.abs(zs - c[2])))
        rel_p = vpi[box] / np.median(vpi[box]) - 1
        rel_s = vsi[box] / np.median(vsi[box]) - 1
        lo = -0.12
        vols = {}
        for name, rel in (('vp', rel_p), ('vs', rel_s)):
            vols[name] = np.round(np.clip(rel / lo, 0, 1) * 255).astype(np.uint8)
        np.savez_compressed(run.dir / 'volumes.npz', **vols)
        final = stages_out[-1]['score']
        run.save({
            'stations': {'sources': src.round(2).tolist(), 'receivers': rec.round(2).tolist()},
            'stages': stages_out, 'score': final,
            'slices': {'x': xs.round(2), 'y': ys.round(2), 'z': zs.round(2),
                       'xz_vp': rel_p[:, j, :].T.round(4), 'xy_vp': rel_p[:, :, k].T.round(4),
                       'xz_vs': rel_s[:, j, :].T.round(4), 'xy_vs': rel_s[:, :, k].T.round(4)},
            'volume': {'shape': list(vols['vs'].shape), 'origin': [float(xs[0]), float(ys[0]), float(zs[0])], 'spacing': g.spacing,
                       'range': [lo, 0.0], 'units': 'fraction of background'},
            'finding': (f"From uniform rock at the travel-time P speed and a textbook S speed, {len(stages_out[-1]['history']) - 1} "
                        f"iterations fitting the whole records bring the misfit to {stages_out[-1]['history'][-1]['misfit']:.2f} "
                        f"(the noise alone: {stages_out[-1]['noise_misfit']:.2f}) and leave the chamber at "
                        f"{100 * final['vs']['chamber_mean_rel']:+.0f}% in S speed "
                        f"({final['vs']['contrast_to_noise']:.1f}x the scatter around it) and {100 * final['vp']['chamber_mean_rel']:+.0f}% "
                        f"in P speed; the background S speed moves from {vs0:.0f} to {final['vs']['background_m_s']:.0f} m/s "
                        f"(true 1830)."),
        })


if __name__ == '__main__':
    main()
