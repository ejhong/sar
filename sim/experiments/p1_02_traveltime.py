"""P1-02 · First-arrival travel-time tomography over one chamber.

    uv run python experiments/p1_02_traveltime.py [--reuse]

The ground is the one-chamber bench (a 6 m air-filled cube, top 12 m down,
in bare Mokattam limestone). Every record comes from the 3-D elastic solver;
noise is added; first breaks are picked; the picks are inverted on a coarser
grid by a code that never sees the model. Two surveys:

    crosshole  four boreholes around the chamber, sources and receivers down them
    surface    a hammer grid and a geophone grid on the ground

`--reuse` skips the wave simulations if their records are cached.
"""
import sys
import time

import numpy as np

from katabasis.compose import Grid, load_site, voxelise
from katabasis.runs import RESULTS, Run
from katabasis.seismic.arrays import Survey, borehole, surface_grid
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, ricker
from katabasis.seismic.picking import add_noise, pick_all
from katabasis.seismic.traveltime import TomoGrid, invert

SITE = 'bench-void'
H_SIM = 0.5
F0 = 120.0                  # surface hammer
F0_CROSS = 250.0            # borehole sparker: first arrivals are P, resolved at ten points per wavelength
T_REC = 0.034
T_REC_CROSS = 0.022
SNR_DB = 30.0
H_INV = 1.0
BOX = (-18.0, 18.0)        # inversion box, x and y: holds every station
ZBOT = -32.0
CACHE = RESULTS / 'cache' / 'p1_02'


def build():
    site = load_site(SITE)
    g = Grid.covering((-28.0, 28.0), (-28.0, 28.0), (-42.0, 2.0), H_SIM)
    model = voxelise(site, g)
    return site, g, model, Medium.from_model(model)


def surveys(g, med):
    holes = [(-12.0, 0.0), (12.0, 0.0), (0.0, -12.0), (0.0, 12.0)]
    src = np.concatenate([borehole(g, med.solid, xy, 5.0, 29.0, 3.0) for xy in holes])
    rec = np.concatenate([borehole(g, med.solid, xy, 5.0, 30.0, 1.0) for xy in holes])
    cross = Survey('crosshole', 'Four boreholes 12 m from the chamber; sparker every 3 m, hydrophones every 1 m, 5-30 m deep',
                   src, rec, 'explosion', {'boreholes': [list(h) for h in holes], 'f0_hz': F0_CROSS})
    surf = Survey('surface', 'A 5 x 5 hammer grid (8 m) over a 9 x 9 geophone grid (4 m), vertical impacts',
                  surface_grid(g, med.solid, 16.0, 8.0), surface_grid(g, med.solid, 16.0, 4.0), 'force', {'f0_hz': F0})
    return [cross, surf]


def simulate(sv: Survey, med: Medium, reuse: bool):
    path = CACHE / f'{sv.name}.npz'
    if reuse and path.exists():
        d = np.load(path)
        return d['traces'], float(d['dt']), d['wavelet']
    f0 = sv.meta['f0_hz']
    sim = Simulation(med, pml_width=16, f0=f0)
    nt = int((T_REC_CROSS if sv.name == 'crosshole' else T_REC) / sim.dt)
    t = np.arange(nt) * sim.dt
    w = ricker(t, f0) * 1e9
    out = np.zeros((len(sv.sources), len(sv.receivers), 3, nt), np.float32)
    t0 = time.time()
    for n, s in enumerate(sv.sources):
        sim.reset()
        res = sim.run([Source(tuple(s), w, sv.source_kind, (0, 0, -1))], Receivers(sv.receivers), nt)
        out[n] = res.traces
        if n % 5 == 0:
            print(f'  {sv.name}: shot {n + 1}/{len(sv.sources)} ({time.time() - t0:.0f} s)', flush=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, traces=out, dt=sim.dt, wavelet=w)
    return out, sim.dt, w


def pairs_for(sv: Survey) -> np.ndarray:
    pr = []
    for a, s in enumerate(sv.sources):
        for b, r in enumerate(sv.receivers):
            same_hole = sv.name == 'crosshole' and np.allclose(s[:2], r[:2])
            if same_hole or np.linalg.norm(s - r) < 3.0:
                continue
            pr.append((a, b))
    return np.array(pr)


def score(inv, model_grid: Grid, site) -> dict:
    tg = inv.grid.grid
    v = 1 / inv.slowness
    X, Y, Z = np.meshgrid(tg.x, tg.y, tg.z, indexing='ij')
    ch = site.feature('chamber')['shape']
    c, half = np.array(ch['centre']), np.array(ch['size']) / 2
    inside = (np.abs(X - c[0]) <= half[0]) & (np.abs(Y - c[1]) <= half[1]) & (np.abs(Z - c[2]) <= half[2]) & inv.grid.active
    ring = (~inside) & inv.grid.active & (np.abs(X - c[0]) <= 12) & (np.abs(Y - c[1]) <= 12) & (np.abs(Z - c[2]) <= 12)
    bg = np.median(v[inv.grid.active])
    rel = v / bg - 1
    cov = inv.coverage
    return {'background_m_s': float(bg),
            'chamber_mean_rel': float(rel[inside].mean()), 'chamber_min_rel': float(rel[inside].min()),
            'surroundings_std_rel': float(rel[ring].std()),
            'contrast_to_noise': float(-rel[inside].mean() / max(rel[ring].std(), 1e-6)),
            'slowest_cell': [float(X[np.unravel_index(np.argmin(rel), rel.shape)]), float(Y[np.unravel_index(np.argmin(rel), rel.shape)]),
                             float(Z[np.unravel_index(np.argmin(rel), rel.shape)])],
            'chamber_ray_length_m': float(cov[inside].sum()) if cov is not None else None,
            'true_air_rel': float(343.0 / bg - 1)}


def slices(inv, site) -> dict:
    """Signed relative speed on the vertical plane through the chamber and the horizontal plane at its centre."""
    tg = inv.grid.grid
    v = 1 / inv.slowness
    bg = np.median(v[inv.grid.active])
    rel = v / bg - 1
    c = site.feature('chamber')['shape']['centre']
    j = int(np.argmin(np.abs(tg.y - c[1])))
    k = int(np.argmin(np.abs(tg.z - c[2])))
    cov = inv.coverage
    return {'x': tg.x.round(2), 'y': tg.y.round(2), 'z': tg.z.round(2),
            'xz_rel': rel[:, j, :].T.round(4), 'xy_rel': rel[:, :, k].T.round(4),
            'xz_coverage_m': cov[:, j, :].T.round(2), 'plane_y': float(tg.y[j]), 'plane_z': float(tg.z[k])}


def volume(inv, name: str) -> dict:
    """The tomogram as a relative-velocity grid for the viewer (slow anomalies bright)."""
    tg = inv.grid.grid
    v = 1 / inv.slowness
    bg = np.median(v[inv.grid.active])
    rel = np.where(inv.grid.active, v / bg - 1, 0)
    lo = -0.12
    q = np.clip(rel / lo, 0, 1)                 # 0 = background or faster, 1 = 12% slow or more
    data = np.round(q * 255).astype(np.uint8)
    return {'name': name, 'data': data, 'grid': tg, 'range': [lo, 0.0], 'quantity': 'relative P-wave speed',
            'units': 'fraction of background', 'coverage': inv.coverage}


def main():
    reuse = '--reuse' in sys.argv
    site, g, model, med = build()
    print(f'grid {g.shape} ({np.prod(g.shape) / 1e6:.2f} M cells), vs_min {med.vs_min:.0f} m/s, '
          f'{med.vs_min / (2.5 * F0) / H_SIM:.1f} points per S wavelength at 2.5 f0')
    params = {'site': SITE, 'sim_spacing_m': H_SIM, 'f0_hz': {'surface': F0, 'crosshole': F0_CROSS},
              'record_s': {'surface': T_REC, 'crosshole': T_REC_CROSS}, 'snr_db': SNR_DB,
              'inversion_spacing_m': H_INV, 'inversion_box': {'x': BOX, 'y': BOX, 'z': [ZBOT, 0.0]}}
    with Run('p1_02_traveltime', 'First-arrival travel-time tomography over one chamber', params) as run:
        rng = np.random.default_rng(2)
        out = {'surveys': {}}
        vols = []
        for sv in surveys(g, med):
            traces, dt, w = simulate(sv, med, reuse)
            ns, nr = traces.shape[:2]
            noisy = add_noise(traces.reshape(ns * nr, 3, -1), SNR_DB, rng)
            picks = pick_all(noisy, 0.5 * dt, dt, w, period=1.0 / sv.meta['f0_hz']).reshape(ns, nr)
            pr = pairs_for(sv)
            t_obs = picks[pr[:, 0], pr[:, 1]]
            inv_grid = Grid.covering(BOX, BOX, (ZBOT, 0.0), H_INV)
            tg = TomoGrid(inv_grid, np.ones(inv_grid.shape, bool))
            # starting model: the straight-ray average speed of the survey's own picks
            dist = np.linalg.norm(sv.sources[pr[:, 0]] - sv.receivers[pr[:, 1]], axis=1)
            v0 = float(np.nanmedian(dist / t_obs))
            t0 = time.time()
            inv = invert(tg, sv.sources, sv.receivers, pr, t_obs, 1 / v0, iterations=6, smooth=3.0, damp=0.05)
            sc = score(inv, g, site)
            print(f'  {sv.name}: v0 {v0:.0f} m/s, rms {inv.history[0]["rms_ms"]:.3f} -> {inv.residual_rms * 1e3:.3f} ms, '
                  f'chamber {100 * sc["chamber_mean_rel"]:+.1f}% (min {100 * sc["chamber_min_rel"]:+.1f}%), '
                  f'CNR {sc["contrast_to_noise"]:.1f} ({time.time() - t0:.0f} s)')
            straight = dist / 3300.0
            out['surveys'][sv.name] = {
                **sv.summary(), 'pairs': len(pr), 'picked': int(np.isfinite(t_obs).sum()),
                'starting_speed_m_s': v0, 'history': inv.history, 'score': sc,
                'pick_minus_straight_ms': {'median': float(np.nanmedian(t_obs - straight) * 1e3),
                                           'p95_abs': float(np.nanpercentile(np.abs(t_obs - straight), 95) * 1e3)},
                'stations': {'sources': sv.sources.round(2).tolist(), 'receivers': sv.receivers.round(2).tolist()}}
            out['surveys'][sv.name]['slices'] = slices(inv, site)
            vols.append(volume(inv, sv.name))
        np.savez_compressed(run.dir / 'volumes.npz', **{v['name']: v['data'] for v in vols},
                            **{v['name'] + '_coverage': v['coverage'] for v in vols})
        out['volumes'] = [{'name': v['name'], 'shape': list(v['data'].shape), 'origin': list(v['grid'].origin),
                           'spacing': v['grid'].spacing, 'range': v['range'], 'quantity': v['quantity'],
                           'units': v['units']} for v in vols]
        cs, ss = out['surveys']['crosshole']['score'], out['surveys']['surface']['score']
        out['finding'] = (f"Crosshole rays through the chamber recover it as a slow patch of {100 * cs['chamber_mean_rel']:+.0f}% "
                          f"(air would be {100 * cs['true_air_rel']:.0f}%): first arrivals bend around a void, so travel-time "
                          f"tomography sees its shadow, not its emptiness. From the surface the same method recovers "
                          f"{100 * ss['chamber_mean_rel']:+.1f}%: no first-arrival ray reaches 12 m in uniform rock.")
        run.save(out)


if __name__ == '__main__':
    main()
