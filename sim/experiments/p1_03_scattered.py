"""P1-03 · The chamber's echo: imaging the scattered wave from the surface.

    uv run python experiments/p1_03_scattered.py [--reuse]

The same hammer and geophone grids as P1-02's surface survey, over the same
chamber. Records are made twice, with and without the chamber; their
difference is the chamber's echo. Noise is added at the level of P1-02
(30 dB below each record's peak) and the echo is migrated back through the
rock by reverse-time migration. The rock without the chamber is known
exactly, so this is the most favourable case for any waveform method.

Also written: the ground's vertical motion for the central shot, with the
chamber and as its echo alone, for the viewer.
"""
import copy
import sys
import time

import numpy as np

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, ricker, surface_index
from katabasis.seismic.imaging import laplace_filter, rtm_shot
from katabasis.seismic.picking import add_noise

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from p1_02_traveltime import CACHE as CACHE02, F0, H_SIM, SITE, SNR_DB, T_REC, build, surveys  # noqa: E402

CACHE = RESULTS / 'cache' / 'p1_03'
BOX = ((-14.0, 14.0), (-14.0, 14.0), (-30.0, -1.0))


def background_medium(site, g):
    raw = copy.deepcopy(site.raw)
    raw['features'] = []
    bg = parse_site(raw, site.directory)
    return Medium.from_model(voxelise(bg, g))


def records(sv, med, name, reuse):
    path = CACHE / f'{name}.npz'
    if reuse and path.exists():
        d = np.load(path)
        return d['traces'], float(d['dt']), d['wavelet']
    sim = Simulation(med, pml_width=16, f0=F0)
    nt = int(T_REC / sim.dt)
    w = ricker(np.arange(nt) * sim.dt, F0) * 1e9
    out = np.zeros((len(sv.sources), len(sv.receivers), 3, nt), np.float32)
    t0 = time.time()
    for n, s in enumerate(sv.sources):
        sim.reset()
        out[n] = sim.run([Source(tuple(s), w, 'force', (0, 0, -1))], Receivers(sv.receivers), nt).traces
        if n % 5 == 0:
            print(f'  {name}: shot {n + 1}/{len(sv.sources)} ({time.time() - t0:.0f} s)', flush=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, traces=out, dt=sim.dt, wavelet=w)
    return out, sim.dt, w


def index_box(g: Grid):
    x0, y0, ztop = g.origin
    h = g.spacing
    ix = slice(int(round((BOX[0][0] - x0) / h)), int(round((BOX[0][1] - x0) / h)) + 1)
    iy = slice(int(round((BOX[1][0] - y0) / h)), int(round((BOX[1][1] - y0) / h)) + 1)
    iz = slice(int(round((ztop - BOX[2][1]) / h)), int(round((ztop - BOX[2][0]) / h)) + 1)
    return ix, iy, iz


def surface_movie(med_full, med_bg, sv, g):
    """Vertical ground velocity for the central shot, every few steps: with the chamber, and the echo alone."""
    centre = int(np.argmin(np.linalg.norm(sv.sources[:, :2], axis=1)))
    s = sv.sources[centre]
    frames = {}
    for name, med in (('with', med_full), ('without', med_bg)):
        sim = Simulation(med, pml_width=16, f0=F0)
        nt = int(T_REC / sim.dt)
        w = ricker(np.arange(nt) * sim.dt, F0) * 1e9
        res = sim.run([Source(tuple(s), w, 'force', (0, 0, -1))], Receivers(sv.receivers[:1]), nt,
                      snapshot_every=6, surface_index=surface_index(med))
        frames[name] = res.snapshots
        t = res.snap_t
    return frames['with'], frames['with'] - frames['without'], t, s


def main():
    reuse = '--reuse' in sys.argv
    site, g, model, med = build()
    sv = surveys(g, med)[1]                             # the surface survey
    med_bg = background_medium(site, g)
    params = {'site': SITE, 'sim_spacing_m': H_SIM, 'f0_hz': F0, 'record_s': T_REC, 'snr_db': SNR_DB,
              'image_box': {'x': BOX[0], 'y': BOX[1], 'z': BOX[2]}, 'background': 'exact (the site without its chamber)'}
    with Run('p1_03_scattered', "The chamber's echo, imaged from the surface", params) as run:
        d02 = CACHE02 / 'surface.npz'
        if d02.exists():
            z = np.load(d02)
            full, dt, w = z['traces'], float(z['dt']), z['wavelet']
        else:
            full, dt, w = records(sv, med, 'with', reuse)
        bg, _, _ = records(sv, med_bg, 'without', reuse)
        ns, nr = full.shape[:2]
        scat = full - bg
        rng = np.random.default_rng(3)
        noisy_full = add_noise(full.reshape(ns * nr, 3, -1), SNR_DB, rng).reshape(full.shape)
        noise = noisy_full - full
        scat_noisy = scat + noise
        e_scat = np.sqrt((scat.astype(np.float64) ** 2).sum((2, 3)))
        e_full = np.sqrt((full.astype(np.float64) ** 2).sum((2, 3)))
        e_noise = np.sqrt((noise.astype(np.float64) ** 2).sum((2, 3)))
        ratio_db = 20 * np.log10(np.median(e_scat / e_full))
        echo_to_noise_db = 20 * np.log10(np.median(e_scat / e_noise))
        print(f'  echo / record {ratio_db:.1f} dB, echo / noise {echo_to_noise_db:.1f} dB')

        box = index_box(g)
        sim = Simulation(med_bg, pml_width=16, f0=F0)
        nt = full.shape[-1]
        image = None
        illum = None
        t0 = time.time()
        for n, s in enumerate(sv.sources):
            img, ill = rtm_shot(sim, Source(tuple(s), w, 'force', (0, 0, -1)), sv.receivers, scat_noisy[n], nt, box)
            image = img if image is None else image + img
            illum = ill if illum is None else illum + ill
            if n % 5 == 0:
                print(f'  rtm: shot {n + 1}/{ns} ({time.time() - t0:.0f} s)', flush=True)
        img = image / (illum + 1e-3 * illum.max())
        img = laplace_filter(img, g.spacing)
        env = np.abs(img)
        env /= np.percentile(env, 99.9)
        xs, ys, zs = g.x[box[0]], g.y[box[1]], g.z[box[2]]
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
        ch = site.feature('chamber')['shape']
        c, half = np.array(ch['centre']), np.array(ch['size']) / 2
        near = (np.abs(X - c[0]) <= half[0] + 2) & (np.abs(Y - c[1]) <= half[1] + 2) & (np.abs(Z - c[2]) <= half[2] + 2)
        shallow = Z > -3.0                              # the acquisition footprint right under the stations
        far = ~near & ~shallow
        peak = np.unravel_index(np.argmax(np.where(shallow, 0, env)), env.shape)
        score = {'echo_over_record_db': float(ratio_db), 'echo_over_noise_db': float(echo_to_noise_db),
                 'image_near_chamber_mean': float(env[near].mean()), 'image_elsewhere_mean': float(env[far].mean()),
                 'image_elsewhere_p99': float(np.percentile(env[far], 99)),
                 'contrast': float(env[near].mean() / env[far].mean()),
                 'brightest_below_3m': [float(X[peak]), float(Y[peak]), float(Z[peak])],
                 'brightest_offset_from_chamber_m': float(np.linalg.norm(np.array([X[peak], Y[peak], Z[peak]]) - c))}
        print('  score', score)
        j = int(np.argmin(np.abs(ys - c[1])))
        k = int(np.argmin(np.abs(zs - c[2])))
        q = np.round(np.clip(env / 1.0, 0, 1) * 255).astype(np.uint8)
        movie, echo, mt, shot = surface_movie(med, med_bg, sv, g)
        np.savez_compressed(run.dir / 'volumes.npz', rtm=q, movie=movie, echo=echo, movie_t=mt)
        # decimate slices to the 1 m grid of the page figures
        run.save({
            'survey': sv.summary() | {'stations': {'sources': sv.sources.round(2).tolist(), 'receivers': sv.receivers.round(2).tolist()}},
            'score': score,
            'slices': {'x': xs[::2].round(2), 'y': ys[::2].round(2), 'z': zs[::2].round(2),
                       'xz_image': env[::2, j, ::2].T.round(3), 'xy_image': env[::2, ::2, k].T.round(3),
                       'plane_y': float(ys[j]), 'plane_z': float(zs[k])},
            'volume': {'shape': list(q.shape), 'origin': [float(xs[0]), float(ys[0]), float(zs[0])], 'spacing': g.spacing,
                       'quantity': 'migrated echo amplitude', 'units': 'normalised'},
            'movie': {'shot': shot.round(2).tolist(), 'frames': int(movie.shape[0]), 'grid_origin': list(g.origin),
                      'spacing': g.spacing, 'shape': list(movie.shape[1:]), 'dt_ms': float((mt[1] - mt[0]) * 1e3)},
            'finding': (f"The chamber's echo is {abs(ratio_db):.0f} dB below the records and {echo_to_noise_db:+.0f} dB against "
                        f"their noise; migrated back through the known rock it lights up the chamber at {score['contrast']:.1f}× the "
                        f"image elsewhere, brightest {score['brightest_offset_from_chamber_m']:.1f} m from its centre."),
        })


if __name__ == '__main__':
    main()
