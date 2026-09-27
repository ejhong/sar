"""P1-04 · Ambient noise: can natural vibration reveal the chamber?

    uv run python experiments/p1_04_ambient.py [--reuse]

The one-chamber bench, the same 9 x 9 geophone grid as the surface surveys,
and no hammer: 48 noise sources at random places on the surface (people
walking, traffic), each an independent stream of band-limited vertical
force. Every geophone records the vertical ground velocity.

1. Retrieval. Cross-correlating two geophones' noise gives the wave a
   source at one would send to the other (seismic interferometry). With the
   noise sources spread over the ground rather than on a boundary enclosing
   the geophones, it is minus the lag derivative of the correlation that
   carries the wave (Snieder 2004, Phys. Rev. E 69, 046610). The ensemble-average
   correlation is computed exactly from each noise source's response; the
   retrieved wave is checked against a direct simulation with a real source
   at the first geophone.
2. Passive echo imaging. Correlations turn 16 geophones into virtual
   sources. As in P1-03, the rock without the chamber is known exactly (the
   most favourable case): its predicted virtual gathers are subtracted, and
   what remains is migrated back.
3. Record length. A finite record adds fluctuations to the correlations
   whose variance falls as 1/T. Migration is linear, so the image of one
   fluctuation realisation at a reference length gives the image noise at
   any length: the record needed for the chamber to stand out follows.

Each source's response is recorded long enough (T_RESP) for its waves to
have crossed the whole array, so the correlations are complete over the lag
window used (T_LAG, less the wavelet delay).
"""
import copy
import sys
import time
from pathlib import Path

import numpy as np

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run
from katabasis.seismic.arrays import snap_to_ground, surface_grid
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, ricker
from katabasis.seismic.imaging import laplace_filter, rtm_shot

SITE = 'bench-void'
H = 0.5
F0 = 100.0                 # centre of the noise band
T_RESP = 0.052             # each noise source's response: its surface wave crosses the array in time
T_LAG = 0.037              # correlation lags kept; the first T0 hold the virtual source's wavelet
T0 = 1.5 / F0              # the virtual source's wavelet (noise autocorrelation) is centred here
N_SOURCES = 48
SOURCE_RADIUS = 24.0       # inside the absorbing boundary (which begins 26 m out)
N_VIRTUAL = 16
SEG_REF = 100              # correlation segments in the reference fluctuation realisation
CACHE = RESULTS / 'cache' / 'p1_04'
BOX = ((-14.0, 14.0), (-14.0, 14.0), (-30.0, -1.0))


def models():
    site = load_site(SITE)
    g = Grid.covering((-34.0, 34.0), (-34.0, 34.0), (-42.0, 2.0), H)
    full = Medium.from_model(voxelise(site, g))
    raw = copy.deepcopy(site.raw)
    raw['features'] = []
    bg = Medium.from_model(voxelise(parse_site(raw, site.directory), g))
    return site, g, full, bg


def responses(med, g, sources, receivers, name, reuse):
    """Vertical velocity at every geophone from a unit band-limited vertical force at each noise source."""
    path = CACHE / f'{name}.npz'
    if reuse and path.exists():
        d = np.load(path)
        return d['z'], float(d['dt']), d['wavelet']
    sim = Simulation(med, pml_width=16, f0=F0)
    nt = int(T_RESP / sim.dt)
    w = ricker(np.arange(nt) * sim.dt, F0) * 1e9
    out = np.zeros((len(sources), len(receivers), nt), np.float32)
    nte = int(0.004 / sim.dt)
    taper = np.ones(nt)
    taper[nt - nte:] = 0.5 * (1 + np.cos(np.linspace(0, np.pi, nte)))   # no hard cut where a late wave is still passing
    t0 = time.time()
    for n, s in enumerate(sources):
        sim.reset()
        out[n] = sim.run([Source(tuple(s), w, 'force', (0, 0, -1))], Receivers(receivers), nt).traces[:, 2] * taper
        if n % 6 == 0:
            print(f'  {name}: source {n + 1}/{len(sources)} ({time.time() - t0:.0f} s)', flush=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, z=out, dt=sim.dt, wavelet=w)
    return out, sim.dt, w


def correlations(R: np.ndarray, virtual: list[int], nt: int, rng=None, segments: int = 0) -> np.ndarray:
    """Virtual gathers (n_virtual, n_rec, nt lags): the symmetrised causal correlation of each virtual
    source's record with every geophone's. With `segments`, one finite-record realisation (random source
    phases, averaged over that many segments, each 2 x the response length); without, the ensemble average."""
    nf = 2 * R.shape[-1]
    F = np.fft.rfft(R, n=nf, axis=-1)                     # (ns, nr, nfreq)
    if segments:
        acc = np.zeros((len(virtual), R.shape[1], F.shape[-1]), np.complex128)
        for _ in range(segments):
            N = (rng.standard_normal((F.shape[0], F.shape[-1])) + 1j * rng.standard_normal((F.shape[0], F.shape[-1]))) / np.sqrt(2)
            U = np.einsum('srf,sf->rf', F, N)             # every geophone's record in this segment
            acc += U[virtual][:, None, :] * np.conj(U[None, :, :])
        C = acc / segments
    else:
        C = np.einsum('saf,srf->arf', F[:, virtual], np.conj(F))
    c = np.fft.irfft(C, n=nf, axis=-1)
    causal = c[..., :nt]
    acausal = np.concatenate([c[..., :1], c[..., :nf - nt:-1]], axis=-1)[..., :nt]
    return 0.5 * (np.conj(causal) + acausal) if np.iscomplexobj(c) else 0.5 * (causal + acausal)


def main():
    reuse = '--reuse' in sys.argv
    site, g, full, bg = models()
    rng = np.random.default_rng(4)
    r = SOURCE_RADIUS * np.sqrt(rng.random(N_SOURCES))
    th = 2 * np.pi * rng.random(N_SOURCES)
    sources = snap_to_ground(g, full.solid, np.stack([r * np.cos(th), r * np.sin(th)], 1))
    receivers = surface_grid(g, full.solid, 16.0, 4.0)
    virtual = list(np.linspace(0, len(receivers) - 1, N_VIRTUAL).round().astype(int))
    params = {'site': SITE, 'sim_spacing_m': H, 'noise_band_centre_hz': F0, 'response_length_s': T_RESP,
              'correlation_lags_s': T_LAG, 'virtual_wavelet_delay_s': T0,
              'noise_sources': N_SOURCES, 'noise_source_area': f'random over the ground within {SOURCE_RADIUS:.0f} m of the centre',
              'geophones': len(receivers), 'virtual_sources': N_VIRTUAL, 'background': 'exact (the site without its chamber)',
              'fluctuation_reference_segments': SEG_REF}
    with Run('p1_04_ambient', 'Ambient noise over one chamber', params) as run:
        Rf, dt, w = responses(full, g, sources, receivers, 'with', reuse)
        Rb, _, _ = responses(bg, g, sources, receivers, 'without', reuse)
        nt = int(T_LAG / dt)
        seg_s = 2 * Rf.shape[-1] * dt

        # 1 · retrieval: the virtual gather from the first virtual source against a real source there
        C_full = correlations(Rf, virtual, nt)
        C_bg = correlations(Rb, virtual, nt)
        a = virtual[0]
        # the virtual source's wavelet: the noise wavelet's autocorrelation, zero lag delayed by t0
        full_ac = np.correlate(w, w, 'full')              # zero lag at index len(w) - 1
        k0 = int(round(T0 / dt))
        src_w = full_ac[len(w) - 1 - k0: len(w) - 1 - k0 + nt].astype(np.float64)
        src_w *= 1e9 / np.abs(src_w).max()
        # the virtual gather: minus the lag derivative of the correlation, zero lag moved to t0 like the wavelet
        gather = lambda c: np.concatenate([np.zeros(c.shape[:-1] + (k0,)), -np.gradient(c, dt, axis=-1)[..., : c.shape[-1] - k0]], axis=-1)
        sim = Simulation(bg, pml_width=16, f0=F0)
        assert abs(sim.dt - dt) < 1e-12 * dt
        direct = sim.run([Source(tuple(receivers[a]), src_w, 'force', (0, 0, -1))], Receivers(receivers), nt).traces[:, 2]
        # compare at a geophone about 16 m away
        d_ab = np.linalg.norm(receivers - receivers[a], axis=1)
        b = int(np.argmin(np.abs(d_ab - 16.0)))
        cor = gather(C_bg[0, b])
        dirb = direct[b].astype(np.float64)
        ncc = float(np.dot(cor, dirb) / np.sqrt(np.dot(cor, cor) * np.dot(dirb, dirb)))
        step = max(1, nt // 220)
        retrieval = {'geophones': [int(a), int(b)], 'distance_m': float(d_ab[b]), 'ncc': ncc,
                     't_ms': (np.arange(nt)[::step] * dt * 1e3).round(3),
                     'correlation': (cor[::step] / np.abs(cor).max()).round(4),
                     'direct': (dirb[::step] / np.abs(dirb).max()).round(4)}
        print(f'  retrieval: geophones {a}-{b}, {d_ab[b]:.1f} m, correlation vs direct ncc {ncc:.3f}')

        # 2 · passive echo imaging (ensemble average) and 3 · one finite-record fluctuation
        scat_mean = gather(C_full - C_bg)
        C_fin = correlations(Rf, virtual, nt, np.random.default_rng(9), SEG_REF)
        fluct = gather(C_fin - C_full)
        data_scale = 1e9 / np.abs(scat_mean).max()          # one scale for both images, so they compare linearly
        echo_rms = float(np.sqrt(np.mean(scat_mean ** 2)))
        fl_rms = float(np.sqrt(np.mean(fluct ** 2)))
        print(f'  echo {echo_rms:.3g} vs fluctuation at {SEG_REF * seg_s:.1f} s {fl_rms:.3g}')

        def image(data):
            box = (slice(int(round((BOX[0][0] - g.origin[0]) / H)), int(round((BOX[0][1] - g.origin[0]) / H)) + 1),
                   slice(int(round((BOX[1][0] - g.origin[1]) / H)), int(round((BOX[1][1] - g.origin[1]) / H)) + 1),
                   slice(int(round((g.origin[2] - BOX[2][1]) / H)), int(round((g.origin[2] - BOX[2][0]) / H)) + 1))
            img, ill = None, None
            for k, va in enumerate(virtual):
                d3 = np.zeros((len(receivers), 3, nt))
                d3[:, 2] = data[k]
                i1, l1 = rtm_shot(sim, Source(tuple(receivers[va]), src_w, 'force', (0, 0, -1)),
                                  receivers, d3, nt, box, data_scale=data_scale)
                img = i1 if img is None else img + i1
                ill = l1 if ill is None else ill + l1
            out = laplace_filter(img / (ill + 1e-3 * ill.max()), H)
            return out, box

        t_img = time.time()
        img_mean, box = image(scat_mean)
        print(f'  mean image ({time.time() - t_img:.0f} s)')
        img_fl, _ = image(fluct)
        print(f'  fluctuation image ({time.time() - t_img:.0f} s)')

        xs, ys, zs = g.x[box[0]], g.y[box[1]], g.z[box[2]]
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
        c = np.array(site.feature('chamber')['shape']['centre'])
        near = (np.abs(X - c[0]) <= 5) & (np.abs(Y - c[1]) <= 5) & (np.abs(Z - c[2]) <= 5)
        far = ~near & (Z < -3)
        env = np.abs(img_mean)
        signal = float(env[near].mean())
        systematic = float(env[far].mean())
        fl_ref = float(np.abs(img_fl)[far].mean())
        T_ref = SEG_REF * seg_s
        lengths = np.array([1, 10, 60, 600, 3600, 36000, 360000], float)
        contrast = signal / (systematic + fl_ref * np.sqrt(T_ref / lengths))
        target = 3.0
        # the record at which the contrast reaches the target; none if the systematic image alone forbids it
        need = float(T_ref * (fl_ref / (signal / target - systematic)) ** 2) if signal / target > systematic else None
        peak = np.unravel_index(np.argmax(np.where(Z < -3, env, 0)), env.shape)
        score = {'image_contrast_ensemble': signal / systematic,
                 'brightest_offset_from_chamber_m': float(np.linalg.norm(np.array([X[peak], Y[peak], Z[peak]]) - c)),
                 'echo_over_fluctuation_at_ref': echo_rms / fl_rms, 'reference_record_s': T_ref,
                 'record_for_contrast_3_s': need,
                 'contrast_vs_record': [{'record_s': float(L), 'contrast': float(cn)} for L, cn in zip(lengths, contrast)]}
        print('  score', {k: v for k, v in score.items() if k != 'contrast_vs_record'})
        j = int(np.argmin(np.abs(ys - c[1])))
        q = np.round(np.clip(env / np.percentile(env, 99.9), 0, 1) * 255).astype(np.uint8)
        np.savez_compressed(run.dir / 'volumes.npz', ambient=q)
        e2 = env / np.percentile(env, 99.9)
        run.save({
            'stations': {'noise_sources': sources.round(2).tolist(), 'receivers': receivers.round(2).tolist(),
                         'virtual': [int(v) for v in virtual]},
            'retrieval': retrieval, 'score': score,
            'slices': {'x': xs[::2].round(2), 'z': zs[::2].round(2), 'xz_image': e2[::2, j, ::2].T.round(3)},
            'volume': {'shape': list(q.shape), 'origin': [float(xs[0]), float(ys[0]), float(zs[0])], 'spacing': H,
                       'quantity': 'migrated passive echo amplitude', 'units': 'normalised'},
            'finding': (f"Noise correlations reproduce the wave between two geophones (correlation {ncc:.2f} with a real source). "
                        f"With the rock known exactly and an unlimited record, passive echo imaging lights the chamber at "
                        f"{signal / systematic:.1f}x the image elsewhere; "
                        + (f"reaching 3x takes about {need:,.0f} s of noise." if need else "3x is not reached at any record length, "
                           "because surface noise sends too little energy down to the chamber and back.")),
        })


if __name__ == '__main__':
    main()
