"""P2-25 · What one image can hold about the chamber, whatever reads it.

    uv run python experiments/p2_25_information_bound.py

P2-03 to P2-24 test methods by running them. This bounds them all at once, with every method not yet written. A radar
image is made of the echo from what the wave reaches, the surface, and the surface's motion during the pass; every
method is a function of the image, so none can hold more about the chamber than the image does, and for ordinary
ground the image's statistics say how much that is (sarsim.information): any test's detection rate can exceed its
false-alarm rate by at most sqrt(F) / 2, F the Fisher information the image holds about the chamber's imprint.

1. The routes. The echo from the chamber: P2-11 puts the bench chamber's roof at least 347 dB down at X band; noted,
   not pursued. The surface: what any photograph shows. The motion: bounded here, generously to the claim; the method
   is handed the ambient motion itself (which only raises the information), receiver noise is left out (likewise), and
   the rock has no attenuation.
2. The bound, checked. The exact sum against a brute-force covariance on a tiny image; the continuum formula used on the
   real geometry against the exact sum; shared motion giving exactly zero; and the bound attained: the locally most
   powerful detector, on the model's speckle and on the lab's synthesizer, scatters by sqrt(F) without the imprint and
   shifts by a F with it, so the bound is what the image holds, not a loose ceiling over it.
3. The chamber, the one-chamber bench (6 m, 15 m down), on the real Giza dwell (P2-01). Below 8 Hz its imprint is
   P2-04's static one (the kernels recomputed if absent), averaged over the waves' directions, under Giza's measured
   microseisms, the noisiest stations on Earth, and the measured 1-3 and 3-8 Hz levels. Above 6 Hz it is P2-26's
   dynamic imprint: within 32 m exactly, beyond as the scattered wave carrying its energy out across the whole
   5 km scene (bounded), under the FTA's urban background over 8-100 Hz and a truck over a bump 15 m away for the whole
   pass, at whichever frequency the chamber answers most.
4. Bright points and the genie: a point as bright as the brightest on open plateau (12 dB, P2-23) and a corner reflector
   (50 dB) over the imprint's peak; and the method handed the ground's reflectivity (which no single image has), every
   pixel 30 dB over the receiver noise (the whole span an image has, P2-11).
5. The claimed deep structure, by P2-04's scaling.

Reported for each: F, the most any method's detection rate can exceed its false-alarm rate, and how many times harder
the ground would have to shake for a detection at 95% with 5% false alarms (d' = 3.29).
"""
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from katabasis.runs import RESULTS, Run, load
from sarsim import synthesize
from sarsim import information as inf
from sarsim.acquisition import DwellGeometry
from sarsim.scene import Scatterers

RID = 'p2_25_information_bound'
SITES = Path(__file__).resolve().parents[2] / 'sites'
KERNELS = RESULTS / 'cache' / 'p2_04_kernels.npz'
TAPER_STATIC = (30.0, 39.0)          # m, as P2-07: the static kernels tapered to zero before their grid's edge
TAPER_DYNAMIC = (32.0, 39.0)         # m, P2-26's maps (recorded to 39.5 m)
FAR_FROM = 32.0                      # m: the scattered wave counted from here out by the far-field bound
H = 0.5                              # m, the radar-frame grid
HALF = 45.0                          # m, its half-width
D_RELIABLE = 3.29                    # d' for 95% detection at 5% false alarms
SNR_GENIE_DB = 30.0
MC_REAL = 600
MC_SYNTH = 80


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def summary(F):
    F = float(F)
    return {'fisher': F, 'edge': float(inf.edge(F)), 'shake_factor_for_reliable': D_RELIABLE / np.sqrt(F) if F > 0 else None}


# ------------------------------------------------------------------ 2. checks

def check_covariance(g):
    """The exact double sum against tr(C+ dC C+ dC) built scatterer by scatterer, on a 32 x 8 image."""
    rng = np.random.default_rng(0)
    shape = (32, 8)
    na, nr = shape
    N = na * nr
    in_a, in_r, nu = inf.band_masks(g, shape)
    kr = np.fft.fftfreq(nr, d=g.dr)
    band = (in_a[:, None] & in_r[None, :]).ravel()
    t = g.nu_to_time(nu)
    Xa, Xr = np.meshgrid(np.arange(na) * g.dx, np.arange(nr) * g.dr, indexing='ij')
    K1, K2 = rng.standard_normal(shape), rng.standard_normal(shape)
    fm = 900.0                                        # fast enough to vary across so small a band
    comps = [(K1, lambda tt: np.cos(2 * np.pi * fm * tt + 0.3)), (K2, lambda tt: np.sin(2 * np.pi * 0.37 * fm * tt))]
    F_sum = inf.fisher_exact(comps, g, shape)
    FA, FR = np.meshgrid(nu, kr, indexing='ij')
    E = np.exp(-2j * np.pi * (np.stack([FA.ravel(), FR.ravel()], 1) @ np.stack([Xa.ravel(), Xr.ravel()], 1).T)) / np.sqrt(N)
    Phi = np.repeat(np.stack([sum(K * tau(t[i]) for K, tau in comps).ravel() for i in range(na)]), nr, axis=0)
    B = E * band[:, None]
    dB = -1j * Phi * B
    C0, dC = B @ B.conj().T, dB @ B.conj().T + B @ dB.conj().T
    Cp = np.linalg.pinv(C0, rcond=1e-10, hermitian=True)
    F_cov = float(np.trace(Cp @ dC @ Cp @ dC).real)
    shared = inf.fisher_exact([(np.full(shape, 0.7), lambda tt: np.cos(2 * np.pi * fm * tt))], g, shape)
    return {'image': list(shape), 'fisher_sum': F_sum, 'fisher_covariance': F_cov,
            'relative_difference': abs(F_sum - F_cov) / F_cov, 'shared_motion_fisher': shared}


def bump(g, shape, L, amp=1e-3):
    """A Gaussian bump of line-of-sight displacement, width L (m), on an image's own pixel grid (and its ground grid)."""
    na, nr = shape
    xa = (np.arange(na) - na // 2) * g.dx
    xg = (np.arange(nr) - nr // 2) * g.dr / np.sin(g.theta)
    Xa, Xg = np.meshgrid(xa, xg, indexing='ij')
    return amp * np.exp(-(Xa ** 2 + Xg ** 2) / (2 * L ** 2))


def check_continuum(g):
    """The continuum formula (used on the real geometry) against the exact sum on the image's own grid, with each
    row's zero-Doppler time, for a bump in the paired-echo and in the stretch regimes."""
    rows = []
    k0 = 4 * np.pi / g.lam
    shape = (1024, 48)
    dg = g.dr / np.sin(g.theta)
    xa = (np.arange(shape[0]) - shape[0] // 2) * g.dx
    for L, fm in ((1.5, 2.0), (1.5, 10.0), (4.0, 0.5)):
        K = bump(g, shape, L)
        Fe = 0.0
        for psi in (0.0, np.pi / 2):                     # the mean over psi of A + Re(B e^{2 i psi})
            al = (2 * np.pi * fm * xa / g.V + psi)[:, None]
            comps = [(k0 * K * np.cos(al), lambda t, fm=fm: np.cos(2 * np.pi * fm * t)),
                     (-k0 * K * np.sin(al), lambda t, fm=fm: np.sin(2 * np.pi * fm * t))]
            Fe += 0.5 * inf.fisher_exact(comps, g, shape)
        Fc = inf.fisher_sinusoid(K, (g.dx, dg), fm, g)
        rows.append({'width_m': L, 'f_hz': fm, 'paired_echo_m': fm * g.V / g.Ka, 'fisher_exact': Fe,
                     'fisher_continuum': Fc, 'ratio': Fc / Fe, 'bound': inf.fisher_bound(K, (g.dx, dg), g),
                     'slow_limit': inf.fisher_slow(2 * np.pi * fm * K, (g.dx, dg), g)})
        print(f"  continuum check L={L} f={fm}: ratio {Fc / Fe:.4f}", flush=True)
    return rows


def scatterer_scene(g, shape, rng, per_cell=2):
    Lx, Lr = shape[0] * g.dx, shape[1] * g.dr
    n = per_cell * int((Lx / g.resolution) * (Lr / (0.886 / g.kr_band)))
    x, r = rng.uniform(-Lx / 2, Lx / 2, n), rng.uniform(-Lr / 2, Lr / 2, n)
    z = np.zeros(n)
    return Scatterers(x=x, y=r / np.sin(g.theta), z=z, amp=np.ones(n), phase=rng.uniform(0, 2 * np.pi, n),
                      iso=np.ones(n), flash=z, nu0=z, sig_nu=np.ones(n), vib_amp=z, vib_freq=z, vib_phase=z,
                      label=np.zeros(n, int))


def check_attained(g):
    """The detector T = (2 / s2) Im sum K z conj(z_tau) on images with and without a boosted imprint: its spread without
    and its shift with, against sqrt(F) and a F. Case A: a 1.5 m bump at 2 Hz (paired echoes); B: an 8 m bump at 2 Hz
    (the stretch, the regime of the chamber's microseism imprint); C: the whole scene moving together."""
    k0 = 4 * np.pi / g.lam
    out = {}
    cases = {'A_compact': ((1024, 48), 1.5, 2.0, 0.25), 'B_smooth': ((1024, 112), 8.0, 2.0, 0.19)}
    for name, (shape, L, fm, a) in cases.items():
        Kp = k0 * bump(g, shape, L)                                  # phase, radians
        tau = lambda t, fm=fm: np.cos(2 * np.pi * fm * t)
        t0 = time.time()
        F = inf.fisher_exact([(Kp, tau)], g, shape)
        rng = np.random.default_rng(5 if name.startswith('A') else 6)
        T0 = np.array([inf.score(inf.speckle_image(g, shape, rng), Kp, tau, g) for _ in range(MC_REAL)])
        T1 = np.array([inf.score(inf.speckle_image(g, shape, rng, Kp, tau, a), Kp, tau, g) for _ in range(MC_REAL)])
        row = {'image': list(shape), 'width_m': L, 'f_hz': fm, 'a': a, 'peak_phase_rad': float(a * Kp.max()),
               'fisher_at_a1': F, 'deflection_expected': float(a * np.sqrt(F)),
               'spread_over_sqrtF': float(T0.std() / np.sqrt(F)), 'spread_se': float(np.sqrt(0.5 / MC_REAL)),
               'shift_over_aF': float((T1.mean() - T0.mean()) / (a * F)),
               'shift_se': float(np.sqrt(2 * F / MC_REAL) / (a * F)),
               'estimator_mean_over_a': float(T1.mean() / F / a), 'estimator_var_times_F': float(np.var(T1 / F) * F),
               'realisations': MC_REAL,
               'slow_limit': inf.fisher_slow(2 * np.pi * fm * bump(g, shape, L), (g.dx, g.dr / np.sin(g.theta)), g)}
        if name.startswith('A'):
            # the lab's synthesizer: scatterers placed at random, two per cell, moved through its own motion callback
            amp = 1e-3
            d_amp = lambda x, y, L=L: amp * np.exp(-(x ** 2 + y ** 2) / (2 * L ** 2))
            motion = lambda x, y, z, t, a=a, fm=fm: a * d_amp(x, y)[None, :] * np.cos(2 * np.pi * fm * t)[:, None]
            S0, S1 = [], []
            for k in range(MC_SYNTH):
                sc = scatterer_scene(g, shape, np.random.default_rng(1000 + k))
                S0.append(inf.score(synthesize(sc, g, shape, dtype=np.complex128), Kp, tau, g))
                S1.append(inf.score(synthesize(sc, g, shape, motion=motion, dtype=np.complex128), Kp, tau, g))
            S0, S1 = np.array(S0), np.array(S1)
            row['synthesizer'] = {'realisations': MC_SYNTH, 'spread_over_sqrtF': float(S0.std() / np.sqrt(F)),
                                  'spread_se': float(np.sqrt(0.5 / MC_SYNTH)),
                                  'paired_shift_over_aF': float(np.mean(S1 - S0) / (a * F)),
                                  'paired_shift_se': float(np.std(S1 - S0) / np.sqrt(MC_SYNTH) / (a * F))}
        row['runtime_s'] = round(time.time() - t0, 1)
        out[name] = row
        print(f"  attained {name}: spread {row['spread_over_sqrtF']:.3f}, shift {row['shift_over_aF']:.3f}", flush=True)
    shape = (1024, 48)
    rng = np.random.default_rng(7)
    tau = lambda t: np.cos(2 * np.pi * 2.0 * t)
    Kc = np.full(shape, 0.3)
    Tsh = [inf.score(inf.speckle_image(g, shape, rng, Kc, tau, 1.0), Kc, tau, g) for _ in range(20)]
    out['C_shared'] = {'image': list(shape), 'phase_rad': 0.3, 'fisher': inf.fisher_exact([(Kc, tau)], g, shape),
                       'largest_score': float(np.max(np.abs(Tsh)))}
    return out


# ------------------------------------------------------------------ 3. the chamber

def static_kernels():
    """P2-04's line-of-sight-ready kernels (ENU displacement per unit strain for eps_ee, eps_nn, eps_en) on its 2 m grid:
    read from P2-04's run if present, else from the cache, else recomputed with P2-04's own code (three relaxations)."""
    for path in (RESULTS / 'p2_04_chamber_imprint' / 'kernels.npz', KERNELS):
        if path.is_file():
            K = np.load(path)
            return {k: K[k] for k in ('xy', 'Kxx', 'Kyy', 'Kxy')}, str(path.relative_to(RESULTS))
    p204 = module('p204', Path(__file__).with_name('p2_04_chamber_imprint.py'))
    import json as _json
    from katabasis.seismic.arrays import surface_grid
    site, grid, full, _ = p204.media(1.0, 60.0)
    rec = surface_grid(grid, full.solid, p204.HALF, p204.SPACING)
    mat = _json.loads((SITES / 'materials.json').read_text())['materials']
    rock = next(m for m in (mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()])
                if m['id'] == 'limestone-mokattam')
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    vp, vs, rho = val('vp'), val('vs'), val('rho')
    mu, lam = rho * vs ** 2, rho * vp ** 2 - 2 * rho * vs ** 2
    nu, E = lam / (2 * (lam + mu)), mu * (3 * lam + 2 * mu) / (lam + mu)
    raw = {}
    for comp in ('xx', 'yy', 'xy'):
        raw[comp], info = p204.kernel(full, rec, {comp: p204.STRESS}, 1.0)
        print(f"  static kernel s_{comp}: {info['runtime_s']:.0f} s", flush=True)
    Kxx = E / (1 - nu ** 2) * (raw['xx'] + nu * raw['yy']) / p204.STRESS
    Kyy = E / (1 - nu ** 2) * (raw['yy'] + nu * raw['xx']) / p204.STRESS
    Kxy = 2 * mu * raw['xy'] / p204.STRESS
    KERNELS.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(KERNELS, xy=rec[:, :2], Kxx=Kxx, Kyy=Kyy, Kxy=Kxy)
    return {'xy': rec[:, :2], 'Kxx': Kxx, 'Kyy': Kyy, 'Kxy': Kxy}, 'recomputed'


def to_radar(M, xs, ys, g):
    """A map on a site grid (x east, y north, M[ix, iy], real or complex) resampled onto the radar frame: azimuth along
    the track, ground range away from the radar, H apart within HALF of the axis; zero off the map."""
    a, r = g.along_track_en, g.ground_range_en
    u = np.arange(-HALF, HALF + H / 2, H)
    A, R = np.meshgrid(u, u, indexing='ij')
    pts = np.stack([(A * a[0] + R * r[0]).ravel(), (A * a[1] + R * r[1]).ravel()], 1)
    f = lambda Z: RegularGridInterpolator((xs, ys), Z, bounds_error=False, fill_value=0.0)(pts).reshape(A.shape)
    return f(M.real) + 1j * f(M.imag) if np.iscomplexobj(M) else f(M)


def taper(xs, ys, r0, r1):
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    return np.clip((r1 - np.hypot(X, Y)) / (r1 - r0), 0, 1)


def static_maps(kern, g, azimuths):
    """LOS displacement per unit strain (slant-range increase, m) for a Rayleigh wave travelling towards each azimuth,
    tapered and resampled onto the radar frame: P2-04's rayleigh_los."""
    xy = kern['xy']
    xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
    ix, iy = np.searchsorted(xs, xy[:, 0]), np.searchsorted(ys, xy[:, 1])
    los = np.asarray(g.los_enu)
    tp = taper(xs, ys, *TAPER_STATIC)
    maps = []
    for phi in azimuths:
        s, c = np.sin(phi), np.cos(phi)
        u = kern['Kxx'] * s * s + kern['Kyy'] * c * c + kern['Kxy'] * s * c
        M = np.zeros((len(xs), len(ys)))
        M[ix, iy] = -(u @ los)
        maps.append(to_radar(M * tp, xs, ys, g))
    return maps, (xs, ys, ix, iy)


def static_cases(g, kern, p204, amb):
    """Each ambient band below 8 Hz: the imprint d = eps0 K_phi(x) cos(...), eps0 = hv V_amp / c with V_amp = sqrt(2) x
    the band's rms, F averaged over the waves' directions (every direction at the whole band's power: waves from all
    sides with independent phases add their informations)."""
    az = np.deg2rad(np.arange(0, 180, 7.5))           # P2-04's directions (a direction and its reverse give the same)
    maps, _ = static_maps(kern, g, az)
    hv, c_rock = p204['hv'], p204['rayleigh_speed_m_s']
    reg = amb['regional']
    bands = [(row['case'], row['vertical_velocity_m_s'], row['frequency_hz'], row['phase_speed_m_s'], row['source'])
             for row in p204['cases'][:3]]
    bands.append(('Giza 3-8 Hz, measured', reg['vertical_3_8_hz']['value'], 5.0, c_rock,
                  "M1-01 median; the rock's own Rayleigh speed, as P2-04 takes for 1-3 Hz"))
    rows = []
    for label, V, f, c, src in bands:
        eps0 = hv * np.sqrt(2) * V / c
        F = np.mean([inf.fisher_sinusoid(eps0 * M, (H, H), f, g) for M in maps])
        bound = np.mean([inf.fisher_bound(eps0 * M, (H, H), g) for M in maps])
        slow = np.mean([inf.fisher_slow(2 * np.pi * f * eps0 * M, (H, H), g) for M in maps])
        peak_phase = 4 * np.pi / g.lam * eps0 * max(np.abs(M).max() for M in maps)
        # the largest stretch of the texture along track: R / V_s times the velocity amplitude's along-track gradient
        stretch = g.R0 / g.V_platform * max(np.abs(np.gradient(2 * np.pi * f * eps0 * M, H, axis=0)).max() for M in maps)
        rows.append({'case': label, 'rms_velocity_m_s': V, 'f_hz': f, 'phase_speed_m_s': c, 'strain_amplitude': eps0,
                     'imprint_peak_phase_rad': peak_phase, 'stretch_peak': stretch, 'fisher_bound_any_frequency': bound,
                     'fisher_slow_limit': slow, 'source': src, **summary(F)})
        print(f"  {label}: F {F:.3g}, edge {inf.edge(F):.2g}", flush=True)
    return rows, maps


def dynamic_information(g, p226):
    """F per unit incident vertical velocity amplitude, per frequency, for each P2-26 case: the tapered map within 39 m
    exactly (averaged over its directions), and the scattered wave beyond FAR_FROM by the bound 4 N <Phi^2>, its energy
    per metre of radius (the far-field cross-width) taken as the largest over the outer annuli and carried unattenuated
    to the scene's edge."""
    maps = np.load(RESULTS / 'p2_26_imprint_spectrum' / 'maps.npz')
    src = json.loads((SITES / 'acquisitions' / 'giza-20250827.json').read_text())['source']
    rows_, cols_ = src['shape']
    scene_radius = 0.5 * min(rows_ * g.dx, cols_ * g.dr / np.sin(g.theta))
    k0 = 4 * np.pi / g.lam
    out = {}
    for name in ('surface', 'surface_favourable', 'p_below', 's_below'):
        xs, ys, fs, Hm = maps[f'{name}_x'], maps[f'{name}_y'], maps[f'{name}_f'], maps[f'{name}_H']
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        R = np.hypot(X, Y)
        dA = float(xs[1] - xs[0]) ** 2
        tp = taper(xs, ys, *TAPER_DYNAMIC)
        near, far, sigma = [], [], []
        for k, f in enumerate(fs):
            Kd = Hm[..., k] / (2j * np.pi * f)                       # displacement per unit incident velocity amplitude
            near.append(np.mean([inf.fisher_sinusoid(to_radar(Kd[d] * tp, xs, ys, g), (H, H), f, g)
                                 for d in range(Kd.shape[0])]))
            E = [float((np.abs(Hm[:, (R >= r) & (R < r + 5), k]) ** 2).sum(axis=1).mean()) * dA / 5.0
                 for r in (25.0, 30.0, 35.0)]                        # m^2 of |H|^2 per metre of radius
            s = max(E)
            sigma.append(s)
            far.append(2 * inf.cells_per_m2(g) * k0 ** 2 * s * (scene_radius - FAR_FROM) / (2 * np.pi * f) ** 2)
        out[name] = {'f_hz': fs, 'fisher_near_per_v2': np.array(near), 'fisher_far_bound_per_v2': np.array(far),
                     'cross_width_m': np.array(sigma)}
        print(f"  dynamic {name}: done", flush=True)
    return out, scene_radius


def main():
    g = DwellGeometry.from_record('giza-20250827')
    amb = json.loads((SITES / 'ambient.json').read_text())
    p204, p226, p211, p223 = (load(r) for r in ('p2_04_chamber_imprint', 'p2_26_imprint_spectrum',
                                                'p2_11_every_way_in', 'p2_23_every_reader'))
    params = {'acquisition': 'giza-20250827', 'grid_m': H, 'half_width_m': HALF, 'taper_static_m': TAPER_STATIC,
              'taper_dynamic_m': TAPER_DYNAMIC, 'far_field_from_m': FAR_FROM, 'reliable_dprime': D_RELIABLE,
              'genie_snr_db': SNR_GENIE_DB, 'monte_carlo': {'model': MC_REAL, 'synthesizer': MC_SYNTH},
              'noise_free': True, 'attenuation': 'none'}
    with Run(RID, 'What one image can hold about the chamber, whatever reads it', params) as run:
        checks = {'covariance': check_covariance(g)}
        print(f"  covariance check: {checks['covariance']['relative_difference']:.1e}", flush=True)
        checks['continuum'] = check_continuum(g)
        checks['attained'] = check_attained(g)

        kern, kern_src = static_kernels()
        maps204, (xs, ys, ix, iy) = static_maps(kern, g, np.deg2rad(np.arange(0, 180, 7.5)))
        # the recomputed kernels against P2-04's published rms map (same directions, untapered)
        los = np.asarray(g.los_enu)
        rms = np.sqrt(np.mean([(-(kern['Kxx'] * np.sin(p) ** 2 + kern['Kyy'] * np.cos(p) ** 2
                                  + kern['Kxy'] * np.sin(p) * np.cos(p)) @ los) ** 2
                               for p in np.deg2rad(np.arange(0, 180, 7.5))], axis=0))
        pub = np.asarray(p204['imprint_map_los_per_strain_m']).ravel()
        kernel_check = {'source': kern_src, 'peak_m_per_strain': float(rms.max()),
                        'published_peak_m_per_strain': float(pub.max()),
                        'largest_difference_over_peak': float(np.max(np.abs(rms - pub)) / pub.max())}
        static, _ = static_cases(g, kern, p204, amb)

        dyn, scene_radius = dynamic_information(g, p226)
        cul = amb['cultural']
        v_urban, band = cul['urban_background']['value'], cul['band_hz']['value']
        v_truck = cul['bus_or_truck_over_bump']['value']
        dynamic = []
        for name, d in dyn.items():
            f = np.asarray(d['f_hz'])
            tot = d['fisher_near_per_v2'] + d['fisher_far_bound_per_v2']
            sel = (f >= band[0]) & (f <= band[1])
            df = float(np.mean(np.diff(f)))
            psd = v_urban ** 2 / (band[1] - band[0])                  # (m/s)^2 / Hz, spread evenly over the band
            F_urban = float(np.sum(tot[sel] * 2 * psd * df))
            F_urban_near = float(np.sum(d['fisher_near_per_v2'][sel] * 2 * psd * df))
            k = int(np.argmax(np.where(sel, tot, 0)))
            F_truck = float(tot[k] * 2 * v_truck ** 2)
            F_truck_near = float(np.max(np.where(sel, d['fisher_near_per_v2'], 0)) * 2 * v_truck ** 2)
            dynamic.append({'case': name, 'urban': {**summary(F_urban), 'fisher_near_only': F_urban_near,
                                                    'rms_velocity_m_s': v_urban, 'band_hz': band},
                            'truck': {**summary(F_truck), 'fisher_near_only': F_truck_near, 'worst_f_hz': float(f[k]),
                                      'rms_velocity_m_s': v_truck},
                            'f_hz': f, 'fisher_per_v2': tot, 'fisher_near_per_v2': d['fisher_near_per_v2'],
                            'cross_width_m': d['cross_width_m']})
            print(f"  {name}: urban edge {inf.edge(F_urban):.2g}, truck edge {inf.edge(F_truck):.2g} "
                  f"at {f[k]:.0f} Hz", flush=True)

        # 4. bright points over the imprint's peak, and the genie
        T = g.aperture_time
        t = np.linspace(-T / 2, T / 2, 2001)
        scr_nat = 10 ** (p223['satellite']['natural']['brightest_db'] / 10)
        scr_cr = 10 ** (p223['manifest']['params']['reflector_scr_db'] / 10)
        k0 = 4 * np.pi / g.lam
        micro = static[0]
        d_peak_micro = micro['imprint_peak_phase_rad'] / k0
        sb = next(d for d in dynamic if d['case'] == 'surface')
        maps226 = np.load(RESULTS / 'p2_26_imprint_spectrum' / 'maps.npz')
        f226 = maps226['surface_f']
        kt = int(np.argmin(np.abs(f226 - sb['truck']['worst_f_hz'])))
        H_peak = float(np.abs(maps226['surface_H'][..., kt]).max())
        d_peak_truck = H_peak * np.sqrt(2) * v_truck / (2 * np.pi * f226[kt])
        points = []
        for label, scr in (('the brightest point on open plateau (P2-23)', scr_nat), ('a corner reflector (P2-23)', scr_cr)):
            for case, dpk, fq in (('microseisms', d_peak_micro, micro['f_hz']),
                                  ('truck over a bump, 15 m', d_peak_truck, float(f226[kt]))):
                F = np.mean([inf.fisher_point(scr, k0 * dpk * np.cos(2 * np.pi * fq * t + psi), t)
                             for psi in (0.0, np.pi / 2)])
                points.append({'point': label, 'scr_db': 10 * np.log10(scr), 'shaking': case, 'f_hz': fq,
                               'peak_displacement_m': dpk, **summary(F)})
        snr = 10 ** (SNR_GENIE_DB / 10)
        genie = [{'case': r['case'], **summary(0.5 * snr * r['fisher_bound_any_frequency'])} for r in static]
        # the truck at its worst frequency: the bound 4 N <Phi^2> on the near map plus the far-field bound, times SNR / 2
        Kd = maps226['surface_H'][..., kt] / (2j * np.pi * f226[kt])
        xs6, ys6 = maps226['surface_x'], maps226['surface_y']
        tp = taper(xs6, ys6, *TAPER_DYNAMIC)
        near_bound = np.mean([inf.fisher_bound(to_radar(Kd[q] * tp, xs6, ys6, g), (H, H), g) for q in range(len(Kd))])
        far_bound = sb['fisher_per_v2'][kt] - sb['fisher_near_per_v2'][kt]
        genie.append({'case': 'truck over a bump, 15 m (surface wave)',
                      **summary(0.5 * snr * (near_bound + far_bound) * 2 * v_truck ** 2)})

        # 5. the claimed deep structure: the static imprint scales as volume / depth^2 (P2-04), and in the stretch regime
        # the information depends on the peak alone (sum of |grad v|^2 over the ground does not change with width)
        rel = p204['claimed'][0]['relative_to_bench']
        claim = {'structure': p204['claimed'][0]['structure'], 'relative_imprint': rel,
                 'microseisms': summary(micro['fisher'] * rel ** 2),
                 'noisiest': summary(static[1]['fisher'] * rel ** 2)}
        pen = p211['A_penetration']['rows'][0]['bands']['X']['two_way_loss_db_if_sand']
        # the air inside the room: the solver's walls are traction-free; air resonating in the room pushes on them with a
        # pressure its impedance allows, a fraction (rho c)_air / (rho c)_rock of the rock's own stress per unit of its Q
        mats = json.loads((SITES / 'materials.json').read_text())['materials']
        mats = {m['id']: m for m in (mats if isinstance(mats, list) else [dict(v, id=k) for k, v in mats.items()])}
        z = lambda m: mats[m]['rho']['value'] * mats[m]['vp']['value']
        air = {'impedance_air_rayl': z('air'), 'impedance_rock_rayl': z('limestone-mokattam'),
               'ratio': z('air') / z('limestone-mokattam'), 'q_for_one_percent': 0.01 / (z('air') / z('limestone-mokattam')),
               'room_air_resonance_hz': mats['air']['vp']['value'] / (2 * 6.0)}

        ub = next(d for d in dynamic if d['case'] == 'surface')
        fav = next(d for d in dynamic if d['case'] == 'surface_favourable')
        at = checks['attained']
        worst = max([r['edge'] for r in static] + [d['urban']['edge'] for d in dynamic])
        finding = (
            f"However one image of ordinary ground is read, it holds almost nothing of the chamber. The chamber reaches the "
            f"image only through its echo, at least {pen:.0f} dB down at X band (P2-11), or through the surface's motion "
            f"during the pass; every method is a function of the image, and for ordinary ground (fully developed speckle) "
            f"the image's statistics fix how much it holds. Motion shared "
            f"by the whole scene holds exactly none: the texture's statistics do not change. What survives is the stretch "
            f"along track that the imprint's gradient makes, and it sets a ceiling on any method's detection rate over its "
            f"false-alarm rate. Under Giza's measured microseisms that ceiling is {static[0]['edge']:.0e}; at the "
            f"noisiest stations on Earth {static[1]['edge']:.0e}; at 1-3 and 3-8 Hz as measured {static[2]['edge']:.0e} "
            f"and {static[3]['edge']:.0e}. Above 6 Hz the chamber scatters the waves (P2-26), and under the FTA's urban "
            f"background for the whole pass, the scattered wave counted across the whole 5 km scene, the ceiling is "
            f"{ub['urban']['edge']:.0e} ({fav['urban']['edge']:.0e} for a room whose roof is 5 m down). Only a truck "
            f"over a bump 15 m from the chamber for the whole 25 s, at the frequency it answers most "
            f"({ub['truck']['worst_f_hz']:.0f} Hz), lifts it to {ub['truck']['edge']:.0e}. A reliable detection "
            f"(95% at 5% false alarms) under the microseisms needs the ground shaking "
            f"{static[0]['shake_factor_for_reliable']:.0e} times harder. The bound is what the image holds, not a loose "
            f"ceiling: the best detector, run on simulated speckle and on the lab's synthesizer, scatters by "
            f"{at['A_compact']['spread_over_sqrtF']:.2f} and {at['A_compact']['synthesizer']['spread_over_sqrtF']:.2f} "
            f"of the bound's sqrt(F) and moves by {at['A_compact']['shift_over_aF']:.2f} and "
            f"{at['A_compact']['synthesizer']['paired_shift_over_aF']:.2f} of its a F.")
        run.save({'checks': checks, 'kernels': kernel_check, 'static': static, 'dynamic': dynamic,
                  'scene_radius_m': scene_radius, 'points': points, 'genie': genie, 'claimed': claim,
                  'penetration_db_x_band': pen, 'air_in_the_room': air, 'worst_ambient_edge': worst, 'cells_per_m2': inf.cells_per_m2(g),
                  'stretch_seconds': g.R0 / g.V_platform, 'paired_echo_m_per_hz': g.V / g.Ka, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
