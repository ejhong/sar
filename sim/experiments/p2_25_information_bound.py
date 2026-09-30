"""P2-25 · What one image can hold about the chamber, whatever reads it: a bound within a stated model.

    uv run python experiments/p2_25_information_bound.py

P2-03 to P2-24 test methods by running them. This bounds them all at once, with every method not yet written, within a
stated measurement model. A radar image is made of the echo from what the wave reaches, the surface, and the surface's
motion during the pass; every method is a function of the image, so none can hold more about the chamber than the image
does. For ordinary ground modelled as fully developed speckle the image's statistics fix how much that is
(sarsim.information): for any test, detection rate - false-alarm rate <= TV <= min(sqrt(KL / 2), sqrt(1 - exp(-KL)))
(for two equally likely cases the best accuracy is 1/2 + TV/2), and KL for the finite change from no chamber to the
chamber is bounded rigorously, with an explicit remainder, by the Fisher information F.

1. The routes. The echo from the chamber: P2-11 puts the bench chamber's roof at least 347 dB down at X band; noted,
   not pursued. The surface: what any photograph shows. The motion: bounded here, generously; the method is told the
   ambient motion itself (which only raises the information), receiver noise is left out (likewise), the rock has no
   attenuation.
2. The bound, checked. The double sum against a brute-force covariance (tiny image); the FFT form on the image's own
   pixels against the double sum; the finite change's exact KL against a^2 F / 2 and against the rigorous bound;
   detection: the exact likelihood-ratio test and the score test, their detection curves and their largest
   detection-minus-false-alarm rate, beside the bound and beside what theory says the best test reaches (0.40 sqrt F for
   a weak signal against the bound's 0.50 sqrt F); estimation, separately: the score's spread and shift against sqrt F
   and a F on the model's speckle and on the lab's synthesizer; shared motion, exactly zero; and the model's
   Doppler-to-time relation (motion multiplied into the spectrum at each bin's stationary time) against the exact
   physics (motion multiplied into each pulse), on a short aperture at the real chirp rate.
3. The chamber, the one-chamber bench (6 m, 15 m down), on the real Giza dwell (P2-01), computed exactly on the image's
   own pixels. Below 8 Hz its imprint is P2-04's static one (the kernels recomputed if absent), averaged over the
   waves' directions, under Giza's regional microseism level (measured 67 km east; the level at the pyramids is
   unmeasured), the noisiest stations on Earth, and the measured 1-3 and 3-8 Hz levels. Above 6 Hz it is P2-26's
   dynamic imprint: within 39 m exactly, beyond as the scattered wave carrying its energy out across the whole 5 km
   scene (bounded, unattenuated), under the FTA's urban background over 8-100 Hz and a truck over a bump 15 m away for
   the whole pass, at whichever frequency and from whichever side the chamber answers most.
4. Bright points and the genie: a point as bright as the brightest on open plateau (12 dB, P2-23) and a corner reflector
   (50 dB) over the imprint's peak; and the method handed the ground's reflectivity (which no single image has), every
   pixel 30 dB over the receiver noise (the whole span an image has, P2-11).
5. The claimed deep structure, by P2-04's scaling.

Reported for each: F, the rigorous bound on KL, the ceiling on detection rate - false-alarm rate (and on accuracy above a
coin toss, half of it), what the best test reaches for a weak signal, and the local extrapolation factor: how many times
stronger the imprint would have to be for the linearised information to reach a reliable detection (d' = 3.29). The last
is an extrapolation within the model, not a physical shaking requirement: long before it, small-phase and linear-rock
assumptions fail.
"""
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.fft import next_fast_len

from katabasis.runs import RESULTS, Run, load, memo
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
HALF = 45.0                          # m: the image's own pixels within this of the axis (ground metres)
D_RELIABLE = 3.29                    # d' for 95% detection at 5% false alarms
SNR_GENIE_DB = 30.0
MC_REAL = 600
MC_SYNTH = 80
MC_DETECT = 4000


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def summary(F, kl=None):
    """The quantities reported for one case: F, the KL bound (given, or F / 2 to first order where the remainder is
    negligible and has been checked), the ceiling on detection - false alarm, half of it on accuracy, the best test's
    weak-signal value, and the local extrapolation factor."""
    F = float(F)
    kl = float(inf.kl_upper(F)) if kl is None else float(kl)
    tv = float(inf.tv_upper(kl))
    return {'fisher': F, 'kl_upper': kl, 'ceiling': tv, 'accuracy_above_half': tv / 2,
            'best_test_weak_signal': float(inf.best_test_tv(F)),
            'local_extrapolation_factor': D_RELIABLE / np.sqrt(F) if F > 0 else None}


def ambient_kl(parts, eps):
    """KL bound for a random superposition (independent phases) of components j, each with mean information A_j and
    largest information over its phase M_j: KL(mixture) <= E KL <= (E F + 2 eps sqrt(E F) + eps^2) / (2 (1 - x)),
    E F = sum A_j, x = sum sqrt(M_j) + eps bounding |E|_F in every realisation (triangle inequality)."""
    EF = float(sum(a for a, _ in parts))
    x = float(sum(np.sqrt(m) for _, m in parts)) + eps
    if x >= 1:
        return EF, np.inf, x
    return EF, (EF + 2 * eps * np.sqrt(EF) + eps ** 2) / (2 * (1 - x)), x


def eps_from_norms(l2, kmax, g):
    """The remainder bound from a map's norms alone (phase units on the image's pixels): with |K| <= kmax,
    sqrt(share sum |K|^(2n)) <= kmax^(n-1) sqrt(share) |K|_2, so eps <= sqrt(share) |K|_2 (exp(2 kmax) - 1 - 2 kmax) / kmax."""
    share = g.band_frac * g.band_frac_r
    if kmax <= 0:
        return 0.0
    return float(np.sqrt(share) * l2 * (np.expm1(2 * kmax) - 2 * kmax) / kmax)


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


def check_grid(g):
    """The FFT form on the image's own pixels (used for every number below) against the double sum."""
    shape = (64, 10)
    xa = (np.arange(shape[0]) - shape[0] / 2) * g.dx
    Xa, Xr = np.meshgrid(xa, (np.arange(shape[1]) - shape[1] / 2) * g.dr, indexing='ij')
    K = 1e-3 * np.exp(-(Xa ** 2 + (Xr / np.sin(g.theta)) ** 2) / 0.5) * (1 + 0.3j * np.sin(Xr))
    k0 = 4 * np.pi / g.lam
    rows = []
    for fm in (2.0, 40.0, 400.0):
        A, B = inf.fisher_grid(K, g, fm, pad=shape)
        Fs = []
        for psi in np.linspace(0, np.pi, 8, endpoint=False):
            Kp = 0.5 * k0 * K * np.exp(2j * np.pi * fm * xa / g.V)[:, None] * np.exp(1j * psi)
            Fs.append(inf.fisher_exact([(Kp, lambda t, fm=fm: np.exp(2j * np.pi * fm * t)),
                                        (np.conj(Kp), lambda t, fm=fm: np.exp(-2j * np.pi * fm * t))], g, shape))
        rows.append({'f_hz': fm, 'mean_over_psi_sum': float(np.mean(Fs)), 'grid_A': A,
                     'largest_sum': float(max(Fs)), 'grid_A_plus_B': A + B,
                     'relative_difference_mean': abs(np.mean(Fs) - A) / A})
    return rows


def small_problem(g, L, fm, shape=(128, 16)):
    na, nr = shape
    xa = (np.arange(na) - na // 2) * g.dx
    xg = (np.arange(nr) - nr // 2) * g.dr / np.sin(g.theta)
    Xa, Xg = np.meshgrid(xa, xg, indexing='ij')
    K = 4 * np.pi / g.lam * 1e-3 * np.exp(-(Xa ** 2 + Xg ** 2) / (2 * L ** 2))       # phase at a = 1
    return shape, K, [(K, lambda t, fm=fm: np.cos(2 * np.pi * fm * t))]


def check_finite_kl(g):
    """The finite change's exact KL (eigenvalues of the in-band covariance) against a^2 F / 2 and the rigorous bound."""
    rows = []
    for L, fm in ((0.4, 1.0), (0.8, 1.5)):
        shape, K, comps = small_problem(g, L, fm)
        F = inf.fisher_exact(comps, g, shape)
        k0 = 4 * np.pi / g.lam
        for a in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
            B = inf.band_operator(comps, g, shape, a)
            kl = inf.kl_exact(B @ B.conj().T)
            eps = inf.remainder(K / k0, g, shape, a=a)
            rows.append({'width_m': L, 'f_hz': fm, 'a': a, 'peak_phase_rad': float(a * K.max()), 'half_a2F': a * a * F / 2,
                         'kl_exact': kl, 'remainder_eps': eps, 'kl_rigorous_bound': float(inf.kl_upper(a * a * F, eps)),
                         'bound_holds': bool(kl <= float(inf.kl_upper(a * a * F, eps)))})
    return rows


def check_detection(g):
    """Detection, not estimation: on a small image, the exact likelihood-ratio test (the best test for the finite change)
    and the score test, each over MC_DETECT images per hypothesis; their detection curves and largest
    detection-minus-false-alarm rate, beside Pinsker's bound from the exact KL, the rigorous bound from F, and the
    weak-signal theory 2 Phi(d/2) - 1."""
    shape, K, comps = small_problem(g, 0.4, 1.0)
    F = inf.fisher_exact(comps, g, shape)
    k0 = 4 * np.pi / g.lam
    out = []
    rng = np.random.default_rng(31)
    for d_target in (0.5, 1.5):
        a = d_target / np.sqrt(F)
        B1 = inf.band_operator(comps, g, shape, a)
        C1 = B1 @ B1.conj().T
        kl = inf.kl_exact(C1)
        n_b = C1.shape[0]
        C1i = np.linalg.inv(C1)
        _, logdet = np.linalg.slogdet(C1)
        d = 1e-4                                            # the score's quadratic form, C'(0), by a central difference
        Bp, Bm = inf.band_operator(comps, g, shape, d), inf.band_operator(comps, g, shape, -d)
        Cd = (Bp @ Bp.conj().T - Bm @ Bm.conj().T) / (2 * d)
        Y0 = (rng.standard_normal((n_b, MC_DETECT)) + 1j * rng.standard_normal((n_b, MC_DETECT))) / np.sqrt(2)
        S = (rng.standard_normal((B1.shape[1], MC_DETECT)) + 1j * rng.standard_normal((B1.shape[1], MC_DETECT))) / np.sqrt(2)
        Y1 = B1 @ S
        lrt = lambda Y: -np.real(np.sum(np.conj(Y) * ((C1i - np.eye(n_b)) @ Y), axis=0)) - logdet
        scr = lambda Y: np.real(np.sum(np.conj(Y) * (Cd @ Y), axis=0))
        row = {'a': float(a), 'deflection': float(a * np.sqrt(F)), 'peak_phase_rad': float(a * K.max()), 'kl_exact': kl,
               'pinsker_from_exact_kl': float(inf.tv_upper(kl)),
               'rigorous_ceiling': float(inf.ceiling(a * a * F, inf.remainder(K / k0, g, shape, a=a))),
               'weak_signal_best_test': float(inf.best_test_tv(a * a * F)), 'realisations': MC_DETECT}
        for name, stat in (('likelihood_ratio', lrt), ('score', scr)):
            s0, s1 = stat(Y0), stat(Y1)
            thr = np.quantile(np.r_[s0, s1], np.linspace(0, 1, 401))
            pfa = np.array([(s0 > t).mean() for t in thr])
            pd = np.array([(s1 > t).mean() for t in thr])
            row[name] = {'largest_pd_minus_pfa': float(np.max(pd - pfa)), 'roc_pfa': pfa[::8], 'roc_pd': pd[::8]}
        row['se'] = float(np.sqrt(2 / MC_DETECT))
        out.append(row)
        print(f"  detection d={row['deflection']:.2f}: LRT {row['likelihood_ratio']['largest_pd_minus_pfa']:.3f}, "
              f"score {row['score']['largest_pd_minus_pfa']:.3f}, Pinsker(exact KL) {row['pinsker_from_exact_kl']:.3f}, "
              f"weak-signal best {row['weak_signal_best_test']:.3f}", flush=True)
    return out


def pulse_fisher(Ka, V, lam, T, fm, L, width, per_cell=2, cnr=10.0, band_only=True):
    """F about a sinusoidal slant-range motion K(u) cos(2 pi fm t + psi) of a strip of scatterers, computed two ways for
    the same scatterers and noise: the physics (motion multiplied into each pulse of chirp_k(t) = exp(-i pi Ka (t -
    tau_k)^2)) and the model (the chirp's spectrum multiplied at each bin's stationary time tau_k - f / Ka). Both are
    kept to the processed band |f| <= Ka T / 2, as the product is. Woodbury in scatterer space; mean over psi."""
    res = 0.886 * V / (Ka * T)
    u = np.arange(-L / 2, L / 2, res / per_cell)
    tau = u / V
    k0 = 4 * np.pi / lam
    Kmot = 1e-3 * np.exp(-u ** 2 / (2 * width ** 2))
    prf = 1.3 * (Ka * T + Ka * L / V + 4 * fm) + 50
    n = int(np.ceil(T * prf / 2)) * 2
    t = (np.arange(n) - n / 2) / prf
    f = np.fft.fftfreq(n, 1 / prf)
    band = np.abs(f) <= Ka * T / 2 if band_only else np.ones(n, bool)
    chirp = np.exp(-1j * np.pi * Ka * (t[:, None] - tau[None, :]) ** 2)
    Cf = np.fft.fft(chirp, axis=0) / np.sqrt(n)
    sig_n2 = per_cell * band.sum() / cnr
    out = {}
    for label in ('exact', 'model'):
        Fsum = 0.0
        for psi in (0.0, np.pi / 2):
            if label == 'exact':
                gt = k0 * Kmot[None, :] * np.cos(2 * np.pi * fm * t[:, None] + psi)
                dAf = np.fft.fft(-1j * gt * chirp, axis=0) / np.sqrt(n)
            else:
                ts = tau[None, :] - f[:, None] / Ka
                dAf = -1j * k0 * Kmot[None, :] * np.cos(2 * np.pi * fm * ts + psi) * Cf
            A, dA = Cf[band], dAf[band]
            G = A.conj().T @ A
            Q = np.linalg.inv(sig_n2 * np.eye(len(u)) + G)
            AdA = A.conj().T @ dA
            Gt, Ht = Q @ G, Q @ AdA
            Jt = (dA.conj().T @ dA - AdA.conj().T @ Q @ AdA) / sig_n2
            Fsum += 0.5 * float(2 * np.real(np.trace(Ht @ Ht)) + 2 * np.real(np.trace(Gt @ Jt)))
        out[label] = Fsum
    return out['exact'], out['model'], n, len(u)


def check_pulse_domain(g):
    """The model's Doppler-to-time relation against pulse-by-pulse physics. (a) The real chirp rate on a 0.5 s aperture
    over 0.5-150 Hz: the stationary-phase approximation is local, on the scale 1 / sqrt(Ka) = 13 ms. (b) Five cycles of
    motion across the aperture, as the microseisms make across the real dwell, at growing time-bandwidth product
    Ka T^2 (the real dwell's is 3.4e6, out of reach pulse by pulse): whether what differs is an edge effect that fades."""
    Ka = abs(g.Ka_signed)
    V = float(np.sqrt(Ka * g.lam * g.R0 / 2))
    rows_a = []
    for fm in (0.5, 5.0, 20.0, 60.0, 150.0):
        Fe, Fm, n, m = pulse_fisher(Ka, V, g.lam, 0.5, fm, 30.0, 5.0)
        rows_a.append({'f_hz': fm, 'fisher_exact_pulses': Fe, 'fisher_model': Fm, 'ratio_exact_over_model': Fe / Fm,
                       'pulses': n, 'scatterers': m})
        print(f"  pulse domain (0.5 s) f={fm}: exact/model {Fe / Fm:.4f}", flush=True)
    rows_b = []
    for T, Ka_s in ((1.0, 1000.0), (2.0, 1000.0), (4.0, 1000.0), (4.0, 4000.0)):
        fm = 5.0 / T
        res = 0.886 * V / (Ka_s * T)
        Fe, Fm, n, m = pulse_fisher(Ka_s, V, g.lam, T, fm, 24 * res, 3 * res)
        rows_b.append({'aperture_s': T, 'ka_hz_s': Ka_s, 'time_bandwidth': Ka_s * T * T, 'f_hz': fm,
                       'fisher_exact_pulses': Fe, 'fisher_model': Fm, 'ratio_exact_over_model': Fe / Fm, 'pulses': n,
                       'scatterers': m})
        print(f"  pulse domain (5 cycles) T={T} Ka={Ka_s}: exact/model {Fe / Fm:.4f}", flush=True)
    return {'real_rate_short_aperture': rows_a, 'five_cycles': rows_b,
            'stationary_time_scale_s': float(1 / np.sqrt(Ka)), 'real_time_bandwidth': float(Ka * g.aperture_time ** 2)}


def scatterer_scene(g, shape, rng, per_cell=2):
    Lx, Lr = shape[0] * g.dx, shape[1] * g.dr
    n = per_cell * int((Lx / g.resolution) * (Lr / (0.886 / g.kr_band)))
    x, r = rng.uniform(-Lx / 2, Lx / 2, n), rng.uniform(-Lr / 2, Lr / 2, n)
    z = np.zeros(n)
    return Scatterers(x=x, y=r / np.sin(g.theta), z=z, amp=np.ones(n), phase=rng.uniform(0, 2 * np.pi, n),
                      iso=np.ones(n), flash=z, nu0=z, sig_nu=np.ones(n), vib_amp=z, vib_freq=z, vib_phase=z,
                      label=np.zeros(n, int))


def check_estimator(g):
    """The score T = (2 / s2) Im sum K z conj(z_tau) on images with and without a boosted imprint: its spread without
    and its shift with, against sqrt(F) and a F, and the estimate T / F against the Cramer-Rao bound. This checks the
    Fisher information as an estimation quantity; detection is checked separately (check_detection). Case A: a 1.5 m bump
    at 2 Hz (paired echoes); B: an 8 m bump at 2 Hz (the stretch, the regime of the chamber's microseism imprint);
    C: the whole scene moving together."""
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
        print(f"  estimator {name}: spread {row['spread_over_sqrtF']:.3f}, shift {row['shift_over_aF']:.3f}", flush=True)
    shape = (1024, 48)
    rng = np.random.default_rng(7)
    tau = lambda t: np.cos(2 * np.pi * 2.0 * t)
    Kc = np.full(shape, 0.3)
    Tsh = [inf.score(inf.speckle_image(g, shape, rng, Kc, tau, 1.0), Kc, tau, g) for _ in range(20)]
    out['C_shared'] = {'image': list(shape), 'phase_rad': 0.3, 'fisher': inf.fisher_exact([(Kc, tau)], g, shape),
                       'largest_score': float(np.max(np.abs(Tsh)))}
    return out


# ------------------------------------------------------------------ 3. the chamber

def pixel_axes(g, half=HALF):
    """The image's own pixels within `half` ground metres of the axis: rows g.dx apart along track, columns g.dr apart in
    slant range (g.dr / sin(theta) of ground)."""
    dg = g.dr / np.sin(g.theta)
    na, nr = int(2 * half / g.dx), int(2 * half / dg)
    return (np.arange(na) - na / 2) * g.dx, (np.arange(nr) - nr / 2) * dg


def to_pixels(M, xs, ys, g, axes):
    """A site-frame map (x east, y north, M[ix, iy] on a regular grid) resampled by cubic splines onto the image's pixels
    (azimuth along the track, ground range away from the radar); zero off the map."""
    a, r = g.along_track_en, g.ground_range_en
    A, R = np.meshgrid(*axes, indexing='ij')
    E, N = A * a[0] + R * r[0], A * a[1] + R * r[1]
    ci = (E - xs[0]) / (xs[1] - xs[0])
    cj = (N - ys[0]) / (ys[1] - ys[0])
    f = lambda Z: map_coordinates(Z, [ci, cj], order=3, mode='constant', cval=0.0)
    return f(M.real) + 1j * f(M.imag) if np.iscomplexobj(M) else f(M)


def grid_pad(g, na, nr, f_m):
    """Padding for the exact form: the imprint and its paired echoes, displaced f_m V_g / |Ka| metres each way."""
    D = f_m * g.V / g.Ka
    return next_fast_len(int(na + 2.4 * D / g.dx) + 16), next_fast_len(2 * nr)


def taper(xs, ys, r0, r1):
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    return np.clip((r1 - np.hypot(X, Y)) / (r1 - r0), 0, 1)




def static_maps(kern, g, azimuths, axes):
    """LOS displacement per unit strain (slant-range increase, m) for a Rayleigh wave travelling towards each azimuth,
    tapered and resampled onto the image's pixels: P2-04's rayleigh_los."""
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
        maps.append(to_pixels(M * tp, xs, ys, g, axes))
    return maps


def static_cases(g, kern, p204, amb):
    """Each ambient band below 8 Hz: d = eps0 K_phi(x) cos(...), eps0 = hv V_amp / c with V_amp = sqrt(2) x the band's
    rms, the waves from every direction with independent phases, each carrying 1/n of the power. Exact on the pixels."""
    az = np.deg2rad(np.arange(0, 180, 7.5))           # P2-04's directions (a direction and its reverse give the same)
    axes = pixel_axes(g)
    maps = static_maps(kern, g, az, axes)
    hv, c_rock = p204['hv'], p204['rayleigh_speed_m_s']
    reg = amb['regional']
    bands = [(row['case'], row['vertical_velocity_m_s'], row['frequency_hz'], row['phase_speed_m_s'], row['source'])
             for row in p204['cases'][:3]]
    bands.append(('Giza 3-8 Hz, measured', reg['vertical_3_8_hz']['value'], 5.0, c_rock,
                  "M1-01 median; the rock's own Rayleigh speed, as P2-04 takes for 1-3 Hz"))
    rows = []
    n = len(maps)
    dg = g.dr / np.sin(g.theta)
    for label, V, f, c, src in bands:
        eps0 = hv * np.sqrt(2) * V / c
        pad = grid_pad(g, *maps[0].shape, f)
        parts, worst = [], np.zeros_like(maps[0])
        for M in maps:
            A, B = inf.fisher_grid(eps0 * M, g, f, pad=pad)
            parts.append((A / n, (A + B) / n))
            worst += np.abs(eps0 * M) / np.sqrt(n)
        eps = inf.remainder(worst, g, pad)
        EF, kl, x = ambient_kl(parts, eps)
        cont = np.mean([inf.fisher_sinusoid(eps0 * M, (g.dx, dg), f, g) for M in maps])
        bound = np.mean([inf.fisher_bound(eps0 * M, (g.dx, dg), g) for M in maps])
        slow = np.mean([inf.fisher_slow(2 * np.pi * f * eps0 * M, (g.dx, dg), g) for M in maps])
        peak_phase = 4 * np.pi / g.lam * eps0 * max(np.abs(M).max() for M in maps)
        stretch = g.R0 / g.V_platform * max(np.abs(np.gradient(2 * np.pi * f * eps0 * M, g.dx, axis=0)).max()
                                           for M in maps)
        rows.append({'case': label, 'rms_velocity_m_s': V, 'f_hz': f, 'phase_speed_m_s': c, 'strain_amplitude': eps0,
                     'imprint_peak_phase_rad': peak_phase, 'stretch_peak': stretch, 'remainder_eps': eps,
                     'norm_bound_x': x, 'fisher_continuum': cont, 'fisher_bound_any_frequency': bound,
                     'fisher_slow_limit': slow, 'source': src, **summary(EF, kl)})
        print(f"  {label}: F {EF:.3g} (continuum {cont:.3g}), ceiling {rows[-1]['ceiling']:.2g}", flush=True)
    return rows


def dynamic_information(g, maps_path, p226):
    """Per P2-26 case, per frequency and direction, per unit incident vertical velocity amplitude: the tapered map within
    39 m exactly on the image's pixels (A and the largest over phase, A + |B|), and the scattered wave beyond FAR_FROM by
    the bound 4 N <Phi^2> (its energy per metre of radius, the far-field cross-width, the largest over the outer annuli,
    carried unattenuated to the farthest corner of the scene). Combined per component by the triangle inequality on
    sqrt F."""
    maps = np.load(maps_path)
    src = json.loads((SITES / 'acquisitions' / 'giza-20250827.json').read_text())['source']
    rows_, cols_ = src['shape']
    scene_radius = 0.5 * float(np.hypot(rows_ * g.dx, cols_ * g.dr / np.sin(g.theta)))   # covers the corners
    k0 = 4 * np.pi / g.lam
    axes = pixel_axes(g)
    out = {}
    for name in ('surface', 'surface_favourable', 'p_below', 's_below'):
        t0 = time.time()
        xs, ys, fs, Hm = maps[f'{name}_x'], maps[f'{name}_y'], maps[f'{name}_f'], maps[f'{name}_H']
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        R = np.hypot(X, Y)
        dA = float(xs[1] - xs[0]) ** 2
        tp = taper(xs, ys, *TAPER_DYNAMIC)
        A_near, M_near, far, sigma, l2, peak, l2_far, p_far = [], [], [], [], [], [], [], []
        dpx = g.dx * g.dr / np.sin(g.theta)                           # ground area of one pixel
        for k, f in enumerate(fs):
            Kd = Hm[..., k] / (2j * np.pi * f)                      # displacement per unit incident velocity amplitude
            pad = None
            a_d, m_d, l_d, p_d = [], [], [], []
            for d in range(Kd.shape[0]):
                P = to_pixels(Kd[d] * tp, xs, ys, g, axes)
                pad = pad or grid_pad(g, *P.shape, f)
                A, B = inf.fisher_grid(P, g, f, pad=pad)
                a_d.append(A)
                m_d.append(A + B)
                l_d.append(k0 * float(np.sqrt(np.sum(np.abs(P) ** 2))))
                p_d.append(k0 * float(np.abs(P).max()))
            A_near.append(a_d)
            M_near.append(m_d)
            l2.append(l_d)
            peak.append(p_d)
            # each direction's own far field (an average over directions cannot support a worst-direction claim):
            # its energy per metre of radius, the largest over the outer annuli
            s = np.max([(np.abs(Hm[:, (R >= r) & (R < r + 5), k]) ** 2).sum(axis=1) * dA / 5.0
                        for r in (25.0, 30.0, 35.0)], axis=0)                       # [direction]
            sigma.append(s)
            far.append(2 * inf.cells_per_m2(g) * k0 ** 2 * s * (scene_radius - FAR_FROM) / (2 * np.pi * f) ** 2)
            # the far field's norms on the pixels, for the remainder: |K|^2 = sigma / (2 pi r) / (2 pi f)^2 in metres
            l2_far.append(k0 * np.sqrt(s * (scene_radius - FAR_FROM) / dpx) / (2 * np.pi * f))
            p_far.append(k0 * np.sqrt(s / (2 * np.pi * FAR_FROM)) / (2 * np.pi * f))
        out[name] = {'f_hz': fs, 'A_near': np.array(A_near), 'M_near': np.array(M_near), 'far_bound': np.array(far),
                     'cross_width_m': np.array(sigma), 'directions': Hm.shape[0], 'l2': np.array(l2),
                     'peak': np.array(peak), 'l2_far': np.array(l2_far), 'peak_far': np.array(p_far)}
        print(f"  dynamic {name}: {time.time() - t0:.0f} s", flush=True)
    return out, scene_radius


def dynamic_cases(dyn, g, amb, maps_path):
    """The urban background (a random spectrum spread evenly over its band, waves from every direction with independent
    phases) and a truck over a bump (one sinusoid, from its worst side, at its worst frequency, its phase random)."""
    cul = amb['cultural']
    v_urban, band = cul['urban_background']['value'], cul['band_hz']['value']
    v_truck = cul['bus_or_truck_over_bump']['value']
    maps = np.load(maps_path)
    rows = []
    for name, d in dyn.items():
        f = np.asarray(d['f_hz'])
        sel = (f >= band[0]) & (f <= band[1])
        df = float(np.mean(np.diff(f)))
        nd = d['directions']
        psd = v_urban ** 2 / (band[1] - band[0])
        parts, parts_near = [], []
        for k in np.nonzero(sel)[0]:
            amp2 = 2 * psd * df / nd                                 # each direction's share of the bin's power
            for j in range(nd):
                a = (np.sqrt(d['A_near'][k][j]) + np.sqrt(d['far_bound'][k][j])) ** 2 * amp2
                m = (np.sqrt(d['M_near'][k][j]) + np.sqrt(d['far_bound'][k][j])) ** 2 * amp2
                parts.append((a, m))
                parts_near.append((d['A_near'][k][j] * amp2, d['M_near'][k][j] * amp2))
        # the remainder: every component's amplitude added at every pixel, near maps and far field alike (safe)
        amps = np.sqrt(2 * psd * df / nd)
        l2_sum = amps * float(d['l2'][sel].sum() + d['l2_far'][sel].sum())
        k_sum = amps * float(d['peak'][sel].sum() + d['peak_far'][sel].sum())
        eps = eps_from_norms(l2_sum, k_sum, g)
        EF, kl, x = ambient_kl(parts, eps)
        EFn, kln, _ = ambient_kl(parts_near, eps)
        # truck: the worst frequency and side, the mean over its phase (random) and the largest (reported)
        tot = (np.sqrt(np.asarray(d['A_near'])) + np.sqrt(np.asarray(d['far_bound']))) ** 2          # [f, direction]
        k, jw = np.unravel_index(np.argmax(np.where(sel[:, None], tot, 0)), tot.shape)
        amp2 = 2 * v_truck ** 2
        eps_t = eps_from_norms(np.sqrt(amp2) * (d['l2'][k][jw] + d['l2_far'][k][jw]),
                               np.sqrt(amp2) * (d['peak'][k][jw] + d['peak_far'][k][jw]), g)
        Ft = float(tot[k, jw]) * amp2
        Ft_max = (np.sqrt(d['M_near'][k][jw]) + np.sqrt(d['far_bound'][k][jw])) ** 2 * amp2
        Ft_near = d['A_near'][k][jw] * amp2
        rows.append({'case': name,
                     'urban': {**summary(EF, kl), 'remainder_eps': eps, 'fisher_near_only': EFn,
                               'ceiling_near_only': float(inf.tv_upper(kln)),
                               'rms_velocity_m_s': v_urban, 'band_hz': band, 'norm_bound_x': x},
                     'truck': {**summary(Ft, float(inf.kl_upper(Ft, eps_t))), 'remainder_eps': eps_t,
                               'fisher_largest_over_phase': Ft_max,
                               'ceiling_largest_over_phase': float(inf.ceiling(Ft_max, eps_t)),
                               'fisher_near_only': Ft_near, 'ceiling_near_only': float(inf.ceiling(Ft_near)),
                               'worst_f_hz': float(f[k]), 'rms_velocity_m_s': v_truck},
                     'f_hz': f, 'fisher_near_per_v2': np.max(d['A_near'], axis=1), 'far_bound_per_v2': d['far_bound'],   # [f, direction]
                     'cross_width_m': d['cross_width_m']})
        print(f"  {name}: urban ceiling {rows[-1]['urban']['ceiling']:.2g}, truck {rows[-1]['truck']['ceiling']:.2g} "
              f"at {f[k]:.0f} Hz", flush=True)
    return rows


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




def main():
    g = DwellGeometry.from_record('giza-20250827')
    amb = json.loads((SITES / 'ambient.json').read_text())
    p204, p226, p211, p223 = (load(r) for r in ('p2_04_chamber_imprint', 'p2_26_imprint_spectrum',
                                                'p2_11_every_way_in', 'p2_23_every_reader'))
    maps_path = RESULTS / 'p2_26_imprint_spectrum' / 'maps.npz'
    params = {'acquisition': 'giza-20250827', 'half_width_m': HALF, 'grid': "the image's own pixels",
              'taper_static_m': TAPER_STATIC, 'taper_dynamic_m': TAPER_DYNAMIC, 'far_field_from_m': FAR_FROM,
              'reliable_dprime': D_RELIABLE, 'genie_snr_db': SNR_GENIE_DB,
              'monte_carlo': {'estimator': MC_REAL, 'synthesizer': MC_SYNTH, 'detection': MC_DETECT},
              'noise_free': True, 'attenuation': 'none', 'scattering': 'fully developed speckle, no reference image'}
    with Run(RID, 'What one image can hold about the chamber, whatever reads it', params) as run:
        checks = {'covariance': check_covariance(g)}
        print(f"  covariance check: {checks['covariance']['relative_difference']:.1e}", flush=True)
        checks['grid'] = check_grid(g)
        checks['continuum'] = check_continuum(g)
        checks['finite_kl'] = check_finite_kl(g)
        checks['detection'] = memo(RID, 'detection', lambda: check_detection(g), __file__, version='checks-v1')
        checks['pulse_domain'] = memo(RID, 'pulse_domain', lambda: check_pulse_domain(g), __file__, version='checks-v1')
        checks['estimator'] = memo(RID, 'estimator', lambda: check_estimator(g), __file__, version='checks-v1')

        kern, kern_src = static_kernels()
        los = np.asarray(g.los_enu)
        rms = np.sqrt(np.mean([(-(kern['Kxx'] * np.sin(p) ** 2 + kern['Kyy'] * np.cos(p) ** 2
                                  + kern['Kxy'] * np.sin(p) * np.cos(p)) @ los) ** 2
                               for p in np.deg2rad(np.arange(0, 180, 7.5))], axis=0))
        pub = np.asarray(p204['imprint_map_los_per_strain_m']).ravel()
        kernel_check = {'source': kern_src, 'peak_m_per_strain': float(rms.max()),
                        'published_peak_m_per_strain': float(pub.max()),
                        'largest_difference_over_peak': float(np.max(np.abs(rms - pub)) / pub.max())}
        static = memo(RID, 'static', lambda: static_cases(g, kern, p204, amb), __file__, version='checks-v1')

        dyn, scene_radius = memo(RID, 'dynamic', lambda: dynamic_information(g, maps_path, p226), __file__,
                                 version='dynamic-v2')   # v2: each direction's own far field
        dynamic = dynamic_cases(dyn, g, amb, maps_path)

        # 4. bright points over the imprint's peak: KL <= SCR <Phi^2> (the point's phase and position handed over, the
        # history unprojected: |exp(-i t) - 1| <= |t|), and the genie: KL <= F_genie / 2
        T = g.aperture_time
        t = np.linspace(-T / 2, T / 2, 2001)
        scr_nat = 10 ** (p223['satellite']['natural']['brightest_db'] / 10)
        scr_cr = 10 ** (p223['manifest']['params']['reflector_scr_db'] / 10)
        k0 = 4 * np.pi / g.lam
        micro = static[0]
        d_peak_micro = micro['imprint_peak_phase_rad'] / k0
        sb = next(d for d in dynamic if d['case'] == 'surface')
        maps226 = np.load(maps_path)
        f226 = maps226['surface_f']
        kt = int(np.argmin(np.abs(f226 - sb['truck']['worst_f_hz'])))
        v_truck = amb['cultural']['bus_or_truck_over_bump']['value']
        H_peak = float(np.abs(maps226['surface_H'][..., kt]).max())
        d_peak_truck = H_peak * np.sqrt(2) * v_truck / (2 * np.pi * f226[kt])
        points = []
        for label, scr in (('the brightest point on open plateau (P2-23)', scr_nat), ('a corner reflector (P2-23)', scr_cr)):
            for case, dpk, fq in (('microseisms', d_peak_micro, micro['f_hz']),
                                  ('truck over a bump, 15 m', d_peak_truck, float(f226[kt]))):
                F = np.mean([inf.fisher_point(scr, k0 * dpk * np.cos(2 * np.pi * fq * t + psi), t)
                             for psi in (0.0, np.pi / 2)])
                kl = scr * (k0 * dpk) ** 2 / 2
                points.append({'point': label, 'scr_db': 10 * np.log10(scr), 'shaking': case, 'f_hz': fq,
                               'peak_displacement_m': dpk, **summary(F, kl)})
        snr = 10 ** (SNR_GENIE_DB / 10)
        genie = [{'case': r['case'], **summary(0.5 * snr * r['fisher_bound_any_frequency'],
                                               0.25 * snr * r['fisher_bound_any_frequency'])} for r in static[:1]]
        xs6, ys6 = maps226['surface_x'], maps226['surface_y']
        tp = taper(xs6, ys6, *TAPER_DYNAMIC)
        Kd = maps226['surface_H'][..., kt] / (2j * np.pi * f226[kt])
        dg = g.dr / np.sin(g.theta)
        axes = pixel_axes(g)
        Fg = 0.5 * snr * max((np.sqrt(inf.fisher_bound(to_pixels(Kd[q] * tp, xs6, ys6, g, axes), (g.dx, dg), g))
                              + np.sqrt(float(dyn['surface']['far_bound'][kt][q]))) ** 2
                             for q in range(len(Kd))) * 2 * v_truck ** 2
        genie.append({'case': 'truck over a bump, 15 m (surface wave)', **summary(Fg, Fg / 2)})

        # 5. the claimed deep structure: the static imprint scales as volume / depth^2 (P2-04), and in the stretch regime
        # the information depends on the peak alone (the sum of |grad v|^2 over the ground does not change with width)
        rel = p204['claimed'][0]['relative_to_bench']
        claim = {'structure': p204['claimed'][0]['structure'], 'relative_imprint': rel,
                 'microseisms': summary(micro['fisher'] * rel ** 2),
                 'noisiest': summary(static[1]['fisher'] * rel ** 2), 'note': 'an estimate by scaling, not computed'}
        pen = p211['A_penetration']['rows'][0]['bands']['X']['two_way_loss_db_if_sand']

        ub = next(d for d in dynamic if d['case'] == 'surface')
        fav = next(d for d in dynamic if d['case'] == 'surface_favourable')
        det = checks['detection']
        pd_short = checks['pulse_domain']['real_rate_short_aperture']
        pd_cycles = checks['pulse_domain']['five_cycles']
        short_dev = max(abs(r['ratio_exact_over_model'] - 1) for r in pd_short if r['f_hz'] >= 5.0)
        slowest = min(pd_short, key=lambda r: r['f_hz'])
        over = [det[i]['pinsker_from_exact_kl'] / det[i]['likelihood_ratio']['largest_pd_minus_pfa'] for i in (0, 1)]
        finding = (
            f"Within its model, one image of ordinary ground holds almost nothing of the chamber, however it is read. The "
            f"chamber reaches the image only through its echo, at least {pen:.0f} dB down at X band (P2-11), or through the "
            f"surface's motion during the pass; every method is a function of the image. For fully developed speckle, "
            f"motion shared by the whole scene holds exactly none, and a rigorous bound on the finite change (the Fisher "
            f"information with an explicit remainder) caps any test's detection rate above its false-alarm rate: "
            f"{static[0]['ceiling']:.1e} under Giza's regional microseism level, {static[1]['ceiling']:.1e} at the "
            f"noisiest stations on Earth, {static[2]['ceiling']:.1e} and {static[3]['ceiling']:.1e} at the measured "
            f"1-3 and 3-8 Hz levels (accuracy above a coin toss: half of each). With P2-26's dynamic imprint and the "
            f"scattered wave counted unattenuated across the whole 5 km scene, the FTA's urban background for the whole "
            f"pass gives {ub['urban']['ceiling']:.1e} ({fav['urban']['ceiling']:.1e} over a room whose roof is 5 m down); "
            f"a truck over a bump 15 m away for the whole pass, at its worst frequency ({ub['truck']['worst_f_hz']:.0f} Hz), "
            f"{ub['truck']['ceiling']:.1e} ({fav['truck']['ceiling']:.1e}). On a small image the exact likelihood-ratio "
            f"test reaches {det[0]['likelihood_ratio']['largest_pd_minus_pfa']:.2f} and "
            f"{det[1]['likelihood_ratio']['largest_pd_minus_pfa']:.2f} where Pinsker's bound from the exact KL allows "
            f"{det[0]['pinsker_from_exact_kl']:.2f} and {det[1]['pinsker_from_exact_kl']:.2f} and weak-signal theory "
            f"predicts {det[0]['weak_signal_best_test']:.2f} and {det[1]['weak_signal_best_test']:.2f}: the ceiling sits "
            f"{over[0]:.2f} and {over[1]:.2f} times above what the best test achieves, a bound, not an attained value. "
            f"The model's Doppler-to-time relation (motion multiplied into the spectrum at each bin's stationary time) "
            f"carries the same information as pulse-by-pulse physics to within {100 * short_dev:.1f}% from 5 to 150 Hz "
            f"on a 0.5 s aperture at the real chirp rate (at {slowest['f_hz']:.1f} Hz, a quarter cycle across it, the "
            f"model holds {1 / slowest['ratio_exact_over_model']:.2f} times more); with five cycles of motion across the "
            f"aperture the model's excess falls from {100 * (1 / pd_cycles[0]['ratio_exact_over_model'] - 1):.1f}% to "
            f"{100 * (1 / pd_cycles[-1]['ratio_exact_over_model'] - 1):.2f}% as the time-bandwidth product rises from "
            f"{pd_cycles[0]['time_bandwidth']:.0e} to {pd_cycles[-1]['time_bandwidth']:.0e} (the real dwell's is "
            f"{checks['pulse_domain']['real_time_bandwidth']:.0e}): an edge effect that fades, on the side of more "
            f"information. The local extrapolation factor under the "
            f"microseisms, {static[0]['local_extrapolation_factor']:.0e}, is a scale within the model, not a physical "
            f"shaking requirement.")
        run.save({'checks': checks, 'kernels': kernel_check, 'static': static, 'dynamic': dynamic,
                  'scene_radius_m': scene_radius, 'points': points, 'genie': genie, 'claimed': claim,
                  'penetration_db_x_band': pen, 'cells_per_m2': inf.cells_per_m2(g),
                  'stretch_seconds': g.R0 / g.V_platform, 'paired_echo_m_per_hz': g.V / g.Ka, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
