"""What one image can hold about a chamber, whatever reads it (Phase 2, P2-25, P2-27).

The chamber can reach a radar image only through what the image is made of: the echo from what the radar wave
reaches (the surface and its first centimetres; P2-11 puts the bench chamber's roof at least 347 dB down at X band),
the surface itself (which any photograph shows), and the surface's motion during the pass. Every method is a function
of the image, so none can hold more about the chamber than the image does (the data-processing inequality). This
module bounds what the motion leaves in the image, for any method, within a stated measurement model.

The model. For an ideal focusing processor a scatterer with reflectivity s at x, whose slant range grows by d(x, t)
during the pass, adds s exp(-2 pi i f.x) exp(-i phi) to the image's spectrum at f = (f_a, f_r), with phi = 4 pi d /
lambda taken at t = f_a / Ka, the slow time at which the pass sees that Doppler bin (sarsim.synth; the relation is the
stationary-phase one, checked against pulse-by-pulse echoes in P2-25). Ordinary ground is fully developed speckle: s
white and circular Gaussian on the image's own pixel grid, so the image is a zero-mean circular Gaussian vector and its
covariance C holds everything it says about the ground. No receiver noise (which could only lower what follows), no
reference image of the same ground (a genie that supplies one is bounded separately).

The bound, finite and rigorous within the model. Let the imprint be phi = a Phi(x, t), a = 1 with the chamber and 0
without. Whiten the covariance change: E = C0^-1/2 (C_a - C0) C0^-1/2. Then, exactly, KL(P_a || P_0) = sum over the
eigenvalues l of I + E of (l - 1 - ln l), and since l - 1 - ln l <= (l - 1)^2 / (2 min(l, 1)) and every |l - 1| is
at most |E|_F,

    KL <= |E|_F^2 / (2 (1 - |E|_F))         while |E|_F < 1.

E is a Fourier transform of exp(-i a dPhi) - 1 (dPhi the imprint's phase difference between two Doppler bins); its
first-order part has Frobenius norm a sqrt(F), F the Fisher information below. The rest is a power series in a dPhi,
|dPhi| <= 2 |Phi|, whose n-th term has Frobenius norm at most (2a)^n / n! sqrt((N_b / N) sum_x |Phi(x)|^(2n)) by
Parseval (N_b / N the share of the image's frequency bins inside the processed band); summed from n = 2 it is eps
(`remainder`). So |E|_F <= a sqrt(F) + eps, no approximation. Any test, however built, has

    detection rate - false-alarm rate <= TV <= min( sqrt(KL / 2), sqrt(1 - exp(-KL)) )

(Pinsker; Bretagnolle-Huber); for two equally likely cases the best accuracy is 1/2 + TV/2. For a weak signal the
log-likelihood ratio is Gaussian with deflection d = a sqrt(F), the best test reaches TV = 2 Phi(d/2) - 1 = 0.40 d, and
the bound gives 0.50 d: tight to a factor 1.25 there (checked with the exact likelihood-ratio test on a small image).

The Fisher information. For a zero-mean circular Gaussian F = tr(C^-1 C' C^-1 C'), and with the bins independent,

    F = sum_{f, f' in band} |Phi~(f - f'; f) - Phi~(f - f'; f')|^2,

Phi~(k; f) the discrete spatial Fourier coefficient (over the image's pixels) of Phi(., t(f)) at k.

1. Motion shared by the whole scene holds no information at all: Phi~ is then nonzero only at k = 0, where the
   difference vanishes. Not little: none.
2. For any motion F <= 4 N <Phi^2>, N the independent resolution cells and <Phi^2> the imprint's mean-square phase
   over them and over the pass.
3. For slow, smooth motion F = (R / V_s)^2 sum_cells <(dv/dx)^2>: the image is locally stretched along track by the
   line-of-sight velocity's gradient times R / V_s (79 s at Giza), and only that stretch survives the speckle.

For a sinusoid, d = Re[K(x) exp(i(2 pi f_m (t + x_a / V_g) + psi))] (x_a / V_g: each row's own zero-Doppler time), F
depends on psi as F(psi) = A + Re(B exp(2 i psi)): `fisher_grid` computes A (the mean over psi) and |B| exactly on the
image's pixel grid; `fisher_sinusoid` is the continuum approximation of A, for quick use. Handing the method the
ground's reflectivity (a genie) raises the information to at most 2 sum_cells SNR <Phi^2>, whatever the scene; a lone
bright point with signal-to-clutter ratio SCR holds F = 2 SCR <Phi_perp^2> about its own motion (Phi_perp: the phase
history less a constant and a linear trend, its unknown phase and position).

Arrays of the imprint are [azimuth, range] on a regular grid (metres), azimuth along the radar's track: the image's own
pixels (dx, dr in slant range) for the exact functions, ground range for the continuum ones.
"""
from __future__ import annotations

import math

import numpy as np


def cells_per_m2(g):
    """Independent resolution cells per square metre of ground: the processed bandwidths, azimuth times ground range."""
    return g.nu_band * g.kr_band * np.sin(g.theta)


def kl_upper(F, eps=0.0):
    """A rigorous upper bound on KL(P_a || P_0) from the Fisher information F (at a = 1) and the remainder bound eps:
    x^2 / (2 (1 - x)), x = sqrt F + eps; infinite where x >= 1, where the bound says nothing."""
    x = np.sqrt(np.asarray(F, float)) + np.asarray(eps, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x < 1, x ** 2 / (2 * (1 - x)), np.inf)


def tv_upper(kl):
    """The most any test's detection rate can exceed its false-alarm rate, given KL: min of Pinsker's sqrt(KL / 2) and
    Bretagnolle-Huber's sqrt(1 - exp(-KL)), at most 1."""
    kl = np.asarray(kl, float)
    with np.errstate(over='ignore', invalid='ignore'):
        return np.minimum(np.minimum(np.sqrt(kl / 2), np.sqrt(-np.expm1(-kl))), 1.0)


def ceiling(F, eps=0.0):
    """detection rate - false-alarm rate <= this, for any test, within the model: tv_upper(kl_upper(F, eps)). For two
    equally likely cases the best accuracy is 1/2 plus half of it."""
    return tv_upper(kl_upper(F, eps))


def best_test_tv(F):
    """What the best test reaches for a weak signal (the log-likelihood ratio Gaussian with deflection sqrt F):
    2 Phi(sqrt(F) / 2) - 1. Asymptotic, for comparison with the ceiling."""
    from scipy.special import erf
    return erf(np.sqrt(np.asarray(F, float)) / (2 * np.sqrt(2)))


def remainder(K, g, shape=None, a=1.0, terms=24):
    """eps, a rigorous bound on the Frobenius norm of the whitened covariance change beyond first order, for
    d = a Re[K exp(i ...)] with K [n_a, n_r] metres on the image's pixels. exp(-i a dPhi) - 1 + i a dPhi is a power
    series in a dPhi, dPhi = Re[K+ dtau+] with |dtau| <= 2; term n has Frobenius norm at most
    (2a)^n / n! sqrt((N_b / N) sum_x |4 pi K / lambda|^(2n)) (Parseval over the image, N_b in-band bins). Summed from
    n = 2; the tail beyond `terms` is negligible for the phases here and is bounded by the last term's ratio test."""
    shape = np.shape(K) if shape is None else shape
    in_a, in_r, _ = band_masks(g, shape)
    share = in_a.mean() * in_r.mean()
    k = 4 * np.pi / g.lam * np.abs(np.asarray(K, complex)).ravel()
    total, last = 0.0, 0.0
    for n in range(2, terms + 1):
        last = (2 * a) ** n / math.factorial(n) * np.sqrt(share * np.sum(k ** (2 * n)))
        total += last
    ratio = 2 * a * k.max() / (terms + 1)
    return float(total + (last * ratio / (1 - ratio) if ratio < 1 else np.inf))


def _spectrum(K, spacing, pad):
    """Continuous Fourier transform of a gridded field, zero-padded to `pad` samples, with its frequencies (cycles/m)."""
    da, dg = spacing
    S = np.fft.fft2(K, s=pad) * (da * dg)
    return S, np.fft.fftfreq(pad[0], da), np.fft.fftfreq(pad[1], dg)


def fisher_sinusoid(K, spacing, f_m, g, pad=None, swing=False):
    """Expected Fisher information (over the phase psi) one image of speckle holds about a slant-range motion
    d(x, t) = K(x) cos(2 pi f_m (t + x_a / V_g) + psi), K [n_a, n_g] in metres on a grid with `spacing` (m).
    Noise-free: the most the image can hold. K may be complex, d = Re[K exp(i(2 pi f_m (t + x_a / V_g) + psi))],
    for a travelling or scattered wave. `pad` must leave room for the paired echoes, which the motion displaces by
    f_m V_g / |Ka| metres along track (1.24 m per hertz at Giza): the default pads by four times that.
    With `swing`, also an upper bound on |B| (F(psi) = A + Re(B exp(2 i psi))): the continuum of fisher_grid's
    geometric sum, |sum over the bins' overlap of exp(4 pi i f_m t)| <= its length times |sinc(2 f_m T_overlap)|, and
    2 |S+ S-| <= |S+|^2 + |S-|^2; returns (A, B_upper)."""
    K = np.asarray(K)
    da, dg = spacing
    if pad is None:
        need = 4 * f_m * g.V / g.Ka
        pad = (int(2 ** np.ceil(np.log2(K.shape[0] + need / da))), int(2 ** np.ceil(np.log2(2 * K.shape[1]))))
    xa = (np.arange(K.shape[0]) - K.shape[0] / 2) * da
    k0 = 4 * np.pi / g.lam
    total, total_b = 0.0, 0.0
    T_p = g.nu_band * g.V / g.Ka                              # the processed span of slow time
    for sgn in (1, -1):
        Kc = k0 * (K if sgn > 0 else np.conj(K)) * np.exp(sgn * 2j * np.pi * f_m * xa / g.V)[:, None]
        S, ka, kg = _spectrum(Kc, spacing, pad)
        rho_a = np.clip(1 - np.abs(ka) / g.nu_band, 0, None)
        kg_band = g.kr_band * np.sin(g.theta)
        rho_g = np.clip(1 - np.abs(kg) / kg_band, 0, None)
        wa = np.sin(np.pi * f_m * g.V * ka / g.Ka) ** 2 * rho_a
        dk = (ka[1] - ka[0]) * (kg[1] - kg[0])
        P = np.abs(S) ** 2 * rho_g[None, :]
        total += float(np.sum(P * wa[:, None]) * dk)
        if swing:
            total_b += float(np.sum(P * (wa * np.abs(np.sinc(2 * f_m * T_p * rho_a)))[:, None]) * dk)
    if swing:
        return cells_per_m2(g) * total, cells_per_m2(g) * total_b
    return cells_per_m2(g) * total


def fisher_grid(K, g, f_m, pad=None):
    """Exact F(psi) = A + Re(B exp(2 i psi)) for d = Re[K exp(i(2 pi f_m (t + x_a / V_g) + psi))] with scatterers on the
    image's own pixel grid (noise-free speckle): K complex [n_a, n_r] metres at the image's spacing (g.dx along track,
    g.dr in slant range), zero-padded to `pad` (default: room for the paired echoes, f_m V_g / |Ka| metres each way).
    Returns (A, |B|): the mean over psi, and the amplitude of its swing (the largest F over psi is A + |B|)."""
    K = np.asarray(K, complex)
    na, nr = K.shape
    if pad is None:
        need = 4 * f_m * g.V / g.Ka
        pad = (int(2 ** np.ceil(np.log2(na + need / g.dx + 8))), int(2 ** np.ceil(np.log2(nr + 8))))
    Na, Nr = pad
    N = Na * Nr
    k0 = 4 * np.pi / g.lam
    xa = (np.arange(na) - na / 2) * g.dx
    ph = np.exp(2j * np.pi * f_m * xa / g.V)[:, None]
    Kp = np.zeros((Na, Nr), complex)
    Km = np.zeros((Na, Nr), complex)
    Kp[:na, :nr] = k0 * K * ph
    Km[:na, :nr] = k0 * np.conj(K) * np.conj(ph)
    Fp, Fm = np.fft.fft2(Kp) / N, np.fft.fft2(Km) / N
    in_a, in_r, _ = band_masks(g, pad)
    ka = np.rint(np.fft.fftfreq(Na) * Na).astype(int)[in_a]
    kr = np.rint(np.fft.fftfreq(Nr) * Nr).astype(int)[in_r]
    lo_a, hi_a, lo_r, hi_r = ka.min(), ka.max(), kr.min(), kr.max()
    Ma, Mr = hi_a - lo_a + 1, hi_r - lo_r + 1                 # contiguous bands
    c = 2 * np.pi * f_m * g.V / (Na * g.dx * g.Ka_signed)     # phase of the motion per azimuth bin
    da = np.arange(-(Ma - 1), Ma)
    dr = np.arange(-(Mr - 1), Mr)
    s2 = np.sin(c * da / 2) ** 2
    Ra, Rr = Ma - np.abs(da), Mr - np.abs(dr)
    # G(da): sum over bins k with k and k - da both in band of exp(i c (2k - da))
    k_start = np.maximum(lo_a, lo_a + da)
    count = np.minimum(hi_a, hi_a + da) - k_start + 1
    q = np.exp(2j * c)
    with np.errstate(divide='ignore', invalid='ignore'):
        geo = np.where(np.abs(1 - q) > 1e-12, (1 - q ** count) / (1 - q), count)
    G = np.exp(-1j * c * da) * np.exp(2j * c * k_start) * geo
    Wa = np.zeros(Na)
    WBa = np.zeros(Na, complex)
    np.add.at(Wa, da % Na, Ra * s2)
    np.add.at(WBa, da % Na, s2 * G)
    Wr = np.zeros(Nr)
    np.add.at(Wr, dr % Nr, Rr)
    A = float(np.sum((np.abs(Fp) ** 2 + np.abs(Fm) ** 2) * Wa[:, None] * Wr[None, :]))
    B = -2 * np.sum(Fp * np.conj(Fm) * WBa[:, None] * Wr[None, :])
    return A, float(np.abs(B))


def fisher_bound(K, spacing, g):
    """F <= 4 N <Phi^2> for the sinusoid above: 2 cells_per_m2 (4 pi / lambda)^2 integral |K|^2 dA. Any frequency."""
    da, dg = spacing
    return 2 * cells_per_m2(g) * (4 * np.pi / g.lam) ** 2 * float(np.sum(np.abs(np.asarray(K)) ** 2)) * da * dg


def fisher_slow(V, spacing, g):
    """The slow limit, F = (R / V_s)^2 sum_cells <(dv/dx_a)^2>, for a sinusoidal line-of-sight velocity of amplitude
    V [n_a, n_g] (m/s): the image's local stretch along track, squared and summed over its cells."""
    da, dg = spacing
    dv = np.gradient(np.asarray(V, float), da, axis=0)
    return (g.R0 / g.V_platform) ** 2 * cells_per_m2(g) * 0.5 * float(np.sum(dv ** 2)) * da * dg


def fisher_genie(K, spacing, g, snr):
    """Fisher information when the method is handed the ground's reflectivity, for the same sinusoid: 2 SNR
    sum_cells <Phi^2>, SNR the pixels' signal-to-noise ratio (power). Holds whatever the scene."""
    return 0.5 * snr * fisher_bound(K, spacing, g)


def fisher_point(scr, phase, t):
    """Fisher information a lone bright point holds about the amplitude of its own phase history `phase` (radians,
    sampled at slow times t across the pass): 2 SCR <phase_perp^2>, the history's part that neither its unknown phase
    (a constant) nor its unknown position (a trend in t) can absorb."""
    t = np.asarray(t, float)
    X = np.stack([np.ones_like(t), t], 1)
    r = phase - X @ np.linalg.lstsq(X, phase, rcond=None)[0]
    return 2 * scr * float(np.mean(r ** 2))


# ---------------------------------------------------------------- exact, on an image's own pixel grid (checks)

def band_masks(g, shape):
    """In-band masks of an image's azimuth and range frequency bins, in numpy's fft order (as sarsim.synth)."""
    nu = np.fft.fftfreq(shape[0], d=g.dx)
    kr = np.fft.fftfreq(shape[1], d=g.dr)
    return np.abs(nu) <= g.nu_band / 2, np.abs(kr) <= g.kr_band / 2, nu


def fisher_exact(components, g, shape, w=None):
    """The exact double sum over an image's own frequency bins, for scatterers on its pixel grid (white speckle):
    `components` is a list of (Kj [n_a, n_r] phase per unit a, tau_j(t)) with Phi(x, t) = sum_j Kj(x) tau_j(t); t is
    slow time relative to each scatterer's zero-Doppler crossing (the synthesizer's). `w` optional per-bin weights
    [n_a, n_r] (default: 1 in band, noise-free). O(n^2): for small images."""
    in_a, in_r, nu = band_masks(g, shape)
    if w is None:
        w = (in_a[:, None] & in_r[None, :]).astype(float)
    t = g.nu_to_time(nu)
    N = shape[0] * shape[1]
    Kt = [np.fft.fft2(K) / N for K, _ in components]          # Phi~_j(k), k in fft order
    taus = [np.asarray(tau(t)) for _, tau in components]
    ia, ir = np.nonzero(w)
    F = 0.0
    for a, r in zip(ia, ir):
        da = (a - np.arange(shape[0])) % shape[0]             # f - f' for every f' (azimuth)
        dr = (r - np.arange(shape[1])) % shape[1]
        diff = sum(Kt_j[np.ix_(da, dr)] * (tau_j[a] - tau_j[:, None]) for Kt_j, tau_j in zip(Kt, taus))
        F += w[a, r] * float(np.sum(w * np.abs(diff) ** 2))
    return F


def band_operator(components, g, shape, a=1.0):
    """B(a) [n_band, N]: the in-band image spectrum (unitary DFT) of scatterers on an image's pixel grid, moved by
    a Phi with Phi(x, t) = sum_j Kj(x) tau_j(t) (radians). For unit white speckle the image's in-band covariance is
    C(a) = B B^H, and B(0) has orthonormal rows, so C(0) = I. Exact, for small images."""
    in_a, in_r, nu = band_masks(g, shape)
    na, nr = shape
    t = g.nu_to_time(nu)
    ia, ir = np.nonzero(in_a[:, None] & in_r[None, :])
    xa, xr = np.meshgrid(np.arange(na), np.arange(nr), indexing='ij')
    E = np.exp(-2j * np.pi * (np.outer(np.fft.fftfreq(na)[ia], xa.ravel()) + np.outer(np.fft.fftfreq(nr)[ir], xr.ravel())))
    E /= np.sqrt(na * nr)
    if a == 0.0 or not components:
        return E
    Phi = sum(np.outer(np.asarray(tau(t))[ia], np.ravel(Kj)) for Kj, tau in components)
    return E * np.exp(-1j * a * Phi)


def kl_exact(C1):
    """KL(CN(0, C1) || CN(0, I)) = sum over eigenvalues l of C1 of (l - 1 - ln l): the finite change, exactly."""
    lam = np.linalg.eigvalsh(C1)
    return float(np.sum(lam - 1 - np.log(lam)))


def score(img, K, tau, g):
    """The locally most powerful statistic for a small real imprint a K(x) tau(t) in one noise-free speckle image:
    T = (2 / s2) Im sum_x K(x) z(x) conj(z_tau(x)), z_tau the image re-formed with each Doppler bin weighted by
    tau(t). E[T] = a F + O(a^2) and var T = F without the imprint, so T / F estimates a at the Cramer-Rao bound. For
    motion shared by the scene (K constant) T is identically zero."""
    Y = np.fft.fft2(img)
    in_a, in_r, nu = band_masks(g, img.shape)
    band = in_a[:, None] & in_r[None, :]
    s2 = float(np.mean(np.abs(Y[band]) ** 2)) / img.size
    zt = np.fft.ifft2(np.asarray(tau(g.nu_to_time(nu)))[:, None] * Y * band)
    z = np.fft.ifft2(Y * band)
    return 2.0 / s2 * float(np.imag(np.sum(K * z * np.conj(zt))))


def speckle_image(g, shape, rng, K=None, tau=None, a=0.0, terms=12):
    """A noise-free speckle image with scatterers on its pixel grid, moved by a K(x) tau(t) (phase, radians): the
    model the Fisher information is computed for. exp(-i a K tau) is expanded in powers of a K (exact to rounding for
    |a K tau| < 0.5), each power one FFT."""
    s = (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)
    in_a, in_r, nu = band_masks(g, shape)
    band = in_a[:, None] & in_r[None, :]
    if K is None or a == 0.0:
        return np.fft.ifft2(np.fft.fft2(s) * band)
    if np.max(np.abs(a * K)) > 0.5:
        raise ValueError('the series needs |a K| < 0.5')
    tt = np.asarray(tau(g.nu_to_time(nu)))[:, None]
    Y = np.zeros(shape, complex)
    term = s.astype(complex)
    for n in range(terms):
        Y += (-1j * a * tt) ** n / math.factorial(n) * np.fft.fft2(term)
        term = term * K
    return np.fft.ifft2(Y * band)
