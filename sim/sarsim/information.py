"""What one image can hold about a chamber, whatever reads it (Phase 2, P2-25).

The chamber can reach a radar image only through what the image is made of: the echo from what the radar wave
reaches (the surface and its first centimetres; P2-11 puts the bench chamber's roof at least 347 dB down at X band),
the surface itself (which any photograph shows), and the surface's motion during the pass. Every method is a function
of the image, so none can hold more about the chamber than the image does (the data-processing inequality). This
module bounds what the motion leaves in the image, for any method.

The image. For an ideal focusing processor a scatterer with reflectivity s at x, whose slant range grows by d(x, t)
during the pass, adds s exp(-2 pi i f.x) exp(-i phi) to the image's spectrum at f = (f_a, f_r), with phi = 4 pi d /
lambda taken at t = f_a / Ka, the slow time at which the pass sees that Doppler bin (sarsim.synth). Ordinary ground is
fully developed speckle: s is white and circular Gaussian, so the image is a zero-mean circular Gaussian vector and
its covariance C holds everything it says about the ground.

The bound. Let the chamber's imprint be phi = a Phi(x, t): a = 1 with the chamber, a = 0 without it. Any test, however
built, has P(detect) - P(false alarm) <= TV <= sqrt(KL / 2) (Pinsker); for a small imprint KL = F a^2 / 2 with F the
Fisher information, so the edge any method can have over a guess is at most sqrt(F) / 2. For a zero-mean circular
Gaussian F = tr(C^-1 C' C^-1 C'), and with speckle white (C diagonal in frequency, band window W, receiver noise),

    F = sum_{f, f'} w(f) w(f') |Phi~(f - f'; f) - Phi~(f - f'; f')|^2,      w = |W|^2 s2 / (|W|^2 s2 + noise),

Phi~(k; f) the spatial Fourier coefficient of Phi(., t(f)) at k. Without noise w is 1 in the band, whatever W is.

1. Motion shared by the whole scene holds no information at all: Phi~ is then nonzero only at k = 0, where the
   difference vanishes. Not little: none.
2. For any motion F <= 4 N <Phi^2>, N the independent resolution cells and <Phi^2> the imprint's mean-square phase
   over them and over the pass.
3. For slow, smooth motion F = (R / V_s)^2 sum_cells <(dv/dx)^2>: the image is locally stretched along track by the
   line-of-sight velocity's gradient times R / V_s (79 s at Giza), and only that stretch survives the speckle.

For a sinusoid, d = K(x) cos(2 pi f_m (t + x_a / V_g) + psi) (x_a / V_g: each row's own zero-Doppler time), the
expectation over psi is, in the continuum,

    F = nu_band kg_band  integral (|K+^(k)|^2 + |K-^(k)|^2) sin^2(pi f_m V_g k_a / Ka) rho_a(k_a) rho_g(k_g) dk,

K+- = (4 pi / lambda) K exp(+-2 pi i f_m x_a / V_g), ^ the continuous Fourier transform, rho the band's normalised
autocorrelation. Handing the method the ground's reflectivity (a genie; no single image has it) only raises the
information, to F = 2 sum_cells SNR <Phi^2> whatever the scene; a lone bright point with signal-to-clutter ratio SCR
holds F = 2 SCR <Phi_perp^2> about its own motion (Phi_perp: the phase history less a constant and a linear trend,
its unknown phase and position).

Arrays of the imprint are [azimuth, ground range] on a regular grid (metres), azimuth along the radar's track.
"""
from __future__ import annotations

import math

import numpy as np


def cells_per_m2(g):
    """Independent resolution cells per square metre of ground: the processed bandwidths, azimuth times ground range."""
    return g.nu_band * g.kr_band * np.sin(g.theta)


def edge(F):
    """The most any method's detection rate can exceed its false-alarm rate: sqrt(KL / 2) with KL = F / 2."""
    return 0.5 * np.sqrt(np.asarray(F, float))


def _spectrum(K, spacing, pad):
    """Continuous Fourier transform of a gridded field, zero-padded to `pad` samples, with its frequencies (cycles/m)."""
    da, dg = spacing
    S = np.fft.fft2(K, s=pad) * (da * dg)
    return S, np.fft.fftfreq(pad[0], da), np.fft.fftfreq(pad[1], dg)


def fisher_sinusoid(K, spacing, f_m, g, pad=None):
    """Expected Fisher information (over the phase psi) one image of speckle holds about a slant-range motion
    d(x, t) = K(x) cos(2 pi f_m (t + x_a / V_g) + psi), K [n_a, n_g] in metres on a grid with `spacing` (m).
    Noise-free: the most the image can hold. K may be complex, d = Re[K exp(i(2 pi f_m (t + x_a / V_g) + psi))],
    for a travelling or scattered wave. `pad` must leave room for the paired echoes, which the motion displaces by
    f_m V_g / |Ka| metres along track (1.24 m per hertz at Giza): the default pads by four times that."""
    K = np.asarray(K)
    da, dg = spacing
    if pad is None:
        need = 4 * f_m * g.V / g.Ka
        pad = (int(2 ** np.ceil(np.log2(K.shape[0] + need / da))), int(2 ** np.ceil(np.log2(2 * K.shape[1]))))
    xa = (np.arange(K.shape[0]) - K.shape[0] / 2) * da
    k0 = 4 * np.pi / g.lam
    total = 0.0
    for sgn in (1, -1):
        Kc = k0 * (K if sgn > 0 else np.conj(K)) * np.exp(sgn * 2j * np.pi * f_m * xa / g.V)[:, None]
        S, ka, kg = _spectrum(Kc, spacing, pad)
        rho_a = np.clip(1 - np.abs(ka) / g.nu_band, 0, None)
        kg_band = g.kr_band * np.sin(g.theta)
        rho_g = np.clip(1 - np.abs(kg) / kg_band, 0, None)
        wa = np.sin(np.pi * f_m * g.V * ka / g.Ka) ** 2 * rho_a
        dk = (ka[1] - ka[0]) * (kg[1] - kg[0])
        total += float(np.sum(np.abs(S) ** 2 * wa[:, None] * rho_g[None, :]) * dk)
    return cells_per_m2(g) * total


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
