"""What one dwell can see of the ground's motion (Phase 2, P2-03).

A look is the image formed from part of the processed Doppler band. In a dwell the band is slow time
(t = f / Ka), so a look W seconds long sees the scene as it was during those W seconds, and any motion is
averaged over them: its gain at frequency f is sinc(f W).

Three facts follow, each exact.

1. Two looks cut from one image share the spectrum they overlap, sample for sample. Their complex coherence is
   sum_overlap |Y|^2 / sqrt(sum_A |Y|^2 sum_B |Y|^2), a function of the image's power spectrum alone: 1 - delta/W
   for a flat spectrum and rectangular looks, whatever the ground does.
2. For the same reason the cross-spectrum of two such looks, Y H_A (Y H_B)^* = |Y|^2 H_A H_B, is real and
   non-negative, so their complex cross-correlation peaks at zero lag whatever the ground does. A tracker that
   correlates complex looks of one image is blind to motion by construction.
3. Motion shared by a stretch of ordinary ground (fully developed speckle) leaves the image's statistics
   unchanged: the Doppler bins are independent circular Gaussians and multiplying them by unit phasors changes
   nothing. What motion does leave is the position of each look's image of an isolated bright scatterer (its
   envelope moves), and the apparent brightness that a displacement gradient produces.

Frequencies here are hertz on the image's own azimuth sampling (the processing PRF), band centred at zero.
"""
from __future__ import annotations

import numpy as np


def frequencies(g, n):
    """Azimuth frequency (Hz) of each row of an n-row image's fftshifted spectrum."""
    return np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / g.prf))


def look_masks(g, n, centres_s, width_s):
    """Rectangular look passbands [K, n] on the fftshifted grid: look k spans width_s seconds of slow time centred
    at centres_s[k] (seconds from the zero-Doppler crossing, the product's convention f = Ka t)."""
    f = frequencies(g, n)[None, :]
    c = g.Ka_signed * np.asarray(centres_s, float)[:, None]
    half = g.Ka * width_s / 2
    m = ((f >= c - half) & (f < c + half)).astype(np.float64)
    inside = np.abs(frequencies(g, n)) <= g.doppler_band_hz / 2
    m *= inside[None, :]
    if not m.any(axis=1).all():
        raise ValueError('a look lies outside the processed band')
    return m


def spectrum(img):
    return np.fft.fftshift(np.fft.fft(img, axis=0), axes=0)


def looks(img, masks):
    S = spectrum(img)
    return np.stack([np.fft.ifft(np.fft.ifftshift(S * m[:, None], axes=0), axis=0) for m in masks])


def coherence(a, b):
    """Complex coherence of two images over their whole extent."""
    return abs(np.vdot(b, a)) / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real)


def coherence_from_power(img, mask_a, mask_b):
    """The same coherence from the image's azimuth power spectrum alone (Parseval): fact 1."""
    P = (np.abs(spectrum(img)) ** 2).sum(axis=1)
    return (P * mask_a * mask_b).sum() / np.sqrt((P * mask_a).sum() * (P * mask_b).sum())


def region_shift(a, b, rows, cols, envelope=False):
    """Azimuth shift (pixels) of b relative to a over a region. The cross-spectrum is summed over the region's range
    columns; the integer peak of the cross-correlation gives a coarse shift, and the slope of the remaining
    cross-spectral phase, weighted by its magnitude, refines it (spectral diversity). Picking the interpolated peak
    instead is biased toward zero by the window's edges, by about sigma^2 / n pixels for a look whose correlation
    width is sigma, which matters for narrow looks. `envelope` correlates magnitudes (mean removed)."""
    A, B = a[rows][:, cols], b[rows][:, cols]
    if envelope:
        A = np.abs(A) - np.abs(A).mean()
        B = np.abs(B) - np.abs(B).mean()
    n = A.shape[0]
    X = (np.fft.fft(B, axis=0) * np.conj(np.fft.fft(A, axis=0))).sum(axis=1)
    cc = np.abs(np.fft.ifft(X))
    coarse = int(np.argmax(cc))
    coarse = coarse - n if coarse > n / 2 else coarse
    k = np.fft.fftfreq(n)
    Xr = X * np.exp(2j * np.pi * k * coarse)
    w = np.abs(Xr)
    keep = w > 0.05 * w.max()
    if keep.sum() < 3:
        return float(coarse)
    ph = np.angle(Xr[keep] * np.conj(Xr[keep][np.argmax(w[keep])]))
    slope = np.polyfit(k[keep], ph, 1, w=w[keep])[0]
    return coarse - slope / (2 * np.pi)


def velocity_series(L, g, rows, cols, envelope=False, reference=None):
    """Line-of-sight velocity (m/s) of a region in each look relative to the reference look, as a single dwell
    allows: every look is registered against the reference look of the same image."""
    c = len(L) // 2 if reference is None else reference
    return np.array([g.shift_to_velocity(region_shift(L[c], L[k], rows, cols, envelope) * g.dx)
                     for k in range(len(L))])


def paired_velocity(L_still, L_moved, g, rows, cols):
    """Velocity of a region in each look, read by comparing each look with the same look of the motionless scene:
    possible only in simulation, it shows what motion the data hold."""
    return np.array([g.shift_to_velocity(region_shift(a, b, rows, cols) * g.dx)
                     for a, b in zip(L_still, L_moved)])


def inject_region_motion(slc, g, weight, displacement):
    """Move part of a real image along the line of sight: rows are azimuth, `weight` [n_rows] (0..1) selects the
    part that moves, `displacement(t)` is its slant-range increase (m) at slow time t. The selected content's
    azimuth spectrum is multiplied by exp(-4 pi j d(t) / lambda) with t = f / Ka, the rest is left as it was."""
    w = np.asarray(weight, float)[:, None]
    f = frequencies(g, slc.shape[0])
    t = f / g.Ka_signed
    phase = np.exp(-4j * np.pi * displacement(t) / g.lam)[:, None]
    moved = np.fft.ifft(np.fft.ifftshift(spectrum(slc * w) * phase, axes=0), axis=0)
    return (slc * (1 - w) + moved).astype(np.complex64)


def look_gain(f_hz, width_s):
    """A look averages velocity over its duration: gain sinc(f W) at frequency f."""
    return np.sinc(np.asarray(f_hz, float) * width_s)

def motionless_twin(img, rng, smooth=(9, 5)):
    """A copy of an image in which nothing moves and nothing is inside: the image's own brightness pattern (smoothed over
    about 3 x 3 resolution cells) and mean spectral envelopes in both axes, with fresh complex Gaussian speckle."""
    from scipy.ndimage import uniform_filter
    img = np.asarray(img)
    I = uniform_filter(np.abs(img) ** 2, smooth)
    z = (rng.standard_normal(img.shape) + 1j * rng.standard_normal(img.shape)).astype(np.complex64) / np.sqrt(2)
    t = np.sqrt(I).astype(np.float32) * z
    Fa = np.sqrt(np.mean(np.abs(np.fft.fft(img, axis=0)) ** 2, axis=1) / np.mean(np.abs(np.fft.fft(t, axis=0)) ** 2, axis=1))
    t = np.fft.ifft(np.fft.fft(t, axis=0) * Fa[:, None], axis=0)
    Fr = np.sqrt(np.mean(np.abs(np.fft.fft(img, axis=1)) ** 2, axis=0) / np.mean(np.abs(np.fft.fft(t, axis=1)) ** 2, axis=0))
    t = np.fft.ifft(np.fft.fft(t, axis=1) * Fr[None, :], axis=1)
    return t.astype(np.complex64)
