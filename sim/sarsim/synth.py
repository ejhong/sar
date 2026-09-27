"""Spectral-domain synthesis of a focused SLC image from point scatterers.

For an ideal focusing processor a stationary point scatterer at (x_i, r_i) contributes
    S(k_r, nu) = w_i g_i(nu) exp(-j 2 pi (k_r r_i + nu x_i))
inside the processed (rectangular) band. Aspect-dependent gain g_i(nu) and line-of-sight
vibration (a phase modulation along nu, because nu maps to slow time) are applied per bin.
The image is the inverse 2-D FFT. Homogeneous ground clutter is complex white noise on the
in-band bins, which is exactly a fully developed speckle field with the processor's spectral
support. Array convention: SLC[azimuth, range].
"""
import numpy as np


def pixel_axes(geom, shape, x_centre=0.0, r_centre=0.0):
    Nx, Nr = shape
    x = x_centre + (np.arange(Nx) - Nx // 2) * geom.dx
    r = r_centre + (np.arange(Nr) - Nr // 2) * geom.dr
    return x, r


def image_coords(scat, geom, x_centre=0.0, r_centre=0.0):
    """Azimuth and slant-range image coordinates (m, relative to the grid centre)."""
    r = scat.y * np.sin(geom.theta) - scat.z * np.cos(geom.theta) - r_centre
    x = scat.x - x_centre
    return x, r


def synthesize(scat, geom, shape, x_centre=0.0, r_centre=0.0, ground_intensity=0.0,
               seed=0, chunk=4096, verbose=False):
    Nx, Nr = shape
    nu = np.fft.fftfreq(Nx, d=geom.dx)
    kr = np.fft.fftfreq(Nr, d=geom.dr)
    in_x = np.abs(nu) <= geom.nu_band / 2
    in_r = np.abs(kr) <= geom.kr_band / 2
    nu_in = nu[in_x]
    kr_in = kr[in_r]
    S = np.zeros((in_x.sum(), in_r.sum()), np.complex64)
    if scat is not None and scat.n > 0:
        x_img, r_img = image_coords(scat, geom, x_centre, r_centre)
        x_img = x_img + (Nx // 2) * geom.dx
        r_img = r_img + (Nr // 2) * geom.dr
        w = scat.amp / geom.band_frac ** 2 * np.exp(1j * scat.phase)
        t = geom.nu_to_time(nu_in)[:, None]
        for i in range(0, scat.n, chunk):
            sl = slice(i, i + chunk)
            Ex = np.exp(-2j * np.pi * np.outer(nu_in, x_img[sl]))
            g = scat.iso[sl] + scat.flash[sl] * np.exp(
                -(nu_in[:, None] - scat.nu0[sl]) ** 2 / (2 * scat.sig_nu[sl] ** 2))
            Ex *= g * w[sl]
            va = scat.vib_amp[sl]
            if np.any(va != 0):
                d = va * np.sin(2 * np.pi * scat.vib_freq[sl] * t + scat.vib_phase[sl])
                Ex *= np.exp(-4j * np.pi * d / geom.lam)
            Er = np.exp(-2j * np.pi * np.outer(r_img[sl], kr_in))
            S += Ex.astype(np.complex64) @ Er.astype(np.complex64)
            if verbose and (i // chunk) % 10 == 0:
                print(f"  synth {i}/{scat.n}", flush=True)
    if ground_intensity > 0:
        rng = np.random.default_rng(seed)
        sigma = np.sqrt(ground_intensity * Nx * Nr / geom.band_frac ** 2 / 2)
        S += (rng.normal(0, sigma, S.shape) + 1j * rng.normal(0, sigma, S.shape)).astype(np.complex64)
    full = np.zeros((Nx, Nr), np.complex64)
    ix = np.where(in_x)[0]
    ir = np.where(in_r)[0]
    full[np.ix_(ix, ir)] = S
    return np.fft.ifft2(full).astype(np.complex64)
