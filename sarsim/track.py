"""Batched complex cross-correlation of image patches with sub-pixel refinement.

`patch_shifts` measures, at every requested pixel, the displacement of image O relative to
image R using a patch x patch window centred on the pixel: the (optionally normalised) cross
power spectrum is inverted for a coarse integer peak, then refined by successive stages of local
matrix-DFT upsampling (1/10, 1/100, ... down to 1/upsample px) and a final parabolic interpolation, following Guizar-Sicairos et al. (2008).
Returned shifts are (row, col) displacements in pixels: O(m, n) ~ R(m - d_row, n - d_col).
"""
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view


def _local_dft_matrix(N, offsets):
    k = np.fft.fftfreq(N) * N  # signed integer frequencies
    return np.exp(2j * np.pi * np.outer(offsets, k) / N)


def _peak2d(cc):
    P, M, _ = cc.shape
    idx = np.argmax(np.abs(cc).reshape(P, -1), axis=1)
    return idx // M, idx % M


def _parabolic(cc, a, b):
    """Sub-sample refinement of a 2-D peak at integer (a, b) by separable parabolic fits of |cc|."""
    P, M, _ = cc.shape
    A = np.abs(cc)
    idx = np.arange(P)
    a0 = np.clip(a, 1, M - 2)
    b0 = np.clip(b, 1, M - 2)
    fm, f0, fp = A[idx, a0 - 1, b0], A[idx, a0, b0], A[idx, a0 + 1, b0]
    den = fm - 2 * f0 + fp          # negative at a maximum
    da = np.where(np.abs(den) > 1e-300, 0.5 * (fm - fp) / np.where(den == 0, 1, den), 0.0)
    gm, g0, gp = A[idx, a0, b0 - 1], A[idx, a0, b0], A[idx, a0, b0 + 1]
    den = gm - 2 * g0 + gp
    db = np.where(np.abs(den) > 1e-300, 0.5 * (gm - gp) / np.where(den == 0, 1, den), 0.0)
    return np.clip(da, -0.5, 0.5), np.clip(db, -0.5, 0.5)


def patch_shifts(R, O, rows, cols, patch=32, upsample=1000, normalize=False, chunk=6144,
                 parabolic=True):
    rows = np.asarray(rows, int)
    cols = np.asarray(cols, int)
    N = patch
    h = N // 2
    H, Wd = R.shape
    if rows.min() - h < 0 or cols.min() - h < 0 or rows.max() - h > H - N or cols.max() - h > Wd - N:
        raise ValueError("target pixels too close to the image edge for this patch size")
    Rv = sliding_window_view(R, (N, N))
    Ov = sliding_window_view(O, (N, N))
    k = np.fft.fftfreq(N) * N
    # refinement stages: +-1 px at 0.1, then +-0.1 at 0.01, ... until the step is <= 1/upsample
    stages = []
    step = 0.1
    while True:
        mu = np.arange(-10, 11) * step
        stages.append((mu, _local_dft_matrix(N, mu)))
        if step <= 1.0 / upsample + 1e-12:
            break
        step /= 10
    d_row = np.empty(rows.size)
    d_col = np.empty(rows.size)
    for i in range(0, rows.size, chunk):
        sl = slice(i, i + chunk)
        PR = Rv[rows[sl] - h, cols[sl] - h]
        PO = Ov[rows[sl] - h, cols[sl] - h]
        X = np.fft.fft2(PR) * np.conj(np.fft.fft2(PO))
        if normalize:
            X = X / (np.abs(X) + 1e-12 * np.abs(X).max())
        cc = np.fft.ifft2(X)
        pr, pc = _peak2d(cc)
        pr = ((pr + N // 2) % N) - N // 2
        pc = ((pc + N // 2) % N) - N // 2
        peak_r = pr.astype(float)
        peak_c = pc.astype(float)
        for si, (mu, Wm) in enumerate(stages):
            ramp = np.exp(2j * np.pi * (k[None, :, None] * peak_r[:, None, None]
                                        + k[None, None, :] * peak_c[:, None, None]) / N)
            ccs = Wm @ (X * ramp) @ Wm.T
            a, b = _peak2d(ccs)
            if parabolic and si == len(stages) - 1:
                da, db = _parabolic(ccs, a, b)
                stp = mu[1] - mu[0]
                peak_r = peak_r + mu[np.clip(a, 1, mu.size - 2)] + da * stp
                peak_c = peak_c + mu[np.clip(b, 1, mu.size - 2)] + db * stp
            else:
                peak_r = peak_r + mu[a]
                peak_c = peak_c + mu[b]
        d_row[sl] = -peak_r
        d_col[sl] = -peak_c
    return d_row, d_col


def shift_image(img, d_row, d_col):
    """Translate an image by a fractional (d_row, d_col) via a spectral phase ramp (for tests)."""
    H, Wd = img.shape
    kr = np.fft.fftfreq(H)[:, None]
    kc = np.fft.fftfreq(Wd)[None, :]
    F = np.fft.fft2(img) * np.exp(-2j * np.pi * (kr * d_row + kc * d_col))
    return np.fft.ifft2(F)
