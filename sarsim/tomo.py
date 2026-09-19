"""Trajectory -> 'depth' focusing, in the two published variants.

paper variant (Biondi & Malanga 2022, eqs 21-24): the K-sample complex trajectory
Y_k = d_az,k + j d_rg,k is projected on steering vectors a_k(z) = exp(j Kz_k z):
h(z) = a(z)^H Y, tomogram = |h(z)|^2. This is a non-uniform DFT of the trajectory; z is the
frequency axis relabelled with Kz.

windows variant (2026 replication protocol v1.5, branch A): for each contiguous window of W
samples the real two-component trajectory q_k is least-squares fitted to
[cos(Kz_k z), sin(Kz_k z)] and scored by adjusted R^2 (n = 2W, p = 4); the pixel's score at
depth z is the maximum over windows.
"""
import numpy as np


def kz_for_bank(bank, geom, lam_s, Nx=None):
    return geom.kz(bank.bands(geom, Nx)["nu_c"], lam_s)


def nyquist_depth(kz):
    dk = np.median(np.diff(kz))
    return np.pi / abs(dk)


def default_depths(kz, n=160, zmax=None):
    """Depth grid from just above zero to 97% of the Nyquist depth (the steering basis is
    exactly degenerate at the Nyquist depth: cos and sin columns become proportional)."""
    zmax = 0.97 * nyquist_depth(kz) if zmax is None else zmax
    return np.linspace(zmax / n, zmax, n)


def focus_paper(Yc, kz, z):
    """|a(z)^H Y|^2 / K^2 for complex trajectories Yc [P, K]."""
    A = np.exp(-1j * np.outer(kz, z))
    H = Yc @ A
    return (np.abs(H) ** 2) / Yc.shape[1] ** 2


def focus_windows(q, kz, z, W=25, chunk=4096):
    """Max-over-windows adjusted R^2 of the steering fit; q [P, K, 2] real. Returns (score [P, Z], best window [P, Z])."""
    P, K, _ = q.shape
    nW = K - W + 1
    Z = z.size
    n, p = 2 * W, 4
    best = np.full((P, Z), -np.inf)
    best_w = np.zeros((P, Z), np.int16)
    for w in range(nW):
        kzw = kz[w:w + W]
        X = np.stack([np.cos(np.outer(z, kzw)), np.sin(np.outer(z, kzw))], axis=2)  # [Z, W, 2]
        G = np.einsum('zki,zkj->zij', X, X)                                        # [Z, 2, 2]
        Greg = G + 1e-9 * np.trace(G, axis1=1, axis2=2)[:, None, None] * np.eye(2)[None]
        Pinv = np.linalg.solve(Greg, np.transpose(X, (0, 2, 1)))                  # [Z, 2, W]
        for i in range(0, P, chunk):
            Qw = q[i:i + chunk, w:w + W, :]                                        # [p, W, 2]
            sstot = ((Qw - Qw.mean(axis=1, keepdims=True)) ** 2).sum(axis=(1, 2))  # centred energy
            ssq = (Qw ** 2).sum(axis=(1, 2))
            B = np.einsum('zik,pkj->pzij', Pinv, Qw)                               # [p, Z, 2, 2]
            fitE = np.einsum('pzij,zik,pzkj->pz', B, G, B)                          # ||X B||^2
            ssres = np.clip(ssq[:, None] - fitE, 0, None)
            r2 = 1 - ssres / np.maximum(sstot[:, None], 1e-30)
            adj = 1 - (1 - r2) * (n - 1) / (n - p - 1)
            better = adj > best[i:i + chunk]
            best[i:i + chunk] = np.where(better, adj, best[i:i + chunk])
            best_w[i:i + chunk] = np.where(better, w, best_w[i:i + chunk])
    return best, best_w


def permute(q, rng):
    """Shuffle the sub-aperture order (same permutation for all pixels)."""
    perm = rng.permutation(q.shape[1])
    return q[:, perm, :]
