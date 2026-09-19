"""End-to-end reconstructed pipeline: SLC -> sub-aperture pairs -> patch shifts -> trajectories -> depth focus."""
import time
import numpy as np
from .subap import SubapBank
from .track import patch_shifts
from .tomo import kz_for_bank, default_depths, focus_paper, focus_windows


def grid_targets(r0, r1, c0, c1, stride):
    rows = np.arange(r0, r1, stride)
    cols = np.arange(c0, c1, stride)
    RR, CC = np.meshgrid(rows, cols, indexing="ij")
    return RR.ravel(), CC.ravel(), (rows, cols)


def line_targets(r0, c0, r1, c1, n=None):
    n = int(max(abs(r1 - r0), abs(c1 - c0))) + 1 if n is None else n
    rows = np.round(np.linspace(r0, r1, n)).astype(int)
    cols = np.round(np.linspace(c0, c1, n)).astype(int)
    return rows, cols


def trajectories(slc, geom, bank, rows, cols, patch=32, upsample=1000, normalize=False, verbose=True):
    """Per-pixel sub-aperture displacement trajectories Y [P, K, 2] (azimuth px, range px)."""
    spec = np.fft.fft(slc, axis=0)
    masks = bank.masks(geom, slc.shape[0])
    P = rows.size
    Y = np.zeros((P, bank.K, 2))
    t0 = time.time()
    for k in range(bank.K):
        R, O = bank.pair(spec, k, masks)
        dr, dc = patch_shifts(R, O, rows, cols, patch=patch, upsample=upsample, normalize=normalize)
        Y[:, k, 0] = dr
        Y[:, k, 1] = dc
        if verbose and (k % 10 == 0 or k == bank.K - 1):
            print(f"  pair {k + 1}/{bank.K}  ({time.time() - t0:.0f}s)", flush=True)
    return Y


def common_mode(Y, mode="global"):
    if mode == "none":
        return np.zeros_like(Y)
    if mode == "global":
        return np.broadcast_to(np.median(Y, axis=0, keepdims=True), Y.shape)
    raise ValueError(mode)


def run_pipeline(slc, geom, rows, cols, bank=None, patch=32, upsample=1000, lam_s=0.48,
                 z=None, common="global", normalize=False, W=25, modes=("paper", "windows"),
                 verbose=True):
    bank = SubapBank() if bank is None else bank
    Y = trajectories(slc, geom, bank, rows, cols, patch, upsample, normalize, verbose)
    q = Y - common_mode(Y, common)
    kz = kz_for_bank(bank, geom, lam_s, slc.shape[0])
    z = default_depths(kz) if z is None else z
    out = {"Y": Y, "q": q, "kz": kz, "z": z, "rows": rows, "cols": cols,
           "bands": bank.bands(geom, slc.shape[0]), "lam_s": lam_s, "patch": patch}
    if "paper" in modes:
        out["tomo_paper"] = focus_paper(q[..., 0] + 1j * q[..., 1], kz, z)
    if "windows" in modes:
        out["score_w"], out["best_w"] = focus_windows(q, kz, z, W=W)
    return out
