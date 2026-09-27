"""First-arrival travel-time tomography in three dimensions.

Forward: the eikonal equation solved by second-order fast marching
(scikit-fmm) from every source. Rays: steepest descent on the travel-time
field from each receiver back to its source; their lengths in each cell form
the Fréchet matrix. Inverse: regularised Gauss-Newton on slowness, each step
solved by LSQR with smoothing (a discrete Laplacian) and damping toward the
starting model. Nothing here knows the true model.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numba as nb
import numpy as np
import scipy.sparse as sp
import skfmm
from scipy.ndimage import map_coordinates
from scipy.sparse.linalg import lsqr

from ..compose.grid import Grid


@dataclass
class TomoGrid:
    """A regular inversion grid: cell centres, x fastest in (nx, ny, nz); z index grows downward."""
    grid: Grid
    active: np.ndarray            # bool (nx, ny, nz): cells below the ground

    @property
    def shape(self):
        return self.grid.shape

    def to_index(self, p: np.ndarray) -> np.ndarray:
        """Site coordinates (n, 3) -> fractional cell indices (n, 3)."""
        x0, y0, ztop = self.grid.origin
        h = self.grid.spacing
        return np.stack([(p[:, 0] - x0) / h, (p[:, 1] - y0) / h, (ztop - p[:, 2]) / h], 1)


def travel_time_field(slowness: np.ndarray, active: np.ndarray, h: float, src_idx: np.ndarray) -> np.ndarray:
    """Travel time from a point source (fractional index) to every cell centre."""
    I, J, K = np.meshgrid(*(np.arange(n) for n in slowness.shape), indexing='ij', sparse=True)
    d = np.sqrt((I - src_idx[0]) ** 2 + (J - src_idx[1]) ** 2 + (K - src_idx[2]) ** 2) * h
    r0 = 0.9 * h
    phi = np.ma.MaskedArray(d - r0, mask=~active)
    speed = np.ma.MaskedArray(1.0 / slowness, mask=~active)
    t = skfmm.travel_time(phi, speed, dx=h, order=2)
    s0 = map_coordinates(slowness, src_idx[:, None], order=1, mode='nearest')[0]
    t = np.ma.filled(t, np.nan) + r0 * s0
    inside = (d <= r0) & active
    t[inside] = d[inside] * s0
    return t


@nb.njit(cache=True)
def _trace(T, gx, gy, gz, start, src, h, max_steps, lengths_idx, lengths_val):
    """Steepest descent from start to src (indices); accumulate path length per cell."""
    nx, ny, nz = T.shape
    p0, p1, p2 = start[0], start[1], start[2]
    ds = 0.3
    n = 0
    for _ in range(max_steps):
        r2 = (p0 - src[0]) ** 2 + (p1 - src[1]) ** 2 + (p2 - src[2]) ** 2
        if r2 < 1.0:
            # finish with a straight segment to the source
            steps = int(np.sqrt(r2) / ds) + 1
            for s in range(steps):
                f = (s + 0.5) / steps
                q0 = p0 + (src[0] - p0) * f
                q1 = p1 + (src[1] - p1) * f
                q2 = p2 + (src[2] - p2) * f
                i = min(max(int(q0 + 0.5), 0), nx - 1)
                j = min(max(int(q1 + 0.5), 0), ny - 1)
                k = min(max(int(q2 + 0.5), 0), nz - 1)
                if n < lengths_idx.shape[0]:
                    lengths_idx[n] = (i * ny + j) * nz + k
                    lengths_val[n] = np.sqrt(r2) / steps * h
                    n += 1
            return n, True
        i = min(max(int(p0 + 0.5), 1), nx - 2)
        j = min(max(int(p1 + 0.5), 1), ny - 2)
        k = min(max(int(p2 + 0.5), 1), nz - 2)
        a0, a1, a2 = gx[i, j, k], gy[i, j, k], gz[i, j, k]
        norm = np.sqrt(a0 * a0 + a1 * a1 + a2 * a2)
        if not np.isfinite(norm) or norm == 0:
            return n, False
        p0 -= ds * a0 / norm
        p1 -= ds * a1 / norm
        p2 -= ds * a2 / norm
        if n < lengths_idx.shape[0]:
            lengths_idx[n] = (i * ny + j) * nz + k
            lengths_val[n] = ds * h
            n += 1
    return n, False


def _gradients(T: np.ndarray):
    Tf = np.where(np.isfinite(T), T, np.nanmax(T) * 2)
    return np.gradient(Tf)


@dataclass
class Inversion:
    grid: TomoGrid
    slowness: np.ndarray
    history: list = field(default_factory=list)
    coverage: np.ndarray | None = None
    residual_rms: float = np.nan


def forward(slowness: np.ndarray, tg: TomoGrid, sources: np.ndarray, receivers: np.ndarray,
            pairs: np.ndarray, rays: bool = True):
    """Predicted times for (source, receiver) pairs, and the sparse Fréchet matrix."""
    h = tg.grid.spacing
    sidx = tg.to_index(sources)
    ridx = tg.to_index(receivers)
    ncell = slowness.size
    t_pred = np.full(len(pairs), np.nan)
    rows, cols, vals = [], [], []
    buf_i = np.zeros(20000, np.int64)
    buf_v = np.zeros(20000, np.float64)
    for s in np.unique(pairs[:, 0]):
        T = travel_time_field(slowness, tg.active, h, sidx[s])
        sel = np.where(pairs[:, 0] == s)[0]
        rr = ridx[pairs[sel, 1]]
        Tf = np.where(np.isfinite(T), T, 0)
        t_pred[sel] = map_coordinates(Tf, rr.T, order=1, mode='nearest')
        if not rays:
            continue
        gx, gy, gz = _gradients(T)
        for n, row in zip(rr, sel):
            m, ok = _trace(T, gx, gy, gz, n.astype(np.float64), sidx[s].astype(np.float64), h, 20000, buf_i, buf_v)
            if not ok:
                continue
            rows.append(np.full(m, row))
            cols.append(buf_i[:m].copy())
            vals.append(buf_v[:m].copy())
    G = None
    if rays:
        if rows:
            G = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(len(pairs), ncell))
        else:
            G = sp.csr_matrix((len(pairs), ncell))
    return t_pred, G


def laplacian(shape: tuple[int, int, int], active: np.ndarray) -> sp.csr_matrix:
    """A discrete Laplacian over the active cells (smoothness)."""
    nx, ny, nz = shape
    idx = np.arange(nx * ny * nz).reshape(shape)
    rows, cols, vals = [], [], []
    r = 0
    for axis in range(3):
        a = np.moveaxis(idx, axis, 0)
        act = np.moveaxis(active, axis, 0)
        c, lo, hi = a[1:-1], a[:-2], a[2:]
        ok = act[1:-1] & act[:-2] & act[2:]
        c, lo, hi = c[ok], lo[ok], hi[ok]
        n = len(c)
        rr = np.arange(r, r + n)
        rows += [rr, rr, rr]
        cols += [c, lo, hi]
        vals += [np.full(n, -2.0), np.ones(n), np.ones(n)]
        r += n
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(r, nx * ny * nz))


def invert(tg: TomoGrid, sources: np.ndarray, receivers: np.ndarray, pairs: np.ndarray, t_obs: np.ndarray,
           s0: float, iterations: int = 6, smooth: float = 2.0, damp: float = 0.05,
           bounds: tuple[float, float] = (1 / 6000, 1 / 300)) -> Inversion:
    """Regularised Gauss-Newton from a uniform starting slowness s0."""
    s = np.full(tg.shape, s0)
    s[~tg.active] = s0
    L = laplacian(tg.shape, tg.active)
    act = tg.active.ravel()
    Dm = sp.diags(act.astype(float))
    hist = []
    ok = np.isfinite(t_obs)
    cov = None
    for it in range(iterations):
        t_pred, G = forward(s, tg, sources, receivers, pairs[ok])
        r = t_obs[ok] - t_pred
        rms = float(np.sqrt(np.nanmean(r ** 2)))
        hist.append({'iteration': it, 'rms_ms': rms * 1e3})
        # scale so the regularisation weights are in units of the data misfit
        scale = np.sqrt(G.multiply(G).sum() / max(1, act.sum())) or 1.0
        A = sp.vstack([G, smooth * scale * L, damp * scale * Dm]).tocsr()
        dev = (s.ravel() - s0)
        b = np.concatenate([r, -smooth * scale * (L @ dev), -damp * scale * (Dm @ dev)])
        ds = lsqr(A, b, atol=1e-6, btol=1e-6, iter_lim=300)[0]
        s = np.clip(s + ds.reshape(tg.shape) * tg.active, *bounds)
        cov = np.asarray(G.sum(0)).reshape(tg.shape)
    t_pred, _ = forward(s, tg, sources, receivers, pairs[ok], rays=False)
    rms = float(np.sqrt(np.nanmean((t_obs[ok] - t_pred) ** 2)))
    hist.append({'iteration': iterations, 'rms_ms': rms * 1e3})
    return Inversion(tg, s, hist, cov, rms)
