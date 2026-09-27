"""Three-dimensional isotropic elastic waves: velocity-stress finite differences.

Staggered grid (Virieux 1986; Graves 1996), fourth order in space and second
order in time. Arrays are indexed [ix, iy, iz] with iz growing downward, on
the composition Grid (iz = 0 is the top layer of cells).

    vx (i+½, j, k)   vy (i, j+½, k)   vz (i, j, k+½)
    sxx syy szz (i, j, k)
    sxy (i+½, j+½, k)   sxz (i+½, j, k+½)   syz (i, j+½, k+½)

Voids, the free surface and terrain are all handled by the vacuum
formulation: air cells carry λ = μ = ρ = 0; buoyancy at a velocity point is
the inverse of the arithmetic mean of its two densities (zero between two
air cells); shear moduli at shear points are harmonic means of the four
surrounding cells, so any contact with air is traction-free. Absorbing
boundaries are convolutional PMLs (Komatitsch & Martin 2007) whose memory
variables live only in the boundary slabs.

Sign conventions: +vz is downward (along +iz). Sources and receivers are
given in site coordinates (x east, y north, z up) and converted here.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numba as nb
import numpy as np

from ..compose.grid import Grid

C1 = 9.0 / 8.0
C2 = -1.0 / 24.0
CFL_LIMIT = 6.0 / (7.0 * math.sqrt(3.0))       # O(2,4) in 3-D


# ---------------------------------------------------------------- model --

@dataclass
class Medium:
    """Material parameters staggered for the scheme (float32)."""
    grid: Grid
    lam: np.ndarray           # at nodes
    lam2mu: np.ndarray        # λ + 2μ at nodes
    mu_xy: np.ndarray
    mu_xz: np.ndarray
    mu_yz: np.ndarray
    bx: np.ndarray            # buoyancy at vx points
    by: np.ndarray
    bz: np.ndarray
    vp_max: float
    vs_min: float             # smallest non-zero S speed
    solid: np.ndarray         # bool at nodes
    near: np.ndarray          # bool: air within two cells (second-order stencils there)

    @classmethod
    def from_arrays(cls, grid: Grid, vp: np.ndarray, vs: np.ndarray, rho: np.ndarray,
                    air: np.ndarray | None = None) -> 'Medium':
        vp = vp.astype(np.float64)
        vs = vs.astype(np.float64)
        rho = rho.astype(np.float64).copy()
        air = (vs <= 0) & (vp < 400) if air is None else air
        rho[air] = 0.0
        vp = np.where(air, 0.0, vp)
        vs = np.where(air, 0.0, vs)
        mu = rho * vs ** 2
        lam = rho * vp ** 2 - 2 * mu

        def avg_rho(axis):
            r = rho
            nxt = np.roll(r, -1, axis=axis)
            m = 0.5 * (r + nxt)
            b = np.where(m > 0, 1.0 / np.where(m > 0, m, 1.0), 0.0)
            sl = [slice(None)] * 3
            sl[axis] = -1
            b[tuple(sl)] = 0.0                        # no wrap-around at the far edge
            return b.astype(np.float32)

        def harm(a, axes):
            vals = [a]
            for ax in axes:
                vals = vals + [np.roll(v, -1, axis=ax) for v in vals]
            stack = np.stack(vals)
            with np.errstate(divide='ignore'):
                inv = np.where(stack > 0, 1.0 / np.where(stack > 0, stack, 1.0), np.inf)
            h = len(vals) / inv.sum(0)
            h[~np.isfinite(h)] = 0.0
            for ax in axes:
                sl = [slice(None)] * 3
                sl[ax] = -1
                h[tuple(sl)] = 0.0
            return h.astype(np.float32)

        solid = ~air
        near = np.zeros(air.shape, bool)
        if air.any():
            from scipy.ndimage import binary_dilation
            near = binary_dilation(air, structure=np.ones((5, 5, 5), bool)) & solid
        return cls(grid, lam.astype(np.float32), (lam + 2 * mu).astype(np.float32),
                   harm(mu, (0, 1)), harm(mu, (0, 2)), harm(mu, (1, 2)),
                   avg_rho(0), avg_rho(1), avg_rho(2),
                   float(vp[solid].max()), float(vs[solid & (vs > 0)].min()) if np.any(solid & (vs > 0)) else 0.0,
                   solid, near)

    @classmethod
    def from_model(cls, model) -> 'Medium':
        """From a voxelised composition (katabasis.compose.voxel.Model)."""
        return cls.from_arrays(model.grid, model.vp, model.vs, model.rho, model.air)


# ------------------------------------------------------------------ PML --

@dataclass
class PML:
    width: int
    faces: tuple[bool, bool, bool, bool, bool, bool]     # xlo xhi ylo yhi ztop zbottom
    # per axis: (a, b) at nodes and at half points, length n
    ax: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    ay: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    az: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]

    @staticmethod
    def profile(n: int, width: int, lo: bool, hi: bool, dx: float, dt: float, vp: float, f0: float,
                reflection: float = 1e-4, power: float = 2.0):
        L = width * dx
        d0 = -(power + 1) * vp * math.log(reflection) / (2 * L)
        amax = math.pi * f0
        out = []
        for shift in (0.0, 0.5):
            pos = (np.arange(n) + shift) * dx
            dist = np.zeros(n)
            if lo:
                dist = np.maximum(dist, (L - pos) / L)
            if hi:
                dist = np.maximum(dist, (pos - (n - 1) * dx + L) / L)
            dist = np.clip(dist, 0, 1)
            d = d0 * dist ** power
            alpha = amax * (1 - dist) * (dist > 0)
            b = np.exp(-(d + alpha) * dt)
            with np.errstate(invalid='ignore', divide='ignore'):
                a = np.where(d > 0, d * (b - 1) / (d + alpha), 0.0)
            out += [a.astype(np.float32), b.astype(np.float32)]
        return tuple(out)      # a_node, b_node, a_half, b_half


# --------------------------------------------------------------- kernels --

@nb.njit(parallel=True, fastmath=True, cache=True)
def _update_velocity(vx, vy, vz, sxx, syy, szz, sxy, sxz, syz, bx, by, bz, near, dt, inv_dx,
                     axn_a, axn_b, axh_a, axh_b, ayn_a, ayn_b, ayh_a, ayh_b, azn_a, azn_b, azh_a, azh_b,
                     px_sxx, py_sxy, pz_sxz, px_sxy, py_syy, pz_syz, px_sxz, py_syz, pz_szz,
                     xmap, ymap, zmap):
    nx, ny, nz = vx.shape
    for i in nb.prange(2, nx - 2):
        mi = xmap[i]
        for j in range(2, ny - 2):
            mj = ymap[j]
            for k in range(2, nz - 2):
                mk = zmap[k]
                if near[i, j, k]:
                    C1 = 1.0
                    C2 = 0.0
                else:
                    C1 = 9.0 / 8.0
                    C2 = -1.0 / 24.0
                # vx at (i+1/2, j, k)
                d1 = (C1 * (sxx[i + 1, j, k] - sxx[i, j, k]) + C2 * (sxx[i + 2, j, k] - sxx[i - 1, j, k])) * inv_dx
                d2 = (C1 * (sxy[i, j, k] - sxy[i, j - 1, k]) + C2 * (sxy[i, j + 1, k] - sxy[i, j - 2, k])) * inv_dx
                d3 = (C1 * (sxz[i, j, k] - sxz[i, j, k - 1]) + C2 * (sxz[i, j, k + 1] - sxz[i, j, k - 2])) * inv_dx
                if mi >= 0:
                    px_sxx[mi, j, k] = axh_b[i] * px_sxx[mi, j, k] + axh_a[i] * d1
                    d1 += px_sxx[mi, j, k]
                if mj >= 0:
                    py_sxy[i, mj, k] = ayn_b[j] * py_sxy[i, mj, k] + ayn_a[j] * d2
                    d2 += py_sxy[i, mj, k]
                if mk >= 0:
                    pz_sxz[i, j, mk] = azn_b[k] * pz_sxz[i, j, mk] + azn_a[k] * d3
                    d3 += pz_sxz[i, j, mk]
                vx[i, j, k] += dt * bx[i, j, k] * (d1 + d2 + d3)
                # vy at (i, j+1/2, k)
                d1 = (C1 * (sxy[i, j, k] - sxy[i - 1, j, k]) + C2 * (sxy[i + 1, j, k] - sxy[i - 2, j, k])) * inv_dx
                d2 = (C1 * (syy[i, j + 1, k] - syy[i, j, k]) + C2 * (syy[i, j + 2, k] - syy[i, j - 1, k])) * inv_dx
                d3 = (C1 * (syz[i, j, k] - syz[i, j, k - 1]) + C2 * (syz[i, j, k + 1] - syz[i, j, k - 2])) * inv_dx
                if mi >= 0:
                    px_sxy[mi, j, k] = axn_b[i] * px_sxy[mi, j, k] + axn_a[i] * d1
                    d1 += px_sxy[mi, j, k]
                if mj >= 0:
                    py_syy[i, mj, k] = ayh_b[j] * py_syy[i, mj, k] + ayh_a[j] * d2
                    d2 += py_syy[i, mj, k]
                if mk >= 0:
                    pz_syz[i, j, mk] = azn_b[k] * pz_syz[i, j, mk] + azn_a[k] * d3
                    d3 += pz_syz[i, j, mk]
                vy[i, j, k] += dt * by[i, j, k] * (d1 + d2 + d3)
                # vz at (i, j, k+1/2)
                d1 = (C1 * (sxz[i, j, k] - sxz[i - 1, j, k]) + C2 * (sxz[i + 1, j, k] - sxz[i - 2, j, k])) * inv_dx
                d2 = (C1 * (syz[i, j, k] - syz[i, j - 1, k]) + C2 * (syz[i, j + 1, k] - syz[i, j - 2, k])) * inv_dx
                d3 = (C1 * (szz[i, j, k + 1] - szz[i, j, k]) + C2 * (szz[i, j, k + 2] - szz[i, j, k - 1])) * inv_dx
                if mi >= 0:
                    px_sxz[mi, j, k] = axn_b[i] * px_sxz[mi, j, k] + axn_a[i] * d1
                    d1 += px_sxz[mi, j, k]
                if mj >= 0:
                    py_syz[i, mj, k] = ayn_b[j] * py_syz[i, mj, k] + ayn_a[j] * d2
                    d2 += py_syz[i, mj, k]
                if mk >= 0:
                    pz_szz[i, j, mk] = azh_b[k] * pz_szz[i, j, mk] + azh_a[k] * d3
                    d3 += pz_szz[i, j, mk]
                vz[i, j, k] += dt * bz[i, j, k] * (d1 + d2 + d3)


@nb.njit(parallel=True, fastmath=True, cache=True)
def _update_stress(vx, vy, vz, sxx, syy, szz, sxy, sxz, syz, lam, lam2mu, mxy, mxz, myz, near, dt, inv_dx,
                   axn_a, axn_b, axh_a, axh_b, ayn_a, ayn_b, ayh_a, ayh_b, azn_a, azn_b, azh_a, azh_b,
                   px_vx, py_vy, pz_vz, py_vx, px_vy, pz_vx, px_vz, pz_vy, py_vz,
                   xmap, ymap, zmap):
    nx, ny, nz = vx.shape
    for i in nb.prange(2, nx - 2):
        mi = xmap[i]
        for j in range(2, ny - 2):
            mj = ymap[j]
            for k in range(2, nz - 2):
                mk = zmap[k]
                if near[i, j, k]:
                    C1 = 1.0
                    C2 = 0.0
                else:
                    C1 = 9.0 / 8.0
                    C2 = -1.0 / 24.0
                # normal stresses at nodes
                exx = (C1 * (vx[i, j, k] - vx[i - 1, j, k]) + C2 * (vx[i + 1, j, k] - vx[i - 2, j, k])) * inv_dx
                eyy = (C1 * (vy[i, j, k] - vy[i, j - 1, k]) + C2 * (vy[i, j + 1, k] - vy[i, j - 2, k])) * inv_dx
                ezz = (C1 * (vz[i, j, k] - vz[i, j, k - 1]) + C2 * (vz[i, j, k + 1] - vz[i, j, k - 2])) * inv_dx
                if mi >= 0:
                    px_vx[mi, j, k] = axn_b[i] * px_vx[mi, j, k] + axn_a[i] * exx
                    exx += px_vx[mi, j, k]
                if mj >= 0:
                    py_vy[i, mj, k] = ayn_b[j] * py_vy[i, mj, k] + ayn_a[j] * eyy
                    eyy += py_vy[i, mj, k]
                if mk >= 0:
                    pz_vz[i, j, mk] = azn_b[k] * pz_vz[i, j, mk] + azn_a[k] * ezz
                    ezz += pz_vz[i, j, mk]
                l = lam[i, j, k]
                l2 = lam2mu[i, j, k]
                sxx[i, j, k] += dt * (l2 * exx + l * (eyy + ezz))
                syy[i, j, k] += dt * (l2 * eyy + l * (exx + ezz))
                szz[i, j, k] += dt * (l2 * ezz + l * (exx + eyy))
                # sxy at (i+1/2, j+1/2, k)
                m = mxy[i, j, k]
                if m != 0.0:
                    d1 = (C1 * (vx[i, j + 1, k] - vx[i, j, k]) + C2 * (vx[i, j + 2, k] - vx[i, j - 1, k])) * inv_dx
                    d2 = (C1 * (vy[i + 1, j, k] - vy[i, j, k]) + C2 * (vy[i + 2, j, k] - vy[i - 1, j, k])) * inv_dx
                    if mj >= 0:
                        py_vx[i, mj, k] = ayh_b[j] * py_vx[i, mj, k] + ayh_a[j] * d1
                        d1 += py_vx[i, mj, k]
                    if mi >= 0:
                        px_vy[mi, j, k] = axh_b[i] * px_vy[mi, j, k] + axh_a[i] * d2
                        d2 += px_vy[mi, j, k]
                    sxy[i, j, k] += dt * m * (d1 + d2)
                # sxz at (i+1/2, j, k+1/2)
                m = mxz[i, j, k]
                if m != 0.0:
                    d1 = (C1 * (vx[i, j, k + 1] - vx[i, j, k]) + C2 * (vx[i, j, k + 2] - vx[i, j, k - 1])) * inv_dx
                    d2 = (C1 * (vz[i + 1, j, k] - vz[i, j, k]) + C2 * (vz[i + 2, j, k] - vz[i - 1, j, k])) * inv_dx
                    if mk >= 0:
                        pz_vx[i, j, mk] = azh_b[k] * pz_vx[i, j, mk] + azh_a[k] * d1
                        d1 += pz_vx[i, j, mk]
                    if mi >= 0:
                        px_vz[mi, j, k] = axh_b[i] * px_vz[mi, j, k] + axh_a[i] * d2
                        d2 += px_vz[mi, j, k]
                    sxz[i, j, k] += dt * m * (d1 + d2)
                # syz at (i, j+1/2, k+1/2)
                m = myz[i, j, k]
                if m != 0.0:
                    d1 = (C1 * (vy[i, j, k + 1] - vy[i, j, k]) + C2 * (vy[i, j, k + 2] - vy[i, j, k - 1])) * inv_dx
                    d2 = (C1 * (vz[i, j + 1, k] - vz[i, j, k]) + C2 * (vz[i, j + 2, k] - vz[i, j - 1, k])) * inv_dx
                    if mk >= 0:
                        pz_vy[i, j, mk] = azh_b[k] * pz_vy[i, j, mk] + azh_a[k] * d1
                        d1 += pz_vy[i, j, mk]
                    if mj >= 0:
                        py_vz[i, mj, k] = ayh_b[j] * py_vz[i, mj, k] + ayh_a[j] * d2
                        d2 += py_vz[i, mj, k]
                    syz[i, j, k] += dt * m * (d1 + d2)


# ------------------------------------------------------ sources/receivers --

def ricker(t: np.ndarray, f0: float, t0: float | None = None) -> np.ndarray:
    t0 = 1.2 / f0 if t0 is None else t0
    a = (math.pi * f0 * (t - t0)) ** 2
    return (1 - 2 * a) * np.exp(-a)


@dataclass
class Source:
    """A point force (N) or an explosion (moment, N·m) at a site-coordinate position."""
    position: tuple[float, float, float]
    wavelet: np.ndarray                  # one value per time step
    kind: str = 'force'                  # 'force' | 'explosion'
    direction: tuple[float, float, float] = (0.0, 0.0, -1.0)    # site coords; default: pushing down


@dataclass
class Receivers:
    """Three-component velocity recorders at site-coordinate positions."""
    positions: np.ndarray                # (n, 3)
    names: list[str] = field(default_factory=list)


def _nearest(grid: Grid, p, half: tuple[float, float, float] = (0, 0, 0)) -> tuple[int, int, int]:
    """Index of the grid point nearest p for a field staggered by `half` cells (x, y, depth)."""
    x0, y0, ztop = grid.origin
    h = grid.spacing
    i = int(round((p[0] - x0) / h - half[0]))
    j = int(round((p[1] - y0) / h - half[1]))
    k = int(round((ztop - p[2]) / h - half[2]))
    return i, j, k


# ----------------------------------------------------------------- solver --

@dataclass
class Result:
    t: np.ndarray                         # (nt_rec,)
    traces: np.ndarray                    # (n_receivers, 3, nt_rec), site axes: east, north, up
    snapshots: np.ndarray | None          # (n_snap, nx, ny) up velocity at the ground
    snap_t: np.ndarray | None
    dt: float
    meta: dict


class Simulation:
    def __init__(self, medium: Medium, dt: float | None = None, pml_width: int = 12,
                 faces: tuple[bool, ...] = (True, True, True, True, False, True), f0: float = 30.0,
                 cfl: float = 0.42, second_order_near_air: bool = True, pml_vp: float | None = None):
        """`pml_vp` tunes the absorbing layers (default: the medium's fastest speed); hold it fixed when
        the medium changes between runs that must compare, as in an inversion."""
        self.m = medium
        g = medium.grid
        self.grid = g
        h = g.spacing
        self.dt = dt if dt is not None else cfl * h / medium.vp_max
        if medium.vp_max * self.dt / h > CFL_LIMIT:
            raise ValueError(f'unstable: CFL {medium.vp_max * self.dt / h:.3f} > {CFL_LIMIT:.3f}')
        nx, ny, nz = g.shape
        self.shape = g.shape
        if not faces[4] and medium.solid.any():
            # the updates skip the outer two layers: the velocity just above the ground must be one they reach
            top = np.argmax(medium.solid, axis=2)[medium.solid.any(axis=2)]
            if top.min() < 3:
                raise ValueError('the ground needs at least three layers of air cells above it')
        self.near = medium.near if second_order_near_air else np.zeros(g.shape, bool)
        self.pml = pml_width
        self.faces = tuple(faces)
        vref = medium.vp_max if pml_vp is None else pml_vp
        prof = lambda n, lo, hi: PML.profile(n, pml_width, lo, hi, h, self.dt, vref, f0)
        self.ax = prof(nx, faces[0], faces[1])
        self.ay = prof(ny, faces[2], faces[3])
        self.az = prof(nz, faces[4], faces[5])
        self.xmap, nxm = self._map(nx, faces[0], faces[1])
        self.ymap, nym = self._map(ny, faces[2], faces[3])
        self.zmap, nzm = self._map(nz, faces[4], faces[5])
        f32 = np.float32
        z = lambda *s: np.zeros(s, f32)
        self.v = [z(nx, ny, nz) for _ in range(3)]
        self.s = [z(nx, ny, nz) for _ in range(6)]
        self.pv = [z(max(nxm, 1), ny, nz), z(nx, max(nym, 1), nz), z(nx, ny, max(nzm, 1)),
                   z(max(nxm, 1), ny, nz), z(nx, max(nym, 1), nz), z(nx, ny, max(nzm, 1)),
                   z(max(nxm, 1), ny, nz), z(nx, max(nym, 1), nz), z(nx, ny, max(nzm, 1))]
        self.ps = [z(max(nxm, 1), ny, nz), z(nx, max(nym, 1), nz), z(nx, ny, max(nzm, 1)),
                   z(nx, max(nym, 1), nz), z(max(nxm, 1), ny, nz), z(nx, ny, max(nzm, 1)),
                   z(max(nxm, 1), ny, nz), z(nx, ny, max(nzm, 1)), z(nx, max(nym, 1), nz)]

    def _map(self, n: int, lo: bool, hi: bool) -> tuple[np.ndarray, int]:
        w = self.pml
        m = -np.ones(n, np.int64)
        c = 0
        for i in range(n):
            if (lo and i < w) or (hi and i >= n - w):
                m[i] = c
                c += 1
        return m, c

    def points_per_wavelength(self, fmax: float) -> float:
        return self.m.vs_min / fmax / self.grid.spacing if self.m.vs_min > 0 else math.inf

    def run(self, sources: list[Source], receivers: Receivers, nt: int, record_every: int = 1,
            snapshot_every: int = 0, surface_index: np.ndarray | None = None, progress: bool = False,
            on_step=None, on_step_every: int = 1) -> Result:
        """Advance nt steps. `on_step(it, sim)` is called after step `it` every `on_step_every` steps."""
        g = self.grid
        h = g.spacing
        dt = np.float32(self.dt)
        inv_dx = np.float32(1.0 / h)
        vx, vy, vz = self.v
        sxx, syy, szz, sxy, sxz, syz = self.s
        m = self.m
        # sources: velocity points for forces, nodes for explosions
        src = []
        for s in sources:
            if len(s.wavelet) < nt:
                raise ValueError('source wavelet shorter than the run')
            if s.kind == 'force':
                # Split between the two staggered points around the nearest node: the exact
                # adjoint of how receivers read the grid, so sources and receivers commute.
                fx, fy, fzup = s.direction
                i, j, k = _nearest(g, s.position)
                for comp, w in ((0, fx), (1, fy), (2, -fzup)):
                    if w == 0:
                        continue
                    for off in (0, -1):
                        ii, jj, kk = (i + off, j, k) if comp == 0 else (i, j + off, k) if comp == 1 else (i, j, k + off)
                        b = [m.bx, m.by, m.bz][comp][ii, jj, kk]
                        src.append(('v', comp, ii, jj, kk, np.float32(0.5 * w) * b / h ** 3, s.wavelet))
            else:
                i, j, k = _nearest(g, s.position)
                src.append(('s', -1, i, j, k, np.float32(-1.0 / h ** 3), s.wavelet))
        # receivers: average the two staggered samples around each node
        rec = [(_nearest(g, p)) for p in receivers.positions]
        nrec = len(rec)
        nt_rec = (nt + record_every - 1) // record_every
        traces = np.zeros((nrec, 3, nt_rec), np.float32)
        ri = np.array([r[0] for r in rec]), np.array([r[1] for r in rec]), np.array([r[2] for r in rec])
        snaps, snap_t = [], []
        if snapshot_every and surface_index is None:
            raise ValueError('snapshots need the surface index (first solid cell under each column)')
        pa = (*self.ax, *self.ay, *self.az)
        for it in range(nt):
            for kind, comp, i, j, k, amp, w in src:
                if kind == 'v':
                    self.v[comp][i, j, k] += dt * amp * w[it]
            _update_velocity(vx, vy, vz, sxx, syy, szz, sxy, sxz, syz, m.bx, m.by, m.bz, self.near, dt, inv_dx, *pa,
                             *self.pv, self.xmap, self.ymap, self.zmap)
            for kind, comp, i, j, k, amp, w in src:
                if kind == 's':
                    e = dt * amp * w[it]
                    sxx[i, j, k] += e
                    syy[i, j, k] += e
                    szz[i, j, k] += e
            _update_stress(vx, vy, vz, sxx, syy, szz, sxy, sxz, syz, m.lam, m.lam2mu, m.mu_xy, m.mu_xz, m.mu_yz,
                           self.near, dt, inv_dx, *pa, *self.ps, self.xmap, self.ymap, self.zmap)
            if it % record_every == 0:
                n = it // record_every
                I, J, K = ri
                traces[:, 0, n] = 0.5 * (vx[I, J, K] + vx[I - 1, J, K])
                traces[:, 1, n] = 0.5 * (vy[I, J, K] + vy[I, J - 1, K])
                traces[:, 2, n] = -0.5 * (vz[I, J, K] + vz[I, J, K - 1])       # up is positive
            if snapshot_every and it % snapshot_every == 0:
                K = surface_index
                I, J = np.meshgrid(np.arange(g.shape[0]), np.arange(g.shape[1]), indexing='ij')
                snaps.append(-0.5 * (vz[I, J, K] + vz[I, J, np.maximum(K - 1, 0)]))
                snap_t.append((it + 1) * self.dt)
            if on_step is not None and it % on_step_every == 0:
                on_step(it, self)
            if progress and it % max(1, nt // 10) == 0:
                print(f'  step {it}/{nt}', flush=True)
        t = (np.arange(nt_rec) * record_every + 1) * self.dt      # velocities live at the half step; ignore
        return Result(t, traces, np.array(snaps, np.float32) if snaps else None,
                      np.array(snap_t) if snap_t else None, self.dt,
                      {'dx': h, 'dt': self.dt, 'nt': nt, 'shape': list(g.shape), 'pml': self.pml,
                       'faces': list(self.faces), 'vp_max': m.vp_max, 'vs_min': m.vs_min})

    def energy(self) -> float:
        """Kinetic energy proxy: sum of v² over the solid cells (for decay checks)."""
        return float(sum(np.sum(v.astype(np.float64) ** 2) for v in self.v))

    def reset(self):
        for a in (*self.v, *self.s, *self.pv, *self.ps):
            a[...] = 0


def surface_index(medium: Medium) -> np.ndarray:
    """Index of the first solid cell from the top in each column."""
    solid = medium.solid
    k = np.argmax(solid, axis=2)
    return k.astype(np.int64)
