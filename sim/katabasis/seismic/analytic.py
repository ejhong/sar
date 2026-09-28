"""Exact solutions the solver is checked against."""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq


def stokes_displacement(t: np.ndarray, force: np.ndarray, r_vec: np.ndarray, j: int,
                        vp: float, vs: float, rho: float) -> np.ndarray:
    """Displacement u_i(t), i = 0..2, from a point force F(t) along axis j in a whole space.

    Aki & Richards (2002), equation 4.23: the near-field term and the far-field
    P and S terms. `force` is sampled on `t` (uniform, starting at 0).
    """
    r = float(np.linalg.norm(r_vec))
    g = r_vec / r
    dt = t[1] - t[0]
    out = np.zeros((3, len(t)))
    F = lambda tt: np.interp(tt, t, force, left=0.0, right=0.0)
    taus = np.arange(r / vp, r / vs + dt / 2, dt / 4)
    for n, tn in enumerate(t):
        near = np.trapezoid(taus * F(tn - taus), taus) if len(taus) > 1 else 0.0
        for i in range(3):
            dij = 1.0 if i == j else 0.0
            out[i, n] = (1 / (4 * math.pi * rho)) * (3 * g[i] * g[j] - dij) / r ** 3 * near \
                + g[i] * g[j] / (4 * math.pi * rho * vp ** 2 * r) * F(tn - r / vp) \
                - (g[i] * g[j] - dij) / (4 * math.pi * rho * vs ** 2 * r) * F(tn - r / vs)
    return out


def rayleigh_speed(vp: float, vs: float) -> float:
    """Root of the Rayleigh equation for a homogeneous half-space."""
    k = (vs / vp) ** 2

    def f(x):          # x = (c / vs)^2
        return x ** 3 - 8 * x ** 2 + (24 - 16 * k) * x - 16 * (1 - k)
    return vs * math.sqrt(brentq(f, 1e-6, 1 - 1e-9))


# ------------------------------------------------ a small void in a strained half-space --

def void_moment(stress: np.ndarray, volume: float, nu: float) -> np.ndarray:
    """Moment tensor (N·m) of the far-field perturbation a spherical void makes in a remote stress (Pa).

    Eshelby (1957): a void is equivalent to an inclusion with eigenstrain e* = (I - S)^-1 e_remote, S the
    Eshelby tensor of a sphere, isotropic with volumetric and deviatoric parts a = (1 + nu) / (3 (1 - nu)) and
    b = 2 (4 - 5 nu) / (15 (1 - nu)). Its moment is V C e*, and since C and S commute for a sphere,
    M = V [ s_vol / (1 - a) + s_dev / (1 - b) ]: about 2.4 and 1.9 times the stress the void removes.
    """
    s = np.asarray(stress, float)
    vol = np.trace(s) / 3 * np.eye(3)
    a = (1 + nu) / (3 * (1 - nu))
    b = 2 * (4 - 5 * nu) / (15 * (1 - nu))
    return volume * (vol / (1 - a) + (s - vol) / (1 - b))


def okada_tensile_point(x: np.ndarray, y: np.ndarray, depth: float, dip_deg: float, potency: float,
                        lam: float, mu: float) -> np.ndarray:
    """Surface displacement (east-like x, north-like y, up) of a tensile point source: Okada (1985), BSSA 75,
    1135, eqs. 25-27 with the I terms of eq. 28 (point source, strike along x, dip dip_deg, opening times area
    `potency`, m^3). Returns [3, ...]."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    d = depth
    dl = np.deg2rad(dip_deg)
    sd, cd = np.sin(dl), np.cos(dl)
    p = y * cd + d * sd
    q = y * sd - d * cd
    R = np.sqrt(x ** 2 + y ** 2 + d ** 2)
    m = mu / (lam + mu)
    I1 = m * y * (1 / (R * (R + d) ** 2) - x ** 2 * (3 * R + d) / (R ** 3 * (R + d) ** 3))
    I2 = m * x * (1 / (R * (R + d) ** 2) - y ** 2 * (3 * R + d) / (R ** 3 * (R + d) ** 3))
    I3 = m * x / R ** 3 - I2
    I5 = m * (1 / (R * (R + d)) - x ** 2 * (2 * R + d) / (R ** 3 * (R + d) ** 2))
    c = potency / (2 * np.pi)
    return np.stack([c * (3 * x * q ** 2 / R ** 5 - I3 * sd ** 2),
                     c * (3 * y * q ** 2 / R ** 5 - I1 * sd ** 2),
                     c * (3 * d * q ** 2 / R ** 5 - I5 * sd ** 2)])


def moment_surface_displacement(M: np.ndarray, depth: float, x: np.ndarray, y: np.ndarray,
                                lam: float, mu: float) -> np.ndarray:
    """Surface displacement [3, ...] (x, y, up) of a buried point source with symmetric moment tensor M in the
    (x, y, up) frame, at `depth` below the surface under (0, 0). M is resolved on its principal axes into three
    orthogonal tensile sources (a tensile source of potency P and unit normal n has moment P (lam I + 2 mu n n)),
    each evaluated with Okada's point source in a frame whose strike is horizontal."""
    w, V = np.linalg.eigh(np.asarray(M, float))
    S = w.sum() / (3 * lam + 2 * mu)
    out = np.zeros((3,) + np.shape(x))
    for k in range(3):
        P = (w[k] - lam * S) / (2 * mu)
        n = V[:, k]
        if n[2] < 0:
            n = -n
        dip = np.rad2deg(np.arccos(np.clip(n[2], -1, 1)))          # a plane's dip equals its normal's tilt
        h = n[:2]
        if np.linalg.norm(h) < 1e-12:
            h = np.array([0.0, -1.0])                              # horizontal plane: any strike will do
        h = h / np.linalg.norm(h)
        # Okada's frame: strike along x', and the normal's horizontal part along -y' (plane dips towards +y')
        ey = -h
        ex = np.array([-ey[1], ey[0]])
        if np.cross(np.r_[ex, 0], np.r_[ey, 0])[2] < 0:
            ex = -ex
        xs = x * ex[0] + y * ex[1]
        ys = x * ey[0] + y * ey[1]
        u = okada_tensile_point(xs, ys, depth, dip, P, lam, mu)
        out[0] += u[0] * ex[0] + u[1] * ey[0]
        out[1] += u[0] * ex[1] + u[1] * ey[1]
        out[2] += u[2]
    return out
