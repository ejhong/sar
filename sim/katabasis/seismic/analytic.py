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


def _mandel(T):
    """A minor-symmetric fourth-order tensor [3, 3, 3, 3] as a 6 x 6 matrix in Mandel's notation (11, 22, 33, 23, 13, 12,
    shear entries scaled by sqrt 2), where composition and inversion are matrix products and inverses."""
    idx = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    w = [1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2)]
    return np.array([[T[a, b, c, d] * w[m] * w[n] for n, (c, d) in enumerate(idx)] for m, (a, b) in enumerate(idx)])


def _vec(e):
    r2 = np.sqrt(2)
    return np.array([e[0, 0], e[1, 1], e[2, 2], r2 * e[1, 2], r2 * e[0, 2], r2 * e[0, 1]])


def _ten(v):
    r2 = np.sqrt(2)
    return np.array([[v[0], v[5] / r2, v[4] / r2], [v[5] / r2, v[1], v[3] / r2], [v[4] / r2, v[3] / r2, v[2]]])


def eshelby_tensor(shape: str, nu: float) -> np.ndarray:
    """Eshelby's tensor S [3, 3, 3, 3] for a sphere, or for a circular cylinder along z (infinitely long; Mura 1987,
    eq. 11.22): S_1111 = S_2222 = (5 - 4 nu) / (8 (1 - nu)), S_1122 = S_2211 = (4 nu - 1) / (8 (1 - nu)),
    S_1133 = S_2233 = nu / (2 (1 - nu)), S_1212 = (3 - 4 nu) / (8 (1 - nu)), S_1313 = S_2323 = 1 / 4, S_33kl = 0."""
    S = np.zeros((3, 3, 3, 3))
    d = np.eye(3)
    if shape == 'sphere':
        a = (1 + nu) / (3 * (1 - nu))
        b = 2 * (4 - 5 * nu) / (15 * (1 - nu))
        for i, j, k, l in np.ndindex(3, 3, 3, 3):
            S[i, j, k, l] = a / 3 * d[i, j] * d[k, l] + b * (0.5 * (d[i, k] * d[j, l] + d[i, l] * d[j, k]) - d[i, j] * d[k, l] / 3)
        return S
    if shape != 'cylinder':
        raise ValueError(shape)
    c = 1 / (8 * (1 - nu))
    S[0, 0, 0, 0] = S[1, 1, 1, 1] = (5 - 4 * nu) * c
    S[0, 0, 1, 1] = S[1, 1, 0, 0] = (4 * nu - 1) * c
    S[0, 0, 2, 2] = S[1, 1, 2, 2] = nu / (2 * (1 - nu))
    for (i, j) in ((0, 1),):
        S[i, j, i, j] = S[j, i, j, i] = S[i, j, j, i] = S[j, i, i, j] = (3 - 4 * nu) * c
    for (i, j) in ((0, 2), (1, 2)):
        S[i, j, i, j] = S[j, i, j, i] = S[i, j, j, i] = S[j, i, i, j] = 0.25
    return S


def void_moment_from_strain(strain: np.ndarray, measure: float, lam: float, mu: float, shape: str) -> np.ndarray:
    """Moment tensor of a void in a remote strain: Eshelby's equivalent inclusion, eigenstrain e* = (I - S)^-1 e (the
    stress inside the void vanishes), moment measure x C e*, `measure` the volume for a sphere and the cross-section's
    area for a cylinder along z (a moment per unit length). For a sphere it equals void_moment(C e, volume, nu)."""
    nu = lam / (2 * (lam + mu))
    Sm = _mandel(eshelby_tensor(shape, nu))
    e_star = _ten(np.linalg.solve(np.eye(6) - Sm, _vec(np.asarray(strain, float))))
    return measure * (lam * np.trace(e_star) * np.eye(3) + 2 * mu * e_star)


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
