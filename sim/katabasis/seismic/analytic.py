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
