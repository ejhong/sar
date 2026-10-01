"""A buried point moment in a homogeneous half-space, at one frequency: the surface displacement by wavenumber integration.

P2-38 needs the claimed structure's imprint at the microseisms' frequency, where the shear wavelength (about 8 km at
0.2 Hz) is only some seven times the depth (1,220 m), so the quasi-static field (Okada's point source,
`analytic.moment_surface_displacement`) is an approximation whose accuracy has to be shown, and the field radiated to
the image's edge, which a static solution leaves out, has to be counted.

Method. Horizontal Fourier transform (u(x) = (2 pi)^-2 int u^(k) exp(i k.x) d^2k), z downward. In the transform the
elastic equations are a first-order system in depth for v = (u_x, u_y, u_z, t_xz, t_yz, t_zz), dv/dz = A v, whose six
modes are three decaying downward and three growing (Re of the eigenvalue below or above zero; the medium's damping, a
quality factor Q on both wave speeds, separates them, and the time convention is exp(-i omega t)). Above the source the
field holds both families, each referred to the depth where it is largest, below it only the downward family; the free
surface carries no traction; across the source v jumps by the moment's equivalent jump,
    [u_z] = M_zz / (lambda + 2 mu),   [t_xz] = i k_x (M_xx - lambda / (lambda + 2 mu) M_zz) + i k_y M_xy,
    [t_yz] = i k_x M_xy + i k_y (M_yy - lambda / (lambda + 2 mu) M_zz),   [u_x] = [u_y] = [t_zz] = 0,
for a moment with no xz or yz part (the force -M_ij d_j delta, integrated across the source plane). Nine linear equations
per wavenumber, each exponential at most one in size, so the solve is stable at every wavenumber.

The surface field then follows from the wavevector's azimuth phi: the response for a wavevector along x is computed for
the moment rotated into that frame, so u^(k, phi) is a trigonometric polynomial of order at most three in phi, and with
u^ = sum_n c_n(k) exp(i n phi), u(r, theta) = (2 pi)^-1 sum_n i^n exp(i n theta) int_0^inf c_n(k) J_n(k r) k dk. The
radial integrals are trapezoidal on a grid refined about the Rayleigh pole and the two branch points.

Checked (tests): the low-frequency limit against Okada's point source for every moment component.
"""
from __future__ import annotations

import numpy as np
from scipy.special import jv

ORDERS = np.arange(-3, 4)
N_PHI = 8


def _system(kx, ky, omega, lam, mu, rho):
    """A(k) for v = (u_x, u_y, u_z, t_xz, t_yz, t_zz), batched over kx, ky (z downward; d/dx -> i k_x)."""
    kx, ky = np.broadcast_arrays(np.asarray(kx, complex), np.asarray(ky, complex))
    g = lam + 2 * mu
    c1 = 4 * mu * (lam + mu) / g
    c2 = 2 * mu * lam / g
    w2 = rho * omega ** 2
    A = np.zeros(kx.shape + (6, 6), complex)
    A[..., 0, 3] = 1 / mu
    A[..., 0, 2] = -1j * kx
    A[..., 1, 4] = 1 / mu
    A[..., 1, 2] = -1j * ky
    A[..., 2, 0] = -1j * lam * kx / g
    A[..., 2, 1] = -1j * lam * ky / g
    A[..., 2, 5] = 1 / g
    A[..., 3, 0] = c1 * kx ** 2 + mu * ky ** 2 - w2
    A[..., 3, 1] = (c2 + mu) * kx * ky
    A[..., 3, 5] = -1j * lam / g * kx
    A[..., 4, 0] = (c2 + mu) * kx * ky
    A[..., 4, 1] = c1 * ky ** 2 + mu * kx ** 2 - w2
    A[..., 4, 5] = -1j * lam / g * ky
    A[..., 5, 2] = -w2
    A[..., 5, 3] = -1j * kx
    A[..., 5, 4] = -1j * ky
    return A


def _jump(kx, ky, M, lam, mu):
    g = lam + 2 * mu
    kx, ky = np.broadcast_arrays(np.asarray(kx, complex), np.asarray(ky, complex))
    Mxx, Myy, Mzz, Mxy = M[..., 0, 0], M[..., 1, 1], M[..., 2, 2], M[..., 0, 1]
    s = np.zeros(kx.shape + (6,), complex)
    s[..., 2] = Mzz / g
    s[..., 3] = 1j * kx * (Mxx - lam / g * Mzz) + 1j * ky * Mxy
    s[..., 4] = 1j * kx * Mxy + 1j * ky * (Myy - lam / g * Mzz)
    return s


def modes(kx, ky, omega, lam, mu, rho):
    """The six modes of A(k), batched: eigenvalues and eigenvectors, the three decaying downward (Re < 0) first."""
    A = _system(kx, ky, omega, lam, mu, rho)
    lamb, W = np.linalg.eig(A)
    order = np.argsort(lamb.real, axis=-1)
    lamb = np.take_along_axis(lamb, order, axis=-1)
    W = np.take_along_axis(W, order[..., None, :], axis=-1)
    if np.any(lamb[..., :3].real >= 0) or np.any(lamb[..., 3:].real <= 0):
        raise FloatingPointError('modes not separated: add damping')
    return lamb, W


def surface_from_modes(md, depth, s):
    """u^ at the surface (x, y, down) for a source jump s [..., 6] at `depth`, from modes() (one solve per wavenumber)."""
    lamb, W = md
    ld, lu = lamb[..., :3], lamb[..., 3:]
    Wd, Wu = W[..., :, :3], W[..., :, 3:]
    Ed = np.exp(ld * depth)                                    # downward family referred to the surface, at the source
    Eu = np.exp(-lu * depth)                                   # upward family referred to the source, at the surface
    n = lamb.shape[:-1]
    S = np.zeros(n + (9, 9), complex)
    S[..., 0:3, 0:3] = Wd[..., 3:6, :]
    S[..., 0:3, 3:6] = Wu[..., 3:6, :] * Eu[..., None, :]
    S[..., 3:9, 0:3] = -Wd * Ed[..., None, :]
    S[..., 3:9, 3:6] = -Wu
    S[..., 3:9, 6:9] = Wd
    rhs = np.zeros(n + (9,), complex)
    rhs[..., 3:9] = s
    x = np.linalg.solve(S, rhs[..., None])[..., 0]
    cd, cu = x[..., 0:3], x[..., 3:6]
    return np.einsum('...ij,...j->...i', Wd[..., 0:3, :], cd) + np.einsum('...ij,...j->...i', Wu[..., 0:3, :], Eu * cu)


def _check_moment(M):
    M = np.asarray(M, complex)
    if np.any(np.abs(M[..., 0, 2]) + np.abs(M[..., 1, 2]) > 1e-12 * max(np.abs(M).max(), 1e-300)):
        raise ValueError('moments with xz or yz parts are not implemented')
    return M


def surface_response_k(kx, ky, M, depth, omega, lam, mu, rho):
    """u^ at the surface (x, y, down) for the wavevectors (kx, ky) (batched) and moments M [..., 3, 3] (z down or up:
    only xx, yy, zz, xy enter; xz and yz must vanish). lam and mu may be complex (damped)."""
    M = _check_moment(M)
    return surface_from_modes(modes(kx, ky, omega, lam, mu, rho), depth, _jump(kx, ky, M, lam, mu))


def k_grid(omega, vp, vs, vr, depth, r_max, Q, k_top=40.0, n_base=6000, n_fine=600):
    """Wavenumbers for the radial integrals: a base grid to k_top / depth fine enough for J_n(k r_max), refined about the
    Rayleigh pole and the two branch points, where the damped integrand varies over a width k / (2 Q)."""
    kmax = k_top / depth
    dk = min(kmax / n_base, 0.25 / max(r_max, depth))
    ks = [np.arange(dk / 2, kmax, dk)]
    for c in (vp, vs, vr):
        k0 = omega / c
        w = k0 / (2 * Q)
        for span, n in ((60 * w, n_fine), (8 * w, n_fine)):
            ks.append(np.linspace(max(k0 - span, dk / 4), k0 + span, n))
    return np.unique(np.concatenate(ks))


def radial_coefficients(M, depth, omega, vp, vs, rho, Q, ks, taper_m=None):
    """c_n(k) [orders, 3, k] for a moment M (x, y, up frame) at `depth`, or for a column of sources: M [n, 3, 3] at
    depths [n], summed (one modal decomposition serves every depth). `taper_m` multiplies by exp(-(k taper_m)^2), a
    Gaussian smoothing of the surface field over about taper_m that lets shallow sources be integrated on a bounded
    grid (used only where the field is read at distances much larger than taper_m)."""
    vpc, vsc = vp * (1 - 0.5j / Q), vs * (1 - 0.5j / Q)
    mu = rho * vsc ** 2
    lam = rho * vpc ** 2 - 2 * mu
    phis = 2 * np.pi * np.arange(N_PHI) / N_PHI
    Ms = np.asarray(M, float).reshape(-1, 3, 3)
    depths = np.atleast_1d(np.asarray(depth, float))
    if len(depths) != len(Ms):
        raise ValueError('one depth per moment')
    out = np.zeros((N_PHI, 3, len(ks)), complex)
    ky = np.zeros_like(ks)
    md = modes(ks, ky, omega, lam, mu, rho)                    # a wavevector along x; the moments rotate instead
    for m, p in enumerate(phis):
        c, s = np.cos(p), np.sin(p)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
        u = np.zeros((len(ks), 3), complex)
        for Mj, zj in zip(Ms, depths):
            Mp = _check_moment(R.T @ Mj @ R)                   # the moment in the frame whose x is along the wavevector
            u += surface_from_modes(md, zj, _jump(ks, ky, np.broadcast_to(Mp, ks.shape + (3, 3)), lam, mu))
        u = u @ R.T                                            # back to (x, y, down)
        u[:, 2] *= -1                                          # down to up
        out[m] = u.T
    if taper_m:
        out *= np.exp(-(ks * taper_m) ** 2)[None, None, :]
    # u^(k, phi) = sum_n c_n exp(i n phi)
    E = np.exp(-1j * np.outer(ORDERS, phis)) / N_PHI
    return np.einsum('nm,mck->nck', E, out)


def radial_profiles(cn, ks, rs):
    """F_n(r) = int c_n(k) J_n(k r) k dk [orders, 3, r] (trapezoid on the nonuniform grid); cn one set [orders, 3, k] or
    a list of sets on the same grid, which share the Bessel functions (a list back)."""
    many = isinstance(cn, (list, tuple))
    sets = list(cn) if many else [cn]
    wts = np.zeros_like(ks)
    d = np.diff(ks)
    wts[:-1] += d / 2
    wts[1:] += d / 2
    outs = [np.zeros((len(ORDERS), 3, len(rs)), complex) for _ in sets]
    for i, n in enumerate(ORDERS):
        J = jv(abs(n), np.outer(rs, ks)) * (1 if n >= 0 or n % 2 == 0 else -1)
        for c, o in zip(sets, outs):
            o[i] = np.einsum('ck,rk->cr', c[i] * (ks * wts)[None, :], J)
    return outs if many else outs[0]


def field_from_profiles(F, rs, x, y):
    """The surface field [3, ...] at (x, y) from radial profiles F [orders, 3, r] on rs (linear in r between samples)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    out = np.zeros((3,) + x.shape, complex)
    for i, n in enumerate(ORDERS):
        ph = (1j ** n) * np.exp(1j * n * th) / (2 * np.pi)
        for c in range(3):
            Fr = np.interp(r, rs, F[i, c].real) + 1j * np.interp(r, rs, F[i, c].imag)
            out[c] += ph * Fr
    return out


def surface_displacement(M, depth, x, y, omega, vp, vs, rho, Q=100.0, ks=None, r_max=None):
    """Surface displacement [3, ...] (x, y, up), complex amplitude at frequency omega (exp(-i omega t)), of a point moment
    M (x, y, up frame; no xz or yz part) at `depth` below (0, 0), in a half-space damped with quality factor Q."""
    from .analytic import rayleigh_speed
    x, y = np.asarray(x, float), np.asarray(y, float)
    rmax = float(np.hypot(x, y).max()) if r_max is None else r_max
    if ks is None:
        ks = k_grid(omega, vp, vs, rayleigh_speed(vp, vs), depth, rmax, Q)
    cn = radial_coefficients(M, depth, omega, vp, vs, rho, Q, ks)
    rs = np.linspace(0, rmax * 1.0001, max(400, int(rmax / max(depth, 1) * 200)))
    return field_from_profiles(radial_profiles(cn, ks, rs), rs, x, y)
