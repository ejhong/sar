"""Pulse-level SAR simulation: echoes pulse by pulse, focused by time-domain back-projection (Phase 2, P2-18).

The lab's fast synthesizer (`synthesize`) writes a focused image directly in the azimuth-frequency domain, mapping
each Doppler frequency to a moment of the pass (t = f / Ka) and applying motion there as a phase. That is the
stationary-phase approximation, efficient and usually accurate, but it is not an independent calculation of what a
radar records. This module is one: it computes, for every transmitted pulse, the range-compressed echo of every
scatterer from its exact distance to the platform at that moment (the scatterer displaced by its own motion), and
forms the image by back-projection, summing every pulse's echo at every pixel's own range history. Range migration,
the exact phase history and motion in both directions of the slant plane come out of the calculation rather than
being assumed.

Geometry: the equivalent straight-line model of a spaceborne pass. In the slant plane the platform moves along u at
V_eff = sqrt(V_s V_g), which reproduces the product's Doppler rate Ka = -2 V_eff^2 / (lambda R0); a scatterer at
(u, r) has its zero-Doppler time at t0 = u / V_eff and closest range r, and its image row advances at V_g per second of
t0, so the pixel grid is the product's: rows row_dt apart in zero-Doppler time, columns dr apart in slant range.
Scatterer motion is a displacement (du along track, dr along the line of sight at closest approach) as a function of
absolute time. The range-compressed pulse is the ideal response of the chirp's bandwidth, sinc((R' - R)/rho) with
rho = c / 2B, stored in a window that follows the scene centre's range (a target's range changes by kilometres over a
dwell, a small scene's by tens of metres). Pulses only have to sample the scene's differential Doppler and the motion:
at PRF p the azimuth ambiguities fall p lambda R / (2 V_eff) (V_g / V_eff) metres away, far outside a small scene.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

C = 299_792_458.0


@dataclass
class EchoSetup:
    lam: float                 # wavelength, m
    R0: float                  # slant range of the scene centre at closest approach, m
    V_eff: float               # equivalent platform speed in the slant plane, m/s
    V_g: float                 # ground-sweep speed of the image rows, m/s
    bandwidth_hz: float        # chirp bandwidth
    dwell_s: float             # duration of the pass
    prf: float                 # simulated pulse rate (only needs to sample the scene's differential Doppler)
    row_dt: float              # image row spacing in zero-Doppler time, s
    dr: float                  # image column spacing in slant range, m
    oversample: int = 4        # range samples per image column in the stored echoes

    @classmethod
    def from_geometry(cls, g, prf: float = 1000.0, dwell_s: float | None = None, oversample: int = 4):
        """From a DwellGeometry: V_eff from the product's own Doppler rate, the chirp bandwidth from its range
        band, the dwell from its processed Doppler band."""
        V_eff = float(np.sqrt(abs(g.Ka_signed) * g.lam * g.R0 / 2))
        bw = C * g.kr_band / 2                                   # kr_band cycles per metre of slant range
        dwell = g.doppler_band_hz / abs(g.Ka_signed) if dwell_s is None else dwell_s
        return cls(lam=g.lam, R0=g.R0, V_eff=V_eff, V_g=g.V, bandwidth_hz=bw, dwell_s=dwell, prf=prf,
                   row_dt=g.row_dt, dr=g.dr, oversample=oversample)

    @property
    def rho(self) -> float:
        return C / (2 * self.bandwidth_hz)

    @property
    def times(self) -> np.ndarray:
        n = int(round(self.dwell_s * self.prf))
        return (np.arange(n) - (n - 1) / 2) / self.prf


@dataclass
class PointScene:
    """Scatterers in the slant plane: u along track (m, zero-Doppler time u / V_eff), r slant range offset from the
    scene centre (m), complex amplitude a; each may vibrate: displacement du(t) = au sin(2 pi f t + phi) along track and
    dr(t) = ar sin(2 pi f t + phi) along the line of sight (increase in range), t absolute time from the pass centre."""
    u: np.ndarray
    r: np.ndarray
    a: np.ndarray
    au: np.ndarray
    ar: np.ndarray
    f: np.ndarray
    phi: np.ndarray

    @classmethod
    def still(cls, u, r, a):
        u, r, a = np.asarray(u, float), np.asarray(r, float), np.asarray(a, complex)
        z = np.zeros_like(u)
        return cls(u=u, r=r, a=a, au=z.copy(), ar=z.copy(), f=z.copy(), phi=z.copy())

    def moving(self, select, ar=0.0, au=0.0, f=0.0, phi=0.0):
        """A copy in which the selected scatterers vibrate."""
        s = PointScene(**{k: np.array(v, copy=True) for k, v in self.__dict__.items()})
        s.ar[select], s.au[select], s.f[select], s.phi[select] = ar, au, f, phi
        return s


@njit(parallel=True, fastmath=False, cache=True)
def _echoes(t, V_eff, R0, u, r, a_re, a_im, au, ar, f, phi, lam, rho, r_lo, dr_s, n_s, half):
    """Range-compressed echoes [pulses, samples] in a window starting r_lo[n] metres at pulse n."""
    n_p = t.shape[0]
    out = np.zeros((n_p, n_s), np.complex128)
    k4 = 4.0 * np.pi / lam
    for n in prange(n_p):
        tn = t[n]
        x_p = V_eff * tn
        for k in range(u.shape[0]):
            s = np.sin(2.0 * np.pi * f[k] * tn + phi[k])
            du = au[k] * s
            dR = ar[k] * s
            dx = x_p - (u[k] + du)
            rr = R0 + r[k] + dR
            R = np.sqrt(rr * rr + dx * dx)
            ph = -k4 * R
            c_re = a_re[k] * np.cos(ph) - a_im[k] * np.sin(ph)
            c_im = a_re[k] * np.sin(ph) + a_im[k] * np.cos(ph)
            pos = (R - r_lo[n]) / dr_s
            i0 = int(np.floor(pos))
            for i in range(i0 - half, i0 + half + 1):
                if i < 0 or i >= n_s:
                    continue
                x = (i - pos) * dr_s / rho
                w = 1.0 if abs(x) < 1e-12 else np.sin(np.pi * x) / (np.pi * x)
                out[n, i] += complex(c_re * w, c_im * w)
    return out


@njit(parallel=True, fastmath=False, cache=True)
def _backproject(echo, t, V_eff, R0, rows_u, cols_r, lam, r_lo, dr_s):
    """Sum every pulse at every pixel's own range history, I = sum_n s_n(R_n) exp(+j 4 pi (R_n - r0) / lambda), r0 the
    pixel's closest range: a baseband image in which a scatterer carries the phase -4 pi r / lambda of its own range, as
    in a product's SLC (without the -r0 term the carrier would ride on the range response and alias on the grid)."""
    n_p, n_s = echo.shape
    ny, nx = rows_u.shape[0], cols_r.shape[0]
    img = np.zeros((ny, nx), np.complex128)
    k4 = 4.0 * np.pi / lam
    for iy in prange(ny):
        uk = rows_u[iy]
        for ix in range(nx):
            rr = R0 + cols_r[ix]
            acc_re = 0.0
            acc_im = 0.0
            for n in range(n_p):
                dx = V_eff * t[n] - uk
                R = np.sqrt(rr * rr + dx * dx)
                pos = (R - r_lo[n]) / dr_s
                i0 = int(np.floor(pos))
                if i0 < 0 or i0 + 1 >= n_s:
                    continue
                w = pos - i0
                s = echo[n, i0] * (1.0 - w) + echo[n, i0 + 1] * w
                ph = k4 * (R - rr)
                cr, ci = np.cos(ph), np.sin(ph)
                acc_re += s.real * cr - s.imag * ci
                acc_im += s.real * ci + s.imag * cr
            img[iy, ix] = complex(acc_re, acc_im)
    return img


def window(setup: EchoSetup, t: np.ndarray, margin_m: float, u_span: float, r_span: float):
    """The stored range window at each pulse: centred on the scene centre's range, wide enough for the scene's
    differential migration (along-track extent times the look's squint) and its range extent."""
    x = setup.V_eff * t
    R_c = np.sqrt(setup.R0 ** 2 + x ** 2)
    squint = np.abs(x) / R_c
    half = r_span / 2 + u_span / 2 * squint.max() + margin_m
    dr_s = setup.dr / setup.oversample
    n_s = int(np.ceil(2 * half / dr_s)) + 1
    return R_c - half, dr_s, n_s


def simulate(setup: EchoSetup, scene: PointScene, rows: int, cols: int, margin_m: float = 8.0, sinc_lobes: float = 4.0,
             pulses: tuple[float, float] | None = None):
    """Echoes of `scene` over the pass, back-projected onto a rows x cols image centred on the scene centre, rows
    row_dt apart in zero-Doppler time and columns dr apart in slant range, like the product's own grid. With `pulses`
    (t_lo, t_hi), only the pulses sent in that stretch of the pass are used."""
    t = setup.times
    if pulses is not None:
        t = t[(t >= pulses[0]) & (t < pulses[1])]
    rows_t0 = (np.arange(rows) - rows // 2) * setup.row_dt
    rows_u = setup.V_eff * rows_t0
    cols_r = (np.arange(cols) - cols // 2) * setup.dr
    u_span = max(np.ptp(rows_u), np.ptp(scene.u) if scene.u.size else 0.0)
    r_span = max(np.ptp(cols_r), np.ptp(scene.r) if scene.r.size else 0.0)
    r_lo, dr_s, n_s = window(setup, t, margin_m, u_span, r_span)
    half_taps = int(np.ceil(sinc_lobes * setup.rho / dr_s))           # the range sinc kept to +-sinc_lobes widths
    e = _echoes(t, setup.V_eff, setup.R0, scene.u, scene.r, scene.a.real.copy(), scene.a.imag.copy(), scene.au,
                scene.ar, scene.f, scene.phi, setup.lam, setup.rho, r_lo, dr_s, n_s, half_taps)
    img = _backproject(e, t, setup.V_eff, setup.R0, rows_u, cols_r, setup.lam, r_lo, dr_s)
    return img / len(t)
