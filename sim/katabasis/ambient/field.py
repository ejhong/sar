"""Ambient surface motion as a superposition of plane Rayleigh waves.

Microseisms are surface waves from distant ocean storms: at any site they
arrive as plane waves from many directions. A stationary random field with
a target vertical velocity spectrum P_z(f) is built from waves on a
frequency grid, several per bin, each with a random azimuth and phase and
an amplitude that makes the bin's power P_z(f) Δf. On a uniform half-space
each wave's surface motion is a retrograde ellipse whose horizontal to
vertical ratio follows from the P and S speeds (below).

Surface displacement of a Rayleigh wave, z down, e^{i(kx - ωt)}: with
q = sqrt(1 - c²/Vp²), s = sqrt(1 - c²/Vs²) and the traction-free surface,

    u_x / u_z = i (1 + s² - 2 q s) / (q (1 - s²))

so |H/V| = (1 + s² - 2qs) / (q (1 - s²)); 0.681 for a Poisson solid.
Waves this long (kilometres) sample far deeper than the site model; the
phase speed is a parameter, and the uniform half-space is a stated
simplification.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..seismic.analytic import rayleigh_speed


def rayleigh_hv(vp: float, vs: float) -> float:
    """Horizontal over vertical surface amplitude of a Rayleigh wave on a uniform half-space."""
    c = rayleigh_speed(vp, vs)
    q = np.sqrt(1 - (c / vp) ** 2)
    s = np.sqrt(1 - (c / vs) ** 2)
    return float((1 + s ** 2 - 2 * q * s) / (q * (1 - s ** 2)))


@dataclass
class PlaneWaveField:
    """Surface velocity from a set of plane Rayleigh waves (site axes: x east, y north, z up)."""
    freq: np.ndarray           # (n,) Hz
    amp: np.ndarray            # (n,) vertical velocity amplitude, m/s
    azimuth: np.ndarray        # (n,) propagation direction, radians clockwise from north
    phase: np.ndarray          # (n,) radians
    speed: np.ndarray          # (n,) phase speed, m/s
    hv: float                  # horizontal / vertical amplitude ratio

    def velocity(self, points: np.ndarray, t: np.ndarray) -> np.ndarray:
        """(npoints, 3, nt) velocity: east, north, up."""
        p = np.atleast_2d(points)[:, :2]
        kx = 2 * np.pi * self.freq / self.speed * np.sin(self.azimuth)       # wavenumber, east
        ky = 2 * np.pi * self.freq / self.speed * np.cos(self.azimuth)       # north
        out = np.zeros((len(p), 3, len(t)))
        w = 2 * np.pi * self.freq
        for n in range(len(self.freq)):
            arg = (p[:, 0:1] * kx[n] + p[:, 1:2] * ky[n]) - w[n] * t[None, :] + self.phase[n]
            vz = self.amp[n] * np.cos(arg)
            # retrograde ellipse: horizontal leads vertical by a quarter period along the propagation direction
            vh = self.hv * self.amp[n] * np.sin(arg)
            out[:, 0] += vh * np.sin(self.azimuth[n])
            out[:, 1] += vh * np.cos(self.azimuth[n])
            out[:, 2] += vz
        return out


def microseisms(freq: np.ndarray, psd_z: np.ndarray, rng: np.random.Generator, per_bin: int = 8,
                speed: float | np.ndarray = 3000.0, hv: float = 0.681) -> PlaneWaveField:
    """A field whose vertical velocity has the one-sided PSD psd_z ((m/s)²/Hz) on the frequency grid
    `freq` (uniformly spaced); `per_bin` waves per frequency bin, isotropic in azimuth. The default speed
    is sites/ambient.json regional.microseism_phase_speed (northern Egypt, from Adam et al. 2025)."""
    f = np.asarray(freq, float)
    df = np.gradient(f)
    F = np.repeat(f, per_bin) + (rng.random(len(f) * per_bin) - 0.5) * np.repeat(df, per_bin)
    A = np.repeat(np.sqrt(2 * psd_z * df / per_bin), per_bin)     # each: power A²/2 = P Δf / per_bin
    n = len(F)
    spd = np.broadcast_to(np.asarray(speed, float), (len(f),))
    return PlaneWaveField(F, A, rng.uniform(0, 2 * np.pi, n), rng.uniform(0, 2 * np.pi, n),
                          np.repeat(spd, per_bin), hv)
