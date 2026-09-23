"""Orbit geometry for dwell products: WGS84, state-vector interpolation, virtual baselines.

The single-image "tomography" method treats the platform's own along-track motion during one
aperture as if it were a multi-baseline tomographic stack. This module computes that virtual
baseline the way the method needs it, but from the product's real state vectors rather than
from a scalar aperture-span approximation. Nothing here asserts that an along-track baseline
carries height or depth information; it reproduces the published geometry faithfully so that
the resulting depth axis can be tested.

Reference for the steering wavenumber: Biondi & Malanga 2022 eq. 22,
Kz = 4 pi B_perp / (lambda_s r sin theta), with lambda_s a declared "sound" wavelength.
"""
from __future__ import annotations

import numpy as np

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


def lla_to_ecef(latitude_deg, longitude_deg, height_m):
    """WGS84 geodetic to ECEF metres. Height is ellipsoidal, not orthometric."""
    lat = np.deg2rad(np.asarray(latitude_deg, float))
    lon = np.deg2rad(np.asarray(longitude_deg, float))
    h = np.asarray(height_m, float)
    if not (np.all(np.isfinite(lat)) and np.all(np.isfinite(lon)) and np.all(np.isfinite(h))):
        raise ValueError("coordinates must be finite")
    if np.any(np.abs(latitude_deg) > 90) or np.any(np.abs(longitude_deg) > 180):
        raise ValueError("coordinates must be WGS84 degrees")
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    return np.stack([(n + h) * np.cos(lat) * np.cos(lon),
                     (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1.0 - WGS84_E2) + h) * np.sin(lat)], axis=-1)


class Ephemeris:
    """Cubic-interpolated state vectors, evaluated in seconds from a declared epoch.

    Times must be strictly increasing and the requested epoch must lie inside the supplied
    span; extrapolation is refused rather than silently continued.
    """

    def __init__(self, times_s, positions_m, velocities_m_s):
        self.t = np.asarray(times_s, float)
        self.p = np.asarray(positions_m, float)
        self.v = np.asarray(velocities_m_s, float)
        if self.t.ndim != 1 or self.t.size < 4:
            raise ValueError("need at least four state vectors")
        if self.p.shape != (self.t.size, 3) or self.v.shape != (self.t.size, 3):
            raise ValueError("positions and velocities must be [N, 3] matching the times")
        if not (np.all(np.isfinite(self.t)) and np.all(np.isfinite(self.p)) and np.all(np.isfinite(self.v))):
            raise ValueError("state vectors must be finite")
        if np.any(np.diff(self.t) <= 0):
            raise ValueError("state-vector times must be strictly increasing")

    @property
    def span(self):
        return float(self.t[0]), float(self.t[-1])

    def at(self, times):
        """Position and velocity at the requested times, by component-wise cubic interpolation."""
        q = np.atleast_1d(np.asarray(times, float))
        if not np.all(np.isfinite(q)):
            raise ValueError("requested times must be finite")
        if q.min() < self.t[0] or q.max() > self.t[-1]:
            raise ValueError(f"requested time outside state-vector span {self.span}")
        from scipy.interpolate import CubicSpline
        pos = np.stack([CubicSpline(self.t, self.p[:, i])(q) for i in range(3)], axis=-1)
        vel = np.stack([CubicSpline(self.t, self.v[:, i])(q) for i in range(3)], axis=-1)
        return pos, vel


def line_of_sight(sensor_pos, target_pos):
    """Unit line-of-sight vector from sensor to target, and the slant range."""
    d = np.asarray(target_pos, float) - np.asarray(sensor_pos, float)
    r = np.linalg.norm(d, axis=-1)
    if np.any(r <= 0):
        raise ValueError("sensor and target coincide")
    return d / r[..., None], r


def perpendicular_baselines(ephemeris, times_s, target_ecef, reference_time_s):
    """Along-track virtual baselines projected perpendicular to the line of sight.

    For each sub-aperture time the platform offset from the aperture-centre state is projected
    onto the direction that is along-track and orthogonal to the reference line of sight. This
    is the quantity the published method substitutes for a cross-track tomographic baseline.
    Returns the baselines, the reference slant range, and the incidence angle in degrees.
    """
    target = np.asarray(target_ecef, float)
    if target.shape != (3,):
        raise ValueError("target must be a single ECEF position")
    ref_pos, ref_vel = ephemeris.at(reference_time_s)
    ref_pos, ref_vel = ref_pos[0], ref_vel[0]
    los, slant_range = line_of_sight(ref_pos, target)
    along = ref_vel - np.dot(ref_vel, los) * los
    norm = np.linalg.norm(along)
    if norm <= 0:
        raise ValueError("velocity is parallel to the line of sight; no along-track direction")
    along = along / norm
    pos, _ = ephemeris.at(times_s)
    baselines = (pos - ref_pos) @ along
    # Incidence at the target: angle between the up direction at the target and the target-to-sensor ray.
    up = target / np.linalg.norm(target)
    incidence = np.rad2deg(np.arccos(np.clip(np.dot(up, -los), -1.0, 1.0)))
    return baselines, float(slant_range), float(incidence)


def steering_wavenumbers(baselines_m, slant_range_m, incidence_deg, lam_s_m):
    """Kz = 4 pi B_perp / (lambda_s r sin theta), the published depth-steering wavenumber."""
    b = np.asarray(baselines_m, float)
    if b.ndim != 1 or b.size < 2 or not np.all(np.isfinite(b)):
        raise ValueError("need at least two finite baselines")
    if not np.isfinite([slant_range_m, incidence_deg, lam_s_m]).all():
        raise ValueError("geometry must be finite")
    if slant_range_m <= 0 or lam_s_m <= 0 or not 0 < incidence_deg < 90:
        raise ValueError("invalid slant range, wavelength model, or incidence angle")
    return 4.0 * np.pi * b / (lam_s_m * slant_range_m * np.sin(np.deg2rad(incidence_deg)))
