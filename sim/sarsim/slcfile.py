"""Synthetic single-look complex products that read like a real ICEYE file (Phase 2, benchmark).

A simulated image is only useful to a reconstruction written for real products if that reconstruction can read it
unchanged. `write_synthetic_slc` writes an HDF5 file with the datasets the ICEYE SLC format uses for what such code
reads (the complex image as s_i, s_q; the azimuth time interval; the zero-Doppler start time; state vectors in ECEF;
the range sampling rate; the Doppler-rate polynomial), marks it as synthetic in its attributes and in a `synthetic`
group recording the scene, the motion and the simulator, and keeps timing, orbit and Doppler consistent with the
pulse-level simulation that made the image:

* the orbit is the simulation's equivalent straight line, flown at V_eff past the scene centre at the product's
  slant range and incidence, so the Doppler rate the geometry implies is the one written;
* row k is zero-Doppler time t0 + k row_dt, the scene centre's row at the pass centre;
* a column's fast time is its slant range over c/2, and the range sampling rate is c / (2 dr).

Nothing in the file is real imagery; its name and attributes say so.
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .echo import C
from .orbit import lla_to_ecef


def _enu_basis(lat_deg, lon_deg):
    la, lo = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    E = np.array([-np.sin(lo), np.cos(lo), 0.0])
    N = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    U = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    return E, N, U


def straight_orbit(setup, lat_deg, lon_deg, height_m, heading_deg, theta_deg, look='right', times=None):
    """State vectors (ECEF) of the equivalent straight-line pass: speed V_eff along the heading, closest approach to the
    scene centre at slant range R0 and incidence theta at t = 0 (flat local geometry about the scene centre)."""
    E, N, U = _enu_basis(lat_deg, lon_deg)
    h = np.deg2rad(heading_deg)
    along = np.sin(h) * E + np.cos(h) * N
    right = np.cos(h) * E - np.sin(h) * N                   # 90 degrees clockwise from the heading, in the plane
    side = right if look == 'right' else -right
    th = np.deg2rad(theta_deg)
    target = lla_to_ecef(lat_deg, lon_deg, height_m)
    # the satellite sits up and to the side opposite the look: the line of sight points down and toward `side`
    sat0 = target + setup.R0 * (np.cos(th) * U - np.sin(th) * side)
    t = np.linspace(-setup.dwell_s / 2 - 2.0, setup.dwell_s / 2 + 2.0, 41) if times is None else np.asarray(times, float)
    pos = sat0[None, :] + setup.V_eff * t[:, None] * along[None, :]
    vel = np.repeat((setup.V_eff * along)[None, :], len(t), axis=0)
    return t, pos, vel, target


def write_synthetic_slc(path, image, setup, Ka_signed, lat_deg, lon_deg, height_m, heading_deg, theta_deg,
                        epoch=datetime(2030, 1, 1, tzinfo=timezone.utc), scene: dict | None = None, look='right'):
    """Write `image` (rows x cols complex, rows the product's zero-Doppler rows, scene centre at the middle row and
    column) as an ICEYE-like SLC. Returns the file's path."""
    import h5py
    path = Path(path)
    rows, cols = image.shape
    t_sv, pos, vel, target = straight_orbit(setup, lat_deg, lon_deg, height_m, heading_deg, theta_deg, look)
    iso = lambda s: (epoch + timedelta(seconds=float(s))).strftime('%Y-%m-%dT%H:%M:%S.%f')
    zd_start = -(rows // 2) * setup.row_dt
    near_range = setup.R0 - (cols // 2) * setup.dr
    with h5py.File(path, 'w') as f:
        f.attrs['synthetic'] = 'yes: simulated imagery, not a satellite acquisition'
        f.attrs['simulator'] = 'sarsim.echo (pulse-level, time-domain back-projection)'
        f['s_i'] = image.real.astype(np.float32)
        f['s_q'] = image.imag.astype(np.float32)
        f['azimuth_time_interval'] = float(setup.row_dt)
        f['zerodoppler_start_utc'] = iso(zd_start)
        f['zerodoppler_end_utc'] = iso(zd_start + (rows - 1) * setup.row_dt)
        f['state_vector_time_utc'] = np.array([iso(s) for s in t_sv], dtype='S32')
        for k, axis in enumerate('XYZ'):
            f[f'pos{axis}'] = pos[:, k]
            f[f'vel{axis}'] = vel[:, k]
        f['range_sampling_rate'] = float(C / (2 * setup.dr))
        f['slant_range_to_first_pixel'] = float(near_range)
        f['doppler_rate_coeffs'] = np.array([float(Ka_signed)])
        f['carrier_frequency'] = float(C / setup.lam)
        f['chirp_bandwidth'] = float(setup.bandwidth_hz)
        f['heading'] = float(heading_deg)
        f['look_side'] = look
        f['coord_center'] = np.array([cols // 2 + 1, rows // 2 + 1, lat_deg, lon_deg], float)
        g = f.create_group('synthetic')
        g.attrs['scene'] = __import__('json').dumps(scene or {})
        g['scene_centre_ecef'] = target
        g['V_eff'] = setup.V_eff
        g['prf_simulated'] = setup.prf
        g['dwell_s'] = setup.dwell_s
    return path
