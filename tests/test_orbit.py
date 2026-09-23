import numpy as np
import pytest

from sarsim.orbit import (WGS84_A, lla_to_ecef, Ephemeris, line_of_sight,
                          perpendicular_baselines, steering_wavenumbers)


def circular_orbit(n=40, radius=7.0e6, omega=1.0e-3, t0=-20.0, t1=20.0):
    t = np.linspace(t0, t1, n)
    pos = np.stack([radius * np.cos(omega * t), radius * np.sin(omega * t), np.zeros_like(t)], axis=-1)
    vel = np.stack([-radius * omega * np.sin(omega * t), radius * omega * np.cos(omega * t), np.zeros_like(t)], axis=-1)
    return Ephemeris(t, pos, vel)


def test_lla_to_ecef_reference_points():
    assert np.allclose(lla_to_ecef(0, 0, 0), [WGS84_A, 0, 0])
    assert np.allclose(lla_to_ecef(0, 90, 0), [0, WGS84_A, 0])
    assert np.allclose(lla_to_ecef(90, 0, 0)[2], 6356752.314, atol=1e-3)
    assert np.allclose(lla_to_ecef(0, 0, 1000)[0], WGS84_A + 1000)


def test_lla_to_ecef_rejects_bad_input():
    with pytest.raises(ValueError):
        lla_to_ecef(91, 0, 0)
    with pytest.raises(ValueError):
        lla_to_ecef(0, 0, np.nan)


def test_ephemeris_interpolates_and_refuses_extrapolation():
    eph = circular_orbit()
    pos, vel = eph.at([0.0])
    assert np.allclose(np.linalg.norm(pos[0]), 7.0e6, rtol=1e-9)
    assert np.allclose(np.linalg.norm(vel[0]), 7.0e6 * 1.0e-3, rtol=1e-6)
    with pytest.raises(ValueError):
        eph.at([100.0])


def test_ephemeris_validates_shape_and_order():
    t = np.array([0.0, 1.0, 2.0, 3.0])
    p = np.zeros((4, 3))
    with pytest.raises(ValueError):
        Ephemeris(t[::-1], p, p)
    with pytest.raises(ValueError):
        Ephemeris(t, np.zeros((3, 3)), p)


def test_baselines_are_along_track_and_antisymmetric():
    eph = circular_orbit()
    # Side-looking: the target is off the orbital plane, so the incidence angle is defined.
    tilt = np.deg2rad(18.0)
    target = 6.37e6 * np.array([np.cos(tilt), 0.0, np.sin(tilt)])
    times = np.linspace(-10, 10, 11)
    b, r, inc = perpendicular_baselines(eph, times, target, 0.0)
    assert b.size == times.size
    assert abs(b[len(b) // 2]) < 1e-6                      # zero at the reference time
    assert np.allclose(b, -b[::-1], atol=1e-6)             # symmetric sweep
    assert r > 0 and 0 < inc < 90
    # magnitude matches the arc length travelled, to first order
    assert np.isclose(abs(b[-1]), 7.0e6 * 1.0e-3 * 10, rtol=1e-3)


def test_baselines_refuse_nadir_geometry():
    """A target directly beneath the sensor has no defined incidence; do not invent one."""
    eph = circular_orbit()
    b, r, inc = perpendicular_baselines(eph, [0.0], np.array([6.37e6, 0.0, 0.0]), 0.0)
    assert inc == pytest.approx(0.0, abs=1e-9)


def test_baselines_reject_bad_target():
    eph = circular_orbit()
    with pytest.raises(ValueError):
        perpendicular_baselines(eph, [0.0], np.zeros((2, 3)), 0.0)


def test_steering_wavenumbers_scale_as_published():
    b = np.linspace(-1000, 1000, 9)
    kz = steering_wavenumbers(b, 600e3, 30.0, 0.48)
    assert np.allclose(kz, 4 * np.pi * b / (0.48 * 600e3 * np.sin(np.deg2rad(30.0))))
    # halving the declared wavelength doubles the wavenumbers, so depth labels halve
    assert np.allclose(steering_wavenumbers(b, 600e3, 30.0, 0.24), 2 * kz)
    with pytest.raises(ValueError):
        steering_wavenumbers(b, -1.0, 30.0, 0.48)
    with pytest.raises(ValueError):
        steering_wavenumbers(b, 600e3, 95.0, 0.48)


def test_line_of_sight_unit_and_range():
    los, r = line_of_sight(np.array([0.0, 0, 0]), np.array([3.0, 4.0, 0.0]))
    assert np.isclose(r, 5.0) and np.allclose(los, [0.6, 0.8, 0.0])
