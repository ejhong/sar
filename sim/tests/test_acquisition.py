"""P2-01: the real dwell as the simulator's geometry, from the committed acquisition records."""
import numpy as np
import pytest

from sarsim import point_targets, synthesize
from sarsim.acquisition import DwellGeometry, RECORDS

NAMES = ['giza-20250827', 'sacsayhuaman-20250822']


@pytest.fixture(params=NAMES)
def g(request):
    return DwellGeometry.from_record(request.param)


def test_records_round_trip(g):
    h = DwellGeometry.from_record(g.record())
    assert h == g


def test_three_speeds_obey_the_geometric_mean(g):
    # the Doppler rate's effective speed is the geometric mean of the platform and ground-sweep speeds
    assert g.V_eff ** 2 / (g.V_platform * g.V) == pytest.approx(1, abs=1e-4)
    assert g.V < g.V_eff < g.V_platform


def test_giza_numbers_match_the_product():
    g = DwellGeometry.from_record('giza-20250827')
    assert g.aperture_time == pytest.approx(24.49, abs=0.01)
    assert g.lam == pytest.approx(0.031228, abs=1e-6)
    assert np.rad2deg(np.arccos(g.los_enu[2])) == pytest.approx(g.theta_deg, abs=0.01)
    assert g.source['window_azimuth'] == 'NONE' and g.source['window_range'] == 'NONE'


def test_slow_time_follows_the_product_convention(g):
    # f = Ka (t - t_zd) with Ka < 0: positive Doppler is earlier
    assert g.Ka_signed < 0
    assert g.nu_to_time(g.nu_band / 4) < 0 < g.nu_to_time(-g.nu_band / 4)
    assert g.nu_to_time(g.nu_band / 2) - g.nu_to_time(-g.nu_band / 2) == pytest.approx(-g.aperture_time)


def test_resolution_times_duration_is_fixed(g):
    for d in (0.3, 1.0, 3.0):
        assert g.look_resolution(d) * d == pytest.approx(0.886 * g.V / g.Ka)
        assert g.look_duration(g.look_resolution(d)) == pytest.approx(d)


def test_a_receding_target_shifts_as_predicted(g):
    one = lambda: point_targets(np.zeros(1), np.zeros(1), np.zeros(1), np.ones(1))
    shape = (256, 64)
    v = 1.5e-3
    a = synthesize(one(), g, shape)
    b = synthesize(one(), g, shape, motion=lambda x, y, z, t: v * t[:, None] * np.ones(len(x)))
    # the shift is a phase ramp across the azimuth spectrum; read it from the cross-spectrum's slope
    j = np.argmax(np.abs(a).max(axis=0))
    A, B = np.fft.fft(a[:, j]), np.fft.fft(b[:, j])
    nu = np.fft.fftfreq(shape[0], d=g.dx)
    inb = np.abs(nu) <= g.nu_band / 2.2
    ph = np.unwrap(np.angle(B[inb] * np.conj(A[inb]))[np.argsort(nu[inb])])
    slope = np.polyfit(np.sort(nu[inb]), ph, 1)[0]
    shift = -slope / (2 * np.pi)
    assert shift == pytest.approx(2 * g.V * v / (g.lam * g.Ka_signed), rel=1e-3)
    assert g.shift_to_velocity(shift) == pytest.approx(v, rel=1e-3)
    assert abs(g.shift_to_velocity(1.0)) == pytest.approx(g.V_platform / g.R0, rel=1e-4)


def test_absolute_time_includes_each_row_zero_doppler_time():
    g = DwellGeometry.from_record('giza-20250827')
    seen = {}

    def field(east, north, t_abs):
        seen['t'] = t_abs
        return np.zeros(t_abs.shape + (3,))
    motion = g.absolute_motion(field)
    xr = np.array([0.0, 70.136])                 # 70 m further along track: 10 ms later
    motion(xr, np.zeros(2), np.zeros(2), np.array([0.0, 1.0]))
    assert seen['t'][:, 1] - seen['t'][:, 0] == pytest.approx(70.136 / g.V)
    assert seen['t'][1, 0] - seen['t'][0, 0] == pytest.approx(1.0)


def test_site_frame_round_trip_and_line_of_sight():
    g = DwellGeometry.from_record('giza-20250827')
    a, r = g.along_track_en, g.ground_range_en
    assert abs(a @ r) < 0.02                       # along track and ground range are nearly perpendicular
    up = np.array([0.0, 0.0, 1e-3])
    assert g.slant_range_increase(up) == pytest.approx(-1e-3 * np.cos(np.deg2rad(g.theta_deg)), rel=1e-3)
    assert (RECORDS / 'giza-20250827.json').is_file()
