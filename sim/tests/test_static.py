"""P2-04: static relaxation with the elastic solver, and Eshelby's void with Okada's point sources."""
import copy

import numpy as np
import pytest

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.seismic.analytic import moment_surface_displacement, void_moment
from katabasis.seismic.arrays import snap_to_ground
from katabasis.seismic.elastic3d import Medium
from katabasis.seismic.static import relax


def test_three_equal_tensile_sources_are_mogis_source():
    lam, mu, d = 1.3, 1.0, 12.0
    nu = lam / (2 * (lam + mu))
    X, Y = np.meshgrid(np.linspace(-30, 30, 7), [0.0, 5.0], indexing='ij')
    u = moment_surface_displacement(np.eye(3) * (3 * lam + 2 * mu), d, X, Y, lam, mu)
    dV = (3 * lam + 2 * mu) / (lam + 2 * mu)             # the cavity volume change of this moment (Mogi)
    R3 = (X ** 2 + Y ** 2 + d ** 2) ** 1.5
    assert u == pytest.approx((1 - nu) * dV / np.pi * np.stack([X / R3, Y / R3, d / R3]), rel=1e-10)


def test_a_void_carries_about_twice_the_stress_it_removes():
    M = void_moment(np.eye(3), 1.0, 0.25)
    assert np.diag(M) == pytest.approx([3 * (1 - 0.25) / (2 * (1 - 2 * 0.25))] * 3)      # 1 / (1 - a) = 2.25
    D = void_moment(np.diag([1.0, -1.0, 0.0]), 1.0, 0.25)
    assert D[0, 0] == pytest.approx(15 * (1 - 0.25) / (7 - 5 * 0.25))                  # 1 / (1 - b)


@pytest.fixture(scope='module')
def bench():
    site = load_site('bench-void')
    g = Grid.covering((-36.0, 36.0), (-36.0, 36.0), (-44.0, 6.0), 2.0)      # three air cells above the ground
    full = Medium.from_model(voxelise(site, g, heterogeneity=False))
    raw = copy.deepcopy(site.raw)
    raw['features'] = []
    bg = Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))
    rec = snap_to_ground(g, full.solid, np.array([[0.0, 0.0], [20.0, 0.0]]))
    return g, full, bg, rec


def test_intact_rock_under_horizontal_stress_does_not_move(bench):
    g, full, bg, rec = bench
    u, _, _, _ = relax(bg, {'xx': 1e6, 'yy': 5e5, 'xy': 2e5}, rec, 0.2, pml_width=6)
    assert np.abs(u).max() == 0.0


def test_the_surface_sinks_over_a_chamber_in_horizontal_tension(bench):
    g, full, bg, rec = bench
    u, hist, _, _ = relax(full, {'xx': 1e6, 'yy': 1e6}, rec, 0.4, pml_width=6)
    lam, mu = full.lam.max(), full.mu_xy.max()
    nu = lam / (2 * (lam + mu))
    ua = moment_surface_displacement(void_moment(np.diag([1e6, 1e6, 0.0]), 216.0, nu), 15.0,
                                     rec[:, 0], rec[:, 1], lam, mu)
    assert u[0, 2] < 0 < u[1, 2]                          # subsidence above, a ring of uplift beyond
    assert u[0, 2] == pytest.approx(ua[2, 0], rel=0.5)    # a 6 m cube on a 2 m grid against a point sphere
    assert abs(hist['u'][0, 2, -1] - hist['u'][0, 2, -2]) < 0.01 * abs(u[0, 2])
