"""The half-space at one frequency reduces to Okada's point source at low frequency; Eshelby's cylinder matches Lamé."""
import unittest

import numpy as np

from katabasis.seismic.analytic import (moment_surface_displacement, rayleigh_speed, void_moment,
                                        void_moment_from_strain)
from katabasis.seismic.halfspace import k_grid, surface_displacement

VP, VS, RHO = 3300.0, 1690.0, 2400.0
MU = RHO * VS ** 2
LAM = RHO * VP ** 2 - 2 * MU


class HalfspaceTests(unittest.TestCase):
    def test_low_frequency_is_okada(self):
        depth = 100.0
        xs = np.linspace(-300, 300, 13)
        X, Y = np.meshgrid(xs, xs, indexing='ij')
        om = 0.02 * VS / depth                       # k_s depth = 0.02
        ks = k_grid(om, VP, VS, rayleigh_speed(VP, VS), depth, 430.0, 100.0, n_base=1500, n_fine=200)
        for M in (np.diag([1e9, 0, 0]), np.diag([0, 0, 1e9]), np.array([[0, 1e9, 0], [1e9, 0, 0], [0, 0, 0]])):
            ud = surface_displacement(M, depth, X, Y, om, VP, VS, RHO, Q=100.0, ks=ks)
            us = moment_surface_displacement(M, depth, X, Y, LAM, MU)
            self.assertLess(np.abs(ud.real - us).max() / np.abs(us).max(), 3e-3)
            # the damping's phase, about 1 / Q, and no more
            self.assertLess(np.abs(ud.imag).max() / np.abs(us).max(), 0.02)


class EshelbyTests(unittest.TestCase):
    def test_sphere_matches_void_moment(self):
        nu = LAM / (2 * (LAM + MU))
        e = np.diag([1e-6, 2e-6, -0.5e-6])
        C = LAM * np.trace(e) * np.eye(3) + 2 * MU * e
        np.testing.assert_allclose(void_moment_from_strain(e, 10.0, LAM, MU, 'sphere'), void_moment(C, 10.0, nu),
                                   rtol=1e-12)

    def test_cylinder_matches_lame_and_kirsch(self):
        nu = LAM / (2 * (LAM + MU))
        a, p = 5.0, 1e5
        # a hole under remote in-plane pressure p (plane strain): u_r = p a^2 / (2 mu r); a line of moment m per unit
        # length gives u_r = m (1 - 2 nu) / (4 pi mu (1 - nu) r)
        M = void_moment_from_strain(np.diag([p / (2 * (LAM + MU))] * 2 + [0.0]), np.pi * a * a, LAM, MU, 'cylinder')
        self.assertAlmostEqual(M[0, 0] * (1 - 2 * nu) / (4 * np.pi * MU * (1 - nu)), p * a * a / (2 * MU), delta=1e-12)
        # pure shear: the deviatoric moment 4 (1 - nu) pi a^2 p
        M = void_moment_from_strain(np.diag([p / (2 * MU), -p / (2 * MU), 0.0]), np.pi * a * a, LAM, MU, 'cylinder')
        self.assertAlmostEqual(M[0, 0] / (4 * (1 - nu) * np.pi * a * a * p), 1.0, places=12)


if __name__ == '__main__':
    unittest.main()
