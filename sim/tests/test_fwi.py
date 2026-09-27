"""Full-waveform inversion: the adjoint-state gradient against finite differences of the misfit."""
import numpy as np
import pytest

from katabasis.compose import Grid
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, ricker
from katabasis.seismic.fwi import Problem, gaussian_sigma


@pytest.fixture(scope='module')
def problem():
    g = Grid.covering((-10, 10), (-10, 10), (-10, 2), 0.5)
    X, Y, Z = np.meshgrid(g.x, g.y, g.z, indexing='ij')
    air = Z > 0
    rho = np.where(air, 0.0, 2400.0)
    vp0 = np.where(air, 0.0, 3000.0)
    vs0 = np.where(air, 0.0, 1700.0)
    blob = np.exp(-((X - 0.5) ** 2 + (Y + 0.5) ** 2 + (Z + 4) ** 2) / 2.0)
    vp_true, vs_true = vp0 * (1 - 0.15 * blob), vs0 * (1 - 0.15 * blob)
    zs = -0.25
    sources = np.array([[-3.0, 0.5, zs], [3.0, -1.0, zs]])
    receivers = np.array([[x, y, zs] for x in (-4.0, 0.0, 4.0) for y in (-4.0, 0.0, 4.0)])
    f0 = 150.0
    sim = Simulation(Medium.from_arrays(g, vp_true, vs_true, rho, air), pml_width=10, f0=f0)
    nt = int(0.02 / sim.dt)
    w = ricker(np.arange(nt) * sim.dt, f0) * 1e9
    obs = np.stack([sim.run([Source(tuple(s), w, 'force', (0, 0, -1))], Receivers(receivers), nt).traces
                    for s in sources]).astype(np.float64)
    ix = slice(int(round((-4 - g.origin[0]) / g.spacing)), int(round((4 - g.origin[0]) / g.spacing)) + 1)
    iz = slice(int(round((g.origin[2] + 1) / g.spacing)), int(round((g.origin[2] + 7) / g.spacing)) + 1)
    box = (ix, ix, iz)
    p = Problem(g, rho, air, vp0, vs0, box, sim.dt, sources, receivers, obs, w, f0, pml=10, every=2)
    Xb, Yb, Zb = X[box], Y[box], Z[box]
    bump = np.exp(-((Xb + 1) ** 2 + (Yb - 1) ** 2 + (Zb + 3.5) ** 2) / 1.5)
    return p, bump


@pytest.mark.parametrize('which', ['vp', 'vs'])
def test_adjoint_gradient_matches_finite_differences(problem, which):
    """Along the gradient itself (the direction the optimiser steps) for both speeds, and along a smooth
    bump for S. For records dominated by the surface wave the P-speed kernel is thousands of times
    weaker and nearly cancels over a bump; it is only held to its own direction."""
    p, bump = problem
    zero = np.zeros_like(bump)
    J, g_vp, g_vs, _ = p.evaluate(zero, zero)
    g = g_vp if which == 'vp' else g_vs
    for d in (g / np.abs(g).max(),) + ((bump,) if which == 'vs' else ()):
        eps = 0.01
        pert = lambda e: (e * d, zero) if which == 'vp' else (zero, e * d)
        fd = (p.evaluate(*pert(eps))[0] - p.evaluate(*pert(-eps))[0]) / (2 * eps)
        adj = float(np.sum(g * d))
        assert adj / fd == pytest.approx(1.0, abs=0.1), f'{which}: adjoint {adj:.4g} vs finite difference {fd:.4g}'


def test_the_low_passed_misfit_has_the_right_gradient(problem):
    """Synthetics and records filtered alike inside the misfit; the adjoint source filtered once more."""
    p0, bump = problem
    p = Problem(p0.grid, p0.rho, p0.air, p0.vp0, p0.vs0, p0.box, p0.dt, p0.sources, p0.receivers, p0.obs.copy(),
                p0.wavelet, p0.f0, pml=p0.pml, every=2, sigma_s=gaussian_sigma(90.0))
    zero = np.zeros_like(bump)
    _, _, g, _ = p.evaluate(zero, zero)
    for d in (g / np.abs(g).max(), bump):
        eps = 0.01
        fd = (p.evaluate(zero, eps * d)[0] - p.evaluate(zero, -eps * d)[0]) / (2 * eps)
        assert float(np.sum(g * d)) / fd == pytest.approx(1.0, abs=0.1)
