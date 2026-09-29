"""Pulse-level simulation (sarsim.echo): focus, a Doppler slice as a stretch of the pass, and a steady velocity's shift."""
import numpy as np

from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup, PointScene, simulate
from sarsim.looks import look_masks, looks

G = DwellGeometry.from_record('giza-20220715')
ST = EchoSetup.from_geometry(G, prf=500.0)


def test_a_point_focuses_where_it_is_with_its_range_phase():
    u = 2.0 * ST.V_eff / ST.V_g
    img = simulate(ST, PointScene.still([u], [1.0], [1.0]), 256, 32)
    a = np.abs(img)
    iy, ix = np.unravel_index(a.argmax(), a.shape)
    assert (iy, ix) == (128 + round(2 / G.dx), 16 + round(1 / G.dr))
    # a scatterer on the pixel's own range carries -4 pi r / lambda (global R0 aside): compare two on-grid points
    one = simulate(ST, PointScene.still([0.0], [0.0], [1.0]), 64, 16)[32, 8]
    two = simulate(ST, PointScene.still([0.0], [G.dr], [1.0]), 64, 16)[32, 9]
    assert abs(np.angle(two / one) - np.angle(np.exp(-4j * np.pi * G.dr / G.lam))) < 0.05


def test_a_doppler_slice_is_the_image_of_that_stretch_of_pulses():
    rng = np.random.default_rng(3)
    n = 400
    sc = PointScene.still(rng.uniform(-4, 4, n) * ST.V_eff / ST.V_g, rng.uniform(-2, 2, n),
                          rng.standard_normal(n) + 1j * rng.standard_normal(n))
    full = simulate(ST, sc, 512, 16)
    part = simulate(ST, sc, 512, 16, pulses=(4.0, 8.0))
    same = looks(full, look_masks(G, 512, [6.0], 4.0))[0]
    mirror = looks(full, look_masks(G, 512, [-6.0], 4.0))[0]
    c = lambda a, b: abs(np.vdot(a.ravel(), b.ravel())) / np.linalg.norm(a) / np.linalg.norm(b)
    assert c(part, same) > 0.97
    assert c(part, mirror) < 0.05


def test_a_steady_recession_shifts_the_image_by_r_v_over_vs():
    v, f = 5e-3, 1e-4
    still = simulate(ST, PointScene.still([0.0], [0.0], [1.0]), 512, 8)
    moved = simulate(ST, PointScene.still([0.0], [0.0], [1.0]).moving(np.array([True]), ar=v / (2 * np.pi * f), f=f), 512, 8)

    def peak(im):
        F = np.fft.fft(np.abs(im[:, 4]))
        Z = np.zeros(len(F) * 64, complex)
        Z[:256], Z[-256:] = F[:256], F[-256:]
        return np.abs(np.fft.ifft(Z)).argmax() / 64
    shift = (peak(moved) - peak(still)) * G.dx
    theory = v * 2 * G.V / (G.lam * G.Ka_signed)
    assert abs(shift - theory) < 0.02 * abs(theory)
