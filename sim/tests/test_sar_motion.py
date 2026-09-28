"""Phase 2 groundwork: the SLC synthesizer moves scatterers by any motion field."""
import numpy as np

from sarsim import Geometry, point_targets, synthesize


def _scene(vib_amp=0.0):
    x = np.array([-20.0, 5.0, 30.0])
    y = np.array([10.0, -15.0, 25.0])
    z = np.zeros(3)
    return point_targets(x, y, z, np.ones(3), vib_amp=vib_amp, vib_freq=0.7, vib_phase=0.3)


def test_a_motion_field_reproduces_the_sinusoid_path():
    g = Geometry()
    shape = (256, 256)
    a = 2e-3
    ref = synthesize(_scene(vib_amp=a), g, shape)
    moved = synthesize(_scene(), g, shape,
                       motion=lambda x, y, z, t: a * np.sin(2 * np.pi * 0.7 * t[:, None] + 0.3) * np.ones(len(x)))
    assert np.allclose(moved, ref, rtol=0, atol=1e-6 * np.abs(ref).max())


def test_no_motion_leaves_the_image_unchanged():
    g = Geometry()
    still = synthesize(_scene(), g, (128, 128))
    zero = synthesize(_scene(), g, (128, 128), motion=lambda x, y, z, t: np.zeros((len(t), len(x))))
    assert np.array_equal(still, zero)


def test_double_precision_keeps_a_billionth_of_a_radian():
    g = Geometry()
    shape = (128, 128)
    lam = g.lam
    d = 1e-9 * lam / (4 * np.pi)                     # a slant-range change worth 1e-9 rad of phase
    still = synthesize(_scene(), g, shape, dtype=np.complex128)
    moved = synthesize(_scene(), g, shape, dtype=np.complex128,
                       motion=lambda x, y, z, t: d * np.ones((len(t), len(x))))
    assert still.dtype == np.complex128
    ratio = np.vdot(still, moved) / np.vdot(still, still)
    assert abs(np.angle(ratio) + 1e-9) < 1e-12        # the phase change comes back to a thousandth of itself
