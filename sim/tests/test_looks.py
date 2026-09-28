"""P2-03: what one dwell can see of the ground's motion."""
import numpy as np
import pytest

from sarsim import synthesize
from sarsim.acquisition import DwellGeometry
from sarsim.looks import (coherence, coherence_from_power, inject_region_motion, look_gain, look_masks, looks,
                          paired_velocity, region_shift, velocity_series)
from sarsim.scene import Scatterers

G = DwellGeometry.from_record('giza-20250827')


def speckle(n_rows, n_cols, rng):
    return (rng.normal(size=(n_rows, n_cols)) + 1j * rng.normal(size=(n_rows, n_cols))).astype(np.complex64)


def test_look_coherence_is_a_function_of_the_power_spectrum_alone():
    rng = np.random.default_rng(0)
    img = speckle(1024, 32, rng)
    m = look_masks(G, 1024, [0.0, 0.5, 1.2], 1.96)
    L = looks(img, m)
    for k in (1, 2):
        assert coherence(L[0], L[k]) == pytest.approx(coherence_from_power(img, m[0], m[k]), abs=1e-6)
    # a flat spectrum gives 1 - delta / W
    flat = np.fft.ifft(np.exp(2j * np.pi * rng.random((1024, 1))), axis=0) * np.ones((1, 4))
    assert coherence_from_power(flat, m[0], m[1]) == pytest.approx(1 - 0.5 / 1.96, abs=0.01)


def smooth_block(n, lo, hi, ramp=40):
    x = np.arange(n, dtype=float)
    return np.clip(np.minimum(x - lo, hi - x) / ramp, 0, 1) ** 2 * (3 - 2 * np.clip(np.minimum(x - lo, hi - x) / ramp, 0, 1))


def test_complex_looks_of_one_image_do_not_see_ground_move():
    rng = np.random.default_rng(1)
    img = speckle(2048, 32, rng)
    w = smooth_block(2048, 800, 1250)
    a = 1e-3                                                                 # a strong acceleration of one block
    moved = inject_region_motion(img, G, w, lambda t: a * t ** 2 / 2)
    centres = np.linspace(-0.9, 0.9, 7)
    m = look_masks(G, 2048, centres, 1.96)
    rows, cols = slice(900, 1150), slice(0, 32)
    single = velocity_series(looks(moved, m), G, rows, cols) - velocity_series(looks(img, m), G, rows, cols)
    paired = paired_velocity(looks(img, m), looks(moved, m), G, rows, cols)
    injected = a * centres
    assert paired == pytest.approx(injected, abs=0.05 * np.abs(injected).max())      # the motion is in the data
    assert np.abs(single).max() < 0.05 * np.abs(injected - injected[3]).max()         # one image cannot see it


def test_the_moved_block_is_displaced_in_each_look():
    rng = np.random.default_rng(2)
    img = speckle(2048, 32, rng)
    w = smooth_block(2048, 800, 1250)
    v = 1.5e-3
    moved = inject_region_motion(img, G, w, lambda t: v * t)
    m = look_masks(G, 2048, [-0.5, 0.5], 1.96)
    got = paired_velocity(looks(img, m), looks(moved, m), G, slice(900, 1150), slice(0, 32))
    assert got == pytest.approx([v, v], rel=0.06)             # the region estimator reads a few per cent short on a narrow look
    far = moved[:600] - img[:600]                  # only the fractional shift's interpolation tails reach this far
    assert np.sqrt(np.mean(np.abs(far) ** 2)) < 0.02 * np.sqrt(np.mean(np.abs(img) ** 2))


def point_scene(extra_amp):
    rng = np.random.default_rng(3)
    n = 20000
    L = 512 * G.dx
    x = np.r_[rng.uniform(-L / 2, L / 2, n), 0.0]
    r = np.r_[rng.uniform(-3, 3, n), 0.0]
    N = n + 1
    return Scatterers(x=x, y=r / np.sin(G.theta), z=np.zeros(N), amp=np.r_[np.ones(n), extra_amp],
                      phase=rng.uniform(0, 2 * np.pi, N), iso=np.ones(N), flash=np.zeros(N), nu0=np.zeros(N),
                      sig_nu=np.ones(N), vib_amp=np.zeros(N), vib_freq=np.zeros(N), vib_phase=np.zeros(N),
                      label=np.zeros(N, int))


def test_only_the_envelope_of_a_bright_point_follows_its_motion():
    shape = (512, 40)
    a = 3e-3
    motion = lambda x, y, z, t: a * t[:, None] ** 2 / 2 * (x == 0.0)[None, :]     # constant acceleration a
    still = synthesize(point_scene(200.0), G, shape)
    moved = synthesize(point_scene(200.0), G, shape, motion=motion)
    centres = np.array([-0.8, 0.0, 0.8])
    m = look_masks(G, 512, centres, 1.2)
    rows, cols = slice(200, 312), slice(14, 26)
    env = velocity_series(looks(moved, m), G, rows, cols, envelope=True) - \
        velocity_series(looks(still, m), G, rows, cols, envelope=True)
    cpx = velocity_series(looks(moved, m), G, rows, cols) - velocity_series(looks(still, m), G, rows, cols)
    expected = a * (centres - centres[1])
    assert env == pytest.approx(expected, abs=0.1 * np.abs(expected).max())
    # a windowed complex correlation catches only part of a bright point's motion, through the window's edges
    assert np.abs(cpx).max() < 0.5 * np.abs(expected).max()


def test_a_look_averages_motion_over_its_duration():
    assert look_gain(0.0, 1.96) == 1.0
    assert abs(look_gain(1.0, 1.96)) < 0.03 and abs(look_gain(2.0, 1.96)) < 0.03


def test_motionless_twin_keeps_brightness_and_spectrum_not_speckle():
    from sarsim.looks import motionless_twin
    rng = np.random.default_rng(3)
    base = (rng.standard_normal((256, 128)) + 1j * rng.standard_normal((256, 128))).astype(np.complex64)
    base[:, 64:] *= 3.0                                   # a brighter half
    t = motionless_twin(base, np.random.default_rng(4))
    ratio = np.mean(np.abs(t[:, 64:]) ** 2) / np.mean(np.abs(t[:, :64]) ** 2)
    assert 7.0 < ratio < 11.0                             # the brighter half stays about 9 times brighter
    inner = (slice(8, -8), slice(72, 120))                # inside the bright half: fresh speckle, not the same
    assert abs(np.corrcoef(np.abs(t[inner]).ravel(), np.abs(base[inner]).ravel())[0, 1]) < 0.1
