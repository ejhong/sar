import numpy as np
from sarsim import Geometry, point_targets, synthesize, SubapBank, patch_shifts, focus_paper, focus_windows, kz_for_bank
from sarsim.track import shift_image
from sarsim.tomo import default_depths


def test_point_target_position_and_amplitude():
    g = Geometry()
    Nx, Nr = 256, 192
    # on-pixel positions: x = 3.0 m -> 12 px, y=10 m, z=4 m -> r = 10 sin - 4 cos
    x, y, z = 3.0, 10.0, 4.0
    sc = point_targets([x], [y], [z], [1.0], phase=[0.0])
    slc = synthesize(sc, g, (Nx, Nr))
    m = np.unravel_index(np.argmax(np.abs(slc)), slc.shape)
    r = y * np.sin(g.theta) - z * np.cos(g.theta)
    assert m[0] == Nx // 2 + round(x / g.dx)
    assert abs((m[1] - Nr // 2) * g.dr - r) <= g.dr  # within one pixel (r is off-grid)
    assert 0.6 < np.abs(slc[m]) <= 1.05


def test_ground_intensity_calibration():
    g = Geometry()
    slc = synthesize(None, g, (256, 256), ground_intensity=0.01, seed=1)
    assert abs(np.mean(np.abs(slc) ** 2) / 0.01 - 1) < 0.05


def test_patch_shifts_recover_subpixel_translation():
    g = Geometry()
    img = synthesize(None, g, (256, 256), ground_intensity=1.0, seed=2)
    d = (0.37, -0.21)
    img2 = shift_image(img, *d)
    rows = np.arange(40, 220, 20)
    cols = np.arange(40, 220, 20)
    RR, CC = np.meshgrid(rows, cols, indexing="ij")
    dr, dc = patch_shifts(img, img2, RR.ravel(), CC.ravel(), patch=32, upsample=1000)
    assert np.abs(np.median(dr) - d[0]) < 0.02  # ~3% window-edge bias at large fractional shifts
    assert np.abs(np.median(dc) - d[1]) < 0.02
    assert np.std(dr) < 0.05


def test_focus_recovers_injected_depth():
    g = Geometry()
    bank = SubapBank()
    kz = kz_for_bank(bank, g, 0.48)
    z = default_depths(kz, n=400)
    z0 = 40.0
    Yc = np.exp(1j * kz * z0)[None, :]
    t = focus_paper(Yc, kz, z)
    assert abs(z[np.argmax(t[0])] - z0) < 1.0
    q = np.stack([0.02 * np.cos(kz * z0), 0.01 * np.sin(kz * z0)], axis=1)[None, :, :]
    s, _ = focus_windows(q, kz, z)
    assert abs(z[np.argmax(s[0])] - z0) < 1.0
    assert s[0].max() > 0.99


def test_subaperture_masks_inside_band():
    g = Geometry()
    bank = SubapBank(K=10)
    mR, mO = bank.masks(g, 512)
    nu = np.fft.fftfreq(512, d=g.dx)
    assert np.all(np.abs(nu[mR.any(axis=0)]) <= g.nu_band / 2 + 1e-9)
    assert np.all(mR.sum(axis=1) > 0) and np.all(mO.sum(axis=1) > 0)


def test_patch_shifts_small_shift_is_precise():
    g = Geometry()
    img = synthesize(None, g, (256, 256), ground_intensity=1.0, seed=5)
    img2 = shift_image(img, 0.012, -0.007)
    rows = np.arange(40, 220, 20)
    cols = np.arange(40, 220, 20)
    RR, CC = np.meshgrid(rows, cols, indexing="ij")
    dr, dc = patch_shifts(img, img2, RR.ravel(), CC.ravel(), patch=32, upsample=1000)
    assert np.abs(np.median(dr) - 0.012) < 0.001 and np.abs(np.median(dc) + 0.007) < 0.001
    assert dr.std() < 0.001
