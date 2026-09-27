"""Unit tests for the pure helpers of the real-data chapter.

The experiments themselves need 10 GB products, so they are not exercised here. Everything
that turns cached arrays into a published number is.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from real_common import target_grid, patch_stats, GIZA_PATCHES, SACSAYHUAMAN_PATCHES, SITES, BANKS  # noqa: E402
from r03_split_dwell import common_profiles  # noqa: E402
from r04_velocity_floor import spectrum  # noqa: E402
from r05_learned_null import blocked_auc  # noqa: E402


def fake_run(pixels=64, depths=40, seed=0, shift=0.0, z_max=12.0):
    rng = np.random.default_rng(seed)
    z = np.linspace(z_max / depths, z_max, depths)
    tomo = rng.random((pixels, depths)) + shift
    return {"z": z, "tomo": tomo.astype(np.float32),
            "score": rng.random((pixels, depths)).astype(np.float32),
            "q": rng.normal(0, 0.003, (pixels, 50, 2)).astype(np.float32),
            "amplitude": rng.lognormal(3, 1, pixels).astype(np.float32),
            "velocity": rng.normal(0, 1e-4, (pixels, 50)).astype(np.float32),
            "meta": {"patch": "x", "label": "X", "kind": "control", "pixels": pixels,
                     "nyquist_depth_m": z_max, "velocity_kind": "line_of_sight_velocity",
                     "sub_aperture_times_s": list(np.linspace(10, -10, 50))}}


def test_sites_are_well_formed():
    assert set(SITES) == {"giza", "sacsayhuaman"}
    for patches in SITES.values():
        for name, label, lat, lon, height, kind in patches:
            assert isinstance(name, str) and name
            assert -90 <= lat <= 90 and -180 <= lon <= 180
            assert kind in {"monument", "control", "urban", "tombs"}
        assert len({p[0] for p in patches}) == len(patches)     # unique keys
    # every site needs at least one control to compare against
    for patches in SITES.values():
        assert any(p[5] == "control" for p in patches)


def test_target_grid_respects_the_margin():
    rows, cols, (r_axis, c_axis) = target_grid((4096, 1536), margin=160, stride_az=12, stride_rg=8)
    assert rows.min() >= 160 and cols.min() >= 160
    assert rows.max() < 4096 - 160 and cols.max() < 1536 - 160
    assert rows.size == r_axis.size * c_axis.size
    assert np.all(np.diff(r_axis) == 12) and np.all(np.diff(c_axis) == 8)


def test_target_grid_refuses_an_impossible_margin():
    rows, cols, _ = target_grid((100, 100), margin=60)
    assert rows.size == 0


def test_patch_stats_reports_finite_summaries():
    s = patch_stats(fake_run())
    for key in ("best_score_mean", "depth_median_m", "trajectory_rms_px", "velocity_rms_mm_s"):
        assert np.isfinite(s[key])
    assert 0 <= s["frac_score_gt_0p5"] <= 1
    assert s["kind"] == "control"


def test_common_profiles_uses_the_shorter_axis_and_normalises():
    a = fake_run(seed=1, z_max=12.0)
    b = fake_run(seed=2, z_max=9.0)
    z, pa, pb = common_profiles(a, b)
    assert z[-1] == pytest.approx(9.0)
    assert pa.shape == pb.shape == (64, 160)
    assert np.allclose(pa.sum(axis=1), 1.0)
    assert np.allclose(pb.sum(axis=1), 1.0)


def test_spectrum_finds_a_planted_tone():
    fs, n = 2.2, 50
    t = np.arange(n) / fs
    v = np.tile(1e-4 * np.sin(2 * np.pi * 0.5 * t), (20, 1))
    f, p = spectrum(v, fs)
    assert f[np.argmax(p)] == pytest.approx(0.5, abs=0.12)
    flat = np.tile(np.zeros(n), (20, 1))
    assert spectrum(flat, fs)[1].max() == 0.0


def test_blocked_auc_separates_and_returns_chance_when_it_should():
    rng = np.random.default_rng(0)
    n = 400
    blocks = np.repeat(np.arange(10), n // 10)
    y = (blocks % 2).astype(float)
    separable = (y[:, None] * 3 + rng.normal(size=(n, 2)))
    mean, std = blocked_auc(separable, y, blocks)
    assert mean > 0.9 and std >= 0
    noise = rng.normal(size=(n, 2))
    mean_null, _ = blocked_auc(noise, y, blocks)
    assert 0.2 < mean_null < 0.8


def test_bank_registry_covers_the_published_and_split_designs():
    assert {"reference", "paper", "v17", "split_a", "split_b"} <= set(BANKS)
    assert BANKS["reference"].shift_hz is None            # common-reference mode
    assert BANKS["paper"].shift_hz == 88.0
    assert BANKS["v17"].shift_hz == 404.0
    assert BANKS["split_a"].start_frac == 0.0 and BANKS["split_b"].start_frac == 0.45
    assert BANKS["split_a"].sweep_frac == BANKS["split_b"].sweep_frac
