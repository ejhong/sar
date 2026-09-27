"""Independent null, paired-response and gate checks for the new T7."""
import numpy as np
import pytest

from experiments.t07_surface_controls import paired_region_metrics, gate_counts


def test_paired_power_uses_the_same_region_and_keeps_absolute_scale():
    baseline = np.arange(1, 13, dtype=float).reshape(2, 2, 3)
    region = np.array([[True, True], [False, False]])
    changed = baseline.copy()
    changed[region] *= 0.25
    changed[~region] *= 100
    result = paired_region_metrics(changed, baseline, region)
    assert result["power_ratio_to_baseline"] == pytest.approx(0.25)
    assert result["fraction_depth_bins_below_baseline"] == 1.0
    assert result["mean_raw_power"] == pytest.approx(result["baseline_mean_raw_power"] / 4)


def test_identical_inputs_have_no_paired_suppression():
    baseline = np.arange(1, 13, dtype=float).reshape(2, 2, 3)
    region = np.ones((2, 2), dtype=bool)
    result = paired_region_metrics(baseline.copy(), baseline, region)
    assert result["power_ratio_to_baseline"] == 1.0
    assert result["fraction_depth_bins_below_baseline"] == 0.0
    with pytest.raises(ValueError, match="no grid targets"):
        paired_region_metrics(baseline, baseline, ~region)


def test_full_gate_distinguishes_resolved_ellipse_small_ellipse_and_zero():
    phase = 2 * np.pi * np.arange(50) / 25
    ellipse = np.column_stack((0.02 * np.cos(phase), 0.01 * np.sin(phase)))
    counts = gate_counts(np.stack((ellipse, ellipse / 10, np.zeros_like(ellipse))))
    assert counts["shape_gate_count"] == 2
    assert counts["full_gate_count"] == 1
    assert counts["full_gate_fraction"] == pytest.approx(1 / 3)
