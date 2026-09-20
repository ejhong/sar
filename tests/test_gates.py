"""Selection-order regression and independent invariants for the corrected branch."""
import numpy as np
import pytest
from sarsim.gates import ellipse_windows, ellipse_gate, two_arc_correction
from sarsim.tomo import focus_branch_b, focus_windows


def test_a_better_rejected_window_does_not_hide_an_accepted_one():
    rng = np.random.default_rng(4)
    t = 2*np.pi*np.arange(25)/25
    ellipse = np.column_stack((.02*np.cos(t), .01*np.sin(t)))
    ellipse += rng.normal(0, .0003, ellipse.shape)
    line = np.column_stack((.02*np.cos(t), np.zeros(25)))
    q = np.concatenate((ellipse, line))[None]
    gates = ellipse_windows(q)
    assert gates.full[0, 0]
    assert not gates.full[0, -1]
    assert gates.adjusted_r2[0, -1] > gates.adjusted_r2[0, 0]
    assert ellipse_gate(q)[-1][0]


def test_rejected_windows_never_supply_a_branch_b_depth():
    t = 2*np.pi*np.arange(50)/25
    q = np.column_stack((.02*np.cos(t), .01*np.sin(t)))[None]
    kz = np.arange(50)*.01
    z = np.linspace(1, 150, 80)
    score, winning, gates = focus_branch_b(q, kz, z)
    for winner in winning[0]:
        assert gates.full[0, winner]
    null_score, null_w, _ = focus_branch_b(q/100, kz, z)
    assert np.isnan(null_score).all()
    assert (null_w == -1).all()
    a, aw = focus_windows(q, kz, z)
    all_windows, all_w = focus_windows(q, kz, z, accepted=np.ones((1, 26), bool))
    np.testing.assert_allclose(a, all_windows)
    np.testing.assert_array_equal(aw, all_w)


def test_two_arcs_keep_a_target_signal_and_resist_one_outlier():
    rng = np.random.default_rng(93)
    common = rng.normal(size=(50, 2))
    plus = np.broadcast_to(common, (11, 50, 2)).copy()
    minus = plus.copy()
    plus[0] += 100
    target = common[None] + .1
    q, qc = two_arc_correction(target, plus, minus)
    np.testing.assert_allclose(q, .1, atol=1e-8)
    np.testing.assert_allclose(q+qc['common'], target)
    assert qc['disagreement'].max() < 1e-8
    with pytest.raises(ValueError):
        ellipse_windows(np.zeros((1, 20, 2)))
