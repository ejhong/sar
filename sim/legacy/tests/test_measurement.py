"""Acquisition invariants and failure cases, not claims of field accuracy."""
import numpy as np
import pytest
from sarsim.geometry import Geometry
from sarsim.measurement import ReflectorAcquisition, WaveDictionary, centred, wilson
from sarsim.waves import WaveModel, surface_wave


def test_complex_slc_preserves_known_phase_motion_and_sign():
    acq = ReflectorAcquisition(Geometry(), shape=(512, 256))
    positions = [-40., 0., 40.]
    motion = np.column_stack([.0003*np.sin(2*np.pi*f*acq.time+.2) for f in (1., 2., 4.)])
    channels = acq.histories(motion, np.inf, np.random.default_rng(1))
    # Pure motion modulates phase, not the isolated pulse's amplitude.
    np.testing.assert_allclose(np.abs(channels), 1, atol=1e-14)
    measured = acq.read_phase(acq.separate(acq.focus(channels, positions), positions))
    np.testing.assert_allclose(measured, centred(motion), atol=1e-14)
    with pytest.raises(ValueError, match='not resolved'):
        acq.separate(acq.focus(channels, [0., 0., 0.]), [0., 0., 0.])


def test_unknown_positive_gain_does_not_change_dictionary_choice():
    t = np.linspace(0, 1, 100)
    waves = np.array([np.column_stack((np.sin(2*np.pi*f*t), np.cos(2*np.pi*f*t))) for f in [1, 2, 3]])
    estimator = WaveDictionary(waves, [0, 100, 200], [2000]*3)
    for gain in [.1, 1., 10.]:
        assert estimator.predict(gain*waves[1])['depth'][0] == 100
    assert wilson(0, 16)[1] > .1  # zero in a small sample is not a zero risk bound


def test_wave_arrival_respects_finite_propagation_and_free_cavity():
    model = WaveModel(duration=.28, spacing=5., receiver_step=150.)
    time, positions, wave = surface_wave(model)
    centre = np.flatnonzero(positions == 0)[0]
    assert np.max(np.abs(wave[time < .02, centre])) < 1e-25
    assert np.max(np.abs(wave[time > .08, centre])) > 1e-12
    _, _, cavity = surface_wave(WaveModel(depth=120., duration=.28, spacing=5., receiver_step=150.))
    assert np.linalg.norm(cavity-wave) > 0
    with pytest.raises(ValueError, match='wholly below'):
        surface_wave(WaveModel(depth=10.))
