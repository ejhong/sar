"""The ambient library traces: levels match their decibels and the measured run they cite."""
import json
from pathlib import Path

import numpy as np

from katabasis.ambient import peterson
from katabasis.runs import RESULTS

LIB = json.loads((Path(__file__).parents[2] / 'sites' / 'ambient.json').read_text())


def test_every_value_names_a_known_source():
    for group in ('bounds', 'regional', 'site_frequencies', 'cultural'):
        for k, v in LIB[group].items():
            if isinstance(v, dict):
                assert v.get('source') in LIB['sources'], f'{group}.{k}'


def test_cultural_levels_match_their_decibels():
    for k, v in LIB['cultural'].items():
        if isinstance(v, dict) and 'vdb' in v:
            assert abs(v['value'] / (10 ** (v['vdb'] / 20) * 2.54e-8) - 1) < 0.02, k


def test_measured_levels_match_the_run():
    run = json.loads((RESULTS / 'm1_01_ambient_levels' / 'summary.json').read_text())
    z = run['components']['BHZ']['band_rms_um_s']
    for key, band in (('microseism_vertical_0.1_0.3_hz', '0.1-0.3 Hz'), ('vertical_1_3_hz', '1.0-3.0 Hz'),
                      ('vertical_3_8_hz', '3.0-8.0 Hz')):
        v = LIB['regional'][key]
        assert abs(v['value'] / (z[band]['all_p50'] * 1e-6) - 1) < 0.03, key
        assert abs(v['range'][0] / (z[band]['all_p10'] * 1e-6) - 1) < 0.03, key
        assert abs(v['range'][1] / (z[band]['all_p90'] * 1e-6) - 1) < 0.03, key


def test_peterson_transcription_hits_the_published_end_points():
    # NLNM and NHNM at 0.1 s, from Tables 3 and 4: A + B log10(0.1)
    assert np.isclose(peterson.acceleration_db(0.1, 'low')[0], -162.36 - 5.64)
    assert np.isclose(peterson.acceleration_db(0.1, 'high')[0], -108.73 + 17.23)
    assert np.isnan(peterson.acceleration_db(0.05, 'low')[0])


def test_rayleigh_ellipse_of_a_poisson_solid():
    from katabasis.ambient.field import rayleigh_hv
    assert abs(rayleigh_hv(np.sqrt(3.0) * 1000, 1000.0) - 0.681) < 0.002


def test_plane_wave_field_has_the_target_spectrum():
    from scipy.signal import welch
    from katabasis.ambient.field import microseisms
    f = np.linspace(0.1, 0.5, 81)
    target = 1e-15 * np.exp(-((f - 0.2) / 0.07) ** 2)          # (m/s)²/Hz, a microseism-like peak
    fld = microseisms(f, target, np.random.default_rng(1), per_bin=6, hv=0.7)
    fs, T = 4.0, 4000.0
    t = np.arange(0, T, 1 / fs)
    v = fld.velocity(np.array([[0.0, 0.0, 0.0]]), t)[0]
    fw, pz = welch(v[2], fs=fs, nperseg=4096)
    band = (fw > 0.12) & (fw < 0.4)
    assert abs(np.trapezoid(pz[band], fw[band]) / np.trapezoid(target, f) - 1) < 0.15
    h = np.sqrt((v[0] ** 2 + v[1] ** 2).mean() / (v[2] ** 2).mean())
    assert abs(h - 0.7) < 0.05
