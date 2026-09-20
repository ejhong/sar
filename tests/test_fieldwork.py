"""Real-file processing invariants and honest publication of field checks."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from fieldwork.processing import (bilinear_geoid, doppler_rate, full_row_spectrum,
                                  intensity_preview, timing_audit, translation_trials)
from fieldwork.publishing import obj_model, validate_model

ROOT = Path(__file__).resolve().parents[1]


def test_preview_streaming_keeps_partial_edges_and_hashes_original_values():
    class BoundedArray:
        def __init__(self, values):
            self.values, self.shape, self.reads = values, values.shape, []
        def __getitem__(self, selection):
            result = self.values[selection]
            self.reads.append(result.shape)
            assert result.shape[0] <= 4
            return result
    i = np.arange(77, dtype=np.float32).reshape(11,7)
    q = np.ones_like(i)*3
    real, imag = BoundedArray(i), BoundedArray(q)
    result, hashes = intensity_preview(real,imag,[1,11],[1,7],[4,4],max_chunk_bytes=4*6*4*5)
    expected = np.array([[np.mean(i[a:min(a+4,11),b:min(b+4,7)].astype(float)**2+9)
                          for b in (1,5)] for a in (1,5,9)])
    np.testing.assert_allclose(result,expected)
    assert real.reads == [(4,6),(4,6),(2,6)]
    assert hashes['i_samples_sha256'] == hashlib.sha256(i[1:11,1:7].astype('<f4').tobytes()).hexdigest()
    assert hashes['q_samples_sha256'] == hashlib.sha256(q[1:11,1:7].astype('<f4').tobytes()).hexdigest()
    with pytest.raises(ValueError,match='budget'):
        intensity_preview(real,imag,[0,11],[0,7],[10,2],max_chunk_bytes=100)
    i[-1,-1] = np.nan
    with pytest.raises(ValueError,match='Nonfinite'):
        intensity_preview(real,imag,[1,11],[1,7],[4,4],max_chunk_bytes=4*6*4*5)


def test_doppler_rate_uses_range_offset_not_absolute_delay():
    coefficients = [-5662.497777974439,1396578.9641082047,-344448059.8265285,84953953226.37804]
    rate = doppler_rate(coefficients,.004048625065198723,975879241.8017589,11562,[0,5781,11561])
    assert rate[1] == coefficients[0]
    assert -5671 < rate[0] < -5670
    assert -5655 < rate[2] < -5654
    # Changing the absolute range origin cannot move the polynomial reference.
    np.testing.assert_array_equal(rate,doppler_rate(coefficients,.04,975879241.8017589,11562,[0,5781,11561]))
    with pytest.raises(ValueError,match='coordinates'):
        doppler_rate(coefficients,.004,1e9,11562,[-1])


def test_collection_duration_does_not_become_image_row_duration():
    meta=json.loads((ROOT/'results/field/giza/report.json').read_text())['source_metadata']
    timing=timing_audit(meta)
    assert 24 < timing['collection_duration_s'] < 25
    assert .7 < timing['zero_doppler_image_span_s'] < .8
    assert 24 < timing['nominal_bandwidth_over_rate_s'][1] < 25
    # Published endpoint timestamps and N-1 intervals differ by ~1.1 output rows.
    # Keep that metadata difference; neither duration is the ~25 s collection.
    assert abs(timing['image_span_from_sampling_s']-timing['zero_doppler_image_span_s']) < 2*meta['azimuth_time_interval']
    assert timing['aperture_time_validated'] is False
    with pytest.raises(ValueError,match='timing'):
        timing_audit(dict(meta,processing_prf=meta['acquisition_prf']))


def test_full_row_spectrum_retains_complex_sign_and_outside_power():
    n=1024;rows=np.arange(n)
    # A negative-frequency complex tone and a positive tone outside the band.
    samples=(np.exp(-2j*np.pi*100*rows/n)+.5*np.exp(2j*np.pi*400*rows/n))[:,None]
    result=full_row_spectrum(samples,1024,400,bins=1024)
    assert result['outside_zero_centered_metadata_band_fraction'] == pytest.approx(.2)
    peak=max(result['curve'],key=lambda r:r['mean_power'])
    assert peak['image_frequency_hz'] == -100
    with pytest.raises(ValueError,match='complex'):
        full_row_spectrum(samples.real,1024,400)


def test_geoid_conversion_is_interpolated_and_cannot_extrapolate():
    tile=dict(west_deg=31.,north_deg=30.,spacing_deg=.25,values_m=[[15.,16.],[17.,18.]])
    N=bilinear_geoid(tile,29.875,31.125)
    assert N == 16.5
    assert 58+N == 74.5  # h = H + N; not a zero altitude from a KML point.
    with pytest.raises(ValueError,match='outside'):
        bilinear_geoid(tile,29.,31.1)


def test_real_texture_control_keeps_signed_inputs_and_zero_trial():
    rng=np.random.default_rng(812)
    sample=rng.normal(size=(192,192))+1j*rng.normal(size=(192,192))
    trials=translation_trials(sample,[[0,0],[.005,0],[0,-.005]],patch=32)
    assert len(trials) == 75
    zero=[t for t in trials if t['injected_row_px']==t['injected_col_px']==0]
    assert max(t['error_px'] for t in zero) < 1e-10
    row=[t for t in trials if t['injected_row_px']>0]
    col=[t for t in trials if t['injected_col_px']<0]
    assert .003 < np.median([t['measured_row_px'] for t in row]) < .007
    assert -.007 < np.median([t['measured_col_px'] for t in col]) < -.003


def test_reference_meshes_preserve_datum_and_do_not_claim_detection():
    models=json.loads((ROOT/'catalog/giza/surveys.json').read_text())['models']
    for model in models:
        validate_model(model)
        mesh=obj_model(model)
        vertices=np.array([[float(v) for v in line.split()[1:]] for line in mesh.splitlines() if line.startswith('v ')])
        assert vertices[:,2].max()==0
        assert vertices[:,2].min()==-model['depth_m']
        assert model['native_target_mask'] is None
        with pytest.raises(ValueError,match='radar result'):
            validate_model(dict(model,field_detection=True))
    hetepheres=next(m for m in models if m['id']=='hetepheres')
    assert hetepheres['parts'][0]['bottom_m']==-27.42
    assert hetepheres['parts'][1]['top_m']==-25.50
    assert hetepheres['parts'][1]['bottom_m']==-27.45


@pytest.mark.parametrize('site',['giza','sacsayhuaman'])
def test_field_reports_keep_failures_and_every_trial_reproduces_tolerance(site):
    path=ROOT/'results/field'/site
    report=json.loads((path/'report.json').read_text())
    with (path/'translation_trials.csv').open() as f:
        trials=list(csv.DictReader(f))
    assert len(trials)==200
    for group in report['translation']['groups']:
        selected=[r for r in trials if float(r['injected_row_px'])==group['injected_row_px'] and float(r['injected_col_px'])==group['injected_col_px']]
        error=np.array([np.hypot(float(r['measured_row_px'])-float(r['injected_row_px']),float(r['measured_col_px'])-float(r['injected_col_px'])) for r in selected])
        assert len(error)==group['patches']==25
        assert np.percentile(error,95)==pytest.approx(group['p95_error_px'])
        assert (group['p95_error_px']<=report['translation']['max_p95_error_px'])==group['meets_0_02px_p95_tolerance']
    assert any(not g['meets_0_02px_p95_tolerance'] for g in report['translation']['groups'])
    assert report['native_registration_validated'] is report['aperture_time_validated'] is False
    assert report['field_depth_results']==0
    for filename,digest in report['output_sha256'].items():
        assert hashlib.sha256((path/filename).read_bytes()).hexdigest()==digest
    published=ROOT/'docs/data/field'/site/'report.json'
    assert json.loads(published.read_text())==report
