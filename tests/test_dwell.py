import numpy as np
import pytest
import h5py

from sarsim.dwell import DwellProduct, DopplerBank, occupied_band, azimuth_shift_to_velocity

PRF = 160000.0
BAND = 138000.0
N_AZ, N_RG = 512, 64


def write_product(path, **overrides):
    """A minimal but structurally faithful ICEYE dwell product."""
    rng = np.random.default_rng(0)
    fields = dict(
        processing_prf=PRF, azimuth_time_interval=1.0 / PRF, acquisition_prf=6000.0,
        carrier_frequency=9.6e9, total_processed_bandwidth_azimuth=BAND,
        number_of_azimuth_samples=N_AZ, number_of_range_samples=N_RG,
        slant_range_to_first_pixel=600000.0, slant_range_spacing=0.15,
        azimuth_ground_spacing=0.044, incidence_center=20.8, first_pixel_time=0.004,
        range_sampling_rate=9.7e8, avg_scene_height=70.0, mean_earth_radius=6371000.0)
    fields.update(overrides)
    text = dict(acquisition_start_utc='2025-08-27T20:26:54.000000',
                acquisition_end_utc='2025-08-27T20:27:19.000000',
                zerodoppler_start_utc='2025-08-27T20:27:06.000000',
                zerodoppler_end_utc='2025-08-27T20:27:06.713000',
                product_type='SpotlightDwellFine', polarization='VV',
                satellite_name='TEST-X', look_side='right', orbit_direction='ASCENDING')
    with h5py.File(path, 'w') as f:
        for k, v in fields.items():
            f[k] = v
        for k, v in text.items():
            f[k] = np.bytes_(v)
        f['doppler_rate_coeffs'] = np.array([-5662.5, 0.0, 0.0, 0.0])
        f['s_i'] = rng.normal(0, 30, (N_AZ, N_RG)).astype(np.int16)
        f['s_q'] = rng.normal(0, 30, (N_AZ, N_RG)).astype(np.int16)
        t = np.linspace(-13.0, 13.0, 27)
        f['state_vector_time_utc'] = np.array(
            [[np.bytes_((np.datetime64('2025-08-27T20:27:06.500000') + np.timedelta64(int(x * 1e6), 'us')).item().isoformat())]
             for x in t])
        r, w = 7.0e6, 1.0e-3
        f['posX'], f['posY'], f['posZ'] = r * np.cos(w * t), r * np.sin(w * t), np.zeros_like(t)
        f['velX'], f['velY'], f['velZ'] = -r * w * np.sin(w * t), r * w * np.cos(w * t), np.zeros_like(t)
        f['coa_pos'] = np.array([r, 0.0, 0.0])
        f['coa_vel'] = np.array([0.0, r * w, 0.0])
        for name in ('LINE', 'SAMP'):
            f[f'RPC/{name}_OFF'] = 0.0
            f[f'RPC/{name}_SCALE'] = 1.0
            f[f'RPC/{name}_NUM_COEFF'] = np.zeros(20)
            f[f'RPC/{name}_DEN_COEFF'] = np.eye(20)[0]
        for name in ('LAT', 'LONG', 'HEIGHT'):
            f[f'RPC/{name}_OFF'] = 0.0
            f[f'RPC/{name}_SCALE'] = 1.0
    return path


@pytest.fixture
def product(tmp_path):
    return DwellProduct(write_product(tmp_path / 'p.h5'))


def test_acquisition_facts_and_aperture(product):
    a = product.acq
    assert a.shape == (N_AZ, N_RG)
    assert a.collection_duration_s == pytest.approx(25.0)
    assert a.image_span_s == pytest.approx(0.713)
    assert a.wavelength_m == pytest.approx(0.031228, rel=1e-4)
    # the aperture is bandwidth over Doppler rate, not the image row span
    assert product.aperture_time_s(N_RG // 2) == pytest.approx(BAND / 5662.5, rel=1e-6)
    assert product.aperture_time_s(N_RG // 2) > 30 * a.image_span_s


def test_ephemeris_matches_declared_centre_of_aperture(product):
    pos, vel = product.ephemeris().at([0.0])
    with h5py.File(product.path, 'r') as f:
        assert np.allclose(pos[0], f['coa_pos'][()], atol=1e-3)
        assert np.allclose(vel[0], f['coa_vel'][()], atol=1e-3)


def test_missing_field_is_refused(tmp_path):
    path = write_product(tmp_path / 'bad.h5')
    with h5py.File(path, 'a') as f:
        del f['total_processed_bandwidth_azimuth']
    with pytest.raises(ValueError, match='missing required fields'):
        DwellProduct(path)


def test_inconsistent_timing_is_refused(tmp_path):
    path = write_product(tmp_path / 'bad2.h5', azimuth_time_interval=1.0 / (PRF * 2))
    with pytest.raises(ValueError, match='inconsistent'):
        DwellProduct(path)


def test_positive_doppler_rate_is_refused(tmp_path):
    path = write_product(tmp_path / 'bad3.h5')
    with h5py.File(path, 'a') as f:
        del f['doppler_rate_coeffs']
        f['doppler_rate_coeffs'] = np.array([+5662.5, 0.0, 0.0, 0.0])
    with pytest.raises(ValueError, match='Doppler rate'):
        DwellProduct(path).doppler_rate_hz_s(0)


def test_crop_is_bounded_and_centred(product):
    data, (r0, c0) = product.crop(100, 32, 64, 16)
    assert data.shape == (64, 16) and data.dtype == np.complex64
    assert r0 == 68 and c0 == 24
    with pytest.raises(ValueError, match='budget'):
        product.crop(100, 32, 512, 64, max_bytes=1024)
    with pytest.raises(ValueError, match='larger than the image'):
        product.crop(100, 32, 4096, 16)


def test_bank_plan_maps_frequency_to_slow_time(product):
    bank = DopplerBank(K=50, width_frac=0.10)
    plan = bank.plan(product.acq, N_AZ, -5662.5)
    assert plan['sub_aperture_s'] == pytest.approx(plan['width_hz'] / 5662.5)
    assert plan['span_s'] > 20.0                         # a dwell really does span seconds
    assert plan['sampling_hz'] == pytest.approx(49 / plan['span_s'])
    assert plan['pair_lag_s'] == 0.0                     # common-reference mode
    assert plan['overlap_within_pair'] == 1.0
    assert np.all(np.diff(plan['times_s']) < 0)          # negative rate reverses the time order


def test_bank_pair_mode_overlap_is_nearly_total(product):
    plan = DopplerBank(K=20, width_frac=0.50, shift_hz=404.0).plan(product.acq, N_AZ, -5662.5)
    assert plan['overlap_within_pair'] > 0.98
    assert plan['pair_lag_s'] == pytest.approx(plan['shift_hz'] / 5662.5)


def test_bank_refuses_impossible_designs(product):
    with pytest.raises(ValueError):
        DopplerBank(K=1).plan(product.acq, N_AZ, -5662.5)
    with pytest.raises(ValueError):
        DopplerBank(K=50, width_frac=0.99, shift_hz=5000.0).plan(product.acq, N_AZ, -5662.5)
    with pytest.raises(ValueError, match='Doppler rate'):
        DopplerBank().plan(product.acq, N_AZ, +5662.5)


def test_masks_cover_the_declared_support(product):
    bank = DopplerBank(K=10, width_frac=0.2)
    plan = bank.plan(product.acq, N_AZ, -5662.5)
    freqs = np.fft.fftshift(np.fft.fftfreq(N_AZ, d=1.0 / PRF))
    ref, off = bank.masks(plan, freqs)
    assert ref.shape == (10, N_AZ)
    assert np.all(ref.sum(axis=1) > 0)
    assert np.all(np.abs(freqs[ref[0] > 0]) <= BAND / 2 + plan['bin_hz'])


def test_occupied_band_recovers_a_planted_support(product):
    rng = np.random.default_rng(1)
    n = 4096
    spec = np.zeros((n, 4), complex)
    freqs = np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / PRF))
    inside = np.abs(freqs) < 30000
    spec[inside] = rng.normal(size=(inside.sum(), 4)) + 1j * rng.normal(size=(inside.sum(), 4))
    crop = np.fft.ifft(np.fft.ifftshift(spec, axes=0), axis=0)
    band = occupied_band(crop, PRF)
    assert band['width_hz'] == pytest.approx(60000, rel=0.05)
    assert abs(band['centroid_hz']) < 2000


def test_velocity_conversion_is_the_published_lever(product):
    # one azimuth pixel of shift at R/V = 78 s corresponds to a sub-mm/s velocity
    v = azimuth_shift_to_velocity(1.0, product.acq, 600000.0, 7661.0)
    assert v == pytest.approx(0.044 * 7661.0 / 600000.0)
    assert 0.0005 < v < 0.001
    with pytest.raises(ValueError):
        azimuth_shift_to_velocity(1.0, product.acq, -1.0, 7661.0)


def test_split_aperture_halves_have_matching_sweeps(product):
    a = DopplerBank(K=25, width_frac=0.10, start_frac=0.0, sweep_frac=0.45).plan(product.acq, N_AZ, -5662.5)
    b = DopplerBank(K=25, width_frac=0.10, start_frac=0.45, sweep_frac=0.45).plan(product.acq, N_AZ, -5662.5)
    assert a['span_s'] == pytest.approx(b['span_s'], rel=1e-9)
    assert a['step_hz'] == pytest.approx(b['step_hz'])
    # the two halves cover disjoint stretches of slow time
    assert min(a['times_s']) > max(b['times_s']) or min(b['times_s']) > max(a['times_s'])


def test_sweep_that_does_not_fit_is_refused(product):
    with pytest.raises(ValueError, match='does not fit'):
        DopplerBank(K=10, width_frac=0.5, start_frac=0.4, sweep_frac=0.9).plan(product.acq, N_AZ, -5662.5)
