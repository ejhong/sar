import h5py
import numpy as np
import pytest

from sarsim.geolocation import audit_iceye_geometry, project_rpc, CONTROL_POINTS


def camera():
    # Deliberately different latitude/longitude scales expose axis swaps.
    rpc = dict(LAT_OFF=30., LONG_OFF=31., HEIGHT_OFF=100.,
               LAT_SCALE=2., LONG_SCALE=4., HEIGHT_SCALE=20.,
               LINE_OFF=200., SAMP_OFF=400., LINE_SCALE=100., SAMP_SCALE=50.)
    for axis in ('LINE', 'SAMP'):
        rpc[f'{axis}_NUM_COEFF'] = np.zeros(20)
        rpc[f'{axis}_DEN_COEFF'] = np.r_[1., np.zeros(19)]
    # row = 200 + 100 * P; sample = 400 + 50 * (L-H).
    rpc['LINE_NUM_COEFF'][2] = 1.
    rpc['SAMP_NUM_COEFF'][1] = 1.
    rpc['SAMP_NUM_COEFF'][3] = -1.
    return rpc


def test_forward_projection_preserves_axes_height_and_pixel_centres():
    rpc = camera()
    np.testing.assert_allclose(project_rpc(rpc, 31., 33., 105.), [250., 412.5])
    np.testing.assert_allclose(project_rpc(rpc, [30., 31.], 33., [100., 110.]),
                               [[200., 425.], [250., 400.]])
    # Nonlinear cross term P*H^2 and rational denominator (1+L).
    rpc['LINE_NUM_COEFF'][16] = 2.
    rpc['LINE_DEN_COEFF'][1] = 1.
    np.testing.assert_allclose(project_rpc(rpc, 31., 33., 110.), [250., 400.])
    rpc['LINE_SCALE'] = -100.
    np.testing.assert_allclose(project_rpc(rpc, 31., 33., 110.), [150., 400.])


def test_rpc_rejects_missing_invalid_and_singular_metadata():
    rpc = camera()
    del rpc['HEIGHT_OFF']
    with pytest.raises(ValueError, match='Missing RPC'):
        project_rpc(rpc, 30., 31., 100.)
    rpc = camera()
    rpc['LAT_SCALE'] = 0.
    with pytest.raises(ValueError, match='nonzero'):
        project_rpc(rpc, 30., 31., 100.)
    rpc = camera()
    rpc['LINE_DEN_COEFF'][:] = 0.
    with pytest.raises(ValueError, match='Singular'):
        project_rpc(rpc, 30., 31., 100.)
    with pytest.raises(ValueError, match='extrapolation'):
        project_rpc(camera(), 40., 31., 100.)


def product(path):
    with h5py.File(path, 'w') as f:
        # Several GB if read in full; no chunks are allocated.
        for name in ('s_i', 's_q'):
            f.create_dataset(name, shape=(50000, 40000), dtype='float32', chunks=(64, 64))
        f['number_of_azimuth_samples'] = 50000
        f['number_of_range_samples'] = 40000
        f['data_orientation'] = 'native'
        f['geo_ref_system'] = 'WGS84'
        for key, value in camera().items():
            f[f'RPC/{key}'] = value
        for name in CONTROL_POINTS:
            f[name] = [401., 201., 30., 31.]
        f['avg_scene_height'] = 100.
        f['incidence_center'] = 30.
        f['slant_range_spacing'] = .2


def test_audit_reads_no_pixels_and_does_not_promote_internal_agreement(tmp_path, monkeypatch):
    path = tmp_path/'native.h5'
    product(path)
    original_read = h5py.Dataset.__getitem__
    def guarded_read(self, key):
        assert self.name not in ('/s_i', '/s_q'), 'Raster read during metadata audit'
        return original_read(self, key)
    monkeypatch.setattr(h5py.Dataset, '__getitem__', guarded_read)
    sidecar = tmp_path/'sicd.xml'
    sidecar.write_text('<SICD xmlns="urn:SICD:1.3.0"><ImageData><NumRows>40000</NumRows>'
                       '<NumCols>49998</NumCols></ImageData></SICD>')
    report = audit_iceye_geometry(path, sidecar)
    np.testing.assert_allclose(report['control_point_checks'][0]['residual_to_zero_based_metadata_px'], [0., 0.])
    np.testing.assert_allclose(report['control_point_checks'][0]['height_plus_10m_shift_px'], [0., -25.])
    assert report['raster_samples_read'] == 0
    assert not report['native_registration_validated']
    assert not report['depth_analysis_ready']
    assert not report['companion_sicd']['transposed_shape_matches_native']
    assert len(report['source_metadata_sha256']) == 64


def test_audit_refuses_axis_and_dimension_assumptions(tmp_path):
    path = tmp_path/'native.h5'
    product(path)
    with h5py.File(path, 'r+') as f:
        f['data_orientation'][()] = 'shadows_down'
    with pytest.raises(ValueError, match='native ICEYE orientation'):
        audit_iceye_geometry(path)
    with h5py.File(path, 'r+') as f:
        f['data_orientation'][()] = 'native'
        f['number_of_range_samples'][()] = 39999
    with pytest.raises(ValueError, match='dimensions disagree'):
        audit_iceye_geometry(path)
