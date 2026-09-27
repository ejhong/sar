import h5py
import numpy as np
import pytest
import tifffile
from sarsim.realdata import inspect_product, read_crop, crop_bounds


def test_large_hdf5_can_be_inspected_and_cropped_without_full_load(tmp_path):
    path = tmp_path/'large.h5'
    # Huge sparse datasets: allocating their full arrays would use many GB.
    with h5py.File(path, 'w') as f:
        i = f.create_dataset('s_i', shape=(50000, 40000), dtype='int16', chunks=(64, 64), fillvalue=3)
        q = f.create_dataset('s_q', shape=i.shape, dtype='int16', chunks=(64, 64), fillvalue=-2)
        f['radar_frequency'] = 9.65e9
        i[100:108, 200:208] = 7
    report = inspect_product(path)
    assert report['shape'] == [50000, 40000]
    assert report['complex_samples'] and not report['depth_analysis_ready']
    crop, saved = read_crop(path, [100, 108], [200, 208])
    np.testing.assert_array_equal(crop, 7-2j)
    assert saved['crop']['shape'] == [8, 8]
    with pytest.raises(ValueError, match='memory budget'):
        read_crop(path, [0, 50000], [0, 40000])
    with pytest.raises(ValueError, match='explicit'):
        crop_bounds([100, 100], None, None)


def test_complex_tiff_preserves_phase_and_intensity_tiff_is_rejected(tmp_path):
    complex_path = tmp_path/'complex.tif'
    intensity_path = tmp_path/'amplitude.tif'
    rng = np.random.default_rng(10)
    image = (rng.normal(size=(32, 48))+1j*rng.normal(size=(32, 48))).astype('complex64')
    tifffile.imwrite(complex_path, image)
    tifffile.imwrite(intensity_path, np.abs(image))
    crop, report = read_crop(complex_path, [3, 10], [12, 20])
    np.testing.assert_array_equal(crop, image[3:10, 12:20])
    assert report['complex_samples']
    with pytest.raises(ValueError, match='intensity images'):
        read_crop(intensity_path, [3, 10], [12, 20])
