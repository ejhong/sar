"""The lab's command places its lines on a synthetic product where the product's own geometry puts the site."""
import json

import h5py
import numpy as np

from katabasis.lab import flat_project, synthetic_mapping
from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup
from sarsim.slcfile import write_synthetic_slc


def test_lab_places_site_points_on_a_synthetic_product(tmp_path):
    g = DwellGeometry.from_record('giza-20220715')
    st = EchoSetup.from_geometry(g)
    lat, lon = g.source['scene_centre']['latitude'], g.source['scene_centre']['longitude']
    rows, cols = 64, 32
    path = write_synthetic_slc(tmp_path / 'SYNTHETIC_test.h5', np.zeros((rows, cols), np.complex64), st, g.Ka_signed,
                               lat, lon, 75.0, g.heading_deg, g.theta_deg)
    with h5py.File(path, 'a') as h:
        h['synthetic'].attrs['site_to_pixel'] = json.dumps({'centre_row': rows // 2, 'centre_col': cols // 2, 'dx': g.dx,
                                                             'dr': g.dr, 'theta_deg': g.theta_deg, 'heading_deg': g.heading_deg})
    m = synthetic_mapping(path)
    assert abs(m['height_m'] - 75.0) < 1e-3
    assert abs(m['origin']['latitude'] - lat) < 1e-12 and abs(m['origin']['longitude'] - lon) < 1e-12
    x, y = np.array([10.0, -5.0, 0.0]), np.array([3.0, 20.0, 0.0])
    r, c = flat_project(m, x, y, np.zeros(3))
    a = x * g.along_track_en[0] + y * g.along_track_en[1]
    gr = x * g.ground_range_en[0] + y * g.ground_range_en[1]
    np.testing.assert_allclose(r, rows // 2 + a / g.dx, atol=1e-3)
    np.testing.assert_allclose(c, cols // 2 + gr * np.sin(g.theta) / g.dr, atol=1e-3)


def test_a_real_product_has_no_synthetic_mapping(tmp_path):
    path = tmp_path / 'plain.h5'
    with h5py.File(path, 'w') as h:
        h['s_i'] = np.zeros((2, 2), np.float32)
    assert synthetic_mapping(path) is None
