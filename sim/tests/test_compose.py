"""The composition model: materials, sites, shapes and the voxeliser."""
import json

import numpy as np
import pytest

from katabasis.compose import Grid, list_sites, load_materials, load_site, voxelise
from katabasis.compose import shapes
from katabasis.compose.site import parse_site


def test_every_material_is_sourced_derived_or_assumed():
    mats = load_materials()
    assert {'air', 'water', 'limestone-mokattam', 'dry-sand-giza'} <= set(mats)
    for m in mats.values():
        for p in (m.vp, m.vs, m.rho):
            assert p.status in {'sourced', 'derived', 'assumed', 'definition'}
            if p.status == 'sourced':
                assert p.source and p.quote


def test_giza_bedrock_is_the_measured_mokattam_range():
    m = load_materials()['limestone-mokattam']
    assert m.vp.range == (2800, 3800) and m.vp.source == 'sharafeldin2019'
    assert m.vs.status == 'derived'


def test_every_site_loads_and_every_feature_has_a_source():
    ids = list_sites()
    assert {'bench-void', 'bench-shafts', 'bench-khafre-claim', 'giza'} <= set(ids)
    for sid in ids:
        s = load_site(sid)
        assert all(f['source'] for f in s.features)


@pytest.mark.parametrize('shape,inside,outside', [
    ({'type': 'box', 'centre': [0, 0, 0], 'size': [2, 4, 6]}, [0.9, 1.9, 2.9], [1.1, 0, 0]),
    ({'type': 'box', 'centre': [0, 0, 0], 'size': [2, 10, 2], 'pitch_deg': 45}, [0, 3, 3], [0, 3, -3]),
    ({'type': 'box', 'centre': [0, 0, 0], 'size': [10, 2, 2], 'yaw_deg': 90}, [0, 4.5, 0], [4.5, 0, 0]),
    ({'type': 'cylinder', 'centre': [0, 0, 0], 'radius': 1, 'height': 4}, [0.7, 0.7, 1.9], [0.8, 0.8, 0]),
    ({'type': 'sphere', 'centre': [1, 1, 1], 'radius': 1}, [1.5, 1.5, 1.5], [0, 0, 0]),
    ({'type': 'prism', 'polygon': [[0, 0], [4, 0], [0, 4]], 'bottom': -1, 'top': 1}, [1, 1, 0], [3, 3, 0]),
    ({'type': 'pyramid', 'centre': [0, 0, 0], 'base': 10, 'height': 10}, [0, 0, 9], [4, 4, 5]),
])
def test_shapes_contain_what_they_should(shape, inside, outside):
    shapes.validate(shape)
    got = shapes.contains(shape, np.array([inside, outside], float))
    assert got.tolist() == [True, False]
    lo, hi = shapes.bounds(shape)
    assert np.all(lo <= np.array(inside)) and np.all(np.array(inside) <= hi)


def test_a_passage_pitched_up_to_the_north_rises_northward():
    box = {'type': 'box', 'centre': [0, 0, 0], 'size': [1, 20, 1], 'pitch_deg': 30}
    assert shapes.contains(box, np.array([[0, 8, 8 * np.tan(np.radians(30))]]))[0]
    assert not shapes.contains(box, np.array([[0, 8, -8 * np.tan(np.radians(30))]]))[0]


def test_voxelised_chamber_has_the_right_volume_and_is_air():
    s = load_site('bench-void')
    g = Grid.covering(tuple(s.extent['x']), tuple(s.extent['y']), (-30.0, 4.0), 0.5)
    m = voxelise(s, g)
    cells = m.cells('chamber')
    assert cells.sum() * g.spacing ** 3 == pytest.approx(216.0, rel=0.02)
    assert np.all(m.vs[cells] == 0) and np.all(m.material[cells] == m.legend.index('air'))
    # sand cover over limestone, air above the ground
    col = m.material[g.shape[0] // 4, g.shape[1] // 4]
    names = [m.legend[c] for c in col]
    iz0 = int(np.argmin(np.abs(g.z - 0.25)))
    assert names[0] == 'air' and names[iz0 + 1] == 'dry-sand-giza' and names[-1] == 'limestone-mokattam'


def test_heterogeneity_perturbs_only_its_stratum_with_the_stated_sigma():
    s = load_site('bench-shafts')
    g = Grid.covering((-60.0, 60.0), (-60.0, 60.0), (-48.0, 4.0), 1.0)
    m = voxelise(s, g)
    mok = (m.material == m.legend.index('limestone-mokattam')) & (m.feature < 0)
    rel = m.vp[mok] / load_materials()['limestone-mokattam'].vp.value - 1
    assert abs(rel.std() - 0.03) < 0.01
    sand = m.material == m.legend.index('dry-sand-giza')
    assert np.allclose(m.vp[sand], load_materials()['dry-sand-giza'].vp.value)


def test_water_table_saturates_what_can_be_saturated():
    raw = json.loads((load_site('bench-void').directory / 'site.json').read_text())
    raw['cover'] = [{'material': 'dry-sand-giza', 'thickness': 10.0}]
    raw['water_table'] = {'z': -5.0}
    s = parse_site(raw)
    g = Grid.covering((-10.0, 10.0), (-10.0, 10.0), (-20.0, 2.0), 1.0)
    m = voxelise(s, g)
    col = [m.legend[c] for c in m.material[5, 5]]
    z = g.z
    assert col[int(np.argmin(np.abs(z + 2.5)))] == 'dry-sand-giza'
    assert col[int(np.argmin(np.abs(z + 7.5)))] == 'wet-sand'


def test_claimed_features_stay_claimed():
    s = load_site('bench-khafre-claim')
    assert {f['status'] for f in s.features} == {'claimed'}
    assert len([f for f in s.features if f['kind'] == 'shaft']) == 8


def test_bad_sites_are_refused():
    raw = json.loads((load_site('bench-void').directory / 'site.json').read_text())
    bad = json.loads(json.dumps(raw))
    bad['features'][0]['source'] = ''
    with pytest.raises(ValueError, match='source'):
        parse_site(bad)
    bad = json.loads(json.dumps(raw))
    bad['strata'][0]['material'] = 'unobtainium'
    with pytest.raises(ValueError):
        parse_site(bad)
