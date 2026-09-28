import numpy as np

from sarsim.ortho import frame_to_lla, grid_height, predicted_brightness, pyramid_height, register, sample, surface


def test_frame_to_lla_scales():
    lat, lon = frame_to_lla({'latitude': 30.0, 'longitude': 31.0}, np.array([0.0, 1000.0]), np.array([1000.0, 0.0]))
    assert abs((lat[0] - 30.0) * 111_320.0 - 1000.0) < 1e-6
    assert abs((lon[1] - 31.0) * 111_320.0 * np.cos(np.deg2rad(30.0)) - 1000.0) < 1e-6


def test_pyramid_height_apex_and_edge():
    shape = {'centre': [10.0, -5.0, 60.0], 'base': 200.0, 'height': 120.0}
    assert pyramid_height(shape, 10.0, -5.0) == 180.0
    assert abs(pyramid_height(shape, 60.0, -5.0) - 120.0) < 1e-9          # half way to the edge: half the height
    assert pyramid_height(shape, 120.0, -5.0) == -np.inf                   # outside the base


def test_surface_takes_the_higher_of_terrain_and_structure():
    grid = {'x0': -100.0, 'y0': -100.0, 'dx': 100.0, 'nx': 3, 'ny': 3, 'z': [50.0] * 9}
    scene = {'terrain': {'kind': 'grid', **grid},
             'structures': [{'shape': {'type': 'pyramid', 'centre': [0.0, 0.0, 50.0], 'base': 100.0, 'height': 80.0}}]}
    assert grid_height(grid, 12.3, -45.6) == 50.0
    z = surface(scene, np.array([0.0, 90.0]), np.array([0.0, 0.0]))
    assert z[0] == 130.0 and z[1] == 50.0


def test_predicted_brightness_faces_the_radar():
    scene = {'terrain': {'kind': 'flat', 'z': 0.0},
             'structures': [{'shape': {'type': 'pyramid', 'centre': [0.0, 0.0, 0.0], 'base': 100.0, 'height': 60.0}}]}
    x = np.arange(-60, 60.1, 1.0)
    X, Y = np.meshgrid(x, x)
    los = np.array([-0.5, 0.0, np.sqrt(0.75)])                             # the radar sits to the west
    b, _ = predicted_brightness(scene, X, Y, los, 1.0)
    west = b[60, 30]                                                       # (x, y) = (-30, 0): the west face
    east = b[60, 90]
    flat = b[0, 0]
    assert west > flat > east


def test_register_recovers_a_known_offset():
    rng = np.random.default_rng(1)
    base = rng.random((120, 90)) + 0.1
    shifted = np.roll(np.roll(base, 7, axis=0), -4, axis=1)
    r, c = np.meshgrid(np.arange(120.0), np.arange(90.0), indexing='ij')
    (dr, dc), corr = register(shifted, (0, 0), 1, 1, r, c, base, max_shift_px=(30, 30))
    assert (dr, dc) == (7.0, -4.0) and corr > 0.9


def test_sample_is_bilinear_and_nan_outside():
    obs = np.arange(12, dtype=float).reshape(3, 4)
    v = sample(obs, (0, 0), 1, 1, np.array([1.0, 0.5, 10.0]), np.array([1.0, 0.5, 1.0]))
    assert v[0] == obs[1, 1] and abs(v[1] - np.mean(obs[:2, :2])) < 1e-12 and np.isnan(v[2])
    w = sample(obs.repeat(2, 0).repeat(2, 1)[::2, ::2], (0, 0), 2, 2, np.array([0.5]), np.array([0.5]))
    assert w[0] == obs[0, 0]                                               # a 2 x 2 cell's centre sits at pixel 0.5
