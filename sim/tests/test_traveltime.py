"""Travel-time tomography: rays and times in uniform rock."""
import numpy as np

from katabasis.compose import Grid
from katabasis.seismic.traveltime import TomoGrid, forward


def test_rays_between_stations_on_the_ground_have_the_straight_length():
    g = Grid.covering((-18, 18), (-18, 18), (-20, 0), 1.0)
    tg = TomoGrid(g, np.ones(g.shape, bool))
    src = np.array([[-12.0, 0.0, -0.25], [4.0, -8.0, -0.25]])
    rec = np.array([[12.0, 0.0, -0.25], [-8.0, 8.0, -0.25], [16.0, 16.0, -0.25]])
    pairs = np.array([[s, r] for s in range(2) for r in range(3)])
    v = 3300.0
    t, G = forward(np.full(g.shape, 1 / v), tg, src, rec, pairs)
    d = np.linalg.norm(src[pairs[:, 0]] - rec[pairs[:, 1]], axis=1)
    lengths = np.asarray(G.sum(1)).ravel()
    assert np.all(lengths > 0), 'every ray must reach its source'
    assert np.allclose(lengths, d, rtol=0.04)
    assert np.allclose(t, d / v, rtol=0.03)


def test_rays_between_stations_on_the_ground_run_in_the_top_cells():
    g = Grid.covering((-18, 18), (-18, 18), (-20, 0), 1.0)
    tg = TomoGrid(g, np.ones(g.shape, bool))
    src = np.array([[-12.0, 0.0, -0.25]])
    rec = np.array([[12.0, 0.0, -0.25]])
    _, G = forward(np.full(g.shape, 1 / 3300.0), tg, src, rec, np.array([[0, 0]]))
    k = np.unravel_index(G.indices, g.shape)[2]
    w = G.data
    assert np.sum(w[k == 0]) / np.sum(w) > 0.95, 'the ray the times come from must be the ray the matrix holds'
