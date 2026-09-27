"""The 3-D elastic solver against exact solutions, and its basic symmetries.

Kept small enough for CI; experiments/p1_01_forward_validation.py runs the
same checks at full size and publishes them.
"""
import numpy as np
import pytest
from scipy.signal import correlate

from katabasis.compose.grid import Grid
from katabasis.seismic.analytic import rayleigh_speed, stokes_displacement
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, ricker

VP, VS, RHO = 2000.0, 2000.0 / np.sqrt(3), 2000.0


def uniform(g, air=None):
    shape = g.shape
    air = np.zeros(shape, bool) if air is None else air
    return Medium.from_arrays(g, np.where(air, 0, VP), np.where(air, 0, VS), np.where(air, 0, RHO), air=air)


def ncc(a, b):
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


@pytest.fixture(scope='module')
def fullspace():
    g = Grid.covering((-60, 60), (-60, 60), (-60, 60), 2.0)
    sim = Simulation(uniform(g), pml_width=10, faces=(True,) * 6, f0=25)
    nt = 330
    t = np.arange(nt) * sim.dt
    w = ricker(t, 25.0) * 1e9
    src = (1.0, 1.0, 1.0)                        # a node: the force is split about it
    recs = np.array([[1.0, 1.0, 33.0], [33.0, 1.0, 1.0], [23.0, 23.0, 23.0]])
    res = sim.run([Source(src, w, 'force', (0, 0, 1))], Receivers(recs), nt)
    return sim, t, w, np.array(src), recs, res


def test_point_force_in_a_whole_space_matches_stokes(fullspace):
    sim, t, w, src, recs, res = fullspace
    for n, r in enumerate(recs):
        u = stokes_displacement(t + sim.dt / 2, w, r - src, 2, VP, VS, RHO)
        v = np.gradient(u, t, axis=1)
        for c in range(3):
            if np.abs(v[c]).max() < 1e-3 * np.abs(v).max():
                continue
            assert ncc(v[c], res.traces[n, c]) > 0.995
            assert np.abs(res.traces[n, c]).max() / np.abs(v[c]).max() == pytest.approx(1.0, abs=0.02)


def test_absorbing_boundaries_let_the_energy_out(fullspace):
    sim, t, w, src, recs, res = fullspace
    e_end_of_source = sim.energy()
    quiet = Source(tuple(src), np.zeros(400), 'force', (0, 0, 1))
    sim.run([quiet], Receivers(recs), 400)
    assert sim.energy() < 2e-3 * e_end_of_source


def test_rayleigh_wave_travels_at_the_rayleigh_speed():
    h = 2.0
    cr = rayleigh_speed(VP, VS)
    fmax = VS / (10 * h)                         # ten points per minimum S wavelength
    f0 = fmax / 2.5
    lam = cr / f0
    g = Grid.covering((-30, 9 * lam + 30), (-40, 40), (-1.5 * lam, 6), h)
    air = np.broadcast_to(g.z[None, None, :] > 0, g.shape).copy()
    sim = Simulation(uniform(g, air), pml_width=14, f0=f0)
    nt = int((9 * lam / cr + 3 / f0) / sim.dt)
    t = np.arange(nt) * sim.dt
    k0 = int(np.argmax(sim.m.solid[0, 0]))
    zs = g.z[k0]
    xs = np.round(np.arange(3, 9) * lam / h) * h + g.origin[0] % h
    res = sim.run([Source((0.0, 0.0, zs - h / 2), ricker(t, f0) * 1e9, 'force', (0, 0, -1))],
                  Receivers(np.array([[x, 0.0, zs] for x in xs])), nt)
    tr = res.traces[:, 2].astype(float)
    up = 8
    fine = lambda s: np.interp(np.arange(0, nt, 1 / up), np.arange(nt), s)
    lags = [(np.argmax(correlate(fine(tr[a]), fine(tr[0]))) - (nt * up - 1)) / up * sim.dt for a in range(1, len(xs))]
    c = 1 / np.polyfit(xs[1:] - xs[0], lags, 1)[0]
    assert c == pytest.approx(cr, rel=0.01)


def test_reciprocity_holds_across_a_void():
    g = Grid.covering((-50, 50), (-50, 50), (-50, 50), 2.0)
    air = np.zeros(g.shape, bool)
    X, Y, Z = np.meshgrid(g.x, g.y, g.z, indexing='ij')
    air[(np.abs(X) < 6) & (np.abs(Y) < 6) & (np.abs(Z) < 6)] = True
    sim = Simulation(uniform(g, air), pml_width=10, faces=(True,) * 6, f0=25)
    nt = 420
    t = np.arange(nt) * sim.dt
    w = ricker(t, 25.0) * 1e9
    A, B = (-17.0, -7.0, -11.0), (19.0, 9.0, 7.0)      # both outside the absorbing layers
    # force along east at A, record everything at B; then force along up at B, record at A
    r1 = sim.run([Source(A, w, 'force', (1, 0, 0))], Receivers(np.array([B])), nt)
    sim.reset()
    r2 = sim.run([Source(B, w, 'force', (0, 0, 1))], Receivers(np.array([A])), nt)
    g_zx = r1.traces[0, 2]          # up velocity at B from east force at A
    g_xz = r2.traces[0, 0]          # east velocity at A from up force at B
    assert ncc(g_zx, g_xz) > 0.999
    assert np.abs(g_zx).max() / np.abs(g_xz).max() == pytest.approx(1.0, abs=0.01)


def test_unstable_time_steps_are_refused():
    g = Grid.covering((-20, 20), (-20, 20), (-20, 20), 2.0)
    with pytest.raises(ValueError, match='unstable'):
        Simulation(uniform(g), dt=1e-3)
