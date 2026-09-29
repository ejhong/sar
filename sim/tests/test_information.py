"""The information one image of speckle holds about a moving imprint (sarsim.information, P2-25)."""
import numpy as np
import pytest

from sarsim import information as inf
from sarsim.acquisition import DwellGeometry


@pytest.fixture(scope='module')
def g():
    return DwellGeometry.from_record('giza-20250827')


def bump(g, shape, L, amp=1e-3):
    na, nr = shape
    xa = (np.arange(na) - na // 2) * g.dx
    xg = (np.arange(nr) - nr // 2) * g.dr / np.sin(g.theta)
    Xa, Xg = np.meshgrid(xa, xg, indexing='ij')
    return amp * np.exp(-(Xa ** 2 + Xg ** 2) / (2 * L ** 2))


def test_shared_motion_holds_nothing(g):
    shape = (64, 12)
    tau = lambda t: np.cos(2 * np.pi * 3.0 * t)
    K = np.full(shape, 0.4)
    assert inf.fisher_exact([(K, tau)], g, shape) == 0.0
    img = inf.speckle_image(g, shape, np.random.default_rng(1), K, tau, 1.0)
    assert abs(inf.score(img, K, tau, g)) < 1e-9


def test_sum_matches_covariance(g):
    """The double sum over frequency bins equals tr(C+ dC C+ dC) built scatterer by scatterer."""
    rng = np.random.default_rng(2)
    shape = (16, 6)
    na, nr = shape
    N = na * nr
    in_a, in_r, nu = inf.band_masks(g, shape)
    kr = np.fft.fftfreq(nr, d=g.dr)
    band = (in_a[:, None] & in_r[None, :]).ravel()
    t = g.nu_to_time(nu)
    K = rng.standard_normal(shape)
    tau = lambda tt: np.cos(2 * np.pi * 1500.0 * tt + 0.2)
    Xa, Xr = np.meshgrid(np.arange(na) * g.dx, np.arange(nr) * g.dr, indexing='ij')
    FA, FR = np.meshgrid(nu, kr, indexing='ij')
    E = np.exp(-2j * np.pi * (np.stack([FA.ravel(), FR.ravel()], 1) @ np.stack([Xa.ravel(), Xr.ravel()], 1).T)) / np.sqrt(N)
    Phi = np.repeat(np.stack([(K * tau(t[i])).ravel() for i in range(na)]), nr, axis=0)
    B = E * band[:, None]
    dB = -1j * Phi * B
    C0, dC = B @ B.conj().T, dB @ B.conj().T + B @ dB.conj().T
    Cp = np.linalg.pinv(C0, rcond=1e-10, hermitian=True)
    F_cov = np.trace(Cp @ dC @ Cp @ dC).real
    assert inf.fisher_exact([(K, tau)], g, shape) == pytest.approx(F_cov, rel=1e-9)


def test_continuum_matches_sum(g):
    """The continuum formula (with each row's zero-Doppler time) against the exact sum on the image's own grid."""
    shape = (512, 24)
    L, fm = 0.8, 2.0
    K = bump(g, shape, L)
    k0 = 4 * np.pi / g.lam
    xa = ((np.arange(shape[0]) - shape[0] // 2) * g.dx)[:, None]
    Fe = 0.0
    for psi in (0.0, np.pi / 2):
        al = 2 * np.pi * fm * xa / g.V + psi
        Fe += 0.5 * inf.fisher_exact([(k0 * K * np.cos(al), lambda t: np.cos(2 * np.pi * fm * t)),
                                      (-k0 * K * np.sin(al), lambda t: np.sin(2 * np.pi * fm * t))], g, shape)
    Fc = inf.fisher_sinusoid(K, (g.dx, g.dr / np.sin(g.theta)), fm, g)
    assert Fc == pytest.approx(Fe, rel=0.05)


def test_bound_and_slow_limit(g):
    shape = (512, 48)
    sp = (g.dx, g.dr / np.sin(g.theta))
    for L, fm in ((0.8, 2.0), (0.8, 20.0), (3.0, 0.25)):
        K = bump(g, shape, L)
        assert inf.fisher_sinusoid(K, sp, fm, g) <= inf.fisher_bound(K, sp, g)
    K = bump(g, shape, 3.0)
    assert inf.fisher_sinusoid(K, sp, 0.25, g) == pytest.approx(inf.fisher_slow(2 * np.pi * 0.25 * K, sp, g), rel=0.05)
    # the slow limit depends on the velocity's gradient alone: the same peak velocity, twice as wide, the same F
    K2 = bump(g, (1024, 96), 6.0)
    assert inf.fisher_slow(K2, sp, g) == pytest.approx(inf.fisher_slow(K, sp, g), rel=0.02)


def test_score_attains_fisher(g):
    """The locally most powerful detector scatters by sqrt(F) without the imprint and moves by a F with it."""
    shape = (256, 24)
    k0 = 4 * np.pi / g.lam
    Kp = k0 * bump(g, shape, 0.8)
    tau = lambda t: np.cos(2 * np.pi * 2.0 * t)
    F = inf.fisher_exact([(Kp, tau)], g, shape)
    a = 1.2 / np.sqrt(F)
    assert a * Kp.max() < 0.3
    rng = np.random.default_rng(3)
    n = 1000
    T0 = np.array([inf.score(inf.speckle_image(g, shape, rng), Kp, tau, g) for _ in range(n)])
    T1 = np.array([inf.score(inf.speckle_image(g, shape, rng, Kp, tau, a), Kp, tau, g) for _ in range(n)])
    assert T0.var() / F == pytest.approx(1.0, abs=4 * np.sqrt(2 / n))
    assert (T1.mean() - T0.mean()) / (a * F) == pytest.approx(1.0, abs=4 * np.sqrt(2 / n) / (a * np.sqrt(F)))


def test_point_and_edge():
    t = np.linspace(-12, 12, 4001)
    phase = 1e-3 * np.cos(2 * np.pi * 0.5 * t)
    assert inf.fisher_point(100.0, phase, t) == pytest.approx(100.0 * 1e-6, rel=0.02)
    assert inf.fisher_point(100.0, 1e-3 + 2e-4 * t, t) == pytest.approx(0.0, abs=1e-20)   # phase and position absorb it
    assert inf.edge(4e-18) == pytest.approx(1e-9)
