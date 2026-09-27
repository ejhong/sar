"""Full-waveform inversion: whole records fitted by adjoint-state gradients.

The misfit is J = ½ Σ ∫ |v_syn − v_obs|² dt over every shot, receiver and
component of the velocity records. Its gradient with respect to the Lamé
parameters is the zero-lag correlation of forward and adjoint strains
(Tromp, Tape & Liu 2005, Geophys. J. Int. 160, 195):

    ∂J/∂λ(x) = −∫ tr ε(x, t) · tr ε†(x, T − t) dt
    ∂J/∂μ(x) = −∫ 2 ε(x, t) : ε†(x, T − t) dt

In the velocity-stress scheme the stresses of a run are those of the
displacement its force drives, so strains follow from stresses through the
compliance: tr ε = tr σ / (3λ + 2μ) and the deviator ε' = σ' / 2μ. For
velocity data the displacement residual is −∂t(v_syn − v_obs); reversed in
time, it is the adjoint force at each receiver. Chained to the P and S
speeds at fixed density:

    ∂J/∂Vp = 2ρVp ∂J/∂λ        ∂J/∂Vs = 2ρVs (∂J/∂μ − 2 ∂J/∂λ)

Multiscale inversion filters the synthetics and the records alike inside
the misfit (never the source, whose filtered tail would not fit the record):
a symmetric Gaussian kernel with zero padding, so the adjoint source is the
filtered residual filtered once more.

The optimiser is limited-memory BFGS with the illumination of the forward
wavefield as a fixed diagonal preconditioner, the gradient muted within a
few metres of every station (where the kernels are singular), bounds by
clipping, and a backtracking line search. Nothing here sees the true model:
only records.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from ..compose.grid import Grid
from .elastic3d import Medium, Receivers, Simulation, Source

DIRS = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))     # east, north, up: the trace components


def lowpass(x: np.ndarray, sigma_s: float, dt: float) -> np.ndarray:
    """Gaussian smoothing in time (half amplitude at sqrt(ln 2 / 2) / (π sigma_s)); symmetric, so self-adjoint."""
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(x, sigma_s / dt, axis=-1, mode='constant', cval=0.0, truncate=4.0)


def gaussian_sigma(f_half: float) -> float:
    """The Gaussian time width whose response is one half at f_half."""
    return float(np.sqrt(np.log(2) / 2) / (np.pi * f_half))


def shot_gradient(sim: Simulation, source: Source, receivers: np.ndarray, obs: np.ndarray,
                  box: tuple[slice, slice, slice], every: int = 4, sigma_s: float = 0.0):
    """Misfit, ∂J/∂λ, ∂J/∂μ and forward illumination inside `box` for one shot; obs is (nrec, 3, nt),
    already filtered when sigma_s > 0."""
    nt = obs.shape[-1]
    dt = sim.dt
    frames: list[np.ndarray] = []

    def keep(it, s: Simulation):
        frames.append(np.stack([a[box] for a in s.s]))

    sim.reset()
    syn = sim.run([source], Receivers(receivers), nt, on_step=keep, on_step_every=every).traces
    filt = (lambda x: lowpass(x, sigma_s, dt)) if sigma_s > 0 else (lambda x: x)
    r = filt(syn.astype(np.float64)) - obs
    J = 0.5 * float(np.sum(r ** 2)) * dt
    r = filt(r)                                              # the filter's transpose (itself)
    f = np.gradient(r[:, :, ::-1], dt, axis=-1)              # −∂t r, reversed in time
    adj = [Source(tuple(p), np.ascontiguousarray(f[i, c]), 'force', d)
           for i, p in enumerate(receivers) for c, d in enumerate(DIRS)]
    m = sim.m
    lam = m.lam[box].astype(np.float64)
    mu = 0.5 * (m.lam2mu[box].astype(np.float64) - lam)
    k3 = 3 * lam + 2 * mu
    safe_mu = np.where(mu > 0, mu, np.inf)
    safe_k3 = np.where(k3 > 0, k3, np.inf)
    g_lam = np.zeros(lam.shape)
    g_mu = np.zeros(lam.shape)
    illum = np.zeros(lam.shape)
    w = dt * every * sim.grid.spacing ** 3                 # the kernels are densities: × cell volume per cell

    def correlate(it, s: Simulation):
        fs = nt - 1 - it
        if fs % every:
            return
        k = fs // every
        if k >= len(frames):
            return
        sf = frames[k]
        sa = [a[box] for a in s.s]
        trf = sf[0] + sf[1] + sf[2]
        tra = sa[0] + sa[1] + sa[2]
        tt = trf * tra
        dev = ((sf[0] - trf / 3) * (sa[0] - tra / 3) + (sf[1] - trf / 3) * (sa[1] - tra / 3)
               + (sf[2] - trf / 3) * (sa[2] - tra / 3) + 2 * (sf[3] * sa[3] + sf[4] * sa[4] + sf[5] * sa[5]))
        tr_eps = tt / safe_k3 ** 2                         # tr ε · tr ε†
        g_lam[...] -= w * tr_eps
        g_mu[...] -= w * 2 * (dev / (4 * safe_mu ** 2) + tr_eps / 3)
        illum[...] += dt * every * (sf.astype(np.float64) ** 2).sum(0)

    sim.reset()
    sim.run(adj, Receivers(receivers[:1]), nt, on_step=correlate, on_step_every=1)
    return J, g_lam, g_mu, illum, syn


@dataclass
class Problem:
    """Records to fit, over a model whose P and S speeds may change inside `box`."""
    grid: Grid
    rho: np.ndarray                  # density, fixed (full grid)
    air: np.ndarray                  # air cells, fixed (full grid)
    vp0: np.ndarray                  # starting speeds (full grid)
    vs0: np.ndarray
    box: tuple[slice, slice, slice]
    dt: float
    sources: np.ndarray              # (ns, 3) site coordinates
    receivers: np.ndarray            # (nr, 3)
    obs: np.ndarray                  # (ns, nr, 3, nt) velocity records
    wavelet: np.ndarray              # the source force history (vertical)
    f0: float
    pml: int = 16
    every: int = 4
    vp_bounds: tuple[float, float] = (0.3, 1.12)     # relative to the start
    vs_bounds: tuple[float, float] = (0.2, 1.12)
    max_vp_vs: float = 1.3                           # keep the bulk modulus positive: Vs ≤ Vp / 1.3
    sigma_s: float = 0.0                             # Gaussian low-pass of synthetics and records (0: none)
    obs_filtered: bool = False                       # the records were low-passed already (before decimation)
    mute_m: tuple[float, float] = (1.5, 3.0)         # gradient zero within the first radius of a station, whole beyond the second
    log: list = field(default_factory=list)

    def __post_init__(self):
        if self.sigma_s > 0 and not self.obs_filtered:
            self.obs = lowpass(self.obs, self.sigma_s, self.dt)
        g, b = self.grid, self.box
        X, Y, Z = np.meshgrid(g.x[b[0]], g.y[b[1]], g.z[b[2]], indexing='ij')
        d = np.full(X.shape, np.inf)
        for p in np.concatenate([self.sources, self.receivers]):
            d = np.minimum(d, np.sqrt((X - p[0]) ** 2 + (Y - p[1]) ** 2 + (Z - p[2]) ** 2))
        r0, r1 = self.mute_m
        self.mute = np.clip((d - r0) / (r1 - r0), 0, 1)
        self.mute = 0.5 - 0.5 * np.cos(np.pi * self.mute)

    def speeds(self, dvp: np.ndarray, dvs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        vp, vs = self.vp0.copy(), self.vs0.copy()
        vp[self.box] *= 1 + dvp
        vs[self.box] *= 1 + dvs
        vs = np.minimum(vs, vp / self.max_vp_vs)
        return vp, vs

    def medium(self, vp, vs) -> Medium:
        return Medium.from_arrays(self.grid, vp, vs, self.rho, self.air)

    def evaluate(self, dvp: np.ndarray, dvs: np.ndarray):
        """Misfit and its gradient with respect to the relative changes dvp, dvs (box arrays)."""
        vp, vs = self.speeds(dvp, dvs)
        sim = Simulation(self.medium(vp, vs), dt=self.dt, pml_width=self.pml, f0=self.f0, pml_vp=float(self.vp0.max()))
        J = 0.0
        gl = gm = il = None
        for n, s in enumerate(self.sources):
            j, a, b, c, _ = shot_gradient(sim, Source(tuple(s), self.wavelet, 'force', (0, 0, -1)), self.receivers,
                                          self.obs[n], self.box, self.every, self.sigma_s)
            J += j
            gl, gm, il = (a, b, c) if gl is None else (gl + a, gm + b, il + c)
        rho, vpb, vsb = self.rho[self.box], vp[self.box], vs[self.box]
        g_vp = 2 * rho * vpb * gl
        g_vs = 2 * rho * vsb * (gm - 2 * gl)
        # chain to the relative changes: v = v0 (1 + d)
        return J, g_vp * self.vp0[self.box], g_vs * self.vs0[self.box], il


def lbfgs(problem: Problem, iterations: int, first_step: float = 0.05, memory: int = 5,
          precondition_eps: float = 0.02, smooth_cells: float = 1.0, start=None, label: str = '') -> dict:
    """Minimise the misfit over (dvp, dvs). The first step moves the largest cell by `first_step`."""
    from scipy.ndimage import gaussian_filter
    shape = problem.vp0[problem.box].shape
    n = int(np.prod(shape))
    lo = np.concatenate([np.full(n, problem.vp_bounds[0] - 1), np.full(n, problem.vs_bounds[0] - 1)])
    hi = np.concatenate([np.full(n, problem.vp_bounds[1] - 1), np.full(n, problem.vs_bounds[1] - 1)])
    x = np.zeros(2 * n) if start is None else np.concatenate([start[0].ravel(), start[1].ravel()])
    split = lambda v: (v[:n].reshape(shape), v[n:].reshape(shape))
    P = None

    def fg(v):
        nonlocal P
        t0 = time.time()
        J, gp, gs, il = problem.evaluate(*split(v))
        gp, gs = gp * problem.mute, gs * problem.mute
        if P is None:                                   # the preconditioner, fixed from the start
            ilm = gaussian_filter(il, 2.0)
            P = 1.0 / (ilm / ilm.max() + precondition_eps)
            P /= P.max()
        gp = gaussian_filter(gp * P, smooth_cells)
        gs = gaussian_filter(gs * P, smooth_cells)
        return J, np.concatenate([gp.ravel(), gs.ravel()]), time.time() - t0

    J, g, sec = fg(x)
    J0 = J
    hist = [{'iteration': 0, 'misfit': 1.0, 'evaluations': 1, 'seconds': round(sec, 1)}]
    print(f'  {label} iteration 0: misfit 1 ({sec:.0f} s)', flush=True)
    S, Y = [], []
    evals = 1
    for it in range(1, iterations + 1):
        # two-loop recursion
        q = -g.copy()
        alph = []
        for s, y in reversed(list(zip(S, Y))):
            a = np.dot(s, q) / np.dot(y, s)
            alph.append(a)
            q -= a * y
        if S:
            q *= np.dot(S[-1], Y[-1]) / np.dot(Y[-1], Y[-1])
        for (s, y), a in zip(zip(S, Y), reversed(alph)):
            b = np.dot(y, q) / np.dot(y, s)
            q += (a - b) * s
        d = q
        if np.dot(d, g) >= 0:                           # not a descent direction: restart
            S, Y = [], []
            d = -g
        step = first_step / np.abs(d).max() if not S else 1.0
        accepted = False
        for _ in range(4):
            xn = np.clip(x + step * d, lo, hi)
            Jn, gn, sec = fg(xn)
            evals += 1
            if Jn < J:
                accepted = True
                break
            step *= 0.4
            print(f'  {label} iteration {it}: step rejected (misfit {Jn / J0:.4f}), shrinking', flush=True)
        if not accepted:
            print(f'  {label} iteration {it}: no decrease; stopping', flush=True)
            break
        s, y = xn - x, gn - g
        if np.dot(s, y) > 1e-12 * np.dot(y, y):
            S.append(s)
            Y.append(y)
            S, Y = S[-memory:], Y[-memory:]
        x, J, g = xn, Jn, gn
        hist.append({'iteration': it, 'misfit': J / J0, 'evaluations': evals, 'seconds': round(sec, 1)})
        print(f'  {label} iteration {it}: misfit {J / J0:.4f} ({evals} evaluations, {sec:.0f} s each)', flush=True)
    dvp, dvs = split(x)
    return {'dvp': dvp, 'dvs': dvs, 'history': hist, 'misfit_start': J0, 'misfit_end': J}
