"""Imaging the scattered wave: reverse-time migration, the first step of full-waveform inversion.

For each shot, the source wavefield is computed forward through the
background model and stored inside an imaging box; the scattered records
(with the target minus without) are injected, time-reversed, at the
receivers and propagated back through the same background. Where the two
wavefields meet in time, the image lights up:

    I(x) = Σ_shots Σ_t v_s(x, t) · v_r(x, t)  /  Σ_shots Σ_t |v_s(x, t)|²

(the velocity cross-correlation imaging condition with source-illumination
compensation). With the background exact, this is the best case for any
waveform method: the gradient FWI would take first.
"""
from __future__ import annotations

import numpy as np

from .elastic3d import Receivers, Simulation, Source


class BoxStore:
    """Velocity snapshots inside an index box, every `every` steps."""

    def __init__(self, box: tuple[slice, slice, slice], nt: int, every: int):
        self.box = box
        self.every = every
        self.frames: list[np.ndarray] = []

    def __call__(self, it, sim: Simulation):
        b = self.box
        vx, vy, vz = sim.v
        self.frames.append(np.stack([vx[b], vy[b], vz[b]]).astype(np.float16))


def rtm_shot(sim: Simulation, source: Source, receivers: np.ndarray, scattered: np.ndarray, nt: int,
             box: tuple[slice, slice, slice], every: int = 3, data_scale: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Image and illumination for one shot. `scattered` is (nrec, 3, nt) in site axes (east, north, up).

    The adjoint sources are the data times `data_scale` (by default, scaled so this shot's largest sample
    is 1e9); pass one scale for every call whose images must add or compare linearly.
    """
    sim.reset()
    store = BoxStore(box, nt, every)
    sim.run([source], Receivers(receivers[:1]), nt, on_step=store, on_step_every=every)
    frames = store.frames
    # adjoint: time-reversed scattered records as forces at the receivers, one per component
    adj = []
    rev = scattered[:, :, ::-1].astype(np.float64)
    scale = data_scale if data_scale is not None else 1e9 / max(np.abs(rev).max(), 1e-30)
    for r, tr in zip(receivers, rev):
        for c, d in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
            adj.append(Source(tuple(r), np.ascontiguousarray(tr[c] * scale), 'force', d))
    sim.reset()
    shape = frames[0].shape[1:]
    image = np.zeros(shape, np.float64)
    illum = np.zeros(shape, np.float64)
    nframes = len(frames)

    def correlate(it, s: Simulation):
        # adjoint step it corresponds to forward step nt-1-it; forward frames were stored at multiples of `every`
        f = nt - 1 - it
        if f % every:
            return
        k = f // every
        if k >= nframes:
            return
        vs = frames[k].astype(np.float32)
        b = box
        vr = np.stack([s.v[0][b], s.v[1][b], s.v[2][b]])
        image[...] += (vs * vr).sum(0)
        illum[...] += (vs * vs).sum(0)

    sim.run(adj, Receivers(receivers[:1]), nt, on_step=correlate, on_step_every=1)
    return image, illum


def laplace_filter(img: np.ndarray, h: float) -> np.ndarray:
    """Suppress the low-wavenumber backscatter RTM leaves along ray paths."""
    from scipy.ndimage import laplace
    return -laplace(img) / h ** 2
