"""Static response by relaxation: the elastic solver run until the rock is at rest.

A wave far longer than a chamber and its depth strains the rock around it almost uniformly and almost
statically, so the chamber's imprint on the surface is the static answer to one question: what does opening
a void do to rock under this uniform stress? The solver answers it directly. Start from the uniform stress,
with the chamber already open so its walls carry traction nothing balances; a horizontal stress
(s_zz = s_xz = s_yz = 0) is an exact equilibrium of intact rock under a free surface, so the chamber is the
only thing out of balance. A light global damping of the velocity lets the waves this releases die away
without changing the equilibrium they die away to, and the velocity at the surface is integrated until the
ground is still: the displacement made is the imprint.

The convolutional PML leaves a uniform stress alone (it acts on spatial derivatives), and at zero frequency it
stretches the domain's edge towards infinity (its coordinate stretch 1 + d / alpha grows without bound where
alpha falls to zero at the outer edge), so the answer approaches the half-space's; the experiment checks it
against domain size and grid spacing.
"""
from __future__ import annotations

import math

import numpy as np

from .elastic3d import Medium, Receivers, Simulation


def relax(medium: Medium, stress: dict, receivers: np.ndarray, duration: float, damping: float = 30.0,
          pml_width: int = 12, f0: float = 5.0, checkpoints: int = 10):
    """Displacement (n, 3; east, north, up) at `receivers` when the rock, loaded with the uniform horizontal
    `stress` {'xx', 'yy', 'xy'} (Pa) in every solid cell, comes to rest. Returns the displacement and its value at
    `checkpoints` evenly spaced times, to judge convergence."""
    sim = Simulation(medium, pml_width=pml_width, f0=f0)
    sxx, syy, szz, sxy, sxz, syz = sim.s
    solid = medium.lam2mu > 0
    sxx[...] = np.where(solid, stress.get('xx', 0.0), 0.0)
    syy[...] = np.where(solid, stress.get('yy', 0.0), 0.0)
    sxy[...] = np.where(medium.mu_xy > 0, stress.get('xy', 0.0), 0.0)
    nt = int(math.ceil(duration / sim.dt))
    decay = np.float32(math.exp(-damping * sim.dt))

    def damp(it, s):
        for v in s.v:
            v *= decay
    res = sim.run([], Receivers(np.asarray(receivers, float)), nt, on_step=damp)
    u = np.cumsum(res.traces.astype(np.float64), axis=2) * sim.dt          # (n, 3, nt)
    marks = np.linspace(nt // checkpoints, nt, checkpoints).astype(int) - 1
    return u[:, :, -1], {'t': (marks + 1) * sim.dt, 'u': u[:, :, marks]}, sim.dt, nt
