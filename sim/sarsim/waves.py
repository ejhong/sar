"""A deliberately limited physical forward model for controlled depth benchmarks.

2-D anti-plane shear displacement in homogeneous rock: u_tt = c² div(grad u).
Free-traction faces at the surface and a circular cavity. The hole represents an
infinite horizontal cylinder, not a vertical column. A Ricker body force has
known timing and position. No intrinsic attenuation, geology or ambient source
model is supplied. Absolute force is arbitrary; a *shared cavity-free reference*
sets the displacement scale before targets are evaluated.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import json
import numpy as np


@dataclass(frozen=True)
class WaveModel:
    depth: float | None = None
    radius: float = 20.0
    frequency: float = 25.0
    speed: float = 2000.0
    spacing: float = 5.0
    duration: float = 0.85
    source_x: float = -150.0
    source_z: float = 10.0
    receiver_step: float = 15.0


def surface_wave(model=WaveModel(), cache=None):
    """Return time, surface positions, displacement; cache includes source hash."""
    if model.depth is not None and model.depth <= model.radius:
        raise ValueError("the cavity must be wholly below the surface")
    if min(model.speed, model.spacing, model.frequency, model.duration) <= 0:
        raise ValueError("wave parameters must be positive")
    key = hashlib.sha256(Path(__file__).read_bytes()+json.dumps(asdict(model), sort_keys=True).encode()).hexdigest()[:20]
    path = Path(cache) / f"wave_{key}.npz" if cache else None
    if path and path.exists():
        with np.load(path) as stored:
            return stored['time'], stored['positions'], stored['displacement']
    dx, c = model.spacing, model.speed
    x = np.arange(-800., 800.+dx/2, dx)
    z = np.arange(0., 800.+dx/2, dx)
    dt = .4*dx/c
    time = np.arange(int(np.ceil(model.duration/dt)))*dt
    X, Z = np.meshgrid(x, z)
    rock = np.ones(X.shape, bool)
    if model.depth is not None:
        rock = X*X + (Z-model.depth)**2 > model.radius**2
    faces_x = rock[:, 1:] & rock[:, :-1]
    faces_z = rock[1:] & rock[:-1]
    edge = np.maximum(np.maximum((np.abs(X)-650)/150, 0), (Z-650)/150)
    damping = 250*np.maximum(edge, 0)**2
    source = (int(round(model.source_z/dx)), int(round((model.source_x-x[0])/dx)))
    if not rock[source]:
        raise ValueError("source is inside the void")
    requested = np.arange(-150., 150.+model.receiver_step/2, model.receiver_step)
    receivers = np.rint((requested-x[0])/dx).astype(int)
    old = np.zeros(X.shape)
    field = np.zeros(X.shape)
    trace = np.zeros((len(time), len(receivers)))
    coefficient = (c*dt/dx)**2
    for k, t in enumerate(time):
        lap = np.zeros_like(field)
        flux_x = (field[:, 1:]-field[:, :-1])*faces_x
        flux_z = (field[1:]-field[:-1])*faces_z
        lap[:, :-1] += flux_x
        lap[:, 1:] -= flux_x
        lap[:-1] += flux_z
        lap[1:] -= flux_z
        new = (2*field-(1-damping*dt/2)*old+coefficient*lap)/(1+damping*dt/2)
        a = np.pi*model.frequency*(t-1.5/model.frequency)
        new[source] += dt**2/dx**2*(1-2*a*a)*np.exp(-a*a)
        new *= rock
        old, field = field, new
        trace[k] = field[0, receivers]
    if not np.isfinite(trace).all():
        raise FloatingPointError("nonfinite wavefield")
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, time=time, positions=x[receivers], displacement=trace)
    return time, x[receivers], trace


def sample_wave(time, displacement, sample_time, start=-.45):
    """Place the active event within an aperture; taper its finite record end.

    This specified acquisition operator is shared by calibration and evaluation.
    It is not a model of ongoing ambient excitation.
    """
    taper = np.minimum(1., np.maximum(0., (time[-1]-time)/.06))
    return np.column_stack([np.interp(sample_time-start, time, d*taper, left=0., right=0.)
                            for d in displacement.T])
