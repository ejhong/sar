"""First-arrival picking, as a processor would do it on crosshole or refraction records.

Noise is added first (a fraction of each record's peak). Each record's first
energy is found where its three-component envelope first rises above both
the noise measured before the shot and a small fraction of the record's
peak. The first arrival is then timed by waveform picking: each component
is cross-correlated with the source's far-field P pulse (the wavelet's time
derivative, which is what a velocity sensor records from a point force or an
explosion), over lags from a period before that first energy to a quarter
period after it; the best match gives the travel time directly, free of the
onset bias a threshold picker has on a smooth wavelet.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import hilbert


def add_noise(traces: np.ndarray, snr_db: float, rng: np.random.Generator) -> np.ndarray:
    """Gaussian noise at a set signal-to-noise ratio against each record's peak (vector) amplitude."""
    peak = np.sqrt((traces.astype(np.float64) ** 2).sum(1)).max(-1)          # (nrec,)
    sigma = peak / 10 ** (snr_db / 20)
    return traces + rng.standard_normal(traces.shape) * sigma[:, None, None]


def envelope(trace: np.ndarray) -> np.ndarray:
    """Three-component envelope of a (components, nt) record (tapered against wrap-around)."""
    x = trace.astype(np.float64)
    n = x.shape[-1]
    pad = np.zeros(x.shape[:-1] + (n,))
    a = np.abs(hilbert(np.concatenate([x, pad], -1), axis=-1))[..., :n]
    return np.sqrt((a ** 2).sum(0))


def first_energy(trace: np.ndarray, quiet: int, fraction: float = 0.04, noise_factor: float = 6.0) -> int:
    env = envelope(trace)
    floor = noise_factor * np.sqrt(np.mean(env[:max(quiet, 4)] ** 2))
    thr = max(floor, fraction * env.max())
    above = np.nonzero(env > thr)[0]
    return int(above[0]) if len(above) else -1


def pick(trace: np.ndarray, dt: float, pulse: np.ndarray, period: float, quiet: int) -> float:
    """Travel time of one (components, nt) record by correlation with the reference pulse.

    `pulse` is the far-field P pulse sampled on the same dt, starting at the
    shot; the returned time is the lag at which it best matches.
    """
    i0 = first_energy(trace, quiet)
    if i0 < 0:
        return np.nan
    n = trace.shape[-1]
    # the pulse's own front, by the same rule: a lag L puts that front at L + front
    front = first_energy(pulse[None, :], 4, noise_factor=0.0)
    w = int(period / dt)
    Ls = np.arange(i0 - front - w, i0 - front + w // 2 + 1)
    best, bestL = -1.0, 0
    x = trace.astype(np.float64)
    m = len(pulse)
    for L in Ls:
        a0, a1 = max(0, L), min(n, L + m)
        if a1 - a0 < m // 3:
            continue
        seg = x[:, a0:a1]
        p = pulse[a0 - L:a1 - L]
        c = np.abs(seg @ p) / (np.sqrt((seg ** 2).sum(1) * (p @ p)) + 1e-30)
        cm = float(c.max())
        if cm > best:
            best, bestL = cm, L
    return bestL * dt


def far_field_pulse(wavelet: np.ndarray, dt: float, period: float, lobes: float = 0.75) -> np.ndarray:
    """The wavelet's derivative from the shot up to `lobes` periods past its front.

    Only the leading part of the pulse is matched: what follows it on a real
    record is soon mixed with ghosts from the free surface and with S waves.
    """
    d = np.gradient(wavelet.astype(np.float64), dt)
    front = first_energy(d[None, :], 4, noise_factor=0.0)
    return d[: front + int(lobes * period / dt) + 1]


def pick_all(traces: np.ndarray, t0: float, dt: float, wavelet: np.ndarray, period: float,
             quiet_s: float | None = None) -> np.ndarray:
    """Travel times for every record (components, nt), `t0` the time of the first sample."""
    pulse = far_field_pulse(wavelet, dt, period)
    quiet = int((quiet_s if quiet_s is not None else 0.25 * period) / dt)
    return np.array([pick(tr, dt, pulse, period, quiet) for tr in traces]) + t0
