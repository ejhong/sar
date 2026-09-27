"""Peterson's new low- and high-noise models (NLNM, NHNM).

Source: Peterson, J.R., 1993, Observations and modeling of seismic background
noise: U.S. Geological Survey Open-File Report 93-322, 94 p.,
https://doi.org/10.3133/ofr93322. Tables 3 and 4 (report pp. 33-34), "Line
parameters for constructing the NLNM [NHNM] curve given the period (P)",
with "NLNM_acc = A + B log10(P) dB referred to 1 (m/sec^2)^2/Hz" and
"NLNM_vel = NLNM_acc + 20.0 log10(P/2π)". The NLNM is "a hypothetical
background spectrum that is unlikely to be duplicated at any single location
on Earth"; the NHNM is "a spectrum of average high background noise power in
the network" (p. 32). Both are defined for periods 0.1 to 100,000 s, i.e. up
to 10 Hz: they say nothing about the band above.
"""
from __future__ import annotations

import numpy as np

# (period from which the segment applies, A, B), transcribed from Tables 3 and 4
NLNM = [(0.10, -162.36, 5.64), (0.17, -166.7, 0.00), (0.40, -170.00, -8.30), (0.80, -166.40, 28.90),
        (1.24, -168.60, 52.48), (2.40, -159.98, 29.81), (4.30, -141.10, 0.00), (5.00, -71.36, -99.77),
        (6.00, -97.26, -66.49), (10.00, -132.18, -31.57), (12.00, -205.27, 36.16), (15.60, -37.65, -104.33),
        (21.90, -114.37, -47.10), (31.60, -160.58, -16.28), (45.00, -187.50, 0.00), (70.00, -216.47, 15.70),
        (101.00, -185.00, 0.00), (154.00, -168.34, -7.61), (328.00, -217.43, 11.90), (600.00, -258.28, 26.60),
        (10000.00, -346.88, 48.75)]
NHNM = [(0.10, -108.73, -17.23), (0.22, -150.34, -80.50), (0.32, -122.31, -23.87), (0.80, -116.85, 32.51),
        (3.80, -108.48, 18.08), (4.60, -74.66, -32.95), (6.30, 0.66, -127.18), (7.90, -93.37, -22.42),
        (15.40, 73.54, -162.98), (20.00, -151.52, 10.01), (354.80, -206.66, 31.63)]
P_MIN, P_MAX = 0.10, 100000.0


def acceleration_db(period: np.ndarray, model: str = 'low') -> np.ndarray:
    """Acceleration power, dB re 1 (m/s^2)^2/Hz, at the given periods (s); NaN outside 0.1-100,000 s."""
    table = NLNM if model == 'low' else NHNM
    p = np.atleast_1d(np.asarray(period, float))
    out = np.full(p.shape, np.nan)
    starts = np.array([t[0] for t in table])
    ok = (p >= P_MIN) & (p <= P_MAX)
    k = np.searchsorted(starts, p[ok], side='right') - 1
    A = np.array([t[1] for t in table])[k]
    B = np.array([t[2] for t in table])[k]
    out[ok] = A + B * np.log10(p[ok])
    return out


def velocity_psd(freq: np.ndarray, model: str = 'low') -> np.ndarray:
    """Velocity power spectral density, (m/s)^2/Hz, at the given frequencies."""
    f = np.asarray(freq, float)
    return 10 ** (acceleration_db(1 / f, model) / 10) / (2 * np.pi * f) ** 2


def band_rms_velocity(f_lo: float, f_hi: float, model: str = 'low', n: int = 2000) -> float:
    """RMS ground velocity (m/s) in a band, from the model's velocity PSD."""
    f = np.geomspace(f_lo, f_hi, n)
    return float(np.sqrt(np.trapezoid(velocity_psd(f, model), f)))
