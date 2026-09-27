"""Controlled complex radar measurement and an optimistic phase readout.

Straight-line stationary-phase approximation, known resolved point reflectors,
rectangular processed band, no autofocus, atmospheric phase or range migration.
Each reflector has independent circular complex receiver noise before focusing.
Focusing is an invertible FFT. A known range-basis separates the reflectors on
readout; this is an optimistic benchmark, not a distributed-ground SAR method.
"""
from dataclasses import dataclass
import numpy as np
from .geometry import Geometry


@dataclass
class ReflectorAcquisition:
    geometry: Geometry
    shape: tuple = (1024, 640)
    bearing_deg: float = 45.

    def __post_init__(self):
        g = self.geometry
        self.nu = np.fft.fftfreq(self.shape[0], g.dx)
        self.bins = np.flatnonzero(np.abs(self.nu) <= g.nu_band/2)
        self.bins = self.bins[np.argsort(self.nu[self.bins])]
        self.time = g.nu_to_time(self.nu[self.bins])
        # SH polarization is horizontal, perpendicular to the surface transect.
        self.projection = np.sin(g.theta)*np.cos(np.deg2rad(self.bearing_deg))

    def pixel_locations(self, positions):
        angle = np.deg2rad(self.bearing_deg)
        x = np.asarray(positions)*np.cos(angle)
        r = np.asarray(positions)*np.sin(angle)*np.sin(self.geometry.theta)
        return (self.shape[0]//2+x/self.geometry.dx,
                self.shape[1]//2+r/self.geometry.dr)

    def histories(self, displacement, snr_db, rng):
        """Ground SH displacement [slow time, reflector] -> independent channels."""
        if displacement.shape[0] != len(self.time):
            raise ValueError("displacement must use the acquisition slow-time axis")
        signal = np.exp(-4j*np.pi*self.projection*displacement/self.geometry.lam)
        if np.isfinite(snr_db):
            sigma = 10**(-snr_db/20)/np.sqrt(2)
            signal += sigma*(rng.normal(size=signal.shape)+1j*rng.normal(size=signal.shape))
        return signal

    def range_basis(self, positions):
        _, cols = self.pixel_locations(positions)
        freq = np.fft.fftfreq(self.shape[1], self.geometry.dr)
        keep = np.abs(freq) <= self.geometry.kr_band/2
        spectrum = np.exp(-2j*np.pi*np.outer(cols*self.geometry.dr, freq))*keep
        return np.fft.ifft(spectrum, axis=1).T

    def focus(self, channels, positions):
        rows, _ = self.pixel_locations(positions)
        carrier = np.exp(-2j*np.pi*np.outer(self.nu[self.bins], rows*self.geometry.dx))
        spectrum = np.zeros(self.shape, complex)
        spectrum[self.bins] = (channels*carrier) @ self.range_basis(positions).T
        return np.fft.ifft(spectrum, axis=0)

    def separate(self, slc, positions):
        """Known resolved range supports: deliberately favorable measurement."""
        rows, _ = self.pixel_locations(positions)
        basis = self.range_basis(positions)
        # Range supports must actually be separable. No assumed underground depth.
        if np.linalg.cond(basis) > 1e6:
            raise ValueError("reflectors are not resolved by this range model")
        mixed = np.fft.fft(slc, axis=0)[self.bins]
        channels = mixed @ np.linalg.pinv(basis.T)
        return channels*np.exp(2j*np.pi*np.outer(self.nu[self.bins], rows*self.geometry.dx))

    def read_phase(self, channels):
        phase = np.unwrap(np.angle(channels), axis=0)
        displacement = -phase*self.geometry.lam/(4*np.pi*self.projection)
        # A constant phase/reference range is unknown; do not use a static image.
        return displacement-displacement.mean(axis=0, keepdims=True)


def centred(values):
    return np.asarray(values)-np.mean(values, axis=-2, keepdims=True)


class WaveDictionary:
    """Profile a positive unknown source gain for each candidate physical model.

    No trajectory steering formula appears in this estimator. Dictionary arrays
    are coarse-grid wave predictions; evaluation uses off-grid depths/finer mesh.
    """
    def __init__(self, waves, depths, speeds):
        self.depths = np.asarray(depths, float)
        self.speeds = np.asarray(speeds, float)
        matrix = centred(np.asarray(waves)).reshape(len(waves), -1)
        norms = np.linalg.norm(matrix, axis=1)
        if np.any(norms == 0):
            raise ValueError("a dictionary trace is empty")
        self.unit = matrix / norms[:, None]
        self.norms = norms
        self.null = self.depths == 0
        if not self.null.any():
            raise ValueError("include a cavity-free candidate")

    def predict(self, observation):
        observation = np.asarray(observation)
        single = observation.ndim == 2
        if single:
            observation = observation[None]
        y = centred(observation).reshape(len(observation), -1)
        explained = np.maximum(y @ self.unit.T, 0)**2
        best = explained.argmax(axis=1)
        null = explained[:, self.null].max(axis=1)
        score = np.maximum(0., explained.max(axis=1)-null)
        residual = np.maximum(0., (y*y).sum(axis=1)[:, None]-explained)
        return {"depth": self.depths[best], "speed": self.speeds[best],
                "improvement": score, "residual": residual, "index": best}


def wilson(successes, trials):
    """95% Wilson interval for independent simulated noise trials."""
    if trials == 0:
        return [None, None]
    z = 1.959963984540054
    p = successes/trials
    den = 1+z*z/trials
    centre = (p+z*z/(2*trials))/den
    half = z*np.sqrt(p*(1-p)/trials+z*z/(4*trials*trials))/den
    return [max(0., float(centre-half)), min(1., float(centre+half))]
