"""Terrain: the ground surface height z(x, y) in site coordinates."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass
class Terrain:
    kind: str
    x: np.ndarray | None = None
    y: np.ndarray | None = None
    z: np.ndarray | None = None          # heights on (x, y), shape (nx, ny)
    flat_z: float = 0.0
    meta: dict | None = None

    def height(self, x, y) -> np.ndarray:
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        if self.kind == 'flat':
            return np.full(np.broadcast(x, y).shape, self.flat_z)
        f = RegularGridInterpolator((self.x, self.y), self.z, bounds_error=False, fill_value=None)
        xx, yy = np.broadcast_arrays(x, y)
        return f(np.stack([xx.ravel(), yy.ravel()], 1)).reshape(xx.shape)

    def sample(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Heights on the grid x (nx) by y (ny)."""
        X, Y = np.meshgrid(x, y, indexing='ij')
        return self.height(X, Y)


def build_terrain(spec: dict, extent: dict, site_dir: Path) -> Terrain:
    kind = spec.get('type', 'flat')
    if kind == 'flat':
        return Terrain('flat', flat_z=float(spec.get('z', 0.0)))
    if kind == 'relief':
        # Smooth generated relief: filtered noise with a set amplitude and wavelength.
        rng = np.random.default_rng(spec.get('seed', 0))
        step = spec.get('wavelength', 40.0) / 8
        x = np.arange(extent['x'][0], extent['x'][1] + step, step)
        y = np.arange(extent['y'][0], extent['y'][1] + step, step)
        noise = rng.standard_normal((len(x), len(y)))
        kx = np.fft.fftfreq(len(x), step)[:, None]
        ky = np.fft.fftfreq(len(y), step)[None, :]
        k = np.hypot(kx, ky) * spec.get('wavelength', 40.0)
        field = np.real(np.fft.ifft2(np.fft.fft2(noise) * np.exp(-(k ** 2))))
        field *= spec.get('amplitude', 1.0) / (np.abs(field).max() + 1e-12)
        return Terrain('relief', x, y, spec.get('z', 0.0) + field, meta=spec)
    if kind == 'dem':
        d = json.loads((site_dir / spec['file']).read_text())
        x = d['x0'] + d['dx'] * np.arange(d['nx'])
        y = d['y0'] + d['dx'] * np.arange(d['ny'])
        z = np.asarray(d['z'], float) - spec.get('datum_m', 0.0)
        return Terrain('dem', x, y, z, meta={k: v for k, v in d.items() if k != 'z'} | {'datum_m': spec.get('datum_m', 0.0)})
    raise ValueError(f'unknown terrain type {kind!r}')
