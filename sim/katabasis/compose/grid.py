"""A regular grid in site coordinates (x east, y north, z up; metres)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Grid:
    origin: tuple[float, float, float]   # centre of cell (0, 0, 0)
    spacing: float
    shape: tuple[int, int, int]          # (nx, ny, nz); z index grows downward from the top

    @classmethod
    def covering(cls, x: tuple[float, float], y: tuple[float, float], z: tuple[float, float],
                 spacing: float) -> 'Grid':
        """Cells centred inside the box; z[0] is the top (highest) layer."""
        nx = int(round((x[1] - x[0]) / spacing))
        ny = int(round((y[1] - y[0]) / spacing))
        nz = int(round((z[1] - z[0]) / spacing))
        return cls((x[0] + spacing / 2, y[0] + spacing / 2, z[1] - spacing / 2), spacing, (nx, ny, nz))

    @property
    def x(self) -> np.ndarray:
        return self.origin[0] + self.spacing * np.arange(self.shape[0])

    @property
    def y(self) -> np.ndarray:
        return self.origin[1] + self.spacing * np.arange(self.shape[1])

    @property
    def z(self) -> np.ndarray:
        """Cell-centre heights, top first (descending)."""
        return self.origin[2] - self.spacing * np.arange(self.shape[2])

    @property
    def extent(self) -> dict:
        h = self.spacing / 2
        return {'x': [float(self.x[0] - h), float(self.x[-1] + h)],
                'y': [float(self.y[0] - h), float(self.y[-1] + h)],
                'z': [float(self.z[-1] - h), float(self.z[0] + h)]}

    def points(self) -> np.ndarray:
        """All cell centres as (nx*ny*nz, 3), C order over (ix, iy, iz)."""
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    def index_box(self, lo: np.ndarray, hi: np.ndarray) -> tuple[slice, slice, slice]:
        """Index slices of cells whose centres may lie in the axis-aligned box [lo, hi]."""
        sx = self._axis(self.x, lo[0], hi[0])
        sy = self._axis(self.y, lo[1], hi[1])
        z = self.z
        iz0 = int(np.searchsorted(-z, -hi[2], side='left'))
        iz1 = int(np.searchsorted(-z, -lo[2], side='right'))
        return sx, sy, slice(max(iz0, 0), min(iz1, len(z)))

    @staticmethod
    def _axis(a: np.ndarray, lo: float, hi: float) -> slice:
        i0 = int(np.searchsorted(a, lo, side='left'))
        i1 = int(np.searchsorted(a, hi, side='right'))
        return slice(max(i0, 0), min(i1, len(a)))
