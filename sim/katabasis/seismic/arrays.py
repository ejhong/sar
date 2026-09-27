"""Acquisition geometries: where the geophones and the sources go.

Positions are in site coordinates. Surface stations are snapped to the first
solid cell under the terrain, as a spiked geophone or a hammer plate is.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..compose.grid import Grid


@dataclass
class Survey:
    name: str
    description: str
    sources: np.ndarray             # (ns, 3)
    receivers: np.ndarray           # (nr, 3)
    source_kind: str = 'force'      # 'force' (hammer, vertical) | 'explosion' (borehole sparker)
    meta: dict = field(default_factory=dict)

    def summary(self) -> dict:
        return {'name': self.name, 'description': self.description, 'sources': len(self.sources),
                'receivers': len(self.receivers), 'source_kind': self.source_kind, **self.meta}


def snap_to_ground(grid: Grid, solid: np.ndarray, xy: np.ndarray) -> np.ndarray:
    """(n, 2) positions -> (n, 3) at the centre of the first solid cell in each column."""
    x0, y0, ztop = grid.origin
    h = grid.spacing
    out = []
    for x, y in xy:
        i = int(round((x - x0) / h))
        j = int(round((y - y0) / h))
        k = int(np.argmax(solid[i, j]))
        out.append([x0 + i * h, y0 + j * h, ztop - k * h])
    return np.array(out, float)


def surface_grid(grid: Grid, solid: np.ndarray, half_width: float, spacing: float,
                 centre: tuple[float, float] = (0.0, 0.0)) -> np.ndarray:
    ax = np.arange(-half_width, half_width + spacing / 2, spacing)
    X, Y = np.meshgrid(ax + centre[0], ax + centre[1], indexing='ij')
    return snap_to_ground(grid, solid, np.stack([X.ravel(), Y.ravel()], 1))


def borehole(grid: Grid, solid: np.ndarray, xy: tuple[float, float], top_depth: float, bottom_depth: float,
             spacing: float) -> np.ndarray:
    """Stations down a vertical borehole, depths measured from the ground at its collar."""
    collar = snap_to_ground(grid, solid, np.array([xy]))[0]
    d = np.arange(top_depth, bottom_depth + spacing / 2, spacing)
    h = grid.spacing
    z = collar[2] - np.round(d / h) * h
    return np.stack([np.full_like(z, collar[0]), np.full_like(z, collar[1]), z], 1)
