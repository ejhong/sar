"""Site frames: WGS84 latitude/longitude to local east/north metres and back."""
from __future__ import annotations

import numpy as np

from sarsim.orbit import lla_to_ecef


def enu_matrix(lat0: float, lon0: float) -> np.ndarray:
    la, lo = np.radians(lat0), np.radians(lon0)
    return np.array([[-np.sin(lo), np.cos(lo), 0.0],
                     [-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)],
                     [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]])


def to_local(lat, lon, lat0: float, lon0: float, h=0.0, h0: float = 0.0) -> np.ndarray:
    """East, north, up (m) of points relative to the origin (lat0, lon0, h0)."""
    p = lla_to_ecef(np.asarray(lat, float), np.asarray(lon, float), np.asarray(h, float) + 0 * np.asarray(lat, float))
    o = lla_to_ecef(lat0, lon0, h0)
    return (p - o) @ enu_matrix(lat0, lon0).T
