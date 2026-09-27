"""Copernicus GLO-30 DEM crops for real sites.

The tile is fetched once from the public bucket and cached outside the
repository; only the small local crop is written into the site folder.
"""
from __future__ import annotations

import json
import math
import urllib.request
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata

from .geo import to_local

CACHE = Path.home() / '.cache' / 'katabasis'
BUCKET = 'https://copernicus-dem-30m.s3.amazonaws.com'


def tile_name(lat: float, lon: float) -> str:
    la, lo = math.floor(lat), math.floor(lon)
    ns, ew = ('N' if la >= 0 else 'S'), ('E' if lo >= 0 else 'W')
    return f'Copernicus_DSM_COG_10_{ns}{abs(la):02d}_00_{ew}{abs(lo):03d}_00_DEM'


def fetch_tile(lat: float, lon: float) -> Path:
    name = tile_name(lat, lon)
    path = CACHE / f'{name}.tif'
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f'{BUCKET}/{name}/{name}.tif', path)
    return path


def read_tile(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Heights (rows north->south), and the latitude/longitude of pixel centres."""
    import tifffile
    with tifffile.TiffFile(path) as t:
        page = t.pages[0]
        z = page.asarray().astype(np.float64)
        sx, sy = page.tags['ModelPixelScaleTag'].value[:2]
        tie = page.tags['ModelTiepointTag'].value
    lon = tie[3] + sx * (np.arange(z.shape[1]) + 0.5)
    lat = tie[4] - sy * (np.arange(z.shape[0]) + 0.5)
    return z, lat, lon


def local_crop(lat0: float, lon0: float, x: tuple[float, float], y: tuple[float, float], dx: float,
               inpaint: list[dict] | None = None) -> dict:
    """Resample the DSM onto a local east/north grid (bilinear in lat/lon).

    `inpaint` lists footprints ({'centre': [x, y], 'half': h}) whose DSM
    heights are replaced by an interpolation of the ground around them, so a
    monument drawn from its survey does not sit on its own 30 m blur.
    """
    z, lat, lon = read_tile(fetch_tile(lat0, lon0))
    f = RegularGridInterpolator((lat[::-1], lon), z[::-1], bounds_error=True)
    xs = np.arange(x[0], x[1] + dx / 2, dx)
    ys = np.arange(y[0], y[1] + dx / 2, dx)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    # invert the local frame by Newton steps on the small-area linearisation
    m_lat = to_local(lat0 + 1e-3, lon0, lat0, lon0)[1] / 1e-3
    m_lon = to_local(lat0, lon0 + 1e-3, lat0, lon0)[0] / 1e-3
    LAT, LON = lat0 + Y / m_lat, lon0 + X / m_lon
    for _ in range(2):
        e = to_local(LAT, LON, lat0, lon0)
        LAT += (Y - e[..., 1]) / m_lat
        LON += (X - e[..., 0]) / m_lon
    H = f(np.stack([LAT.ravel(), LON.ravel()], 1)).reshape(X.shape)
    masked = np.zeros(X.shape, bool)
    for fp in inpaint or []:
        cx, cy = fp['centre']
        masked |= (np.abs(X - cx) <= fp['half']) & (np.abs(Y - cy) <= fp['half'])
    if masked.any():
        H[masked] = griddata((X[~masked], Y[~masked]), H[~masked], (X[masked], Y[masked]), method='linear')
    return {'x0': float(xs[0]), 'y0': float(ys[0]), 'dx': float(dx), 'nx': len(xs), 'ny': len(ys),
            'units': 'metres above EGM2008 (orthometric)',
            'source': f'Copernicus GLO-30 DSM, tile {tile_name(lat0, lon0)}, resampled bilinearly from 1 arc-second',
            'inpainted': inpaint or [],
            'z': np.round(H, 2).tolist()}


def write_crop(path: Path, crop: dict) -> None:
    path.write_text(json.dumps(crop, separators=(',', ':')))
