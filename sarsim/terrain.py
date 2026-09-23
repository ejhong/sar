"""Terrain heights from a Copernicus DEM tile, for placing targets in a radar image.

An RPC projection needs an ellipsoidal height. Published monument elevations are orthometric
and a single assumed value costs accuracy: at Giza a ten-metre height error moves the projected
column by tens of metres on the ground, which is larger than the chambers a positive control
would have to resolve. This module reads the free Copernicus 30 m DEM, which is a
GeoTIFF in WGS84 with EGM2008 orthometric heights, and converts to ellipsoidal with a declared
geoid undulation.

Source: https://copernicus-dem-30m.s3.amazonaws.com/ (Copernicus DEM GLO-30, open licence).
The undulation must be supplied; this module will not invent one.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile


class DemTile:
    """One north-up geographic GeoTIFF tile, sampled bilinearly."""

    def __init__(self, path):
        self.path = Path(path)
        with tifffile.TiffFile(self.path) as f:
            page = f.pages[0]
            tags = {t.name: t.value for t in page.tags.values()}
            if 'ModelTiepointTag' not in tags or 'ModelPixelScaleTag' not in tags:
                raise ValueError('tile is not a georeferenced north-up GeoTIFF')
            tie = tags['ModelTiepointTag']
            scale = tags['ModelPixelScaleTag']
            if len(tie) < 6 or len(scale) < 2 or scale[0] <= 0 or scale[1] <= 0:
                raise ValueError('unsupported tiepoint or pixel scale')
            self.lon0, self.lat0 = float(tie[3]), float(tie[4])
            self.dlon, self.dlat = float(scale[0]), float(scale[1])
            self.data = page.asarray().astype(np.float32)
        if self.data.ndim != 2:
            raise ValueError('DEM tile must be a single band')
        n_lat, n_lon = self.data.shape
        self.bounds = (self.lon0, self.lat0 - n_lat * self.dlat,
                       self.lon0 + n_lon * self.dlon, self.lat0)

    def contains(self, latitude, longitude):
        w, s, e, n = self.bounds
        return bool(np.all((longitude >= w) & (longitude <= e) & (latitude >= s) & (latitude <= n)))

    def orthometric(self, latitude, longitude):
        """Bilinear DEM height in metres above the tile's vertical datum."""
        lat = np.atleast_1d(np.asarray(latitude, float))
        lon = np.atleast_1d(np.asarray(longitude, float))
        if not self.contains(lat, lon):
            raise ValueError(f'coordinate outside DEM tile bounds {self.bounds}')
        x = (lon - self.lon0) / self.dlon
        y = (self.lat0 - lat) / self.dlat
        n_lat, n_lon = self.data.shape
        x0 = np.clip(np.floor(x).astype(int), 0, n_lon - 2)
        y0 = np.clip(np.floor(y).astype(int), 0, n_lat - 2)
        fx, fy = x - x0, y - y0
        d = self.data
        top = d[y0, x0] * (1 - fx) + d[y0, x0 + 1] * fx
        bottom = d[y0 + 1, x0] * (1 - fx) + d[y0 + 1, x0 + 1] * fx
        out = top * (1 - fy) + bottom * fy
        if not np.all(np.isfinite(out)):
            raise ValueError('DEM returned nonfinite heights')
        return out if np.ndim(latitude) else float(out[0])

    def ellipsoidal(self, latitude, longitude, geoid_undulation_m):
        """h = H + N, with N supplied explicitly; no default undulation is assumed."""
        if not np.isfinite(geoid_undulation_m):
            raise ValueError('geoid undulation must be a finite number')
        return self.orthometric(latitude, longitude) + float(geoid_undulation_m)
