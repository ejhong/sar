"""A real product laid onto a site: where each patch of ground falls in the image, and the image resampled onto the
ground (orthorectified), so the viewer can drape it on the terrain.

The product's RPC model maps latitude, longitude and ellipsoidal height to (line, sample). A site is a local frame (x
east, y north, metres, orthometric z) around an origin; over a site a kilometre or two across a flat mapping from the
frame to latitude and longitude is good to centimetres. The surface is the site's terrain plus its structures (a digital
surface model), because the radar sees the top of a pyramid, not the rock beneath it.

The RPC is not exact and the geoid is only approximately known, so a residual offset remains. It is fitted once per
product and site by correlating the image with what the surface predicts: each surface cell returns the wave in
proportion to the cosine of its local incidence, and cells that fall into one image cell (layover) add up.
"""
from __future__ import annotations

import numpy as np

from .geolocation import project_rpc

M_PER_DEG_LAT = 111_320.0


def frame_to_lla(origin: dict, x, y):
    """Latitude and longitude (degrees) of site-frame points (x east, y north, metres)."""
    lat0, lon0 = origin['latitude'], origin['longitude']
    lat = lat0 + np.asarray(y, float) / M_PER_DEG_LAT
    lon = lon0 + np.asarray(x, float) / (M_PER_DEG_LAT * np.cos(np.deg2rad(lat0)))
    return lat, lon


def grid_height(grid: dict, x, y):
    """Bilinear height from a grid block {x0, y0, dx, nx, ny, z (x fastest)}, clamped at the edges."""
    z = np.asarray(grid['z'], float).reshape(grid['ny'], grid['nx'])
    fx = np.clip((np.asarray(x, float) - grid['x0']) / grid['dx'], 0, grid['nx'] - 1)
    fy = np.clip((np.asarray(y, float) - grid['y0']) / grid['dx'], 0, grid['ny'] - 1)
    ix = np.minimum(np.floor(fx).astype(int), grid['nx'] - 2)
    iy = np.minimum(np.floor(fy).astype(int), grid['ny'] - 2)
    tx, ty = fx - ix, fy - iy
    return (z[iy, ix] * (1 - tx) * (1 - ty) + z[iy, ix + 1] * tx * (1 - ty)
            + z[iy + 1, ix] * (1 - tx) * ty + z[iy + 1, ix + 1] * tx * ty)


def pyramid_height(shape: dict, x, y):
    """Height of a square pyramid's surface above sea level at (x, y), or -inf outside its base."""
    cx, cy, zb = shape['centre']
    a = np.deg2rad(shape.get('yaw_deg', 0.0))
    u = (np.asarray(x) - cx) * np.cos(a) + (np.asarray(y) - cy) * np.sin(a)
    v = -(np.asarray(x) - cx) * np.sin(a) + (np.asarray(y) - cy) * np.cos(a)
    half = shape['base'] / 2
    r = np.maximum(np.abs(u), np.abs(v))
    return np.where(r < half, zb + shape['height'] * (1 - r / half), -np.inf)


def surface(scene: dict, x, y):
    """The site's surface: terrain, raised by any pyramid standing on it (orthometric metres)."""
    t = scene['terrain']
    z = np.full(np.shape(x), float(t['z'])) if t['kind'] == 'flat' else grid_height(t, x, y)
    for st in scene.get('structures', []):
        if st['shape']['type'] == 'pyramid':
            z = np.maximum(z, pyramid_height(st['shape'], x, y))
    return z


def project(product, scene: dict, x, y, z, geoid_m: float, shift=(0.0, 0.0)):
    """(line, sample) of site points in the product; z orthometric, geoid_m the undulation to ellipsoidal."""
    lat, lon = frame_to_lla(scene['frame']['origin'], x, y)
    rc = project_rpc(product._rpc, np.ravel(lat), np.ravel(lon), np.ravel(np.asarray(z, float) + geoid_m))
    return (rc[:, 0] + shift[0]).reshape(np.shape(x)), (rc[:, 1] + shift[1]).reshape(np.shape(x))


def multilook(product, r0: int, r1: int, c0: int, c1: int, la: int, lr: int, block: int = 4096):
    """Mean intensity over la x lr pixel cells of the product rows [r0, r1) and columns [c0, c1), read in blocks."""
    import h5py
    r0, c0 = max(r0, 0), max(c0, 0)
    r1, c1 = min(r1, product.shape[0]), min(c1, product.shape[1])
    nr, nc = (r1 - r0) // la, (c1 - c0) // lr
    out = np.zeros((nr, nc), np.float32)
    step = (block // la) * la
    with h5py.File(product.path, 'r') as f:
        for a in range(0, nr * la, step):
            b = min(a + step, nr * la)
            i = f['s_i'][r0 + a:r0 + b, c0:c0 + nc * lr].astype(np.float32)
            q = f['s_q'][r0 + a:r0 + b, c0:c0 + nc * lr].astype(np.float32)
            p = i * i + q * q
            out[a // la:b // la] = p.reshape((b - a) // la, la, nc, lr).mean(axis=(1, 3))
    return out, (r0, c0)


def predicted_brightness(scene: dict, x, y, los_enu, cell: float):
    """Relative backscatter of each surface cell: the cosine of its local incidence (zero when facing away)."""
    z = surface(scene, x, y)
    gy, gx = np.gradient(z, cell, cell)
    n = np.stack([-gx, -gy, np.ones_like(z)], axis=-1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return np.clip(n @ np.asarray(los_enu, float), 0.0, None), z


def register(obs: np.ndarray, origin: tuple, la: int, lr: int, rows, cols, brightness,
             max_shift_px: tuple = (1800, 500)):
    """Image-space offset (lines, samples) that best aligns the predicted brightness with the observed intensity.

    The prediction is splatted into the observed cells (layover adds), both are log-scaled, and their cross-
    correlation is searched within max_shift_px. Returns the offset and the normalised correlation at the peak."""
    nr, nc = obs.shape
    pr = np.floor((np.ravel(rows) - origin[0] + 0.5) / la).astype(int)     # pixel k's centre is at coordinate k
    pc = np.floor((np.ravel(cols) - origin[1] + 0.5) / lr).astype(int)
    ok = (pr >= 0) & (pr < nr) & (pc >= 0) & (pc < nc)
    pred = np.zeros((nr, nc))
    np.add.at(pred, (pr[ok], pc[ok]), np.ravel(brightness)[ok])
    a = np.log10(obs + 1e-12)
    a = a - a.mean()
    b = np.log10(pred + 1e-3)
    b = b - b.mean()
    F = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real
    ii = np.fft.fftfreq(nr) * nr
    jj = np.fft.fftfreq(nc) * nc
    win = (np.abs(ii)[:, None] * la <= max_shift_px[0]) & (np.abs(jj)[None, :] * lr <= max_shift_px[1])
    F = np.where(win, F, -np.inf)
    i, j = np.unravel_index(np.argmax(F), F.shape)
    di, dj = ii[i], jj[j]
    corr = float(F[i, j] / (np.linalg.norm(a) * np.linalg.norm(b)))
    return (float(di * la), float(dj * lr)), corr


def sample(obs: np.ndarray, origin: tuple, la: int, lr: int, rows, cols):
    """Bilinear sample of the multilooked image at product (line, sample) positions; NaN outside it."""
    fr = (np.asarray(rows, float) - origin[0] - (la - 1) / 2) / la            # a cell's centre, in pixel coordinates
    fc = (np.asarray(cols, float) - origin[1] - (lr - 1) / 2) / lr
    nr, nc = obs.shape
    inside = (fr >= 0) & (fr <= nr - 1) & (fc >= 0) & (fc <= nc - 1)
    fr, fc = np.clip(fr, 0, nr - 1.001), np.clip(fc, 0, nc - 1.001)
    i, j = np.floor(fr).astype(int), np.floor(fc).astype(int)
    tr, tc = fr - i, fc - j
    v = (obs[i, j] * (1 - tr) * (1 - tc) + obs[i + 1, j] * tr * (1 - tc)
         + obs[i, j + 1] * (1 - tr) * tc + obs[i + 1, j + 1] * tr * tc)
    return np.where(inside, v, np.nan)
