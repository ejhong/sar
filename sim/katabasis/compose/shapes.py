"""Solid primitives in site coordinates, with inside tests and bounding boxes.

Every shape is a dict as written in a site file:

    box       centre [x,y,z], size [sx,sy,sz], yaw_deg, pitch_deg
              yaw turns the box about z (counter-clockwise from east); pitch
              then tilts its local y axis upward (a passage rising northward
              has positive pitch, one descending northward negative)
    cylinder  centre [x,y,z], radius, height (vertical axis)
    sphere    centre, radius
    prism     polygon [[x,y],...], bottom, top (vertical extrusion)
    pyramid   centre (of the base) [x,y,z], base (side), height, yaw_deg
"""
from __future__ import annotations

import math

import numpy as np

KINDS = ('box', 'cylinder', 'sphere', 'prism', 'pyramid')


def _rot(yaw_deg: float, pitch_deg: float) -> np.ndarray:
    """Rotation taking local (box) coordinates to site coordinates."""
    cy, sy = math.cos(math.radians(yaw_deg)), math.sin(math.radians(yaw_deg))
    cp, sp = math.cos(math.radians(pitch_deg)), math.sin(math.radians(pitch_deg))
    yaw = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]])
    pitch = np.array([[1, 0, 0], [0, cp, -sp], [0, sp, cp]])   # about local x: +pitch lifts +y
    return yaw @ pitch


def validate(shape: dict) -> None:
    t = shape.get('type')
    if t not in KINDS:
        raise ValueError(f'unknown shape type {t!r}')
    need = {'box': ('centre', 'size'), 'cylinder': ('centre', 'radius', 'height'),
            'sphere': ('centre', 'radius'), 'prism': ('polygon', 'bottom', 'top'),
            'pyramid': ('centre', 'base', 'height')}[t]
    for k in need:
        if k not in shape:
            raise ValueError(f'{t} needs {k!r}')
    if t == 'box' and min(shape['size']) <= 0:
        raise ValueError('box size must be positive')
    if t == 'prism' and (len(shape['polygon']) < 3 or shape['top'] <= shape['bottom']):
        raise ValueError('prism needs 3+ vertices and top > bottom')


def bounds(shape: dict) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned bounding box (lo, hi)."""
    t = shape['type']
    if t == 'box':
        c = np.asarray(shape['centre'], float)
        h = np.asarray(shape['size'], float) / 2
        corners = np.array([[sx, sy, sz] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]) * h
        pts = corners @ _rot(shape.get('yaw_deg', 0), shape.get('pitch_deg', 0)).T + c
        return pts.min(0), pts.max(0)
    if t == 'cylinder':
        c = np.asarray(shape['centre'], float)
        r, h = shape['radius'], shape['height'] / 2
        return c - [r, r, h], c + [r, r, h]
    if t == 'sphere':
        c = np.asarray(shape['centre'], float)
        r = shape['radius']
        return c - r, c + r
    if t == 'prism':
        p = np.asarray(shape['polygon'], float)
        return (np.array([*p.min(0), shape['bottom']]), np.array([*p.max(0), shape['top']]))
    if t == 'pyramid':
        c = np.asarray(shape['centre'], float)
        r = shape['base'] / 2 * math.sqrt(2)
        return c - [r, r, 0], c + [r, r, shape['height']]
    raise ValueError(t)


def contains(shape: dict, pts: np.ndarray) -> np.ndarray:
    """Boolean mask of points (N, 3) inside the shape."""
    t = shape['type']
    pts = np.asarray(pts, float)
    if t == 'box':
        R = _rot(shape.get('yaw_deg', 0), shape.get('pitch_deg', 0))
        local = (pts - np.asarray(shape['centre'], float)) @ R        # R orthonormal: inverse = transpose
        h = np.asarray(shape['size'], float) / 2
        return np.all(np.abs(local) <= h, axis=1)
    if t == 'cylinder':
        d = pts - np.asarray(shape['centre'], float)
        return (d[:, 0] ** 2 + d[:, 1] ** 2 <= shape['radius'] ** 2) & (np.abs(d[:, 2]) <= shape['height'] / 2)
    if t == 'sphere':
        d = pts - np.asarray(shape['centre'], float)
        return np.einsum('ij,ij->i', d, d) <= shape['radius'] ** 2
    if t == 'prism':
        z = pts[:, 2]
        inside = (z >= shape['bottom']) & (z <= shape['top'])
        return inside & _in_polygon(pts[:, 0], pts[:, 1], np.asarray(shape['polygon'], float))
    if t == 'pyramid':
        c = np.asarray(shape['centre'], float)
        yaw = math.radians(shape.get('yaw_deg', 0))
        d = pts - c
        u = d[:, 0] * math.cos(yaw) + d[:, 1] * math.sin(yaw)
        v = -d[:, 0] * math.sin(yaw) + d[:, 1] * math.cos(yaw)
        frac = 1 - d[:, 2] / shape['height']
        half = shape['base'] / 2 * frac
        return (d[:, 2] >= 0) & (frac >= 0) & (np.abs(u) <= half) & (np.abs(v) <= half)
    raise ValueError(t)


def _in_polygon(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Even-odd rule, vectorised over points."""
    inside = np.zeros(x.shape, bool)
    xj, yj = poly[-1]
    for xi, yi in poly:
        crosses = ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi + 1e-300) + xi)
        inside ^= crosses
        xj, yj = xi, yi
    return inside


def volume(shape: dict) -> float:
    t = shape['type']
    if t == 'box':
        return float(np.prod(shape['size']))
    if t == 'cylinder':
        return math.pi * shape['radius'] ** 2 * shape['height']
    if t == 'sphere':
        return 4 / 3 * math.pi * shape['radius'] ** 3
    if t == 'prism':
        p = np.asarray(shape['polygon'], float)
        area = 0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], 1)) - np.dot(p[:, 1], np.roll(p[:, 0], 1)))
        return float(area * (shape['top'] - shape['bottom']))
    if t == 'pyramid':
        return shape['base'] ** 2 * shape['height'] / 3
    raise ValueError(t)


def smallest_dimension(shape: dict) -> float:
    t = shape['type']
    if t == 'box':
        return float(min(shape['size']))
    if t == 'cylinder':
        return float(min(2 * shape['radius'], shape['height']))
    if t == 'sphere':
        return float(2 * shape['radius'])
    if t == 'prism':
        lo, hi = bounds(shape)
        return float(min(hi - lo))
    if t == 'pyramid':
        return float(min(shape['base'], shape['height']))
    raise ValueError(t)
