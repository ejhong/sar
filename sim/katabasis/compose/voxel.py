"""Composition -> grid: materials, elastic properties and feature masks."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import shapes
from .grid import Grid
from .hetero import von_karman
from .site import Site


@dataclass
class Model:
    grid: Grid
    legend: list[str]                   # material id for each code in `material`
    material: np.ndarray                # uint8 (nx, ny, nz)
    vp: np.ndarray                      # float32 m/s
    vs: np.ndarray                      # float32 m/s
    rho: np.ndarray                     # float32 kg/m3
    feature: np.ndarray                 # int16, index into feature_ids, -1 for none
    feature_ids: list[str]
    surface: np.ndarray                 # terrain height (nx, ny)
    stats: dict = field(default_factory=dict)

    @property
    def air(self) -> np.ndarray:
        return self.material == self.legend.index('air')

    def cells(self, fid: str) -> np.ndarray:
        return self.feature == self.feature_ids.index(fid)


def _column_z(grid: Grid) -> np.ndarray:
    return grid.z[None, None, :]


def voxelise(site: Site, grid: Grid, heterogeneity: bool = True) -> Model:
    mats = site.materials
    legend = site.used_materials()
    code = {m: i for i, m in enumerate(legend)}
    nx, ny, nz = grid.shape
    Z = _column_z(grid)
    surface = site.terrain.sample(grid.x, grid.y)
    S = surface[:, :, None]

    material = np.full(grid.shape, code['air'], np.uint8)
    stratum = np.full(grid.shape, -1, np.int8)
    X, Y = np.meshgrid(grid.x, grid.y, indexing='ij')
    for i, s in enumerate(site.strata):
        top = np.minimum(s.top_height(X, Y, site.terrain), surface)[:, :, None]
        below = np.broadcast_to(Z <= top, grid.shape)
        material[below] = code[s.material]
        stratum[below] = i
    above = np.broadcast_to(Z > S, grid.shape)
    material[above] = code['air']
    stratum[above] = -1

    depth_below_surface = S - Z
    t0 = 0.0
    for c in site.cover:
        t1 = t0 + c['thickness']
        m = np.broadcast_to((depth_below_surface >= t0) & (depth_below_surface < t1), grid.shape)
        material[m] = code[c['material']]
        stratum[m] = -2                                  # cover: no heterogeneity
        t0 = t1

    if site.water_table is not None:
        wet = np.broadcast_to(Z <= site.water_table, grid.shape)
        for mid in list(code):
            sat = mats[mid].saturated_as
            if sat:
                material[wet & (material == code[mid])] = code[sat]

    vp_t = np.array([mats[m].vp.value for m in legend], np.float32)
    vs_t = np.array([mats[m].vs.value for m in legend], np.float32)
    rho_t = np.array([mats[m].rho.value for m in legend], np.float32)
    vp, vs, rho = vp_t[material], vs_t[material], rho_t[material]

    if heterogeneity:
        for i, s in enumerate(site.strata):
            h = s.heterogeneity
            if not h:
                continue
            f = von_karman(grid.shape, grid.spacing, tuple(h['corr_m']), h.get('hurst', 0.5),
                           h['sigma'], h.get('seed', i)).astype(np.float32)
            m = stratum == i
            vp[m] *= 1 + f[m]
            vs[m] *= 1 + f[m]

    for st in site.structures:
        sl = grid.index_box(*shapes.bounds(st['shape']))
        sub = _points(grid, sl)
        inside = shapes.contains(st['shape'], sub).reshape(_shape(sl, grid))
        m = code[st['material']]
        _assign(material, sl, inside, m)
        for arr, table in ((vp, vp_t), (vs, vs_t), (rho, rho_t)):
            _assign(arr, sl, inside, table[m])

    feature = np.full(grid.shape, -1, np.int16)
    fids = [f['id'] for f in site.features]
    counts = {}
    for k, f in enumerate(site.features):
        sl = grid.index_box(*shapes.bounds(f['shape']))
        sub = _points(grid, sl)
        if len(sub) == 0:
            counts[f['id']] = 0
            continue
        inside = shapes.contains(f['shape'], sub).reshape(_shape(sl, grid))
        m = code[f['fill']]
        _assign(material, sl, inside, m)
        for arr, table in ((vp, vp_t), (vs, vs_t), (rho, rho_t)):
            _assign(arr, sl, inside, table[m])
        _assign(feature, sl, inside, k)
        counts[f['id']] = int(inside.sum())

    stats = {'cells': int(material.size), 'air_cells': int((material == code['air']).sum()),
             'feature_cells': counts,
             'material_cells': {m: int((material == code[m]).sum()) for m in legend}}
    return Model(grid, legend, material, vp, vs, rho, feature, fids, surface, stats)


def _shape(sl, grid: Grid) -> tuple[int, int, int]:
    return tuple(len(range(*s.indices(n))) for s, n in zip(sl, grid.shape))


def _points(grid: Grid, sl) -> np.ndarray:
    X, Y, Z = np.meshgrid(grid.x[sl[0]], grid.y[sl[1]], grid.z[sl[2]], indexing='ij')
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)


def _assign(arr: np.ndarray, sl, mask: np.ndarray, value) -> None:
    view = arr[sl]
    view[mask] = value
