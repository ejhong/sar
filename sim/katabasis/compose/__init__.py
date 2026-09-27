"""Compositions: what a site is made of, as data, and how it becomes a grid.

A site is a stack of materials under a terrain, with features cut into it:

    terrain       flat, generated relief, or a DEM crop
    cover         thin layers measured down from the terrain (sand sheet, soil)
    strata        layers whose tops are planes, dipping planes or the terrain
    water_table   materials below it take their saturated form
    heterogeneity a von Karman random medium per stratum
    structures    solids above ground (a pyramid of masonry)
    features      voids and fills cut into everything else

Nothing here knows about waves or radar; the voxeliser turns a composition
into whichever property grid a simulation asks for.
"""
from .materials import Material, Property, load_materials
from .site import Site, load_site, list_sites, SITES_DIR
from .grid import Grid
from .voxel import Model, voxelise

__all__ = ['Material', 'Property', 'load_materials', 'Site', 'load_site', 'list_sites', 'SITES_DIR',
           'Grid', 'Model', 'voxelise']
