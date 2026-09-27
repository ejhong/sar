"""sarsim: a small, self-contained toolkit for testing single-image "SAR Doppler tomography".

Modules
-------
geometry  radar geometry and the Doppler <-> slow-time <-> baseline <-> aspect mappings
scene     scatterer scenes (stepped pyramid, desert rocks, vibrating patches) with occlusion
synth     spectral-domain synthesis of a focused single-look complex (SLC) image
subap     Doppler sub-aperture filter bank (reference / offset pairs)
track     batched complex phase cross-correlation (patch pixel tracking, sub-pixel)
tomo      trajectory -> "depth" focusing (paper variant and replication variant)
pipeline  end-to-end reconstructed pipeline
viz       figure styling helpers
"""
from .geometry import Geometry
from .scene import Scatterers, stepped_pyramid, desert_rocks, vibrating_patch, point_targets, concat
from .synth import synthesize
from .subap import SubapBank
from .track import patch_shifts
from .tomo import focus_paper, focus_windows, kz_for_bank
from .pipeline import run_pipeline, grid_targets, line_targets
