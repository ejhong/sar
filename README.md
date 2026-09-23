# sar — testing single-image "SAR Doppler tomography"

Code and data behind **[Can radar reveal what lies below?](https://ejhong.github.io/sar/)**

The question is whether one radar image can map structure under the ground. The claim under
test is the single-SLC "Doppler tomography" of Biondi & Malanga (2022, retracted 2026) and the
2025 Khafre "underground city" announcement. The repository answers it in two halves.

**Simulation.** Reimplement the published pipeline exactly as described, then feed it scenes
whose contents are known: a motionless stepped pyramid, empty desert, two rocks, a wall of
courses, ground that really vibrates. Whatever it draws at depth in those scenes cannot be
geology. A separate coupled benchmark gives depth recovery a fair chance to succeed, using a
wave-equation forward model and a physical-template inverse, and measures when it does.

**Real acquisitions.** Run the same chain on two ICEYE Spotlight Dwell Fine products, one of
the Giza plateau and one of Sacsayhuamán. A dwell is the most favourable case the method can be
given: its Doppler axis spans about 24.5 s of real aperture, so sub-apertures are separated by
seconds and genuine ground motion is at least in principle within reach. Depth geometry uses
the products' own orbit state vectors rather than a scalar baseline approximation.

## Layout

```
sarsim/           the toolkit
  geometry.py     synthetic radar geometry; Doppler <-> slow time <-> baseline <-> aspect
  scene.py        scatterer scenes: stepped pyramid with layover/occlusion/flash, rocks, vibration
  synth.py        spectral-domain synthesis of a focused SLC (+ speckle clutter)
  subap.py        reference/offset sub-aperture bank in spatial frequency (synthetic side)
  track.py        batched complex cross-correlation with upsampled-DFT sub-pixel refinement
  tomo.py         trajectory -> depth: paper variant and the derivative protocol's branches
  gates.py        the derivative protocol's ellipse acceptance gates
  pipeline.py     end-to-end synthetic run
  waves.py        2-D SH cavity forward model for the coupled benchmark
  measurement.py  reflector acquisition model and the physical-template inverse
  dwell.py        REAL DATA: ICEYE dwell product adapter and Doppler bank in hertz
  orbit.py        REAL DATA: WGS84, state-vector ephemeris, virtual baselines, steering Kz
  geolocation.py  strict RPC projection and geometry audit
  realdata.py     bounded, read-only product inspection
  viz.py          figure styling
experiments/
  t01..t07        simulation chapter: what the method computes and why
  validation.py   coupled physical benchmark and its controls
  feasibility.py  phase readout, motion detection, AI controls
  robustness.py   stress tests of the simulation examples
  real_common.py  REAL DATA: patch orchestration and caching
  r01..r05        REAL DATA: within-image controls, depth axis, split dwell,
                  velocity floor, learned null
fieldwork/        real-image checks: spectra, timing audit, planted translations, previews
catalog/          site definitions and survey references
research/         known-void inventory, geometry audits
scripts/          product inspection and audit entry points
site/, docs/      page template and the built GitHub Pages site
tests/            pytest suites for every module above
```

## Run

```bash
python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest tests

# simulation chapter
cd experiments && ../.venv/bin/python run_all.py

# real chapter (needs the two ICEYE products; paths in experiments/real_common.py)
cd experiments && ../.venv/bin/python r01_giza_controls.py
../.venv/bin/python r02_depth_axis.py
../.venv/bin/python r03_split_dwell.py
../.venv/bin/python r04_velocity_floor.py
../.venv/bin/python r05_learned_null.py

cd .. && .venv/bin/python build_site.py
```

The 10 GB products stay outside the repository. Nothing reads a full raster; every real-data
step works on a bounded, geolocated crop and records the crop origin in its report.

Adding an experiment: write `experiments/<id>.py` that saves figures under `results/<id>/figs`
and a `summary.json` with `id, order, tag, eyebrow, title, question, finding, limitations,
method, figures, metrics`, then rebuild the site.

## Status

Simulation chapter complete. Real chapter: Giza processed at four monuments and two controls,
with the depth-axis, split-dwell, velocity-floor and learned-null tests reported. Open: the
positive control against surveyed chambers, and Sacsayhuamán as a second site.

`METHOD_AUDIT.md` is the running ledger of what is and is not established, including
corrections to earlier versions of this work.
