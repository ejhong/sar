# SAR Depth, Tested

Source for [the open research notebook](https://ejhong.github.io/sar/). The site
separates three questions: can SAR measure surface motion, does that motion
constrain underground depth, and do the Khafre observations establish columns?

The new physical suite connects a buried cavity to a wavefield, complex radar
measurement and depth estimator. It finds conditional recovery under favorable
assumptions, poor recovery at lower SNR, and false alarms under model mismatch.
It also evaluates the reconstructed patch-shift method and surface artifacts.
**Published depth results are still simulations.** The new
[field atlas](https://ejhong.github.io/sar/field.html) contains actual Giza and
Sacsayhuaman radar imagery, dimensioned historical chamber models, full-array
finiteness checks, survey-height sensitivity, full-row image spectra and known
translation controls on real complex texture. All calibration failures are
retained. No field depth reconstruction has been validated. See [METHOD_AUDIT.md](METHOD_AUDIT.md) for the distinction
between original publications, the public derivative protocol and this code.

## Run

```bash
python3.11 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python build_site.py
.venv/bin/python -m pytest tests -q
python3 -m http.server 4175 --directory docs
```

Committed reports and figures rebuild both pages without private inputs or cached
arrays. Fonts and their licenses are served locally. The 3D viewer has no external
JavaScript dependency. Frontend source lives in `site/assets/`; `build_site.py`
copies it into `docs/`. Edit the source, not generated publishing output.

```bash
# Physical benchmarks: a few minutes on the development machine.
.venv/bin/python experiments/validation.py
.venv/bin/python build_site.py

# Corrected legacy gate studies: longer full scene processing.
.venv/bin/python experiments/t06_ordering.py
.venv/bin/python experiments/t07_surface_controls.py

# All active studies, supporting models and site build.
.venv/bin/python experiments/run_all.py
```

Set `MPLCONFIGDIR` to a writable location if needed. Computational results use
fixed seeds; timing depends on hardware. The source hashes in the physical report
identify the implementation. Large wave and scene caches remain ignored.

## New physical studies

`experiments/validation.py` writes `results/validation/design.json` before running.
It commits aggregate metrics, per-trial CSVs and publication-sized PNGs:

- **Motion:** known 1/2/4 Hz displacement, separate null calibration, phase versus
  intensity, and an explicit focused-SLC round trip.
- **Depth:** independent wave physics relative to the steering formula; a finer
  forward mesh and depths excluded from the physical-template library. The
  inverse shares a wave-equation family with the forward model; this is not
  field validation or a second independently implemented wave solver.
- **SNR and mismatch:** detection, localization and false alarms across supplied
  signal strengths, wave speeds, source positions and cavity radii. The SNR
  follow-up is disclosed as exploratory; null thresholds remain separate.
- **Depth versus frequency:** independently varied physical burial depths and
  excitation frequencies feed the reconstructed patch tracker at a fixed
  wavelength. Rejected diagnostic peaks are kept separate from detections.
- **Source ambiguity:** cavity-free alternatives fit free complex source
  strengths at seven receivers, then face fourteen held-out receivers.
- **Independent nulls:** 32 stationary image realizations, independent flank
  samples and per-window full gating, preserving the correlated filter bank.

The model is deliberately limited: 2D homogeneous SH waves, a horizontal
cylindrical void, known active source, seven resolved coherent reflectors,
independent receiver noise and ideal focusing. Physical displacement amplitudes
and SNR are supplied assumptions, not measurements at Giza. These trials do not
simulate vertical columns, spiral ramps, heterogeneous geology or ambient
excitation. Noise repeats are not new geological sites.

## Processing audit

`sarsim/gates.py` selects the best harmonic **within each window**, then applies
shape and amplitude criteria. This fixes the earlier global-window selection
error. `focus_branch_b` admits only fully accepted windows to the real
least-squares depth fit. It returns missing output if none passes.

New independent scene controls use two-arc geometric-median correction; legacy
T1–T7 retain their declared global median. Earlier permutation plots are
historical diagnostics, because arbitrary shuffles change overlapping-look
correlation. The physical SH model tests one polarization; gate rejection of
those signals is not a universal result about all elastic-wave polarizations.

## Large real SAR files

The current originals are stored at `~/tmp/sar/`, outside the repository and
verified as excluded from Time Machine. This home-directory folder is distinct
from the system's temporary `/tmp`. Measured Giza shaft/chamber references and
the proposed comparison are recorded in
[Giza surveyed-void candidates](research/giza_known_voids.md).

Native geolocation now has a separate, metadata-only
[audit command](scripts/audit_iceye_geometry.py) and reproducible reports for
[Giza](research/giza_geometry_audit.json) and
[Sacsayhuaman](research/sacsayhuaman_geometry_audit.json). These expose terrain
height sensitivity and metadata alignment differences; they do not certify
target placement or compute depth. The [field workflow](fieldwork/README.md) now records actual image checks and
unfitted survey projections. The next milestone is accepted independent surface
correspondences with held-out error, a qualified survey datum and a validated
aperture-time adapter, followed by a fixed known-void comparison.

Keep files where they are. Metadata inspection does not load the raster:

```bash
.venv/bin/python scripts/inspect_sar.py /path/to/scene.h5
.venv/bin/python scripts/inspect_sar.py /path/to/scene.h5 \
  --rows 1000 1512 --cols 2000 2512 --out data/inspection
```

Crops require explicit bounds and have a working-memory limit. Outputs default
to ignored `data/inspection/`. HDF5 I/Q and uncompressed complex TIFF are
supported. SICD/NITF uses optional `sarpy`; compressed TIFF needs a windowed
reader. The inspector reports missing information and supplies no simulated
geometry. A crop's spectrum is a QC preview, not a validated aperture-time or
depth axis. Appropriate support is chosen after acquisition inspection.

The machine-readable [field plan](site/real_data_plan.json) requires native
coordinates, metadata, source/model calibration, known-positive and negative
controls, and independent acquisitions/surveys. Unknown ground is not a negative
label. The older `sarsim/iceye.py` is experimental groundwork and should not be
used as a validated real-acquisition geometry adapter.

## Repository

- `sarsim/`: geometry, scenes, image synthesis, registration, selection, wave
  physics, controlled measurement and bounded real-product inspection.
- `experiments/t01_*`–`t07_surface_controls.py`: original artifact, sensitivity
  and paired surface-control studies; retained as supporting records.
- `experiments/feasibility.py`: earlier ZIP-inspired, uncoupled wave/receiver
  demonstrations and artificial AI controls, preserved for provenance.
- `experiments/robustness.py`: repeated surface-spacing examples and algebraic
  statistic checks, distinct from physical depth validation.
- `experiments/t07_wells.py`: archived staged-column illustration. Removed from
  the page's visual argument and the active all-studies runner.
- `catalog/`: site/acquisition registry, explicit field-check designs and source-linked
  survey geometry; structured for additional sites and acquisitions.
- `fieldwork/`: bounded field measurements, exportable charts, model validation and
  publishing. See its README for reproducibility and current scientific limits.
- `results/field/`: small real-image products, provenance and every translation trial.
- `site/page.html`, `site/field.html`, `site/assets/`, `site_content.py`,
  `build_site.py`: canonical templates/assets, data-derived conclusions and renderer.
- `docs/`: generated GitHub Pages site, figures and downloadable data. Private
  local-file inspection output is never automatically published.
- `tests/`: numerical invariants, meaningful gate regressions, sparse huge-file
  reads, report integrity, and optional desktop/mobile browser checks.

Run `node tests/browser_smoke.mjs` with the local page served and headless Chrome
on debugging port 9231. It checks the data-driven scenario explorer, keyboard
controls, deep links, figure dialogs, downloads and layouts from 320 to 1440 px.
It also checks site switching, image zoom, surface-point overlays, both survey
models, section/depth controls, OBJ exports and all field-calibration results.

Research and simulations by ejhong. Original simulation work was developed with
Claude Fable 5.1; subsequent review, coupled physical benchmarks and page revisions
were developed with Codex. Archived comparison images remain copyrighted by
their authors; see `results/t07_wells/published/CREDITS.md`.
