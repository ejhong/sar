# Doppler Tomography, Tested

Source for [the research notebook](https://ejhong.github.io/sar/): a reconstruction of single-image SAR Doppler tomography, tested on known synthetic scenes. The experiments show how stationary surface scattering can produce depth-like patterns. They do not establish the cause of a specific published feature, reproduce the original authors' data, or demonstrate field detection of cavities.

## Run

```bash
python3.11 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest tests
.venv/bin/python build_site.py
python3 -m http.server 4175 --directory docs
```

Committed reports and figures are sufficient to rebuild the page. The page's fonts are served locally; their licenses are in `docs/fonts/`.

Run the small supporting studies independently:

```bash
.venv/bin/python experiments/feasibility.py
.venv/bin/python experiments/robustness.py
.venv/bin/python build_site.py
```

Run the paired T7 surface-target comparison (13 scenes; several minutes):

```bash
.venv/bin/python experiments/t07_surface_controls.py
.venv/bin/python build_site.py
```

Run every scene experiment, supporting control, overview figure, and the site build:

```bash
.venv/bin/python experiments/run_all.py
```

The two full scene analyses take several minutes each. Other controls are smaller; runtime depends on hardware. Set `MPLCONFIGDIR` to a writable directory if matplotlib cannot write its usual cache.

## Contents

- `sarsim/`: radar geometry, scenes, focused-image synthesis, sub-apertures, patch registration, and two focusing formulations.
- `experiments/t01_*` through `t06_*`: stationary pyramid, desert, spacing mechanism, parameter sensitivity, imposed vibration, and selection gates.
- `experiments/t07_surface_controls.py`: one fixed pyramid with no added points, the original eight-point layout, and three seeded random perimeter layouts. Three prescribed strength settings give 13 scenes. Paired measurements use unsmoothed power at identical pixels, a fixed depth scale, shared display references, and the T6 acceptance gate on every grid trajectory. Per-point design and measurements, scene summaries, and all settings are committed and downloadable.
- `experiments/t07_wells.py`: the former T7, now an explicitly illustrative appendix. Hand-placed points, depth rescaling, smoothing and chosen isosurfaces produce renderings; these do not explain the Khafre observations. The legacy filename and page anchor remain usable. The all-experiments runner also regenerates this appendix.
- `experiments/feasibility.py`: ideal pulse receiver, shear-wave cavity model, AI site-confounding controls, and a same-input phase-versus-translation check. Adapted from the supplied `SAR-Voids-Quick-Experiments.zip`; the script is self-contained and does not require private inputs.
- `experiments/robustness.py`: 40 two-scatterer runs, wavelength rescaling, and positive/null controls for the permutation statistic.
- `experiments/overview.py`: opening figure from the scene arrays.
- `experiments/reporting.py`: descriptions that keep the displayed claims consistent with the computed metrics.
- `results/`: portable reports, figures and per-run CSVs. Large arrays remain in ignored `results/cache/`.
- `site/page.html`, `build_site.py`, `docs/style.css`, `docs/site.js`: page template, renderer, presentation and progressive enhancement. `docs/index.html` and `docs/data/` are generated and committed for GitHub Pages.
- `tests/`: numerical regressions and static-site checks.

The page offers expandable experiments, full-size figures, downloadable JSON/CSV data, an algebraic depth-scale explorer, a glossary and a roadmap for independent field validation. Native details and image links work without JavaScript.

## Reproducibility and limits

Scene cache keys include scatterer data, geometry, image content and numerical source files. Zero-energy patches are rejected; constant trajectories have no eligible fit window (score 0, window index -1). Numerical checks include independent steering calibration, an end-to-end two-scatterer example, wavelength invariance, and cache invalidation.

The cavity-wave and radar-receiver models are **not coupled**. The wave amplitude is uncalibrated; the receiver imposes its own amplitudes. The AI controls are artificial datasets, not trained satellite classifiers. No real satellite scene or surveyed void labels were processed in these new controls.

In T6, the shape-only gate passes some static pixels, while the implemented full gate accepts none at the tested scale. T5's response relative to a known static baseline is not a calibrated detection limit. The page reports these distinctions explicitly.

T7 is a paired sensitivity study of this reconstruction, not an independent replication of Khafre. Its 13 scenes share one background; strength variants reuse the same points and phases. The three random layouts and their relative amplitudes are chosen by fixed seeds before inspecting outputs. Added points are ideal isotropic targets without terrain shadowing or an object-specific scattering model. The measured power response is separate from acceptance by the implemented gate. No published column count, depth, or shape is predicted by these controls.

To add an experiment, save `results/<id>/summary.json` with `id`, `order`, `title`, `question`, `finding`, `limitations`, `method`, `figures`, and `metrics`. Figure paths are relative to the repository. Add its short label and leading figures to `EXPERIMENTS` in `build_site.py`. Update report language in `experiments/reporting.py` where needed, then rebuild.

Research and simulations by ejhong. Original simulation work was developed with Claude Fable 5.1; subsequent numerical review, supporting controls and page revisions were developed with Codex. Published comparison images remain copyrighted by their authors; see `results/t07_wells/published/CREDITS.md`.
