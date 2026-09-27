# Field acquisitions and survey references

The [field atlas](https://ejhong.github.io/sar/field.html) separates actual SAR
intensity, external historical geometry, and measured processing checks.
There are **no field depth reconstructions**. A 3D survey model is a reference
that a later prediction must explain; it is never an output of the radar code.

## Layout

| Location | Purpose |
| --- | --- |
| `catalog/sites.json` | Sites and their individually identified acquisitions |
| `catalog/<site>/design.json` | Explicit native bounds and fixed check settings |
| `catalog/<site>/surveys.json` | Local metric reference geometry, citations, datums and simplifications |
| `catalog/giza/landmarks.json` | Selected published surface points and four original EGM96 grid nodes |
| `fieldwork/processing.py` | Bounded raster reads, intensity reduction, spectrum, timing and known shifts |
| `scripts/run_field_checks.py` | Run checks against an explicitly supplied private acquisition |
| `results/field/<site>/` | Small derived WebP images, report, complete trial CSV |
| `fieldwork/publishing.py` | Verify provenance and export the atlas and reference OBJ meshes |
| `site/field.html`, `site/assets/` | Page templates and canonical frontend source |
| `docs/` | Portable published output; no private raster required for a rebuild |
| `tests/test_fieldwork.py`, `tests/browser_smoke.mjs` | Numerical and desktop/mobile interaction checks |

The existing physical and mechanism studies remain in `experiments/`. Field
checks have identifiers F01–F03 on their own page; they are not additional
simulated depth experiments. Existing T1–T7 retain their scope and history.

## Reproduce the field checks

Keep the full originals in their existing backup-excluded folder. The command
opens them read-only, writes only small derived outputs, and saves no large
intermediate arrays. Full intensity overviews read every I/Q value in chunks;
they do not load a complete raster into memory. Complex analyses use explicit
bounded selections. Overview images never enter phase analysis.

```bash
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python scripts/run_field_checks.py \
  --site giza --input ~/tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python scripts/run_field_checks.py \
  --site sacsayhuaman --input ~/tmp/sar/sacsayhuaman/ICEYE_X35_SLC_SLEDF_5907293_20250822T152433.h5
.venv/bin/python build_site.py
.venv/bin/python -m pytest tests -q
```

Without the private data, the last two commands verify and publish the saved
records. Code, design and derived-artifact hashes reject stale results during
the build. Original sample hashes identify values in declared selections,
separately for the full I and Q overview arrays. They are not vendor checksums
of the HDF5 container. Full-array finiteness does not certify calibration,
geolocation or the absence of other corruption.

F01 projects independently published coordinates **without fitting a correction**.
EGM96 interpolation supplies a provisional `h = H + N` conversion. The GPMP
vertical datum tie is not established, KML altitude zero is not terrain height,
and no observed native correspondences or held-out errors are accepted. The
page's ±5 m adjustment is a sensitivity scenario, not a confidence interval.

F02 distinguishes collection duration, zero-Doppler image span and nominal
bandwidth/rate duration. Its spectra use all rows from three 24-column strips.
No taper, deramp or tuned band recentering is applied. Fourier frequencies of
the focused image are **not vibration frequencies**. Native HDF5/SICD image
plane differences and a dwell-specific frequency-to-time mapping remain open.

F03 injects known Fourier translations into one 512 × 512 complex crop per
acquisition. Twenty-five fixed 32 × 32 patches are evaluated for eight inputs,
with both signs and axes and a zero input. All 200 trials per acquisition are
published. Six groups meet the declared 0.02 px 95th-percentile tolerance; two
larger-shift groups fail in each scene. Finite windows and texture cause bias.
No independent noise was added, and these are dependent patches, not independent
field trials. Recovering an image translation is not evidence of ground motion.
The 0.02 px engineering tolerance is looser than the legacy 0.005 px gate;
inspect the actual errors instead of treating a pass as sensitivity validation.

## Reference geometry

Models use local metres, x east / y north / z up, with z = 0 at the stated
historical shaft mouth. The viewer uses orthographic projection and true
proportions. Plan, section, orbit and chamber focus change the camera only.
The depth slider moves a ruler plane; it does not reveal data or interpolate
an underground radar volume. JSON and OBJ exports retain the source and
simplifications. Separate prism envelopes can overlap at junctions and are
not a watertight engineering reconstruction.

Idu's reported shaft depth is 8.05 m; its model roof depth of 5.89 m is derived
from floor alignment and the measured chamber height. Hetepheres' observed
roof is 25.50 m below the mouth, measured room height is 1.95 m, and reported
pit bottom is 27.42 m. Their 3 cm datum discrepancy is retained. Neither model
represents modern fill, cover or cavity conditions in August 2025. Osiris is
listed as a source-backed candidate without an invented 3D arrangement.

## Add a site or an acquisition

1. Create a site record with a stable ID, descriptive text, survey registry and
   acquisition records in `catalog/sites.json`. An empty `models` list is valid.
2. Give every acquisition a unique ID and its exact filename, acquisition date,
   platform, geometry-audit name, design path and report path. For another scan
   of a site, use e.g. `giza/2026-01-01/design.json` and
   `giza/2026-01-01/report.json`; never overwrite a previous acquisition.
3. Audit native geometry first. The current processor supports native ICEYE
   I/Q products only; other sensors need an explicit adapter. Add the matching
   audit under `research/`, referenced by that acquisition’s `geometry_audit` field.
   The CLI validates dimensions, RPC values, filename and byte count.
4. Declare explicit display crops, full-row spectrum strips and translation
   controls. Bounds, signs, tolerance and selection history belong in the
   design before inspecting the resulting measurements.
5. Add only source-supported geometry. State local axes, datum, measured versus
   derived dimensions, simplifications, and acquisition-date condition. Keep
   `native_target_mask` and `field_detection` null until separately validated.
6. Run the field command, build and tests. Review the page at mobile and desktop
   sizes. The atlas discovers sites and acquisitions from the registry.

The first Giza comparison is exploratory because the surveys informed design.
Subsequent sites can be withheld from tuning. The next scientific deliverable
is accepted surface correspondences with held-out residuals, qualified survey
labels and a validated dwell-product motion adapter. Only then can a fixed
processing method be scored on chamber location, roof depth, extent, misses and
false alarms. Unknown ground must not be labelled a confirmed negative.
