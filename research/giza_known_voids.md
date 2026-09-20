# Giza: surveyed-void candidates

Reviewed 20 September 2026. This is a reference inventory and proposed field
comparison, **not a SAR detection result**. Historical excavation measurements
must be checked against the site's condition at the acquisition date.

## Available acquisition

The user's originals are now in `~/tmp/sar/`, an ordinary home-directory folder.
`tmutil isexcluded` reported **Excluded** for that folder and both HDF5 files.
Keep full rasters and large intermediates outside this repository. Reference
PDFs and working previews for this review are in
`/private/tmp/sar-giza-survey-review-20260920/`; they are temporary, reproducible
downloads, not the only copies of the SAR data.

Giza product: `ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5`.
The native float32 I/Q arrays are 113,692 × 11,562 samples. Acquisition start/end
are 2025-08-27 20:26:54.327037 / 20:27:19.089845 UTC. Metadata and small complex
sample reads at the start, centre and end passed preliminary checks. These
checks do not verify every sample or validate a depth estimator.

The supplied quicklook shows the pyramid plateau and surrounding cemeteries.
Exact target-to-pixel registration has not been established. The companion
SICD XML describes 11,558 × 113,620 pixels, whereas the HDF5 array is
113,692 × 11,562: **transposing SICD coordinates is insufficient to map them
onto this HDF5 product**. Use the HDF5/SLC metadata and verify registration on
surface landmarks before defining target pixels. A quicklook overlay is not a
native-image registration solution.

The reproducible [Giza geometry audit](giza_geometry_audit.json) now evaluates
the native RPC metadata without reading SAR pixels. At the centre, adding 10 m
to ellipsoidal height changes the projected range coordinate by −60.85 pixels
(about −26.31 m using a flat-ground conversion). With the scene-average height,
the RPC centre differs from the zero-based `coord_center` index by +3.65 rows
and −30.15 columns. The four corner differences are about +3.61–3.63 rows and
+4.32–5.82 columns. These are metadata consistency checks using an unverified
height assumption, **not measured geolocation errors against surveyed points**.
No offset has been fitted away and no target pixels have been certified.

The [Sacsayhuaman audit](sacsayhuaman_geometry_audit.json) provides the same
diagnostics for that acquisition. Its negative `LINE_SCALE` is retained as
supplied rather than silently changed to a positive number.

Reproduce these small reports without copying or loading the full raster:

```bash
.venv/bin/python scripts/audit_iceye_geometry.py \
  ~/tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5 \
  --sicd-xml /path/to/ICEYE_X33_SICD_SLEDF_951562307_20250827T202654.xml \
  --out /tmp/giza_geometry_audit.json
```

The sidecar argument is optional. Heights supplied to the RPC are WGS84
ellipsoidal heights, not chamber depths or unconverted elevations above sea
level. Source-code hashes and selected-metadata hashes accompany each report;
the large rasters have not been checksummed in full.

## Measured reference candidates

| Candidate | Measured evidence | Reference and locator | Remaining qualification |
|---|---|---|---|
| Osiris Shaft beneath Khafre's causeway | Three connected levels with dimensioned chamber plans and sections. Shaft A is reported as 9.62 m long; these individual shaft lengths are not interchangeable with roof depth below the modern surface. | Hawass (2007), [excavation report](https://giza.fas.harvard.edu/pubdocs/476/full/), pp. 381–388, especially Fig. 2 on p. 382 and Figs. 3, 4, 6. | Establish the survey datum and ground registration. The report describes water on the lowest level; determine water/fill conditions at the 2025 acquisition. |
| Hetepheres I, G 7000 X | Pit depth 27.42 m; chamber extends 5.22 m south of the shaft, is 2.67–2.77 m wide and had a measured height of 1.95 m. | Reisner & Smith (1955), [Giza Necropolis II](https://giza.fas.harvard.edu/pubdocs/131/full/), printed pp. 14–15 (PDF pages 41–42), Figs. 19 and 21. | Pit depth is not a direct roof-depth label. Confirm the modern cover, fill and entrance treatment. The book distinguishes measured geometry from a proposed intended room height; use the measured geometry. |
| Idu burial, G 7102 C; adjacent Qar complex as additional candidates | Idu's shaft is recorded as 8.05 m deep: 5.71 m in rock plus 2.34 m of upper lining. Its chamber dimensions and nearby shafts are separately documented. | Simpson (1976), [Giza Mastabas 2](https://giza.fas.harvard.edu/pubdocs/141/full/), printed p. 29 (PDF page 39), Fig. 13; Qar plans and sections are in Figs. 2–8. | Separate the underground burial chamber from open courts and exposed architecture. Confirm current cavity state and the depth datum. Surface architecture is an important competing explanation. |

These are archaeological surveys, not verified contemporary georeferenced
point clouds. Their usefulness does not depend on constructing decorative 3D
models. A measured plan and section can supply an evaluation geometry once its
coordinate transform, uncertainty and present-day applicability are established.

The [general plan of the Eastern Cemetery](https://giza.fas.harvard.edu/mapsandplans/5613/full/)
(HUMFA_EG002032; *Giza Necropolis I*, map 3) shows G 7000 X, G 7101 and G 7102
relative to the queens' pyramids and adjoining mastabas. It includes a 50 m
scale and points north to the right. The plan was visually inspected on
20 September 2026. It is a useful registration source, but its scanned drawing
does not by itself supply a modern WGS84 transform or per-landmark elevations.

The next deliverable is a survey-to-native-image overlay with independent
landmark residuals and an uncertainty envelope. Do not freeze chamber-sized
target masks until that uncertainty is small enough for the proposed comparison.

## Proposed comparison

1. Establish native SAR geometry and registration from acquisition metadata and
   independently identifiable surface landmarks. Preserve sufficient aperture
   support; do not equate the duration of a spatial crop with observation time.
2. Define survey coordinates, roof/floor elevations, surface elevation, datum,
   uncertainty and cavity condition. Keep uncertain or filled regions labelled
   as such instead of calling them air-filled positives.
3. Check motion sensitivity and registration with controlled inputs. Freeze
   processing, physical assumptions, regions, acceptance criteria and depth
   resolution before inspecting the resulting underground maps.
4. Save predictions before overlaying surveyed geometry. Evaluate location,
   depth, extent and misses separately; retain rejected outputs. Test agreement
   against displaced reference maps and an explanation based on visible surface
   structures. Do not tune depth scales to make a selected shaft match.
5. Reserve other tombs or acquisitions from tuning. These sources have already
   been inspected during design, so the initial comparison is exploratory, not
   a genuinely blinded trial. Independent curation can provide stronger held-out
   evaluation later.
6. Use investigated comparison areas for false-alarm claims. An unmarked area
   on an incomplete excavation map is not established cavity-free ground.

Success at shallow mapped chambers would be evidence within those conditions;
it would not validate much deeper columns. Failure would constrain the tested
method and acquisition, not establish the absence of known cavities or rule out
every possible radar-based inference.

## Reference download provenance

PDFs were retrieved through Harvard Digital Giza's download links and read
locally; only this inventory is stored in the repository.

| Local reference name | SHA-256 |
|---|---|
| `osiris_shaft.pdf` | `6432a210d78050a9effec1b6032d059e4ba0fbc1fb31c1db524d1bbab33a3bf6` |
| `giza_necropolis_2.pdf` | `7bc04870ab103c82f60d5194aceba11cb685a997cc2ec389659ca2a5204fc44b` |
| `giza_mastabas_2.pdf` | `d966a6206bc1c3db3c1a9dc44c4d17849311fc3773c6326c21a38b55d717de29` |
