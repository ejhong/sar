# Where this project is, and what comes next

Written so the work can be picked up cold. Style rules are in `STYLE.md`; what is and is not
established is in `METHOD_AUDIT.md`.

## The goal

Decide whether one radar image can map structure under the ground, as claimed by Biondi &
Malanga (2022, retracted August 2026) and the 2025 Khafre "underground city" announcement.
Either build the best feasible version of that method, or show precisely why it cannot work,
with reproducible evidence rather than argument.

## Data in hand

Two ICEYE Spotlight Dwell Fine products, about 10.5 GB each, kept **outside** the repository at
`~/tmp/sar/`:

| | Giza | Sacsayhuamán |
|---|---|---|
| product | ICEYE X33, 2025-08-27 | ICEYE X35, 2025-08-22 |
| shape | 113,692 × 11,562 | 114,606 × 11,482 |
| aperture | 24.5 s | 19.9 s |
| incidence | 20.8° | 20.5° |

Plus a Copernicus 30 m DEM tile for Giza (free, `copernicus-dem-30m.s3.amazonaws.com`).

## The three facts everything rests on

1. **The Doppler axis is slow time.** The image spans 0.7 s of row time but 24.5 s of aperture,
   linked by the Doppler rate. Never use the image row span as aperture time.
2. **The virtual baseline is enormous.** The platform moves 168 km during one dwell, so the
   depth axis of the published steering formula folds back at a few metres and repeats exactly
   at twice that. Published depths need either 1% of the aperture or a chosen wavelength.
3. **Time diversity and coherence are the same axis.** Looks separated by seconds are separated
   by degrees; coherence halves at 0.98 s and is a tenth by 1.76 s. There is no bank design that
   buys slow-motion sensitivity without destroying coherence.

## Settled

Each of these is measured on the real acquisitions, not argued.

- **The budget (R14).** Four independent requirements, each missed: radar reach by 1,277x,
  record length by 81x, motion sensitivity by 8x, lateral resolution by 97x. Repairing any one
  changes none of the others.
- **Positive control (R13).** Hundreds of surveyed burial shafts in the Eastern and Western
  Cemeteries, 5 to 30 m deep, produce no signal at their own depths against bare plateau,
  before or after matching on surface brightness.
- **Within-image controls (R1, R6).** Monuments and open ground score the same, at Giza and at
  Sacsayhuamán, on two satellites and two continents.
- **The obstruction (R9, R10).** Slow time and look angle are the same axis in a dwell.
  Coherence halves at 0.98 s. The brightest thousandth of pixels behave exactly like the
  darkest half, so "a bright stable scatterer would do better" is not supported here.
- **The array idea (R11).** Treating the image as a dense seismic array and cross-correlating
  pixel pairs finds nothing: a peak correlation of 0.003 against 0.45 for an injected wavefield
  put through the identical stack, which the stack recovers on every informative bin.
- **Independent code (R12).** The derivative protocol v1.7 was installed unchanged and run on a
  real Khafre crop. Its registration returns only the integer correlation peak, so it produces
  one distinct vector value across the whole target track, and its shipped depth grid runs to
  100 m where its real ambiguity limit for that geometry is 50 m.
- The depths do not survive a change of look angle within one acquisition (R3).
- A classifier finds nothing in the depth profiles beyond surface brightness (R5).
- Rendered like the published slides, real data reproduces their appearance everywhere (R8).

## Open, in priority order

1. **Single-chamber positive control.** R13 answers the question at cemetery scale, which needs
   placement only to about 50 m. A single named chamber still needs metre-level registration:
   fit a two-parameter offset on monument base corners surveyed by the Glen Dash Foundation,
   validate on held-out corners, then project the King's and Queen's Chambers, the Grand
   Gallery, the Osiris Shaft, Hetepheres G 7000 X and Idu G 7102 C.
2. **Known voids in the viewer.** Once placement is solved, add the surveyed chamber meshes as a
   comparison layer beside the computed voxels. `catalog/giza/surveys.json` already holds
   dimensioned models; they carry no georeferenced position yet.
3. **The honest instrument.** With coherent looks the dwell yields a roughly two-second velocity
   series at about 8 Hz with a floor near 75 µm/s. Characterising that properly, including on
   strong point scatterers, is publishable independently of the pyramid.
4. **Send the derivative-protocol findings to Seyfzadeh.** Its example geometry assumes a 1 km
   aperture span where the real one for that bank is 5 km; its reference and offset passbands
   still overlap 96%; its depth grid runs past its own ambiguity limit.

## Mistakes already made, do not repeat

- The ellipse gate must pick the best harmonic **within** each window before gating, not jointly
  over windows and modes.
- `nyquist_depth` returns π/ΔKz, the range a real two-component fit can use. The complex output
  repeats at **twice** that. Say which convention a number uses.
- A bank swept across the whole aperture produces looks with zero coherence; any "floor"
  measured that way is the search window, not a precision. Confine banks to the R9 window.
- In `sarsim.viz.CMAP_SEQ`, light means low.
- Patch size does not cap the depth range; scatterers far outside the patch alias onto the axis.
