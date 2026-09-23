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

- The method does not distinguish monuments from open desert, at Giza or Sacsayhuamán.
- The depths do not survive a change of look angle within one acquisition.
- A classifier finds nothing in the depth profiles beyond surface brightness.
- Rendered like the published slides, real data reproduces their appearance everywhere.

## Open, in priority order

1. **Positive control against surveyed chambers.** The decisive remaining test. Blocked on
   placement: a 10 m height error moves the projection further than a chamber is wide, and a
   30 m DEM cannot resolve a pyramid. Path: fit a two-parameter offset on monument base corners
   surveyed by the Glen Dash Foundation, validate on held-out corners, then project the King's
   and Queen's Chambers, the Grand Gallery, the Osiris Shaft, Hetepheres G 7000 X and Idu
   G 7102 C and score against chance.
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
