# Where Katabasis is, and what comes next

Written so the work can be picked up cold. The look is `web/DESIGN.md`; what is and is not established is
`METHOD_AUDIT.md`; the first investigation's record is `sim/legacy/ROADMAP.md` (its measured findings stand).

## The goal

Decide whether one radar image can map structure under the ground, as claimed by Biondi & Malanga (2022, retracted
August 2026) and the 2025 Khafre "underground city" announcement. Either build the best feasible version of that
method, or show precisely why it cannot work, with reproducible evidence rather than argument.

The lab does it by simulation over ground whose contents are known, in two phases:

1. **Seismic tomography, simulated** — the reference ceiling: what the instruments built for the job recover.
2. **The satellite over the same ground** — an ICEYE dwell simulated over the same shaking ground, compared with the
   geophones, then applied to the real products.

## Built

- **Repository** — `sites/` (compositions as data), `sim/` (the Python lab, uv), `web/` (Astro + Three.js, built into
  `docs/`, which Pages serves from `main`). The first investigation is intact in `sim/legacy/` and at `/sar/archive/`.
- **Composition model** (`sim/katabasis/compose`) — a materials library where every property is sourced, derived or
  marked assumed; terrain (flat, generated relief, Copernicus DEM), cover, strata with dipping or surface-following
  tops, water table, von Kármán heterogeneity, structures and features; one voxeliser to any grid.
- **Sites** — `bench-void` (one chamber), `bench-shafts` (a shaft-tomb field), `bench-khafre-claim` (the claim placed as
  stated), `giza` (DEM, three pyramids, surveyed chambers with approximate placement).
- **The underworld viewer** — block diagrams with strata patterns on the faces, a section plane, glowing chambers,
  claimed structures dashed, recovered volumes ray-marched, survey stations, and wavefield playback.
- **3-D elastic solver** (`sim/katabasis/seismic/elastic3d.py`) — validated (P1-01): Stokes point force within 0.3%,
  Rayleigh speed within 0.5% at ten points per wavelength, reciprocity across a void, PML energy to 10⁻⁹.
- **Travel-time tomography** (P1-02) on the one-chamber bench: crosshole recovers the chamber as a −11% shadow (3.5× the
  scatter); the surface survey recovers nothing.
- **Echo imaging** (P1-03): the chamber's scattered wave, 21 dB below the records and 6 dB below their noise, migrated
  back through the exactly known rock images the chamber from the surface (6× the image elsewhere, brightest 2.8 m from
  its centre, on its top face). The hammer blow and the echo alone play in the viewer (Show: Waves).
- **Giza from Petrie (1883)**: Khufu's chambers and passages placed from his summary of interior positions (sec. 64);
  Khafre's passage axis and Belzoni's Chamber's east-west span and size from ch. 9; the rest still approximate.

## Next, in order

1. **Ambient-noise tomography** (P1-04): random surface sources, cross-correlation, Green's-function retrieval checked
   against direct simulation by reciprocity; then imaging with virtual sources, and how record length and noise level
   set what it can see. The seismic cousin of the radar claim.
2. **Full-waveform inversion** beyond the first gradient: iterated, multiscale, starting from a model a practitioner
   would have rather than the exact background.
3. **Detectability sweeps**: chamber size × depth × array aperture and spacing × frequency, per method.
4. **The shaft field and Giza**: the positive control at cemetery scale; Khufu's and Khafre's chambers under realistic
   arrays (Khufu is placed from Petrie; Khafre's lower chamber and passages still need positions).
5. **The claimed underworld**: would ordinary seismology see eight 10 m shafts to 640 m? Needs a coarse grid at low
   frequency.
6. **Phase 2**: see `web/src/pages/satellite.astro` for the plan and the known gaps in `sim/sarsim` (ICEYE geometry,
   arbitrary scatterer motion, decorrelation, per-row zero-Doppler time).

## Open decisions

- **Loose sand cover** — S waves near 300 m/s need a grid about three times finer; the benches are bare rock for now.
- **Densities** are assumed throughout until an open source is found.
- **Deploy** — the site builds into `docs/` and is pushed to `main`; a GitHub Actions deploy is possible once Pages is
  switched to "GitHub Actions".

## Mistakes already made, do not repeat

- Put reciprocity test points outside the absorbing layers; PMLs are not reciprocal.
- Threshold and AIC pickers on a smooth wavelet are biased by the S wave and free-surface ghosts that follow the P;
  time first arrivals by correlation with the leading lobe of the far-field pulse.
- A 6 m void delays a first arrival by about 0.1–0.3 ms: picks must be better than that, which is why crosshole work
  uses high frequencies.
- Keep every station inside the travel-time inversion grid.
- A vertical hammer's direct P along a free surface is some 30 dB below its surface wave: at 30 dB noise a picker times
  the surface wave. Surface refraction picks use 16 stacked blows (50 dB), a low trigger and quality control whose
  threshold never drops below a quarter period (so a real anomaly's delay is not thrown away).
- From the first investigation: never use the image row span as aperture time; `nyquist_depth` returns π/ΔKz; confine
  sub-aperture banks to the coherent window.
