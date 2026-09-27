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
  scatter); the surface survey recovers nothing (+0.0%): its first arrivals run along the top cells and never reach
  12 m. (Corrected: until the ray fix the surface inversion did not move at all.)
- **Echo imaging** (P1-03): the chamber's scattered wave, 21 dB below the records and 6 dB below their noise, migrated
  back through the exactly known rock images the chamber from the surface (6× the image elsewhere, brightest 2.8 m from
  its centre, on its top face). The hammer blow and the echo alone play in the viewer (Show: Waves).
- **Ambient noise** (P1-04): 48 noise sources on the ground, the same geophones. Minus the lag derivative of the noise
  correlation retrieves the wave between two geophones (0.90 against a real source). Used as virtual sources for echo
  imaging with the rock known exactly and an unlimited record, the image is only 2.2× brighter around the chamber than
  elsewhere (the hammer: 6×), its brightest point 7 m off; no record length reaches 3×. Fluctuations of a finite record
  are imaged from one realisation and scaled as 1/√T.
- **Full-waveform inversion** (`sim/katabasis/seismic/fwi.py`): adjoint-state gradients of P and S speed from stress
  correlations, L-BFGS with an illumination preconditioner; the gradient is checked against finite differences.
- **Giza from Petrie (1883)**: Khufu's chambers and passages placed from his summary of interior positions (sec. 64);
  Khafre's passage axis and Belzoni's Chamber's east-west span and size from ch. 9; the rest still approximate.

## Next, in order

1. **Full-waveform inversion** (P1-05, running): the surface survey's noisy records fitted whole from a practitioner's
   start (Vp from the surface travel times, Vs = Vp/√3, 3% fast), 60 Hz on a 1 m grid then 120 Hz on 0.5 m. The
   full-band stage commits the inverse crime (same grid and solver as the data); say so wherever it is shown. Then the
   ambient-noise dispersion inversion practitioners actually run (a 1-D shear profile), for completeness.
2. **Make the near-air stencils symmetric**: surface-to-depth reciprocity is off by 1–2% with them (exact without), which
   limits adjoint gradients near the free surface.
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
  the surface wave. Surface refraction picks use 16 stacked blows (50 dB) and a low trigger. Crosshole quality control
  keeps a floor of a quarter period (so the chamber's real delay is not thrown away); the surface survey's must not,
  because at short offsets the surface wave arrives within a quarter period of the P.
- Rays between stations on the ground must stay inside the grid and be credited to the cells they pass through (the
  top row), the same cells the predicted times come from. Until that was fixed the surface inversion never moved, and
  its "nothing" was an accident; an intermediate version drew a false chamber from smoothing alone.
- The solver's updates skip the outer two index layers: keep at least three air cells above the ground (the solver now
  refuses fewer). With two, the free surface is wrong and reciprocity fails for surface stations.
- Hold the absorbing layers' tuning (`pml_vp`) fixed across runs that must compare: otherwise raising a speed anywhere
  retunes every boundary, and finite-difference gradient checks measure the boundary, not the model.
- FWI kernels are densities: the derivative for a cell is the kernel times the cell volume.
- Noise sources spread over the ground retrieve the wave as minus the lag derivative of the correlation, not the
  correlation itself (Snieder 2004).
- From the first investigation: never use the image row span as aperture time; `nyquist_depth` returns π/ΔKz; confine
  sub-aperture banks to the coherent window.
