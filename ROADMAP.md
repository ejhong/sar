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
- **Scattered-wave imaging** (P1-03, `sim/experiments/p1_03_scattered.py`, `katabasis/seismic/imaging.py`): written,
  not yet run. Run it next (about 20 minutes), then `katabasis.export` publishes the image and the wavefield movies.

## Next, in order

1. **Ambient-noise tomography** (P1-04): random surface sources, cross-correlation, Green's-function retrieval checked
   against direct simulation by reciprocity; then imaging with virtual sources, and how record length and noise level
   set what it can see. The seismic cousin of the radar claim.
2. **Full-waveform inversion** beyond the first gradient: iterated, multiscale, starting from a model a practitioner
   would have rather than the exact background.
3. **Detectability sweeps**: chamber size × depth × array aperture and spacing × frequency, per method.
4. **The shaft field and Giza**: the positive control at cemetery scale; Khufu's and Khafre's chambers under realistic
   arrays. Re-derive the chamber placements from Petrie (1883) first (the legacy table's Khufu y axis ran south).
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
- On the surface survey the picker times the surface wave (1,570 m/s), not the weak direct P of a vertical hammer:
  trigger on the P window (radial component, lower threshold) before surface refraction results are trusted.
- From the first investigation: never use the image row span as aperture time; `nyquist_depth` returns π/ΔKz; confine
  sub-aperture banks to the coherent window.
