# Where Underworld is, and what comes next

Written so the work can be picked up cold. The look is `web/DESIGN.md`; what is and is not established is
`METHOD_AUDIT.md`; the first investigation's record is `sim/legacy/ROADMAP.md` (its measured findings stand).

## State (28 September 2026): Phase 2 has answered the question

The desktop session took over from the cloud session and ran Phase 2 on the real Giza geometry. The answer is no,
for three independent reasons, each sufficient on its own (METHOD_AUDIT.md, "Phase 2"):

1. **One dwell is nearly blind to ordinary ground moving (P2-03).** Registering one image's looks against each
   other, as every single-image method does, recovers 1 ± 1% of a block of real ground's motion by complex
   correlation: the looks share their spectrum sample for sample, so the correlation locks at zero.
2. **The chamber's imprint is 10^4 to 10^8 below what a dwell can detect (P2-04, P2-05).** Under Giza's
   microseisms the imprint on the line-of-sight velocity is 3.2e-6 um/s; the smallest a dwell could detect is
   900 um/s on the plateau's texture and 62 um/s on a bright scatterer directly above the chamber.
3. **The published depth is the frequency of the ground's motion (P2-06):** 6.76 m per hertz on this geometry,
   whatever is underground.

Run whole (P2-07), over simulated desert shaking at Giza's level with the chamber's imprint, imaged on the real dwell
and put through the published pipeline unchanged, the depth sections with and without the chamber differ by 1e-8.
The site is now **Underworld**: the front page leads with the answer, `/pictures/` walks through it for anyone,
`/satellite/` holds the evidence.

P2-03 also corrects the first investigation's record: R9's coherence curve is an identity of the looks, R4's floor is
the scatter of an estimator that does not respond to ground motion, and R7 never recovered its injected motion.

Since then (still 28 September):

- **P2-08** reproduces the published pictures' pillars, bands and blocks from the method's arithmetic on the real
  Khafre image (Parseval, the axis's resolution, and the repeat depth where every steering phase coincides).
- **P2-09** puts a wave into the image and reads 3,264 virtual sensors: they are deaf (0% by complex correlation), and
  where magnitudes hear, the chamber would need 2e19 sensors (R11's control had been added after tracking).
- **P2-11** closes every route, not only the published one: penetration, vibration within an image, deformation
  between images under tides, air pressure, heat and the water table, and the surface's own material. X band reaches
  about 30 cm into the driest sand; only P band (Biomass) reaches a few metres, at 25 m resolution. The owner's wish,
  a viewer that looks underground from radar data, is therefore possible only for shallow ground at low frequency, not
  from ICEYE's X band, and not for anything like the claim.
- **The 2026 micro-motion papers** (Vattulainen et al., IEEE Access; Lotti et al., SCHM), the "vibrations from
  satellites work" results, measure isolated bright targets (corner reflectors, a steel bridge) moving at 0.7 mm/s and
  up, with errors near 1 mm/s: what P2-03 says works, 2e8 times the chamber's imprint. Recorded in METHOD_AUDIT.
- **A third product**: the owner supplied ICEYE X13's 2022 Spotlight Dwell over Khufu (`~/tmp/sar/giza2/`, record
  `sites/acquisitions/giza-20220715.json`: 35.0 degrees incidence, 24.56 s, looking east-north-east, so Khufu's east face
  is at 86 degrees local incidence, nearly edge-on and mostly shadowed step by step).
- **P2-13** runs a public gated reconstruction of the method unchanged (its source hashed) on the 2022 image along
  seven loops round Khafre's faces: 5 positions pass on the real image, 13 and 8 on two motionless twins, none over open
  desert, 9 with 2 mm/s vibrations planted (none at a planted frequency). Its pairs' images are 6 s slices sharing
  98.75% of their spectrum, its modes stand for 1.07 to 10.7 Hz where such an image keeps at most 4.7%, a vibration's
  trace is at least 3,000 times as long as it is wide where the gate accepts at most 10, and 31 of 35 passing depths are
  the mode times 6.17 m or its mirror. On `/satellite/#gates`, in general terms.
- **P2-14** runs the implementation's own frozen Khafre geometry (the base corners at 70.5 m above sea level, as sent
  and corrected for the geoid) with the same controls and supports 1 to 5: at its own support of four or five nothing
  passes on the corrected lines in any image; loosened to one, the real image passes 69, the twins 28 and 29, open
  desert 53. Its own steering record gives 6.00 s of every pixel's pass per image. The Great Pyramid's eleven lines use
  a newer format the published code does not read; they lie on the east face, nearly edge-on to the pass.
- **The viewer's satellite data** is exported (`sim/katabasis/export/radar.py`: the track, the draped image, the
  virtual sensors, the method's volumes with and without the chamber over one whole period of its depth axis, and the
  Khafre volume); the viewer's Satellite mode that draws them is next.

- **Start here:** `CLAUDE.md`, this file, `web/DESIGN.md`. `cd sim && uv sync` then `uv run pytest` (≈35 s here);
  `cd web && npm ci && npm run check && npm test && npm run build`. On this Mac, scikit-fmm builds only with
  `CPLUS_INCLUDE_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1`.
- **Not in git** (regenerated by rerunning): `sim/results/cache/`, each run's `*.npz` (P2-04's `kernels.npz`
  included), `sim/data/`. The ICEYE products live only on the desktop in
  `~/tmp/sar/`; their acquisition records are in `sites/acquisitions/`, so everything but P2-01's reading and
  P2-03's real-image part reruns without them (P2-03 then carries its real-image numbers forward, saying so).
- **Machine time** on the desktop (16 cores): a static kernel on the 1 m bench grid 20 s, on a 0.5 m grid 4 min;
  P2-03 about 1 min; P2-05 40 s. In the cloud (4 cores) P1-04 took ~40 min of responses.
- **Stale:** `sim/legacy/results/r07_injected_motion` holds the 23 September run while its script has changed
  since; P1-05's finding says "to the noise level" where the run stopped at its iteration cap.

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
- **Ambient motion, measured** (M1-01, `/ambient/`): MedNet KEG at Kottamya (67 km east of Giza, open data at INGV),
  twelve days of 1996, response removed, McNamara–Buland spectra: median vertical 0.049 µm/s at 0.1–0.3 Hz, 0.029 µm/s at
  1–3 Hz; bracketed by Peterson (1993, Tables 3–4, in `sim/katabasis/ambient/peterson.py`). This replaced the first
  investigation's unsourced "0.1–10 µm/s": the radar's 77 µm/s floor is at least 22× (Peterson high) and about 1,600×
  (measured) above the microseisms. Giza's ground rings near 0.6 Hz and Khufu's pyramid near 2.3 Hz (ELGabry et al.
  2026, Sci. Rep. 16, 14032). No open absolute measurement on the plateau itself has been found.
- **Full-waveform inversion** (P1-05, `sim/katabasis/seismic/fwi.py`): adjoint-state gradients of P and S speed from
  stress correlations, L-BFGS with an illumination preconditioner; gradients checked against finite differences. On
  P1-03's noisy surface records, from uniform rock at the travel-time Vp and Vs = Vp/√3 (3% fast), traces within 4 m
  of the hammer excluded, ten iterations bring the misfit to 0.371 (noise alone: 0.327). The chamber comes back as
  −15.8% in S speed, 5.7× the scatter, slowest cell 1.8 m from its centre; −2.4% in P speed; background Vs 1,881 →
  1,866 m/s (true 1,830). Same solver and grid as the records (the inverse crime). A 60 Hz first stage was tried and
  dropped: this source has too little energy there, and it fitted the near field and noise.
- **Giza from Petrie (1883)**: Khufu's chambers and passages placed from his summary of interior positions (sec. 64);
  Khafre's passage axis and Belzoni's Chamber's east-west span and size from ch. 9; the rest still approximate.

- **Phase 2 on the real geometry** (`sim/sarsim/acquisition.py`, `sim/sarsim/looks.py`,
  `sim/katabasis/seismic/static.py`): P2-01 the dwell's geometry from its products (three speeds, the product's
  slow-time sign, virtual baselines from the state vectors); P2-02 arbitrary motion in the synthesizer; P2-03 what one
  dwell can see of motion; P2-04 the chamber's imprint by static relaxation, checked against Eshelby and Okada;
  P2-05 the budget; P2-06 the depth-frequency identity; P2-07 the whole chain end to end (its figures are rendered by
  the run into `figs/` and exported with it).
- **The real pass on the ground** (P2-12, `#giza/satellite`): the 2025 image resampled onto the Giza terrain, the
  method's volumes from it at Khafre, Khufu, Menkaure and open plateau south of Khafre; a pyramid's picture differs from
  open plateau's by no more than two stretches of open plateau differ from each other.
- **The viewer's Satellite mode** (`web/src/viewer/scene/Radar.ts`, key `s`, `#bench-void/satellite`): the ICEYE
  satellite and its spotlight beam over the site (not to scale), the simulated image draped on the ground, the image's
  virtual sensors each carrying the true motion (ochre) beside what the image reports (cinnabar, complex or magnitudes),
  and the method's volume with and without the chamber in cinnabar. Radar volumes are listed apart in Recovered mode;
  the Khafre claim carries the method's real-image volume. P2-08, P2-09 and P2-11 are on `/satellite/`.

## Next, in order

1. **P2-10, Khufu on two real passes** (run locally, in `sim/private/`, not yet published): lines across the east
   face on the 2022 and 2025 images, with a motionless twin of each. Pre-stated tests: the two passes must agree, and
   the real image must differ from its twin.
2. **A viewer that takes real radar data** (the owner's end goal, bounded by P2-11). Built for the 2025 Giza pass
   (P2-12; the Giza site's Satellite mode): the pass from the product, the image resampled onto the terrain
   (`sarsim.ortho`: RPC over terrain and pyramids, residual offset fitted on predicted brightness), the X-band reach,
   and the published method's volumes at the pyramids and open plateau (the first investigation's patch runner). Any
   product goes through the same steps; next, a command that takes a product path and a site. Then, for low-frequency
   data over dry ground (Biomass P band, NISAR or ALOS-2 L band), a shallow mode that shows buried features within the
   top metres, the one place radar does see underground.
3. **Optional confirmations on the real image**, none of which can change the result: the Khufu known-chamber
   positive control (metre-level registration on surveyed base corners, then the King's and Queen's Chambers and the
   Grand Gallery under the published chain).
4. **Tidy the first investigation's record**: rerun or retire R7 (its script and results disagree); mark R9 and R10
   on the archive as identities of the looks rather than measurements of the ground.
5. Phase 1, if it continues for its own sake: full-waveform inversion beyond the inverse crime; the ambient field in
   every site; symmetric near-air stencils; detectability sweeps; the shaft field and Giza; the claimed underworld
   at low frequency. None bears on the radar question.

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
- Two looks cut from one image share their overlapping spectrum sample for sample: their coherence is an identity of
  the spectrum, and their complex correlation locks at zero whatever the ground does. Test any single-image tracker
  against a region that moves relative to still surroundings, read against a motionless twin; motion common to the
  whole image is undetectable in speckle and proves nothing.
- A look averages motion over its duration (gain sinc(f W)): 1.96 s looks pass little above 0.3 Hz. Space looks
  closer than the Nyquist interval of the fastest tested motion, or the fit aliases.
- Picking the interpolated peak of a narrow look's cross-correlation reads short by about sigma^2 / n pixels;
  `sarsim.looks.region_shift` refines by the cross-spectral phase slope instead.
- Keep a bright test scatterer out of a region meant to test ordinary ground: it dominates the correlation.
- Convergence runs on smaller domains must keep their receivers out of the absorbing layers.
