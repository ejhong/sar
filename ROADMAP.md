# Where Underworld is, and what comes next

Written so the work can be picked up cold. The look is `web/DESIGN.md`; what is and is not established is
`METHOD_AUDIT.md`; the first investigation's record is `sim/legacy/ROADMAP.md` (its measured findings stand).

## State (29 September 2026): the methods tested do not establish subsurface imaging; the claims are graded

Revised after an independent review (METHOD_AUDIT.md, "Independent review"): the site no longer says "cannot" or
"every route is closed". It grades each claim on `/satellite/#sure` as established, measured with limits, modelled or
open, lists this lab's own corrections, and says what would settle the question (`/satellite/#settle`). The reasons,
in order of strength:

1. **The published depth is the frequency of the registration series (P2-06):** 6.76 m per hertz on this geometry
   with the paper's 0.48 m (the review notes the paper also implies 0.24 m), with a repeat and a mirror, whatever is
   underground. Algebraic; no validated model turns it into a physical depth.
2. **Correlating one image's looks registers little of ordinary ground's motion (P2-03):** 1 ± 1% of a block of real
   ground's planted motion by complex correlation. Measured for that estimator; the shared-spectrum identity is exact
   for whole images and a measured bias for patches. P2-15 measures the six-second pairs of the gated reconstruction.
3. **In our model, the chamber's imprint is 10^4 to 10^8 below what a dwell can detect (P2-04, P2-05):** a small
   chamber in uniform rock, quasi-static, under the regional microseism level (not measured on the plateau during a
   pass).
4. **A bound on every way of reading one image (P2-25), with the imprint at every frequency (P2-26):** for ordinary
   ground modelled as fully developed speckle, whatever reads one image raises its detection rate above its false-alarm
   rate on the bench chamber by at most 1.3e-9 under Giza's microseisms. Under a truck close by a whole-image reference
   fails, but one line's does not: in the benchmark (P2-36, revision 2) one image lets no method find a shallow room,
   its place or its shape reliably, beside a bouncing lorry; at the image's calibrated noise level (12.7 dB per cell)
   not even an oracle told every reflectivity reaches 95% at 5%. A theorem under stated models; a bound, not what a detector achieves (on small
   test images the best test comes within 1.3 to 1.4 times of it).

Run whole (P2-07), over simulated desert shaking at Giza's level with the chamber's imprint, imaged on the real dwell
and put through the published pipeline unchanged, the depth sections with and without the chamber differ by 1e-8.
The site is now **Underworld**: the front page leads with what we found, `/pictures/` walks through it for anyone,
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
- **P2-15** measures what the gated reconstruction's six-second pairs register of a vibration, through its own masks
  and registration on simulated images: each image moves by 0.06 to 45 times the boxcar value depending on frequency,
  and a pair catches 3-46% of the shift between its images for a bright point, 0.6-2.8% for texture; at 2 mm/s the
  pair shifts change by at most about one 0.01 px step. It replaces the boxcar numbers the site had quoted.

Then (29 September):

- **P2-16, P2-17** hang the gated reconstruction's fit scores in 3-D across the Great Pyramid (95 lines 3 m apart,
  9,595 positions): at a support of one the real image passes 1,355, a motionless copy 660, the same raster over open
  plateau 1,059; at its own support of four, 4, 4 and 0. The scores are no higher inside the surveyed chambers than at
  the same depths elsewhere. The image's texture, not the monument, sets how much passes.
- **P2-18** checks the radar simulator against a pulse-by-pulse one (`sarsim/echo.py`): a Doppler slice is a stretch of
  the pass, shown by a calculation that never assumes it; the fast injection agrees where shifts are resolvable.
- **P2-19** writes synthetic products in the ICEYE layout (`sarsim/slcfile.py`) that the gated reconstruction reads
  unchanged: a corner reflector vibrating at 5 to 20 mm/s is read by a standard tracker and not registered by the
  pipeline, which passes nothing at the reflector.
- **P2-20** runs the gated reconstruction unchanged over the one-chamber bench on synthetic products shaken as Giza
  shakes, inside the band it reads: with the chamber at its real level, boosted ten thousand and a hundred million
  times, and beside a random perturbation of the imprint's size, the changes never gather over the chamber (6, 6 and
  9 of 178, 242 and 300 changed positions lie within 12 m, where chance puts 5 to 9; the null, 4 of 182).
- **P2-21**: across Khafre, run from the 2022 pass (lines east-west and north-south) and the 2025 pass, what the gated
  reconstruction draws does not survive a change of pass beyond chance (1.33 times chance in positions, 27% of depths
  agreeing against 25%), while two layouts on one image agree strongly: each picture belongs to its image.
- **P2-22, P2-23, the decisive benchmark**: one vibrator beside the bench (10.69 Hz, 30 m), solved with and without the
  chamber, read over the same 24.6 s by geophones and by the satellite. Geophones find the chamber with 11 kN; the
  satellite, with 21 corner reflectors and a detector that knows both answers, needs 2.5 x 10^10 N (a floor: at that
  force its tracker keeps 4% of the motion); natural ground offers nothing to track. P2-24 runs the gated reconstruction
  on the same shaking at that force: the chamber changes its picture exactly as random noise of the chamber's size does
  (646 and 648 positions changed, 32 and 32 within 12 m).
- **The site made concise** (the owner's request): one overview with the assessment (the pictures and claim pages folded
  in), six tabs, the hero's headline under the site's name, the lab's hero without mouse capture.
- **The lab's command** runs the gated reconstruction over any area: across Khafre (`lab_khafre`) the real image passes
  1,218 of 9,595 positions at a support of one and none at four, its motionless copy 539 and none, Khufu's pattern.
- **The lab reorganised by instrument** (the owner found "instruments" beside "satellite" confusing): the ground,
  geophones, satellite; methods as cards with their controls one click away; colour by instrument.
- **Sacsayhuamán, a second real site (P2-33 and three lab runs)** (the owner's request): the ICEYE X35 pass of
  22 August 2025 laid on Copernicus terrain with the site's own EGM2008 undulation (46.2 m; the lab command now reads
  each site's), places from OpenStreetMap, and six areas chosen on the laid image before any method ran: the zigzag
  walls, the Rodadero, fields, three stretches of houses. The paper-style pipeline draws the houses alike (0.98-1.00) and
  gives them the strongest columns (8-12% of their pixels among the noisiest 5%); the walls draw apart (0.84-0.89), the
  Rodadero near the houses (0.93-0.97). The areas' motionless copies keep neither the shared profile nor the walls'
  difference, so both come from something the real image holds beyond its brightness pattern (X band reaches about
  30 cm); what, is open. The gated reconstruction passes 463, 2,411 and 1,126 of 9,595 positions over the walls, the
  city grid and the fields at a support of one (its copies 293, 1,743 and 657).
- **The bands drawn, the two levels explained** (the owner's questions): the paper-style volumes are now exported in
  2.5 m depth steps, smoothed only across the ground (the old 4 m depth smoothing and 5 m steps kept a fifth of the 13 m
  bands); P2-34 reads the gated reconstruction's depth scale from each run: 6.17 m a turn and a 308 m repeat on the 2022
  pass, 3.80 m and 190 m on the 2025 passes, where each fit stands twice, near the surface and at its mirror.
- **Rungs and rooms (P2-35)**: as the paper states the method its bands stay beads (0.48) and its repeat block is faint
  (1.2 times the mean); neighbouring pixels sharing their shifts make rungs (up to 0.96), and 20 closely spaced pairs make
  a room (4.4 times). The site's "draws the same shapes" is corrected on the satellite page, the overview and the lab.

Then (29 September, evening; the cloud session, on branch `claude/upbeat-pascal-swk910`, to merge):

- **P2-25, the proof the site lacked** (the owner: "a proof saying the SAR images can't have enough info regardless of
  the method"): `sim/sarsim/information.py` computes the Fisher information one image of speckle holds about a moving
  imprint, exactly and in the continuum; Pinsker turns it into a ceiling on any test. Motion shared by the scene holds
  exactly nothing; slow motion leaves only the texture's stretch along track (R/V_s = 79 s times the velocity's
  gradient). Checked against a brute-force covariance (6e-16) and reached by the locally most powerful detector on model
  speckle and on the lab's synthesizer. The bench chamber on the real dwell: 1.3e-9 over a coin toss under Giza's
  microseisms, 9e-8 at the noisiest stations, 1e-3 under an urban hum all pass, 2e-2 with a truck bouncing 15 m away
  all pass (8e-2 over a room under a 5 m roof); a genie told the reflectivity, with the truck, 0.37, the one case left
  open. Penetration is noted from P2-11, not pursued. On the overview as step 10 (the routes, and the ceiling per kind
  of shaking) and on `/satellite/#bound` (the derivation, its checks, every case).
- **P2-26, the imprint at every frequency** (the resonance escape, and the resonance the method's proponents invoke): the elastic
  solver with and without the chamber under a surface wave from the side and P and S waves from below, 6-120 Hz (240 Hz
  at 0.5 m on a smaller model). No room rings (every response broad, Q under 2); the bench room's imprint rises to 0.16
  of the motion near 80 Hz; a room under a 5 m roof moves the ground over it up to twice the motion near 60 Hz. The
  dynamic imprint meets P2-04's static one at 7 Hz. The air inside a room has 1/19,000 of the rock's impedance.

- **Start here:** `CLAUDE.md`, this file, `web/DESIGN.md`. `cd sim && uv sync` then `uv run pytest` (≈35 s here);
  `cd web && npm ci && npm run check && npm test && npm run build`. On this Mac, scikit-fmm builds only with
  `CPLUS_INCLUDE_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1`.
- **Not in git** (regenerated by rerunning): `sim/results/cache/`, each run's `*.npz` (P2-04's `kernels.npz`
  included; P2-25 recomputes it into the cache when absent; P2-26's `maps.npz` feeds P2-25), `sim/data/`. A full
  `katabasis.export` on a machine without the desktop's local volumes drops them from `web/public/data/sites/`: export
  only the runs there, or restore the sites afterwards. The ICEYE products live only on the desktop in
  `~/tmp/sar/`; their acquisition records are in `sites/acquisitions/`, so everything but P2-01's reading and
  P2-03's real-image part reruns without them (P2-03 then carries its real-image numbers forward, saying so).
- **Machine time** on the desktop (16 cores): a static kernel on the 1 m bench grid 20 s, on a 0.5 m grid 4 min;
  P2-03 about 1 min; P2-05 40 s. In the cloud (4 cores) P1-04 took ~40 min of responses, P2-26 61 min (sixteen solver
  runs), P2-25 23 min (most of it the synthesizer's Monte Carlo).
- **Stale:** `sim/legacy/results/r07_injected_motion` holds the 23 September run while its script has changed
  since; P1-05's finding says "to the noise level" where the run stopped at its iteration cap.

- **The finite proof, integrated** (29 September, night): an independent derivation (Theorems A–E: processing
  contracts TV; the exact complex Gaussian divergence; a certificate straight from a phase-energy envelope with the common
  background's covariance floor; the oracle; two-world depth risk) is ported as `sim/sarsim/finite.py` with its checks and
  applied to the lab's own runs by P2-30, beside P2-25's tight bound and P2-29's oracle. The site's `/proof/` sets out the
  three layers; the open regime (strong shaking close by, known exactly; the thin-roof room with a lorry) is kept in view.
  P2-25 now carries the scattered wave to the scene's farthest corner (the derivation's correction). A container restart
  cost one P2-26 run; P2-26 now keeps each solver run on disk (`results/cache/p2_26/`) and resumes.

- **The organising question, widened** (the owner, 29 September night): does one complex SAR acquisition hold
  reproducible information about actual underground presence, location or geometry that survives controls and predicts
  withheld ground, even if depth stays uncertain? Presence, horizontal location, shape, relative vertical geometry and
  absolute depth are scored apart; an aliased or uncalibrated depth scale no longer dismisses a method by itself. The site
  now reads the alias tests as a calibration ambiguity. P2-31 (the shape test) ran: an isolated room, an L-shaped and a
  branching tunnel and no cavity, generated through the forward model, read by the published method and a reference
  detector at the real level and at labelled amplified diagnostics. At 20 rad the reference detector names every layout
  and places it within 0.6 m; the published blind map names 4 of 9 (p = 0.37, rerun in double precision on three
  independent grounds) and centres no closer than the controls; at the real level every reader is at chance (audit item
  24). P2-32 bounds shape discrimination as detection is.
- **Third review** (30 September): P2-27's smooth-pattern curve is not a universal requirement (a pattern with more
  structure in the same envelope holds far more); the requirement for every pattern is on energy (Theorem C), and
  under a truck close by it rests on how far the scattered wave carries. Large upper bounds are not working regimes;
  shape pairs go through empty ground; the direction guard and the far field by direction are fixed (audit, third
  review). Next, in the review's order: quantify the physical assumptions as declared ranges (excitation, material,
  geometry, coupling, scattering, the far field's attenuation, spatial tails, numerical error), then a shape proof
  over that declared range for two meaningfully different layouts. After it: a
  physically strong positive control (a vibrator or a lorry beside the layouts, not a boosted imprint), multi-level
  layouts for relative vertical geometry, the gated reconstruction on the same products (desktop), and real withheld
  ground (`BLIND_TEST.md`).
- **Fourth review** (30 September): under strong shaking the ground without a room already moves with the source, so
  the speckle certificate and tight bound have no valid reference there (now labelled surrogates); the oracle with
  receiver noise, where common motion cancels, is the valid bound. Tails bounded by the analytic point source, the 39 m
  disk untapered, sampled sides stated, grid allowances applied, outcomes classed (near chance, 95/5 excluded,
  unresolved), the blind test's statistics corrected. **Next milestone:** one fully specified, auditable comparison:
  the same source in both worlds, a valid reference covariance or calibrated noise floor, bounded tails and numerical
  error, a stated scattering and noise model, and bounds for presence, location and shape; parameter ranges after it.
- **Fifth review: the benchmark** (30 September, P2-36, `BENCHMARK.md`): four worlds identical but for the cavity (none,
  a room under a 5 m roof, an L tunnel, the room 6 m east), one force 15 m west (the FTA's truck over a bump, held at
  69 Hz all pass), the whole image, the tail an assumed surface-wave envelope. The certificate applied line by line
  gives each world a valid reference (floor 0.95 beside the lorry; the whole image gives none): one image reaches at
  most 0.03 to 0.06 above chance for presence, a 6 m location or shape (0.14 under every tabulated allowance); only
  the oracle told the texture reaches 0.9, above 33 to 38 dB (30 dB assumed). Checked pulse by pulse on fixed scenes
  and against the exact divergence on lines; the score test implemented on synthesised images reaches its prediction.
  Code-level corrections adopted with it (audit, fifth review). **Next:** the Khafre claim on the real image's powers
  (needs a power map from the desktop), a calibrated noise floor for the real product, the published method on the
  benchmark's strong-case images, then declared parameter ranges around the benchmark.
- **Sixth review: the noise floor traced, the benchmark frozen** (30 September; audit, sixth review). P2-29's −26.7 dB
  NESZ had no source; ICEYE's Table 2-11 gives −18 to −15 dB for Dwell and neither product carries a noise field, so
  the SNR per cell is calibrated from the image (`sarsim.radiometry`: median ground −5.3 dB over −18 dB, 12.7 dB) and
  every oracle row moves with it: the oracle told every reflectivity now excludes 95% at 5% for every case computed.
  The benchmark gained the 10 m room, a correct average for the quiet case, the microseisms' worst-case envelope, a
  Monte Carlo with one noise law and separate calibration grounds, and wording that keeps the envelope beyond 70 m,
  the oracle's crossing and "achieved" in their place; frozen as revision 2. P2-38 computes the claimed deep structure;
  P2-37 reads the benchmark's images with the published method and a told detector. **Next:** the real acquisition, a
  diffuse-plus-bright-scatterer model informed by the image, with Khafre's geometry and a 2 km map of summed power
  (desktop); the claimed shafts; the scattered field measured beyond 70 m.
- **Seventh and eighth reviews: revision 3, the claim computed** (1 October; audit, seventh and eighth reviews). The
  noise floor was misread: ICEYE's documentation 6.0.8 gives Dwell Fine (the 2025 product) −23.7 to −12.6 dB, 6.0.0
  −18 to −15 dB, and neither is a measurement of this acquisition; the SNR is now a set of conditional scenarios with
  the most favourable, 18.4 dB, as the headline, and every oracle row moved with it (beside the lorry 0.13 to 0.42,
  still excluded). Revision 2's other errors (the score test's evaluation grounds inside its calibration range, four
  stale verification rows, the worst realisation's pass mean, no noise fallback in the covariance floor) are repaired
  in revision 3, its worlds unchanged. P2-37 reads the benchmark's images: the told detector names every world at the
  positive control, the published method's blind map is at chance at every level. P2-38 computes the claim, shafts and
  deep void together at 0.2 Hz through a damped half-space (`katabasis.seismic.halfspace`), on both passes, the oracle
  weighted by the image's own brightness over 2 km: any reader of one image at most 2.4 × 10⁻⁵ above chance, the oracle
  6.3 × 10⁻⁴. The one open corner of the strong case is bright persistent ground (the oracle at 40 dB 0.95 to 1.0).
  **Next:** the real acquisition, a diffuse-plus-bright-scatterer model fitted to the image with Khafre's geometry (the
  open corner); the applicable noise from ICEYE; the pyramid in the claim's model; the scattered field measured beyond
  70 m.

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
- **The lab** (`/underworld/`, `web/src/viewer/`): arranged as the lab works. The site, then what to look with: the
  ground as it is, the geophones on it, or the satellite over it (keys t, g, s). Each instrument lists its methods as
  cards built from the exported scene (`web/src/viewer/catalogue.ts`, tested against the exports); the open card picks
  where the method ran, what went in and its control (controls dashed), and says why its picture looks as it does. The
  stage names the picture in one line, the legend shows only what is drawn, colour follows the instrument (verdigris for
  the geophones, cinnabar for the satellite, whichever method), and the badge says what the voxels hold. Old links
  (`recovered`, `waves`) still open. The satellite's pass, the draped image and the image's virtual sensors
  (`scene/Radar.ts`) sit under the satellite; the waves sit under the geophones. On a real site every place a method
  was run shows at once (the owner's plateau "underworld"), made the same way (method, pass, lines, input, support),
  each labelled on the stage and outlined on the ground so that ground outside every outline reads as not yet
  computed; a label or a place chip flies there. The pass and the way the lines are laid are choices wherever they
  were run: things a real feature should survive (P2-21 measures whether it does).
- **The lab's command** (`python -m katabasis.lab gated`, `sim/katabasis/lab.py`): the gated reconstruction, unchanged,
  over any area of a real product (placed through its RPC) or a synthetic one (through the mapping its file records),
  with a motionless twin; each run writes `sim/results/lab_<name>/` and the export adds it to the lab. The method page's
  run table lists every published run and lab run with its record and the command that repeats it.

## Next, in order

1. **The end-to-end benchmark** (the independent review's plan): one elastic simulation shakes the bench, the same
   shaking read by geophones and by the satellite (pulse-level simulator, synthetic products), then every method,
   geophone and satellite, scored on the same truth, with a favourable physical-template detector for the satellite
   beside the gated reconstruction and the paper-style pipeline, and a sweep of chamber size, depth and shaking to map
   where each succeeds. P2-18 (the simulator against pulses), P2-19 (a known vibration) and P2-20 (the known bench) are
   its first pieces.
2. **The lab's own reading**: the best single-image estimator this lab can build (magnitude tracking of bright
   scatterers, the template detector), in the lab beside the published ones, scored the same way.
3. **The paper-style pipeline in the lab's command**, so its pictures too can be made over any area of any product.
4. **The whole plateau** with the gated reconstruction, area by area, each beside its motionless copy.
5. **P2-10, Khufu on two real passes** (run locally, in `sim/private/`, not yet published): lines across the east
   face on the 2022 and 2025 images, with a motionless twin of each. Pre-stated tests: the two passes must agree, and
   the real image must differ from its twin.
6. **Known-chamber tests on the real image**, each with its outcome stated in advance: the Khufu known-chamber
   positive control (metre-level registration on surveyed base corners, then the King's and Queen's Chambers and the
   Grand Gallery under the published chain).
7. **Low-frequency data over dry ground** (Biomass P band, NISAR or ALOS-2 L band): a shallow mode that shows buried
   features within the top metres, the one place radar does see underground.
8. **Tidy the first investigation's record**: rerun or retire R7; mark R9 and R10 on the archive as identities of the
   looks rather than measurements of the ground.

## Long term: the owner's ideas, kept for later (29 September 2026)

Not being built now; recorded so they are designed for, not lost.

1. **An atlas of the underworld.** A 3-D map of the real world, terrain and imagery on top, with what is known of the
   ground beneath: published models, borehole logs, seismic and other surveys, each a layer carrying its source, method,
   date, resolution and uncertainty, added to as surveys are done. Design: every subsurface object is a statement with
   provenance, not a fact; surveys that disagree are shown side by side and scored where a truth exists, never merged
   silently; "no data" stays distinct from "nothing there"; interpolation is visible when switched on. Ingest from the
   usual formats (GeoTIFF, SEG-Y, LAS, GeoJSON, glTF) into each site's frame (x east, y north, z up); large areas stream
   as tiles. The site compositions already hold the seed of this: features with a source and a status (surveyed,
   representative, claimed, truth).
2. **Every subsurface method in the simulation lab.** Ground-penetrating radar, electrical resistivity, microgravity,
   magnetics, muon tomography, InSAR time series, low-frequency radar, alongside the seismic methods and the satellite
   methods already here: each run over the same known ground and scored the same way, the claim's method one among many.
3. **One common layer for every reader:** the evidence each gives for a void at each place (a detector's likelihood),
   drawn in one colour, so that different instruments can finally be compared on the same quantity.
4. **The whole plateau** with the gated reconstruction, area by area, filling the mosaic; and every area from both passes
   and both line layouts, so P2-21's test runs everywhere.
5. **A known-target test at the Osiris Shaft**, once its position and levels are taken from a published survey: the
   frozen recipe beside the shaft, beside ground with no shaft, and on a motionless copy.
6. **The success boundary in full:** P2-22/23 swept over chamber depth and size, source frequency and distance, and
   geophone record length (P2-23 matches the satellite's 24.6 s; longer records help geophones as the square root).
7. **A colour-scale lesson in the lab:** the same volume in a perceptually ordered scale and in a rainbow scale, to show
   how hue bands make structure appear.
8. **The mechanism the method's proponents propose** (an air column resonating in the chamber), once its excitation,
   damping and coupling are specified well enough to compute the surface motion. Not modelled yet: the solver leaves the
   air out (the earlier impedance argument is withdrawn). P2-27 states what any such mechanism would have to deliver:
   the room-specific motion on the ground a detection in one image needs, under stated conditions.
9. **Measured validation** (kept apart, unfinished): ground motion recorded on the plateau during a pass, or shaker and
   corner-reflector data with instruments on them, read alongside the radar; everything so far is simulated.
10. **The frozen, blinded calibration test** (`BLIND_TEST.md`): the synthetic products and sealed key here, the
   published focusing in the cloud, the gated reconstruction on the desktop; the two-depth and hollow-against-solid test
   through the gated reconstruction too (P2-28 runs the published chain only, the gated code being desktop-only).

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
