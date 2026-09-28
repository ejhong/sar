# What is being tested

This notebook evaluates conditional feasibility and failure mechanisms. It does
not reproduce the original Khafre data or claim that a failed reconstruction
settles the existence of underground columns.

Three methods must remain distinguishable:

1. Biondi and coauthors' published descriptions and original implementation.
2. Public **reconstructions** of the method (a derivative protocol, v1.5, inspected at a fixed revision), which
   acknowledge missing original details and add empirical choices of their own.
3. This repository's specified simulation and processing choices.

## Corrected selection and focusing

The previous ellipse implementation maximized adjusted R² over *both windows
and modes* before gating. A perfect one-dimensional fit could defeat a slightly
noisier accepted ellipse in another window. The regression in `tests/test_gates.py`
constructs that counterexample.

`sarsim/gates.py` now chooses the best harmonic separately in each W25 window,
then tests adjusted R² ≥ 0.25, minor/major axis ≥ 0.1 and minor semiaxis ≥ 0.005 px.
`focus_branch_b` in `sarsim/tomo.py` restricts real two-component least-squares
fitting to those accepted windows. No accepted window yields NaN scores and
window index -1. The ungated branch remains available as a diagnostic/reference.
Legacy complex projection is not presented as the derivative selective branch.

T6 and T7 have been rerun with per-window selection. Gate counts refer to
dependent target samples, not independent field detections.

## Nuisance and geometry

New independent-scene controls use two separate ten-sample flanks at least 48 px
from the target strip. Each side is aggregated by a two-component geometric
median; their mean is subtracted. Side disagreement and target identities are
recorded. The older T1–T7 figures retain their declared global-median correction.
The physical reflector comparison assumes favorable zero nuisance rather than
pretending to implement a terrain-based flank correction on seven isolated points.

Our registration uses complex 32×32 patches and Fourier refinement. It does not
implement every context-window refinement described in derivative documentation.
Our synthetic geometry is straight-line and flat-earth, with an explicit
frequency-to-effective-time mapping. Real-product state vectors, azimuth
resampling, support, autofocus and acquisition timing have not been validated.

## Physical benchmark

The new chain is:

`2D shear-wave cavity model → surface displacement → phase-modulated radar channels → focused SLC → favorable phase readout → physical template inverse`.

The forward wave solver never uses the Doppler-to-depth steering formula.
Forward evaluation uses 2.5 m spacing and depths 70, 110, 150, 190 and 230 m.
The prediction library uses 5 m spacing and a disjoint 20 m depth grid. Both
share an SH wave-equation family: this is not an independent field validation or
an independently implemented elastodynamic solver.

The cavity is an infinite **horizontal** cylindrical hole in a 2D anti-plane
model, with radius 20 m, a known active Ricker source, uniform rock and no intrinsic
attenuation. A shared cavity-free reference sets the force/displacement scale;
targets are never normalized individually. The source is not a model of ambient
excitation. Three-dimensional columns, spiral geometry and kilometre-scale
depth are outside the simulated domain.

The radar is a stationary-phase spectral approximation with seven separated
reflectors, a specified SH-to-line-of-sight projection, independent complex
receiver noise and invertible azimuth focusing. Readout uses known range supports
and removes constant phase. These favorable conditions avoid claiming that a
practical distributed-ground extraction method has been built. Noise trials use
the algebraically equivalent separated-channel readout; explicit focused-image
round trips verify the equality. SNR is per channel/aperture sample, not a measured
image SNR at Giza.

The reference inverse profiles a positive unknown source gain against physical
wave templates. Its candidates and null thresholds are fixed before their
evaluation observations are scored. Evaluation depths are withheld from the
library; source location and timing are known in the favorable case. This is a
held-out automated benchmark, not formal human blinding or a preregistered trial.
The 20/40/60 dB sweep was added after low recovery in the initial 20 dB study;
it uses separate seeds and new null calibration. The entire sweep, including
failures, is reported.

## Controls and uncertainty

Each detection threshold uses independent cavity-free calibration noise.
Evaluation nulls include finer-mesh discrepancy. The notebook reports misses,
localization, and false alarms together, including cases where model mismatch
breaks the nominal 1% threshold. Wilson intervals summarize noise-run variation
under specified models; repeated noise draws are not new geological sites.

Independent stationary image realizations preserve correlations induced by the
sub-aperture bank. Arbitrary look shuffles in legacy T6 do not preserve those
correlations and are retained only as historical diagnostics. Shape-only passes
are not false positives of the full rule.

The source-ambiguity test fits unrestricted complex forcing at one frequency.
Seven weights can match seven receivers, and fourteen held-out receivers test
generalization. This conditional non-uniqueness example does not establish that
such forcing occurs at Khafre or that broadband, constrained inference is impossible.

## Real observations

`scripts/inspect_sar.py` provides metadata-only inspection by default and explicit,
bounded native-image crops. It supplies no default simulation geometry and makes
no depth inference. HDF5 I/Q and memory-mappable complex TIFF are supported;
SICD/NITF uses optional sarpy. Compressed TIFFs require a separate windowed reader.
Published depth results remain synthetic. Real Giza and Sacsayhuaman acquisitions
have passed preliminary metadata and small-sample read checks. The separate
`scripts/audit_iceye_geometry.py` evaluates native RPC geolocation metadata,
reports height sensitivity and metadata-coordinate residuals, and checks the
companion SICD dimensions when supplied. It reads no SAR pixels and establishes
neither independent landmark registration nor depth accuracy. Historical
[Giza survey candidates](research/giza_known_voids.md) have been identified;
qualified field labels and evaluated chamber detections remain pending.

Before field interpretation: verify product metadata and physical aperture
support, freeze regions and processing, include known motion and static controls,
compare independent acquisitions, and evaluate surveyed cavities plus matched
intact ground at sites excluded from calibration. Unknown subsurface labels
cannot serve as confirmed negatives. Repeated features alone do not exclude
stable surface artifacts.

## Real acquisitions — 23 September 2026

The real chapter processes two ICEYE Spotlight Dwell Fine products. What follows records what
that processing does and does not establish.

### What is validated

The dwell adapter reproduces each product's own declared centre-of-aperture state exactly: the
cubic state-vector interpolation evaluated at the midpoint of the collection returns `coa_pos`
and `coa_vel` to within a millimetre. The RPC projection places Khafre, Khufu, Menkaure and the
Sphinx inside the Giza frame with no fitted offset, and the Khafre crop shows the laid-over
near face and the shadow where the geometry requires them. The registration estimator was
already checked on this actual radar texture with planted translations down to 0.005 px, with a
worst-case error of 0.12 px.

The steering geometry now uses the product's real orbit state vectors. Both the 2022 paper and
the derivative protocol v1.7 substitute a scalar aperture-span approximation; v1.7's own README
records that as an outstanding deficiency.

### What is a consistency check, not a validation

Slow time is recovered from Doppler frequency through the product's own Doppler-rate
polynomial, evaluated in the convention used by NGA's sarpy ICEYE reader. The aperture duration
this implies, 24.49 s, agrees with the declared collection duration of 24.76 s to about one
percent. That agreement is a local consistency check on the adapter. A vendor-confirmed
frequency-to-time model for dwell products has not been obtained, and the sign and zero point of
the Doppler centroid are taken from the measured spectrum of each crop rather than from the
`dc_estimate_coeffs` polynomial.

### What is not established

No published figure has been reproduced. The bank parameters, patch positions, nuisance
geometry and wavelength model behind the Khafre images were never disclosed, so the runs here
reproduce the method as described rather than a specific result. Patch centres and heights are
our choices; heights use published elevations plus a nominal geoid undulation, which is good
enough to land on the right structure but not a surveyed tie.

A scene-median common mode is removed from every patch, so any motion common to a whole patch
is removed with it. The velocity floor in R4 is therefore a bound on differential motion, and
it mixes measurement noise with real clutter decorrelation; it is an upper bound on achievable
precision, not a fundamental limit. Natural desert scattering is the worst case; corner
reflectors would do better.

Negative results on our six patches do not establish the absence of anything underground. They
establish that this processing chain does not distinguish those patches from one another.

### The core geometric objection, stated precisely

The published steering wavenumber is `Kz = 4 pi B_perp / (lambda_s r sin theta)`, the standard
multi-baseline tomographic form. In multi-baseline tomography `B_perp` is a cross-track
separation between passes, which is what resolves height. Here it is the platform's own
along-track displacement inside a single aperture, projected perpendicular to the line of
sight. Substituting one for the other is the step that converts an oscillation rate across look
angle into a depth. R2 measures the consequence with real state vectors: the along-track
baseline spans about 168 km, so the depth axis repeats every few metres unless most of the
aperture is discarded.

### Correction, 23 September 2026: the velocity floor

An earlier version of R4 and R7 swept sub-apertures across the whole processed band. R9 then
measured what that costs: coherence with the band-centre look halves at 0.98 s of separation
and is a tenth by 1.76 s, so the extreme looks of a full-band sweep share no coherence at all.
The registration peak between them is a random position inside the patch, and the "velocity
floor" of about 4,400 µm/s those runs reported was the width of that search, not a measurement
precision. Both experiments now use a bank confined to the coherent window, which gives a floor
near 77 µm/s over a series lasting about two seconds. The earlier number should not be cited.

### The modal branch, examined directly, 23 September 2026

The derivative protocol has a second branch, described in its own documents as the selective
modal or ellipse-gated branch, and often referred to informally as a spring model. R15 runs its
published `modal.py` and `focus.py` unmodified, with their sha256 checked before import, on real
sub-pixel registration vectors from the Giza acquisition. Three findings, in order of weight.

**The depth recurrence is a property of the model, not of our implementation and not of the
data.** The depth stage fits each vector sequence to the two columns `cos(Kz_k z)` and
`sin(Kz_k z)`. When the steering wavenumbers lie on a lattice `Kz_k = K0 + k dKz`, substituting
`z + 2 pi / dKz` adds `2 pi k` to each argument, a whole number of turns, plus one constant
offset common to every k. The two columns are therefore the same pair of columns rotated by a
fixed angle, which spans the same plane, so every least-squares residual is unchanged. For this
geometry the ladder is even to about 0.04 per cent and the measured score correlation one repeat
apart is 0.99947. The protocol states the same result in its section 11.4, gives the unsigned
interval as `pi / |dKz|` by way of its signed-depth equivalence in section 11.3, and carries a
ghost rule warning that repeated peaks are not independent structures. None of this was hidden;
the issue is the size of that interval, 13.7 m here, against the depths being reported.

**The modal gate cannot break the ambiguity, by construction.** It runs before the depth stage
and only decides which windows are admitted to it. The recurrence identity above holds for every
input sequence without exception, so a subset of admitted windows inherits it: gated windows
repeat at 0.99946 against 0.99947 ungated.

**The gate is not selective.** At the frozen Branch B settings it admits nothing on this
acquisition, because its absolute minor-axis threshold of 0.005 px is about 25 times the
displacements actually present. With that size threshold removed so the shape rule can act, it
admits 1.9 per cent on Khafre against 1.2 per cent on open plateau, and the modal score
separates the two at AUC 0.504. A random walk with no structure at all passes 73.9 per cent of
the time, because the fit rewards smoothness along the sweep rather than anything underground.

The harmonic index m is indexed by sub-aperture position, not by time. The protocol's own symbol
table says it "is a harmonic index, not automatically a physical cavity eigenmode". A genuine
mechanical resonance is a different proposal and fails for a different reason: at a shear speed
of 2,200 m/s a chamber 30 m down rings near 18 Hz, while a dwell whose coherence halves in about
a second can carry only 0.08 to 0.51 Hz, short by a factor of 36.

### Retired: the old browser smoke test, 23 September 2026

`tests/browser_smoke.mjs` still drove the controls of the retired essay-style site (`#snr`,
`#scenario`, `#predicted-depth`) and could not have passed since the dashboard replaced it. It
has been rewritten against the current pages and now also checks the field atlas. It skips the
volume-rendering assertions when the browser has no WebGL2, so it needs Chrome started with
`--use-angle=swiftshader` to exercise the viewer.

### Correction of scope, 23 September 2026: "you need minutes"

It was objected that the micro-motion literature, including long-standing work from the
University of Zurich, establishes SAR sensitivity in the hundreds of micrometres from a single
pass, and that measuring phase evolution precisely does not take minutes. Both points are
correct, and the site previously stated the time requirement too broadly.

Our own numbers agree with that literature rather than disputing it. R4's velocity floor on this
dwell is 77 µm/s over a two-second series, about 13 times better than the roughly 1,000 µm/s
reported for single-pass micro-motion on bridges. Phase precision is not what fails here, and
the argument never depended on it failing.

The "minutes" requirement applies to one specific route: using the image as a dense seismometer
array and reconstructing structure from wave propagation, which is the strongest form of the
claim and the one R11 tests. Ambient-noise tomography converges only after many wave periods,
and microseism periods are seconds. It does not apply to measuring a single displacement, which
is immediate. Argument step 3 has been rewritten to say which of the two it means.

R16 then shows precision is not the binding constraint at all, by a bound that contains no
precision term. Reaching a feature of size s through the ground needs an elastic wave no longer
than about 2s, hence a frequency of at least Vs/2s, hence sub-apertures short enough to sample
it. Azimuth resolution is inversely proportional to integration time, so

    rho = rho_full * T * f_s  >=  rho_full * T * Vs / s,

and requiring `rho <= s` gives `s_min = sqrt(rho_full * T * Vs)`. For this acquisition, with
rho_full = 0.050 m and T = 24.76 s, that is 52 m, and it never falls below 27 m across azimuth
resolutions from 0.02 to 0.5 m and shear speeds from 1,500 to 2,500 m/s. No term in it responds
to better phase precision, better coherence or better signal-to-noise.

**What R16 does not exclude.** The largest claimed Khafre structure, at about 79 m, is above the
bound, so this test does not rule it out and must not be cited as if it did. It is excluded
instead by the depth axis, which repeats every 27.4 m and therefore cannot place anything at
648 m. Every chamber actually surveyed at Giza is far below the bound: the Subterranean Chamber
by 43 times, the King's Chamber by 101, a mastaba shaft by 1,209.

## Phase 2: the satellite over simulated ground — 28 September 2026

Phase 2 simulates the real Giza dwell over ground whose contents and motion are known
(`sim/experiments/p2_0*.py`, `sim/sarsim/acquisition.py`, `sim/sarsim/looks.py`,
`sim/katabasis/seismic/static.py`). What it establishes, what it corrects in the first
investigation's record, and what it leaves open follow.

### What is validated

- **Geometry (P2-01).** `DwellGeometry` is read from each product's metadata into
  `sites/acquisitions/`. It keeps the three speeds a dwell has (platform 7,662 m/s, ground sweep of
  the rows 7,014 m/s, 7,330 m/s from the Doppler rate); the effective speed is their geometric mean
  to 9 ppm, and the incidence it implies matches the product's to 0.001°. On this geometry the
  synthesizer's point response is 4.48 × 22.1 cm as predicted, and a target receding at 2 mm/s moves
  −0.159 m in azimuth and is recovered as 2.000 mm/s.
- **Static relaxation (P2-04).** Intact rock under a horizontal stress stays exactly still; with the
  chamber open the result changes by 3.5% (epicentre) and 10% (ring) between 1 m and 0.5 m grids and
  not at all between 100 and 180 m domains. Eshelby's void with Okada's point sources, checked exactly
  against Mogi's source, agrees within the 40% expected of a point sphere against a 6 m cube whose roof
  is 12 m down.
- **The detection threshold (P2-05).** The formula's per-look scatter for a bright scatterer (21 µm/s)
  is checked by simulation (16 µm/s), and a sway at the computed threshold is found in 24 of 24
  realisations at a 5% false-alarm rate.

### What is established

1. **Complex registration of one image's looks is blind to the ground's motion (P2-03).** Two looks
   cut from one image share their overlapping spectrum sample for sample, so their cross-spectrum is
   |Y|² H_A H_B, real and non-negative, and their correlation locks at zero lag whatever the ground
   does. A 20 m block of ground swaying at 2 mm/s comes back at 3 ± 1% on simulated desert and 1 ± 1%
   on real Giza texture, though it is in the data (102% and 101% read against a motionless twin, which
   no real dwell has). For fully developed speckle this is exact in distribution: every Doppler bin is
   an independent circular Gaussian, so a phase history common to a patch changes nothing. Magnitude
   registration catches what is point-like in the ground (10 ± 2% simulated, 21 ± 7% real), and a lone
   bright scatterer's envelope follows its motion in full (100 ± 2%).
2. **The chamber's imprint on the ground's motion is 6.6 × 10⁻⁵ of that motion (P2-04).** A Rayleigh
   wave strains the ground by hv V / c whatever its frequency; the imprint on the line-of-sight
   velocity over the one-chamber bench is 3.2 × 10⁻⁶ µm/s under Giza's measured microseisms and
   2.3 × 10⁻⁴ µm/s at the noisiest stations on Earth.
3. **One dwell falls short of that imprint by 10⁴ to 10⁸ (P2-05).** Allowing every look, a 3-sigma
   threshold and the best tracker: 3 × 10⁸ on the plateau's own texture, 2 × 10⁷ on a bright scatterer
   directly above the chamber, 2 × 10⁶ on a corner reflector no archive image contains; 4 × 10⁶,
   3 × 10⁵ and 3 × 10⁴ at the noisiest stations on Earth.
4. **The published depth is the frequency of the ground's motion (P2-06).** Kz grows linearly with
   slow time on a dwell, so the focusing is a Fourier transform over slow time and puts motion at
   z = f λ_s R sin θ / (2 V_s): 6.76 m per hertz on this geometry with λ_s = 0.48 m. Given the exact
   motion at every pair, it draws the same microseism at the same 1.3 m over a chamber 15 m down, one
   30 m down, and open ground (profiles identical to 10⁻¹⁵).

5. **Run whole, the chain cannot tell the chamber is there (P2-07).** Simulated desert over the one-chamber bench,
   shaking at Giza's measured microseism level with the chamber's imprint (P2-04's kernels), imaged on the real dwell in
   double precision and put through the published pipeline unchanged: the chamber changes the image by 1e-10 and the
   depth section by 1e-8. Its imprint must be multiplied by 1e10 (tens of radians of phase) before the section changes
   visibly, and then the change is spread over every depth (0.68 to 1.25 of its mean): the texture is reshuffled, nothing
   is drawn at the chamber. Over the chamber the dwell's velocity series (R4's reading) wanders 10,000 times more than a
   seismometer there records. With the published lambda_s the axis ends at 13.4 m, above the chamber; lambda_s only
   relabels the same pixels (P2-06).

6. **The published pictures' pillars, bands and blocks are the method's own arithmetic (P2-08).** Averaged over a
   period of the depth axis, the power drawn at a pixel equals that pixel's trajectory energy (Parseval; correlation
   1.000 on the real Khafre image), so a pixel whose registration wanders is bright at every depth and stands as a
   column in any thick-slab rendering; such pixels lie where the surface is darker or changes brightness, and open
   plateau with no monument draws the tallest columns. Along a column the power rises and falls at the axis's
   resolution (13 m relabelled): the bands. At z = 0 and at every repeat depth (648 m and 1,296 m relabelled) all
   steering phases coincide and the power is the pixel's mean offset squared (ratio 0.994): the blocks. Below 648 m
   the picture repeats.
7. **A radar image is not a dense seismic array (P2-09).** With the motion put into the image itself (a 60 m, 0.15 Hz
   test wave), 3,264 virtual sensors read by complex correlation recover 0% of the motion under them and by magnitudes
   −6%. Still sensors do average down as n^−0.50, but a deaf sensor stacks to a precise zero. Where magnitudes hear
   (a fifth, over real point-like ground), the chamber would need about 2 × 10¹⁹ independent sensors over it; there
   are 3,602 there, 7 × 10⁷ in the whole image and 4 × 10¹⁴ on all of Earth's land. The ground's own trembling would
   need 1,186 times more than the image holds.
8. **No route carries a Giza room to an X-band radar, and none carries the claimed structure to any radar (P2-11).**
   An echo depends only on how strongly each scatterer returns the wave and where it is while watched, so a room can
   reach the record by penetration, by vibration within one image, by deformation between images under varying
   loads, or by changing the surface's material, and in no other way. Penetration: X band reaches about 30 cm into
   the driest sand (TanDEM-X), L band 1 to 2 m (SIR-A's buried valleys), P band up to 5 m (Biomass, 2025); limestone
   conducts 50 to 200 times more than dry sand, so the bench chamber's roof is at least 347 dB down at X band and the
   claimed structure 34,000 dB, against about 30 dB an image spans; only P band could reach a shallow room, inside
   its 25 m resolution cell with the surface. Vibration: the smallest motion ever measured from orbit (0.66 mm/s RMS
   on a corner reflector, 2026) is 2 × 10⁸ times the chamber's imprint, which reaches it only under shaking of
   0.6 m/s during the pass, over reflectors placed in advance. Between passes: the largest load, the year's heat
   wave, moves the ground over the chamber by 2.4 µm, which needs about 8 × 10⁵ passes at 0.5 mm each and is 180
   times smaller than the ground's own seasonal breathing; tides, air pressure and the water table give 0.03 to
   0.2 µm. The surface's temperature changes by 0.1 mK. Any processing, gates and learning included, is a function of
   the record and cannot restore what these routes do not carry.

Any one of 1, 3 and 4 is sufficient on its own; they fail for independent reasons (the estimator, the
physics of the signal, and the algebra of the depth axis). R16's resolution bound, which rests on the
same trade between a look's duration and its resolution (1.10 m s at Giza, P2-01), is a fourth. 8 extends
the question beyond the published method to every route a radar record offers.

### Motion measured from orbit, 2026

Two papers from the Strathclyde and Trento groups (one co-author, Clemente, also wrote the 2020 tracking papers with
Biondi) measure vibration from single spaceborne images against synchronous ground truth, and are the strongest
evidence that single-image micro-motion is real:

- Vattulainen et al. (2026), *Assessment of spaceborne SAR micro-motion measurement for vibration-based SHM*, IEEE
  Access 14, 6043–6064, doi:10.1109/ACCESS.2026.3652346. Corner reflectors on shakers, set in open grass away from
  other reflectors, nine Umbra spotlight images (5–16 s), vertical motion at 1–4 Hz tracked by sub-pixel offset
  tracking of the brightest pixel through overlapping sub-apertures. Peak radial velocities from 95.5 down to
  1.4 mm/s; RMS errors from 20.6 mm/s (22%) to 0.9 mm/s (63%); frequencies found in every test, down to 0.10 mm RMS
  displacement (0.66 mm/s RMS velocity).
- Lotti et al. (2026), *Monitoring bridge vibrations via spaceborne SAR micro-Doppler*, Structural Control and Health
  Monitoring, doi:10.1155/stc/3858095. A steel footbridge over the Clyde, chosen because its steel and sharp corners
  are bright; peak radial velocities 0.5–2 mm/s, errors of about 1 mm/s, dominant frequencies found at 0.06 Hz
  resolution.

Both are what P2-03 predicts works: an isolated, bright, point-like target's envelope follows its own motion. Neither
tracks natural ground or speckle, and the smallest motion measured is 13,000 times the ground's whole trembling at
Giza and 2 × 10⁸ times the chamber's imprint on it. Their measured floor is some ten to a hundred times above the
corner-reflector threshold P2-05 computes, so P2-05's thresholds are generous to the claim.

### Corrections to the first investigation's record

- **R9 measured the looks, not the ground.** Its coherence between looks of one image is a function
  of the image's power spectrum alone: 1 − Δ/W for looks W wide and a flat spectrum. With its 1.96 s
  looks that is 0.98 s to half and 1.76 s to a tenth, which is what it reported, identically on the
  pyramid and on desert. Simulated desert with no decorrelation of any kind reproduces its curve to
  ±0.005, and the real crop's own power spectrum reproduces its coherence to 10⁻⁹. The "coherent
  window" is the look width one chooses; no experiment here has measured how stable the ground itself
  is across a dwell. The obstruction the first investigation drew from R9 still holds, for a deeper
  reason: overlapping looks share their data and disjoint looks share no speckle.
- **R4's 77 µm/s is the scatter of an estimator that does not respond to the ground's motion**, not a
  sensitivity to it: its complex common-reference registration recovers 1 ± 1% of real motion on real
  texture (P2-03). For bright point-like targets a magnitude tracker does follow motion, with a
  per-look scatter of the same order, which is what the micro-motion literature measures on bridges.
- **R7 never recovered its injected motion.** Its committed run (23 September: the full-band bank, 0.2,
  0.5 and 1 Hz, amplitudes up to 1 cm) has recovery gains within about ±0.1 of zero at every amplitude
  above the smallest, so its finding's "it works" is not what its numbers show. That is what P2-03
  predicts: the motion was common to the whole crop, which for speckle is undetectable in principle, and
  complex registration of one image's looks locks at zero. The script has since been changed to the
  coherent bank and to 1, 2 and 3 Hz without being rerun; those frequencies sit on the nulls of its
  1.96 s looks, which return 10, 4 and 1% even for a bright point (P2-03). R7 should not be cited as a
  positive control.
- **R10 inherits R9's identity.** It graded pixels by brightness and computed the same complex
  coherence between looks, which is set by the spectrum whatever the pixels are. "Bright targets
  decorrelate as fast as desert" is therefore not established by it; a lone bright scatterer's
  envelope in fact follows its motion across the whole dwell (P2-03).
- **R11's positive control never passed through the tracker.** Its synthetic plane wave was added to the
  velocity series after tracking, so its detection says the stacking works, not that the image's pixels hear a
  wave; and its band, 0.9–3.5 Hz, is averaged away by its own 1.96 s looks (a bright point swaying at 1, 2 and
  3 Hz comes back at 10, 4 and 1%, P2-03). P2-09 puts the wave into the image instead: the sensors are deaf.
- **The platform moves 187.6 km over the processed aperture**, not 168 km; 168.8 km is the span
  between the reference bank's first and last sub-aperture centres (P2-01).
- **The flat-earth simulator of the first chapter put the Doppler rate 8% off** for this geometry.
  Its experiments used generic parameters and none of their conclusions rest on the real rate.

### What is not established

- The imprint is computed for one chamber in uniform bare rock, quasi-statically. Layering, the sand
  cover and the plateau's 0.6 Hz resonance are not modelled; site amplification of a few times would
  not change the shortfall's order of magnitude. The traffic case is only a guide, since at 15 Hz the
  waves are not long compared with the depth.
- The magnitude tracker is the best of those tried, not a proven optimum. The argument that nothing
  can do much better for ordinary ground is the distributional one in 1: what survives in speckle is
  the apparent brightness a displacement gradient produces, which for this imprint is of order 10⁻¹¹.
- Giza's own ambient level is unmeasured in the open literature; the microseism level is from
  Kottamya, 67 km east (M1-01), bracketed by Peterson's global models.
- The known-chamber test at Khufu on two real passes (2022 and 2025) is run locally and not yet published. Given 1
  and 4 it could only confirm, not change, the result.
- P2-11's penetration bounds use dry sand's measured penetration as an upper bound for limestone, and a limestone
  conductivity from the standard GPR table; no microwave loss measurement of Mokattam limestone itself is used.
