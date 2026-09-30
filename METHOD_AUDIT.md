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

**Superseded in part, 28 September.** The agreement holds for bright, isolated targets, which is what
the literature measures (the 2026 papers: corner reflectors and a steel footbridge at 0.7 mm/s and up).
It does not hold for ordinary ground: P2-03 shows R4's 77 µm/s is the scatter of an estimator that
recovers about 1% of a patch of ground's motion, not a sensitivity to it, and Giza's trembling is some
13,000 times below the smallest motion measured from orbit. Even perfectly measured motion would not
show a chamber: the published depth is the motion's frequency (P2-06), the ground's trembling comes in
waves kilometres long, and passive imaging with perfect instruments and the rock known exactly reaches
only 2.2× with unlimited listening (P1-04).

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

1. **Complex registration of one image's looks registers little of the ground's motion (P2-03).** Two
   looks cut from one image share their overlapping spectrum, so for whole images and real, non-negative
   masks their cross-spectrum is |Y|² H_A H_B and their correlation peaks at zero lag whatever the ground
   does; for patch registration, which mixes frequencies, this is a strong bias to be measured rather than
   an exact identity (an independent check found two of seven cropped static patches off by one 0.01 px
   step). A 20 m block of ground swaying at 2 mm/s comes back at 3 ± 1% on simulated desert and 1 ± 1%
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
   plateau with no monument draws them too (P2-12 counts them). Along a column the power rises and falls at the axis's
   resolution (13 m relabelled): the bands. At z = 0 and at every repeat depth (648 m and 1,296 m relabelled) all
   steering phases coincide and the power is the pixel's mean offset squared (ratio 0.994): the blocks. Below 648 m
   the picture repeats.
7. **A radar image is not a dense seismic array (P2-09).** With the motion put into the image itself (a 60 m, 0.15 Hz
   test wave), 3,264 virtual sensors read by complex correlation recover 0% of the motion under them and by magnitudes
   −6%. Still sensors do average down as n^−0.50, but a deaf sensor stacks to a precise zero. Where magnitudes hear
   (a fifth, over real point-like ground), the chamber would need about 2 × 10¹⁹ independent sensors over it; there
   are 3,602 there, 7 × 10⁷ in the whole image and 4 × 10¹⁴ on all of Earth's land. The ground's own trembling would
   need 1,186 times more than the image holds.
8. **As estimated, no route carries a Giza room to an X-band radar in useful measure, nor the claimed structure to any
   radar (P2-11).** An echo depends on how strongly each scatterer returns the wave and where it is while watched; we
   group the ways a room could reach the record into penetration, vibration within one image, deformation between
   images under varying loads, and changes of the surface's material. These are estimates with stated assumptions,
   not a proof that no route or processing could ever carry anything. Penetration: X band reaches about 30 cm into
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

9. **On the real image, the method's picture is the surface's (P2-12).** Run on the 2025 pass over Khafre, Khufu,
   Menkaure and three stretches of open plateau, the published pipeline's power at every pixel follows that pixel's
   trajectory energy (1.0000 on every patch). Every patch draws columns; the pyramids, whose surfaces are busier, draw
   more of the strongest (5 to 9% of their pixels are among the noisiest 5% of all six patches, against 1.6 to 3.6% on
   open plateau), while their typical pixel wanders as much (medians within 7%). Its mean depth profiles are alike over
   the three pyramids (0.88 to 0.91),
   whose surfaces are alike, and a pyramid's differs from open plateau's (0.61 to 0.83) by no more than two stretches
   of open plateau differ from each other (0.44 to 0.72). The image itself is resampled onto the ground through the
   product's RPC over terrain and pyramids, the residual offset fitted against predicted brightness (1.0 m and 2.6 m;
   the fit's correlation is low, 0.03, and the alignment was checked by eye against the pyramids' bases).

10. **Stricter gates pass motionless images, and a simple vibrating target's trace does not pass them (P2-13).** A public
    implementation of a gated reconstruction (317 pairs of 32,330 Hz masks 404 Hz apart over the central 12 s,
    registration rounded to 0.01 px, arcs subtracted; a window of 50 pairs kept only if its two shift components fit an
    ellipse at modes 1 to 10 with adjusted R² ≥ 0.25, a minor semi-axis ≥ 0.005 px and an axis ratio ≥ 0.1; ten same-mode
    windows among the central starts; four or five contiguous positions at one start; the depth focus only there) is run
    unchanged with its frozen X13 profile, its source hashed. Along seven loops round Khafre's faces, 20 to 95 m above
    the base (2,828 positions), it passes 5 positions (one feature) on the real 2022 image, 13 and 8 (three and two) on
    two motionless twins, none on the same loops over open desert in the same image, and 9 with vibrations of 2 mm/s
    along the line of sight planted in the image at 0.26 Hz and 3.66 Hz: the real image's feature again, and one new
    feature on the 0.26 Hz stretch at mode 3, which stands for 3.2 Hz. Each image of its pairs spans 6.0 s and the two
    share 98.75% of their spectrum: 36% of the shifts are exactly zero and 75% within one step, and the planted
    vibrations move them by 0.006 and 0.002 px rms. (The run also quotes a "full tracker" at 0.063 and 0.016 px and a
    boxcar gain of at most 4.7%, under 10⁻⁵ at modes 5 and 10; both are rectangular-average values, not this
    registration's response, which P2-15 measures instead.) A mode is m cycles in 50 pairs (0.935 s), so the modes stand
    for 1.07 to 10.7 Hz; a 3.66 Hz tone is expected to win mode 3. A simple target moving along the line of sight moves in
    each image along azimuth by its velocity times 84 s and in range by its displacement, so its trace is at least 3,000
    times as long as it is wide; the gate keeps ellipses at most 10 times, so for such a target what passes is shaped by
    the registration's noise; more complicated scenes would need their own forward simulation. 31 of the 35 passing
    positions have their best depth within half a cycle of the mode times 6.17 m or of its mirror (308.7 m less): that
    is what the depth fit computes after the gates keep sinusoids, not a significance test, and the 35 positions form a
    handful of clusters. The twins and the desert are stress tests, not a calibrated false-alarm rate. At the first gate the real image passes 2.5 times the twins' windows and open desert 1.9
    times (36% of shifts exactly zero against 41%); planting changes that count by 0.7%, and passing early does not
    predict passing at the end. The lines are the lab's, on the faces, placed through the product's RPC over the site's
    surface (offset (44, 4) px from predicted brightness, correlation 0.03); the implementation's own lines, round the
    base, are 11.

11. **On the reconstruction's own lines round Khafre's base, nothing passes at its own settings (P2-14).** The
    profile's frozen geometry (Khafre's base corners, Petrie's as placed by Dash, projected through the RPC at
    H = 70.5 m above sea level given as an ellipsoidal height) and its EGM96-corrected version (h = 85.9 m, about 50 px
    nearer in range; the lab's projection of its corners reproduces the file exactly) are run unchanged, with pixel
    support reported from one to five. At the profile's own support of four or five, nothing passes on the corrected
    lines in the real image, two motionless twins, the planted image or the same square over open desert; the
    uncorrected lines pass 4 positions in the real image and in the planted, none in the twins. At support one the
    corrected lines pass 69 positions in the real image, 28 and 29 in the twins, 67 planted and 53 over open desert; at
    two, 12, 4, 0, 8 and 4. Planting adds and removes a few positions (18 start passing and 15 stop, over all supports):
    in the 0.26 Hz stretch at modes 1 to 3, which that frequency cannot produce, and in the 3.66 Hz stretch one position
    at the nearest mode. 352 of the 395 positions passing at support one have their best depth within half a cycle of
    their mode times 6.17 m or its mirror, as a frequency fit must. The pipeline's own steering uses a Doppler rate of
    5,390 Hz/s at Khafre, so each 32,330 Hz image holds 6.00 s of every pixel's pass, and puts the four sides'
    zero-Doppler moments within 0.03 s of each other. The corrected file's projection.height_m is a text note, which
    the steering reads as a number; it is set to the corners' ellipsoidal height.

12. **What the six-second pairs register of a vibration, measured (P2-15).** The reconstruction's own frequency_masks and
    register_one are applied to simulated images on the 2022 geometry (a bright point alone, the point 30 dB over
    clutter, a 20 m block of clutter), with line-of-sight sinusoids of 2 and 20 mm/s at the modes' frequencies and at
    0.1, 0.26, 0.5 and 3.66 Hz, put in through the Doppler-time mapping. Measured through the reconstruction's own masks and registration, on simulated images of the 2022 geometry with motion put in through the Doppler-time mapping, the six-second images do carry a vibration, but not as a boxcar average: read against its motionless twin, a lone bright point's image moves by 0.062 to 45 times the boxcar value, depending on frequency (at 2 mm/s and 0.5 Hz, 0.129 px where the boxcar gives 0.0085 px). The pair registration then catches only part of the shift between its two images: 3% to 46% for the bright point and 0.6% to 2.8% for a 20 m block of texture, where that shift is resolvable. At 2 mm/s, some 40,000 times Giza's microseisms, the pair shifts change by at most 0.011 px for the point, about the 0.01 px rounding step, and 0.0017 px for texture. The moving point passes 0 of the first gate's windows at every frequency and both speeds, its trace a line; with clutter or texture the traces large enough to measure are at most 0.06 as wide as they are long, under the gate's 0.1, while still texture alone passes 24 windows.
    The injection is the lab's usual approximation (a phase history on a focused image, without range migration,
    occlusion or changing reflectivity); a pulse-level simulation or shaker data would test it.

13. **Over the known bench, the gated reconstruction draws the same whether the chamber is there or not (P2-20).**
    Synthetic products in the ICEYE layout, read by the reconstruction's own loader and run unchanged (51 east-west
    lines of 51 positions, 2.4 m apart, over the one-chamber bench), hold speckle shaken as Giza shakes: the measured
    microseisms at 0.2 Hz and the measured 1-3 Hz level placed on the gates' own second mode (2.14 Hz), where the
    method looks. The chamber's imprint (P2-04's map, tapered at its edge, times each band's strain) rides on that
    motion at its real level, boosted ten thousand and a hundred million times (3.6 mm/s over the chamber), beside a
    null: a random perturbation the size of the real-level imprint at the same pixels. Stated before the run: a
    response would gather the changed positions over the chamber (within 12 m of its axis, 3% of the positions) as the
    boost grows. Against the product without the chamber, the fit scores change at 178 positions with the chamber (6
    of them near it), 182 with the null (4), 242 at ten thousand (6) and 300 at a hundred million (9), where chance
    puts 5 to 9 near; 275 to 280 positions pass at support one in every product, none at the profile's own four. The
    changes never gather over the chamber. Changes in the last digits of single precision move a position's fit score
    by up to 0.43 (the chamber) and 0.50 (the null) on the method's scale of 0 to 1. An earlier version of this run
    shook the ground only at 0.2 Hz, below the band the method reads; it was stopped and replaced before publication.

14. **What the gated reconstruction draws across Khafre does not survive a change of pass (P2-21).** The same square
    (9,595 positions, 3 m apart) run three ways with the lab's command, each beside a motionless copy: the 2022 pass
    with its lines east-west and north-south, and the 2025 pass east-west. Stated before the run: a response to what is
    below would agree across the runs beyond chance and beyond an image's agreement with its own motionless copy. The two
    layouts on one image share 4.9 times the chance number of passing positions, with best depths agreeing at 86% against
    26% by chance: one image's texture decides both. The two passes share 1.33 times chance, depths agreeing at 27%
    against 25%; the 2022 image and its own copy, which has nothing under it, share 1.32 times chance and agree in depth
    at 38%. Nothing survives the change of pass beyond chance.

15. **One shaking, read by geophones and by the satellite (P2-22, P2-23).** A vertical vibrator 30 m from the bench's
    chamber, at the gated reconstruction's tenth mode (10.69 Hz), solved to steady state by the elastic solver with and
    without the chamber (and with a larger, shallower one): the chamber changes the surface motion by up to 1.4% (6.5%
    for the larger room), within about 19 m of its axis. Read over the same 24.6 s by 49 geophones in Kottamya's measured
    ambient noise (or Peterson's high-noise model), and by the satellite on synthetic products of the 2022 geometry: 21
    corner reflectors 50 dB over speckle laid across the chamber in ground range, tracked through 20 ms looks by a
    tracker calibrated on reflectors of known motion (noise per quadrature 1.5 mm/s at 50 dB, 4.4 at 40 dB, 9.3 at 30 dB).
    The favourable detector knows both answers. Geophones find the chamber with 11 kN (1.1 MN at a noisy site); the
    satellite needs 2.5 x 10^10 N, 2.2 million times more, the ground over the chamber then moving at 19 cm/s. Open
    plateau in the real 2022 image offers nothing to track (its brightest point 12 dB over the speckle; 33 points of
    10 dB or more where speckle alone gives 32). Run whole, the chain recovers 99% of the vibrator's own motion where the
    reflectors move at up to 50 mm/s (the measurement control) but only 4% at the idealised boundary force, so the
    satellite's boundary is a floor; detecting a 0.3% change on the vibrator's motion would also need the tracker's gain
    known to about 0.1%. Limits: one chamber, one source position and frequency, the speckle around the reflectors held
    still, a flat bench; the record length is the dwell's, and longer geophone records would lower their boundary further.

16. **Where the data hold the chamber, the gated reconstruction does not read it (P2-24).** Synthetic products of the
    bench shaken by P2-22's vibrator, the row of corner reflectors and the real image's bright points in them, the motion
    entered exactly (Bessel terms of the phase history), at the force where P2-23's favourable detector finds the chamber
    (2.5 x 10^10 N) and at a thousandth of it; each with the chamber, without it, and with random noise the size of the
    chamber's change at the same pixels. Strong shaking at the method's own tenth mode makes it pass more positions (633
    at support one, 67 at four, against 275 and none under Giza's own shaking in P2-20). The chamber changes its scores at
    646 positions, 32 within 12 m of the axis; the random noise at 648, also 32. At a thousandth of the force, 454 and 456.
    The chamber changes the method's picture exactly as noise of its size does.

17. **Within its model, one image of ordinary ground holds almost nothing of the chamber, however it is read (P2-25).**
    Items 1 to 16 test methods by running them; this bounds every method at once, within a stated measurement model. An
    image is a function of the echo from what the wave reaches, the surface, and the surface's motion during the pass,
    and a method is a function of the image (the data-processing inequality). For any test, detection rate minus
    false-alarm rate is at most the total variation TV between the image's distributions with and without the chamber;
    with the two equally likely, accuracy exceeds a coin toss by at most TV/2. For fully developed speckle the finite
    change is bounded rigorously, not to first order: whitening the covariance change into E,
    KL <= |E|_F^2 / (2 (1 - |E|_F)) with |E|_F <= sqrt(F) + eps, F the Fisher information on the image's own pixels and
    eps an explicit power-series remainder (`sim/sarsim/information.py`); TV <= min(sqrt(KL/2), sqrt(1 - exp(-KL))).
    Motion shared by the scene holds exactly nothing; slow smooth motion leaves only the texture's stretch along track
    (R/V_s = 79.3 s times the velocity's gradient). Checked: the sum against a brute-force covariance (6e-16); the exact
    finite KL on small images at most 1.006 of F/2, the rigorous bound holding in every case;
    detection measured apart from estimation: on a small image the exact likelihood-ratio test reaches
    0.20 and 0.54 (the score test 0.19 and 0.52) where Pinsker from the exact KL allows
    0.25 and 0.75: the ceiling is a bound, 1.28 and 1.39 times what the best test
    achieves, not an attained value; the score's spread and shift (estimation, not detection) match sqrt(F) and aF; and the
    Doppler-to-time relation against pulse-by-pulse physics, the same information to
    3.6% from 5 to 150 Hz on a short aperture, the model's excess an edge effect falling to
    0.07% as the time-bandwidth product grows (the real dwell's 3e+06), on the side of more information. On the
    real dwell, generously (the method told the ambient shaking, no receiver noise, no attenuation), under Giza's regional
    microseism level (Kottamya, 67 km east; not an upper limit at the pyramids) detection minus false alarm is at most
    1.3e-09; 9.3e-08 at the noisiest stations on Earth; 1.3e-08 and 2.8e-08 at the regional 1-3 and 3-8 Hz levels;
    with P2-26's dynamic imprint and the scattered wave carried unattenuated to the scene's farthest corner
    (3.54 km), 1.6e-03 under the FTA's urban background all pass and 2.2e-02 for a truck over a bump 15 m
    away all pass (9.0e-03 and 0.12 over a room whose roof is 5 m down); a corner reflector over the imprint's
    peak 2.4e-07 under the microseisms; a genie told the reflectivity, 1.9e-06 (with the truck 0.48:
    open). These ceilings take the reference covariance as the identity (the background motion shared); P2-30's looser
    certificate includes the background. The local extrapolation factor (how much stronger the imprint would have to be for
    the linearised information to reach a reliable detection) is a scale within the model, not a physical shaking
    requirement; P2-27 states the requirement properly.

18. **In the rooms modelled, the responses are broad; faster shaking reaches a room better (P2-26).** The lab's elastic
    solver, with and without the chamber, under a surface wave from the side, P and S waves from below, P2-22's favourable
    room (10 m, roof 5 m down) and the bench room 30 m down, in rock without attenuation, 6-120 Hz (240 Hz at 0.5 m on a
    smaller model). The bench room's imprint on the line-of-sight velocity rises to 0.16 of a surface wave's motion
    near 79 Hz; 0.31 (P) and 0.34 (S) from below; 0.043 for the room 30 m down; the thin roof moves the ground
    over it up to 1.98 times the passing motion near 65 Hz (2.15 on the finer grid). The width of the largest
    response's envelope is reported as such, not as a modal Q, and marked where the model's band cuts it (the surface and
    S cases). Because 0.3-0.5 s records cannot exclude a long-lived narrow mode, each room was also struck on its ceiling
    and recorded for 4 s beside intact rock: the largest modal Q resolved from the spectra (0.25 Hz bins) is 2.2, with no
    peak within 30 dB of the largest narrower than a bin; the energy's decay time is set by the model's residual floor
    (intact rock with no room 6 to 16 s; a room at most 1.07 of it), so a mode weaker than that floor is not
    excluded. Over 4 s of a passing surface wave, 5.4e-05 of the imprint's energy over the thin roof (6e-09 over the
    bench room) comes after the wave has passed. At 6-12 Hz the imprint is 0.96 to 1.82 of P2-04's static answer.
    Halving the grid changes it by at most 6% below 120 Hz (11% for the thin roof). Only these rooms,
    excitations and windows are covered: fractured, layered or coupled structures, the air inside a room and the rock's
    nonlinear response are not modelled.

19. **What one image would need: an energy requirement for every pattern, and a curve for smooth ones (P2-27).** The
    bound turned around, for the Giza dwell (and the 2022 acquisition), fully developed speckle without a reference image
    or receiver noise, motion coherent for the whole pass at its worst phase, and a target of 95% found at 5% false
    alarms. For every pattern, from Theorem C: a room-specific motion whose energy, the integral of max_t |u_los|^2 over
    the ground, stays below 9.8e-10 m^4 (Giza's background motion in the reference covariance; 1.7e-8 m^4 with the
    reference taken as the identity) is excluded whatever its shape, frequency or time history, an rms of 1.0 um over
    1,000 m^2. For one smooth (Gaussian) pattern, exactly within the model (each line's KL from its covariance's
    determinant, the range band widened): at least 3.0 mm at 0.2 Hz, where only the gradient counts, and 47 um at 80 Hz over
    279 m^2. The smooth curve binds only that family: at 0.2 Hz the Gaussian holds 4.7e-4 of what the most informative
    pattern inside its envelope could (the same envelope times a cosine of 4 m period holds 39 times as much, of 1 m period
    489 times); at 20 to 80 Hz it holds half, and the pattern matters little. The bench room's modelled static imprint under
    the regional microseisms falls short of the energy requirement by 4.6e5 in amplitude (8.3e8 against the smooth curve).
    Under a truck over a bump 15 m away all pass, with the scattered wave carried unattenuated to the scene's farthest
    corner, the bench room's dynamic imprint falls short by only 2.5 to 3.7, and the room under a 5 m roof not at all
    (0.45 to 0.97); counted within 39 m, 33 to 35 and 5.2 to 8.1: under strong shaking the conclusion rests on how far the
    scattered wave carries. A change of the ambient motion's phase alone moves the ground by at most twice that motion, so
    such a mechanism is excluded only while the area it changes stays below E_min / (2a)^2: 5.5e4 m^2 under the regional
    microseisms (the room's footprint is 128 m^2), 11 m^2 at the noisiest stations, 2,131 m^2 under the truck at 80 Hz. It
    says nothing of routes outside the model or of mechanisms with gain.

20. **The published picture's depth does not follow the room; its horizontal change follows the image's (P2-28).** The
    published chain run whole (P2-07's) over the bench room 15 m down, the same room 30 m down and a granite block 15 m
    down, each imprint by the solver settling under strain (the 15 m room's peak 1.000 of P2-04's), boosted to 0.5
    and 2 rad, with lambda_s = 1.45 m so the axis passes below both rooms. At 2 rad each room puts 16% of its change in
    its own depth band where 15% of the axis lies; both peak near 14 m. The granite block's change is the hollow
    room's with its sign flipped (-0.45), as their marks on the ground are (-0.88). Horizontally the change gathers over the
    room (2.6 times elsewhere) and a random perturbation of the same pixels gathers almost as much (2.4): in a
    paired comparison the picture changes where the image changes. The blind question is P2-31's.

21. **The oracle told everything but the answer (P2-29).** Given the reflectivity, the background motion, the shaking and
    both models, the raw echoes are Gaussian with the receiver noise's covariance, and the best test's detection minus
    false alarm is exactly 2 Phi(Delta/2) - 1; a method told less does no better, whatever the texture (only independent
    uniform phases needed), and the background motion cancels exactly. Pulse by pulse on the real dwell the formula holds
    (0.993 +- 0.025; 1.011 +- 0.030 with the background ten times larger). At 30 dB per cell (above the
    brightest natural ground over ICEYE's best Dwell NESZ, -26.7 dB): 1.6e-06 under the regional microseisms,
    1.2e-04 at the noisiest stations; the imprint would have to grow 8,415 times (117 at the noisiest stations) even after assumed
    allowances (local level 10x the regional, site amplification 3x, 10 dB more SNR). Strong shaking close by, known exactly,
    is not excluded: the oracle's bound reaches 0.41 with a truck bouncing beside the bench room all pass, and 0.99 beside
    a room under a 5 m roof, with the scattered wave carried unattenuated to the scene's corner (0.24 within 32 m alone); 14
    cases are not excluded after the allowances. These are upper bounds (energy bounds, averaged over phases): a large one
    shows that this argument does not exclude detection, not that any detector achieves it. Depth: small rooms 15 m down in any non-negative mix leave 71% of the
    bench room 30 m down's mark unmatched (78% for the bench room 15 m down; two small rooms solved together agree with their
    sum to 3.7%): in uniform rock the mark's shape carries depth, and telling those depths apart takes about
    1.4 times detection's signal. The sought depth obstruction was not found.

22. **Finite certificates, loose but simple, the background included (P2-30).** An independent derivation's theorems
    (`sim/sarsim/finite.py`, its checks ported): processing contracts TV; the exact complex Gaussian divergence; a
    certificate straight from a phase-energy envelope (|C0^-1/2 dA|_F <= q gives |E| <= 2q + q^2), with a proved covariance
    floor for the common background (0.056); the oracle; two-world depth risk. On the lab's runs: at most
    5.017% found at 5% false alarms under the regional microseism envelope with the background included (TV
    1.7e-04, where P2-25's tighter bound gives 1.3e-09); an arbitrary phase modulation within the background envelope
    on 52,600 pixels at most 9.55%, over the whole image undecided; the thin-roof room with a truck-level
    harmonic at most 27.3% with the reference covariance taken as the identity, still excluded with 11% more phase (89.1%
    before the far field was taken direction by direction; the old bound took four times the average), and not excluded
    with Giza's background motion in the reference.

23. **Telling one layout from another is bounded as detection is (P2-32).** For any two layouts A and B,
    TV(A, B) <= TV(A, none) + TV(none, B), and the closest admissible pair of depth variants is no easier than the pair at
    the same depth, so the same-depth bound covers every allowed vertical transformation. Each layout's imprint is the
    lab's solver (statically for the ambient levels, dynamically for a lorry). The tight bound and the certificate take each
    pair through empty ground (subtracting two layouts' motions does not remove the reference covariance); the oracle uses
    the difference itself, exact there. Under Giza's regional microseisms the oracle tells an L-shaped from a straight
    tunnel at most 1.5e-06 and a column with a spiral from a plain column at most 2.0e-06 (tight 2.6e-09 and 9.1e-09; the
    certificate with the background at most 5.110% at 5% for every pair); discrimination needs 1.2 to 4.1 times detection's
    signal. With a lorry bouncing beside the layouts all pass, known exactly, the oracle's upper bound reaches 0.82 for an
    L-shaped against a branching tunnel, still excluding 95% at 5% until the signal grows 1.2 times: unresolved, as
    detection is. Only these shallow bench-scale layouts in uniform rock are covered.

24. **Where the data hold a layout, a reference detector names and places it; the published picture does neither (P2-31).**
    A room, an L-shaped and a branching tunnel, each imprint from the lab's solver, the image from the synthesizer, the
    unrelated ground and shaking matched across layouts and controls (no cavity, a motionless copy), three realisations,
    read at the real level and at labelled amplified diagnostics (0.3, 2 and 20 rad). At 20 rad (the positive control) a
    reference detector told the shaking names 100% of layouts and, scanning every placement, puts its peak a median 0.4 m
    from the layout's centre (controls 36 m). At 2 rad it names 78% when told where to look (p = 0.008), but its response
    there stands only 1.9 spreads above the scatter over placements and a scan peaks 34 m away. The published method's plan
    map, read blind, names 56% at 20 rad (p = 0.14 against guessing, nine images) and centres its top places 5.6 m from the
    layout's centre against the controls' 6.2 m; its largest value separates cavities from controls with AUC 0.59. Its
    change against the same ground without the cavity, a diagnostic it never has, names 67% at 0.3 and 2 rad (p = 0.04 each,
    uncorrected for a dozen comparisons) and centres 5.2 m away against the motionless copy's 12.1 m (the picture changes
    near where the image changes, as in P2-28), and 33% at 20 rad. After the predeclared depth family its change peaks
    1.5 m from the true depth at 20 rad and the motionless copy's 1.2 m: with 12 candidates, closeness says little. At
    the real level every reader is at chance. Location was first scored as the distance to the nearest footprint point;
    because the tunnels reach within a few metres of almost any central place (chance 3.4 m), the centroid distance was
    added after the first run and is the one reported. Amplified levels are diagnostics, not physical predictions; nine
    images per level is small, and a physically strong positive control (a lorry or vibrator beside the layouts) is next.

These reasons differ in kind and strength. 4 is algebraic and holds whatever the ground does: the depth axis is a
frequency axis with a repeat and a mirror, and no validated model turns it into a physical depth; that makes depth
uncalibrated and non-unique, and does not by itself decide horizontal location or shape; P2-31 tests those directly, and the published picture recovers neither where a reference detector recovers both (24). 1 and 12 are
measured, for particular estimators and configurations; 13 is a known-truth test of one of them, on synthetic products. 2, 3 and 8 are models with stated assumptions (a small chamber
in uniform rock, quasi-static strain, the regional microseism level, dry sand's penetration as a bound for limestone).
Together they show that the methods tested do not establish subsurface imaging; they do not prove every conceivable
route impossible. R16's resolution bound (1.10 m s at Giza, P2-01) limits short looks; strong prior information could
in principle detect weak changes below it. 17 is of a different kind: a theorem that bounds every way of reading one
image at once, tested or not, including methods with strong prior information, within a stated measurement model:
ordinary ground as fully developed speckle with no reference image, the Doppler-to-time relation (checked against
pulse-by-pulse physics), bright points bounded separately, and any scene by a genie that knows the reflectivity. Its
inputs, the chamber's imprint (2, 18) and the regional ambient levels, are models. It is a bound: on small test images
the best test comes within a factor of about 1.3 to 1.4 of it, not onto it. 19 turns it into a conditional requirement: on the energy of the motion a room would have to add, for any pattern, and
on its amplitude for smooth ones. None of it covers penetration at longer wavelengths, changes of
reflectivity during the pass, texture between speckle and bright points, or many images; and none of it has been
checked against measured ground motion and radar data taken together.

### Independent review, 28 September 2026: corrections adopted

An independent review of this repository (commit bd76ceb) and of the reconstruction's public v1.8 code reproduced
several component checks and found overstatements on both sides. Adopted:

- **Scope.** Universal impossibility language ("cannot", "every route is closed") is replaced by scoped claims: the
  methods tested do not establish subsurface imaging; the depth model is unvalidated and aliased; the ambient-vibration
  benchmarks are models showing severe sensitivity problems.
- **Geometry.** P2-13 ran the lab's loops up Khafre's faces after misreading the reconstruction's 70.5 m (the base's
  height above sea level) as a height above the base; P2-14 runs its own base lines, with hashes and manifest.
- **The boxcar.** sinc(fT) is a rectangular average's response, not the registration's; P2-03's own point target
  departs from it (21% against 2% at 0.5 Hz with 1.96 s looks), and P2-15 measures the six-second configuration.
- **The shared-spectrum identity** is exact for whole images and a measured bias for patch registration.
- **3.66 Hz** is expected to win mode 3; a mode alone does not show whether a planted motion came back.
- **Twins and desert** are stress tests, not calibrated nulls; open desert is not surveyed ground.
- **Depth following the mode** is what the fit computes after the gates select sinusoids, not a significance result.
- **Models are conditional**: the imprint (6.6 × 10⁻⁵) and the budget assume a small chamber in uniform rock,
  quasi-static strain and a regional ambient level not measured on the plateau during a pass.
- **Injections are approximations**: a phase history on a focused image, reasonable for diagnosis, not a pulse-level
  simulation of a vibrating pyramid.
- **The sound wavelength**: the review notes that the 2022 paper implies 0.48 m in its text (12,500 Hz at 6,000 m/s)
  and 0.24 m in a displayed expression; the depth scale is proportional to it and the identity holds either way.
- **The orange volume** in the viewer is the paper-style pipeline (the 2025 image, 50 pairs, no gates, focused power,
  relabelled to repeat at 648 m, smoothed), not the gated reconstruction; the viewer now says so.

None of the review's corrections was rejected; where it asked for a measurement (the six-second transfer function)
the measurement was made (P2-15).

### The organising question, widened (29 September 2026)

The owner asked that useful underground information, horizontal location and shape included, be the goal even where
absolute depth is ambiguous. Adopted: presence, horizontal location, shape, relative vertical geometry and absolute depth
are separate claims. The depth-alias results (P2-06, P2-08, P2-28) show that the published depth is uncalibrated,
repeats and mirrors, and that changing the sound wavelength rescales it without changing the fit: a calibration
ambiguity and non-unique depth, not by themselves evidence that recovered horizontal positions or shapes are artefacts.
The evidence that the published and gated pictures carry no underground information at the real level rests on the
controls instead: the same picture with and without a known chamber (P2-07, P2-20), patterns in motionless copies and
open ground (P2-13, P2-14), pictures that do not survive a change of pass (P2-21), and a picture that does not follow the
chamber where the data demonstrably hold it (P2-24). The finite bounds (P2-25, P2-29, P2-30) bound presence within their
model, and location and shape with it, since telling two layouts apart is harder than telling either from none; a bound
that excludes 95% at 5% is not by itself a proof of no prospecting value, which is excluded only where the bound on
detection minus false alarm is itself small. The shape test (P2-31, item 24) found that where the data hold a layout
a reference detector names and places it and the published picture does neither.

### What the proof can and cannot claim (29 September 2026)

Three goals, kept apart (after the second review's framing):

- **Limits within a stated measurement and physical model:** achieved, conditionally. Detection (P2-25, P2-29, P2-30) and
  shape discrimination (P2-32: a branching tunnel or none, an L-shaped or a straight tunnel, a column with a spiral or
  without), each finite and algorithm-independent, the same-depth bound covering every allowed vertical transformation.
- **Limits for a defined class of realistic archaeological settings:** a feasible research goal, needing measured ground
  motion at the site, measured coupling over known voids, a qualified observation model and noise floor, and experiments in
  the regime the bounds leave open (strong shaking close by, known exactly).
- **That every conceivable single-image application has no archaeological value:** too broad to claim from this work.

The useful outcome is the boundary: which conditions rule out useful recovery (within the model, the ground's own
trembling) and which deserve experiments (strong, known shaking near a shallow structure; mechanisms with gain once
specified).

### Second review, 29 September 2026: corrections adopted

A second independent review of the branch that added P2-25 and P2-26 (commit d7e9b9e) asked that they not be presented
as a completed proof. Adopted:

- **The finite change.** KL = a^2 F / 2 is local. The bound now controls the finite change from no chamber to the
  chamber rigorously within the model: whitening the covariance change into E, KL <= |E|_F^2 / (2 (1 - |E|_F)) with
  |E|_F <= sqrt(F) + eps, eps an explicit power-series remainder (`sarsim.information.remainder`), checked against the
  exact KL on small images (P2-25).
- **Probability terms.** The bound caps detection rate minus false-alarm rate (the total variation, TV); with two equally
  likely cases the best accuracy exceeds a coin toss by TV / 2. The site had called TV "beating a coin toss by"; corrected
  everywhere.
- **A bound, not an attained value.** The Cramer-Rao check is an estimation check; it does not show that any test reaches
  Pinsker's ceiling. Detection is now measured on its own: the exact likelihood-ratio test and the score test, with
  their detection curves, beside Pinsker's bound from the exact KL. "Reached by the best detector" is withdrawn.
- **Extrapolation labelled.** "1.3 billion times stronger shaking" was the linearised information extrapolated, within
  the model, to a reliable detection. The run labels it so; the site replaces it with the requirement curve (P2-27),
  computed exactly within the model.
- **A conditional requirement, not a universal exclusion.** P2-27 states its acquisition, scattering model, footprint and
  target; it covers mechanisms that act through the surface's motion in one image and says nothing of routes outside the
  model or of mechanisms not specified. The phase-only argument covers motion; changes of reflectivity during the pass
  are outside it.
- **The Doppler-to-time relation** the bound rests on is checked against pulse-by-pulse physics (P2-25): where the two
  differ the model holds more information, by an edge effect that fades as the time-bandwidth product grows.
- **Regional, not local.** The ambient levels are regional (Kottamya, 67 km east) and are not upper limits at the
  pyramids.
- **Resonance.** P2-26's "Q" was the width of the largest response's envelope, not a modal Q; the S-wave case's band was
  cut by the model's upper limit, so its 3.5 was no measurement at all; and records of 0.3-0.5 s cannot exclude a
  long-lived narrow mode (a 60 Hz mode with Q 600 decays by 1/e in 3.2 s). P2-26 now reports the envelope's width as
  such, rings each room on its ceiling for 4 s beside intact rock, resolves modal Q from the spectrum's peaks at 0.25 Hz,
  and follows a passing wave for 4 s; its statement is limited to the rooms, excitations and windows modelled.
- **Shortcuts dropped.** "At strains of 10^-11 nothing else is conceivable" (rock nonlinearity) and the air-impedance
  argument against an air resonance (a quality factor near 190) are withdrawn; the air inside a room and the rock's
  nonlinear response are listed as not modelled.
- **Measured validation has not been done.** P2-22 and P2-23 compare geophones and the satellite in simulation. A test
  against ground motion measured on the plateau during a pass, or shaker and reflector data with instruments on them,
  remains a separate, unfinished task.
- **Next, in the order the review gave:** two depths and hollow against solid (P2-28), then a frozen, blinded calibration
  test (`BLIND_TEST.md`, written before any product for it exists).

### Third review, 30 September 2026: corrections adopted

A third independent review of the branch (commit 579754e) re-ran the finite-proof tests, built its own numerical
checks, and found that some conclusions went beyond what the calculations establish. Adopted:

- **The requirement curve is not universal.** P2-27 computed detectability for one smooth pattern (a Gaussian envelope
  oscillating at a stated frequency) and then said any mechanism must supply motion of that size. What an image holds
  depends on the spatial pattern, not only on amplitude, frequency and area. The review's example, reproduced with the
  lab's own `fisher_grid`: a Gaussian 4 m wide at 0.2 Hz multiplied by a cosine of 4 m period is nowhere larger and holds
  39 times the information (489 times for a 1 m period); in a small phase-only model with the exact finite divergence the
  modulated pattern, with half the energy, is 16 times more distinguishable (`tests/test_finite.py`). P2-27 now labels
  its curves as that family's and states the requirement for every pattern from Theorem C, on the motion's energy
  (item 19). The quiet-ground conclusion survives with a smaller margin (4.6e5 in amplitude instead of 8.3e8); under a
  truck close by, with the scattered wave carried unattenuated across the scene, the margin shrinks to a few times for the
  bench room and vanishes under a thin roof, so there the conclusion rests on how far the scattered wave carries.
- **A large upper bound is not a working regime.** P2-29 called an oracle value of 0.99 "enough" for 95% at 5%. Its
  applications use energy bounds and averages over phases, so the values are upper bounds (now `tv_upper`). A small valid
  bound shows every permitted detector performs poorly; a large one shows only that the argument does not exclude
  detection; an implemented detector succeeding on independent tests would be the evidence that a task is achievable.
  The strong-shaking cases are unresolved, not demonstrated opportunities.
- **The angular guard.** P2-32 sampled wave directions every 15 degrees and enlarged the peak by 1.01; a function
  a + b cos 2 phi + c sin 2 phi sampled so can peak sec(15 deg) = 1.0353 above its samples, now used. (P2-25 and P2-30
  sample every 7.5 degrees, where sec(7.5 deg) = 1.0086 < 1.01, unchanged.)
- **The reference covariance.** P2-32 bounded two layouts through the difference of their displacement fields with the
  empty ground's covariance. In the speckle model the same added phase is more or less distinguishable depending on the
  motion already present, so subtraction does not remove the reference. The tight bound and the certificate now take each
  pair through empty ground, TV(A, B) <= TV(A, none) + TV(none, B); the oracle keeps the difference, where, told the
  reflectivity and the background, the two worlds' mean echoes differ by exp(i k u_A) - exp(i k u_B) whatever the
  background.
- **The far field by direction.** P2-25 (and P2-29, P2-30 and P2-32 after it) averaged the scattered wave's far-field
  energy over directions and combined it with the worst direction's near field. Each direction is now paired with its own
  far field (the truck's bound for the bench room moves from 0.022 to 0.023; P2-30's certificate for the thin-roof room,
  which had bounded the far field by four times the average, tightens from 89.1% to 27.3% at 5% with the reference as the
  identity, and stays unexcluded with the background in it).
- **Phase modulation.** A cavity's influence need not appear as a new frequency, and a change of phase can carry
  information; the linear solver's complex response already includes amplitude and phase. A fixed-amplitude phase change
  moves the ground by at most twice its amplitude, but turning that into an exclusion needs a justified extent: the
  affected area cannot be assumed to be the room's footprint. P2-27 now gives the area such a mechanism would have to
  change to escape the bound (5.5e4 m^2 under the regional microseisms).
- **The defensible conclusion** is that useful recovery is excluded under explicitly bounded conditions, not that
  vibration-based underground mapping from one SAR acquisition can never work.
- **Next, in the order the review gave:** the corrections above; then quantify the physical assumptions (excitation,
  material properties, geometry, coupling, radar scattering, the far field's attenuation, omitted spatial tails and
  numerical error) as declared ranges; then a shape proof organised around two meaningfully different layouts,
  distinguishable after the allowed depth ambiguity, with a conservative bound on their differing echoes throughout the
  declared range, so that every algorithm receiving less information inherits it. One hard pair rules out uniformly
  reliable mapping across a class, not detection of every member.

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

The same tracking was first applied by the method's authors to bridges (Biondi, Addabbo, Ullo, Clemente and Orlando,
2020, *Remote Sensing* 12, 3852) and to the Mosul Dam (*IEEE JSTARS* 13, 6337), without synchronous ground truth; the
2026 assessment notes that gap and fills it for bright targets.

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
- The bound of 17 holds exactly for fully developed speckle and for bright points taken one at a time; real texture
  holds scatterers in between, whose information lies between the two bounds, and the genie bound, which assumes the
  reflectivity known, covers any scene. It is computed noise-free and with the method told the ambient motion, both
  generous. The imprint is P2-04's static one below 8 Hz and P2-26's dynamic one above, for one chamber in uniform bare
  rock without attenuation; the scattered wave is carried unattenuated across the whole scene by a bound.
- Giza's own ambient level is unmeasured in the open literature; the microseism level is from
  Kottamya, 67 km east (M1-01), bracketed by Peterson's global models. None of these levels is an upper limit at the
  pyramids, whose traffic and visitors add their own shaking.
- Nothing here has been validated against measurement taken together: no ground motion recorded on the plateau during a
  pass, no shaker or corner-reflector data with instruments on them read alongside the radar. P2-22 and P2-23 compare
  geophones and the satellite in simulation. This is a separate, unfinished task.
- P2-27's requirement is conditional: on the acquisition, fully developed speckle without a reference, one coherent
  pattern at its worst phase, and the stated target. It does not exclude mechanisms that are not specified, nor routes
  outside the model (reflectivity changing during the pass, the echo from beneath at longer wavelengths, texture between
  speckle and bright points, many images).
- P2-26 covers the rooms, excitations and observation windows it models: fractured, layered or coupled structures, the
  air inside a room and the rock's nonlinear response are not modelled.
- The known-chamber test at Khufu on two real passes (2022 and 2025) is run locally and not yet published. Stated in
  advance: if a frozen method scores the surveyed chambers above matched controls, blind, on both passes, that would
  weaken 1 and the scope of 4; if it does not, it strengthens them. P2-16 makes the same comparison, in the open, for
  the gated reconstruction across the whole pyramid.
- P2-13's loops are the lab's, on Khafre's faces; P2-14 runs the implementation's own base lines. The motionless twin
  keeps the image's smoothed brightness and mean spectra, not its fine texture, and registers a little more steadily
  (41% of shifts exactly zero against 36%), which is why open desert in the same image is run beside it. The
  implementation's eleven Great Pyramid lines use a newer format the published pipeline does not read and are not
  run; they lie on the east face, nearly edge-on to this pass (86 degrees) and as dark as the darkest tenth of the
  scene around the pyramid. Whether their pixels are mostly receiver noise would need a ray check against the terrain;
  darkness alone does not show it.
- P2-11's penetration bounds use dry sand's measured penetration as an upper bound for limestone, and a limestone
  conductivity from the standard GPR table; no microwave loss measurement of Mokattam limestone itself is used.
