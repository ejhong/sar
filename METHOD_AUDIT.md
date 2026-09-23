# What is being tested

This notebook evaluates conditional feasibility and failure mechanisms. It does
not reproduce the original Khafre data or claim that a failed reconstruction
settles the existence of underground columns.

Three methods must remain distinguishable:

1. Biondi and coauthors' published descriptions and original implementation.
2. The public **derivative** protocol v1.5, inspected at revision
   `81451263db5a5ab4624b74d64ee65cb1b611b453` of
   [Replication and Verification](https://github.com/BiondiProtocol/Replication-and-Verification-Biondi-Protocol/tree/81451263db5a5ab4624b74d64ee65cb1b611b453).
   Its README acknowledges missing original details and added empirical choices.
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
