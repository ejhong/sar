# Fourth review: what the proof still needs to be airtight, and what the site needs

Updated 30 September 2026 after the external reviewer's six comments on `439c48b` (R1 to R6 below) and a code-level
verification of those comments and of this review's own first version: eleven claims, each read against the code by
one reader and then by a second reader told to refute the first, with recomputation wherever a number was involved.
Where the first version of this review was wrong it says so. Branch `claude/upbeat-pascal-swk910` at `439c48b` against
`main` at `8cc150b`; on the branch `uv run pytest` (164 tests) and `npm run check && npm test && npm run build` pass.

**The short version.** The reviewer's six points all hold in the code. The largest is R1: every lorry number outside
the oracle is computed against a reference that leaves out the lorry's own wave, and with that wave in the reference
Theorem C proves no floor under P2-30's own convention, so those rows are uncertified either way. The fix is one exact
computation per world (A1). Two things the verification found that neither the reviewer nor this review had seen:
P2-31 synthesises its images in single precision while its real-level imprints are 20 to 40 times below single
precision's rounding, so its "real level" row measures rounding, not the imprint (A7); and 99% of each lorry energy in
P2-27 is a far-field extrapolation that no grid check covers (A3). Two things in the first version were wrong: the
cube's cos 4φ guard (the static kernel is exactly a quadratic form in the wave's direction, so the sec δ guards are
exact for the static rows, A9) and the background's Fisher information (1.5 × 10⁻⁹ from the slow limit, now 3.0 × 10⁻⁹
with the pass's own sampling, A1). The rest stands: the Khafre claim is uncomputed (A2), the envelope needs a declared
product (A8), the overview's step 10 is an essay (B2), and the branch must merge `main` before anything else (B1). The
quiet-ground conclusion survives every correction here with margins of 10⁵ and more; the lorry rows, whose margins are
2 to 4, do not survive them unchanged.

## A. The proof

### A1. One reference per world, computed exactly (R1)

What the code does, verified. P2-30 builds its floor of 0.0559 from the three regional harmonics alone (lines 88 to
97) and reuses it for every dynamic row's `background_*` fields (168 to 174) and, through P2-27's `any_pattern`, for
E_min (9.75 × 10⁻¹⁰ m⁴ = 1.74 × 10⁻⁸ × 0.0559). P2-25's dynamic rows whiten by C₀ = I (`fisher_grid`, `ambient_kl`);
the only justification in its text, "the method is told the ambient motion", tells the mathematics nothing. P2-26
defines the imprint as with-minus-without the room against the incident motion, so in the lorry world the no-room
reference carries the lorry's incident wave, whose line-of-sight envelope is 7.3 times the regional one at 52 Hz.
Under P2-30's own whole-image convention (the incident wave at its 15 m amplitude on every pixel) √Q_bg is 2.7 to 36
for the truck and urban rows, so `baseline_floor_from_phase_energy` returns no floor: Theorem C proves nothing in
those worlds under that convention, and the regional floor does not bound them. With the incident wave confined to a
declared extent a floor reappears (the truck at full amplitude over at most 4.5 × 10⁴ m² at 52 Hz with the regional
harmonics also in the reference, 6.7 × 10³ m² at 20 Hz, 1.1 × 10³ m² at 8 Hz); P2-30 declares no extent or
attenuation, so its dynamic `background_*` rows are unsupported either way. The identity-reference lorry values (thin
roof 27.3% at 5%) are conditional on C₀ = I in a world whose reference is demonstrably not I under that convention, so
they are no better certified than the floor-based ones. Only the oracle's lorry rows are clean: its Δ² cancels the
background exactly, scatterer by scatterer, checked pulse by pulse (0.993 ± 0.025, and 1.011 ± 0.030 with the
background ten times larger). P2-32's lorry rows are oracle-only and need no change.

The reviewer's two alternatives. An invariance that removes the common motion exists for motion shared by the whole
scene: a phase common to every scatterer at each pulse is a diagonal unitary on the raw rows, so C₀ = I exactly (in
the band-limited model as `fisher_grid` references each row to its own zero-Doppler time, nearly so: F ≈ 2 × 10⁻¹⁰
over the image; whether that residual is real or a convention of the model is unresolved between the two readers, and
either way the floor is nearly one). The lorry's wave at 8 to 100 Hz has a Rayleigh wavelength of 17 to 210 m across a
5 km scene, so it is not shared and no invariance removes it; the same holds, less severely, for the 1 to 3 Hz and 3
to 8 Hz rows (wavelengths 845 m and 338 m), which P2-30's background rows cover and the audit should not call shared.
A calibrated receiver-noise floor would bound λ_min(C₀) below for every excitation, but the certificate and tight
layers carry no receiver noise (information.py's docstring: none, "which could only lower what follows") and only the
oracle has it (σ²I at 30 dB per cell). In the model's normalisation (texture power 1 per cell) noise at 30 dB gives a
floor of about 10⁻³, 56 times below the regional 0.0559; E_min scales with the floor and the margins with its square
root, so a noise floor would leave the static rows excluded by about 6 × 10⁴ and certify nothing for the lorry rows.
Say that with the number, or compute it.

Relabel now, before the computation. (i) In P2-30 drop or mark the dynamic rows' `background_*` fields: "with the
lorry's own wave in the reference no floor is proved"; keep the identity values labelled "conditional on C₀ = I, which
holds for motion shared by the whole scene but not for a wave from a source 15 m away"; a corrected background row
needs a declared extent or attenuation for the incident wave. (ii) In P2-27 report the lorry comparisons against
`E_min_identity`, conditional on C₀ = I (bench 10.5 to 15.8, thin roof 1.9 to 4.1; near-only 141 to 150 and 22 to 34),
and say the floor-based E_min applies to the regional-background world only. (iii) Propagate to the P2-27 finding, the
proof, satellite and overview pages and the audit, all of which quote the floor-based lorry numbers.

The computation. Two routes; do the first, keep the second as the analytic backstop.

1. **Exact, per line, with each world's own reference in it.** P2-27's `line_kl` computes KL exactly for one line
   along track with the full range band. In the model the lines are independent exactly (white scatterers; with the
   whole range band each range pixel's scatterers reach only their own line), and the narrower image is a function of
   the wider (Theorem A), so Σ over lines bounds the SLC's KL with no Fisher information, no remainder and no floor.
   Generalise it from the Gaussian family to (a) the imprint's actual map on the line, direction by direction, and (b)
   the world's common phase ψ(x, t) in both covariances, KL = tr(C₀⁻¹C₁) − n_b − log det(C₀⁻¹C₁). Lines the imprint
   does not touch contribute zero. With ψ present the window along track is a truncation: show the sum is insensitive
   to the window's length. Run it for the regional microseisms, for the lorry's direct wave and for the urban level
   (Q_bg 10.5 to 48 under the same convention: no floor there either). It replaces "reference covariance I" everywhere
   with one exact number per world, and the stat row's four numbers with one.
2. **The lemma, with corrected numbers.** λ_min(C_bg) ≥ 1 − √F_bg − ε_bg holds (information.py's derivation is for any
   phase field; checked exactly on 64 × 8 images), and carrying the background into the imprint's bound needs one more
   term: ‖E′‖_F ≤ (√F_imp + ε_imp + (e^{2P_bg} − 1)(√(2F_bound,imp) + ε_imp)) / λ_min(C_bg). The first version of this
   review quoted F_bg = 1.5 × 10⁻⁹ from the slow-limit formula; that omits the pass's own along-track sampling, which
   adds f/V = 2.85 × 10⁻⁵ cycles/m to the wave's spatial frequency. For a 15 km wave at 0.2 Hz the interior value is
   3.0 × 10⁻⁹ travelling one way along track, 4.8 × 10⁻¹⁰ the other way, 2.7 × 10⁻¹⁰ across track (3,535 m disc). Two
   conventions must be kept apart: the model's own covariance over the image's pixels, for a field that fills the
   image, is dominated by the image's edge or periodic wrap (1.6 × 10⁻⁵ exact over the full 5 km strip; hard-edged
   discs on coarse grids give grid-dependent 4 to 8 × 10⁻⁷), an artefact of cutting a field that continues past the
   image, but it is what the model gives; the rigorous route to the interior value treats the image as a window on a
   periodic box that contains it (105 km × 5 km), whose in-band covariance bounds the window's from below, with F_box
   21 times the interior value. Floors, at P2-25's rms amplitude: regional microseisms 0.9997 (box, ε_bg 6.0 × 10⁻⁵)
   or 0.996 (wrap, ε_bg 1.3 × 10⁻⁵), and 0.9932 (wrap) at P2-30's line-of-sight envelope amplitude, 1.7 times larger;
   the noisiest stations 0.68 (box; ε_bg = 0.31 dominates) or 0.65 (wrap); the lorry's direct wave spreading as 1/√r
   from 15 m to the scene's corner at 3,535 m, `fisher_bound` 0.572 and ε = 3 × 10⁻⁵, floor 0.24 (0.44 with the
   travelling-wave F of 0.315). The cross term is negligible for the microseisms (7.6 × 10⁻¹² against √F_imp = 2.6 ×
   10⁻⁹) and a 21% addition at the noisiest stations. Whether a spatially uniform oscillation counts as shared motion
   in the model is unresolved between the two readers: the exact band-operator argument gives C₀ = I for motion shared
   per band bin, while `fisher_grid`, which references each row to its own zero-Doppler time, gives a constant K the
   spatial frequency f/V and F ≈ 1.7 × 10⁻¹⁰ over the image; either way the floor is nearly one, and the audit's
   "shared" should say which it means.

To compute, not assume: E_min for the regional world rises from 9.8 × 10⁻¹⁰ toward 1.7 × 10⁻⁸ m⁴ and the static margin
from 4.6 × 10⁵ toward 1.9 × 10⁶; for the lorry world the lemma's floor is 0.24 to 0.44 and route 1 does better;
whether the thin-roof rows (0.45, 0.97) cross 1 is unknown until it is run.

### A2. The Khafre claim is the one case not computed, and a reader goes there first

P2-25's `claimed` block is "an estimate by scaling, not computed" (relative imprint 0.36 of the bench room's). The
site's question is the claim. Compute it: P2-04's kernel on `bench-khafre-claim` (composed, `listed: false`) or
`giza-deep`, then the tight layer and the oracle for the 80 m cube at 1,220 m, whose imprint spreads over a kilometre.
For the oracle use the real 2025 image's pixel powers over that footprint, Σ_j SNR_j = Σ|z_j|² / NESZ per cell, the
pyramid's bright faces included: a real-data anchor that needs no texture model, and the strongest single number the
site could hold. The product lives only on the desktop; a downsampled power map (a few hundred kB) exported into
`sites/acquisitions/` makes the run reproducible in the cloud. Do the 2022 pass as well.

### A3. The displacement integral needs a declared error budget: tail, taper, far field, directions, grid (R2)

Verified, with numbers, all on the worst direction of P2-27's `imprint_energy` unless said otherwise.

- **The static tail.** The code sums u² over P2-04's 41 × 41 grid at 2 m (r_max 57 m) and adds a 1/r² tail from R_edge
  = 55 m whose constant comes from 21 points that all sit at the four corner azimuths; the annulus 40 < r < 55 m
  outside the square (69% of its area) is counted nowhere (4.5% of the total). The grid's fully sampled annuli decay
  as r^−1.3 to r^−1.7 over 22 to 38 m: pre-asymptotic. The second reader re-ran P2-04's kernel on a 100 m domain with
  receivers to ±84 m (three runs of about 200 s): the vertical component is r^−2.8 to r^−3.3 beyond 40 m and the
  horizontal r^−1.4 to r^−2.0 over 40 to 84 m, so 1/r² is the right line-of-sight asymptote; the true energy beyond 55
  m is 1.23 times the coded tail, the full worst-direction energy 1.064 times the coded, and the static margin moves
  from 4.57 × 10⁵ to 4.43 × 10⁵ (the 100 m domain reproduces the 60 m domain's kernels on every ring to 0.05%). Fix:
  extend P2-04's receiver grid to ±84 m and replace the assumption with a measured disc plus a tail under 5%; or state
  the tail as an assumed asymptote (p = 2, the slowest force-free static decay), take its constant from the 36 to 40 m
  ring's azimuthal mean, start it at 40 m, and report the fitted exponent beside it; or take the tail from the
  analytic point-moment field already in the codebase (`katabasis.seismic.analytic.moment_surface_displacement` with
  `void_moment`, the cavity's moment tensor in a half-space, u_z ∝ cos 2φ / r² because M_xx ≠ M_yy, whose azimuthal
  nodes are why the four corner azimuths under-read the 55 m ring's mean u² by 0.79 to 0.84), which reproduces the
  solver's line-of-sight field beyond 20 m with correlation 0.999 (scale 1.36 at 20 to 40 m falling to 1.22 at 70 to
  84 m) and, scaled at the 36 to 40 m ring, gives 2.17 beyond 55 m and 0.93 for the uncovered 40 to 55 m region, both
  above the measured 1.97 and 0.87, at no computing cost. That is the envelope with its coefficient and angular
  dependence the reviewer asked for: p = 2, the coefficient from the ring, the cos 2φ form stated. Replace "bounded as
  1/r²" and "(generous)" with that.
- **The taper.** `energy_near` weights |H|² by tp², full to 32 m and zero at 39 m (∫tp² over 32 to 39 m is 7/3 of 7
  m). The untapered energy within 39 m is at most 1.07 to 1.20 times the coded value (estimated with the far-field σ
  for the annulus; the exact value needs `maps.npz`), so the near-only margins move from 35.4 and 33.3 (bench) and
  8.05 and 5.17 (thin roof) to 33.8, 30.3, 7.77 and 4.82. The same quantity carries three labels: "within 39 m
  exactly" (P2-25, P2-27), "counted within 39 m only" (proof page, audit) and "within 32 m" (satellite page, overview,
  P2-29's "within 32 m alone"). One phrase, "the solver's field at full weight to 32 m, tapered to zero at 39 m", or
  compute it untapered (`maps.npz`, four runs of about 270 s).
- **The far field, 99% of every lorry energy.** In P2-25 the annulus 32 to 39 m is counted twice (tapered in the near
  field, at full σ in the far bound), an over-count of 0.13% of the far term, conservative, and, given constant σ
  beyond 40 m, the totals under-count nothing. But 98.6 to 99.4% of each dynamic energy in P2-27, and nearly all of
  the dynamic Fisher information, is the far-field bound: the cross-width σ (energy per metre of radius, the largest
  over the 25 to 40 m annuli) carried unattenuated to the scene's corner. That σ never exceeds its 25 to 40 m maximum
  beyond 40 m is the assumption the word "unattenuated" stands for; a Rayleigh part keeps constant energy per metre of
  radius while body-wave parts fall, so it is plausible and untested in the code. Test it on P2-26's own model: raise
  its recorded half-width from 40 m (to 60 m within the present box, whose absorbing layer starts at 90 × 60 m, or
  enlarge the box to reach 84 m) and measure σ on the annuli beyond 40 m at the worst frequencies (four cases of about
  270 s at the present size, more when enlarged); name it as the load-bearing assumption of every lorry row.
- **Directions.** P2-27 takes the largest of the sampled directions with no guard; P2-25's docstrings say "every
  direction", "whichever side" and "worst side", and the proof and satellite pages "each wave direction with its own",
  where the dynamic maps are four quarter-turn copies of one west-wave run (the diagonal incidence is not simulated;
  the site is symmetric, the model box, 220 × 160 × 103 m with the line source at x = −75 m, is not); P2-30 alone says
  "no unsampled-direction claim". With four samples no guard is possible: any angular harmonic of order two or more
  can vanish at all four. In the long-wave limit worst-of-four under-reads the continuous worst direction by about 10%
  in amplitude (the static analogue gives 1.198 in energy). Either simulate 8 or 12 incidence directions at the worst
  frequency (guard sec(mh/2) for angular content of order m; ka is 0.2 to 1.3 for the bench room over 20 to 120 Hz) or
  label the dynamic energies worst-of-four everywhere.
- **The grid.** P2-26's fine grid gives smaller near-field imprints (amplitude 0.94 to 0.99 for the bench, 0.90 to
  0.98 for the thin roof), so the 1 m grid overestimates the near field, conservative in sign; but the check covers
  the near field under a surface wave only (nothing for P or S from below) and not the far-field σ. P2-30's uniform
  11% phase stress is labelled a stress, not an error bar, correctly. Fix: declare one allowance per row with the
  check it rests on. Static rows: P2-04's 1 m against 0.5 m grid (the epicentre 18.04 against 18.67 µm, 3.5%; the
  uplift ring 3.37 against 3.04 µm, 10%) and the domain check (100 m against 60 m, 0.05% rms on every ring). Dynamic
  near field: P2-26's 0.5 m check, surface incidence on the small model only; run it for P and S from below. Far-field
  σ: no check exists; measure it at two grids on the enlarged P2-26 model above. Until each is run, P2-30's 11% is a
  stress, not the allowance the reviewer asked for.

None of this moves the quiet-ground conclusion (4.6 × 10⁵ becomes 4.4 × 10⁵). All of it matters for the lorry rows:
the bench room's any-pattern margins are 3.7 and 2.5, the thin roof's 0.97 and 0.45 (the 20 Hz row within 3% of 1).

### A4. Three states, per layer, and headlines that name the 5 m roof (R3)

Verified: "found ≤ false alarms + TV" is stated once in prose on the proof page and nowhere as a table state. The
boundary table's lorry row folds the bench room under the lorry, which the oracle, the one lorry layer A1 leaves
standing, excludes at 95/5 (found at most 46% at 5%, accuracy at most 70.6%; the certificate's 23% with the
background, 8% with the identity and the tight layer's 0.023 are A1's uncertified numbers and stay labelled
conditional), into "not excluded, unresolved", driven by the thin roof (oracle 0.99); its shape cell says "not
excluded: 0.82" where 0.82 < 0.9 and P2-32's own finding says still excluded until the signal grows 1.2 times. The
proof page's headline and the satellite page's say "strong shaking close by, known exactly, is not excluded" without
naming the thin roof; the overview names it. The satellite page's "open" is a fourth notion (growth after allowances ≤
1) that attaches to rows as low as TV 0.016. In the runs `excludes_095_at_005` is a bare TV < 0.9 test, so TV 0.3 is
flagged True. P2-32 also applies Giza's floor to the noisiest-station rows where P2-30 gives no floor (certificate
cells of 6.1 to 13.8% at 5% for those pairs): drop them or mark them "no floor at this level", and compute the floor
per band as P2-30 does; the page's 5.110% reads only the microseism row and stands. One candidate for the middle state
does not work: the noisiest ground is near chance on every layer (certificate with the identity 2.9 × 10⁻³, tight 9.3
× 10⁻⁸, oracle 1.2 × 10⁻⁴; P2-30 gives that level no background-in certificate, √Q_bg ≈ 55). The example that does
work is the bench room under the lorry.

Fix: declare one rule and apply it per layer (suggested: near chance if TV ≤ 0.01, that is found ≤ false alarms + 1%
and accuracy ≤ 50.5%; "95% at 5% excluded, not near chance" if 0.01 < TV < 0.9; unresolved if TV ≥ 0.9 or the bound is
uninformative, where the bound is silent and shows neither detectability nor usefulness), a row's state from the
smallest valid bound (for the lorry rows, after A1, the oracle's until route 1 has run), every cell printing "found ≤
false alarms + X". Split the lorry row into the bench room (excluded, not near chance) and the thin roof (unresolved);
fix the shape cell; add a state column to the three-layers table; add a near-chance branch (v ≤ 0.01) to the
calculator's three and reword its v < 0.9 branch to "the 95% at 5% target is excluded; a detector could still find up
to X% at 5%"; reword P2-30's finding "still excluded if its phase is 11% larger" to "still excludes 95% at 5% (31.3%
at most)"; name the 5 m roof in every headline that says "not excluded"; rename the satellite's "open" to "not
excluded after allowances"; add a `state` field to the runs beside the boolean, and to the `background_included` rows,
which have no boolean today. For the owner's goal, say plainly which of presence, location and shape each state
applies to: one hard pair rules out uniformly reliable mapping over a class, not the detection of every member.

### A5. Every realisation, every case, and how many images

- **The ambient field's phases.** The static rows are the expectation over the field's phases at a fixed even split of
  the power over 24 directions (12 in P2-32); `ambient_kl` already carries the worst-phase bound in its denominator.
  For every realisation and every split the triangle inequality gives F ≤ Σ_d (A_d + |B_d|) ≤ 2n × mean, so the TV
  ceiling rises by at most √(2n) = 6.9 (4.9 in practice: 1.3 × 10⁻⁹ becomes 6.3 × 10⁻⁹). The page shows the
  expectation, and the proof page's "static imprint, every direction" reads as worst-over-directions. The worst-phase
  case is n aligned phasors, a local level √n louder with probability about e^−n, and the cap grows with finer
  direction sampling, so the expectation is the right ceiling for a measured rms level: say so, and give √(2n) as the
  cost of "every realisation".
- **The scene's phases.** The oracle's Jensen step bounds the average over the scatterers' phases; the plateau is one
  realisation. State a concentration: the raw check's spread (±7% at 802 scatterers) scales as 1/√N; or the Markov
  form: with probability at least 1 − 1/t over the phases Δ² ≤ t E Δ², so at t = 10⁶ the microseism row is still ≤ 1.6
  × 10⁻³.
- **Composite alternatives and several images.** TV(mixture, P₀) ≤ max over members, so the largest-imprint member
  covers every prior over rooms and layouts. N independent images: KL adds, TV ≤ √N times the one-image bound; state
  how many passes the 2025 announcement used.

### A6. Wording that overstates, with the replacements (R4, R5, R6c)

- **"Energy" and "rms".** E = ∫ max_t |u_los|² dA is m⁴, a sum of squares inherited from `finite.py`'s dimensionless
  phase energy, not a mechanical energy; and "an rms of 1.0 µm over 1,000 m²" is √(E_min / A), the area-rms of the
  temporal peak, that is the peak of a motion filling 1,000 m² evenly (the time-rms of such a sinusoid is 0.70 µm;
  P2-30 names the same quantity "peak displacement", correctly). Replace throughout: the P2-27 finding and docstring,
  audit item 19 and the third-review entry, ROADMAP, the proof page's mechanism heading and card, the satellite page,
  and the overview's "carry enough energy" (→ "its peak displacement squared, added up over the ground, must reach a
  stated value"). Keep the key names or rename `energy_m4` to `E_m4`.
- **"4.7 × 10⁻⁴ of what the most informative pattern could".** The denominator is `fisher_bound`, 4N⟨Φ²⟩, an upper
  bound no pattern attains (equality would need every spectral component at the sin² maximum and no out-of-band loss
  at once); the best pattern found inside the same envelope reaches 0.89 of it (a travelling wave of 0.5 m period,
  whose half-period is the 0.25 m paired-echo shift; 0.84 for a sign-modulated one; the band-edge limit is 0.90).
  Replace with "4.7 × 10⁻⁴ of the ceiling 4N⟨Φ²⟩ that bounds every pattern inside its envelope (an upper bound no
  pattern attains; a travelling wave of 0.5 m period inside the same envelope reaches 0.89 of it, the band-edge limit
  0.90, so the Gaussian holds about 5 × 10⁻⁴ of what the best pattern does)". The 39× and 489× examples stand as
  computed.
- **128 m² against 55,000 m².** 128 m² is the count of 2 m cells at or above half the peak of the rms static map (the
  room is 36 m² in plan); the same map is above a tenth of its peak over 3,900 m² and never below 1.7% of it anywhere
  on the 80 m grid, and its extrapolated 1% contour lies at r ≈ 100 to 130 m, 3 to 5 × 10⁴ m², the same order as the
  55,000 m² disc (radius 133 m). The juxtaposition cannot stand in for a stated extent. Drop "the room's own footprint
  is about 128" (hard-coded in the overview) and label 55,000 m² (a disc of radius 133 m) as the least area a
  phase-only change at the full 2a must cover to escape this exclusion: a necessary condition, not a demonstration
  that any mechanism reaches it, under its stated amplitude (2a_los = 1.33 × 10⁻⁷ m from 4.9 × 10⁻⁸ m/s at 0.2 Hz) and
  covariance assumption (the regional floor 0.0559; with the identity reference the same threshold is 992,000 m²,
  radius 562 m, and after A1 it moves with the floor); a change that tapers is judged by its E; the extent is the
  mechanism's to state. Worth adding: a 2a change shaped like the static imprint (effective area 275 to 350 m²) would
  still be 160 to 200 times short in area, so the exclusion survives that shape.
- **"Harder than telling either from none".** Three prose occurrences (overview step 10, the audit's "The organising
  question, widened", BLIND_TEST §4) and the built overview; the code and the proof page use the sum, TV(A, B) ≤ TV(A,
  none) + TV(none, B), at most twice as easy and never harder. It is false in general and the branch's own data show
  it: in P2-32's dynamic layer the L-shaped-against-branching pair (0.82) is easier than branching-against-none
  (0.73), and the proof page prints both. Replace with the sum inequality everywhere and note that for the static
  pairs computed the difference happens to carry 1.2 to 4.1 times less signal. The closest-pair argument in P2-32's
  docstring is valid for uniform (worst-case) reliability over the depth allowance, not for an average criterion over
  a depth prior; say which.

### A7. The shape test and the blind test (R6, and a defect the reviewer did not see)

- **Single precision.** `sarsim.synth.synthesize` defaults to complex64 and its own docstring warns that a chamber's
  imprint (10⁻⁹ rad) needs complex128; P2-07 and P2-28 pass complex128, P2-31's three synthesize calls do not. Its
  real-level unit imprints are 1.5 to 3 × 10⁻⁹ rad, 20 to 40 times below single precision's rounding, and in the
  results the real-level cavity images differ from the no-cavity image by the same 10⁻⁵ as no-cavity differs from the
  motionless copy. So "at the real level every reader is at chance" (audit item 24, ROADMAP, the proof page) measures
  rounding, not the imprint. The conclusion is what the bound predicts anyway, but the row is not evidence for it.
  Rerun P2-31 in complex128; any BLIND_TEST product at the real level must be made that way.
- **Independence and p-values.** Nine images per level are three layouts on three shared speckle seeds under one
  microseism realisation, so the three images of a seed are paired, not independent, and the Binomial(9, ⅓) p-values
  are wrong. The exact within-seed label-permutation test (6³ = 216 arrangements) gives: blind at 20 rad, 5 of 9, p =
  0.17 (not 0.14); paired 6 of 9, p = 0.037; reference at 2 rad, 7 of 9, p = 0.028 (not 0.008); reference at 20 rad, 9
  of 9, p = 0.0046 (not 5 × 10⁻⁵); and when a map does not respond 3 of 9 is forced (the null's support is {3}, not a
  chance rate). Report the exact p, say "three seeds, layouts paired within a seed", and let the display print "1/216"
  at its floor. The paired reader's presence AUC is 0.5 by construction (the no-cavity control's paired map is
  identically zero, the motionless copy's above every cavity): drop or flag it. The controls' presence values are
  near-identical within a seed, so the control sample holds three distinct values, not six.
- **Sample size.** With p = 0.17 for 5 of 9, "the published picture, read blind, does neither" is not shown. Either 30
  scenes, each with its own speckle seed and its own microseism realisation (P2-31 draws one, `rng(72)`, for all 42
  images; the kernels are cached, the synthesizer is the cost), or "no better than chance was shown".
- **BLIND_TEST.md.** (a) "its 95% interval by permuting the labels" asks a null distribution to be a confidence
  interval (the permutation null for 12 against 12 spans 0.27 to 0.74 around 0.5 by construction): ask for a bootstrap
  or DeLong interval and a separate permutation p. (b) The boosted set is one level, 2 rad; at 2 rad even the
  reference detector told the shaking cannot locate a layout by scanning (34 m against the controls' 36 m), though its
  scanning presence statistic already separates cavities from controls (AUC 0.96). So "fails on the boosted set: the
  method does not read the room's motion, however strong" must become "at the tested amplifications", the set must
  include a positive-control level (20 rad), and the outcome must be scored per claim, as §4 says. (c) The "harder
  than" sentence (A6). (d) Evaluation units are independent scenes; shared seeds serve paired comparisons and are not
  independent trials. BLIND_TEST's 24 generated scenes are already drawn independently (§2); the change is to P2-31's
  design and to how §4 and §5 count trials.

### A8. The physical envelope as a declared product, with the numerical budget in it

The third review's next step and the reviewer's milestone: one table on the proof page, the margin that survives it
stated as a product.

| factor | value | status |
|---|---|---|
| local level at the plateau over Kottamya | ×10 | assumed (no open measurement) |
| site amplification | ×3 | assumed |
| cover or layering under the room | ×? | run one layered bench at low frequency rather than assume |
| SNR for the oracle | +10 dB (×3.2 in amplitude) | assumed |
| far-field attenuation of the scattered wave | none | generous by construction; σ beyond 40 m untested (A3) |
| static tail beyond the grid | 3% in amplitude (6% in energy) | measured on the 100 m domain (A3) |
| static kernel grid | 3.5% at the peak, 10% on the uplift ring | measured, P2-04's 1 m against 0.5 m (A3) |
| near-field grid | 6 to 11% | measured, sign conservative; no check for P or S from below (A3) |
| wave direction, dynamic rows | worst of four | about 10% in the long-wave limit; simulate more (A3) |
| wave direction, static rows | sec 7.5° | exact (A9) |
| product, in amplitude | about 120 (P2-29's ×95 for the three assumed allowances, times the numerical factors) | |

Against margins of 4.4 × 10⁵ (or about 2 × 10⁶ after A1) the ground's own trembling survives every declared allowance
with three orders of magnitude to spare; under a lorry (2 to 4, near-only tens) it does not. That sentence is the
airtight form of the boundary, and the table is what a proponent must argue with.

### A9. Small rigour items, revised

- **Withdrawn from the first version: the cube's cos 4φ guard.** The static kernel is Kxx sin²φ + Kyy cos²φ + Kxy sin
  φ cos φ, exactly a + b cos 2φ + c sin 2φ (fit residual 3 × 10⁻¹⁶), so sec 7.5° (P2-25; P2-30's 1.01) and sec 15°
  (P2-32) are exact for the static rows, and the 24-direction mean is the exact continuous mean. P2-27's static
  energy, quadratic in u, does carry 4φ terms; its 24-sample maximum is 0.11% below the continuous one, and a sec 15°
  guard or the trigonometric interpolant closes it. The dynamic rows' gap is A3's.
- **The certificate's constant.** `perturbation_certificate` and the page agree (TV ≤ ρ/(2(1 − ρ))) but both are
  looser than the page's own Theorem B step by 1/(1 − ρ) in KL: valid and conservative, 1.45× at the thin-roof row
  (0.223 printed, 0.185 supported). Tightening it moves q_b from 0.282 to 0.342, E_min by ×1.48 and every
  `short_by_any_pattern` by ×1.21, which takes the thin-roof row from 0.97 to 1.17, across 1; it propagates to P2-27,
  P2-32 and the export. Either tighten code, tests, page and calculator together and regenerate, or change the page's
  Theorem B line to the constant the code uses. The certificate uses Pinsker only where the tight layer takes the
  smaller of Pinsker and Bretagnolle–Huber; the page says both.
- `finite.py` cites `PROOF.md`, which does not exist.
- ROADMAP's State item 4 still says "reached by the best detector in simulation", withdrawn by the second review.
- The claims table's shape row could say that the shape test's layouts are single-level and say nothing about relative
  vertical geometry; its "not shown" for that claim stands (P2-28 tested two depths). ## B. The site

1. **Merge `main` first.** The branch is fifteen commits behind `origin/main` (8cc150b: Sacsayhuamán, P2-33 to P2-35);
   merge `origin/main`, not a local `main` that may sit on an unrelated history. `git merge-tree` shows conflicts in
   `METHOD_AUDIT.md`, `web/src/pages/satellite.astro`, `sim/katabasis/export/runs.py` (both sides extend `PUBLISHED`),
   `web/public/data/runs/index.json` and `docs/*`. Regenerate the index with the export (runs only on a cloud machine,
   per the roadmap's note), rebuild `docs/`, run both suites and the stale-export test before pushing.
2. **Step 10 on the overview is an essay.** 1,116 words in one column with the right column empty, where steps 1 to 9
   run 67 to 281 words beside their pictures; the FAQ under "could anything change it" is 911; the overview grew from
   3,381 to 4,954 words; on a phone step 10 is five screens of prose before any figure. Keep the headline, one
   paragraph (what an image is made of; shared motion leaves almost nothing; the cap), the Routes figure, the Ceiling
   chart, one line on what stays open and the link to `/proof/`. Everything else is already on the proof page.
3. **One headline number, one unit, one state.** The proof page's stat row shows 5.017%, 1.3 × 10⁻⁹, 1.6 × 10⁻⁶ and
   27%: two units, three layers, two of them uncertified for the lorry (A1). After A1 lead with the one exact number
   per world; until then lead with the tight number for the microseisms, the oracle for the lorry, and label the rest.
   The overview says one in 770 million and the proof page says 5.017%; a reader cannot tell they are the same claim.
4. **The three-layers and boundary tables** take the state rule of A4; the certificate notes ("background included" on
   three rows, "reference covariance I" on the fourth) become one phrase after A1.
5. **The Routes figure** has a stray "×" glyph beside "Its echo" and arrows that do not meet the cards; check the
   phone stacking.
6. **The satellite page's status line** says "revised after two independent reviews"; the audit lists three.
7. **Precision theatre.** "5.017%", "5.110%", "0.0559", "1.31 billion pixels": round to what the reader uses, and say
   "0.02 above chance" where that is the point.
8. **The requirement chart** has seven legend entries and its 20 and 80 Hz lines overlap; keep 0.2 Hz, 80 Hz and the
   two any-pattern lines, with direct labels, and relabel the y-axis per A6 (peak displacement, not energy).
9. `SectionMap.astro` declares `fmt` and never uses it (the one hint from `astro check`); the overview hard-codes
   "about 128" where the proof page reads it from the run.
10. **The answer section's** 456-word paragraph and the mechanism card's 300 words are prose where the rest of the
    page is dashboard; the five-theorem cards are the right form, copy it.
11. **The audit** takes a "Fourth review, 30 September 2026: corrections adopted" entry as before; item 17's
    "reference covariance I", item 24's real-level sentence and the third-review entry's "energy" change with A1, A7
    and A6.

## C. The milestone and the order

The reviewer's next milestone is the right organising deliverable: one fully specified, auditable comparison, the same
source in both worlds, a valid reference covariance for that world, bounded displacement tails and numerical errors, a
stated scattering and noise model, and bounds for presence, location and shape; parameter ranges only after that
comparison is sound.

The reviewer's priority for the owner's goal is location and shape with depth left uncertain. In this repository those
bounds are P2-32's: the static pairs at most 5.11% at 5% by the certificate with the background in (the regional
floor, valid there), tight pair ceilings 1.4 to 9.1 × 10⁻⁹ (at most 4.5 × 10⁻⁸ at the worst phase, n = 12), and the
lorry pairs by the oracle alone, 0.41 to 0.82 (the hardest, L-shaped against branching, 0.82, needs 1.2 times more
signal to cross 0.9); their depth allowance is the same-depth-pair argument, valid for uniform reliability only (A6).
The states of A4 apply per pair, and the shape test (A7) is the empirical side of the same priority.

On the reviewer's assessment: agreed for the quiet ground and for the room under a 5 m roof (oracle 0.99, unresolved).
"The stronger-shaking cases remain unresolved" is too broad for the bench room under the lorry: the oracle's bound,
the one lorry number A1 leaves standing, is found ≤ false alarms + 0.41 (accuracy ≤ 70.6%), so 95% at 5% is excluded
there and what stays open is usefulness short of it, and the case after the oracle's ×95 allowances; P and S from
below are the same (0.53, 0.60). No general impossibility claim remains on the site (the satellite page withdrew it
and the audit's scope item says so), and the audit already carries "one hard pair rules out uniformly reliable mapping
across a class, not detection of every member": carry that sentence to the proof page's shape section and to the
overview, where the pages still say "harder than".

In this repository the milestone is:

1. Merge `main` (B1), then the relabelling of A1 so nothing uncertified stays labelled certified while the rest runs.
2. A1 route 1: the exact per-line layer for the regional microseisms, the lorry and the urban level, the presence rows
   and P2-32's static pairs (which use the same regional floor and the same C₀ = I tight layer; its lorry pairs are
   oracle-only and stand). One run; it rewrites the stat row, the layers table, P2-27's any-pattern rows and P2-30's
   floor.
3. The scattering and noise model, stated once on the proof page: fully developed speckle (white circular Gaussian
   scatterers on the pixel grid, 28.1 cells per m²) with no receiver noise for the certificate and tight layers; white
   receiver noise at 30 dB per cell above σ⁰ ≈ 0 dB over ICEYE's best Dwell NESZ of −26.7 dB, +10 dB as an allowance,
   for the oracle. The proof page's open-items table has the row and names what is missing, "a qualified observation
   operator and noise floor for the real product": A2's Σ|z_j|²/NESZ over the footprint fills the noise floor; the
   operator (layover, decorrelation, the real focusing) stays declared, not modelled.
4. A3's error budget: the kernel grid to ±84 m, σ measured beyond 40 m on an enlarged P2-26 model, 8 or 12 dynamic
   directions at the worst frequency, the 39 m disc untapered; each row of A8's table filled from a run, not an
   assumption.
5. A7: P2-31 rerun in complex128 with the exact p-values and thirty seeds; BLIND_TEST.md corrected.
6. A2: the claim on the real image, both passes.
7. A4, A5, A6, A9: states, sentences and constants, an hour or two.
8. B2 and B3: cut step 10, one headline number.
9. The rest of B, the audit entry, publish.

## Appendix: numbers reproduced during the verification (a rerun should match them)

- Regional microseism floor as coded: LOS envelope 7.16 × 10⁻⁸ m, phase 2.88 × 10⁻⁵ rad, Q_bg 0.583 over 7.03 × 10⁸
  in-band bins, floor 0.0559. The lorry's incident wave: LOS envelope 5.2 × 10⁻⁷ m at 52 Hz, Q_bg 31 (16 at 72 Hz, 11
  at 88 Hz): no floor. Any positive floor by Theorem C needs an envelope under 9.4 × 10⁻⁸ m.
- F_bg of a 15 km plane wave at 0.2 Hz, regional level, interior: 3.0 × 10⁻⁹ (along track, one way), 4.8 × 10⁻¹⁰ (the
  other), 2.7 × 10⁻¹⁰ (across track), 3,535 m disc; exact over the image's own pixels with wrap 1.6 × 10⁻⁵; ε_bg 1.3 ×
  10⁻⁵ over the image (6.0 × 10⁻⁵ for the box), at P2-25's rms amplitude. Noisiest stations: ×5,175 in F, ε_bg 0.068
  (interior) to 0.31 (box). Lorry direct wave: `fisher_bound` 0.572, travelling-wave F 0.315, ε 3 × 10⁻⁵.
- Static energy, worst direction (52.5°): grid 18.04 + coded tail 1.59 = 19.63 m⁴ per strain²; measured within 84 m
  20.0 (18.04 on the square, 0.87 in the 40 to 55 m region outside it, 1.09 from 55 to 84 m), plus a 1/r² tail beyond
  84 m of 0.85 to 0.89: 20.88 to 20.90 in all; margin 4.57 × 10⁵ → 4.43 × 10⁵. P2-27's 24-sample maximum 19.631
  against the continuous 19.653.
- Taper deficit within 39 m: ×1.10 and ×1.20 (bench, 20 and 80 Hz), ×1.07 and ×1.15 (thin roof).
- Far-field share of P2-27's dynamic energies: 98.6 to 99.4%.
- Exact within-seed permutation p (216 arrangements): 5/9 → 0.167, 6/9 → 0.037, 7/9 → 0.028, 9/9 → 0.0046.
- Certificate constant: q_b 0.2817 (code) against 0.3423 (derivation); thin-roof row 0.223 against 0.185.
- Ambient worst-phase cap: √(2n) = 6.93 for n = 24; the even-split worst phase 4.86.