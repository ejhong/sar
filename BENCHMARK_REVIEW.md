# Review of the benchmark (P2-36): the certificate holds; two things to fix before it is cited

30 September 2026. Reviewed at `acf24fc` on `claude/upbeat-pascal-swk910`: `BENCHMARK.md`, the generated
`sim/results/p2_36_benchmark/report.md` and `summary.json`, `sarsim.finite.per_line_certificate` and `line_floor` with
their tests, `sim/experiments/p2_36_benchmark.py`, and the audit's fifth-review entry. On the branch `uv run pytest`
passes (166 tests). Beyond the branch's own checks I ran the certificate on 24 random lines with common phases up to
1.2 rad and differences up to 0.5 rad: the floor never exceeded the covariance's smallest eigenvalue and the exact
divergence never exceeded the certificate; where √Q_c ≥ 1 or ρ ≥ 1 the certificate is silent, as designed. The
fast-forward of `main` is clean (52 commits, `origin/main` an ancestor of the branch; `docs/` built with the benchmark).

**Verdict.** The per-line certificate is valid, and it settles the fourth review's first point for the single image:
each line's reference now carries that line's own share of the lorry's wave, world 0's own cavity and the microseisms'
envelope, with floors of 0.95 beside the lorry and 0.99 elsewhere, where the whole-image floor failed. The benchmark
has the form both reviewers asked for. Two things must change before it is cited as closing the strong case, and one
sentence needs a different argument.

## 1. What holds, and why

- **The algebra.** With C₀ = B₀B₀* + σ²I and B₀ = E∘e^{−iφ₀}, ‖C₀^{−1/2}B₀‖ ≤ 1, so the whitened change has
  Frobenius norm at most 2q + q² with q² = ‖B₁ − B₀‖²_F / λ_min(C₀); λ_min(C₀) ≥ σ² + (1 − ‖B₀ − E‖_F)² by Weyl; and
  KL ≤ ρ²/(2(1 − ρ)) is the tightened constant. Lines are independent exactly once the range band is widened to every
  bin (white scatterers and noise; each range pixel's scatterers reach only their own line), and the real image is a
  function of the widened one, so Σ over lines is an upper bound. Nothing here is approximate.
- **The branch's own checks agree.** On a 2,048-pixel line through the largest difference the exact divergence is 0.12
  to 0.18 of the certificate at every amplification that leaves it informative; pulse by pulse on fixed scenes the
  coherent echo difference is 0.994 ± 0.006 of the ensemble formula, unchanged with the lorry's wave ten times larger;
  the lines beside the lorry have floors of 0.951 and the verification line 0.995.
- **The single-image bound does not depend on the assumed SNR.** Noise enters L1 only through the floor, which is
  already 0.95 to 0.99; the sweep from 20 to 60 dB leaves every L1 value unchanged. The desktop's calibration and
  noise fields therefore change only the oracle rows (the privileged detector told every reflectivity), whose bound
  crosses 0.9 at 33 to 38 dB. The owner's plan should say so: item 2 calibrates the oracle, not the exclusion.

## 2. Not like for like: the case that was open was a 10 m room under a plane wave

The message says "under the lorry, a single image drops from unresolved to excluded". The case that was unresolved
(P2-26's `surface_favourable`, P2-29's open row, P2-30's 27% at 5%) is a **10 m room whose roof is 5 m down**
(`FAVOURABLE = {'centre': [0, 0, -10], 'size': [10, 10, 10]}` in P2-22 and P2-26) under a **plane surface wave at the
FTA level everywhere**. The benchmark's WA is a **6 m room** with a 5 m roof under a **point force 15 m away**, the FTA
level at the site only. Part of the drop is the smaller room and the source's decay, not the certificate: at 69 Hz the
Rayleigh wavelength is about 25 m, so ka is 0.8 for the 6 m room and 1.3 for the 10 m one, and the scattered amplitude
scales between a² and a³, that is 2.8 to 4.6 times more for the 10 m room; its 5 m roof spans 10 m instead of 6 and
flexes more (P2-26 found it moving the ground over it up to twice the passing wave). The L1 margins are 27 (room) and
14.5 (L tunnel) in amplitude, so the 10 m room is probably still excluded, but that is a guess until it is run.

Fix: add a fifth world, WC = the 10 m room with the 5 m roof, and a fifth pair W0 against WC, and let the finding say
"the case the earlier runs left open, run in the benchmark's frame, gives L1 = x". Until then the sentence should read
"for the benchmark's 6 m room". The same wording change applies to the proof page's headline, the roadmap's State list
(item 4: "one image allows at most 0.06 for a room under a 5 m roof") and the audit's item 28.

## 3. The strong-case bound is 97% an assumed far field

Within the measured 70 m the single-image bound is 0.006 to 0.009 and the oracle's 0.09 to 0.14; the tail (an
unattenuated surface wave carried to the image's edge at 2.5 km) supplies 97 to 98% of the phase energy. That is the
conservative direction for an exclusion, and the exclusion survives the tail at ×5 (L1 at most 0.138), so the
conclusion does not rest on the envelope's exact form. Three things to say and one to do:

- Say in the finding that the headline number is 97% envelope, and give the within-70 m number beside it; a reader
  who sees 0.033 and then learns 97% of it is assumed will feel misled if the finding does not say so first.
- The measured σ(r) falls over 40 to 70 m (slopes −0.34 to −0.58; σ(67.5)/σ(42.5) = 0.78 to 0.86), so the flat
  envelope already overstates the tail; and "unattenuated to 2.5 km" at 69 Hz is far outside any real limestone
  (Q = 50 gives an attenuation length near 400 m). The Q = 50 rows (L1 0.011 to 0.019, oracle 0.16 to 0.27) are the
  realistic ones; say which is the headline and which is the worst case.
- The frequency rule picks 69 Hz by the oracle on the disc; L1 rises with frequency to the band's edge, and the
  finding rightly reports the band's worst (0.072). The proof page's stat row uses the frozen frequency
  (`bmMax` takes the maximum over pairs at 69 Hz, 30 dB), so it prints 0.062 where the finding says 0.072; make the
  page use the band's worst.
- To do, later: enlarge P2-26's or the benchmark's box so σ(r) is measured to 100 m or more at the worst frequency; a
  falling σ measured beyond 70 m would replace the envelope with a measurement (fourth review, A3).

## 4. Smaller points

- **The quiet case's ensemble statement.** The certificate is evaluated at the expected Q over the field's phases and
  called "an explicitly ensemble-averaged statement". Jensen runs the wrong way for that: the certificate is convex in
  Q, so certificate(E Q) ≤ E certificate(Q), which does not bound the average over realisations. What makes it right
  is that at these values (Q of order 10⁻¹³ per line) the certificate is linear in Q to relative order q, so the
  expectation passes through it with negligible error; say that, and keep the worst realisation's √(2n) beside it as
  the exact statement. (The oracle's ensemble step is fine: its TV is concave in Δ².)
- **`FROZEN['r_obs_m'] = 39` and the docstring** say the observation is a 39 m disc with the tail bounded separately;
  `BENCHMARK.md` and the report say the whole image is observed, with 39 m and 70 m as reported sub-discs. Make the
  docstring match the document.
- **States.** 0.033 prints "near chance" and 0.062 "95% at 5% excluded" for numbers that differ by a factor of two;
  the rule is declared and the report prints the number beside the word, which is right; keep it that way on the
  pages too (found ≤ false alarms + 3.3%).
- **Location.** Le Cam's floor of 3 m × (1 − TV) on the 6 m pair is stated correctly; it says a reader cannot tell
  which of two places 6 m apart the room is in, and nothing finer.
- **P2-31's rerun** in double precision gives the published blind map 4 of 9 (p = 0.37) at 20 rad: honest, and it
  should replace "does neither" on the pages wherever that survives.
- **D3** (the published method on the benchmark's strong-case images, with a positive control) is the one layer still
  missing from the table; the message's item 4 covers it.

## 5. For the desktop's two exports

The oracle needs Σ_j SNR_j over the scatterers inside the claim's footprint, so the power map should be |z|² **summed**
into cells (5 m is fine; the claim's imprint varies over 100 m), not averaged, over at least 2 km around Khafre for
each pass, with the calibration constant that turns |z|² into σ⁰, the NESZ per cell (or the product's noise-equivalent
field), and the cell grid's georeference. With that, the Khafre claim's oracle is a number from the real image, and the
noise calibration replaces the assumed 30 dB in every oracle row of P2-29, P2-32 and P2-36 at once.

## 6. Order

Publish (the fast-forward is clean). Then the wording of §2 and §3 (an hour), the fifth world (one solver run, about
the cost of the other four), the page's band-worst headline, and the docstring; then the desktop's exports feed the
oracle and the claim; the published method on the strong-case images can run meanwhile.
