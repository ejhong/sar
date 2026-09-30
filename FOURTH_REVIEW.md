# Fourth review: what the proof still needs to be airtight, and what the site needs

Reviewed 30 September 2026: branch `claude/upbeat-pascal-swk910` at `439c48b` (P2-25 to P2-32, `/proof/`, the
overview's step 10, the satellite page's `#bound`) against `main` at `8cc150b`. On the branch `uv run pytest` (164
tests) and `npm run check && npm test && npm run build` pass. `sarsim/information.py` and `sarsim/finite.py` were read
whole; the numbers quoted below that are not on the site were computed with them (the snippet is in the appendix).

**The short version.** The theorems are right and the checks are real. What is not yet airtight is the join between the
layers. The page leads with three numbers for one claim (5.017% found at 5%, 1.3 × 10⁻⁹, 1.6 × 10⁻⁶) because the tight
layer takes the reference covariance as the identity while the certificate proves a floor of 0.056 under it, and the
gap between those two is not physics, it is a loose lemma: the background's own information across the whole scene is
about 10⁻⁹, so the reference is within 4 × 10⁻⁵ of the identity. Close that (A1), compute the case the site is actually
about, the Khafre claim on the real image (A3), state the physical envelope as a declared product of factors with the
margin that survives it (A7), and the proof is as tight as a proof about a model can be. The rest is bookkeeping,
sample size and cutting.

## A. The proof

### A1. Put the background into the tight layer: the identity assumption costs 4 × 10⁻⁵, not a factor 18

The tight bound (P2-25) is stated "with the reference covariance taken as the identity"; the certificate (P2-30)
carries the background by proving λ_min(C₀) ≥ (1 − √Q_bg)² = 0.056 from the background's Frobenius phase energy,
which gives every pixel the background's peak phase and so loses a factor 1/√0.056 = 4.2 in TV and 18 in E_min for
nothing. The background is long-wavelength motion nearly shared by the scene, which is exactly what the tight layer
knows holds almost nothing. Two routes, do the first, keep the second as the analytic backstop:

1. **Exact, per line, with the background in both worlds.** P2-27's `line_kl` already computes KL exactly for one
   line along track with the full range band. In the model the lines are independent exactly (white scatterers; with the
   whole range band each range pixel's scatterers reach only their own line), and the narrower image is a function of
   the wider one (Theorem A), so Σ over lines bounds the SLC's KL with no Fisher information, no remainder and no
   floor. Generalise `line_kl` from the Gaussian family to (a) the imprint's actual map on the line, direction by
   direction, and (b) a common phase ψ(x, t) in both worlds, with KL = tr(C₀⁻¹C₁) − n_b − log det(C₀⁻¹C₁). Lines the
   imprint does not touch contribute zero. With ψ present the window along track is a truncation; show the sum is
   insensitive to the window's length (the check a reviewer will ask for). Run it for the microseisms and for the
   lorry's own direct wave (A2). This replaces "reference covariance I" everywhere with one exact number per case, and
   the four-number stat row with one.
2. **The lemma.** ‖C_bg − I‖_F ≤ √F_bg + ε_bg by the tight layer's own bound applied to the background's phase, so
   λ_min(C₀) ≥ 1 − √F_bg − ε_bg. For a plane Rayleigh wave at the regional level (4.9 × 10⁻⁸ m/s, 0.2 Hz, 3,000 m/s)
   across the whole 3.5 km scene the slow-limit formula gives F_bg = 1.5 × 10⁻⁹, √F_bg = 3.8 × 10⁻⁵ (confirm with
   `fisher_sinusoid` on the pixel grid). The noisiest stations: √F_bg = 2.8 × 10⁻³, so that row gets a
   background-included certificate too (the Frobenius floor fails there, √Q_bg > 1; this one does not). The imprint's
   first-order term changes by at most 2ψ_max√F plus a gradient term (the background's phase varies by
   2π × 40 m / 15 km across the footprint), both below 10⁻⁴ relative here.

Consequences to compute, not assume: E_min (Theorem C) rises from 9.8 × 10⁻¹⁰ toward the identity value
1.7 × 10⁻⁸ m⁴; the static margin from 4.6 × 10⁵ toward 1.9 × 10⁶ in amplitude; the thin-roof-under-a-lorry rows
(0.45, 0.97) move by the square root of the floor's improvement, which may or may not take them past 1 once the
lorry's own wave is in the reference (A2).

### A2. The lorry rows carry the microseisms in the reference but not the lorry's own wave

In the lorry world the dominant motion common to both hypotheses is the lorry's direct wave, not the microseisms.
P2-30's `background_included` uses the three regional harmonics only, and P2-25's lorry rows use the identity. So the
static rows and the lorry rows are not the same theorem, and only the oracle's lorry rows are clean (there the
background cancels exactly). The lorry's direct wave at 15 m has a peak phase of 1.7 × 10⁻⁴ rad; carried unattenuated
as 1/√r to the scene's edge its any-frequency cap is √F ≤ 0.76 (0.12 within 100 m), so even the crude lemma gives
λ_min ≥ 0.24, four times the certificate's floor, and A1's exact route does far better. Until it is computed, say on
the page that the tight and certificate lorry rows take the lorry's wave out of the reference.

### A3. The Khafre claim is the one case not computed, and a reader goes there first

P2-25's `claimed` block is "an estimate by scaling, not computed" (relative imprint 0.36 of the bench room's). The
site's question is the claim. Compute it: P2-04's kernel on `bench-khafre-claim` (composed, `listed: false`) or on
`giza-deep`, then the tight layer and the oracle for the 80 m cube at 1,220 m, whose imprint spreads over a
kilometre. For the oracle use the real 2025 image's pixel powers over that footprint, Σ_j SNR_j = Σ|z_j|² / NESZ per
cell, the pyramid's bright faces included: a real-data anchor that needs no texture model, and the strongest single
number the site could hold. The product lives only on the desktop; a downsampled power map (a few hundred kB) exported
into `sites/acquisitions/` makes the run reproducible in the cloud. Do the 2022 pass as well.

### A4. Every realisation, not the average

Two places take an expectation where the statement should hold for every case:

- P2-25's static rows use F = Σ_d A_d / n: waves from every direction with independent phases, each with 1/n of the
  power. That is the expectation over the field's phases and over how its power is split. The triangle inequality on
  ‖E‖_F gives, for every realisation and every split, F ≤ Σ_d (A_d + |B_d|): at most n = 24 times the mean, √24 = 4.9
  in TV. Free at these margins, and it turns "on average over the ambient field" into "for every realisation".
- P2-29's Jensen step bounds the oracle's TV averaged over the scatterers' phases; the plateau is one realisation.
  Δ² is a sum over the scatterers in the imprint's footprint (thousands within its half-peak, far more in its tail), so
  it concentrates: the check's spread of ±7% at 802 scatterers scales as 1/√N. One sentence on the page, either as
  Var/mean² ≈ 1/N_eff or in the Markov form: with probability at least 1 − 1/t over the phases, Δ² ≤ t E Δ²; at
  t = 10⁶ the microseism row is still ≤ 1.6 × 10⁻³.

### A5. Two sentences that pre-empt two objections

- **Composite alternatives.** The bound is computed for one specified room. For any prior over rooms, layouts, depths
  and sizes, TV(mixture, P₀) ≤ max over the members, so the largest-imprint member (the thin-roof room) covers every
  mixture. Say so once.
- **Several images.** KL adds over independent images: N images give TV ≤ √N times the one-image bound (Pinsker), and
  the oracle's Δ² grows by N. State how many passes the 2025 announcement used; with two the bound is 1.4 times.

### A6. Say which layer covers which data, and take the headline off the Doppler-to-time check

The stationary-phase model was checked against pulses to 3.6% from 5 to 150 Hz on a 0.5 s aperture; at 0.5 Hz the
model holds 1.47 times more. "On the side of more information" is an observation, not a proof, and a reviewer will say
so. Either prove the direction or route the universal statement through the oracle, which is pulse by pulse and needs
no such relation, and let the SLC layers carry the check as a note. A three-row table settles what covers what: raw
echoes, any texture → the oracle (Theorem D, at a stated SNR); the SLC under speckle → the tight layer and the
certificate (the SLC is a function of the raw data, so the oracle bounds it too); any processing of either → Theorem A.
The phrase "the method is told the ambient motion" in P2-25 should go: with C₀ = I nothing is told; A1 is what telling
it would mean.

### A7. The physical envelope as a declared product, with the margin that survives it

The third review's next step, and the thing that makes the boundary a boundary. One table on the proof page:

| factor | value | status |
|---|---|---|
| local level at the plateau over Kottamya | ×10 | assumed (no open measurement) |
| site amplification | ×3 | assumed |
| cover or layering under the room | ×? | run one layered bench at low frequency rather than assume |
| SNR for the oracle | +10 dB | assumed |
| far-field attenuation of the scattered wave | none | generous by construction |
| grid error | 6–11% | measured (P2-26) |
| wave direction between samples | sec 15° | see A8 |
| product, in amplitude | about 300 | |

Against margins of 4.6 × 10⁵ (or ~2 × 10⁶ after A1) the ground's own trembling survives every declared allowance
with three orders of magnitude to spare; under a lorry (2 to 4) it does not. That sentence is the airtight form of the
boundary, and the table is what a proponent must argue with.

### A8. Small rigour items

- **The direction guard for a cube.** sec δ covers a + b cos 2φ + c sin 2φ; the bench room is a cube, so its
  imprint's dependence on wave direction has cos 4φ terms, for which the guard is sec 2δ: 1.035 at 7.5° (P2-25,
  P2-30 use 1.01) and 1.155 at 15° (P2-32 uses 1.035). State the harmonic content assumed or use the larger guard.
- **`perturbation_certificate`** uses KL ≤ ρ²/(2(1 − ρ)²); the derivation supports ρ²/(2(1 − ρ)) (Σe² = ‖E‖_F² ≤ ρ²
  and e − log(1 + e) ≤ e²/(2(1 − ρ))). Valid either way; make the code and the page's TV ≤ ρ/(2(1 − ρ)) agree with
  the derivation, or say why the extra (1 − ρ) is kept.
- `finite.py`'s docstring cites `PROOF.md`, which does not exist.
- ROADMAP's State list, item 4, still says "reached by the best detector in simulation", which the second review
  withdrew; the audit and the page say "a bound, not an attained value".
- The relative-vertical-geometry row of the claims table says "not shown"; it is untestable with single-level layouts
  (P2-31 says so); say "not testable yet", not "not shown".

### A9. Nine images cannot carry "the published picture, read blind, does neither"

The headline of the shape test rests on 5 of 9 named against 3 of 9 by chance (p = 0.14) and on centroid distances of
5.6 m against 6.2 m. The reference detector's 9 of 9 is decisive (p = 5 × 10⁻⁵); the published method's is not. Either
run 30 realisations (the kernels are cached; the synthesizer is the cost) or change the headline to "no better than
chance was shown" until then. The same for the "44% at the real level" cell, which is 4 of 9.

## B. The site

1. **Merge `main` first.** The branch is 10 commits behind (Sacsayhuamán, P2-33 to P2-35). `git merge-tree` shows
   conflicts in `METHOD_AUDIT.md`, `web/src/pages/satellite.astro`, `sim/katabasis/export/runs.py` (both sides extend
   `PUBLISHED`), `web/public/data/runs/index.json` and `docs/*`. Regenerate the index with the export (runs only on a
   cloud machine, per the roadmap's note), rebuild `docs/`, run both suites and the stale-export test before pushing.
2. **Step 10 on the overview is an essay.** 1,116 words in one column with the right column empty, where steps 1 to 9
   run 67 to 281 words beside their pictures; the FAQ under "could anything change it" is 911; the overview grew from
   3,381 to 4,954 words. On a phone step 10 is five screens of prose before any figure. Keep the headline, one paragraph
   (what an image is made of; shared motion leaves nothing; the cap), the Routes figure, the Ceiling chart, one line on
   what stays open and the link to `/proof/`. Everything else is already on the proof page.
3. **One headline number, one unit.** The proof page's stat row shows 5.017%, 1.3 × 10⁻⁹, 1.6 × 10⁻⁶ and 27%: two
   units, three layers. After A1 lead with the one exact number; until then lead with the tight number and give the
   certificate as "the loose check". The overview says one in 770 million and the proof page says 5.017%; a reader
   cannot tell they are the same claim.
4. **The three-layers table** notes "background included" on three certificate rows and "reference covariance I" on
   the fourth; after A1 every row says the same thing.
5. **The Routes figure** has a stray "×" glyph beside "Its echo" and arrows that do not meet the cards; check the phone
   stacking.
6. **The satellite page's status line** says "revised after two independent reviews"; the audit lists three.
7. **Precision theatre.** "5.017%", "5.110%", "0.0559", "1.31 billion pixels": round to what the reader uses (5.02%,
   "1.3 billion"), and say "0.02 above chance" where that is the point.
8. **The requirement chart** has seven legend entries and its 20 and 80 Hz lines overlap; keep 0.2 Hz, 80 Hz and the
   two any-pattern lines, with direct labels.
9. `SectionMap.astro` declares `fmt` and never uses it (the one hint from `astro check`).
10. **The answer section's** 456-word sub-paragraph and the mechanism card's 300 words are prose where the rest of the
    page is dashboard; the five-theorem cards are the right form, copy it.
11. **The audit** takes a "Fourth review, 30 September 2026: corrections adopted" entry as before, and item 17's
    "reference covariance I" wording changes with A1.

## C. In what order

1. Merge `main` (B1).
2. A1, route 1: the exact per-line layer with the background, microseisms and the lorry's wave (A2). One run;
   it rewrites the stat row, the three-layers table and the certificate's floor.
3. A3: the claim on the real image, both passes.
4. A7: the declared-range table, with the layered bench run if it can be had in a day; otherwise the row says assumed.
5. A4, A5, A6, A8: sentences and small fixes, an hour.
6. A9: thirty realisations, or the softer headline.
7. B2 and B3: cut step 10, one headline number.
8. The rest of B, then the audit entry, then publish.

## Appendix: the background's own information, reproduced

Run in `sim/` on the branch (`uv run python`):

```python
import numpy as np
from sarsim.acquisition import DwellGeometry
from sarsim import information as inf
g = DwellGeometry.from_record('giza-20250827')
k0, R, A = 4 * np.pi / g.lam, 3535.0, np.pi * 3535.0 ** 2
for label, v0, f, c in [('regional microseisms', 4.9e-8, 0.2, 3000.0), ('noisiest stations', 3.525e-6, 0.2, 3000.0),
                        ('1-3 Hz', 2.9e-8, 2.0, 1690.8), ('3-8 Hz', 2.5e-8, 5.0, 1690.8)]:
    kw = 2 * np.pi * f / c
    F_slow = (g.R0 / g.V_platform) ** 2 * inf.cells_per_m2(g) * 0.5 * (v0 * kw) ** 2 * A   # a plane Rayleigh wave over the scene
    print(f'{label}: F_bg (slow limit) = {F_slow:.2g}, sqrt = {np.sqrt(F_slow):.2g}')
v15, f = 1e-4 * np.sqrt(2), 51.8                       # the lorry's direct wave at 15 m, spreading as 1/sqrt(r), no attenuation
K15 = v15 / (2 * np.pi * f)
for Rmax in (100.0, 3535.0):
    F_cap = 2 * inf.cells_per_m2(g) * k0 ** 2 * 2 * np.pi * 15.0 * K15 ** 2 * (Rmax - 15.0)   # 4 N <Phi^2>, any frequency
    print(f'lorry wave to {Rmax:.0f} m: F <= {F_cap:.2g}, sqrt = {np.sqrt(F_cap):.2g}; peak phase {k0 * K15:.2g} rad')
```

Output on this machine: regional microseisms F_bg = 1.5 × 10⁻⁹ (√ = 3.8 × 10⁻⁵); noisiest stations 7.6 × 10⁻⁶
(2.8 × 10⁻³); 1–3 Hz 1.6 × 10⁻⁷ (4.0 × 10⁻⁴); 3–8 Hz 7.5 × 10⁻⁷ (8.7 × 10⁻⁴); the lorry's wave √F ≤ 0.12 within
100 m and 0.76 to the scene's edge, peak phase 1.7 × 10⁻⁴ rad. The slow limit is the right form at 0.2 Hz over a
15 km wavelength (P2-25's own check: 1.37 against 1.33 exact at 0.5 Hz); confirm on the pixel grid before quoting.
