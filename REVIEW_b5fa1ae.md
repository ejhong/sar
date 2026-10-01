# Review at b5fa1ae (revision 3): good enough to publish, with one wording fix; three things would make it great

1 October 2026. Reviewed on `claude/upbeat-pascal-swk910` at `b5fa1ae`: the front page and satellite page answer
(`ClaimVerdict`), `BENCHMARK.md` revision 3 and its erratum, the regenerated P2-36 report, the audit's seventh and
eighth reviews, `sarsim/radiometry.py`, P2-37, P2-38 with `katabasis/seismic/halfspace.py`, and the rebuilt pages
(served from the branch and read at desktop and phone widths). The simulation suite passes (172 tests), `astro check`
is clean, the 11 web tests pass, and a fresh build reproduces the committed `docs/` except for hashed asset names. The
18 volumes and 2 wavefields of the bench site are byte-identical to `main`'s, and the export now keeps published
entries whose arrays are elsewhere.

**Verdict.** Publish. The answer is supported by three independent legs that a careful reader can follow from the
front page to the runs: the claim's own structures computed at their depths leave any reader of one image within
2.4 × 10⁻⁵ of chance; a privileged reader told every scatterer is 2,100 times short even at the most favourable noise
floor the maker specifies; and the published method, read blind on images of known ground, is at chance even with the
cavities' effect amplified a thousandfold, where a detector told what to look for names every world. "What it rests
on" and "what could still change it" are the right two boxes and say the right things. Every point of my afee349 review
is handled, and handled correctly (the list is in §3). One sentence on the front page must change before it goes live,
and three additions would move it from compelling to great.

## 1. The one fix before publishing

The first stat on the front page and the satellite page says "the most any way of reading one image of Khafre can do
better than a coin toss at telling whether the claimed shafts and void are there: 1 in 42 thousand". The number is
the bound on detection rate minus false-alarm rate (2.4 × 10⁻⁵); the margin over a coin toss is half of it. This is
the slip the second review corrected everywhere else, back on the most-read line of the site. Either print the
overview's wording ("raise its rate of finding them above its rate of false alarms by at most one in 42 thousand") or
keep "coin toss" and print 1 in 84 thousand. The overview's step 10 and the proof page have it right.

## 2. Three things that would make it great

1. **Say the margin after the allowances where the objection will come.** "The trembling measured 67 km away" is the
   first thing a sceptic reads in "what it rests on". The run already answers it: with a local level ten times the
   regional, site amplification three times and the cube's moment at 1.5, applied to the background too, the
   structures are still 96 times short for one image and 46 times for the oracle. Put those two numbers in "what could
   still change it" and the objection is met on the same screen.
2. **Compute the open corner for Khafre itself.** The pyramid is a bright persistent target (the 2 km map's brightest
   cell is 29.5 dB over the median), the shafts lie under it, and buses pass on the plateau road. That is exactly the
   case the bounds leave open, in the claim's own geometry. Everything needed exists: the shaft kernel from P2-38, the
   benchmark's point source with its distance as a parameter, and the measured brightness over 2 km. One run with the
   lorry at the road's nearest approach to Khafre, the shafts' dynamic response at the worst sampled frequency and the
   oracle weighted by the real brightness would turn "the one case these bounds leave open" into a number, open or
   closed. Not a blocker for publishing; the most valuable next computation.
3. **Let the first screen carry the answer.** The hero still reads "The radar methods tested do not show what lies
   below", the scoped statement from the first review; the answer "In our assessment, no" is one screen down, and on a
   phone two. The hero's kicker already asks "Can a satellite see beneath the ground?"; its headline can now answer the
   question readers bring in the verdict's own words, with the scoped sentence as the line beneath.

## 3. What I verified, and how each afee349 point was handled

- **The noise floor.** ICEYE's documentation 6.0.8, Table 2-11, read again today: Dwell −26.7 to −15.6 dB, Dwell Fine
  −23.7 to −12.6 dB, Dwell Precise −20.7 to −11.3 dB, scene-centre values. The 2025 product is Dwell Fine, so the
  headline 18.4 dB is the measured −5.3 dB over −23.7 dB, the most favourable specified for the mode; 6.0.0's −18 dB
  and the measured lower bound stay as scenarios; the 2022 pass carries its own 18.8 dB. The erratum's account (a
  traceable figure misapplied, then wrongly called untraceable) is accurate and the module's docstring says so. The
  mode comes from the product code (SLEDF, SLED); say that in the record, since only ICEYE can confirm it.
- **Revision 3.** Worlds unchanged (the hash says so, so P2-37's images carry over); the score test with named disjoint
  seed streams, 400 calibration and 100 + 100 evaluation grounds, 73% (63 to 81%) at 15% (9 to 24%) false alarms, AUC
  0.91 against 0.93 predicted, chance at the real level; the 15% reported as measured with the Kolmogorov–Smirnov
  comparison rather than tuned away, which is the honest choice; the four stale verification rows recomputed (0.12 to
  0.18 of the certificate, as before); the worst realisation's pass mean corrected; the floor's positive part with the
  noise fallback, tested.
- **P2-38.** The deep void through the damped half-space by wavenumber integration, equal to Okada at low frequency
  within 5 × 10⁻⁴ (tested); its dynamic integral 1.99 times the static one, Q = 100 and 1,000 within 0.5%. The shafts
  by the static solver near the surface (truncation at 100 and 200 m changes the near field by 0.9%, so the top hundred
  metres rule the surface mark, as expected) and Eshelby cylinder columns beyond, calibrated by the larger of the two
  ring ratios (1.47). The shafts dominate (1.0 × 10⁻¹⁰ m against the void's 2.5 × 10⁻¹²): a Kirsch-sized number for
  a 5 m hole under 1.5 × 10⁻¹¹ strain, so the physics is in proportion. Both passes at their own SNR; the oracle
  weighted by the 2 km brightness map and, beyond it, by the brightest cell (a bound five times the at-median value,
  both near chance). The allowances now scale the background too, and the floor falls back on the receiver noise, which
  is why the margins after allowances are 96 and 46. The cube allowance cites the measured 1.22 to 1.36; the loading is
  "at most as at the surface"; the shafts are in the first sentence.
- **P2-37.** The told detector names 100% at the positive control (location 16 of 16, shape 16 of 16, p = 0.004 within
  grounds), 62% at ×200; the published method's blind map is at chance at every level (AUC 0.52 to 0.53, 33%, 8 of 16
  on both pairs). Eight grounds are few, but the control's strength makes the null result decisive. The paired reader's
  presence AUC is written as 0.0 in the summary; the finding says it is not scored, so write it as null.
- **The proof page's statements.** The mixture bound is now the prior's average of its members with the correct reading
  of a hard pair; the Gaussian divergence states equal means and ‖E‖ < 1; both floors use the positive part; √N is the
  KL ceiling's with independence to be justified; the location pair's equal-prior accuracy is 52.3%; detecting a room
  does not imply telling its depth. The Routes figure now shows the benchmark's world WA and ends "found ≤ false
  alarms + 4.2% (P2-36)" on the proof page and the overview; the claim's two rows are in the benchmark's table; the
  bright corner leads the open table with its number.

## 4. Small things, none blocking

- The front page's step 10 says "at most 17% found beside a bouncing lorry" where the proof page's table says 16% for
  the same row (0.116 + 0.05, rounded two ways). Round once.
- The satellite page's index has 22 entries and its lede six lines; it is the evidence page, so that is tolerable, but
  "The claim, answered" and "Where it stands" could be the only two before a fold.
- `katabasis.runs.memo` keyed on inputs is the right fix for the stale rows; the one-line note in the erratum is enough.

## 5. Order

The coin-toss line (minutes). Publish. Then the two numbers in "what could change it", the hero headline, the
Khafre-corner run, the paired AUC as null, and the two rounding and index items.
