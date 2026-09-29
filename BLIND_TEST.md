# A frozen, blinded calibration test

Written 29 September 2026, before any product for it exists. It is the test the lab's own results cannot replace: a
single-image method (the published focusing, the gated reconstruction, or any other) run frozen on products whose
contents neither its operator nor the scorer knows, and scored afterwards against the truth. The lab's bound (P2-25)
and requirement curve (P2-27) predict the outcome within their model; the test is worth running because it does not
depend on that model.

## 0. Gates before it, in order

An independent review (29 September 2026) set the path to a real mapper, whichever way the evidence goes. Each gate is
passed before the next is attempted:

1. **Recover motion before interpreting depth.** A measured moving reflector with synchronised instruments on it: the
   frozen processing must recover the withheld motion's place, timing, frequency, amplitude and phase, beside stationary
   controls and ordinary scatterers. This qualifies the chain, not cavity detection.
2. **Establish a cavity-specific surface response.** Surveyed voids and matched intact ground under measured excitation,
   an array measuring amplitude and phase; a forward model fitted on some sensors and frequencies predicts the rest.
3. **Predict withheld cavities with frozen settings** (this document's test): presence, place, depth and size scored
   separately; any focus-coordinate calibration fixed on development examples and tested for transfer.
4. **Only then a voxel map**, from validated estimates, with its resolution, calibration domain and alternatives shown.

The lab's bounds (P2-25, P2-29, P2-30) say which regimes are worth testing: a regime the certificate leaves open (strong
nearby shaking, known exactly) is where a controlled experiment could be designed.

## 1. Freeze, then make

Before any product is made, publish (commit to this repository) and hash:

- the method: code at a commit, its SHA-256, the profile or settings, the rule that places its lines, the depth
  window it reports, and how it is run;
- the scoring rule below, its thresholds and tolerances, and the success criterion;
- the generator of the products (below) with its seed withheld: only the seed's SHA-256 is published.

Nothing frozen may change after the first product exists. A change starts a new test with a new seed.

## 2. The products

Synthetic, from the lab's pipeline on the real Giza dwell geometry (P2-01), each 101 m by 111 m of open desert:

- 24 scenes. For each, the generator draws: whether there is a room (probability one half); if so, its depth (10, 15,
  20 or 30 m to its centre), size (4, 6 or 10 m), contents (hollow, or a granite block), and position (uniform within
  30 m of the scene's centre).
- Shaking: Giza's measured microseisms and 1-3 Hz level (M1-01), as P2-07 and P2-20 shake. Each scene is made twice:
  at the real level, and with the room's imprint boosted to a stated peak phase (2 rad, where P2-07's and P2-28's
  outputs change visibly), because a test only at the real level cannot tell a method that does not respond from one
  that responds too weakly.
- Controls mixed in and labelled only in the sealed key: a motionless copy of one scene, a scene of pure speckle, and
  a scene with a random perturbation of a room's size at its pixels (P2-20's null).

Real products, once chosen by someone outside the lab from published surveys: known voids (the Khufu chambers, the
Osiris Shaft when its position and levels are taken from a survey) and matched surveyed controls, each cut to the same
size, labels withheld from the operator.

## 3. Predictions

For every product the frozen method writes, before the key is opened: a score for "a room is here", and if it reports
one, its position (east, north) and depth. The predictions file is hashed and committed; then the seed is revealed and
the key regenerated from it.

## 4. Scoring, stated in advance

- Detection: the area under the ROC curve of the scores against presence, with its 95% interval by the permutation
  of labels (1,000 permutations).
- Location: a hit is a reported room within 6 m of the true one horizontally and within 5 m in depth; every other
  report is a false positive. Hits and false positives per product, and the depth error of each hit.
- Depth: whether reported depth rises with true depth (rank correlation over the rooms found), since a method whose depth
  is the shaking's frequency reports the same depth for every room (P2-06, P2-28).
- Contents: whether the scores or maps for hollow rooms differ from those for granite blocks beyond the null.

Success, for the boosted set: AUC at least 0.8 with its interval above 0.5, a median depth error of 5 m or less, and a
positive depth correlation. At the real level the lab's bound predicts an AUC indistinguishable from 0.5 for any method
(P2-25: detection minus false alarm at most about 1e-9 under the microseisms); a method that does better at the real
level would contradict the model (speckle, the Doppler-to-time relation, the imprint), which is why the real level is in
the test.

## 5. What each outcome would mean

- Fails on the boosted set: the method does not read the room's motion, however strong; its pictures come from elsewhere.
- Passes on the boosted set, fails at the real level: the method reads strong motion, and the real ground's motion is
  too weak for it, as the bound says.
- Passes at the real level: the lab's model is wrong somewhere, and the finding is the model's error, to be found.
- On real products: a pass on known voids and not on matched controls, blind and on both passes, would weaken the
  lab's conclusions directly (METHOD_AUDIT.md, "What is not established").

## 6. Who runs what

The lab can make the synthetic products and seal the key; the published focusing runs in the cloud, the gated
reconstruction only on the desktop, where its public code is installed. The method's proponents may run their own
method on the same products; the key stays sealed until every predictions file is committed.
