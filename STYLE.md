# Presentation style

Eugene's stated preference, and the rule for every page in `docs/`.

## The shape of it

**Dashboard, not essay.** The page is an instrument panel: stat tiles, compact result cards,
tables and figures on a grid. A reader should be able to stand back and take in the state of the
investigation without reading a paragraph. Narrative is available, but it is collapsed behind
the numbers rather than wrapped around them.

## Voice

- **Headers are statements, not questions.** "The depth axis repeats every 13.6 m", not
  "Does the depth axis repeat?". A reader scanning only the headings should get the findings.
- **No flourish.** No rhetorical build-up, no "three different claims", no section numbers that
  pose a question and answer it two screens later. Say the thing.
- **Findings as bullets.** A finding is a short claim plus its number. Long narrative paragraphs
  belong in the method details, collapsed.
- **One idea per line.** If a sentence needs a semicolon it is two bullets.

## Typography and density

- **Small type, high density.** Body around 14px, headings that do not dominate the viewport.
  The reader should see several findings per screen without scrolling.
- **No hero-scale display type** except the page title once. No 60px numbers unless the number
  is the single headline of the page.
- **Tight vertical rhythm.** Sections separated by a rule, not by empty screens.
- **Tables over prose for numbers.** Every measurement lives in a table or a stat row.
- Beautiful still matters: restrained palette, real typographic care, generous horizontal
  measure. Dense is not cramped.

## Structure

- **The viewer leads.** The three-dimensional underworld viewer is the most interesting artefact
  on the site and belongs at the top, not linked from a card halfway down.
- **Real data before simulation.** Measurement first; simulated benchmarks are supporting.
- **Collapse the long tail.** Method, metrics and reproducibility go in `<details>`.
- Each result block: statement heading, one-line finding, the number, the figure, limits.
- **Lead with a stat row.** The headline numbers of the whole investigation sit at the top, each
  traceable to one experiment.
- Grid layouts over stacked full-width sections. Two or three columns where the content allows.

## What the reader must always be able to tell apart

1. Published claims by Biondi, Malanga and the Khafre project.
2. The public derivative protocol (Seyfzadeh, v1.5 and v1.7).
3. This repository's reimplementation and its own choices.
4. Simulation versus measurement on real acquisitions.

Never let a rendering flatter the method. If a figure is smoothed or colour-scaled the way the
published slides were, say so in the caption.
