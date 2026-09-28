# Underworld — can a satellite see beneath the ground?

Code and data behind **[Underworld](https://ejhong.github.io/sar/)**, which asks what can see beneath the ground, and
how far down, by building ground whose contents are known and simulating the instruments over it. Its answer for the
satellite claim is shown step by step at [/sar/pictures/](https://ejhong.github.io/sar/pictures/).

The claim under test is the single-image SAR "Doppler tomography" of Biondi & Malanga (2022, retracted 2026) and the
2025 Khafre "underground city" announcement. The lab answers it in two phases:

1. **Seismic tomography, simulated.** Full 3-D elastic waves through composed test sites and Giza, recorded by the arrays
   near-surface surveys deploy, inverted by travel-time, ambient-noise and full-waveform tomography, and scored against the
   truth. The reference ceiling.
2. **The satellite over the same ground.** An ICEYE Spotlight Dwell with the real Giza parameters, simulated over
   ground whose contents and motion are known, and checked against the real product.

**The answer is no, for three independent reasons** (`METHOD_AUDIT.md`, "Phase 2"). Registering one image's looks
against each other, as every single-image method does, is blind to ordinary ground moving: the looks share their
spectrum sample for sample (P2-03). A chamber changes the ground's own motion by 6.6 × 10⁻⁵ of it, 10⁴ to 10⁸ below
what one dwell could detect (P2-04, P2-05). And the published depth is the frequency of the ground's motion, 6.76 m per
hertz on the real geometry, whatever is underground (P2-06). Run whole, over simulated shaking ground, the published
method draws the same picture with and without the chamber (P2-07).

The lab's Python package keeps its first name, `katabasis`. The first investigation (the published pipeline
reimplemented and run on two real ICEYE dwells) is kept intact in
`sim/legacy/` and its site at [/sar/archive/](https://ejhong.github.io/sar/archive/).

## Layout

```
sites/            materials.json (sourced) and one composition per site: bench-void, bench-shafts,
                  bench-khafre-claim, giza
sim/              the Python lab (uv)
  katabasis/
    compose/      materials, site schema, shapes, terrain, heterogeneity, voxeliser, DEM crops
    seismic/      3-D elastic solver, arrays, inversions (Phase 1, in progress)
    export/       everything the site reads, written to web/public/data/
  sarsim/         the SAR toolkit from the first investigation (Phase 2 builds on it)
  legacy/         the first investigation: experiments, fieldwork, results, catalog, tests
  tests/
web/              the site: Astro + Three.js; builds into docs/
docs/             the built site, served by GitHub Pages from main
```

## Run

```bash
cd sim && uv sync && uv run pytest            # simulation and checks
uv run python -m katabasis.export             # regenerate web/public/data/
cd ../web && npm ci && npm run check && npm test && npm run build   # site → docs/
```

`METHOD_AUDIT.md` records what is and is not established; `ROADMAP.md` is the plan.
