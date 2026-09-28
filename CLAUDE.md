# Underworld — notes for working in this repo

A simulation lab for imaging the underworld (the site is "Underworld"; the Python package keeps its first name,
`katabasis`): can seismic arrays, and a single satellite radar image, see chambers under the ground? Phase 1 simulates seismic tomography over known ground; Phase 2 flies the ICEYE dwell over the same
ground. Owner: Eugene Jhong. Read `ROADMAP.md` (state and next steps) and `web/DESIGN.md` (the look) first.

## Principles

- **Known ground, simulated instruments, honest scoring.** Sites are compositions whose contents are known; solvers never
  see inversions and inversions never see the model; every tomogram is scored against the truth.
- **Every number traces.** Material properties carry a source and a quoted line, or are marked derived or assumed
  (`sites/materials.json`). Never a number from memory.
- **Results are generated, never edited.** Each run writes `sim/results/<id>/summary.json` with a manifest (date,
  commit, runtime, parameters); `katabasis.export` is the only path to `web/public/data/`; a test fails if it is stale.
- **Four things kept apart:** published claims, the derivative protocol, this lab's work, and simulated versus measured.
- **Headings are statements**, dashboard not essay, small dense type (the legacy `sim/legacy/STYLE.md` still holds).

## Layout

- `sites/` — `materials.json` and one `site.json` per site (bench-void, bench-shafts, bench-khafre-claim, giza).
- `sim/katabasis/compose` — materials, site schema, shapes, terrain, heterogeneity, voxeliser, DEM crops, benches.
- `sim/katabasis/seismic` — `elastic3d` (the solver), `analytic`, `arrays`, `picking`, `traveltime`, `imaging`.
- `sim/katabasis/export` — sites, runs, volumes and wavefields for the viewer.
- `sim/experiments/p1_*.py` — the published runs. `sim/legacy/` — the first investigation, intact.
- `web/` — Astro + Three.js; `src/viewer/` the underworld (engine, scene: Block, Volume, Wavefield, sections), `src/pages/`
  (the front page leads with the answer; `pictures` is the plain-language walk-through; `satellite` the evidence).
- `docs/` — the built site (committed; Pages serves it from `main`).

## Conventions

- Site coordinates: x east, y north, z up, metres. Solver arrays are [ix, iy, iz] with iz downward; +vz is down inside
  the solver, traces are reported east, north, up.
- Forces are injected as the adjoint of the receiver read-out (split about a node), so sources and receivers commute.
- Voxel files for the viewer are uint8, x fastest, then y, then depth.

## Workflow

- `cd sim && uv run pytest` (≈1 min); experiments: `uv run python experiments/p1_0N_*.py`; then
  `uv run python -m katabasis.export`.
- `cd web && npm run check && npm test && npm run build` (builds into `../docs`); screenshots:
  `SHOT_CHROMIUM=/opt/pw-browsers/chromium npx tsx scripts/shot.ts /sar/underworld/ ../shots/x.png` against
  `npx astro preview`.
- Publish: commit, push the working branch, then fast-forward `main` (the owner asked for live pushes at junctures).
