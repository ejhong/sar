# Design — mist and deep water

An instrument more than a magazine: small dense type, statement headings, and one image allowed to be beautiful — the
underworld. Cool and muted: mist for reading, deep slate-blue where the drama and the instruments are.

- **Mist (reading).** Page `#e8eef2`, sheets and cards `#f5f8fa`, borders `#d3dde5`/`#b7c6d2`, ink `#1c2833`, body
  `#3d4c5a`, muted `#6b7d8e`.
- **Deep water (drama).** The top band and footer `#132030` with a faint lapis glow; the front-page hero and the
  viewer ground `#0f1822` with a lamp-lit vignette; instrument panels `#16222e`/`#1b2836` with text `#dfe6ec`; full-bleed
  feature bands (the measured budget) in deep water with light type.
- **The seal.** The golden bough, Aeneas's passage into the underworld and back (the Virgil in the top band): gold on
  lapis, drawn in `src/components/Seal.astro` and `public/favicon.svg`.
- **Type.** IBM Plex Sans (UI, section headings), IBM Plex Mono (data, labels, captions), Newsreader (display titles and
  prose), all self-hosted through fontsource so every reader sees the same page. Headings are statements.
- **Colour roles — one each, everywhere, 2-D and 3-D.** lapis = voids and chambers; verdigris = what an instrument
  recovered; ochre = sensors and sources; porphyry = claimed structures, always dashed; cinnabar = the radar. Three sets:
  on mist (`#2c5aa0 #16876e #a8741f #7a55a6 #aa4b33`), on deep water (`#5b8fe0 #23a07f #be8a2e #9a78d3 #d46e50`), both
  passing the dataviz six checks; and emissive glows in the viewer (`#82b3ff #74d8bc #e8b95c #c0a6da #f0906f`).
- **Front page.** An introduction, not an instrument: the hero tours the underworld on its own (no controls, never
  catching the scroll), then the claim with its timeline, the ladder of instruments (the organising figure: borehole to
  orbit, each with its verdict), the latest results, and the measured budget in a deep band.
- **The underworld.** Each site a geological block diagram: terrain lit low with contours, strata on the faces in their
  conventional patterns, depth ticks, a section plane with a patterned cut face, chambers glowing through the rock with
  plumb lines, pyramids as survey outlines, recovered volumes ray-marched in verdigris, stations in ochre, the wavefield
  as a sheet of light (the echo alone in lapis). The panel leads with the site and a Show control (truth · recovered ·
  waves); features pin on click; a live scale bar; keyboard shortcuts; the view lives in the URL.
- **Reading pages.** One mist sheet with a quiet index; a scoreboard first where there are results; numbered figures;
  every run with a collapsible details panel linking its data and code; terms from the other community underlined with
  their definition on hover (`src/data/glossary.ts`).
- **Charts.** Hand-made SVG at build time: 2 px lines, hairline grids, legends for two or more series, hover read-outs.
  Tomogram sections diverge about the background: slower in verdigris, faster in umber, true chambers outlined in lapis.
- **Honesty in the interface.** Every figure says simulated, measured, claimed or planned; every run names its date and
  commit; every material number traces to a source or says assumed.

Tokens live in `src/styles/tokens.css` and `src/viewer/engine/theme.ts`; keep them in step.
