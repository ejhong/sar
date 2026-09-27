# Design — limestone and basalt

The craft of ejhong/knots (and The OM Project before it), in a register drawn from the subject: an instrument more than a
magazine, small dense type, and one image allowed to be beautiful — the underworld.

- **Surfaces.** Day is the page: limestone `#f4efe5`, fresh-cut sheets `#faf7f0`, a sandstone band `#e4dccd` with the
  title and a line of Virgil in Georgia italic, 12 px cards with 1 px `#ddd3c2` borders. Night is the underworld: basalt
  `#0f1214` behind the viewer with a soft lamp-lit vignette, slate panels `#1a1f22` with mono labels.
- **Type (system only, no web fonts).** Sans for UI and headings (title 20/700, sections 16/700), Georgia for prose
  (13.5–15 px), SF Mono for labels, data and captions (9.5–11 px). Headings are statements, never questions.
- **Colour roles — one each, everywhere, 2-D and 3-D.**
  lapis `#2c5aa0` / night `#82b3ff` = voids and chambers (the truth or the survey);
  verdigris `#16876e` / `#74d8bc` = what an instrument recovered;
  ochre `#a8741f` / `#e8b95c` = sensors and sources;
  porphyry `#7a55a6` / `#c0a6da` = claimed structures, always dashed;
  cinnabar `#aa4b33` / `#f0906f` = the radar (Phase 2).
  Everything else is stone, ink and paper. The day palette passes the dataviz six checks against `--sheet`; night
  colours are emissive glows in the viewer, not chart fills.
- **Charts.** Hand-made SVG, drawn at build time, 2 px lines, hairline grids, a legend for two or more series, hover
  read-outs. Tomogram sections use a diverging scale about the background: slower in verdigris, faster in umber, a
  neutral grey between, true chambers outlined in lapis, stations in ochre.
- **The underworld.** Each site is a geological block diagram: terrain lit low with contours, strata on the faces in
  their conventional patterns (limestone bedded blocks, sand stipple, marl dashes, igneous crosses, masonry large
  blocks), the water table dashed, depth ticks on each face, a section plane that reveals a patterned cut face, chambers
  glowing through the rock with plumb lines to the ground, pyramids as survey outlines, recovered volumes ray-marched in
  verdigris, stations as ochre points, the wavefield as a sheet of light (the echo alone in lapis). Slow turntable when
  idle; motion stops under `prefers-reduced-motion`.
- **Layout.** The instrument: a basalt panel beside the viewer card. Embedded in a page the viewer never steals the
  scroll (ctrl/⌘ to zoom, two fingers to turn); on phones the viewer is pinned above and the panel scrolls below.
  Reading pages are one limestone sheet with a quiet index, dashboard pieces on a grid: stat rows first, then figures.
- **Honesty in the interface.** Every figure says simulated, measured, claimed or planned; every run names its date and
  commit; every material number traces to a source or says assumed.

Tokens live in `src/styles/tokens.css` and `src/viewer/engine/theme.ts`; keep them in step.
