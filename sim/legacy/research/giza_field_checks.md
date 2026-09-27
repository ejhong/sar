# Real-image checks and survey geometry — 20 September 2026

This advances the known-void inventory to actual native-image checks. It does
not produce a chamber detection or validate the absence of underground columns.
The first Giza comparison remains exploratory: reference surveys were inspected
while designing the work.

## What is measured

Both supplied native ICEYE complex arrays were read in bounded chunks, with no
nonfinite I/Q values encountered. The full I and Q sample arrays have separate
SHA-256 digests. These identify the read values, not a vendor-authenticated file
checksum. Small intensity overviews and detail crops are published. Originals
remain at `~/tmp/sar/`; no full raster or large cache was copied into the repo.

The [Giza report](../results/field/giza/report.json) and
[Sacsayhuaman report](../results/field/sacsayhuaman/report.json) retain acquisition
metadata, code/design hashes, sample selections, display bounds, raw spectrum
curves and complete translation-trial records. They rebuild the page without
access to the original private files.

## Survey placement

The [Glen Dash Foundation survey index](http://www.dashfoundation.org/Surveys-at-Giza.htm)
links downloadable control and feature points. The
[2012 survey KMZ](http://www.dashfoundation.org/downloads/archaeology/survey/GDFS-2012-for-GE-21vii18.kmz)
contains 1,539 point records. **KML point altitudes are zero**; the actual survey
elevations appear in each point's description, identified as height above mean
sea level. Thirteen surface points near the Eastern Cemetery were transcribed
with their original descriptions, including qualifications such as “Maybe”.
The registry preserves the source checksum.

ICEYE's [metadata specification](https://sar.iceye.com/6.0.0/productFormats/metadata/)
uses WGS84 ellipsoidal heights. A four-node extract of the
[NGA EGM96 grid distributed by PROJ](https://cdn.proj.org/us_nga_egm96_15.tif)
gives a local geoid undulation near +15.46 m, using bilinear interpolation.
The conversion is `h = H + N`, following
[GeographicLib's height convention](https://geographiclib.sourceforge.io/C++/doc/geoid.html).
This does **not** establish the exact vertical tie of the GPMP survey to EGM96.

Projecting these points with versus without that conversion changes range
placement by about 94 native pixels, or 41 m in a local flat-ground conversion.
These are **computed projection differences**, not measured registration errors.
The result is larger than the chamber dimensions, illustrating why the datum
cannot be ignored.

The overlay uses the native RPC with no fitted offset. Surface structures are
recognizable at monument scale, but the selected historical ground/casing edges
are not unambiguous radar point reflectors. No manually accepted correspondence,
held-out RMSE, fitted survey-to-radar transform or buried target mask is claimed.
The ±5 m interactive adjustment is a sensitivity scenario, not a confidence
interval. Elevation, layover, modern alterations and landmark identity still
need resolution before chamber-scale scoring.

## Timing and spectral support

Giza collection duration is 24.762808 s. Its supplied zero-Doppler image span is
0.712901 s. The latter tracks where the focused image rows fall and must not
become the aperture duration. The native output intervals and supplied endpoint
timestamps also differ by about 1.1 output rows, which is retained in the report.

The [NGA SarPy ICEYE reader](https://github.com/ngageoint/sarpy/blob/master/sarpy/io/complex/iceye.py)
evaluates Doppler-rate coefficients relative to the mid-range delay. On this
Giza file the centre value is about −5,662.5 Hz/s, yielding a nominal
bandwidth/rate duration near 24.49 s. Evaluating the polynomial at the absolute
range delay would be erroneous. This agreement is a consistency check, not a
validated product-specific time mapping. The generic reader and supplied SICD
also describe image-plane conventions that must be reconciled with native HDF5.

Three 24-column strips per acquisition preserve the full native azimuth extent.
Their Fourier power spectra retain the zero-centred metadata bandwidth as a
diagnostic reference. Roughly 0.6–0.8% of measured power lies outside that band,
depending on scene and strip. No taper, deramp, tuned recentering or physical
vibration-frequency axis is applied. The curves and standalone SVGs are saved.

## Real-texture translation control

The design specifies one 512 × 512 crop per scene, eight known translations,
25 fixed interior patches, the existing 32 × 32 complex tracker and a 0.02 px
95th-percentile error tolerance. Translation is applied independently using
SciPy's Fourier-shift function. Every patch, including failures, is published.

In each scene six input groups meet that tolerance. Inputs of (−0.25, +0.4) px
and (+1, −1) px fail; their group 95th-percentile errors are about 0.03–0.04 px.
The ±0.005 px controls are recovered much more closely. This reveals finite-
window/texture bias and verifies sign/axis handling on real image content. It
does not measure ambient motion, add independent noise, or validate the phase
response to actual surface vibration. Patches are dependent, and the engineering
tolerance is looser than the legacy 0.005 px motion gate.

## Models and remaining evaluation

Idu and Hetepheres have dimensioned local 3D survey models, source citations,
explicit datums and simplified geometry. Their measured versus derived roof,
floor and shaft depths remain separate. The Hetepheres pit is drawn at its
bottom width through its full height, a substantial disclosed simplification.
The models represent envelopes of historically excavated spaces; acquisition-
date cover/fill is unknown. Osiris remains a documented reference without a
fabricated arrangement of its lower chambers. No Sacsayhuaman underground model
is invented in the absence of a registered survey.

Next: accept identifiable independent surface correspondences, quantify held-out
placement error, qualify survey height/fill, and validate the dwell-product
motion adapter. Then freeze target/flank coordinates and physical assumptions,
save unrendered predictions, and score misses, depth/position/shape error and
competing surface explanations. Unknown ground is not a confirmed negative;
repetition of a pattern does not exclude a stable artifact. Shallow known-void
performance would not automatically validate kilometre-scale Khafre claims.
