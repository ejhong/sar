"""Generate docs/index.html from results/*/summary.json. Add a test = add a results folder and rerun."""
import glob, html, json, os, shutil, time

ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
DOCS = os.path.join(ROOT, "docs")
SITE_URL = "https://ejhong.github.io/sar/"
REPO_URL = "https://github.com/ejhong/sar"

TITLE = "Pillars from Nothing"
DECK = ("A reproducible test of the single-image “SAR Doppler tomography” behind the Khafre underground-city claims. "
        "A simulated pyramid that never moves, on ground that contains nothing, produces the same kind of structures at depth. "
        "Every figure here comes from code in the repository and reruns in minutes on a laptop.")

PLANNED = [
    ("07 · Real data", "The ICEYE Dwell scene of Giza", "Run the identical pipeline on a 25-second ICEYE Dwell Fine acquisition (27 Aug 2025): Khafre versus empty plateau versus the Sphinx versus Cairo apartment blocks, with the same statistics as T1 and T2."),
    ("08 · Known chambers", "Positive control inside Khufu", "The King's and Queen's chambers and the Grand Gallery sit at surveyed positions. A method that maps 648 m wells should locate 10 m rooms 40 m up. Score the tomograms against the survey, against chance."),
    ("09 · Cross-geometry", "Ascending versus descending", "Real subsurface architecture is fixed in the ground. Compare tomograms from two headings; T3 predicts the features will move by the projected-spacing law."),
    ("10 · Learned null", "A classifier that cannot tell", "Train a small discriminator on tomograms from simulated static scenes versus the real pyramid; an AUC near 0.5 is the machine-learning form of T2."),
]

SOURCES = [
    ("Biondi & Malanga 2022, <i>Synthetic Aperture Radar Doppler Tomography Reveals Details of Undiscovered High-Resolution Internal Structure of the Great Pyramid of Giza</i>, arXiv:2208.00811 (Remote Sensing 14, 5231; retracted 10 August 2026)", "https://arxiv.org/abs/2208.00811"),
    ("Retraction notice, Remote Sensing 18, 2679 (2026)", "https://www.mdpi.com/2072-4292/18/16/2679"),
    ("Retraction Watch, 31 August 2026: journal retracts paper claiming network of corridors inside the Great Pyramid", "https://retractionwatch.com/2026/08/31/retraction-corridors-inside-great-pyramid-giza/"),
    ("Seyfzadeh, <i>Replication and Verification: Biondi Protocol</i>, GitHub, 17 September 2026 (protocol v1.5, K = 50 sub-apertures, B_shift = 88 Hz, 32 × 32 px u-PCC, W = 25 steering fit)", "https://github.com/BiondiProtocol/Replication-and-Verification-Biondi-Protocol"),
    ("Biondi, <i>Scanning Inside Volcanoes with Synthetic Aperture Radar Echography Tomographic Doppler Imaging</i>, Remote Sensing 14, 3828 (2022)", "https://www.mdpi.com/2072-4292/14/15/3828"),
    ("Guizar-Sicairos, Thurman & Fienup 2008, <i>Efficient subpixel image registration algorithms</i>, Optics Letters 33, 156 (the upsampled-DFT registration used for the patch shifts)", "https://doi.org/10.1364/OL.33.000156"),
    ("Peterson 1993, <i>Observations and modeling of seismic background noise</i>, USGS Open-File Report 93-322 (ambient ground-motion levels)", "https://pubs.usgs.gov/publication/ofr93322"),
    ("ICEYE Dwell imaging modes (25 s spotlight collection; CSI and SAR video products)", "https://www.iceye.com/sar-data/imaging-modes/dwell"),
]


def esc(s):
    return html.escape(str(s), quote=True)


def fmt(v):
    if isinstance(v, float):
        if v == 0:
            return "0"
        if abs(v) >= 1000 or abs(v) < 0.001:
            return f"{v:.3g}"
        return f"{v:.3f}".rstrip("0").rstrip(".")
    return str(v)


def metrics_table(m):
    rows = []
    for k, v in m.items():
        if isinstance(v, (int, float, str)) or v is None:
            rows.append(f"<tr><td>{esc(k)}</td><td class='num'>{esc(fmt(v))}</td></tr>")
        elif isinstance(v, list) and v and isinstance(v[0], dict):
            keys = list(v[0].keys())
            sub = "<table><tr>" + "".join(f"<th>{esc(x)}</th>" for x in keys) + "</tr>"
            for row in v:
                sub += "<tr>" + "".join(f"<td class='num'>{esc(fmt(row.get(x)))}</td>" for x in keys) + "</tr>"
            sub += "</table>"
            rows.append(f"<tr><td>{esc(k)}</td><td>{sub}</td></tr>")
        elif isinstance(v, dict) and all(isinstance(x, (int, float, str)) or x is None for x in v.values()):
            sub = ", ".join(f"{esc(a)}: {esc(fmt(b))}" for a, b in v.items())
            rows.append(f"<tr><td>{esc(k)}</td><td>{sub}</td></tr>")
    return "<table><tr><th>metric</th><th>value</th></tr>" + "".join(rows) + "</table>" if rows else ""


def load_tests():
    tests = []
    for p in sorted(glob.glob(os.path.join(RESULTS, "*", "summary.json"))):
        with open(p) as f:
            t = json.load(f)
        if t.get("id") != "method":
            tests.append(t)
    return sorted(tests, key=lambda t: t.get("order", 99))


def method_figures():
    p = os.path.join(RESULTS, "method", "summary.json")
    if not os.path.exists(p):
        return ""
    with open(p) as f:
        t = json.load(f)
    out = ""
    for fg in t.get("figures", []):
        rel = copy_fig(t, fg)
        out += f"""<figure class="wide"><div class="frame"><img src="{rel}" alt="" loading="lazy"></div><figcaption>{fg.get('caption','')}</figcaption></figure>\n"""
    return out


JPEG_ABOVE_BYTES = 600_000


def copy_fig(test, fig):
    """Copy a figure into docs/figs; speckle-heavy PNGs above 600 KB become quality-88 JPEGs."""
    src = fig["file"]
    dst_dir = os.path.join(DOCS, "figs", test["id"])
    os.makedirs(dst_dir, exist_ok=True)
    base = os.path.basename(src)
    if os.path.getsize(src) > JPEG_ABOVE_BYTES:
        from PIL import Image
        base = os.path.splitext(base)[0] + ".jpg"
        dst = os.path.join(dst_dir, base)
        im = Image.open(src).convert("RGB")
        im.save(dst, "JPEG", quality=88, optimize=True)
        # remove a stale PNG of the same name
        stale = os.path.join(dst_dir, os.path.basename(src))
        if os.path.exists(stale):
            os.remove(stale)
    else:
        dst = os.path.join(dst_dir, base)
        shutil.copyfile(src, dst)
    return f"figs/{test['id']}/{base}"


def card(t):
    return f"""<article class="card" id="card-{esc(t['id'])}">
  <span class="tag {esc(t.get('tag',''))}">{esc(t.get('tag','')) or 'test'}</span>
  <div class="eyebrow" style="margin-top:12px">{esc(t.get('eyebrow',''))}</div>
  <h3><a href="#{esc(t['id'])}">{esc(t['title'])}</a></h3>
  <p class="question">{esc(t['question'])}</p>
  <dl><dt>Finding</dt><dd>{esc(t['finding'])}</dd><dt>Limitations</dt><dd>{esc(t['limitations'])}</dd></dl>
</article>"""


def planned_card(eyebrow, title, text):
    return f"""<article class="card">
  <span class="tag planned">planned</span>
  <div class="eyebrow" style="margin-top:12px">{esc(eyebrow)}</div>
  <h3>{esc(title)}</h3>
  <p class="question">{esc(text)}</p>
</article>"""


def section(t):
    figs = ""
    for fg in t.get("figures", []):
        rel = copy_fig(t, fg)
        figs += f"""<figure class="wide"><div class="frame"><img src="{rel}" alt="{esc(fg.get('alt', ''))}" loading="lazy"></div><figcaption>{fg.get('caption','')}</figcaption></figure>\n"""
    metrics = metrics_table(t.get("metrics", {}))
    return f"""<section class="part" id="{esc(t['id'])}"><div class="wrap">
  <div class="part-head"><div class="n">{esc(t.get('eyebrow',''))}</div><h2>{esc(t['title'])}</h2></div>
  <div class="prose"><p class="lede">{esc(t['question'])}</p></div>
  <div class="finding"><div class="label">Finding</div><div>{esc(t['finding'])}</div></div>
  <div class="prose"><p><b>What was run.</b> {esc(t.get('method',''))}</p></div>
  {figs}
  <details><summary class="kicker" style="cursor:pointer">Numbers behind this test</summary>{metrics}</details>
  <div class="limits" style="margin-top:22px"><div class="label">Limitations</div><div>{esc(t['limitations'])}</div></div>
</div></section>"""


def method_section(geom_summary, params):
    grow = "".join(f"<tr><td>{esc(k)}</td><td class='num'>{esc(fmt(v))}</td></tr>" for k, v in geom_summary.items())
    prow = "".join(f"<tr><td>{esc(k)}</td><td>{esc(v)}</td></tr>" for k, v in params.items())
    return f"""<section class="part" id="method"><div class="wrap">
  <div class="part-head"><div class="n">Method</div><h2>What the published pipeline actually computes</h2></div>
  <div class="prose">
    <p class="lede">Strip away the springs, magma and “sonic imaging”, and the 2022 paper together with the 2026 replication protocol describe a definite computation on one focused radar image. This site reimplements it and feeds it scenes whose contents are known exactly.</p>
    <ol>
      <li><b>One SLC.</b> A single-look complex image, azimuth by slant range. Its 2-D Fourier transform has a rectangular support; the azimuth (Doppler) axis maps one-to-one onto slow time, platform position and look angle.</li>
      <li><b>Sub-aperture pairs.</b> Two band-pass filters, each half the Doppler band wide and offset from one another by 88 Hz, are swept across the spectrum in K = 50 steps. Each step yields a reference image and an offset image with slightly different look angles.</li>
      <li><b>“Micro-motion.”</b> At every pixel of interest a 32 × 32 patch of the reference image is cross-correlated with the same patch of the offset image and the peak located to a thousandth of a pixel. The 50 resulting shift vectors are the pixel's trajectory. The paper calls this vibration; it is the change of the local interference pattern with look angle.</li>
      <li><b>“Tomography.”</b> A steering matrix borrowed from multi-baseline SAR tomography, exp(j K<sub>z</sub> z), is applied along the trajectory with K<sub>z</sub> = 4π B<sub>⊥</sub> / (λ<sub>s</sub> r sin θ), where B<sub>⊥</sub> is the platform's <i>along-track</i> offset at each sub-aperture and λ<sub>s</sub> a declared “sound wavelength” (0.48 m). The paper projects the complex trajectory on the steering vectors; the replication protocol fits cos and sin of K<sub>z</sub>z over 25-sample windows and keeps the best adjusted R². Either way the depth axis is the trajectory's oscillation rate, relabelled.</li>
    </ol>
    <p>Substituting the definitions gives the result that drives everything below: two scatterers a distance Δx apart in azimuth produce a trajectory that oscillates at a rate proportional to Δx, and the focusing step places it at</p>
    <p style="text-align:center"><b>z = Δx · λ<sub>s</sub> sin θ / λ ≈ 8.9 Δx</b> (for X-band at 35° incidence)</p>
    <p>with the constant a free choice. Depth is surface geometry in disguise.</p>
  </div>
  {method_figures()}
  <h3>Simulation and pipeline settings</h3>
  <div class="figrow"><div><table><tr><th>geometry</th><th></th></tr>{grow}</table></div><div><table><tr><th>pipeline</th><th></th></tr>{prow}</table></div></div>
  <div class="prose"><p>The simulator synthesises the focused image in the spectral domain from explicit point scatterers, with layover, occlusion, aspect-dependent “flash” from block facets, fully developed speckle clutter, and optional line-of-sight vibration applied as a phase history. Nothing in the simulator knows about depth. Code: <a href="{REPO_URL}">{REPO_URL}</a>.</p></div>
</div></section>"""


def build():
    tests = load_tests()
    os.makedirs(os.path.join(DOCS, "figs"), exist_ok=True)
    import sys
    sys.path.insert(0, ROOT)
    from sarsim import Geometry
    geom_summary = {k: (float(v) if not isinstance(v, str) else v) for k, v in Geometry().summary().items()}
    params = {"sub-aperture pairs K": "50", "sub-aperture width": "50% of the processed band", "reference/offset shift": "88 Hz (snapped to 4 FFT bins)",
              "patch": "32 × 32 px complex cross-correlation, unnormalised", "sub-pixel refinement": "upsampled DFT to 1/1000 px + parabolic peak",
              "common-mode": "global median trajectory removed", "declared wavelength λs": "0.48 m (paper)", "depth focus": "paper: |a(z)ᴴY|²; replication: max-over-windows adjusted R², W = 25",
              "depth grid": "160 steps to 97% of the Nyquist depth"}
    cards = "\n".join(card(t) for t in tests) + "\n" + "\n".join(planned_card(*p) for p in PLANNED)
    sections = "\n".join(section(t) for t in tests)
    toc = "".join(f'<li><a href="#{esc(t["id"])}">{esc(t.get("eyebrow","").split("·")[-1].strip())}</a></li>' for t in tests)
    sources = "".join(f'<li><a href="{esc(u)}">{s}</a></li>' for s, u in SOURCES)
    updated = time.strftime("%-d %B %Y")
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(TITLE)}</title>
<meta name="description" content="{esc(DECK)}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,300;0,6..72,400;0,6..72,500;0,6..72,600;1,6..72,400;1,6..72,500&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<link rel="stylesheet" href="style.css">
<meta property="og:type" content="website">
<meta property="og:title" content="{esc(TITLE)}">
<meta property="og:description" content="{esc(DECK)}">
<meta property="og:url" content="{SITE_URL}">
<meta property="og:image" content="{SITE_URL}figs/t01_static_pyramid/t01_slc.jpg">
</head>
<body>
<a class="skip-link" href="#main">Skip to content</a>
<header class="masthead"><div class="wrap">
  <div class="eyebrow">Radar, pyramids &amp; a method under test</div>
  <h1 style="margin-top:14px">{esc(TITLE)}</h1>
  <p class="deck">{esc(DECK)}</p>
  <div class="verdict"><b>Where this stands.</b> The 2022 paper was retracted in August 2026 for “serious methodological flaws”. The tests here go one step further than the retraction: they show <i>what</i> the method computes and reproduce its characteristic output, layered bands and repeated units at depth, from a scene that contains nothing and never moves. The next phase runs the same code on a real 25-second ICEYE scene of Giza.</div>
  <div class="meta"><span><b>Updated</b> {esc(updated)}</span><span><a href="{REPO_URL}">Code, experiments &amp; results</a></span><span><b>Phase</b> simulation complete · real data next</span></div>
</div></header>
<nav class="toc" aria-label="On this page"><div class="wrap"><ol>
  <li><a href="#tests">Tests</a></li><li><a href="#method">Method</a></li>{toc}<li><a href="#roadmap">Roadmap</a></li><li><a href="#sources">Sources</a></li>
</ol></div></nav>
<main id="main">
<section class="part" id="tests"><div class="wrap">
  <div class="part-head"><div class="n">Tests</div><h2>Six tests, one verdict</h2></div>
  <div class="prose"><p>Each test keeps its question, finding and limitations together. Tests are numbered in the order they were run; later ones will be added as the real-data phase proceeds and earlier ones stay as run.</p></div>
  <div class="cards">
{cards}
  </div>
</div></section>
{method_section(geom_summary, params)}
{sections}
<section class="part" id="roadmap"><div class="wrap">
  <div class="part-head"><div class="n">Roadmap</div><h2>What comes next</h2></div>
  <div class="prose">
    <p>The simulation phase settles what the method is. The real-data phase asks whether the published Giza results behave like the simulation, using an ICEYE X33 Dwell Fine scene acquired on 27 August 2025 (25 s illumination, 600 MHz bandwidth, SLC of 9.8 GB, with ICEYE's own colour sub-aperture and SAR-video products as a zero-code check of how strongly the pyramid's scattering changes with look angle).</p>
    <ul>
      <li>Same pipeline, same statistics, on Khafre, empty plateau, Sphinx and Cairo housing (T7).</li>
      <li>Positive control against the surveyed chambers of Khufu (T8).</li>
      <li>Two headings, one ground: do the features move as T3 predicts (T9)?</li>
      <li>A discriminator trained to tell real from simulated tomograms (T10).</li>
    </ul>
    <p>Anyone can rerun everything: <code>pip install -r requirements.txt</code>, then <code>python experiments/t01_static_pyramid.py</code> and so on, then <code>python build_site.py</code>.</p>
  </div>
</div></section>
<section class="part" id="sources"><div class="wrap">
  <div class="part-head"><div class="n">Sources</div><h2>Papers, notices and code</h2></div>
  <div class="prose"><ol>{sources}</ol></div>
</div></section>
</main>
<footer><div class="wrap">Built from <code>results/*/summary.json</code> by <code>build_site.py</code>. Figures are generated by the experiment scripts; none are edited by hand. Simulation and analysis code by Eugene Hong with Claude Fable 5.1.</div></footer>
</body>
</html>
"""
    with open(os.path.join(DOCS, "index.html"), "w") as f:
        f.write(page)
    print(f"wrote docs/index.html with {len(tests)} tests")


if __name__ == "__main__":
    build()
