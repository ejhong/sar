"""Build the research notebook from portable experiment reports and figures."""
from datetime import date
from pathlib import Path
import html
import json
import re
import shutil
from string import Template

from PIL import Image
from experiments.reporting import refine_summary

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
REPO_URL = "https://github.com/ejhong/sar"
EXPERIMENTS = {
    "t01_static_pyramid": ("Static scene", "Patterns without buried objects", "t01_transect_az", "t01_maps"),
    "t02_desert_null": ("Null control", "Desert and pyramid scores overlap", "t02_histograms", "t02_maps"),
    "t03_mechanism": ("Mechanism", "Surface spacing predicts the peak", "t03_two_scatterers", "t03_nonlocal"),
    "t04_parameters": ("Parameters", "The depth scale is a processing choice", "t04_lambda", "t04_patch"),
    "t05_motion": ("Positive input", "Motion enters; frequency becomes depth", "t05_sensitivity", "t05_depth_is_frequency"),
    "t06_ordering": ("Selection gates", "The full gate rejects these static pixels", "t06_permutation"),
    "t07_surface_controls": ("Paired controls", "What changes when surface points are added?", "t07_paired_metrics", "t07_controlled_maps"),
}


def esc(value):
    return html.escape(str(value), quote=True)


def plain(value):
    return html.unescape(re.sub(r"<[^>]+>", "", str(value)))


def fmt(value):
    if value is None:
        return "—"
    return f"{value:.4g}" if isinstance(value, float) else str(value)


def resolve_figure(value):
    path = Path(value)
    if path.is_absolute():
        if "results" in path.parts:
            path = ROOT.joinpath(*path.parts[path.parts.index("results"):])
    else:
        path = ROOT / path
    path = path.resolve()
    if not path.is_relative_to(RESULTS.resolve()):
        raise ValueError(f"Figure outside results: {value}")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def copy_figure(group, source):
    source = resolve_figure(source)
    destination = DOCS / "figs" / group
    destination.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as im:
        width, height = im.size
        if source.suffix.lower() == ".png" and source.stat().st_size > 600_000:
            name = source.stem + ".jpg"
            im.convert("RGB").save(destination / name, "JPEG", quality=90, optimize=True)
        else:
            name = source.name
            shutil.copyfile(source, destination / name)
    return f"figs/{group}/{name}", width, height


def figure(group, item, eager=False):
    src, width, height = copy_figure(group, item["file"])
    caption = item.get("caption", "")
    alt = item.get("alt") or plain(caption).split(". ")[0].rstrip(".")
    priority = 'loading="eager" fetchpriority="high"' if eager else 'loading="lazy"'
    return f'''<figure class="research-figure"><a class="figure-link" href="{src}" data-zoom aria-label="Enlarge figure: {esc(alt)}"><img src="{src}" width="{width}" height="{height}" alt="{esc(alt)}" {priority} decoding="async"><span class="zoom-hint" aria-hidden="true">View full size ↗</span></a><figcaption>{caption}</figcaption></figure>'''


def json_report(path):
    return json.loads(path.read_text())


def copy_data(name, value):
    destination = DOCS / "data"
    destination.mkdir(exist_ok=True)
    (destination / name).write_text(json.dumps(value, indent=2, allow_nan=False))
    return f"data/{name}"


def metrics_table(values):
    rows = []
    for key, value in values.items():
        label = key.replace("_", " ")
        if isinstance(value, dict):
            rendered = metrics_table(value)
        elif isinstance(value, list) and value and isinstance(value[0], dict):
            keys = list(value[0])
            heads = "".join(f"<th scope='col'>{esc(k.replace('_', ' '))}</th>" for k in keys)
            cells = "".join("<tr>" + "".join(f"<td>{esc(fmt(row.get(k)))}</td>" for k in keys) + "</tr>" for row in value)
            rendered = f"<div class='table-scroll' tabindex='0' role='region' aria-label='{esc(label)}'><table><thead><tr>{heads}</tr></thead><tbody>{cells}</tbody></table></div>"
        else:
            rendered = esc(fmt(value))
        rows.append(f"<tr><th scope='row'>{esc(label)}</th><td>{rendered}</td></tr>")
    return "<div class='table-scroll' tabindex='0' role='region' aria-label='Experiment metrics'><table class='metrics'><tbody>" + "".join(rows) + "</tbody></table></div>"


def experiment_section(report):
    key = report["id"]
    label, short, *lead = EXPERIMENTS.get(key, ("Experiment", report["title"]))
    figures = report.get("figures", [])
    first = [f for name in lead for f in figures if Path(f["file"]).stem == name]
    rest = [f for f in figures if f not in first]
    if not first and figures:
        first, rest = figures[:1], figures[1:]
    body = "".join(figure(key, f) for f in first)
    if rest:
        body += f"<details class='more-figures'><summary>{len(rest)} additional figures &amp; comparisons <span aria-hidden='true'>+</span></summary><div class='additional-figures'>" + "".join(figure(key, f) for f in rest) + "</div></details>"
    data = copy_data(key + ".json", report)
    extra_downloads = ""
    for filename in report.get("downloads", []):
        destination = key + "_" + filename
        shutil.copyfile(RESULTS / key / filename, DOCS / "data" / destination)
        extra_downloads += f'<a href="data/{destination}" download>{esc(Path(filename).stem.replace("_", " ").capitalize())} ↓</a>'
    source_note = ""
    if key == "t07_surface_controls":
        source_note = "<p class='source-note'>The previous eight-point rendering is retained in the <a href='#t07_wells'>illustrative appendix</a>. It is excluded from the numbered evidence and does not reproduce the Khafre observations.</p>"
        body = gate_table(report["metrics"]["runs"]) + body
    return f'''<details class="experiment" id="{key}"><summary><span class="experiment-number">{report['order']:02d}</span><span class="experiment-heading"><span class="eyebrow">{esc(label)}</span><span class="experiment-title">{esc(short)}</span></span><span class="expand" aria-hidden="true">+</span></summary><div class="experiment-body"><h3>{esc(report['title'])}</h3><p class="experiment-question">{esc(report['question'])}</p><p class="finding">{esc(report['finding'])}</p><div class="scope"><strong>What this leaves open</strong><p>{esc(report['limitations'])}</p></div>{body}{source_note}<details class="technical"><summary>Method, numbers &amp; reproducibility <span aria-hidden="true">+</span></summary><p>{esc(report.get('method', ''))}</p>{metrics_table(report.get('metrics', {}))}<p class="resource-links"><a href="{REPO_URL}/blob/main/experiments/{key}.py">Experiment code ↗</a><a href="{data}" download>Results JSON ↓</a>{extra_downloads}</p></details></div></details>'''


def gate_table(runs):
    rows = []
    for run in runs:
        name = ("Unchanged pyramid" if run["layout"] == "none" else
                "Original layout" if run["layout"] == "original" else
                "Random layout " + run["layout"].split("_")[-1])
        rows.append(f"<tr><th scope='row'>{esc(name)}</th><td>{fmt(run['strength'])}</td><td>{fmt(run['median_point_power_ratio'])}</td><td>{run['shape_gate_count']:,}</td><td>{run['full_gate_count']:,}</td></tr>")
    return f'''<details class="technical"><summary>Check the power response and both acceptance gates <span aria-hidden="true">+</span></summary><p>Every scene uses the same {runs[0]['sampled_pixels']:,} grid locations. Counts are dependent pixels, not independent detections. The full gate adds the 0.005 px minor-axis floor to the shape requirements. Power ratios compare fixed point regions with those same locations in the unchanged scene; no ratio is assigned to the baseline itself.</p><div class="table-scroll" tabindex="0" role="region" aria-label="Paired controls and acceptance gates"><table><thead><tr><th scope="col">Layout</th><th scope="col">Strength</th><th scope="col">Median power ratio</th><th scope="col">Shape-only passes</th><th scope="col">Full-gate passes</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></details>'''


def appendix_section(report):
    key = report["id"]
    data = copy_data(key + ".json", report)
    figures = report.get("figures", [])
    leading = [f for f in figures if Path(f["file"]).stem in ("t07_plan", "t07_isosurface")]
    rest = [f for f in figures if f not in leading]
    body = "".join(figure(key, f) for f in leading)
    if rest:
        body += "<details class='more-figures'><summary>Earlier comparisons and additional renderings <span aria-hidden='true'>+</span></summary><div class='additional-figures'>" + "".join(figure(key, f) for f in rest) + "</div></details>"
    return f'''<section class="section appendix" id="appendix"><div class="wrap"><div class="section-heading"><span class="section-number">APPENDIX / ILLUSTRATION</span><h2>A constructed example,<br>with a narrower role.</h2><p>The former T7 is retained for transparency. Eight bright surface points and presentation settings were deliberately chosen. The eight-point arrangement and the roughly 630 m displayed extent are inputs to this illustration, not independent predictions about Khafre. See <a href="#t07_surface_controls">the paired T7 controls</a> for measured effects before smoothing.</p></div><details class="support-detail" id="{key}"><summary>Illustration: surface points and rendering choices <span aria-hidden="true">+</span></summary><div class="detail-inner"><p class="finding">{esc(report['finding'])}</p><div class="scope"><strong>Illustrative scope</strong><p>{esc(report['limitations'])}</p></div>{body}<p class="source-note">Published comparison images belong to the Khafre Research Project. <a href="{REPO_URL}/blob/main/results/t07_wells/published/CREDITS.md">Image sources &amp; credits</a>. These renderings do not establish how the published features were produced.</p><details class="technical"><summary>Illustration settings and archived numbers <span aria-hidden="true">+</span></summary><p>{esc(report['method'])}</p>{metrics_table(report['metrics'])}<p class="resource-links"><a href="{REPO_URL}/blob/main/experiments/t07_wells.py">Illustration code ↗</a><a href="{data}" download>Illustration results ↓</a></p></details></div></details></div></section>'''


def support_figure(stem, caption, alt):
    return figure("feasibility", {"file": f"results/feasibility/{stem}.png", "caption": caption, "alt": alt})


def controls_section(data):
    phase30 = next(r for r in data["motion"]["rows"] if r["amplitude_um"] == 30)
    null = data["motion"]["rows"][0]
    site = next(r for r in data["ai"]["summary"] if r["experiment"] == "site_identity_only" and r["split"] == "held_out_sites")
    wave = data["waves"]["numerical_check"]
    download = copy_data("feasibility.json", data)
    readout = support_figure("phase_readout", "<b>Two measurements of the same input.</b> A patch receives a uniform phase rotation. Its phase follows the imposed cycle; its translation remains zero. This is an algebraic control, with no SAR focusing or physical scene motion.", "Imposed and recovered phase-equivalent displacement agree, while patch translation stays at zero")
    receiver = support_figure("motion_detection", "<b>An ideal pulse-time receiver.</b> Separate noise trials set the threshold before evaluation. The model contains one coherent scatterer, independent noise and perfect platform compensation; its frequency lies exactly on a Fourier bin.", "Detection rises with imposed vibration amplitude in the ideal phase receiver; pulse intensity stays near false-alarm level")
    waves = support_figure("cavity_wave_response", "<b>A mechanical connection in a controlled model.</b> Difference between cavity and cavity-free surface waveforms under the same prescribed source. These normalized differences are not satellite detection probabilities.", "Surface wave differences for ideal cylindrical cavities centered 60, 120 and 240 metres deep")
    ai = support_figure("ai_controls", "<b>High scores can have the wrong explanation.</b> Left: site recognition supplies a perfect patch-split score. Right: a separately planted feature supports generalization to unseen sites. Error bars show variation across ten synthetic datasets.", "AUC falls from 1 to about 0.51 with whole-site holdouts; a planted transferable signal raises held-out-site AUC")
    return f'''<section class="section support" id="signal"><div class="wrap"><div class="section-heading"><span class="section-number">03 / SIGNAL &amp; INFERENCE</span><h2>A measurable signal is<br>only the first step.</h2><p>Could underground geometry affect surface motion in a way radar can use? These supporting experiments examine separate parts of that chain.</p></div>
<ol class="signal-chain"><li><span>01</span><strong>Underground geometry</strong><small>Cavity, material, source</small></li><li><span>02</span><strong>Surface motion</strong><small>Amplitude and wave pattern</small></li><li><span>03</span><strong>Radar measurement</strong><small>Phase, noise, processing</small></li><li><span>04</span><strong>Validated inference</strong><small>Detection, depth, geometry</small></li></ol>
<p class="chain-note">The wave and receiver models below are independent. No calibrated cavity-to-radar simulation or field detection is demonstrated here.</p>
<article class="support-study" id="phase-control"><div class="study-copy"><span class="eyebrow">S1 / What is measured</span><h3>Phase and translation answer different questions.</h3><p>A correlation peak tells us how far a patch shifts. A phase measurement can respond while that shift remains zero. This explains why T5's response cannot be treated as a universal limit on radar vibration measurements.</p><p>The control applies a phase rotation equivalent to a ±30 μm displacement. Both estimators receive exactly the same patch.</p></div><div>{readout}</div></article>
<details class="support-detail"><summary>How sensitive is an ideal phase receiver? <span aria-hidden="true">+</span></summary><div class="detail-inner"><p>With a 3.1 cm wavelength, 8 s record, 512 pulses/s and 20 dB per-pulse signal-to-noise ratio, the phase detector finds the 30 μm input in <strong>{phase30['phase_detection_rate']:.1%}</strong> of 800 trials. Independent zero-signal trials produce <strong>{null['phase_detection_rate']:.2%}</strong> detections against a nominal 1% threshold. Distributed clutter, autofocus, atmospheric effects and real acquisition geometry are absent.</p>{receiver}<p>Airborne vibration measurements provide a precedent for controlled target vibrometry; they do not validate satellite cavity mapping. <a href="https://epublications.marquette.edu/electric_fac/571/">Wang et al., 2012 ↗</a></p></div></details>
<article class="support-study" id="wave-control"><div class="study-copy"><span class="eyebrow">S2 / The physical connection</span><h3>A cavity can change the surface wavefield.</h3><p>A simple shear-wave model gives different surface responses with and without a buried cavity. The location and timing of the difference change with cavity depth.</p><p>The model represents an infinitely long cylindrical cavity with radius 20 m. It has a known active source, uniform rock and no intrinsic attenuation. Its displacement amplitude is uncalibrated.</p></div><div>{waves}</div></article>
<p class="study-footnote">Numerical check: halving grid spacing gives a cavity-perturbation correlation of {wave['cavity_difference_correlation_coarse_fine']:.3f}; the relative change at 120 m shifts from 24.5% to {wave['fine_grid_relative_waveform_L2_change_depth120']:.1%}. This supports the qualitative model effect, not a depth limit.</p>
<article class="support-study" id="learning-control"><div class="study-copy"><span class="eyebrow">S3 / Learning from the right evidence</span><h3>A model can recognize a site without detecting a cavity.</h3><p>In the site-only control, randomly split patches score AUC 1.00. Holding out whole sites reduces it to <strong>{site['mean_auc']:.2f}</strong>, near chance. No transferable cavity signal was supplied.</p><p>A second control deliberately supplies a transferable feature and gets better held-out performance. Test for added information beyond location, geology and visible surface clues.</p></div><div>{ai}</div></article>
<p class="source-note">Supporting models adapted and rerun from the supplied SAR-Voids-Quick-Experiments package. All data here are synthetic. <a href="{REPO_URL}/blob/main/experiments/feasibility.py">Code ↗</a> · <a href="{download}" download>Results &amp; versions ↓</a></p></div></section>'''


def robustness_section(data):
    src = figure("robustness", {"file": "results/robustness/robustness.png", "caption": "<b>Stress-testing the examples.</b> Left: ten phase realizations at each of four surface separations. Right: the ordering statistic tested on planted noisy sinusoids and independent noise. These are checks of the model and statistic, not field validation.", "alt": "Peak depths across forty two-scatterer runs and permutation rejection rates for ordered and noise controls"})
    ordering = data["ordering_control"]
    url = copy_data("robustness.json", data)
    return f'''<details class="robustness"><summary><span><span class="eyebrow">Additional checks</span><strong>Do the examples survive simple stress tests?</strong></span><span class="expand" aria-hidden="true">+</span></summary><div class="detail-inner"><p>Across {data['runs']} two-scatterer runs, median peak error relative to the spacing prediction is {data['median_peak_error_m']:.2f} m; the largest is {data['max_peak_error_m']:.2f} m. End-to-end reruns at three wavelengths verify the predicted rescaling.</p><p>The permutation test flags {ordering['ordered_rejection_fraction']:.1%} of planted ordered trajectories and {ordering['noise_rejection_fraction']:.1%} of independent-noise trajectories, at a nominal 5% level. It is sensitive to this positive control, which is generated in trajectory space rather than through a radar simulation.</p>{src}<p class="resource-links"><a href="{REPO_URL}/blob/main/experiments/robustness.py">Code ↗</a><a href="{url}" download>Results ↓</a></p></div></details>'''


def build():
    DOCS.mkdir(exist_ok=True)
    reports = [refine_summary(json_report(path)) for path in RESULTS.glob("t*/summary.json")]
    tests = sorted([report for report in reports if report.get("role") != "illustration"], key=lambda x: x["order"])
    illustrations = [report for report in reports if report.get("role") == "illustration"]
    controls = json_report(RESULTS / "feasibility" / "metrics.json")
    robustness = json_report(RESULTS / "robustness" / "metrics.json")
    overview = RESULTS / "overview" / "overview.png"
    if overview.is_file():
        hero_item = {"file": str(overview), "caption": "<b>Known input. Computed output.</b> A stationary simulated pyramid, with no buried objects, and the normalized depth spectrum along one transect. Brightness is relative spectral power, not confidence in an underground structure.", "alt": "Stationary simulated pyramid and its depth-like spectral output along a marked transect"}
        hero = figure("overview", hero_item, eager=True)
        hero_src = copy_figure("overview", str(overview))[0]
    else:
        hero = figure("t03_mechanism", tests[2]["figures"][0], eager=True)
        hero_src = copy_figure("t03_mechanism", tests[2]["figures"][0]["file"])[0]
    sections = "".join(experiment_section(t) for t in tests)
    for filename in list((RESULTS / "feasibility").glob("*.csv")) + list((RESULTS / "robustness").glob("*.csv")):
        shutil.copyfile(filename, DOCS / "data" / filename.name)
    latest = max(t.get("date", "2026-09-19") for t in tests + [controls, robustness])
    updated = date.fromisoformat(latest).strftime("%d %B %Y").lstrip("0")
    versions = " · ".join(f"{esc(k)} {esc(v)}" for k, v in controls.get("versions", {}).items())
    downloads = "".join(f"<a href='data/{t['id']}.json' download>T{t['order']} results ↓</a>" for t in tests)
    template = Template((ROOT / "site" / "page.html").read_text())
    page = template.substitute(updated=updated, hero=hero, hero_src=hero_src, count=len(tests),
                               experiments=sections, controls=controls_section(controls),
                               robustness=robustness_section(robustness), downloads=downloads, versions=versions,
                               appendix="".join(appendix_section(report) for report in illustrations))
    (DOCS / "index.html").write_text(page)
    print(f"Built docs/index.html: {len(tests)} scene experiments, supporting controls, robustness checks")


if __name__ == "__main__":
    build()
