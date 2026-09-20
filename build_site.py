"""Build the research notebook from portable experiment reports and figures."""
from datetime import date
from pathlib import Path
import html
import hashlib
import json
import re
import shutil
from string import Template

from PIL import Image
from experiments.reporting import refine_summary
from site_content import validation_content

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
        source_note = "<p class='source-note'>The previous eight-point rendering is <a href='#t07_wells'>archived separately</a>. These paired power measurements do not reproduce the Khafre observations.</p>"
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
    validation = json_report(RESULTS / "validation" / "metrics.json")
    if validation['design'] != json_report(RESULTS / 'validation' / 'design.json'):
        raise ValueError('Incomplete validation run: design and results disagree; rerun experiments/validation.py')
    for path, digest in validation['source_sha256'].items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Physical benchmark source changed: {path}; rerun experiments/validation.py')
    content = validation_content(validation, figure)
    hero_src = copy_figure("validation", "results/validation/noise_transition.png")[0]
    sections = "".join(experiment_section(t) for t in tests)
    for filename in list((RESULTS / "feasibility").glob("*.csv")) + list((RESULTS / "robustness").glob("*.csv")):
        shutil.copyfile(filename, DOCS / "data" / filename.name)
    for path in (RESULTS / "validation").glob("*"):
        if path.suffix in (".json", ".csv") and not path.name.endswith(".partial.json"):
            shutil.copyfile(path, DOCS / "data" / ("validation_" + path.name))
    for name in ("method_audit.json", "real_data_plan.json"):
        shutil.copyfile(ROOT / "site" / name, DOCS / "data" / name)
    copy_data("feasibility.json", controls)
    for report in illustrations:
        copy_data(report['id'] + '.json', report)
    latest = max(t.get("date", "2026-09-19") for t in tests + [controls, robustness, validation])
    updated = date.fromisoformat(latest).strftime("%d %B %Y").lstrip("0")
    versions = " · ".join(f"{esc(k)} {esc(v)}" for k, v in validation.get("versions", {}).items())
    downloads = "".join(f"<a href='data/{t['id']}.json' download>T{t['order']} results ↓</a>" for t in tests)
    template = Template((ROOT / "site" / "page.html").read_text())
    page = template.substitute(updated=updated, hero_src=hero_src, count=len(tests),
                               experiments=sections, robustness=robustness_section(robustness),
                               downloads=downloads, versions=versions, **content)
    (DOCS / "index.html").write_text(page)
    print(f"Built docs/index.html: coupled physical benchmarks, method audit, {len(tests)} supporting experiments")


if __name__ == "__main__":
    build()
