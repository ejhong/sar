"""The real-acquisition chapter of the site, built from results/r0*/summary.json.

Simulation established what the published method computes. This chapter reports what it does
on the two ICEYE dwell acquisitions we hold, and is written to lead the page: every claim here
rests on measurements of real radar data rather than on a scene we built.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
REPO_URL = "https://github.com/ejhong/sar"
ORDER = ["r01_giza_controls", "r02_depth_axis", "r03_split_dwell",
         "r04_velocity_floor", "r05_learned_null", "r06_second_site", "r07_injected_motion"]


def esc(value):
    return html.escape(str(value), quote=True)


def load_reports():
    reports = {}
    for name in ORDER:
        path = RESULTS / name / "summary.json"
        if path.is_file():
            reports[name] = json.loads(path.read_text())
    return reports


def acquisition_table(acq):
    rows = [
        ("Satellite and mode", f"{acq['satellite']} · {acq['product_type']} · {acq['polarization']}"),
        ("Collected", f"{acq['collection_start'][:19].replace('T', ' ')} UTC"),
        ("Image size", f"{acq['shape'][0]:,} by {acq['shape'][1]:,} complex samples"),
        ("Collection duration", f"{acq['collection_duration_s']:.2f} s"),
        ("Zero-Doppler image span", f"{acq['image_span_s']:.3f} s"),
        ("Processed azimuth bandwidth", f"{acq['bandwidth_hz']:,.0f} Hz"),
        ("Wavelength and incidence", f"{acq['wavelength_m'] * 100:.2f} cm at {acq['incidence_center_deg']:.1f}°"),
        ("Sample spacing", f"{acq['azimuth_spacing_m']:.3f} m azimuth, {acq['range_spacing_m']:.3f} m slant range"),
    ]
    body = "".join(f"<tr><th scope='row'>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in rows)
    return ("<div class='table-scroll' tabindex='0' role='region' aria-label='Acquisition facts'>"
            f"<table class='metrics'><tbody>{body}</tbody></table></div>")


def strip(reports):
    """Three headline numbers, each traceable to one experiment."""
    items = []
    r02 = reports.get("r02_depth_axis")
    if r02:
        depth = r02["metrics"]["nyquist_full_aperture"][1]["nyquist_full_aperture_m"]
        items.append(("Depth axis", f"{depth:.1f} m",
                      "How deep the steering basis reaches before it repeats, using the whole aperture and the "
                      "product's own state vectors. Published depths are hundreds of times larger."))
    r01 = reports.get("r01_giza_controls")
    if r01:
        ks = r01["metrics"]["ks_statistic_khafre_vs_controls"]
        items.append(("Khafre vs desert", f"{ks:.3f}",
                      "Kolmogorov-Smirnov distance between the structure scores on Khafre and on open plateau, "
                      "in the same acquisition. Zero would mean the two are indistinguishable."))
    r04 = reports.get("r04_velocity_floor")
    if r04:
        floor = r04["metrics"]["best_floor_um_s"]
        items.append(("Velocity floor", f"{floor:,.0f} <small>µm/s</small>",
                      "Measured line-of-sight velocity noise on the brightest real targets. Ambient ground motion "
                      "in the microseism band is 0.1 to 10 µm/s."))
    if not items:
        return ""
    cells = "".join(
        f"<article><span class='eyebrow'>{esc(label)}</span><strong>{value}</strong><p>{esc(note)}</p></article>"
        for label, value, note in items)
    return f"<div class='result-strip'>{cells}</div>"


def block(report, figure):
    key = report["id"]
    figs = report.get("figures", [])
    lead = "".join(figure(key, f) for f in figs[:1])
    rest = figs[1:]
    more = ""
    if rest:
        more = (f"<details class='more-figures'><summary>{len(rest)} further figure"
                f"{'s' if len(rest) > 1 else ''} <span aria-hidden='true'>+</span></summary>"
                f"<div class='additional-figures'>" + "".join(figure(key, f) for f in rest) + "</div></details>")
    metrics = json.dumps(report.get("metrics", {}), indent=2, allow_nan=False)
    (ROOT / "docs" / "data").mkdir(parents=True, exist_ok=True)
    (ROOT / "docs" / "data" / f"{key}.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    return f"""<div class="research-block" id="{esc(key)}">
  <span class="eyebrow">{esc(report['eyebrow'])}</span>
  <h3>{esc(report['title'])}</h3>
  <p class="experiment-question">{esc(report['question'])}</p>
  <p class="finding">{esc(report['finding'])}</p>
  {lead}{more}
  <div class="scope"><strong>What this leaves open</strong><p>{esc(report['limitations'])}</p></div>
  <details class="technical"><summary>Method &amp; numbers <span aria-hidden="true">+</span></summary>
    <p>{esc(report.get('method', ''))}</p>
    <p class="resource-links"><a href="{REPO_URL}/blob/main/experiments/{esc(key)}.py">Experiment code ↗</a><a href="data/{esc(key)}.json" download>Results JSON ↓</a></p>
  </details>
</div>"""


def real_section(figure, height_shift_px='', height_shift_m='', field_portal=''):
    reports = load_reports()
    if not reports:
        return ("<section class='section roadmap' id='field'><div class='wrap'>"
                "<div class='section-heading'><span class='section-number'>05 / REAL ACQUISITIONS</span>"
                "<h2>Processing in progress.</h2><p>The real-data chapter is being computed.</p></div>"
                "</div></section>")
    acq = None
    if "r01_giza_controls" in reports:
        acq = reports["r01_giza_controls"]["metrics"].get("acquisition")
    blocks = "".join(block(reports[k], figure) for k in ORDER if k in reports)
    facts = acquisition_table(acq) if acq else ""
    return f"""<section class="section roadmap" id="field"><div class="wrap">
  <div class="section-heading"><span class="section-number">05 / WHAT DO THE REAL SCANS SHOW?</span>
    <h2>The same acquisition,<br>the same answer everywhere.</h2>
    <p>Everything above is simulation. This chapter is measurement. We hold two ICEYE Spotlight Dwell Fine
    acquisitions, one of the Giza plateau and one of Sacsayhuamán, and the published processing chain has now been
    run on the Giza scene at four monuments and two stretches of open desert, all within the same twenty-five
    seconds of radar.</p></div>
  <div class="field-status"><span class="evidence-badge">Real-data status</span><p>A dwell product is the most
  favourable case this method can be given: its Doppler axis spans about 24.5 s of real aperture, so sub-apertures
  are separated by seconds rather than milliseconds, and genuine ground motion is at least in principle within
  reach. The results below use the product's own orbit state vectors for the steering geometry rather than the
  scalar baseline approximation the published protocols use.</p></div>
  {facts}
  {strip(reports)}
  {blocks}
  {field_portal}
  <div class="support-heading"><span class="eyebrow">Still to do</span><h3>Score the method against chambers whose geometry is measured.</h3>
  <p>Giza has surveyed shafts and burial chambers. A fair positive control compares a frozen method with those
  surveys and reports the misses as well as the matches. Placement has to be settled first: changing the assumed
  surface elevation by ten metres moves the metadata projection by about <strong>{height_shift_px} range
  pixels</strong>, roughly {height_shift_m} m on flat ground, which is larger than the chambers themselves.</p></div>
"""
