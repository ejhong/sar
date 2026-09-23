"""Build docs/index.html as a dashboard.

Every result in results/*/summary.json becomes one compact card: statement heading, the finding
split into bullets, its headline numbers, a figure and collapsed method detail. The three
dimensional viewer leads the page. Presentation rules are in STYLE.md.
"""
from datetime import date
from pathlib import Path
import html
import json
import re
import shutil

from PIL import Image

from fieldwork.publishing import build_field_site

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
REPO = "https://github.com/ejhong/sar"

ACQ = "r00_acquisitions"
# Ordered as an argument, not by number: the summary, then the positive control, then what the
# method does and does not distinguish, then why, then whether our own code is to blame.
REAL = ["r14_budget", "r13_known_voids", "r01_giza_controls", "r06_second_site", "r08_sections",
        "r02_depth_axis", "r09_coherence", "r10_scatterer_coherence", "r11_array",
        "r03_split_dwell", "r04_velocity_floor", "r07_injected_motion", "r05_learned_null",
        "r12_independent_implementation"]
SIM = ["t01_static_pyramid", "t02_desert_null", "t03_mechanism", "t04_parameters",
       "t05_motion", "t06_ordering", "t07_surface_controls", "t07_wells"]

# Which numbers each card puts on its footer strip. One place to maintain.
HEADLINE = {
    "r14_budget": [],
    "r13_known_voids": [("tombs", "tomb_band_share_mean", "{:.4f}"),
                        ("plateau", "control_band_share_mean", "{:.4f}"),
                        ("matched effect", "standardised_difference_brightness_matched", "{:+.3f}")],
    "r10_scatterer_coherence": [],
    "r11_array": [("real", "median_peak_correlation_real", "{:.3f}"),
                  ("injected wave", "median_peak_correlation_control", "{:.2f}")],
    "r12_independent_implementation": [],
    "r01_giza_controls": [("monuments", "monument_mean", "{:.4f}"), ("open desert", "control_mean", "{:.4f}"),
                          ("standardised gap", "cohens_d", "{:+.3f}")],
    "r02_depth_axis": [("folds back at", None, None), ("aperture for 648 m", "percent_of_aperture_for_648m_K50", "{:.1f}%")],
    "r09_coherence": [("half coherence", None, None), ("tenth coherence", None, None)],
    "r03_split_dwell": [("same pixel", None, None), ("unrelated pixels", None, None)],
    "r08_sections": [("repeats every", None, None), ("correlation there", None, None)],
    "r04_velocity_floor": [("velocity floor", "best_floor_um_s", "{:,.0f} µm/s")],
    "r07_injected_motion": [("static floor", "static_floor_um_s", "{:,.0f} µm/s"),
                            ("threshold", "detection_threshold_velocity_um_s", "{:,.0f} µm/s")],
    "r05_learned_null": [("depth profiles", None, None), ("brightness", None, None)],
    "r06_second_site": [("cut stone", "monument_mean", "{:.3f}"), ("hillside", "control_mean", "{:.3f}"),
                        ("city", "urban_mean", "{:.3f}")],
}


def esc(v):
    return html.escape(str(v), quote=True)


def plain(v):
    return html.unescape(re.sub(r"<[^>]+>", "", str(v)))


def report(name):
    path = RESULTS / name / "summary.json"
    return json.loads(path.read_text()) if path.is_file() else None


def copy_figure(group, source):
    src = Path(source)
    if not src.is_absolute():
        src = ROOT / src
    src = src.resolve()
    if not src.is_relative_to(RESULTS.resolve()) or not src.is_file():
        raise FileNotFoundError(source)
    dest = DOCS / "figs" / group
    dest.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as im:
        if src.suffix.lower() == ".png" and src.stat().st_size > 500_000:
            name = src.stem + ".jpg"
            im.convert("RGB").save(dest / name, "JPEG", quality=88, optimize=True)
        else:
            name = src.name
            shutil.copyfile(src, dest / name)
    return f"figs/{group}/{name}"


def bullets(text, limit=5):
    parts = [p.strip() for p in re.split(r"(?<=[.!?]) +", str(text)) if p.strip()]
    if not parts:
        return "", []
    lead = parts.pop(0)
    # A one-word verdict is not a paragraph; fold it onto the sentence it answers.
    while parts and len(lead) < 34:
        lead = lead + " " + parts.pop(0)
    return lead, parts[:limit]


def numbers(rep):
    """The footer strip. Entries with a None path are filled from the finding text."""
    spec = HEADLINE.get(rep["id"], [])
    m = rep.get("metrics", {})
    out = []
    for label, key, fmt in spec:
        if key and key in m and m[key] is not None:
            out.append((label, fmt.format(m[key])))
    if rep["id"] == "r14_budget":
        out = [(l["name"].lower(), f"{l['shortfall']:,.0f}×") for l in m["legs"]]
    if rep["id"] == "r10_scatterer_coherence":
        g = m["grades"][-1]
        out = [("brightest 0.1%", f"{g['half_coherence_s']:.2f} s"),
               ("darkest half", f"{m['grades'][0]['half_coherence_s']:.2f} s")]
    if rep["id"] == "r12_independent_implementation":
        ok = [r for r in m["runs"] if r.get("returncode") == 0]
        if ok:
            out = [("distinct vector values", str(ok[0]["vector_values_unique"])),
                   ("its depth grid", "100 m"),
                   ("its real limit", f"{ok[-1]['recurrence_depth_m']:.0f} m")]
    if rep["id"] == "r02_depth_axis":
        out.insert(0, ("folds back at", f"{m['nyquist_full_aperture'][1]['nyquist_full_aperture_m']:.1f} m"))
    if rep["id"] == "r09_coherence":
        p = m["patches"][0]
        out = [("half coherence", f"{p['half_coherence_s']:.2f} s"), ("a tenth", f"{p['tenth_coherence_s']:.2f} s")]
    if rep["id"] == "r03_split_dwell":
        r = m["runs"][0]
        out = [("same pixel", f"{r['profile_correlation_median']:.3f}"),
               ("unrelated", f"{r['shuffled_correlation_median']:.3f}")]
    if rep["id"] == "r08_sections":
        p = m["panels"][0]
        out = [("repeats every", f"{p['repetition_period_m']:.1f} m"),
               ("correlation", f"{p['repetition_correlation']:.4f}")]
    if rep["id"] == "r05_learned_null":
        a = m["auc"]
        out = [("depth profiles", f"{a['depth profile only']['mean']:.2f}"),
               ("brightness", f"{a['amplitude only']['mean']:.2f}"),
               ("chance", "0.50")]
    if rep["id"] == "r07_injected_motion" and m.get("detection_threshold_velocity_um_s") is None:
        out = [("static floor", f"{m['static_floor_um_s']:,.0f} µm/s")]
    return out


def card(rep, thumb=True):
    lead, rest = bullets(rep["finding"])
    figs = rep.get("figures", [])
    img = ""
    if thumb and figs:
        src = copy_figure(rep["id"], figs[0]["file"])
        alt = plain(figs[0].get("caption", ""))[:110]
        img = (f'<figure><a href="{src}"><img src="{src}" alt="{esc(alt)}" loading="lazy" '
               f'decoding="async"></a></figure>')
    gallery = ""
    extras = figs[1:] if thumb else figs
    if extras:
        thumbs = "".join(
            f'<a href="{copy_figure(rep["id"], f["file"])}" title="{esc(plain(f.get("caption", ""))[:90])}">'
            f'<img src="{copy_figure(rep["id"], f["file"])}" alt="{esc(plain(f.get("caption", ""))[:70])}" loading="lazy"></a>'
            for f in extras)
        gallery = f'<div class="db-gallery">{thumbs}</div>'
    items = "".join(f"<li>{esc(b)}</li>" for b in rest)
    nums = "".join(f"<div>{esc(k)}<b>{esc(v)}</b></div>" for k, v in numbers(rep))
    (DOCS / "data").mkdir(exist_ok=True)
    (DOCS / "data" / f"{rep['id']}.json").write_text(json.dumps(rep, indent=2, allow_nan=False))
    return f"""<article class="db-card" id="{esc(rep['id'])}">{img}<div class="db-card-body">
<span class="db-tag {esc(rep.get('tag',''))}">{esc(rep['eyebrow'])}</span>
<h3>{esc(rep['title'])}</h3>
<p class="db-find">{esc(lead)}</p>
{f'<ul>{items}</ul>' if items else ''}
{f'<div class="db-num">{nums}</div>' if nums else ''}
<details><summary>Method, limits &amp; data</summary>
<p><b>Method.</b> {esc(rep.get('method',''))}</p>
<p><b>Leaves open.</b> {esc(rep['limitations'])}</p>
{gallery}
<p><a href="{REPO}/blob/main/experiments/{esc(rep['id'])}.py">Code ↗</a> · <a href="data/{esc(rep['id'])}.json" download>Results JSON ↓</a></p>
</details></div></article>"""


def stat(label, value, note, href):
    return (f'<a class="db-stat" href="{href}"><span>{esc(label)}</span>'
            f'<strong>{value}</strong><p>{esc(note)}</p></a>')


def stats(reps):
    m = {k: v["metrics"] for k, v in reps.items()}
    out = []
    if "r14_budget" in m:
        for leg in m["r14_budget"]["legs"]:
            out.append(stat(leg["name"], f"{leg['shortfall']:,.0f}<small>×</small>",
                            f"Short of what the claim needs. Required {leg['need_v']:,.1f} {leg['unit']}, "
                            f"available {leg['have_v']:,.2f}.", "#r14_budget"))
    if "r13_known_voids" in m:
        d = m["r13_known_voids"]["standardised_difference_brightness_matched"]
        out.append(stat("Known voids", f"{d:+.3f}",
                        "Standardised difference between hundreds of surveyed burial shafts and bare plateau, "
                        "at the shafts' own depths.", "#r13_known_voids"))
    if "r11_array" in m:
        a = m["r11_array"]
        out.append(stat("Array test", f"{a['median_peak_correlation_real']:.3f}",
                        f"Peak cross-correlation between array elements on real data. An injected wavefield "
                        f"gives {a['median_peak_correlation_control']:.2f}.", "#r11_array"))
    if "r02_depth_axis" in m and False:
        d = m["r02_depth_axis"]["nyquist_full_aperture"][1]["nyquist_full_aperture_m"]
        out.append(stat("Depth axis", f"{d:.1f}<small> m</small>",
                        "Where the steering basis folds back on the full aperture. It repeats exactly at twice that.",
                        "#r02_depth_axis"))
    if "r01_giza_controls" in m:
        c = m["r01_giza_controls"]["cohens_d"]
        out.append(stat("Monument vs desert", f"{c:+.3f}",
                        "Standardised difference in structure score across 287,280 pixels of one acquisition.",
                        "#r01_giza_controls"))
    if "r09_coherence" in m:
        h = m["r09_coherence"]["patches"][0]["half_coherence_s"]
        out.append(stat("Coherence half-life", f"{h:.2f}<small> s</small>",
                        "Look separation at which the ground stops looking like itself. Microseism needs seconds.",
                        "#r09_coherence"))
    if "r04_velocity_floor" in m:
        f = m["r04_velocity_floor"]["best_floor_um_s"]
        out.append(stat("Velocity floor", f"{f:,.0f}<small> µm/s</small>",
                        "Measured on coherent looks. Ambient ground motion is 0.1 to 10 µm/s.",
                        "#r04_velocity_floor"))
    if "r03_split_dwell" in m:
        r = m["r03_split_dwell"]["runs"][0]
        out.append(stat("Split-dwell agreement", f"{r['profile_correlation_median']:.3f}",
                        f"Correlation between the two halves of one acquisition. Unrelated pixels give {r['shuffled_correlation_median']:.3f}.",
                        "#r03_split_dwell"))
    if "r05_learned_null" in m:
        a = m["r05_learned_null"]["auc"]["depth profile only"]["mean"]
        out.append(stat("Learned null", f"{a:.2f}",
                        "Area under the curve for a classifier given only the depth profiles. Chance is 0.50.",
                        "#r05_learned_null"))
    return "".join(out[:6])


def supporting_cards():
    """Hand-composed cards for the three results that predate the summary.json schema."""
    out = []
    v = json.loads((RESULTS / "validation" / "metrics.json").read_text())
    rows = {r["snr_db"]: r for r in v["noise_sweep"]["rows"] if r["family"] == "matched_assumptions"}
    hi = max(rows), rows[max(rows)]
    lo = min(rows), rows[min(rows)]
    out.append(dict(
        key="validation", tag="mechanism", eyebrow="Steelman · coupled benchmark",
        title="Given a favourable physical model, depth is recoverable",
        lead=(f"A buried cavity drives a wave-equation forward model, the surface motion modulates a radar "
              f"acquisition, and a physical-template inverse recovers the depth."),
        points=[f"At {hi[0]:.0f} dB channel signal-to-noise, {hi[1]['detected_within_20m_rate']:.0%} of trials "
                f"locate the cavity within 20 m, with a {hi[1]['false_alarm_rate']:.0%} false-alarm rate.",
                f"At {lo[0]:.0f} dB that falls to {lo[1]['detected_within_20m_rate']:.0%}.",
                "The forward model never uses the Doppler-to-depth steering formula, so this is an independent "
                "check that depth information can survive in principle, not a vindication of the published method.",
                "It assumes a known active source, uniform rock, a horizontal cylindrical void and favourable "
                "readout. None of those hold at Giza."],
        nums=[("high SNR", f"{hi[1]['detected_within_20m_rate']:.0%}"), ("low SNR", f"{lo[1]['detected_within_20m_rate']:.0%}")],
        figs=["results/validation/noise_transition.png", "results/validation/depth_recovery.png",
              "results/validation/coupled_measurement.png", "results/validation/source_ambiguity.png",
              "results/validation/assumption_stress.png", "results/validation/frequency_depth.png",
              "results/validation/independent_controls.png", "results/validation/motion_validation.png"],
        code="experiments/validation.py"))
    f = json.loads((RESULTS / "feasibility" / "metrics.json").read_text())
    out.append(dict(
        key="feasibility", tag="mechanism", eyebrow="Steelman · measurement chain",
        title="The phase channel and its units work in a controlled acquisition",
        lead="A known displacement passes through focusing and is recovered from the complex image.",
        points=["Positive controls confirm phase, sign and units survive an idealised focusing operation.",
                "Artificial detection datasets isolate the confounding a learned detector would face.",
                "These are favourable resolved targets and perfect platform compensation, not distributed "
                "pyramid scattering."],
        nums=[], figs=["results/feasibility/phase_readout.png", "results/feasibility/motion_detection.png",
                       "results/feasibility/cavity_wave_response.png", "results/feasibility/ai_controls.png"],
        code="experiments/feasibility.py"))
    r = json.loads((RESULTS / "robustness" / "metrics.json").read_text())
    o = r["ordering_control"]
    out.append(dict(
        key="robustness", tag="null", eyebrow="Controls · stress tests",
        title="The simulated examples survive repetition and reseeding",
        lead=f"Across {r['runs']} two-scatterer runs the median peak sits {r['median_peak_error_m']:.2f} m from the "
             f"spacing prediction, at most {r['max_peak_error_m']:.2f} m.",
        points=[f"The permutation test flags {o['ordered_rejection_fraction']:.0%} of planted ordered trajectories "
                f"and {o['noise_rejection_fraction']:.1%} of pure noise, at a nominal 5 per cent level.",
                "End-to-end reruns at three declared wavelengths reproduce the predicted rescaling."],
        nums=[("median error", f"{r['median_peak_error_m']:.2f} m"), ("worst", f"{r['max_peak_error_m']:.2f} m")],
        figs=["results/robustness/robustness.png"], code="experiments/robustness.py"))
    html_out = []
    (DOCS / "data").mkdir(exist_ok=True)
    for c in out:
        # publish the source metrics so every card on the page has its data behind it
        shutil.copyfile(RESULTS / c["key"] / "metrics.json", DOCS / "data" / f"{c['key']}.json")
        src = copy_figure(c["key"], c["figs"][0])
        extra = "".join(f'<a href="{copy_figure(c["key"], p)}"><img src="{copy_figure(c["key"], p)}" '
                        f'alt="" loading="lazy"></a>' for p in c["figs"][1:])
        pts = "".join(f"<li>{esc(p)}</li>" for p in c["points"])
        nums = "".join(f"<div>{esc(k)}<b>{esc(val)}</b></div>" for k, val in c["nums"])
        html_out.append(f"""<article class="db-card" id="{esc(c['key'])}">
<figure><a href="{src}"><img src="{src}" alt="" loading="lazy"></a></figure><div class="db-card-body">
<span class="db-tag {esc(c['tag'])}">{esc(c['eyebrow'])}</span><h3>{esc(c['title'])}</h3>
<p class="db-find">{esc(c['lead'])}</p><ul>{pts}</ul>
{f'<div class="db-num">{nums}</div>' if nums else ''}
<details><summary>Figures &amp; code</summary>
<div class="db-gallery">{extra}</div>
<p><a href="{REPO}/blob/main/{c['code']}">Code ↗</a> · <a href="data/{esc(c['key'])}.json" download>Results JSON ↓</a></p></details></div></article>""")
    return "".join(html_out)


VIEWER = """<div class="db-viewer">
  <div class="db-viewer-stage">
    <canvas id="db-canvas" aria-label="Three-dimensional depth volume. Drag to rotate, scroll to zoom."></canvas>
    <div class="db-hint" id="db-hint"></div>
    <div class="db-axes">surface at the top · depth downward<br>drag to rotate · scroll to zoom</div>
  </div>
  <div class="db-viewer-side">
    <h4>Patch</h4><div class="db-chips" id="db-chips"></div>
    <h4>Rendering</h4>
    <div class="db-slider"><label for="db-threshold"><span>Threshold</span><span id="db-threshold-v">0.42</span></label>
      <input type="range" id="db-threshold" min="0" max="0.95" step="0.01" value="0.42"></div>
    <div class="db-slider"><label for="db-density"><span>Opacity</span><span id="db-density-v">0.30</span></label>
      <input type="range" id="db-density" min="0.02" max="0.6" step="0.01" value="0.30"></div>
    <div class="db-slider"><label for="db-clip"><span>Depth shown</span><span id="db-clip-v">100%</span></label>
      <input type="range" id="db-clip" min="0.05" max="1" step="0.01" value="1"></div>
    <div class="db-slider"><label for="db-surface"><span>Radar surface</span><span id="db-surface-v">0.85</span></label>
      <input type="range" id="db-surface" min="0" max="1" step="0.05" value="0.85"></div>
    <label class="db-check"><input type="checkbox" id="db-repeat" checked> Mark where the volume repeats</label>
    <label class="db-check"><input type="checkbox" id="db-voids" checked> Mark surveyed void depths</label>
    <label class="db-check"><input type="checkbox" id="db-survey" checked> Show surveyed chambers</label>
    <p class="db-survey-key"><i style="background:#4ce09e"></i>chamber <i style="background:#6bb8fa"></i>passage <i style="background:#ffb340"></i>muon-detected void</p>
    <div class="db-legend"><span>low</span><i></i><span>high</span></div>
    <h4>This volume</h4><div class="db-readout" id="db-readout"></div>
  </div>
</div>
<p class="db-viewer-note"><b>What you are looking at.</b> The grey plane on top is the real radar image of that
patch, at 27 cm sampling. Everything below it is what the published method outputs as depth, rendered as a volume.
The red plane marks where that output starts repeating exactly, because the steering basis has come back to itself.
Switch between a monument and open ground, or between Giza and Sacsayhuamán: if the method were responding
to buried structure, these would not look alike. On Khufu and Khafre the wireframes are the real chambers, passages
and muon-detected voids, placed from excavation surveys and published muon imaging rather than from the radar; the
volume simply repeats past them. On the cemetery patches a green band marks the depth range of hundreds of excavated
burial shafts.
<a href="underworld.html">Open the full viewer ↗</a></p>"""


def acquisitions_band():
    rep = report(ACQ)
    if not rep or not rep.get("figures"):
        return ""
    item = rep["figures"][0]
    src = copy_figure(ACQ, item["file"])
    facts = rep["metrics"]["acquisitions"]
    cells = "".join(
        f'<div><span>{esc(a["label"].split(" · ")[0])}</span>'
        f'{esc(a["satellite"])} · {esc(a["collected"][:10])} · '
        f'{a["shape"][0]:,} × {a["shape"][1]:,} samples · {a["collection_duration_s"]:.1f} s · '
        f'{a["incidence_deg"]:.1f}° incidence · {a["patches"]} patches</div>'
        for a in facts)
    (DOCS / "data").mkdir(exist_ok=True)
    (DOCS / "data" / f"{ACQ}.json").write_text(json.dumps(rep, indent=2, allow_nan=False))
    return (f'<figure class="db-band" id="{esc(ACQ)}"><a href="{src}"><img src="{src}" alt="Ten patches of two '
            f'ICEYE dwell acquisitions" loading="lazy"></a>'
            f'<figcaption>{item.get("caption", "")}</figcaption>'
            f'<div class="db-band-facts">{cells}</div></figure>')


ARGUMENT = [
    ("Radar does not go through rock",
     "At this wavelength the signal dies about half a metre into dry limestone. Nothing underground is lit up, so "
     "any claim about what is down there has to be read off how the surface moves.",
     "r14_budget"),
    ("The surface movement in question is very slow",
     "Ground shaken by distant oceans and cities moves with periods of several seconds. To measure that you have "
     "to watch the same ground for minutes. One radar pass gives under two seconds before the ground stops "
     "looking like itself from a changed angle.",
     "r09_coherence"),
    ("So the movement is too small and too brief to see",
     "We measured what one pass can actually detect on the brightest real targets. It is about eight times larger "
     "than the strongest ground motion it would need to sense, and the recording is eighty times too short.",
     "r04_velocity_floor"),
    ("Where voids are certainly known, nothing appears",
     "Giza's cemeteries hold hundreds of excavated burial shafts five to thirty metres deep, beside bare plateau "
     "with none. At the shafts' own depths, the two are the same.",
     "r13_known_voids"),
    ("And the depth numbers are manufactured",
     "The depth scale is set by a constant nobody measures, and the picture repeats exactly every few metres. Two "
     "halves of the same pass put the same ground at different depths.",
     "r02_depth_axis"),
]


def argument_block():
    items = "".join(
        f'<a class="db-step" href="#{esc(k)}"><span>{i + 1:02d}</span><b>{esc(t)}</b><p>{esc(b)}</p></a>'
        for i, (t, b, k) in enumerate(ARGUMENT))
    return f'<div class="db-argument">{items}</div>'


FALSIFIERS = [
    ("A positive control succeeding",
     "Locate a surveyed void at its known position with the processing settings frozen beforehand, and score the "
     "result against chance. R13 does this at cemetery scale and finds nothing; a single named chamber needs "
     "metre-level registration, which is the next task on the roadmap."),
    ("Depths that reproduce across look angle",
     "Split one acquisition into two disjoint halves, give each an identical sweep, and show that the same ground "
     "is placed at the same depth by both. R3 measures a correlation of 0.067 between halves against 0.074 for "
     "unrelated pixels."),
    ("Depths that survive their own free parameters",
     "Show a feature whose depth does not move when the declared sound wavelength or the sub-aperture width "
     "changes. At present the depth axis scales linearly with a constant no measurement constrains."),
    ("Coherence on a target class we have not measured",
     "A deployed corner reflector, or any target holding complex coherence beyond a few seconds of look "
     "separation, would lengthen the usable record and weaken R9 and R11. Nothing in either scene does."),
    ("Agreement between two acquisitions of the same ground",
     "Two dwells of the same plateau on different dates should place the same structure at the same depth. We hold "
     "one acquisition per site and cannot run this; it is the cleanest test we are missing."),
]

NOTES = [
    ("Measuring vibration from radar is real",
     "Estimating micro-motion from SAR sub-apertures is established and works: on ships, bridges and other "
     "strong, stable structures vibrating fast enough that very short lags suffice. Independent work in 2026 "
     "recovered 1 to 4 Hz motion on controlled reflectors. Nothing here questions that. The subsurface claim is a "
     "different regime in every variable that matters: a distributed natural target, a frequency band a hundred "
     "times lower, and an inference about depth rather than about motion."),
    ("No published figure is reproduced here",
     "The sub-aperture parameters, patch positions, nuisance geometry and wavelength model behind the Khafre "
     "images were never released, so no result on this page reproduces a specific published image. Everything "
     "reproduces the method as described, with settings chosen and published here. If the disagreement lies in "
     "those settings, releasing the originals would resolve it quickly."),
    ("The public derivative protocol, run unmodified",
     "Version 1.7 of the openly published derivative protocol was installed without changes, its own tests passed, "
     "and it was given a real Khafre crop. Its registration returns only the integer peak of the cross-correlation "
     "and produced one distinct displacement value across an entire target track, where real inter-look "
     "displacements are thousandths of a pixel. Its example geometry also assumes a one-kilometre aperture span "
     "where real state vectors give about ten, moving its ambiguity limit from 494 m to 50 m, below the 100 m grid "
     "it searches. Both are fixable."),
    ("Everything is open and rebuilds from source",
     "The two radar products are commercial acquisitions and are not redistributed, but every script, "
     "configuration, intermediate result and figure is in the repository, and this page is generated from the "
     "saved reports. Corrections to earlier versions of this work, including two of our own mistakes, are kept in "
     "the method audit. Disagreement is best expressed as a rerun."),
]


def falsifier_block():
    items = "".join(
        f"<tr><td><b>{esc(t)}</b></td><td>{esc(b)}</td></tr>" for t, b in FALSIFIERS)
    return f'<table class="db-table db-falsify"><tbody>{items}</tbody></table>'


def author_block():
    items = "".join(
        f'<article class="db-note"><h3>{esc(t)}</h3><p>{esc(b)}</p></article>' for t, b in NOTES)
    return f'<div class="db-grid two">{items}</div>'


def method_figure():
    rep = report("method")
    if not rep or not rep.get("figures"):
        return ""
    item = rep["figures"][0]
    src = copy_figure("method", item["file"])
    return (f'<figure style="margin:0 0 18px"><a href="{src}"><img src="{src}" alt="" loading="lazy" '
            f'style="width:100%;border:1px solid var(--line);border-radius:3px"></a>'
            f'<figcaption style="font-size:11.5px;line-height:1.6;color:var(--dim);margin-top:7px;max-width:1000px">'
            f'{item.get("caption", "")}</figcaption></figure>')


def method_block(reps):
    r01 = reps.get("r01_giza_controls")
    acq = r01["metrics"]["acquisition"] if r01 else {}
    rows = [
        ("Acquisitions", "ICEYE X33 Giza 2025-08-27 · ICEYE X35 Sacsayhuamán 2025-08-22, Spotlight Dwell Fine, VV"),
        ("Image", f"{acq.get('shape',[0,0])[0]:,} × {acq.get('shape',[0,0])[1]:,} complex samples, "
                  f"{acq.get('azimuth_spacing_m',0):.3f} m azimuth, {acq.get('range_spacing_m',0):.3f} m slant range"),
        ("Aperture", "24.5 s of slow time carried on the Doppler axis; the image row span is 0.713 s and is not the aperture"),
        ("Sub-apertures", "Reference and offset passbands swept across the processed band; the published designs overlap 96–99% at lags of 14–72 ms"),
        ("Registration", "32 × 32 complex cross-correlation, upsampled DFT to 1/1000 px, validated on this radar texture with planted shifts down to 0.005 px"),
        ("Depth geometry", "Kz = 4π B⊥ / (λs r sin θ) with B⊥ the along-track platform offset from the product's own 101 state vectors, not a scalar approximation"),
        ("Depth convention", "π/ΔKz is the range a real two-component fit can use; the complex output repeats at twice that"),
    ]
    body = "".join(f"<tr><td>{esc(k)}</td><td>{esc(v)}</td></tr>" for k, v in rows)
    return f'<table class="db-table"><tbody>{body}</tbody></table>'


SOURCES = [
    ("Original method", "Biondi &amp; Malanga (2022), SAR Doppler tomography of the Great Pyramid",
     "https://arxiv.org/abs/2208.00811"),
    ("Retraction", "Remote Sensing, 10 August 2026", "https://www.mdpi.com/2072-4292/18/16/2679"),
    ("Derivative protocol", "Replication and Verification, v1.5 and v1.7",
     "https://github.com/BiondiProtocol/Replication-and-Verification-Biondi-Protocol"),
    ("Measured SAR vibration", "Vattulainen et al. (2026), spaceborne SAR micro-motion",
     "https://ieeexplore.ieee.org/Xplore/home.jsp"),
    ("Registration", "Guizar-Sicairos et al. (2008), efficient subpixel image registration",
     "https://doi.org/10.1364/OL.33.000156"),
    ("Ambient ground motion", "Peterson (1993), seismic background noise",
     "https://pubs.usgs.gov/publication/ofr93322"),
    ("Product metadata", "ICEYE HDF5 product format", "https://sar.iceye.com/6.0.0/productFormats/metadata/"),
    ("Terrain", "Copernicus DEM GLO-30", "https://registry.opendata.aws/copernicus-dem/"),
    ("Our audit", "What is and is not established here", f"{REPO}/blob/main/METHOD_AUDIT.md"),
    ("Roadmap", "Project state and what comes next", f"{REPO}/blob/main/ROADMAP.md"),
]


def build():
    DOCS.mkdir(exist_ok=True)
    build_field_site(DOCS)          # the field atlas stays as a linked appendix
    reps = {name: r for name in REAL + SIM if (r := report(name))}
    real = [reps[n] for n in REAL if n in reps]
    sim = [reps[n] for n in SIM if n in reps]
    updated = date.fromisoformat(max(r.get("date", "2026-09-19") for r in reps.values())).strftime("%d %B %Y").lstrip("0")
    src = "".join(f'<li><span>{esc(a)}</span><a href="{esc(c)}">{b} ↗</a></li>' for a, b, c in SOURCES)
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SAR Depth, Tested</title>
<meta name="description" content="Single-image SAR depth tomography tested on two real ICEYE dwell acquisitions of Giza and Sacsayhuamán: fourteen measurements including a positive control against hundreds of surveyed burial shafts, with a three-dimensional viewer for the depth volumes.">
<link rel="stylesheet" href="style.css">
<link rel="stylesheet" href="dashboard.css">
<link rel="icon" href="favicon.svg">
<meta property="og:title" content="SAR Depth, Tested">
<meta property="og:description" content="Nine tests of single-image SAR depth tomography on two real ICEYE dwell acquisitions.">
<meta property="og:image" content="https://ejhong.github.io/sar/figs/r08_sections/r08_sections.jpg">
</head>
<body>
<header class="db-top"><div class="wrap">
  <a class="db-mark" href="./"><span aria-hidden="true">◒</span> SAR / DEPTH, TESTED</a>
  <nav class="db-nav">
    <a href="#argument">Argument</a><a href="#viewer">Viewer</a><a href="#scans">Scans</a><a href="#real">Real data</a><a href="#sim">Simulation</a><a href="#support">Steelman</a>
    <a href="#falsify">Falsification</a><a href="#scope">Scope</a><a href="#method">Method</a><a href="field.html">Field atlas</a><a href="{REPO}">Code ↗</a>
  </nav>
</div></header>

<section class="db-title"><div class="wrap">
  <h1>One radar image cannot map what is under the ground</h1>
  <p>Single-image SAR Doppler tomography, as published for the Great Pyramid in 2022 and announced for Khafre in
  2025, reimplemented and run on two real ICEYE Spotlight Dwell Fine acquisitions of Giza and Sacsayhuamán.
  Fourteen measurements, eight simulations, three supporting benchmarks. Every figure is generated by code in the
  repository.</p>
  <div class="db-verdict"><b>The claim needs four things at once, and each fails on its own.</b> X-band reaches
  about half a metre into dry limestone, so nothing below the surface is illuminated. Resolving the ground motion
  the method invokes needs a record of order 150 seconds; one acquisition gives under two before the looks stop
  resembling each other. That motion is 0.1 to 10 µm/s and the measured floor is 77. And surface waves long enough
  to sample 648 m cannot resolve better than about a kilometre laterally, against wells said to be ten metres
  across. Independently: the method finds no signal at the known depths of hundreds of surveyed burial shafts, it
  does not separate a pyramid from open desert in the same acquisition, and its depth axis is periodic by
  construction with a freely chosen scale.</div>
</div></section>

<section class="db-stats"><div class="wrap" style="display:contents">{stats(reps)}</div></section>

<section class="db-sec" id="argument"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">The argument</span><h2>Five steps, each one measured</h2></div>
  <p>Every step links to the experiment that establishes it. Each is enough on its own.</p></div>
  {argument_block()}
</div></section>

<section class="db-sec" id="viewer"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Interactive · WebGL</span><h2>The depth volumes, in three dimensions</h2></div>
  <p>Real output from the Giza acquisition. Compare a pyramid with empty plateau directly.</p></div>
  {VIEWER}
</div></section>

<section class="db-sec" id="scans"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">The data</span><h2>What the radar recorded</h2></div>
  <p>Ten patches, two acquisitions, no processing. Every later result is a statement about these pixels.</p></div>
  {acquisitions_band()}
</div></section>

<section class="db-sec" id="real"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Measurement · {len(real)} tests</span><h2>What the real acquisitions show</h2></div>
  <p>The budget first, then the positive control, then what the method distinguishes, then why, then whether our
  own code is to blame. Figures enlarge; method and limits are collapsed.</p></div>
  <div class="db-grid">{''.join(card(r) for r in real)}</div>
</div></section>

<section class="db-sec" id="sim"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Simulation · {len(sim)} tests</span><h2>What the method computes, in scenes whose contents are known</h2></div>
  <p>Synthetic scenes establish the mechanism: a motionless pyramid, empty desert, two rocks, ground that really moves.</p></div>
  <div class="db-grid">{''.join(card(r, thumb=False) for r in sim)}</div>
</div></section>

<section class="db-sec" id="support"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Supporting · 3 results</span><h2>Where the method is given every advantage</h2></div>
  <p>A coupled physical benchmark, the measurement chain it rests on, and the stress tests behind the simulated examples.</p></div>
  <div class="db-grid">{supporting_cards()}</div>
</div></section>

<section class="db-sec" id="falsify"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Falsification</span><h2>What would change this conclusion</h2></div>
  <p>Stated in advance and in testable form. Any one of these, demonstrated, would reopen the question.</p></div>
  {falsifier_block()}
</div></section>

<section class="db-sec" id="scope"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Scope</span><h2>What this does and does not claim</h2></div>
  <p>Worth reading before drawing a conclusion from anything above.</p></div>
  {author_block()}
  <div class="db-repro">
    <h3>Reproduce any number on this page</h3>
    <pre><code>git clone https://github.com/ejhong/sar &amp;&amp; cd sar
python3.11 -m venv .venv &amp;&amp; .venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest tests                      # 80 tests

cd experiments &amp;&amp; ../.venv/bin/python run_all.py      # simulation chapter, self-contained
../.venv/bin/python r13_known_voids.py                # needs the ICEYE products
cd .. &amp;&amp; .venv/bin/python build_dashboard.py         # rebuilds this page</code></pre>
    <p>Product paths are in <code>experiments/real_common.py</code>. Nothing reads a full raster; every step works
    on a bounded, geolocated crop and records its origin. What is and is not established, including corrections to
    our own earlier results, is kept in <a href="{REPO}/blob/main/METHOD_AUDIT.md">METHOD_AUDIT.md</a>.</p>
  </div>
</div></section>

<section class="db-sec" id="method"><div class="wrap">
  <div class="db-head"><div><span class="db-kicker">Method</span><h2>What was run</h2></div></div>
  {method_figure()}
  <div class="db-cols">
    <div>{method_block(reps)}</div>
    <div><h4 style="font:9px var(--mono);letter-spacing:.12em;text-transform:uppercase;color:var(--dim);margin:0 0 9px">Sources</h4>
    <ul class="db-sources">{src}</ul></div>
  </div>
</div></section>

<footer class="db-foot"><div class="wrap">
  Updated {esc(updated)} · every figure and number is produced by code in the repository and rebuilt from
  <code>results/*/summary.json</code> · <a href="{REPO}">source</a> ·
  <a href="{REPO}/blob/main/METHOD_AUDIT.md">audit</a> · <a href="{REPO}/blob/main/ROADMAP.md">roadmap</a>
</div></footer>

<script src="underworld.js"></script>
<script>
(function () {{
  const canvas = document.getElementById('db-canvas');
  const hint = document.getElementById('db-hint');
  let view;
  try {{ view = new Underworld(canvas, {{ set textContent(v) {{ hint.textContent = v; }} }}); }}
  catch (e) {{ hint.innerHTML = 'WebGL2 is needed for the volume view. The same volumes appear as sections in <a href="#r08_sections" style="color:#9ad">the section test</a>.'; return; }}
  const bind = (id, key, render) => {{
    const el = document.getElementById(id), out = document.getElementById(id + '-v');
    el.addEventListener('input', () => {{
      const v = parseFloat(el.value); view.settings[key] = v;
      if (out) out.textContent = render ? render(v) : v.toFixed(2);
      view.draw();
    }});
  }};
  bind('db-threshold', 'threshold'); bind('db-density', 'density');
  bind('db-clip', 'depthClip', v => Math.round(v * 100) + '%'); bind('db-surface', 'surfaceMix');
  const keyOf = {{'db-repeat': 'showRepeat', 'db-voids': 'showVoids', 'db-survey': 'showSurvey'}};
  ['db-repeat', 'db-voids', 'db-survey'].forEach(id => {{
    const el = document.getElementById(id), key = keyOf[id];
    el.addEventListener('change', () => {{ view.settings[key] = el.checked; view.draw(); }});
  }});
  const fmt = (v, n) => (v === null || v === undefined ? '—' : Number(v).toFixed(n));
  fetch('data/voxels/index.json').then(r => r.json()).then(index => {{
    const holder = document.getElementById('db-chips'), buttons = [];
    const show = async (bundle, button) => {{
      buttons.forEach(b => b.setAttribute('aria-pressed', String(b === button)));
      let h; try {{ h = await view.load(bundle); }} catch (e) {{ hint.textContent = 'Could not load ' + bundle.patch; return; }}
      const ex = h.extent_m, dz = h.depth_axis;
      hint.innerHTML = '<b>' + h.label + '</b><br>' + Math.round(ex.along_azimuth) + ' × ' +
        Math.round(ex.across_slant_range) + ' m of ground<br>' + Math.round(ex.depth) + ' m of model depth';
      document.getElementById('db-readout').innerHTML = [
        ['Satellite', h.acquisition.satellite], ['Collected', h.acquisition.collection_start.slice(0, 10)],
        ['Voxels', h.shape.join(' × ')],
        ['Voxel size', fmt(h.spacing_m.along_azimuth, 2) + ' × ' + fmt(h.spacing_m.across_slant_range, 2) + ' × ' + fmt(h.spacing_m.depth, 2) + ' m'],
        ['Declared λs', fmt(dz.lambda_s_m, 2) + ' m'], ['Repeats every', fmt(dz.repeat_period_m, 1) + ' m'],
        ['Surveyed voids', h.known_voids ? fmt(h.known_voids.depth_min_m, 0) + '–' + fmt(h.known_voids.depth_max_m, 0) + ' m' : 'none here'],
        ['Surveyed chambers', h.survey && h.survey.features ? h.survey.features.length + ' drawn' : 'none here'],
        ['Terrain', h.terrain_height_m === null ? '—' : fmt(h.terrain_height_m, 0) + ' m'],
      ].map(([k, v]) => '<div><span>' + k + '</span><span>' + v + '</span></div>').join('');
    }};
    const sites = [];
    index.bundles.forEach(bundle => {{
      let g = sites.find(x => x.site === bundle.site);
      if (!g) {{ g = {{ site: bundle.site, items: [] }}; sites.push(g); }}
      g.items.push(bundle);
    }});
    sites.forEach((group, gi) => {{
      if (sites.length > 1) {{
        const h = document.createElement('p'); h.className = 'db-site';
        h.textContent = (index.sites && index.sites[group.site]) || group.site;
        holder.appendChild(h);
      }}
      group.items.forEach((bundle, i) => {{
        const b = document.createElement('button'); b.type = 'button'; b.setAttribute('aria-pressed', 'false');
        b.innerHTML = (bundle.short_label || bundle.label) + '<small>' + bundle.kind + '</small>';
        b.addEventListener('click', () => show(bundle, b));
        holder.appendChild(b); buttons.push(b);
        if (gi === 0 && i === 0) show(bundle, b);
      }});
    }});
  }}).catch(() => {{ hint.textContent = 'No volume bundles found.'; }});
}})();
</script>
</body>
</html>
"""
    (DOCS / "index.html").write_text(page)
    print(f"Built docs/index.html as a dashboard: {len(real)} real tests, {len(sim)} simulation tests")


if __name__ == "__main__":
    build()
