"""Checks for publishing a portable report with intact evidence and downloads."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json
import csv
import pytest

import build_site


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids, self.references, self.images = [], [], []
        self.experiments, self.details = [], {}
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.append(attrs["id"])
            if "experiment" in attrs.get("class", "").split():
                self.experiments.append(attrs["id"])
            if tag == "details":
                self.details[attrs["id"]] = attrs
        for name in ("href", "src"):
            if attrs.get(name):
                self.references.append(attrs[name])
        if tag == "img" and attrs.get("src"):
            self.images.append(attrs)


def test_published_page_has_valid_local_links_and_accessible_figures():
    page = Page((build_site.DOCS / "index.html").read_text())
    assert len(page.ids) == len(set(page.ids))
    for reference in page.references:
        parsed = urlsplit(reference)
        if parsed.scheme or parsed.netloc:
            continue
        if parsed.path:
            path = build_site.DOCS / unquote(parsed.path)
            assert path.is_file(), reference
        elif parsed.fragment:
            assert parsed.fragment in page.ids, reference
    # Check the important evidence directly rather than requiring more pictures.
    assert {"signal", "depth", "ambiguity", "field", "scenario-lab"} <= set(page.ids)
    assert any(image['src'].endswith('noise_transition.png') for image in page.images)
    assert not any('/t07_wells/' in image['src'] for image in page.images)
    for image in page.images:
        assert image.get("alt"), image
        assert int(image["width"]) > 0 and int(image["height"]) > 0


def test_old_absolute_figure_paths_resolve_in_a_new_checkout():
    file = next((build_site.RESULTS / "t03_mechanism" / "figs").glob("*.png"))
    relative = file.relative_to(build_site.ROOT)
    assert build_site.resolve_figure(f"/different/machine/checkout/{relative}") == file


def test_downloads_are_portable_and_match_the_reports():
    for report in build_site.RESULTS.glob("t*/summary.json"):
        source = build_site.refine_summary(json.loads(report.read_text()))
        published = json.loads((build_site.DOCS / "data" / f"{source['id']}.json").read_text())
        assert source["metrics"] == published["metrics"]
        assert source["finding"] == published["finding"]
        assert all(not Path(f["file"]).is_absolute() for f in published["figures"])


def test_site_fonts_are_packaged_with_their_licenses():
    for font, license in (("newsreader.ttf", "Newsreader-LICENSE.txt"), ("ibm-plex-mono.ttf", "IBMPlexMono-LICENSE.txt")):
        assert (build_site.DOCS / "fonts" / font).stat().st_size > 10_000
        assert "SIL OPEN FONT LICENSE" in (build_site.DOCS / "fonts" / license).read_text()


def test_field_atlas_has_intact_links_and_portable_source_assets():
    page=Page((build_site.DOCS/'field.html').read_text())
    assert len(page.ids)==len(set(page.ids))
    assert {'scene','survey','checks','registration-check','aperture-check','translation-check'} <= set(page.ids)
    for reference in page.references:
        parsed=urlsplit(reference)
        if parsed.scheme or parsed.netloc:
            continue
        if parsed.path:
            assert (build_site.DOCS/unquote(parsed.path)).is_file(),reference
        elif parsed.fragment:
            assert parsed.fragment in page.ids,reference
    for source,destination in [('main.css','style.css'),('main.js','site.js'),('field.css','field.css'),('model-viewer.js','model-viewer.js')]:
        assert (build_site.ROOT/'site/assets'/source).read_bytes()==(build_site.DOCS/destination).read_bytes()
    for path in (build_site.DOCS/'data/field').rglob('*.json'):
        assert '/Users/' not in path.read_text(),path


def test_constructed_illustration_is_separate_from_numbered_experiments():
    page = Page((build_site.DOCS / "index.html").read_text())
    assert set(page.experiments) == set(build_site.EXPERIMENTS)
    assert len(page.experiments) == 7
    assert "t07_surface_controls" in page.experiments
    assert "t07_wells" in page.details
    assert "t07_wells" not in page.experiments
    archived = json.loads((build_site.DOCS / "data" / "t07_wells.json").read_text())
    assert archived["role"] == "illustration"


def test_interactive_results_are_exactly_the_committed_trials():
    import re
    text = (build_site.DOCS / "index.html").read_text()
    embedded = re.search(r'<script id="scenario-data" type="application/json">(.*?)</script>', text, re.S)
    metrics = json.loads((build_site.RESULTS / "validation" / "metrics.json").read_text())
    assert json.loads(embedded.group(1)) == metrics['noise_sweep']['rows']
    assert json.loads((build_site.DOCS / "data" / "validation_metrics.json").read_text()) == metrics
    assert metrics['frequency_depth']['accepted_reflectors'] == sum(r['accepted'] for r in metrics['frequency_depth']['rows'])
    assert set(metrics['design']['evaluation_depths_m']).isdisjoint(metrics['design']['dictionary_depths_m'])
    assert metrics['depth']['oracle_and_roundtrip']['focused_slc_roundtrip_error_m'] < 1e-12


@pytest.mark.parametrize('filename,section,number_key', [
    ('snr_trials.csv', 'noise_sweep', 'snr_db'),
    ('depth_trials.csv', 'depth', 'amplitude_um'),
])
def test_trial_records_independently_reproduce_displayed_rates(filename, section, number_key):
    root = build_site.RESULTS / 'validation'
    metrics = json.loads((root / 'metrics.json').read_text())
    with (root / filename).open() as stream:
        rows = list(csv.DictReader(stream))
    groups = metrics[section]['rows' if section == 'noise_sweep' else 'groups']
    for group in groups:
        selected = [r for r in rows if r['family'] == group['family'] and r['estimator'] == group['estimator']
                    and float(r[number_key]) == group[number_key]]
        positive = [r for r in selected if float(r['true_depth_m']) > 0]
        negative = [r for r in selected if float(r['true_depth_m']) == 0]
        assert len(positive) == group['positive_trials']
        assert len(negative) == group['null_trials']
        assert sum(r['accepted'] == 'True' for r in positive)/len(positive) == group['detection_rate']
        assert sum(r['accepted'] == 'True' for r in negative)/len(negative) == group['false_alarm_rate']
        assert sum(r['within_20m'] == 'True' for r in positive)/len(positive) == group['detected_within_20m_rate']
