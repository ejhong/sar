"""Checks that the published dashboard is portable and its evidence intact."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json

import pytest

import build_dashboard

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids, self.references, self.images, self.cards = [], [], [], []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.append(attrs["id"])
            if "db-card" in attrs.get("class", "").split():
                self.cards.append(attrs["id"])
        for key in ("href", "src"):
            if key in attrs:
                self.references.append(attrs[key])
        if tag == "img" and "src" in attrs:
            self.images.append(attrs["src"])


@pytest.fixture(scope="module")
def page():
    assert (DOCS / "index.html").is_file(), "run build_dashboard.py first"
    return Page((DOCS / "index.html").read_text())


def local(reference):
    parts = urlsplit(reference)
    return not parts.scheme and not parts.netloc and parts.path


def test_every_local_link_resolves(page):
    missing = []
    for reference in page.references:
        if not local(reference):
            continue
        target = unquote(urlsplit(reference).path)
        if not (DOCS / target).exists():
            missing.append(reference)
    assert not missing, f"broken local references: {missing}"


def test_every_anchor_has_a_target(page):
    anchors = {r[1:] for r in page.references if r.startswith("#") and len(r) > 1}
    assert anchors <= set(page.ids), f"dangling anchors: {sorted(anchors - set(page.ids))}"


def test_images_exist_and_are_not_empty(page):
    assert page.images, "the page shows no figures"
    for src in page.images:
        path = DOCS / unquote(urlsplit(src).path)
        assert path.is_file() and path.stat().st_size > 0, f"missing figure {src}"


SUPPORTING = {"validation", "feasibility", "robustness"}


def test_every_result_card_publishes_its_data(page):
    assert len(page.cards) >= 12, "expected the full set of result cards"
    for key in page.cards:
        data = DOCS / "data" / f"{key}.json"
        assert data.is_file(), f"no published data for {key}"
        report = json.loads(data.read_text())
        if key in SUPPORTING:
            assert report.get("date"), f"{key} metrics carry no date"
            continue
        for field in ("id", "title", "finding", "limitations", "method", "metrics"):
            assert report.get(field), f"{key} is missing {field}"
        assert report["id"] == key


def test_published_reports_match_the_results_tree(page):
    for key in page.cards:
        published = json.loads((DOCS / "data" / f"{key}.json").read_text())
        name = "metrics.json" if key in SUPPORTING else "summary.json"
        source = ROOT / "results" / key / name
        assert source.is_file(), f"{key} has no source {name}"
        assert published == json.loads(source.read_text()), f"{key} is stale; rebuild the site"


def test_real_results_lead_the_page(page):
    text = (DOCS / "index.html").read_text()
    assert text.index('id="viewer"') < text.index('id="real"') < text.index('id="sim"'), \
        "the viewer and real data must come before the simulation chapter"


def test_headings_are_statements_not_questions(page):
    """STYLE.md: headings state findings."""
    text = (DOCS / "index.html").read_text()
    import re
    headings = re.findall(r"<h[123][^>]*>(.*?)</h[123]>", text, re.S)
    questions = [h for h in headings if "?" in h]
    assert not questions, f"question headings: {questions}"


def test_volume_bundles_are_complete_and_bounded(page):
    index = DOCS / "data" / "voxels" / "index.json"
    assert index.is_file(), "the viewer has no bundle index"
    manifest = json.loads(index.read_text())
    assert manifest["bundles"], "no volume bundles published"
    for bundle in manifest["bundles"]:
        base = DOCS / "data" / "voxels" / bundle["site"] / bundle["patch"]
        header = json.loads(base.with_suffix(".json").read_text())
        raw = base.with_suffix(".bin")
        assert raw.is_file() and raw.stat().st_size == header["bytes"], f"{bundle['patch']} volume is wrong size"
        assert raw.stat().st_size <= 4_000_000, "bundle too large for a web page"
        assert base.with_suffix(".png").is_file(), f"{bundle['patch']} has no surface texture"
        nx, ny, nz = header["shape"]
        assert nx * ny * nz == header["bytes"]
        assert header["depth_axis"]["repeat_period_m"] > 0


def test_fonts_ship_with_their_licences():
    for font in (DOCS / "fonts").glob("*.ttf"):
        stem = font.stem.split("-")[0].replace("ibm", "IBMPlex").replace("_", "")
        licences = list((DOCS / "fonts").glob("*LICENSE*"))
        assert licences, "no font licences published"


def test_field_atlas_still_builds_and_links(page):
    assert (DOCS / "field.html").is_file()
    atlas = Page((DOCS / "field.html").read_text())
    missing = [r for r in atlas.references if local(r) and not (DOCS / unquote(urlsplit(r).path)).exists()]
    assert not missing, f"broken field atlas references: {missing}"
