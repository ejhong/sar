"""Checks for publishing a portable report with intact evidence and downloads."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json

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
    assert len(page.images) >= 30
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


def test_constructed_illustration_is_separate_from_numbered_experiments():
    page = Page((build_site.DOCS / "index.html").read_text())
    assert set(page.experiments) == set(build_site.EXPERIMENTS)
    assert len(page.experiments) == 7
    assert "t07_surface_controls" in page.experiments
    assert "t07_wells" in page.details
    assert "t07_wells" not in page.experiments
    archived = json.loads((build_site.DOCS / "data" / "t07_wells.json").read_text())
    assert archived["role"] == "illustration"
