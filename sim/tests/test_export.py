"""What the website shows must be what the code produces now."""
import json

from katabasis.export.sites import DATA, scene
from katabasis.compose import list_sites, load_site


def test_exported_scenes_are_current():
    for sid in list_sites():
        on_disk = json.loads((DATA / 'sites' / sid / 'scene.json').read_text())
        fresh = json.loads(json.dumps(scene(load_site(sid)), separators=(',', ':'), ensure_ascii=False))
        for k in ('volumes', 'surveys', 'wavefields', 'radar'):
            on_disk.pop(k, None)
            fresh.pop(k, None)
        assert on_disk == fresh, f'{sid}: run `uv run python -m katabasis.export`'
