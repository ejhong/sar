import json

from .radar import export_radar
from .runs import export_runs
from .sites import DATA, export_sites
from .volumes import export_p1_02, export_p1_03, export_p1_04, export_p1_05, export_p2_22

KEPT = ('volumes.json', 'wavefields.json')


def published():
    """Each site's published volume and wavefield lists, read before anything is rewritten."""
    out = {}
    for f in KEPT:
        for p in (DATA / 'sites').glob(f'*/{f}'):
            out[(p.parent.name, f)] = json.loads(p.read_text())
    return out


def keep_published(before):
    """Put back every published entry this run did not regenerate (its arrays are on another machine): an export run
    where some arrays are absent leaves those volumes as they are instead of dropping them."""
    for (site, f), old in before.items():
        p = DATA / 'sites' / site / f
        new = json.loads(p.read_text()) if p.exists() else []
        have = {e.get('id') for e in new}
        lost = [e for e in old if e.get('id') not in have]
        if lost:
            p.write_text(json.dumps(new + lost, indent=1))
            print(f"{site} {f}: kept {len(lost)} published entries whose arrays are not here")


if __name__ == '__main__':
    before = published()
    for fn in (export_p1_02, export_p1_03, export_p1_04, export_p1_05, export_p2_22):
        try:
            for v in fn():
                print(f"volume {v['id']}")
        except FileNotFoundError as e:         # a run's arrays are not committed; its exported volumes stay as they are
            print(f'{fn.__name__}: kept the published volumes ({e.filename} is not here)')
    for r in export_radar():                  # after the tomograms: export_p1_02 rewrites volumes.json
        print(f"radar {r['site']}: {', '.join(r['volumes'])}")
    keep_published(before)
    for s in export_sites():
        print(f"{s['id']:22s} {s['kind']:5s} features {s['features']} volumes {s['volumes']}")
    for r in export_runs():
        print(f'run {r}')
