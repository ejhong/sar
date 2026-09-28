from .radar import export_radar
from .runs import export_runs
from .sites import export_sites
from .volumes import export_p1_02, export_p1_03, export_p1_04, export_p1_05

if __name__ == '__main__':
    for fn in (export_p1_02, export_p1_03, export_p1_04, export_p1_05):
        try:
            for v in fn():
                print(f"volume {v['id']}")
        except FileNotFoundError as e:         # a run's arrays are not committed; its exported volumes stay as they are
            print(f'{fn.__name__}: kept the published volumes ({e.filename} is not here)')
    for r in export_radar():                  # after the tomograms: export_p1_02 rewrites volumes.json
        print(f"radar {r['site']}: {', '.join(r['volumes'])}")
    for s in export_sites():
        print(f"{s['id']:22s} {s['kind']:5s} features {s['features']} volumes {s['volumes']}")
    for r in export_runs():
        print(f'run {r}')
