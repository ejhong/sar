from .runs import export_runs
from .sites import export_sites
from .volumes import export_p1_02, export_p1_03

if __name__ == '__main__':
    for v in export_p1_02() + export_p1_03():
        print(f"volume {v['id']}")
    for s in export_sites():
        print(f"{s['id']:22s} {s['kind']:5s} features {s['features']} volumes {s['volumes']}")
    for r in export_runs():
        print(f'run {r}')
