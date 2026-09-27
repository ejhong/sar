from .sites import export_sites

if __name__ == '__main__':
    for s in export_sites():
        print(f"{s['id']:22s} {s['kind']:5s} features {s['features']} volumes {s['volumes']}")
