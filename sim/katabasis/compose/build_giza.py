"""Regenerate sites/giza/site.json and terrain.json from their sources.

    uv run python -m katabasis.compose.build_giza

Monuments: Khufu and Khafre from the first investigation's catalogue
(sim/legacy/catalog/giza/underworld.json, which cites Petrie 1883 and
Maragioglio & Rinaldi); Menkaure and the Sphinx from their Wikipedia
infoboxes (coordinates and dimensions quoted in the site file).
Terrain: Copernicus GLO-30 DSM with the pyramid footprints inpainted.
"""
from __future__ import annotations

import json
from pathlib import Path

from .dem import local_crop, write_crop
from .geo import to_local

ROOT = Path(__file__).resolve().parents[3]
LEGACY = ROOT / 'sim' / 'legacy' / 'catalog' / 'giza' / 'underworld.json'
OUT = ROOT / 'sites' / 'giza'
ORIGIN = (29.976025, 31.130794)          # Khafre base centre
EXTENT = {'x': [-640.0, 960.0], 'y': [-760.0, 640.0], 'z': [-120.0, 220.0]}
KIND = {'chamber': 'chamber', 'passage': 'passage', 'void': 'void'}


def main() -> None:
    legacy = json.loads(LEGACY.read_text())
    structures, features = [], []
    footprints = []
    for mon in legacy['monuments']:
        e, n, _ = to_local(mon['latitude'], mon['longitude'], *ORIGIN)
        base = mon['base_elevation_m']
        structures.append({
            'id': mon['id'], 'name': mon['name'], 'material': 'masonry',
            'shape': {'type': 'pyramid', 'centre': [round(float(e), 2), round(float(n), 2), base],
                      'base': mon['base_side_m'], 'height': mon['height_m'], 'yaw_deg': 0.0},
            'source': 'sim/legacy/catalog/giza/underworld.json (Petrie 1883; Maragioglio & Rinaldi)'})
        footprints.append({'centre': [round(float(e), 1), round(float(n), 1)], 'half': mon['base_side_m'] / 2 + 30})
        # The legacy Khufu table puts the North Face Corridor at y = -110, i.e. its y runs south;
        # its Khafre table runs north. Both are brought to y = north here.
        flip = -1.0 if mon['id'] == 'khufu' else 1.0
        for f in mon['features']:
            x, y, z = f['centre']
            shape = {'type': 'box', 'centre': [round(float(e) + x, 2), round(float(n) + flip * y, 2), round(base + z, 2)],
                     'size': f['size']}
            if 'incline_deg' in f:
                shape['pitch_deg'] = -f['incline_deg']     # legacy sign is the sense of travel
            features.append({
                'id': f"{mon['id']}-{f['name'].lower().replace(' ', '-').replace(chr(39), '')}",
                'name': f"{f['name']} ({mon['name']})", 'kind': KIND.get(f['kind'], 'void'),
                'status': 'surveyed', 'fill': 'air', 'shape': shape, 'source': f['source'],
                'note': (f.get('note', '') + ' Simplified box; placement approximate until re-derived from Petrie (1883).').strip(),
                'placement': 'approximate'})

    menk = to_local(29.97250, 31.12833, *ORIGIN)
    sphinx = to_local(29.97526, 31.13758, *ORIGIN)
    footprints.append({'centre': [round(float(menk[0]), 1), round(float(menk[1]), 1)], 'half': 104.6 / 2 + 25})
    structures.append({
        'id': 'menkaure', 'name': 'Menkaure', 'material': 'masonry',
        'shape': {'type': 'pyramid', 'centre': [round(float(menk[0]), 2), round(float(menk[1]), 2), None],
                  'base': 103.4, 'height': 61.0, 'yaw_deg': 0.0},
        'source': 'Wikipedia, Pyramid of Menkaure (infobox): "29°58′21″N 31°07′42″E", base "102.2 by 104.6 metres", current height "61 m". Base drawn square at the mean side; base elevation from the terrain.'})

    crop = local_crop(*ORIGIN, tuple(EXTENT['x']), tuple(EXTENT['y']), 20.0, inpaint=footprints)
    OUT.mkdir(parents=True, exist_ok=True)
    write_crop(OUT / 'terrain.json', crop)

    # Menkaure's base elevation: the inpainted ground at its centre
    import numpy as np
    zs = np.asarray(crop['z'])
    ix = int(round((menk[0] - crop['x0']) / crop['dx']))
    iy = int(round((menk[1] - crop['y0']) / crop['dx']))
    structures[-1]['shape']['centre'][2] = round(float(zs[ix, iy]), 2)
    sphinx_z = round(float(zs[int(round((sphinx[0] - crop['x0']) / crop['dx'])), int(round((sphinx[1] - crop['y0']) / crop['dx']))]), 1)

    site = {
        'id': 'giza', 'name': 'Giza plateau', 'kind': 'real',
        'summary': 'Khufu, Khafre and Menkaure on Mokattam limestone, with the chambers and passages already surveyed inside and beneath them.',
        'frame': {'axes': 'x east, y north, z up; metres', 'origin': {'latitude': ORIGIN[0], 'longitude': ORIGIN[1], 'label': 'Khafre base centre'},
                  'vertical': 'z = metres above EGM2008 (orthometric)'},
        'extent': EXTENT,
        'terrain': {'type': 'dem', 'file': 'terrain.json', 'datum_m': 0.0},
        'cover': [],
        'strata': [
            {'name': 'Mokattam Formation', 'material': 'limestone-mokattam', 'top': 'surface'},
            {'name': 'Older limestone (assumed)', 'material': 'limestone', 'top': -60.0},
        ],
        'water_table': {'z': 15.0},
        'structures': structures,
        'features': features,
        'landmarks': [{'id': 'sphinx', 'name': 'Great Sphinx', 'position': [round(float(sphinx[0]), 1), round(float(sphinx[1]), 1), sphinx_z],
                       'source': 'Wikipedia, Great Sphinx of Giza (infobox): "29.97526°N 31.13758°E"; length "73 metres", height "20 metres".'}],
        'notes': [
            'Bedrock: the plateau is Mokattam limestone at or near the surface; its base (-60 m here) and the rock beneath are assumed, not surveyed.',
            'Water table: "The average water table elevation is about +15 m" (Sharafeldin et al. 2019, abstract).',
            'Chambers and passages are simplified boxes from the first investigation; their placement is approximate (the legacy Khufu table ran y southward and is corrected here) and will be re-derived from Petrie (1883) before a Giza simulation relies on them.',
            'Terrain is a DSM: monument footprints are replaced by the surrounding ground and the pyramids are drawn from their survey dimensions.',
        ],
        'sources': {'sharafeldin2019': 'doi:10.5194/gi-8-29-2019', 'copernicus': 'Copernicus GLO-30 DEM, ESA/Airbus, via copernicus-dem-30m.s3.amazonaws.com'},
    }
    (OUT / 'site.json').write_text(json.dumps(site, indent=1))
    print(f'wrote {OUT/"site.json"} ({len(features)} features) and terrain {crop["nx"]}x{crop["ny"]}')


if __name__ == '__main__':
    main()
