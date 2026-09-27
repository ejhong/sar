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

import numpy as np

from .dem import local_crop, write_crop
from .geo import to_local

ROOT = Path(__file__).resolve().parents[3]
LEGACY = ROOT / 'sim' / 'legacy' / 'catalog' / 'giza' / 'underworld.json'
OUT = ROOT / 'sites' / 'giza'
ORIGIN = (29.976025, 31.130794)          # Khafre base centre
EXTENT = {'x': [-640.0, 960.0], 'y': [-760.0, 640.0], 'z': [-120.0, 220.0]}
KIND = {'chamber': 'chamber', 'passage': 'passage', 'void': 'void'}
INCH = 0.0254
PETRIE = 'Petrie (1883), The Pyramids and Temples of Gizeh, sec. 64, summary of interior positions (inches from the centre and above the pavement)'


def khufu_petrie() -> list[dict]:
    """Khufu's chambers and passages from Petrie's summary (x east, y north, z above the pavement; metres).

    Each passage runs between two tabulated points; chambers use tabulated walls, floors and roofs with the
    dimensions Petrie gives in the same chapter. Passage bores are 41.6 x 47.3 inches.
    """
    def pt(n_from_centre, e, z):
        return np.array([e * INCH, n_from_centre * INCH, z * INCH])

    def passage(a, b, bore=(41.6 * INCH, 47.3 * INCH)):
        mid = (a + b) / 2
        d = b - a
        horiz = float(np.hypot(d[0], d[1]))
        length = float(np.linalg.norm(d))
        pitch = float(np.degrees(np.arctan2(d[2], horiz)))
        # the box's local y runs from a to b; pitch lifts +y, so orient it southward-positive when b is south of a
        yaw = float(np.degrees(np.arctan2(-d[0], d[1])))
        return {'type': 'box', 'centre': mid.round(2).tolist(), 'size': [round(bore[0], 2), round(length, 2), round(bore[1], 2)],
                'yaw_deg': round(yaw, 2), 'pitch_deg': round(pitch, 2)}

    # entrance (descending) passage: doorway to its south end (sec. 64 rows 1-2)
    entrance = passage(pt(4010.0, 287.0, 668.2), pt(306.0, 286.4, -1181.0))
    # ascending passage: its beginning to its end (rows 7-8)
    ascending = passage(pt(3016.3, 286.6, 179.9), pt(1626.8, 287.0, 852.6))
    # grand gallery: floor from the end of the ascending passage to its virtual south end (rows 8, 11); 82 in wide, 339 in high
    g0, g1 = pt(1626.8, 284.4, 852.6), pt(-61.7, 284.4, 1689.0)
    gallery = passage(g0, g1, bore=(82 * INCH, 339 * INCH))
    lift = 339 * INCH / 2 * np.cos(np.radians(26.3))
    gallery['centre'][2] = round(gallery['centre'][2] + lift, 2)
    # King's chamber: N wall 330.6 S, width 206.3 (to 537 S); E wall 305.0 E, length 412.6; floor 1692.8, roof 1922.6
    kc = {'type': 'box', 'centre': [round((305.0 - 412.6 / 2) * INCH, 2), round(-(330.6 + 206.3 / 2) * INCH, 2),
                                    round((1692.8 + 1922.6) / 2 * INCH, 2)],
          'size': [round(412.6 * INCH, 2), round(206.3 * INCH, 2), round((1922.6 - 1692.8) * INCH, 2)]}
    # Queen's chamber: E wall 308 E, W end 72 E (length 226.47), ridge on the E-W centre line (width 205.85);
    # floor 834.4, ridge 1078.7 (gabled roof, drawn as a box to the ridge)
    qc = {'type': 'box', 'centre': [round((308 + 72) / 2 * INCH, 2), round(0.3 * INCH, 2), round((834.4 + 1078.7) / 2 * INCH, 2)],
          'size': [round(226.47 * INCH, 2), round(205.85 * INCH, 2), round((1078.7 - 834.4) * INCH, 2)]}
    # subterranean chamber: centre 203 S, 25.9 E, -1056; walls 300 E / 252 W (553 in), N-S about 325 in; height as the legacy table
    sc = {'type': 'box', 'centre': [round(25.9 * INCH, 2), round(-203 * INCH, 2), round(-1056 * INCH, 2)],
          'size': [round(553 * INCH, 2), round(325 * INCH, 2), 5.0]}
    rows = [('kings-chamber', "King's Chamber", 'chamber', kc), ('queens-chamber', "Queen's Chamber", 'chamber', qc),
            ('grand-gallery', 'Grand Gallery', 'passage', gallery), ('ascending-passage', 'Ascending Passage', 'passage', ascending),
            ('descending-passage', 'Descending Passage', 'passage', entrance),
            ('subterranean-chamber', 'Subterranean Chamber', 'chamber', sc)]
    return [{'id': f'khufu-{i}', 'name': f'{n} (Khufu)', 'kind': k, 'shape': sh, 'source': PETRIE} for i, n, k, sh in rows]


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
        if mon['id'] == 'khufu':
            for f in khufu_petrie():
                sh = dict(f['shape'])
                sh['centre'] = [round(float(e) + sh['centre'][0], 2), round(float(n) + sh['centre'][1], 2), round(base + sh['centre'][2], 2)]
                features.append({'id': f['id'], 'name': f['name'], 'kind': f['kind'], 'status': 'surveyed', 'fill': 'air',
                                 'shape': sh, 'source': f['source'],
                                 'note': 'From Petrie\'s tabulated positions; chambers drawn as boxes (the Queen\'s gabled roof to its ridge).',
                                 'placement': 'surveyed'})
        # Khafre from Petrie (1883), secs. 73-75, 78: the passage axis is 490.3 in E of the middle of the N face;
        # Belzoni's chamber is 557.9 x 195.8 in, 206.4 high to the walls, its door 104.3-144.9 in from its E wall,
        # so its E wall is about 614.9 in E and its W wall about 57 in E of the centre.
        petrie_khafre = {
            "Belzoni's Chamber": {'x': (614.9 - 557.9 / 2) * INCH, 'size': [557.9 * INCH, 195.8 * INCH, 206.4 * INCH],
                                  'note': "East-west position and size from Petrie (secs. 73-75); north-south and height from the first investigation's table."},
            'Upper Descending Passage': {'x': 490.3 * INCH, 'note': "Axis from Petrie (sec. 73); the rest from the first investigation's table."},
            'Lower Chamber': {'size': [411.9 * INCH, 123.1 * INCH, None], 'note': "Size from Petrie (sec. 78); position from the first investigation's table."},
        }
        for f in mon['features']:
            if mon['id'] == 'khufu' and not f['name'].startswith(('ScanPyramids', 'North Face')):
                continue                                     # replaced by Petrie above
            f = dict(f)
            ov = petrie_khafre.get(f['name']) if mon['id'] == 'khafre' else None
            if ov:
                f['centre'] = list(f['centre'])
                if 'x' in ov:
                    f['centre'][0] = ov['x']
                if 'size' in ov:
                    f['size'] = [round(v, 2) if v is not None else f['size'][i] for i, v in enumerate(ov['size'])]
                f['source'] = f['source'] + '; Petrie (1883), ch. 9'
                f['note'] = ov['note']
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
            "Khufu's chambers and passages are placed from Petrie (1883, sec. 64); the ScanPyramids voids and all of Khafre's features are still the first investigation's simplified boxes, approximate in placement, and are next to be re-derived.",
            'Terrain is a DSM: monument footprints are replaced by the surrounding ground and the pyramids are drawn from their survey dimensions.',
        ],
        'sources': {'sharafeldin2019': 'doi:10.5194/gi-8-29-2019', 'copernicus': 'Copernicus GLO-30 DEM, ESA/Airbus, via copernicus-dem-30m.s3.amazonaws.com'},
    }
    (OUT / 'site.json').write_text(json.dumps(site, indent=1))
    print(f'wrote {OUT/"site.json"} ({len(features)} features) and terrain {crop["nx"]}x{crop["ny"]}')


if __name__ == '__main__':
    main()
