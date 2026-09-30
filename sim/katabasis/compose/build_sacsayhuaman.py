"""Regenerate sites/sacsayhuaman/site.json and terrain.json from their sources.

    uv run python -m katabasis.compose.build_sacsayhuaman

A second real site for the satellite methods: the Inca walls of Sacsayhuamán above Cusco, the Rodadero outcrop beside
them, and the city's houses below, all inside one ICEYE Spotlight Dwell Fine pass (22 August 2025). Nothing is claimed
beneath it; the site is composed so that what each method draws over megalithic walls can be set beside what it draws
over houses and open hillside.

Frame: the origin at Wikipedia's coordinates for Sacsayhuamán. Places: OpenStreetMap objects (ids below, read through
the Overpass API on 29 September 2026). Terrain: Copernicus GLO-30 DSM (buildings and walls included, blurred at 30 m).
Geoid: EGM2008, the model the DSM's heights are in, through PROJ.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .dem import local_crop, write_crop
from .geo import to_local

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'sites' / 'sacsayhuaman'
ORIGIN = (-13.50778, -71.98222)          # Wikipedia, Sacsayhuamán: "13°30′28″S 71°58′56″W / 13.50778°S 71.98222°W"
EXTENT = {'x': [-700.0, 1000.0], 'y': [-1250.0, 700.0], 'z': [3150.0, 3750.0]}
GEOID_M = 46.2                           # EGM2008 undulation at the origin (PROJ 9.5.1, grid us_nga_egm08_25): 46.19 m
OSM = '© OpenStreetMap contributors (ODbL), read through the Overpass API on 29 September 2026'
# (id, name, latitude, longitude, OpenStreetMap object, what it is)
PLACES = [
    ('walls', 'The zigzag walls', -13.508618, -71.982085, 'ways 431479258 and 431479257 ("Murallas de Sacsayhuamán", '
     'barrier=retaining_wall), the midpoint of their extent', 'three tiers of Inca terrace walls; "The longest of the three '
     'walls is about 400 meters. They are about 6 meters tall" (Wikipedia)'),
    ('muyuqmarka', 'Muyuqmarka', -13.509012, -71.982802, 'way 480267375 (historic=castle), its centre',
     'the round tower\'s foundations on the hilltop behind the walls'),
    ('rodadero', 'Rodadero (Suchuna)', -13.506724, -71.981948, 'node 5857164644 ("Suchuna")',
     'the outcrop across the esplanade from the walls, with carved seats ("Trono del Inka", node 5857165093)'),
    ('chinkana', 'Chinkana Chica', -13.505862, -71.980801, 'node 5857164643 ("Chinkana Chica")',
     'a short passage through the rock north of the Rodadero'),
    ('cristo', 'Cristo Blanco', -13.509566, -71.978088, 'node 4746310711 ("Cristo Blanco")',
     'the statue on Pukamuqu hill, east of the walls'),
    ('san-blas', 'San Blas church', -13.514902, -71.974311, 'way 23799145 ("San Blas", amenity=place_of_worship)',
     'the church at the heart of the San Blas quarter\'s houses'),
    ('plaza', 'Plaza de Armas', -13.516767, -71.978779, 'node 6266294785 ("Plaza de Armas")',
     'the city\'s main square, colonial buildings on Inca foundations round it'),
]


def main() -> None:
    crop = local_crop(*ORIGIN, tuple(EXTENT['x']), tuple(EXTENT['y']), 10.0)
    OUT.mkdir(parents=True, exist_ok=True)
    write_crop(OUT / 'terrain.json', crop)
    zs = np.asarray(crop['z'])

    def ground(x, y):
        ix = int(round((x - crop['x0']) / crop['dx']))
        iy = int(round((y - crop['y0']) / crop['dx']))
        return round(float(zs[ix, iy]), 1)

    landmarks = []
    for pid, name, lat, lon, obj, what in PLACES:
        e, n, _ = to_local(lat, lon, *ORIGIN)
        landmarks.append({'id': pid, 'name': name, 'position': [round(float(e), 1), round(float(n), 1), ground(e, n)],
                          'what': what, 'source': f'OpenStreetMap {obj}; {OSM}'})
    site = {
        'id': 'sacsayhuaman', 'name': 'Sacsayhuamán', 'kind': 'real',
        'summary': 'The Inca walls above Cusco, the Rodadero outcrop and the city\'s houses, in one satellite pass: '
                   'what each method draws over megalithic walls, beside houses and open fields.',
        'frame': {'axes': 'x east, y north, z up; metres',
                  'origin': {'latitude': ORIGIN[0], 'longitude': ORIGIN[1], 'label': 'Sacsayhuamán (Wikipedia\'s coordinates)'},
                  'vertical': 'z = metres above EGM2008 (orthometric)',
                  'geoid_m': GEOID_M,
                  'geoid_source': 'EGM2008 undulation at the origin, 46.19 m, through PROJ 9.5.1 (us_nga_egm08_25) and pyproj 3.7.2'},
        'extent': EXTENT,
        'terrain': {'type': 'dem', 'file': 'terrain.json', 'datum_m': 0.0},
        'cover': [],
        'strata': [{'name': 'Limestone (assumed)', 'material': 'limestone', 'top': 'surface'}],
        'structures': [],
        'features': [],
        'landmarks': landmarks,
        'notes': [
            'No structure is claimed beneath this site. It is composed for the satellite methods, to set megalithic walls '
            'beside houses and open hillside in one image.',
            'The rock beneath is assumed limestone, as the walls\' largest blocks are ("the largest Limestone block", '
            'Wikipedia); it is not surveyed, and no run here uses it.',
            'Terrain is a DSM: the walls, houses and trees are in it, blurred over its 30 m cells.',
        ],
        'sources': {
            'wikipedia': 'Wikipedia, Sacsayhuamán: coordinates "13°30′28″S 71°58′56″W"; altitude "3,701 metres"; "The longest '
                         'of the three walls is about 400 meters. They are about 6 meters tall."',
            'osm': f'OpenStreetMap: {OSM}',
            'copernicus': 'Copernicus GLO-30 DEM, ESA/Airbus, via copernicus-dem-30m.s3.amazonaws.com',
            'egm2008': 'EGM2008 geoid (NGA), grid us_nga_egm08_25 from cdn.proj.org, evaluated with PROJ 9.5.1',
        },
    }
    (OUT / 'site.json').write_text(json.dumps(site, indent=1, ensure_ascii=False))
    print(f'wrote {OUT / "site.json"} ({len(landmarks)} places) and terrain {crop["nx"]}x{crop["ny"]}; '
          f'ground {float(zs.min()):.0f} to {float(zs.max()):.0f} m')


if __name__ == '__main__':
    main()
