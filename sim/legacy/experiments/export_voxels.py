"""Export tomogram volumes, surface textures and georeferencing for the three-dimensional viewer.

One bundle per patch: a uint8 voxel volume, a surface amplitude texture, a terrain height patch
from the Copernicus DEM, and a JSON header carrying the georeferencing and the honest depth
conventions. The viewer reads nothing else, so adding a site means adding bundles.

Volumes are written in the viewer's axis order [depth][across][along] so a WebGL 3-D texture can
be uploaded without transposition.
"""
import json
import os
import sys
import time

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter

from real_common import run_patch, SITES, RESULTS, PRODUCTS
from sarsim.dwell import DwellProduct
from sarsim.tomo import nyquist_depth
from sarsim.terrain import DemTile

OUT = os.path.join(RESULTS, 'voxels')
DEM_PATH = os.environ.get('GIZA_DEM', '')
GEOID_M = 15.46          # EGM96 at Giza, from catalog/giza/landmarks.json
DOWN_AZ, DOWN_RG, N_DEPTH = 2, 2, 64

# Independently known voids, as a depth band rather than a position. Individual shafts cannot be
# placed in the image to better than their own width, but the band they occupy is well
# established from excavation, and it needs no registration to draw.
KNOWN_VOIDS = {
    ('giza', 'cemetery_east'): dict(label='Surveyed burial shafts, Eastern Cemetery',
                                    depth_min_m=5.0, depth_max_m=30.0,
                                    source='Reisner and Simpson mastaba surveys; hundreds of excavated shafts'),
    ('giza', 'cemetery_west'): dict(label='Surveyed burial shafts, Western Cemetery',
                                    depth_min_m=5.0, depth_max_m=30.0,
                                    source='Reisner Western Cemetery surveys; hundreds of excavated shafts'),
}
PERIODS = 2.0
# Monuments get a deeper axis so the surveyed chambers, which lie 100 to 180 m below the apex,
# fall inside the box. The volume simply repeats across that range, which is the comparison.
DEEP_PERIODS, DEEP_DEPTH = 8.0, 192
DEEP_PATCHES = {('giza', 'khufu'), ('giza', 'khafre')}
UNDERWORLD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          'catalog', 'giza', 'underworld.json')


def survey_features(site, key, run, product=None):
    """Surveyed chambers placed in the patch frame: horizontally by the product's own RPC,
    vertically as depth below the apex, which is what the radar surface is over a pyramid.

    Nothing here comes from the radar. Placement is good to a few metres at monument scale and
    is not a registration claim."""
    ident = {'khufu': 'khufu', 'khafre': 'khafre'}.get(key)
    if site != 'giza' or ident is None or not os.path.isfile(UNDERWORLD):
        return None
    cat = json.load(open(UNDERWORLD))
    mon = next((m for m in cat['monuments'] if m['id'] == ident), None)
    if mon is None:
        return None
    meta = run['meta']
    acq = meta['acquisition']
    origin = meta['crop_origin']
    step_az = acq['azimuth_spacing_m'] * meta['stride'][0] * DOWN_AZ
    step_rg = acq['range_spacing_m'] * meta['stride'][1] * DOWN_RG
    rows, cols = run['grid_rows'], run['grid_cols']
    centre_row = origin[0] + (rows[0] + rows[-1]) / 2
    centre_col = origin[1] + (cols[0] + cols[-1]) / 2
    metres_per_deg_lat = 111132.0
    metres_per_deg_lon = 111320.0 * np.cos(np.deg2rad(mon['latitude']))
    out = []
    for f in mon['features']:
        east, north, up = f['centre']
        lat = mon['latitude'] + north / metres_per_deg_lat
        lon = mon['longitude'] + east / metres_per_deg_lon
        height = mon['base_elevation_m'] + up + GEOID_M
        try:
            r, c = product.geolocate(lat, lon, height)
        except Exception:
            continue
        out.append({
            'name': f['name'], 'kind': f['kind'], 'source': f.get('source', ''),
            'note': f.get('note', ''),
            'centre_m': [float((r - centre_row) * acq['azimuth_spacing_m']),
                         float((c - centre_col) * acq['range_spacing_m'])],
            'depth_below_apex_m': float(mon['height_m'] - up),
            'size_m': [float(f['size'][0]), float(f['size'][1]), float(f['size'][2])],
        })
    return {'monument': mon['name'], 'reference': 'depth measured below the apex; horizontal '
            'placement by the product RPC, not fitted to the radar',
            'features': out, 'references': cat['references']}


def volume(run, periods=PERIODS, n_depth=N_DEPTH):
    n_az, n_rg = run['grid_rows'].size, run['grid_cols'].size
    q = run['q'].astype(float).reshape(n_az, n_rg, -1, 2)
    kz = run['kz']
    period = 2.0 * nyquist_depth(np.sort(kz))
    z = np.linspace(period / n_depth, periods * period, n_depth)
    A = np.exp(-1j * np.outer(kz, z))
    flat = q.reshape(-1, q.shape[2], 2)
    power = np.abs((flat[..., 0] + 1j * flat[..., 1]) @ A) ** 2
    vol = np.log10(power.reshape(n_az, n_rg, z.size) + 1e-20)
    vol = gaussian_filter(vol, (1.0, 1.0, 1.2))
    return z, period, vol[::DOWN_AZ, ::DOWN_RG]


def to_uint8(vol):
    lo, hi = np.percentile(vol, [2.0, 99.5])
    scaled = np.clip((vol - lo) / max(hi - lo, 1e-9), 0, 1)
    return (scaled * 255).astype(np.uint8), float(lo), float(hi)


def surface_png(run, path):
    n_az, n_rg = run['grid_rows'].size, run['grid_cols'].size
    amp = run['amplitude'].astype(float).reshape(n_az, n_rg)[::DOWN_AZ, ::DOWN_RG]
    db = 20 * np.log10(np.maximum(amp, 1e-3))
    lo, hi = np.percentile(db, [3, 99.5])
    img = np.clip((db - lo) / max(hi - lo, 1e-9), 0, 1)
    Image.fromarray((img.T * 255).astype(np.uint8), mode='L').save(path, optimize=True)
    return float(lo), float(hi)


def export(site, key, dem=None, product=None):
    run = run_patch(site, key, 'paper', verbose=False)
    meta = run['meta']
    acq = meta['acquisition']
    deep = (site, key) in DEEP_PATCHES
    z, period, vol = volume(run, periods=DEEP_PERIODS if deep else PERIODS,
                            n_depth=DEEP_DEPTH if deep else N_DEPTH)
    data, lo, hi = to_uint8(vol)
    d = os.path.join(OUT, site)
    os.makedirs(d, exist_ok=True)
    # viewer order: depth slowest, then across-range, then along-azimuth
    np.transpose(data, (2, 1, 0)).tofile(os.path.join(d, f'{key}.bin'))
    tex_lo, tex_hi = surface_png(run, os.path.join(d, f'{key}.png'))
    nx, ny, nz = data.shape
    ground_az = acq['azimuth_spacing_m'] * meta['stride'][0] * DOWN_AZ
    # Slant range is kept as measured. Converting to ground range with one incidence angle is
    # wrong over a laid-over monument, where the axis mixes ground distance and height.
    slant_rg = acq['range_spacing_m'] * meta['stride'][1] * DOWN_RG
    header = {
        'site': site, 'patch': key, 'label': meta['label'], 'kind': meta['kind'],
        'shape': [int(nx), int(ny), int(nz)],
        'axis_order': 'binary is [depth][across_range][along_azimuth], uint8',
        'spacing_m': {'along_azimuth': float(ground_az), 'across_slant_range': float(slant_rg),
                      'depth': float(z[1] - z[0])},
        'extent_m': {'along_azimuth': float(nx * ground_az), 'across_slant_range': float(ny * slant_rg),
                     'depth': float(z[-1])},
        'centre_lonlat': [meta['longitude'], meta['latitude']],
        'assumed_height_m': meta['height_m'],
        'terrain_height_m': None if dem is None else float(dem.orthometric(meta['latitude'], meta['longitude'])),
        'depth_axis': {'lambda_s_m': meta['lam_s_m'],
                       'repeat_period_m': float(period),
                       'real_fit_ambiguity_limit_m': float(period / 2),
                       'note': ('Depth is the oscillation rate of the registration trajectory relabelled through '
                                'a declared sound wavelength. The volume repeats exactly every repeat_period_m.')},
        'known_voids': KNOWN_VOIDS.get((site, key)),
        'survey': survey_features(site, key, run, product),
        'value_range_log10': [lo, hi], 'surface_db_range': [tex_lo, tex_hi],
        'acquisition': {k: acq[k] for k in ('satellite', 'product_type', 'collection_start',
                                            'wavelength_m', 'incidence_center_deg', 'polarization')},
        'bytes': int(data.size),
    }
    with open(os.path.join(d, f'{key}.json'), 'w') as f:
        json.dump(header, f, indent=2)
    return header


def main(argv):
    t0 = time.time()
    dem = DemTile(DEM_PATH) if DEM_PATH and os.path.isfile(DEM_PATH) else None
    if dem is None:
        print('  no DEM supplied (set GIZA_DEM); terrain heights omitted', flush=True)
    wanted = argv or ['giza']
    index = []
    for site in wanted:
        for entry in SITES[site]:
            try:
                h = export(site, entry[0], dem if site == 'giza' else None,
                           product=DwellProduct(PRODUCTS[site]))
            except (FileNotFoundError, KeyError, StopIteration):
                print(f'  {site}/{entry[0]}: no cached run, skipped', flush=True)
                continue
            index.append({k: h[k] for k in ('site', 'patch', 'label', 'kind', 'shape', 'extent_m', 'bytes')})
            short = (h['label'].replace(' pyramid', '').replace(' and temples', '')
                     .replace('Open plateau, ', 'Plateau ').replace('Open hillside, ', 'Hillside '))
            index[-1]['short_label'] = short[0].upper() + short[1:]
            index[-1]['known_voids'] = bool(h.get('known_voids'))
            index[-1]['survey_features'] = len((h.get('survey') or {}).get('features', []))
            print(f"  {site}/{entry[0]}: {h['shape']} = {h['bytes'] / 1e6:.1f} MB", flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'index.json'), 'w') as f:
        json.dump({'bundles': index,
                   'sites': {'giza': 'Giza plateau, Egypt', 'sacsayhuaman': 'Sacsayhuaman, Peru'},
                   'generated': time.strftime('%Y-%m-%d'),
                   'geoid_undulation_m': GEOID_M}, f, indent=2)
    print(f'wrote {len(index)} bundles in {time.time() - t0:.0f}s -> {OUT}')


if __name__ == '__main__':
    main(sys.argv[1:])
