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

from real_common import run_patch, SITES, RESULTS
from sarsim.tomo import nyquist_depth
from sarsim.terrain import DemTile

OUT = os.path.join(RESULTS, 'voxels')
DEM_PATH = os.environ.get('GIZA_DEM', '')
GEOID_M = 15.46          # EGM96 at Giza, from catalog/giza/landmarks.json
DOWN_AZ, DOWN_RG, N_DEPTH = 2, 1, 96
PERIODS = 2.0


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


def export(site, key, dem=None):
    run = run_patch(site, key, 'paper', verbose=False)
    meta = run['meta']
    acq = meta['acquisition']
    z, period, vol = volume(run)
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
                h = export(site, entry[0], dem if site == 'giza' else None)
            except (FileNotFoundError, KeyError, StopIteration):
                print(f'  {site}/{entry[0]}: no cached run, skipped', flush=True)
                continue
            index.append({k: h[k] for k in ('site', 'patch', 'label', 'kind', 'shape', 'extent_m', 'bytes')})
            short = (h['label'].replace(' pyramid', '').replace(' and temples', '')
                     .replace('Open plateau, ', 'Plateau ').replace('Open hillside, ', 'Hillside '))
            index[-1]['short_label'] = short[0].upper() + short[1:]
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
