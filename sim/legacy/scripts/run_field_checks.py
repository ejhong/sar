"""Publish bounded image/measurement checks from an explicitly supplied SLC.

Example: python scripts/run_field_checks.py --site giza --input /path/scene.h5
Only derived display images, metrics and provenance enter results/field/.
Originals are read-only. No large intermediate arrays are saved anywhere.
"""
import argparse
import csv
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import h5py
import numpy as np
from PIL import Image

from fieldwork.processing import (bilinear_geoid, check_native, full_row_spectrum,
                                  intensity_preview, read_complex, timing_audit,
                                  translation_trials)
from sarsim.geolocation import project_rpc
from fieldwork.figures import spectrum_figure


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def metadata(file):
    names = ('acquisition_start_utc', 'acquisition_end_utc', 'zerodoppler_start_utc',
             'zerodoppler_end_utc', 'acquisition_prf', 'processing_prf',
             'azimuth_time_interval', 'azimuth_ground_spacing', 'slant_range_spacing',
             'incidence_center', 'number_of_azimuth_samples', 'number_of_range_samples',
             'doppler_rate_coeffs', 'first_pixel_time', 'range_sampling_rate',
             'total_processed_bandwidth_azimuth', 'carrier_frequency', 'avg_scene_height',
             'processor_version', 'product_type', 'polarization')
    result = {}
    for name in names:
        value = file[name][()]
        if isinstance(value, bytes):
            value = value.decode('utf8')
        elif isinstance(value, np.ndarray):
            value = value.tolist()
        elif isinstance(value, np.generic):
            value = value.item()
        result[name] = value
    return result


def render_preview(power, view, meta, output):
    if not np.any(power > 0):
        raise ValueError('Empty intensity preview')
    db = 10*np.log10(np.maximum(power, np.finfo(np.float32).tiny))
    low, high = np.percentile(db, [3, 99.5])
    if high <= low:
        raise ValueError('Preview has no display contrast')
    grey = np.uint8(np.clip((db-low)/(high-low), 0, 1)*255)
    if view['flip_rows']:
        grey = grey[::-1]
    rows, cols = view['rows'], view['cols']
    ground_range = meta['slant_range_spacing']/np.sin(np.deg2rad(meta['incidence_center']))
    aspect = (cols[1]-cols[0])*ground_range/((rows[1]-rows[0])*meta['azimuth_ground_spacing'])
    width = min(power.shape[1], 1400)
    height = round(width/aspect)
    if height > 1400:
        width = round(width*1400/height)
        height = 1400
    im = Image.fromarray(grey).resize((width, height), Image.Resampling.LANCZOS)
    im.save(output, 'WEBP', quality=92, method=6)
    return dict(image=output.name, width=width, height=height, display_db_limits=[float(low), float(high)],
                approximate_ground_width_m=(cols[1]-cols[0])*ground_range,
                approximate_azimuth_height_m=(rows[1]-rows[0])*meta['azimuth_ground_spacing'],
                intensity_block=view['block'], calibrated_backscatter=False)


def survey_projection(site, audit, meta):
    path = ROOT/'catalog'/site/'landmarks.json'
    if not path.exists():
        return None
    survey = json.loads(path.read_text())
    rpc = audit['source_metadata']['rpc']
    points = []
    for p in survey['points']:
        lat, lon, H = p['latitude'], p['longitude'], p['orthometric_height_m']
        N = bilinear_geoid(survey['geoid'], lat, lon)
        correct = project_rpc(rpc, lat, lon, H+N)
        wrong = project_rpc(rpc, lat, lon, H)
        sensitivity = [(project_rpc(rpc, lat, lon, H+N+dh)-correct).tolist() for dh in [-5, 5]]
        points.append(dict(**p, geoid_undulation_m=N, provisional_ellipsoidal_height_m=H+N,
                           projected_native_row_col=correct.tolist(),
                           using_msl_as_ellipsoid_row_col=wrong.tolist(),
                           height_minus_plus_5m_shifts_px=sensitivity))
    shifts = [p['projected_native_row_col'][1]-p['using_msl_as_ellipsoid_row_col'][1] for p in points]
    return dict(source=survey['source'], reference=survey['reference'], points=points,
                mean_datum_range_shift_px=float(np.mean(shifts)),
                approximate_mean_datum_ground_shift_m=float(np.mean(shifts)*meta['slant_range_spacing']/np.sin(np.deg2rad(meta['incidence_center']))),
                observed_landmark_count=0, held_out_rmse_m=None, native_registration_validated=False,
                interpretation='Independent survey coordinates projected without a fitted correction. No radar correspondences accepted: old casing/ground corners are not unambiguous point reflectors. GPMP vertical datum tie, terrain and layover remain unresolved. The ±5 m envelope is a sensitivity scenario, not a confidence interval.')


def run(site, source):
    catalog = json.loads((ROOT/'catalog/sites.json').read_text())
    record = next(s for s in catalog['sites'] if s['id'] == site)
    matches = [a for a in record['acquisitions'] if a['filename'] == source.name]
    if len(matches) != 1:
        raise ValueError('Input filename does not match the explicit acquisition record')
    acquisition = matches[0]
    design_path = ROOT/'catalog'/acquisition['design']
    design = json.loads(design_path.read_text())
    out = (ROOT/'results/field'/acquisition['report']).parent
    out.mkdir(parents=True, exist_ok=True)
    audit_path = ROOT/'research'/acquisition['geometry_audit']
    audit = json.loads(audit_path.read_text())
    if source.stat().st_size != audit['file_bytes']:
        raise ValueError('Input size differs from the inspected acquisition')
    with h5py.File(source, 'r') as file:
        shape = check_native(file)
        if list(shape) != audit['source_metadata']['shape']:
            raise ValueError('Input shape differs from the geolocation audit')
        for key, value in audit['source_metadata']['rpc'].items():
            if not np.array_equal(file['RPC'][key][()],value):
                raise ValueError('Input RPC differs from the geolocation audit')
        meta = metadata(file)
        views = []
        for view in design['views']:
            print(f'{site}: reading {view["id"]} intensity in bounded chunks', flush=True)
            power, hashes = intensity_preview(file['s_i'], file['s_q'], view['rows'], view['cols'], view['block'])
            display = render_preview(power, view, meta, out/f'{view["id"]}.webp')
            views.append(dict(**view, **display, **hashes))
            del power
        timing = timing_audit(meta)
        spectra = []
        for strip in design['spectrum_strips']:
            print(f'{site}: full-row Fourier support at columns {strip["cols"]}', flush=True)
            samples = read_complex(file, strip['rows'], strip['cols'])
            spectrum = full_row_spectrum(samples, meta['processing_prf'], meta['total_processed_bandwidth_azimuth'])
            spectra.append(dict(**strip, **spectrum, complex_samples_sha256=hashlib.sha256(samples.astype('<c8').tobytes()).hexdigest()))
            del samples
        control = design['translation']
        samples = read_complex(file, control['rows'], control['cols'])
        trials = translation_trials(samples, control['shifts_px'], control['patch'], control['upsample'])
        groups = []
        for dr, dc in control['shifts_px']:
            selection = [r for r in trials if r['injected_row_px'] == dr and r['injected_col_px'] == dc]
            errors = [r['error_px'] for r in selection]
            groups.append(dict(injected_row_px=dr, injected_col_px=dc, patches=len(selection),
                               median_measured_row_px=float(np.median([r['measured_row_px'] for r in selection])),
                               median_measured_col_px=float(np.median([r['measured_col_px'] for r in selection])),
                               median_error_px=float(np.median(errors)), p95_error_px=float(np.percentile(errors,95)),
                               max_error_px=max(errors), meets_0_02px_p95_tolerance=bool(np.percentile(errors,95) <= control['max_p95_error_px'])))
        translation = dict(control, groups=groups, trials_csv='translation_trials.csv',
                           complex_samples_sha256=hashlib.sha256(samples.astype('<c8').tobytes()).hexdigest(),
                           interpretation='Known translations of real complex texture. Patches from one crop are dependent; no inferential confidence interval or micromotion claim. The declared 0.02 px tolerance is an engineering check, not a physical sensitivity limit.')
    with (out/'translation_trials.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(trials[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(trials)
    spectrum_figure(spectra,timing,out/'spectrum.svg')
    spectrum_figure(spectra,timing,out/'spectrum-mobile.svg',compact=True)
    dependencies = [Path('scripts/run_field_checks.py'), Path('fieldwork/processing.py'),
                    Path('fieldwork/figures.py'),
                    Path('sarsim/track.py'), Path('sarsim/geolocation.py'),
                    Path('catalog')/acquisition['design'],
                    audit_path.relative_to(ROOT)]
    if (ROOT/'catalog'/site/'landmarks.json').exists():
        dependencies.append(Path(f'catalog/{site}/landmarks.json'))
    report = dict(version=1, site=site, acquisition=acquisition, design_sha256=digest(design_path),
                  versions=dict(python=platform.python_version(),**{name:version(name) for name in ['numpy','scipy','h5py','Pillow','matplotlib']}),
                  file_bytes=source.stat().st_size, native_shape=list(shape), source_metadata=meta,
                  source_metadata_sha256=hashlib.sha256(json.dumps(meta,sort_keys=True).encode()).hexdigest(),
                  hash_scope='Selected metadata and original selected sample values; not the HDF5 container bytes. Full overview checks I/Q finiteness and hashes both complete arrays in bounded chunks.',
                  source_sha256={str(p):digest(ROOT/p) for p in dependencies},
                  output_sha256={p.name:digest(p) for p in out.iterdir() if p.suffix in ('.webp','.csv','.svg')},
                  figures=dict(spectrum='spectrum.svg',spectrum_mobile='spectrum-mobile.svg'),
                  views=views, timing=timing, spectra=spectra, translation=translation,
                  survey_projection=survey_projection(site, audit, meta),
                  field_depth_results=0, native_registration_validated=False, aperture_time_validated=False,
                  interpretation='Real-image quality and measurement controls only. Survey models are external references. No underground detection, depth estimate or absence claim is produced.')
    write_json(ROOT/'results/field'/acquisition['report'], report)
    print(json.dumps(dict(site=site, collection_s=timing['collection_duration_s'],
                          translation_groups_passing=sum(g['meets_0_02px_p95_tolerance'] for g in groups),
                          translation_groups_total=len(groups),
                          spectrum_outside_fractions=[s['outside_zero_centered_metadata_band_fraction'] for s in spectra]), indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sites = json.loads((ROOT/'catalog/sites.json').read_text())['sites']
    parser.add_argument('--site', choices=[s['id'] for s in sites], required=True)
    parser.add_argument('--input', type=Path, required=True)
    args = parser.parse_args()
    run(args.site, args.input)
