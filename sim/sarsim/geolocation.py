"""Strict, metadata-only ICEYE RPC diagnostics, not a depth geometry adapter.

ICEYE specifies the GeoTIFF RPC00B model. Inputs are WGS84 longitude/latitude
and ellipsoidal height, not height above local ground or orthometric elevation.
Outputs retain the continuous pixel-centre convention of the supplied RPC.
No fitted offset, pixel rounding or terrain correction is applied here.

References:
https://sar.iceye.com/6.0.0/productFormats/metadata/
https://gdal.org/en/stable/development/rfc/rfc22_rpc.html
"""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

import numpy as np


RPC_SCALARS = tuple(f'{axis}_{kind}' for axis in
                    ('LINE', 'SAMP', 'LAT', 'LONG', 'HEIGHT')
                    for kind in ('OFF', 'SCALE'))
RPC_COEFFICIENTS = tuple(f'{axis}_{kind}_COEFF' for axis in
                         ('LINE', 'SAMP') for kind in ('NUM', 'DEN'))
CONTROL_POINTS = ('coord_center', 'coord_first_near', 'coord_first_far',
                  'coord_last_near', 'coord_last_far')


def project_rpc(rpc, latitude, longitude, height_m):
    """Return continuous (line, sample) RPC coordinates; do not assume alignment.

    Scalar or broadcastable input arrays are accepted. Coefficients must be
    complete, finite RPC00B values. Signed image scales are preserved as supplied
    (some ICEYE products use a negative LINE_SCALE). Reject extrapolation beyond normalized
    +/-1 coordinates, and singular denominators. Bounds are a conservative
    extrapolation guard, not an accuracy guarantee or a scene-footprint test.
    """
    values = {}
    for name in RPC_SCALARS + RPC_COEFFICIENTS:
        if name not in rpc:
            raise ValueError(f'Missing RPC field: {name}')
        value = np.asarray(rpc[name], dtype=float)
        expected = (20,) if name in RPC_COEFFICIENTS else ()
        if value.shape != expected or not np.all(np.isfinite(value)):
            raise ValueError(f'Invalid RPC field: {name}')
        if name.endswith('_SCALE') and value == 0:
            raise ValueError(f'RPC scale must be nonzero: {name}')
        if name in ('LAT_SCALE', 'LONG_SCALE', 'HEIGHT_SCALE') and value < 0:
            raise ValueError(f'Geographic RPC scale must be positive: {name}')
        values[name] = value
    lat, lon, height = np.broadcast_arrays(latitude, longitude, height_m)
    if not all(np.all(np.isfinite(a)) for a in (lat, lon, height)):
        raise ValueError('Coordinates must be finite')
    if np.any(np.abs(lat) > 90) or np.any(np.abs(lon) > 180):
        raise ValueError('Coordinates must be WGS84 degrees')
    p = (lat-values['LAT_OFF'])/values['LAT_SCALE']
    l = (lon-values['LONG_OFF'])/values['LONG_SCALE']
    h = (height-values['HEIGHT_OFF'])/values['HEIGHT_SCALE']
    if any(np.any(np.abs(a) > 1) for a in (p, l, h)):
        raise ValueError('Coordinate outside normalized RPC domain; extrapolation refused')
    terms = np.stack((np.ones_like(l), l, p, h, l*p, l*h, p*h,
                      l*l, p*p, h*h, l*p*h, l**3, l*p*p, l*h*h,
                      l*l*p, p**3, p*h*h, l*l*h, p*p*h, h**3), axis=-1)
    projected = []
    for axis in ('LINE', 'SAMP'):
        denominator = terms @ values[f'{axis}_DEN_COEFF']
        # A relative test also catches nearly cancelling terms.
        magnitude = np.abs(terms) @ np.abs(values[f'{axis}_DEN_COEFF'])
        if np.any(np.abs(denominator) <= 1e-12*magnitude):
            raise ValueError(f'Singular RPC denominator: {axis}')
        coordinate = (values[f'{axis}_OFF'] + values[f'{axis}_SCALE'] *
                      (terms @ values[f'{axis}_NUM_COEFF']) / denominator)
        if not np.all(np.isfinite(coordinate)):
            raise ValueError(f'Nonfinite RPC projection: {axis}')
        projected.append(coordinate)
    return np.stack(projected, axis=-1)


def _small_dataset(group, name, shape):
    """Validate the schema before any read, including deliberately malformed files."""
    import h5py
    if name not in group or not isinstance(group[name], h5py.Dataset):
        raise ValueError(f'Missing metadata dataset: {name}')
    data = group[name]
    if data.shape != shape or data.dtype.kind not in 'iuf':
        raise ValueError(f'Invalid metadata dataset: {name}')
    value = np.asarray(data[()])
    if not np.all(np.isfinite(value)):
        raise ValueError(f'Nonfinite metadata: {name}')
    return value.tolist()


def audit_iceye_geometry(path, sicd_xml=None):
    """Read small numeric metadata only; leave native registration unvalidated.

    coord_* has [column, row, latitude, longitude] with one-based indices in
    ICEYE's examples. We report the residual after subtracting one from those
    indices to compare with standard zero-based RPC pixel centres. The supplied
    points have no per-point height here; using average scene height is only a
    consistency diagnostic, not independent ground-truth validation.
    """
    import h5py
    path = Path(path)
    with h5py.File(path, 'r') as file:
        for name in ('s_i', 's_q'):
            if name not in file or not isinstance(file[name], h5py.Dataset):
                raise ValueError('Expected native ICEYE s_i and s_q arrays')
        shape = file['s_i'].shape
        if len(shape) != 2 or shape != file['s_q'].shape:
            raise ValueError('I/Q raster shapes disagree or are not two dimensional')
        rows = _small_dataset(file, 'number_of_azimuth_samples', ())
        cols = _small_dataset(file, 'number_of_range_samples', ())
        if tuple(shape) != (rows, cols):
            raise ValueError('Native raster dimensions disagree with metadata')
        # Numeric dimensions alone cannot establish an axis convention.
        if 'data_orientation' not in file or file['data_orientation'].shape != ():
            raise ValueError('Missing native orientation metadata')
        orientation = file['data_orientation'][()]
        if isinstance(orientation, bytes):
            orientation = orientation.decode('ascii')
        if orientation != 'native':
            raise ValueError('Only native ICEYE orientation is supported')
        if 'geo_ref_system' not in file or file['geo_ref_system'].shape != ():
            raise ValueError('Missing WGS84 reference metadata')
        reference = file['geo_ref_system'][()]
        if isinstance(reference, bytes):
            reference = reference.decode('ascii')
        if reference != 'WGS84':
            raise ValueError('Only WGS84 ICEYE metadata is supported')
        if 'RPC' not in file or not isinstance(file['RPC'], h5py.Group):
            raise ValueError('Missing RPC metadata group')
        rpc = {name: _small_dataset(file['RPC'], name, ()) for name in RPC_SCALARS}
        rpc.update({name: _small_dataset(file['RPC'], name, (20,))
                    for name in RPC_COEFFICIENTS})
        points = {name: _small_dataset(file, name, (4,)) for name in CONTROL_POINTS}
        height = _small_dataset(file, 'avg_scene_height', ())
        incidence = _small_dataset(file, 'incidence_center', ())
        spacing = _small_dataset(file, 'slant_range_spacing', ())
        if not 0 < incidence < 90 or spacing <= 0:
            raise ValueError('Invalid incidence angle or slant-range spacing')
    metadata = dict(rpc=rpc, control_points=points, shape=list(shape),
                    avg_scene_height_m=height, incidence_center_deg=incidence,
                    slant_range_spacing_m=spacing, data_orientation=orientation,
                    geo_ref_system=reference)
    checks = []
    for name, (col, row, lat, lon) in points.items():
        projected = project_rpc(rpc, lat, lon, height)
        shifted = project_rpc(rpc, lat, lon, height+10)
        checks.append(dict(name=name, latitude=lat, longitude=lon,
                           metadata_row_col_one_based=[row, col],
                           rpc_line_sample=projected.tolist(),
                           residual_to_zero_based_metadata_px=(projected-[row-1, col-1]).tolist(),
                           height_plus_10m_shift_px=(shifted-projected).tolist()))
    center = checks[0]
    shifts = []
    for delta in (-10, -5, -1, 0, 1, 5, 10):
        projected = project_rpc(rpc, center['latitude'], center['longitude'], height+delta)
        shift = projected-np.array(center['rpc_line_sample'])
        shifts.append(dict(height_change_m=delta, line_shift_px=float(shift[0]),
                           sample_shift_px=float(shift[1]),
                           approximate_flat_ground_range_shift_m=float(
                               shift[1]*spacing/np.sin(np.deg2rad(incidence)))))
    report = dict(
        schema_version=1, filename=path.name, file_bytes=path.stat().st_size,
        raster_samples_read=0, native_registration_validated=False,
        depth_analysis_ready=False, source_metadata=metadata,
        signed_image_scales={key: rpc[key] for key in ('LINE_SCALE', 'SAMP_SCALE') if rpc[key] < 0},
        source_metadata_sha256=hashlib.sha256(json.dumps(metadata, sort_keys=True,
                                                        allow_nan=False).encode()).hexdigest(),
        hash_scope='Selected metadata only; not a checksum of the full raster.',
        coordinate_convention='Unadjusted RPC pixel-centre coordinates; coord_* indices converted from one based for residuals.',
        height_assumption='Average scene ellipsoidal height at every metadata point; per-point terrain height is unverified.',
        control_point_checks=checks, center_height_sensitivity=shifts,
        interpretation='Metadata consistency and projection sensitivity only. No landmark registration, motion measurement or subsurface detection.',
        next_requirements=[
            'Verify independent surface landmarks and their WGS84 ellipsoidal heights.',
            'Resolve residual offsets, terrain assumptions and pixel conventions without fitting to buried targets.',
            'Register survey plans and their depth datums, with uncertainty and acquisition-date cavity state.',
            'Validate aperture timing and support separately before any motion-to-depth calculation.'],
        references=['https://sar.iceye.com/6.0.0/productFormats/metadata/',
                    'https://gdal.org/en/stable/development/rfc/rfc22_rpc.html'])
    if sicd_xml is not None:
        sicd_xml = Path(sicd_xml)
        root = ET.fromstring(sicd_xml.read_bytes())
        def integer(name):
            node = root.find(f'./{{*}}ImageData/{{*}}{name}')
            if node is None:
                raise ValueError(f'Missing SICD ImageData/{name}')
            return int(node.text)
        sicd_shape = [integer('NumRows'), integer('NumCols')]
        report['companion_sicd'] = dict(filename=sicd_xml.name, shape=sicd_shape,
            sha256=hashlib.sha256(sicd_xml.read_bytes()).hexdigest(),
            transposed_shape_matches_native=sicd_shape[::-1] == list(shape),
            interpretation='Dimension check only; even equal dimensions would not establish a coordinate transform.')
    return report
