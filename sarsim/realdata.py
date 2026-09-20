"""Bounded, read-only inspection of local SAR products.

No default radar geometry, no implicit full-raster reads, no depth inference.
ICEYE HDF5 I/Q pairs and complex TIFFs are supported. SICD/NITF uses optional
sarpy. A preview crop is not a validated sub-aperture processing domain.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import h5py
import tifffile


def value(v):
    array = np.asarray(v)
    if array.size > 32:
        return {"shape": list(array.shape), "dtype": str(array.dtype)}
    if array.dtype.kind in 'SUO':
        def decode(x):
            return x.decode('utf-8', errors='replace') if isinstance(x, bytes) else str(x)
        return decode(array.item()) if array.size == 1 else [decode(x) for x in array.ravel()]
    if array.size == 1:
        v = array.item()
        return v if not isinstance(v, float) or np.isfinite(v) else None
    return [v if not isinstance(v, float) or np.isfinite(v) else None for v in array.ravel().tolist()]


def inspect_product(path):
    path = Path(path)
    report = {"filename": path.name, "file_bytes": path.stat().st_size,
              "scope": "metadata only; no raster loaded", "complex_samples": False,
              "depth_analysis_ready": False}
    if h5py.is_hdf5(path):
        metadata, datasets = {}, {}
        with h5py.File(path, 'r') as f:
            def visit(name, item):
                if isinstance(item, h5py.Dataset):
                    datasets[name] = {"shape": list(item.shape), "dtype": str(item.dtype)}
                    # Bound *all* reads, including string arrays and variable-length data.
                    if item.size <= 32 and item.dtype.kind not in 'OV' and item.nbytes <= 8192:
                        metadata[name] = value(item[()])
            f.visititems(visit)
            for key, item in f.attrs.items():
                if np.asarray(item).size <= 32:
                    metadata['@'+key] = value(item)
            if 's_i' in f and 's_q' in f:
                si, sq = f['s_i'], f['s_q']
                if si.shape != sq.shape or si.ndim != 2:
                    raise ValueError('I/Q arrays must have the same two-dimensional shape')
                report.update(format='HDF5 I/Q', shape=list(si.shape), complex_samples=True,
                              datasets=datasets, metadata=metadata, sample_storage=[str(si.dtype), str(sq.dtype)])
            else:
                report.update(format='HDF5 (unmapped arrays)', datasets=datasets, metadata=metadata)
    elif path.suffix.lower() in ('.tif', '.tiff'):
        with tifffile.TiffFile(path) as f:
            series = f.series[0]
            report.update(format='TIFF', shape=list(series.shape), dtype=str(series.dtype),
                          complex_samples=series.dtype.kind == 'c',
                          metadata={"axes": series.axes, "description": (f.pages[0].description or '')[:8192]})
    elif path.suffix.lower() in ('.ntf', '.nitf', '.sicd'):
        try:
            from sarpy.io.complex.converter import open_complex
        except ImportError:
            report.update(format='NITF/SICD (reader needed)', reader='optional sarpy package')
        else:
            reader = open_complex(str(path))
            try:
                report.update(format='SICD', shape=list(reader.get_data_size_as_tuple()[0]), complex_samples=True,
                              metadata=reader.get_sicds_as_tuple()[0].to_dict())
            finally:
                reader.close()
    else:
        raise ValueError('supported inspection: HDF5, TIFF, SICD/NITF')
    report['required_before_depth_processing'] = [
        'Verify native azimuth/range axes and SLC product type',
        'Verify acquisition mode, aperture timing, processed Doppler support and centroid',
        'Use state vectors, geolocation, incidence and slant range at the target',
        'Retain adequate aperture/spatial support; a preview crop is not a new acquisition',
        'Freeze target regions, independent flanks and matched controls before depth output',
        'Run injected-motion, static-surface, window-gate and registration checks',
        'Supply independent surveyed positives and negatives for any detection claim'
    ]
    return report


def crop_bounds(shape, rows, cols, max_bytes=64*1024*1024):
    if len(shape) != 2:
        raise ValueError('a two-dimensional native complex raster is required')
    if rows is None or cols is None:
        raise ValueError('explicit row and column bounds are required; full-image reads are never implicit')
    r0, r1 = rows
    c0, c1 = cols
    if not (0 <= r0 < r1 <= shape[0] and 0 <= c0 < c1 <= shape[1]):
        raise ValueError('crop bounds outside the native raster')
    # Reserve working memory for conversion as well as the complex output.
    if (r1-r0)*(c1-c0)*32 > max_bytes:
        raise ValueError('crop exceeds the configured working-memory budget')
    return slice(r0, r1), slice(c0, c1)


def read_crop(path, rows, cols, max_bytes=64*1024*1024):
    path = Path(path)
    report = inspect_product(path)
    if not report['complex_samples']:
        raise ValueError('complex samples were not verified; intensity images cannot supply phase')
    rr, cc = crop_bounds(report['shape'], rows, cols, max_bytes)
    if report['format'] == 'HDF5 I/Q':
        with h5py.File(path, 'r') as f:
            crop = f['s_i'][rr, cc].astype(np.float32)+1j*f['s_q'][rr, cc].astype(np.float32)
    elif report['format'] == 'TIFF':
        try:
            data = tifffile.memmap(path)
        except ValueError as error:
            raise ValueError('compressed TIFF needs a windowed reader such as rasterio; refusing a full-image decode') from error
        crop = np.asarray(data[rr, cc], dtype=np.complex64).copy()
    else:
        from sarpy.io.complex.converter import open_complex
        reader = open_complex(str(path))
        try:
            crop = np.asarray(reader[rr, cc], dtype=np.complex64)
        finally:
            reader.close()
    if not np.isfinite(crop).all():
        raise ValueError('the requested complex crop contains nonfinite samples')
    report['crop'] = {"rows": list(rows), "cols": list(cols), "shape": list(crop.shape),
                      "sha256_complex64": hashlib.sha256(np.ascontiguousarray(crop, dtype=np.complex64).tobytes()).hexdigest(),
                      "note": "This hashes the selected samples, not the complete acquisition."}
    return crop.astype(np.complex64), report
