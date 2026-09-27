"""Bounded native-SLC measurements. None of these functions estimates depth.

Rows and columns always mean the original HDF5 arrays. Display averaging never
feeds phase analysis. Full-row spectra retain azimuth support, but their Fourier
axis is an image frequency, not a validated observation-time axis.
"""
from datetime import datetime
import hashlib

import numpy as np
from scipy import fft


def bounds(shape, rows, cols):
    for size, pair in zip(shape, (rows, cols)):
        if len(pair) != 2 or any(int(x) != x for x in pair) or not 0 <= pair[0] < pair[1] <= size:
            raise ValueError('Explicit, nonempty native bounds must lie inside the raster')
    return tuple(map(int, rows)), tuple(map(int, cols))


def check_native(file):
    if file['s_i'].shape != file['s_q'].shape or len(file['s_i'].shape) != 2:
        raise ValueError('I/Q arrays must have identical two-dimensional shapes')
    if file['data_orientation'][()].decode() != 'native':
        raise ValueError('Only native ICEYE arrays are supported')
    shape = tuple(int(file[k][()]) for k in ('number_of_azimuth_samples', 'number_of_range_samples'))
    if shape != file['s_i'].shape:
        raise ValueError('Raster dimensions disagree with native metadata')
    return shape


def intensity_preview(i_data, q_data, rows, cols, block, max_chunk_bytes=64*1024**2):
    """Incoherent block mean, including partial edge blocks, in bounded chunks.

    Hashes identify the selected original float samples, in little-endian float32
    row order, separately for I and Q. These are not whole-file checksums.
    """
    if i_data.shape != q_data.shape:
        raise ValueError('I/Q shapes disagree')
    rows, cols = bounds(i_data.shape, rows, cols)
    br, bc = block
    if any(int(x) != x or x <= 0 for x in block):
        raise ValueError('Averaging blocks must be positive integers')
    br, bc = int(br), int(bc)
    nr, nc = rows[1]-rows[0], cols[1]-cols[0]
    # Up to five float32-sized arrays/scratch buffers; no full-raster read.
    bytes_per_row = nc*4*5
    if bytes_per_row*br > max_chunk_bytes:
        raise ValueError('A single averaging block exceeds the chunk memory budget')
    chunk_rows = max(br, (max_chunk_bytes//bytes_per_row//br)*br)
    output = np.empty(((nr+br-1)//br, (nc+bc-1)//bc), np.float32)
    hashes = [hashlib.sha256(), hashlib.sha256()]
    col_starts = np.arange(0, nc, bc)
    counts = np.minimum(bc, nc-col_starts)
    for start in range(rows[0], rows[1], chunk_rows):
        stop = min(start+chunk_rows, rows[1])
        real = np.asarray(i_data[start:stop, cols[0]:cols[1]], dtype='<f4')
        imag = np.asarray(q_data[start:stop, cols[0]:cols[1]], dtype='<f4')
        if not np.isfinite(real).all() or not np.isfinite(imag).all():
            raise ValueError('Nonfinite complex samples in preview region')
        for digest, values in zip(hashes, (real, imag)):
            digest.update(values.tobytes())
        power = real*real+imag*imag
        if not np.isfinite(power).all():
            raise ValueError('Intensity overflow in preview region')
        row_starts = np.arange(0, len(power), br)
        summed = np.add.reduceat(np.add.reduceat(power, row_starts, axis=0,
                                               dtype=np.float64), col_starts, axis=1)
        averaged = summed/(np.minimum(br, len(power)-row_starts)[:, None]*counts)
        output[(start-rows[0])//br:(start-rows[0])//br+len(averaged)] = averaged
    return output, dict(zip(('i_samples_sha256', 'q_samples_sha256'), (h.hexdigest() for h in hashes)))


def read_complex(file, rows, cols, max_bytes=128*1024**2):
    shape = check_native(file)
    rows, cols = bounds(shape, rows, cols)
    # Both components, assembled complex samples and downstream scratch allowance.
    if (rows[1]-rows[0])*(cols[1]-cols[0])*32 > max_bytes:
        raise ValueError('Complex selection exceeds working-memory budget')
    real = np.asarray(file['s_i'][rows[0]:rows[1], cols[0]:cols[1]], np.float32)
    imag = np.asarray(file['s_q'][rows[0]:rows[1], cols[0]:cols[1]], np.float32)
    if not np.isfinite(real).all() or not np.isfinite(imag).all():
        raise ValueError('Nonfinite complex samples')
    return real+1j*imag


def doppler_rate(coefficients, first_pixel_time, sampling_rate, range_samples, columns):
    """ICEYE/SarPy convention: polynomial argument is time relative to mid range.

    See NGA sarpy/io/complex/iceye.py calculate_drate_sf_poly and zd_ref_time.
    Return the supplied sign, with no simulated default or sign repair.
    """
    coefficients = np.asarray(coefficients, float)
    if coefficients.ndim != 1 or not np.isfinite(coefficients).all() or not 1 <= len(coefficients) <= 8:
        raise ValueError('Invalid Doppler-rate polynomial')
    columns = np.asarray(columns, float)
    if (not np.isfinite([first_pixel_time, sampling_rate, range_samples]).all()
            or sampling_rate <= 0 or range_samples <= 0 or not np.isfinite(columns).all()
            or np.any(columns < 0) or np.any(columns >= range_samples)):
        raise ValueError('Invalid range-time coordinates')
    # Subtract analytically to avoid cancellation of the two absolute times.
    return np.polynomial.polynomial.polyval((columns-range_samples/2)/sampling_rate, coefficients)


def timing_audit(metadata):
    elapsed = lambda start, end: (datetime.fromisoformat(metadata[end])-datetime.fromisoformat(metadata[start])).total_seconds()
    collection = elapsed('acquisition_start_utc', 'acquisition_end_utc')
    image_span = elapsed('zerodoppler_start_utc', 'zerodoppler_end_utc')
    prf, dt = metadata['processing_prf'], metadata['azimuth_time_interval']
    if not np.isfinite([prf,dt,metadata['acquisition_prf'],metadata['total_processed_bandwidth_azimuth']]).all() or metadata['acquisition_prf'] <= 0:
        raise ValueError('Nonfinite or nonpositive acquisition sampling')
    if collection <= 0 or image_span <= 0 or prf <= 0 or dt <= 0 or abs(prf*dt-1) > 1e-6:
        raise ValueError('Inconsistent acquisition timing')
    cols = np.array([0, metadata['number_of_range_samples']/2,
                     metadata['number_of_range_samples']-1])
    rates = doppler_rate(metadata['doppler_rate_coeffs'], metadata['first_pixel_time'],
                         metadata['range_sampling_rate'], metadata['number_of_range_samples'], cols)
    if np.any(rates >= 0):
        raise ValueError('Unsupported Doppler-rate sign; resolve processor convention first')
    bandwidth = metadata['total_processed_bandwidth_azimuth']
    if not 0 < bandwidth <= prf:
        raise ValueError('Processed bandwidth must fit the output sampling rate')
    return dict(collection_duration_s=collection, zero_doppler_image_span_s=image_span,
                image_span_from_sampling_s=(metadata['number_of_azimuth_samples']-1)*dt,
                acquisition_prf_hz=metadata['acquisition_prf'], processing_prf_hz=prf,
                doppler_rate_near_center_far_hz_s=rates.tolist(),
                nominal_bandwidth_over_rate_s=(bandwidth/abs(rates)).tolist(),
                processed_bandwidth_hz=bandwidth,
                interpretation='Bandwidth/rate is a nominal local consistency check. Image-row time is not aperture time. A dwell-product frequency-to-time adapter is not yet validated.',
                aperture_time_validated=False)


def full_row_spectrum(samples, processing_prf, metadata_bandwidth, bins=600):
    """Averaged full-row image spectrum; no cropping, taper, deramp or time labels."""
    if (samples.ndim != 2 or not np.iscomplexobj(samples) or not np.isfinite(samples).all()
            or not np.isfinite([processing_prf,metadata_bandwidth]).all()
            or not 0 < metadata_bandwidth <= processing_prf):
        raise ValueError('Provide finite complex full-row samples and a positive PRF')
    frequencies = fft.fftshift(fft.fftfreq(samples.shape[0], 1/processing_prf))
    power = np.mean(np.abs(fft.fftshift(fft.fft(samples, axis=0), axes=0))**2, axis=1, dtype=np.float64)
    if power.sum() <= 0:
        raise ValueError('Empty spectrum')
    outside = float(power[np.abs(frequencies) > metadata_bandwidth/2].sum()/power.sum())
    groups = np.array_split(np.arange(len(power)), min(bins, len(power)))
    curve = [dict(image_frequency_hz=float(frequencies[g].mean()), mean_power=float(power[g].mean())) for g in groups]
    return dict(outside_zero_centered_metadata_band_fraction=outside, curve=curve,
                axis='Fourier frequency of full native image rows, scaled by processing PRF; not vibration frequency.',
                interpretation='A zero-centered metadata band is a diagnostic reference, not an independently validated aperture window.')


def translation_trials(samples, shifts, patch=32, upsample=1000):
    """Known image translations on actual complex texture, independent of depth.

    Fourier shifts act on the larger crop; local tracker windows exclude its
    periodic boundary. Finite windows and texture need not recover shifts exactly.
    """
    from scipy.ndimage import fourier_shift
    from sarsim.track import patch_shifts
    if min(samples.shape) < patch*4:
        raise ValueError('Translation control needs a larger support than its patches')
    axis = [np.linspace(patch, n-patch-1, 5, dtype=int) for n in samples.shape]
    rows, cols = np.meshgrid(*axis, indexing='ij')
    rows, cols = rows.ravel(), cols.ravel()
    source_spectrum = fft.fftn(samples)
    records = []
    for dr, dc in shifts:
        if max(abs(dr), abs(dc)) > patch/4:
            raise ValueError('Injected translation too large for the fixed local patches')
        moved = fft.ifftn(fourier_shift(source_spectrum, (dr, dc)))
        measured = np.column_stack(patch_shifts(samples, moved, rows, cols, patch=patch,
                                               upsample=upsample))
        for row, col, value in zip(rows, cols, measured):
            records.append(dict(local_row=int(row), local_col=int(col),
                                injected_row_px=float(dr), injected_col_px=float(dc),
                                measured_row_px=float(value[0]), measured_col_px=float(value[1]),
                                error_px=float(np.linalg.norm(value-[dr, dc]))))
    return records


def bilinear_geoid(tile, latitude, longitude):
    """Interpolate an explicitly georeferenced small EGM grid tile (pixel points)."""
    values = np.asarray(tile['values_m'], float)
    if (not np.isfinite([latitude,longitude,tile['west_deg'],tile['north_deg'],tile['spacing_deg']]).all()
            or tile['spacing_deg'] <= 0):
        raise ValueError('Invalid geoid coordinates or spacing')
    x = (longitude-tile['west_deg'])/tile['spacing_deg']
    y = (tile['north_deg']-latitude)/tile['spacing_deg']
    if values.ndim != 2 or not np.isfinite(values).all() or not 0 <= x < values.shape[1]-1 or not 0 <= y < values.shape[0]-1:
        raise ValueError('Coordinates outside the finite geoid tile')
    j, i = int(x), int(y)
    a, b = x-j, y-i
    return float((1-b)*((1-a)*values[i,j]+a*values[i,j+1])+
                 b*((1-a)*values[i+1,j]+a*values[i+1,j+1]))
