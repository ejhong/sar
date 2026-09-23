"""Processing adapter for ICEYE spotlight-dwell products, and the Doppler sub-aperture bank.

A dwell product is the case where single-image "Doppler tomography" is least obviously
hopeless: the processed azimuth bandwidth corresponds to roughly 25 s of real aperture time,
so sub-apertures are separated by seconds, not milliseconds. This module makes that mapping
explicit and refuses to guess any acquisition fact.

Two facts govern everything here.

*The image azimuth axis is not time.* An ICEYE dwell SLC spans about 0.7 s of zero-Doppler
row time but 25 s of aperture. Slow time is carried by the Doppler axis through
f = f_dc + Ka (t - t_zd), with Ka the (negative) Doppler rate.

*The azimuth shift of a sub-aperture measures line-of-sight velocity.* A target moving at
v_los acquires a Doppler offset 2 v_los / lambda, which azimuth focusing maps to a position
offset of (R / V) v_los. Tracking sub-apertures against a common reference therefore recovers
a velocity time series, which is a physically defined quantity; tracking a reference against a
slightly shifted offset band, as the published protocol does, measures a short-lag difference
of that series instead.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

from .orbit import Ephemeris, lla_to_ecef, perpendicular_baselines, steering_wavenumbers

C_LIGHT = 299792458.0
_SCALARS = ('processing_prf', 'azimuth_time_interval', 'acquisition_prf', 'carrier_frequency',
            'total_processed_bandwidth_azimuth', 'number_of_azimuth_samples',
            'number_of_range_samples', 'slant_range_to_first_pixel', 'slant_range_spacing',
            'azimuth_ground_spacing', 'incidence_center', 'first_pixel_time',
            'range_sampling_rate', 'avg_scene_height', 'mean_earth_radius')
_TEXT = ('acquisition_start_utc', 'acquisition_end_utc', 'zerodoppler_start_utc',
         'zerodoppler_end_utc', 'product_type', 'polarization', 'satellite_name',
         'look_side', 'orbit_direction')
_RPC = tuple(f'{a}_{k}' for a in ('LINE', 'SAMP', 'LAT', 'LONG', 'HEIGHT') for k in ('OFF', 'SCALE')) \
     + tuple(f'{a}_{k}_COEFF' for a in ('LINE', 'SAMP') for k in ('NUM', 'DEN'))


def _utc(value):
    return datetime.fromisoformat(value.decode() if isinstance(value, bytes) else str(value))


def _text(value):
    return (value.decode() if isinstance(value, bytes) else str(value)).strip()


@dataclass(frozen=True)
class Acquisition:
    """Every acquisition fact the processing needs, read once from the product."""
    path: str
    satellite: str
    product_type: str
    polarization: str
    look_side: str
    orbit_direction: str
    shape: tuple
    prf_hz: float
    azimuth_time_interval_s: float
    bandwidth_hz: float
    wavelength_m: float
    azimuth_spacing_m: float
    range_spacing_m: float
    slant_range_first_m: float
    incidence_center_deg: float
    scene_height_m: float
    collection_start: str
    collection_end: str
    zero_doppler_start: str
    zero_doppler_end: str

    @property
    def collection_duration_s(self):
        return (_utc(self.collection_end) - _utc(self.collection_start)).total_seconds()

    @property
    def image_span_s(self):
        return (_utc(self.zero_doppler_end) - _utc(self.zero_doppler_start)).total_seconds()

    def slant_range_m(self, column):
        return self.slant_range_first_m + np.asarray(column, float) * self.range_spacing_m

    def record(self):
        d = asdict(self)
        d['shape'] = list(self.shape)
        d.update(collection_duration_s=self.collection_duration_s, image_span_s=self.image_span_s)
        return d


class DwellProduct:
    """Read-only access to an ICEYE dwell SLC: metadata, geolocation, ephemeris, bounded crops."""

    def __init__(self, path):
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(self.path)
        with h5py.File(self.path, 'r') as f:
            missing = [k for k in _SCALARS + _TEXT if k not in f]
            if missing:
                raise ValueError(f'product is missing required fields: {missing}')
            if 's_i' not in f or 's_q' not in f:
                raise ValueError('product has no s_i/s_q complex arrays')
            if f['s_i'].shape != f['s_q'].shape or f['s_i'].ndim != 2:
                raise ValueError('I and Q arrays must share one two-dimensional shape')
            m = {k: float(f[k][()]) for k in _SCALARS}
            t = {k: _text(f[k][()]) for k in _TEXT}
            self.shape = tuple(int(v) for v in f['s_i'].shape)
            self._rpc = {k: np.asarray(f['RPC/' + k][()], float) for k in _RPC}
            self._drate = np.asarray(f['doppler_rate_coeffs'][()], float)
            self._sv = (np.asarray([_utc(v[0]) for v in f['state_vector_time_utc'][()]]),
                        np.stack([f['posX'][()], f['posY'][()], f['posZ'][()]], axis=-1),
                        np.stack([f['velX'][()], f['velY'][()], f['velZ'][()]], axis=-1))
            self._coa = (np.asarray(f['coa_pos'][()], float), np.asarray(f['coa_vel'][()], float))
        if int(m['number_of_azimuth_samples']) != self.shape[0] or int(m['number_of_range_samples']) != self.shape[1]:
            raise ValueError('declared sample counts disagree with the array shape')
        if abs(m['processing_prf'] * m['azimuth_time_interval'] - 1) > 1e-6:
            raise ValueError('processing PRF and azimuth time interval are inconsistent')
        if not 0 < m['total_processed_bandwidth_azimuth'] <= m['processing_prf']:
            raise ValueError('processed bandwidth does not fit the output sampling rate')
        self.acq = Acquisition(
            path=str(self.path), satellite=t['satellite_name'], product_type=t['product_type'],
            polarization=t['polarization'], look_side=t['look_side'], orbit_direction=t['orbit_direction'],
            shape=self.shape, prf_hz=m['processing_prf'], azimuth_time_interval_s=m['azimuth_time_interval'],
            bandwidth_hz=m['total_processed_bandwidth_azimuth'],
            wavelength_m=C_LIGHT / m['carrier_frequency'],
            azimuth_spacing_m=m['azimuth_ground_spacing'], range_spacing_m=m['slant_range_spacing'],
            slant_range_first_m=m['slant_range_to_first_pixel'], incidence_center_deg=m['incidence_center'],
            scene_height_m=m['avg_scene_height'], collection_start=t['acquisition_start_utc'],
            collection_end=t['acquisition_end_utc'], zero_doppler_start=t['zerodoppler_start_utc'],
            zero_doppler_end=t['zerodoppler_end_utc'])
        self._meta = m

    # ---------------------------------------------------------------- geometry
    @property
    def epoch(self):
        """Aperture-centre epoch: the midpoint of the declared collection."""
        a, b = _utc(self.acq.collection_start), _utc(self.acq.collection_end)
        return a + (b - a) / 2

    def ephemeris(self):
        times, pos, vel = self._sv
        seconds = np.array([(t - self.epoch).total_seconds() for t in times])
        return Ephemeris(seconds, pos, vel)

    def doppler_rate_hz_s(self, column):
        """Doppler rate at a range column, in the product's own polynomial convention."""
        n = self.shape[1]
        column = np.asarray(column, float)
        if np.any(column < 0) or np.any(column >= n):
            raise ValueError('range column outside the image')
        tau = (column - n / 2) / self._meta['range_sampling_rate']
        rate = np.polynomial.polynomial.polyval(tau, self._drate)
        if np.any(rate >= 0):
            raise ValueError('non-negative Doppler rate; resolve the processor convention first')
        return rate

    def zero_doppler_time_s(self, row):
        """Row zero-Doppler time in seconds from the aperture-centre epoch."""
        start = (_utc(self.acq.zero_doppler_start) - self.epoch).total_seconds()
        return start + np.asarray(row, float) * self.acq.azimuth_time_interval_s

    def aperture_time_s(self, column):
        """Bandwidth divided by Doppler rate: the real aperture duration at this column."""
        return float(self.acq.bandwidth_hz / abs(self.doppler_rate_hz_s(column)))

    def geolocate(self, latitude_deg, longitude_deg, height_m):
        from .geolocation import project_rpc
        line, sample = project_rpc(self._rpc, latitude_deg, longitude_deg, height_m)
        return float(line), float(sample)

    # ---------------------------------------------------------------- data
    def crop(self, row_centre, col_centre, n_rows, n_cols, max_bytes=512 * 1024 ** 2):
        """Bounded complex crop centred as closely as the image edges allow."""
        n_rows, n_cols = int(n_rows), int(n_cols)
        if n_rows <= 0 or n_cols <= 0:
            raise ValueError('crop sizes must be positive')
        if n_rows > self.shape[0] or n_cols > self.shape[1]:
            raise ValueError('crop is larger than the image')
        if n_rows * n_cols * 8 > max_bytes:
            raise ValueError(f'crop of {n_rows * n_cols * 8 / 1e6:.0f} MB exceeds the {max_bytes / 1e6:.0f} MB budget')
        r0 = int(np.clip(int(row_centre) - n_rows // 2, 0, self.shape[0] - n_rows))
        c0 = int(np.clip(int(col_centre) - n_cols // 2, 0, self.shape[1] - n_cols))
        with h5py.File(self.path, 'r') as f:
            i = f['s_i'][r0:r0 + n_rows, c0:c0 + n_cols]
            q = f['s_q'][r0:r0 + n_rows, c0:c0 + n_cols]
        data = (i.astype(np.float32) + 1j * q.astype(np.float32)).astype(np.complex64)
        if not np.isfinite(data).all():
            raise ValueError('crop contains nonfinite samples')
        return data, (r0, c0)


# -------------------------------------------------------------------- bank
@dataclass(frozen=True)
class DopplerBank:
    """K sub-apertures swept across the processed Doppler support.

    `width_frac` is the sub-aperture width as a fraction of the processed bandwidth.
    `shift_hz` is the reference-to-offset shift inside a pair; None selects common-reference
    mode, in which every sub-aperture is registered against the band-centre sub-aperture.
    `span_frac` shrinks the swept range, which is how the published configurations, whose
    sweep is a few hundred hertz per step, are reproduced. `start_frac` moves the start of the
    sweep along the processed band, which selects a portion of the dwell in slow time and makes
    a split-aperture consistency test possible inside one acquisition. `sweep_frac`, when given,
    fixes the swept fraction of the band directly and overrides `span_frac`, so that two halves
    of one aperture can be given identical sweeps.
    """
    K: int = 50
    width_frac: float = 0.10
    shift_hz: float | None = None
    span_frac: float = 1.0
    start_frac: float = 0.0
    sweep_frac: float | None = None

    def plan(self, acq, n_bins, doppler_rate_hz_s, centroid_hz=0.0):
        if self.K < 2 or not 0 < self.width_frac < 1 or not 0 < self.span_frac <= 1:
            raise ValueError('need K >= 2, 0 < width_frac < 1, 0 < span_frac <= 1')
        if not 0 <= self.start_frac < 1:
            raise ValueError('start_frac must lie in [0, 1)')
        if self.sweep_frac is not None and not 0 < self.sweep_frac <= 1:
            raise ValueError('sweep_frac must lie in (0, 1]')
        if n_bins < 16:
            raise ValueError('azimuth crop is too short to resolve the band')
        rate = float(doppler_rate_hz_s)
        if rate >= 0:
            raise ValueError('Doppler rate must be negative in this convention')
        bin_hz = acq.prf_hz / n_bins
        band = acq.bandwidth_hz
        width = max(round(self.width_frac * band / bin_hz), 2) * bin_hz
        shift = 0.0 if self.shift_hz is None else max(round(self.shift_hz / bin_hz), 1) * bin_hz
        start = self.start_frac * band
        free = band - width - shift - start
        if free <= 0:
            raise ValueError('sub-aperture width, offset and start leave no room inside the processed band')
        sweep = self.span_frac * free if self.sweep_frac is None else self.sweep_frac * band
        if sweep > free + 1e-9:
            raise ValueError('requested sweep does not fit inside the processed band')
        step = np.floor(sweep / (self.K - 1) / bin_hz) * bin_hz
        if step < bin_hz:
            raise ValueError('crop is too short for distinct sub-apertures; widen the crop or reduce K')
        lo = centroid_hz - band / 2 + start + np.arange(self.K) * step
        centre = lo + width / 2
        # Slow time relative to the row's zero-Doppler crossing, from f = f_dc + Ka (t - t_zd).
        times = (centre - centroid_hz) / rate
        return {'bin_hz': bin_hz, 'width_hz': width, 'shift_hz': shift, 'step_hz': step,
                'start_hz': start,
                'ref_lo_hz': lo, 'off_lo_hz': lo + shift, 'centre_hz': centre,
                'times_s': times, 'sub_aperture_s': width / abs(rate),
                'span_s': float(times.max() - times.min()),
                'sampling_hz': float((self.K - 1) / (times.max() - times.min())) if self.K > 1 else 0.0,
                'pair_lag_s': shift / abs(rate),
                'overlap_within_pair': float(max(0.0, width - shift) / width)}

    def masks(self, plan, freqs):
        f = np.asarray(freqs, float)[None, :]
        lo_r = plan['ref_lo_hz'][:, None]
        lo_o = plan['off_lo_hz'][:, None]
        w = plan['width_hz']
        ref = ((f >= lo_r) & (f < lo_r + w)).astype(np.float32)
        off = ((f >= lo_o) & (f < lo_o + w)).astype(np.float32)
        if not ref.any(axis=1).all() or not off.any(axis=1).all():
            raise ValueError('a sub-aperture is empty on this frequency grid')
        return ref, off


def occupied_band(crop, prf_hz, threshold_db=-12.0):
    """Measured Doppler support and centroid of a crop, as a check on the declared bandwidth."""
    spectrum = np.abs(np.fft.fft(crop, axis=0)) ** 2
    power = np.fft.fftshift(spectrum.mean(axis=1))
    freqs = np.fft.fftshift(np.fft.fftfreq(crop.shape[0], d=1.0 / prf_hz))
    keep = power > power.max() * 10 ** (threshold_db / 10)
    if not keep.any():
        raise ValueError('no azimuth spectral support above threshold')
    inside = freqs[keep]
    return {'low_hz': float(inside.min()), 'high_hz': float(inside.max()),
            'width_hz': float(inside.max() - inside.min()),
            'centroid_hz': float((freqs * power).sum() / power.sum()),
            'freqs_hz': freqs, 'power': power / power.max()}


def azimuth_shift_to_velocity(shift_px, acq, slant_range_m, platform_speed_m_s):
    """Convert an azimuth registration shift to line-of-sight velocity, v = dx V / R."""
    if slant_range_m <= 0 or platform_speed_m_s <= 0:
        raise ValueError('slant range and platform speed must be positive')
    return np.asarray(shift_px, float) * acq.azimuth_spacing_m * platform_speed_m_s / slant_range_m
