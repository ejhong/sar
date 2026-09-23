"""Shared orchestration for the real-data (dwell) chapter.

One patch of a real ICEYE dwell SLC goes in; sub-aperture trajectories, a depth volume and
the summary statistics used for the within-image comparison come out. Everything that depends
on an acquisition fact is read from the product; nothing is defaulted.
"""
from __future__ import annotations

import json
import os
import sys
import time
import hashlib

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from sarsim.dwell import DwellProduct, DopplerBank, occupied_band, azimuth_shift_to_velocity  # noqa: E402
from sarsim.orbit import lla_to_ecef, perpendicular_baselines, steering_wavenumbers  # noqa: E402
from sarsim.track import patch_shifts  # noqa: E402
from sarsim.tomo import focus_paper, focus_windows, nyquist_depth  # noqa: E402

RESULTS = os.path.join(ROOT, 'results')
CACHE = os.path.join(RESULTS, 'cache')
PRODUCTS = {
    'giza': os.path.expanduser('~/tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5'),
    'sacsayhuaman': os.path.expanduser('~/tmp/sar/sacsayhuaman/ICEYE_X35_SLC_SLEDF_5907293_20250822T152433.h5'),
}

# Ellipsoidal height = published elevation + local geoid undulation (about +15.5 m at Giza,
# about -41 m in the Cusco region). Height only has to be good enough to land the RPC on the
# right structure; every patch is verified against the image before it is used.
GIZA_PATCHES = [
    ('khafre',        'Khafre pyramid',      29.97603, 31.13078, 158 + 15.5, 'monument'),
    ('khufu',         'Khufu pyramid',       29.97917, 31.13417, 146 + 15.5, 'monument'),
    ('menkaure',      'Menkaure pyramid',    29.97250, 31.12833, 125 + 15.5, 'monument'),
    ('sphinx',        'Sphinx and temples',  29.97528, 31.13750,  35 + 15.5, 'monument'),
    ('desert_west',   'Open plateau, west',  29.97600, 31.11800,  80 + 15.5, 'control'),
    ('desert_south',  'Open plateau, south', 29.96700, 31.12600,  75 + 15.5, 'control'),
]

# Sacsayhuaman sits on steep Andean terrain at about 3.6 km, where a single assumed height
# places a patch far less precisely than at flat Giza. Every patch here was checked against the
# image before use: the zigzag terrace walls, the Rodadero outcrop, the Cusco street grid and an
# open hillside are all recognisable at their projected positions.
SACSAYHUAMAN_PATCHES = [
    ('walls',       'Sacsayhuaman terrace walls', -13.5095, -71.9820, 3660, 'monument'),
    ('rodadero',    'Rodadero outcrop',           -13.5085, -71.9845, 3640, 'monument'),
    ('cusco',       'Cusco street grid',          -13.5170, -71.9785, 3360, 'urban'),
    ('slope_north', 'Open hillside, north',       -13.5040, -71.9820, 3680, 'control'),
]
SITES = {'giza': GIZA_PATCHES, 'sacsayhuaman': SACSAYHUAMAN_PATCHES}

N_AZ, N_RG = 4096, 1536
MARGIN = 160
STRIDE_AZ, STRIDE_RG = 12, 8
LAM_S = 0.48
PATCH_PX = 32

BANKS = {
    'reference': DopplerBank(K=50, width_frac=0.10, shift_hz=None),
    'paper': DopplerBank(K=50, width_frac=0.50, shift_hz=88.0),
    'v17': DopplerBank(K=20, width_frac=0.50, shift_hz=404.0, span_frac=0.062),
    # Two disjoint halves of the same dwell, given identical sweeps so their depth axes match.
    'split_a': DopplerBank(K=25, width_frac=0.10, start_frac=0.00, sweep_frac=0.45),
    'split_b': DopplerBank(K=25, width_frac=0.10, start_frac=0.45, sweep_frac=0.45),
}


def digest(*parts):
    return hashlib.sha256('|'.join(str(p) for p in parts).encode()).hexdigest()[:12]


def target_grid(shape, margin=MARGIN, stride_az=STRIDE_AZ, stride_rg=STRIDE_RG):
    rows = np.arange(margin, shape[0] - margin, stride_az)
    cols = np.arange(margin, shape[1] - margin, stride_rg)
    RR, CC = np.meshgrid(rows, cols, indexing='ij')
    return RR.ravel(), CC.ravel(), (rows, cols)


def run_patch(site, key, bank_name='paper', lam_s=LAM_S, n_az=N_AZ, n_rg=N_RG,
              stride_az=STRIDE_AZ, stride_rg=STRIDE_RG, force=False, verbose=True):
    """Process one named patch and cache the trajectories, geometry and depth volume."""
    if site not in SITES:
        raise KeyError(f'unknown site {site}')
    entry = next((p for p in SITES[site] if p[0] == key), None)
    if entry is None:
        raise KeyError(f'unknown patch {site}/{key}')
    name, label, lat, lon, height, kind = entry
    bank = BANKS[bank_name]
    tag = digest(site, key, bank_name, lam_s, n_az, n_rg, stride_az, stride_rg, PATCH_PX)
    path = os.path.join(CACHE, f'real_{site}_{key}_{bank_name}_{tag}.npz')
    if os.path.exists(path) and not force:
        d = np.load(path, allow_pickle=True)
        out = {k: d[k] for k in d.files}
        out['meta'] = json.loads(str(out['meta']))
        return out

    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    product = DwellProduct(PRODUCTS[site])
    acq = product.acq
    row, col = product.geolocate(lat, lon, height)
    crop, (r0, c0) = product.crop(row, col, n_az, n_rg)
    rate = float(product.doppler_rate_hz_s(c0 + n_rg / 2))
    t_zd = float(product.zero_doppler_time_s(r0 + n_az / 2))
    band = occupied_band(crop, acq.prf_hz)
    plan = bank.plan(acq, n_az, rate, centroid_hz=band['centroid_hz'])

    eph = product.ephemeris()
    lo, hi = eph.span
    times = np.clip(t_zd + plan['times_s'], lo, hi)
    target = lla_to_ecef(lat, lon, height)
    bperp, slant, incidence = perpendicular_baselines(eph, times, target, float(np.clip(t_zd, lo, hi)))
    kz = steering_wavenumbers(bperp, slant, incidence, lam_s)
    order = np.argsort(kz)
    z_nyq = nyquist_depth(kz[order])
    z = np.linspace(z_nyq / 160, 0.97 * z_nyq, 160)

    rows, cols, axes = target_grid(crop.shape, stride_az=stride_az, stride_rg=stride_rg)
    freqs = np.fft.fftshift(np.fft.fftfreq(n_az, d=1.0 / acq.prf_hz))
    ref_mask, off_mask = bank.masks(plan, freqs)
    spectrum = np.fft.fftshift(np.fft.fft(crop, axis=0), axes=0)

    K = bank.K
    Y = np.zeros((rows.size, K, 2), np.float32)
    if bank.shift_hz is None:                       # common-reference mode
        centre = K // 2
        base = np.fft.ifft(np.fft.ifftshift(spectrum * ref_mask[centre][:, None], axes=0), axis=0)
    for k in range(K):
        a = np.fft.ifft(np.fft.ifftshift(spectrum * ref_mask[k][:, None], axes=0), axis=0)
        b = base if bank.shift_hz is None else \
            np.fft.ifft(np.fft.ifftshift(spectrum * off_mask[k][:, None], axes=0), axis=0)
        dr, dc = patch_shifts(a, b, rows, cols, patch=PATCH_PX, upsample=1000)
        Y[:, k, 0], Y[:, k, 1] = dr, dc
        if verbose and (k % 10 == 0 or k == K - 1):
            print(f'    {site}/{key}/{bank_name} sub-aperture {k + 1}/{K}  {time.time() - t0:.0f}s', flush=True)

    q = Y - np.median(Y, axis=0, keepdims=True)
    tomo = focus_paper(q[..., 0] + 1j * q[..., 1], kz[order], z).astype(np.float32)
    score, best = focus_windows(q[:, order, :].astype(float), kz[order], z, W=min(25, K - 1))
    amplitude = np.abs(crop[rows, cols])
    speed = float(np.linalg.norm(eph.at([float(np.clip(t_zd, lo, hi))])[1][0]))
    # An azimuth shift maps to line-of-sight velocity only when the two bands share a centre
    # time. In common-reference mode that holds and the series is a velocity; in pair mode the
    # bands are separated by the R/O lag, so the same conversion yields a velocity increment
    # across that lag. The name records which one this is.
    velocity = azimuth_shift_to_velocity(Y[..., 0], acq, slant, speed)
    velocity_kind = 'line_of_sight_velocity' if bank.shift_hz is None else 'velocity_increment_over_pair_lag'

    meta = {
        'site': site, 'patch': name, 'label': label, 'kind': kind,
        'latitude': lat, 'longitude': lon, 'height_m': height,
        'rpc_row': row, 'rpc_col': col, 'crop_origin': [r0, c0], 'crop_shape': list(crop.shape),
        'bank': bank_name, 'K': K, 'width_frac': bank.width_frac, 'shift_hz': bank.shift_hz,
        'span_frac': bank.span_frac, 'start_frac': bank.start_frac, 'sweep_frac': bank.sweep_frac,
        'sub_aperture_times_s': plan['times_s'].tolist(), 'lam_s_m': lam_s, 'patch_px': PATCH_PX,
        'stride': [stride_az, stride_rg], 'pixels': int(rows.size),
        'doppler_rate_hz_s': rate, 'zero_doppler_time_s': t_zd,
        'measured_band_hz': band['width_hz'], 'measured_centroid_hz': band['centroid_hz'],
        'declared_band_hz': acq.bandwidth_hz,
        'sub_aperture_s': plan['sub_aperture_s'], 'centres_span_s': plan['span_s'],
        'sampling_hz': plan['sampling_hz'], 'pair_lag_s': plan['pair_lag_s'],
        'overlap_within_pair': plan['overlap_within_pair'],
        'baseline_span_m': float(bperp.max() - bperp.min()),
        'slant_range_m': slant, 'incidence_deg': incidence, 'platform_speed_m_s': speed,
        'nyquist_depth_m': float(z_nyq), 'runtime_s': time.time() - t0,
        'velocity_kind': velocity_kind,
        'acquisition': acq.record(),
    }
    out = {'Y': Y, 'q': q.astype(np.float32), 'kz': kz[order], 'z': z, 'tomo': tomo,
           'score': score.astype(np.float32), 'best': best, 'rows': rows, 'cols': cols,
           'grid_rows': axes[0], 'grid_cols': axes[1], 'amplitude': amplitude.astype(np.float32),
           'velocity': velocity.astype(np.float32), 'bperp': bperp,
           'band_freqs': band['freqs_hz'].astype(np.float32), 'band_power': band['power'].astype(np.float32),
           'meta': json.dumps(meta)}
    np.savez_compressed(path, **out)
    out['meta'] = meta
    if verbose:
        print(f'  {site}/{key}/{bank_name}: {rows.size:,} px in {time.time() - t0:.0f}s -> {path}', flush=True)
    return out


def patch_stats(out):
    """Summary statistics used for the within-image comparison."""
    m = out['meta']
    score = out['score']
    smax = score.max(axis=1)
    zbest = out['z'][np.argmax(score, axis=1)]
    tomo = out['tomo']
    tn = tomo / np.maximum(tomo.sum(axis=1, keepdims=True), 1e-30)
    amp = out['amplitude']
    return {
        'patch': m['patch'], 'label': m['label'], 'kind': m['kind'], 'pixels': m['pixels'],
        'amplitude_median': float(np.median(amp)), 'amplitude_p99': float(np.percentile(amp, 99)),
        'trajectory_rms_px': float(out['q'][..., 0].std()),
        'velocity_kind': m['velocity_kind'],
        'velocity_rms_mm_s': float(out['velocity'].std() * 1e3),
        'best_score_mean': float(smax.mean()), 'best_score_p95': float(np.percentile(smax, 95)),
        'frac_score_gt_0p5': float((smax > 0.5).mean()),
        'depth_median_m': float(np.median(zbest)),
        'spectrum_peak_to_median': float(tn.mean(axis=0).max() / np.median(tn.mean(axis=0))),
        'nyquist_depth_m': m['nyquist_depth_m'],
    }
