"""P2-01 · The real dwell as the simulator's geometry.

    uv run python experiments/p2_01_geometry.py

Reads the two ICEYE Spotlight Dwell Fine products' metadata (no imagery) and writes their
acquisition records to sites/acquisitions/, from which `sarsim.acquisition.DwellGeometry` builds
the simulator's geometry anywhere, without the products. Where the products are absent (the
cloud), the committed records are used and only the synthetic checks run again.

Checks, each against the product's own numbers:

1. Three speeds, not one: the platform's V_s (state vectors), the ground sweep of the image rows
   V_g (row spacing over row time), and V_eff from the Doppler rate. V_eff^2 = V_s V_g, and the
   flat-earth model's single speed misplaces the Doppler rate by about nine per cent either way.
2. The aperture: processed bandwidth over Doppler rate; the virtual baseline from the state
   vectors, fitted by a cubic; the span of the first investigation's sub-aperture centres (the
   "168 km" it quoted) against the whole processed aperture.
3. Looks: resolution times duration is fixed by the orbit, V_g / |Ka| (x 0.886).
4. The synthesizer on this geometry: a point target's impulse response has the predicted widths,
   and a target receding at constant speed moves in azimuth by 2 V_g v / (lambda Ka), which is
   R v / V_s: the conversion every velocity series in Phase 2 uses.
"""
import json
from pathlib import Path

import numpy as np

from katabasis.runs import Run
from sarsim import point_targets, synthesize
from sarsim.acquisition import DwellGeometry, RECORDS

PRODUCTS = {
    'giza-20250827': Path.home() / 'tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5',
    'sacsayhuaman-20250822': Path.home() / 'tmp/sar/sacsayhuaman/ICEYE_X35_SLC_SLEDF_5907293_20250822T152433.h5',
}
REFERENCE_BANK = {'K': 50, 'width_frac': 0.10}     # the first investigation's reference bank (legacy real_common)


def geometry(name):
    path = PRODUCTS[name]
    if path.is_file():
        g = DwellGeometry.from_product(path, name=name)
        RECORDS.mkdir(parents=True, exist_ok=True)
        rec = g.record()
        rec['about'] = ('Acquisition record for the Phase 2 simulator, read from the product metadata by '
                        'sim/experiments/p2_01_geometry.py. Imagery is not included.')
        (RECORDS / f'{name}.json').write_text(json.dumps(rec, indent=1))
        return g, 'product'
    return DwellGeometry.from_record(name), 'record'


def widths(img, axis, up=16):
    """-3 dB width (pixels) of the brightest point's response along an axis: the complex line is interpolated
    `up` times by zero-padding its spectrum (the slant axis holds only 1.4 samples per resolution cell), then
    the half-power crossings are found by linear interpolation."""
    i, j = np.unravel_index(np.argmax(np.abs(img)), img.shape)
    z = img[:, j] if axis == 0 else img[i, :]
    n = len(z)
    F = np.fft.fftshift(np.fft.fft(z))
    pad = np.zeros(n * up, complex)
    pad[(n * up - n) // 2:(n * up + n) // 2] = F
    line = np.abs(np.fft.ifft(np.fft.ifftshift(pad))) ** 2
    c = int(np.argmax(line))
    half = line[c] / 2
    lo = c
    while lo > 0 and line[lo] > half:
        lo -= 1
    hi = c
    while hi < len(line) - 1 and line[hi] > half:
        hi += 1
    fl = lo + (half - line[lo]) / (line[lo + 1] - line[lo])
    fh = hi - 1 + (line[hi - 1] - half) / (line[hi - 1] - line[hi])
    return (fh - fl) / up


def peak_position(img, up=64):
    """Sub-pixel azimuth position of the brightest response (row index), by zero-padded FFT interpolation."""
    j = np.unravel_index(np.argmax(np.abs(img)), img.shape)[1]
    col = img[:, j]
    n = len(col)
    F = np.fft.fftshift(np.fft.fft(col))
    pad = np.zeros(n * up, complex)
    pad[(n * up - n) // 2:(n * up + n) // 2] = F
    fine = np.abs(np.fft.ifft(np.fft.ifftshift(pad)))
    return np.argmax(fine) / up


def synthetic_checks(g):
    shape = (512, 128)
    still = synthesize(point_targets(np.array([0.0]), np.array([0.0]), np.array([0.0]), np.ones(1)), g, shape)
    w_az = widths(still, 0) * g.dx
    w_rg = widths(still, 1) * g.dr
    v = 2e-3                                   # m/s, receding
    moving = synthesize(point_targets(np.array([0.0]), np.array([0.0]), np.array([0.0]), np.ones(1)), g, shape,
                        motion=lambda x, y, z, t: v * t[:, None] * np.ones(len(x)))
    shift = (peak_position(moving) - peak_position(still)) * g.dx
    predicted = 2 * g.V * v / (g.lam * g.Ka_signed)
    return {'impulse_width_azimuth_m': w_az, 'predicted_azimuth_m': g.resolution,
            'impulse_width_slant_m': w_rg, 'predicted_slant_m': 0.886 / g.kr_band,
            'receding_speed_m_s': v, 'azimuth_shift_m': shift, 'predicted_shift_m': predicted,
            'recovered_speed_m_s': float(g.shift_to_velocity(shift))}


def main():
    params = {'products': {k: p.name for k, p in PRODUCTS.items()}, 'reference_bank': REFERENCE_BANK,
              'look_durations_s': [0.25, 0.5, 1.0, 1.96, 4.0]}
    with Run('p2_01_geometry', 'The real dwell as the simulator\'s geometry', params) as run:
        out = {}
        for name in PRODUCTS:
            g, origin = geometry(name)
            ka_flat_g = 2 * g.V ** 2 / (g.lam * g.R0)
            ka_flat_s = 2 * g.V_platform ** 2 / (g.lam * g.R0)
            T = g.aperture_time
            w = REFERENCE_BANK['width_frac']
            centres = np.linspace(-(1 - w) / 2, (1 - w) / 2, REFERENCE_BANK['K']) * g.nu_band
            b_centres = g.nu_to_baseline(centres)
            s = g.source
            out[name] = {
                'origin': origin, 'satellite': s['satellite'], 'mode': s['mode'],
                'collection_start': s['collection_start'],
                'speeds_m_s': {'platform': g.V_platform, 'ground_sweep': g.V, 'effective': g.V_eff},
                'veff2_over_vs_vg': g.V_eff ** 2 / (g.V_platform * g.V),
                'doppler_rate_hz_s': g.Ka_signed,
                'flat_earth_rate_error': {'with_ground_speed': ka_flat_g / g.Ka - 1,
                                          'with_platform_speed': ka_flat_s / g.Ka - 1},
                'wavelength_m': g.lam, 'slant_range_m': g.R0, 'incidence_deg': g.theta_deg,
                'incidence_from_los_deg': float(np.rad2deg(np.arccos(g.los_enu[2]))),
                'los_enu': list(g.los_enu), 'heading_deg': g.heading_deg,
                'aperture_s': T, 'collection_s': s['collection_duration_s'], 'image_span_s': s['image_span_s'],
                'doppler_band_hz': g.doppler_band_hz, 'aspect_span_deg': g.aspect_span_deg,
                'baseline_span_m': s['baseline_span_m'], 'baseline_fit_residual_m': s['baseline_fit_residual_m'],
                'reference_bank_centre_span_m': float(b_centres.max() - b_centres.min()),
                'azimuth_resolution_m': g.resolution, 'slant_resolution_m': 0.886 / g.kr_band,
                'pixel_m': {'azimuth': g.dx, 'slant': g.dr},
                'windows': {'range': s['window_range'], 'azimuth': s['window_azimuth']},
                'time_space_m_s': g.time_space_m_s,
                'looks': [{'duration_s': d, 'resolution_m': float(g.look_resolution(d))}
                          for d in params['look_durations_s']],
                'velocity_per_azimuth_metre': float(abs(g.shift_to_velocity(1.0))),
                'platform_speed_over_range': g.V_platform / g.R0,
                'synthetic': synthetic_checks(g),
            }
        gz = out['giza-20250827']
        finding = (f"The real dwell has three speeds, not one: the platform's {gz['speeds_m_s']['platform']:,.0f} m/s, "
                   f"the image rows' sweep over the ground at {gz['speeds_m_s']['ground_sweep']:,.0f} m/s, and "
                   f"{gz['speeds_m_s']['effective']:,.0f} m/s from the Doppler rate, their geometric mean to "
                   f"{abs(gz['veff2_over_vs_vg'] - 1) * 1e6:.0f} parts per million. The first investigation's "
                   f"flat-earth simulator used one and put the Doppler rate {abs(gz['flat_earth_rate_error']['with_ground_speed']) * 100:.0f}% off. "
                   f"Over the {gz['aperture_s']:.2f} s processed aperture the platform moves "
                   f"{gz['baseline_span_m'] / 1e3:.1f} km ({gz['reference_bank_centre_span_m'] / 1e3:.0f} km between the "
                   f"reference bank's first and last sub-aperture centres, the figure quoted before). A look's azimuth "
                   f"resolution times its duration is fixed at {0.886 * gz['time_space_m_s']:.2f} m s: a 1.96 s look "
                   f"resolves {gz['looks'][3]['resolution_m']:.2f} m. On this geometry the synthesizer's point response "
                   f"is {gz['synthetic']['impulse_width_azimuth_m'] * 100:.1f} cm by {gz['synthetic']['impulse_width_slant_m'] * 100:.0f} cm, "
                   f"and a target receding at 2 mm/s moves {gz['synthetic']['azimuth_shift_m']:.2f} m in azimuth "
                   f"(predicted {gz['synthetic']['predicted_shift_m']:.2f} m).")
        run.save({'sites': out, 'finding': finding})


if __name__ == '__main__':
    main()
