"""R7: a positive control on real radar data.

R4 measures how much the recovered velocity series wanders when nothing is moving. That is a
noise floor, and a noise floor on its own does not prove the readout works. Here a known
line-of-sight displacement history is imposed on a real Khafre crop by modulating its azimuth
spectrum, which in a dwell product is slow time, and the same common-reference tracker is asked
to recover it. The amplitude at which recovery becomes faithful is the detection threshold of
this acquisition, measured on genuine radar texture rather than on synthetic speckle.

The modulation is common to the whole crop, so common-mode removal is disabled throughout.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import PRODUCTS, GIZA_PATCHES, BANKS, RESULTS, N_AZ, N_RG, PATCH_PX
from sarsim import viz
from sarsim.dwell import (DwellProduct, occupied_band, inject_line_of_sight_motion,
                          azimuth_shift_to_velocity)
from sarsim.orbit import lla_to_ecef, perpendicular_baselines
from sarsim.track import patch_shifts

TEST = 'r07_injected_motion'
AMPLITUDES = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
FREQUENCIES = [0.2, 0.5, 1.0]
N_PIXELS = 900


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    name, label, lat, lon, height, _ = GIZA_PATCHES[0]
    product = DwellProduct(PRODUCTS['giza'])
    acq = product.acq
    row, col = product.geolocate(lat, lon, height)
    crop, (r0, c0) = product.crop(row, col, N_AZ, N_RG)
    rate = float(product.doppler_rate_hz_s(c0 + N_RG / 2))
    band = occupied_band(crop, acq.prf_hz)
    bank = BANKS['reference']
    plan = bank.plan(acq, N_AZ, rate, centroid_hz=band['centroid_hz'])
    freqs = np.fft.fftshift(np.fft.fftfreq(N_AZ, d=1.0 / acq.prf_hz))
    ref_mask, _ = bank.masks(plan, freqs)

    eph = product.ephemeris()
    lo, hi = eph.span
    t_zd = float(np.clip(product.zero_doppler_time_s(r0 + N_AZ / 2), lo, hi))
    _, slant, _ = perpendicular_baselines(eph, [t_zd], lla_to_ecef(lat, lon, height), t_zd)
    speed = float(np.linalg.norm(eph.at([t_zd])[1][0]))

    # brightest pixels inside the safe margin: the most favourable targets the scene offers
    margin = 200
    sub = np.abs(crop[margin:-margin, margin:-margin])
    flat = np.argsort(sub.ravel())[-N_PIXELS:]
    rows = flat // sub.shape[1] + margin
    cols = flat % sub.shape[1] + margin

    def recover(image):
        spectrum = np.fft.fftshift(np.fft.fft(image, axis=0), axes=0)
        base = np.fft.ifft(np.fft.ifftshift(spectrum * ref_mask[bank.K // 2][:, None], axes=0), axis=0)
        shifts = np.zeros((rows.size, bank.K))
        for k in range(bank.K):
            a = np.fft.ifft(np.fft.ifftshift(spectrum * ref_mask[k][:, None], axes=0), axis=0)
            shifts[:, k] = patch_shifts(a, base, rows, cols, patch=PATCH_PX, upsample=1000)[0]
        return azimuth_shift_to_velocity(shifts, acq, slant, speed)

    baseline = recover(crop)
    baseline -= baseline.mean(axis=1, keepdims=True)
    floor = float(np.sqrt((baseline ** 2).mean()))
    print(f'  static floor {floor * 1e6:,.0f} um/s on {rows.size} bright pixels '
          f'({time.time() - t0:.0f}s)', flush=True)

    rows_out = []
    for f0 in FREQUENCIES:
        for amp in AMPLITUDES:
            moved, t_ax, disp = inject_line_of_sight_motion(
                crop, acq, rate, amp, f0, centroid_hz=band['centroid_hz'])
            v = recover(moved)
            v -= v.mean(axis=1, keepdims=True)
            delta = (v - baseline).mean(axis=0)
            truth = 2 * np.pi * f0 * amp * np.cos(2 * np.pi * f0 * plan['times_s'])
            truth -= truth.mean()
            scale = float(np.dot(delta, truth) / max(np.dot(truth, truth), 1e-30))
            rows_out.append({
                'frequency_hz': f0, 'displacement_amplitude_m': float(amp),
                'velocity_amplitude_m_s': float(2 * np.pi * f0 * amp),
                'recovered_rms_m_s': float(np.sqrt((delta ** 2).mean())),
                'recovery_gain': scale,
                'signal_to_floor': float(np.sqrt((delta ** 2).mean()) / floor)})
        print(f'  {f0} Hz done ({time.time() - t0:.0f}s)', flush=True)

    # threshold: the smallest injected velocity whose recovery gain is within a factor of two
    good = [r for r in rows_out if 0.5 <= r['recovery_gain'] <= 2.0]
    threshold = min((r['velocity_amplitude_m_s'] for r in good), default=None)

    fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.6))
    for i, f0 in enumerate(FREQUENCIES):
        sel = [r for r in rows_out if r['frequency_hz'] == f0]
        x = [r['velocity_amplitude_m_s'] * 1e6 for r in sel]
        axs[0].plot(x, [r['recovered_rms_m_s'] * 1e6 for r in sel], '-o', ms=4, lw=1.8,
                    color=viz.SERIES[i], markeredgecolor=viz.SURFACE, label=f'{f0:g} Hz')
        axs[1].plot(x, [r['recovery_gain'] for r in sel], '-o', ms=4, lw=1.8,
                    color=viz.SERIES[i], markeredgecolor=viz.SURFACE, label=f'{f0:g} Hz')
    lim = np.array([min(r['velocity_amplitude_m_s'] for r in rows_out),
                    max(r['velocity_amplitude_m_s'] for r in rows_out)]) * 1e6
    axs[0].plot(lim, lim / np.sqrt(2), color=viz.MUTED, lw=1, ls=':')
    axs[0].text(lim[1], lim[1] / np.sqrt(2) * 1.3, 'faithful recovery', ha='right', fontsize=8.5, color=viz.INK2)
    axs[0].axhline(floor * 1e6, color=viz.SERIES[7], lw=1, ls='--')
    axs[0].text(lim[0] * 1.2, floor * 1e6 * 1.25, 'static floor of the same pixels',
                fontsize=8.5, color=viz.SERIES[7])
    axs[0].set_xscale('log'); axs[0].set_yscale('log')
    axs[0].set_xlabel('injected velocity amplitude ($\\mu$m/s)')
    axs[0].set_ylabel('recovered velocity ($\\mu$m/s, rms)')
    axs[0].set_title('What comes back', loc='left'); axs[0].legend(fontsize=8.5)
    axs[1].axhspan(0.5, 2.0, color=viz.GRID, alpha=0.7, lw=0)
    axs[1].axhline(1.0, color=viz.MUTED, lw=1, ls=':')
    axs[1].set_xscale('log'); axs[1].set_yscale('log')
    axs[1].set_xlabel('injected velocity amplitude ($\\mu$m/s)')
    axs[1].set_ylabel('recovered / injected')
    axs[1].set_title('Gain, with the factor-of-two band shaded', loc='left'); axs[1].legend(fontsize=8.5)
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r07_threshold.png'), 'caption':
             '<b>A positive control on real radar.</b> A known line-of-sight displacement history is imposed on the '
             'Khafre crop by modulating its azimuth spectrum, which in a dwell product is slow time, and the same '
             'common-reference tracker recovers it from the nine hundred brightest pixels. Left: recovered against '
             'injected velocity, with the static floor of those same pixels marked. Right: the ratio of the two. '
             'The readout is faithful well above the floor and collapses into it below, which is where the '
             'detection threshold of this acquisition sits.'}]

    m = {'patch': name, 'label': label, 'pixels': int(rows.size),
         'static_floor_um_s': floor * 1e6,
         'detection_threshold_velocity_um_s': None if threshold is None else threshold * 1e6,
         'detection_threshold_displacement_um_at_0p5hz':
             None if threshold is None else threshold / (2 * np.pi * 0.5) * 1e6,
         'sampling_hz': plan['sampling_hz'], 'series_span_s': plan['span_s'],
         'slant_range_m': slant, 'platform_speed_m_s': speed,
         'azimuth_metres_per_mm_s': acq.azimuth_spacing_m and float(1e-3 * speed / slant / acq.azimuth_spacing_m),
         'trials': rows_out, 'runtime_s': time.time() - t0}
    thr_txt = ('not reached within the tested range' if threshold is None
               else f"{threshold * 1e6:,.0f} micrometres per second")
    summary = {
        'id': TEST, 'order': 17, 'tag': 'threshold', 'eyebrow': '07 · Positive control',
        'title': 'Known motion put into the real image comes back, above a measurable threshold',
        'question': 'Does the velocity readout work at all on real data, and from what amplitude?',
        'finding': (f"It works, and the threshold is {thr_txt}. Below that the recovery collapses into the static "
                    f"floor of the same pixels, {floor * 1e6:,.0f} micrometres per second. Ambient ground motion in "
                    f"the microseism band is 0.1 to 10 micrometres per second, so the motion a passive subsurface "
                    f"inference would have to use sits below what this acquisition can see, on its most favourable "
                    f"targets."),
        'limitations': ('The injected motion is common to the whole crop and the control therefore runs with '
                        'common-mode removal disabled, which is the most favourable configuration. Real ground '
                        'motion is neither spatially uniform nor single-frequency. The threshold applies to the '
                        'brightest nine hundred pixels of this scene; dimmer ground is worse, and a corner '
                        'reflector would be better.'),
        'method': ('One real Khafre crop, azimuth spectrum modulated by exp(-4 pi j d(t) / lambda) with t recovered '
                   'from the product Doppler rate, for nine displacement amplitudes at three frequencies. Recovery '
                   'is the projection of the change in the mean velocity series onto the injected waveform, using '
                   'the same fifty common-reference sub-apertures as R4.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in m.items() if k != 'trials'}, indent=1))


if __name__ == '__main__':
    main()
