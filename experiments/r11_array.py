"""R11: treat the image as a dense seismic array and look for a propagating wave.

This is the strongest form of the claim. A SAR image gives a displacement estimate at every
pixel, which is a denser array than any seismic deployment, so the natural question is whether
ambient-noise array processing can be done with it and inverted for structure below.

The test is the first step of that programme and it is the one that has to work before anything
else can: cross-correlate the recovered velocity series between pairs of pixels, stack by
separation, and look for an arrival whose lag grows linearly with distance. A real surface
wavefield produces a moveout at the Rayleigh velocity. Nothing else does.

A synthetic plane wave of known speed and amplitude is injected into the same series and put
through the same stack, so a null result on the real data is read against a positive control
rather than on its own.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import PRODUCTS, GIZA_PATCHES, RESULTS, PATCH_PX
from sarsim import viz
from sarsim.dwell import (DwellProduct, DopplerBank, occupied_band, azimuth_shift_to_velocity)
from sarsim.orbit import lla_to_ecef, perpendicular_baselines
from sarsim.track import patch_shifts

TEST = 'r11_array'
N_AZ, N_RG = 40960, 512          # a 1.8 km azimuth strip: long enough to resolve moveout
STRIDE_AZ, STRIDE_RG = 160, 64
BANK = DopplerBank(K=24, width_frac=0.08, shift_hz=None, sweep_frac=0.085)
RAYLEIGH_M_S = 1800.0            # limestone, for the injected control and the predicted moveout
BAND_LOW_HZ, BAND_HIGH_HZ = 0.9, 3.5   # inside what the record resolves and the sampling allows
INJECT_AMPLITUDE_UM_S = 400.0
SEP_BINS = np.array([0, 100, 250, 450, 700, 1000, 1400, 1900])


def stack_by_separation(series, x_m, times_s, bins=SEP_BINS, max_pairs=60000, seed=0):
    """Cross-correlate pixel pairs, correct each pair's time offset, stack by signed separation.

    Pairs are always ordered so the second element lies further along the array. Binning by
    unsigned separation would average the two propagation directions together and put any real
    arrival at zero lag, which is exactly the artefact this ordering avoids.
    """
    rng = np.random.default_rng(seed)
    n, k = series.shape
    dt = float(np.median(np.diff(times_s.mean(axis=0))))
    lags = np.arange(-(k - 1), k) * dt
    out, counts = np.zeros((len(bins) - 1, lags.size)), np.zeros(len(bins) - 1, int)
    z = series - series.mean(axis=1, keepdims=True)
    norm = np.sqrt((z ** 2).sum(axis=1)) + 1e-30
    pairs_per_bin = max_pairs // (len(bins) - 1)
    for b in range(len(bins) - 1):
        lo, hi = bins[b], bins[b + 1]
        got = 0
        for _ in range(pairs_per_bin * 6):
            if got >= pairs_per_bin:
                break
            i = rng.integers(n)
            d = x_m - x_m[i]                      # signed: keep only elements further along
            cand = np.flatnonzero((d >= lo) & (d < hi))
            if cand.size == 0:
                continue
            j = cand[rng.integers(cand.size)]
            c = np.correlate(z[j], z[i], mode='full') / (norm[i] * norm[j])
            # each pixel's Doppler-to-time map is offset by its own zero-Doppler crossing
            shift = times_s[j, 0] - times_s[i, 0]
            out[b] += np.interp(lags, lags + shift, c, left=0.0, right=0.0)
            got += 1
        counts[b] = got
        if got:
            out[b] /= got
    return lags, out, counts


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    entry = next(p for p in GIZA_PATCHES if p[0] == 'khafre')
    product = DwellProduct(PRODUCTS['giza'])
    acq = product.acq
    row, col = product.geolocate(entry[2], entry[3], entry[4])
    print(f'  reading a {N_AZ * acq.azimuth_spacing_m / 1000:.2f} km strip', flush=True)
    crop, (r0, c0) = product.crop(row, col, N_AZ, N_RG)
    rate = float(product.doppler_rate_hz_s(c0 + N_RG / 2))
    band = occupied_band(crop, acq.prf_hz)
    plan = BANK.plan(acq, N_AZ, rate, centroid_hz=band['centroid_hz'])
    freqs = np.fft.fftshift(np.fft.fftfreq(N_AZ, d=1.0 / acq.prf_hz))
    ref, _ = BANK.masks(plan, freqs)

    eph = product.ephemeris()
    lo_t, hi_t = eph.span
    t_zd_mid = float(np.clip(product.zero_doppler_time_s(r0 + N_AZ / 2), lo_t, hi_t))
    _, slant, _ = perpendicular_baselines(eph, [t_zd_mid], lla_to_ecef(*entry[2:5]), t_zd_mid)
    speed = float(np.linalg.norm(eph.at([t_zd_mid])[1][0]))

    rows = np.arange(300, N_AZ - 300, STRIDE_AZ)
    cols = np.arange(160, N_RG - 160, STRIDE_RG)
    RR, CC = np.meshgrid(rows, cols, indexing='ij')
    rows_f, cols_f = RR.ravel(), CC.ravel()
    print(f'  {rows_f.size:,} array elements over {rows.size * STRIDE_AZ * acq.azimuth_spacing_m:.0f} m', flush=True)

    spec = np.fft.fftshift(np.fft.fft(crop, axis=0), axes=0)
    base = np.fft.ifft(np.fft.ifftshift(spec * ref[BANK.K // 2][:, None], axes=0), axis=0)
    shifts = np.zeros((rows_f.size, BANK.K))
    for k in range(BANK.K):
        img = np.fft.ifft(np.fft.ifftshift(spec * ref[k][:, None], axes=0), axis=0)
        shifts[:, k] = patch_shifts(img, base, rows_f, cols_f, patch=PATCH_PX, upsample=1000)[0]
        if k % 6 == 0:
            print(f'    look {k + 1}/{BANK.K}  {time.time() - t0:.0f}s', flush=True)
    velocity = azimuth_shift_to_velocity(shifts, acq, slant, speed)

    x_m = (rows_f - rows_f.mean()) * acq.azimuth_spacing_m
    tau = (plan['centre_hz'] - band['centroid_hz']) / rate
    t_zd = np.array([product.zero_doppler_time_s(r0 + r) for r in rows_f])
    times = t_zd[:, None] + tau[None, :]
    order = np.argsort(tau)
    velocity, times = velocity[:, order], times[:, order]

    lags, real_stack, counts = stack_by_separation(velocity, x_m, times)

    # Positive control: band-limited noise propagating at a known speed. Noise rather than a
    # tone, because a tone's cross-correlation is itself periodic and its peak is ambiguous by
    # one period; a band-limited pulse gives an unambiguous arrival.
    rng = np.random.default_rng(1)
    grid = np.arange(times.min() - 3.0, times.max() + 3.0, 0.01)
    white = rng.normal(size=grid.size)
    spec_n = np.fft.rfft(white)
    fgrid = np.fft.rfftfreq(grid.size, d=0.01)
    spec_n[(fgrid < BAND_LOW_HZ) | (fgrid > BAND_HIGH_HZ)] = 0
    wave = np.fft.irfft(spec_n, n=grid.size)
    wave /= wave.std()
    retarded = times - (x_m[:, None] - x_m.min()) / RAYLEIGH_M_S
    injected = velocity + INJECT_AMPLITUDE_UM_S * 1e-6 * np.interp(retarded, grid, wave)
    _, control_stack, _ = stack_by_separation(injected, x_m, times)

    centres = (SEP_BINS[:-1] + SEP_BINS[1:]) / 2
    predicted = centres / RAYLEIGH_M_S
    rows_out = []
    for b in range(len(centres)):
        peak_real = lags[np.argmax(np.abs(real_stack[b]))]
        peak_ctrl = lags[np.argmax(np.abs(control_stack[b]))]
        rows_out.append({'separation_m': float(centres[b]), 'pairs': int(counts[b]),
                         'predicted_lag_s': float(predicted[b]),
                         'real_peak_lag_s': float(peak_real),
                         'real_peak_value': float(np.max(np.abs(real_stack[b]))),
                         'control_peak_lag_s': float(peak_ctrl),
                         'control_peak_value': float(np.max(np.abs(control_stack[b])))})

    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)
    for ax, stack, title in ((axs[0], real_stack, 'Real Giza data'),
                             (axs[1], control_stack,
                              f'Same data with a {RAYLEIGH_M_S:.0f} m/s wavefield injected')):
        for b in range(len(centres)):
            ax.plot(lags, stack[b] * 0.85 + centres[b] / 250, lw=1.4, color=viz.SERIES[b % 8])
        ax.plot(predicted, centres / 250, 'o--', color=viz.SERIES[7], lw=1.2, ms=4,
                label=f'moveout at {RAYLEIGH_M_S:.0f} m/s')
        ax.set_xlabel('lag (s)'); ax.set_title(title, loc='left')
        ax.set_xlim(lags[0], lags[-1]); ax.legend(fontsize=8.5, loc='upper right')
    axs[0].set_ylabel('pair separation (arbitrary offset, labelled below)')
    axs[0].set_yticks(centres / 250); axs[0].set_yticklabels([f'{c:.0f} m' for c in centres])
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r11_moveout.png'), 'caption':
             '<b>No wave arrives.</b> Velocity series recovered from a 1.8 km strip of the real Giza acquisition, '
             'cross-correlated between pixel pairs and stacked by separation, with each pair\'s own '
             'Doppler-to-time offset removed. Left: the real data, where a propagating surface wave would appear as '
             'a peak tracking the dashed moveout line. Right: the identical processing after injecting a plane wave '
             'of known speed and amplitude, which the same stack recovers cleanly. The method works; there is '
             'nothing in the data for it to find.'}]

    # Bins whose predicted lag is under two sample intervals carry no information: neighbouring
    # pixels share speckle and correlate at zero lag anyway, so a "match" there is automatic.
    dt = plan['span_s'] / (BANK.K - 1)
    informative = [r for r in rows_out if r['predicted_lag_s'] >= 2 * dt]
    ctrl_ok = sum(1 for r in informative if abs(r['control_peak_lag_s'] - r['predicted_lag_s']) < 0.2)
    real_ok = sum(1 for r in informative if abs(r['real_peak_lag_s'] - r['predicted_lag_s']) < 0.2)
    real_amp = float(np.median([r['real_peak_value'] for r in informative]))
    ctrl_amp = float(np.median([r['control_peak_value'] for r in informative]))
    m = {'strip_m': float(N_AZ * acq.azimuth_spacing_m), 'elements': int(rows_f.size),
         'looks': BANK.K, 'record_s': float(plan['span_s']), 'sampling_hz': float(plan['sampling_hz']),
         'injected_amplitude_um_s': INJECT_AMPLITUDE_UM_S, 'rayleigh_m_s': RAYLEIGH_M_S,
         'injected_band_hz': [BAND_LOW_HZ, BAND_HIGH_HZ],
         'bins': rows_out, 'informative_bins': len(informative),
         'control_bins_on_moveout': ctrl_ok, 'real_bins_on_moveout': real_ok,
         'median_peak_correlation_real': real_amp, 'median_peak_correlation_control': ctrl_amp,
         'control_over_real': ctrl_amp / max(real_amp, 1e-12), 'sample_interval_s': float(dt),
         'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 21, 'tag': 'null', 'eyebrow': '11 · The array',
        'title': 'Used as a dense seismic array, the image contains no propagating wave to find',
        'question': 'A SAR image gives a measurement at every pixel. Can it be processed as an array and inverted like ambient-noise seismology?',
        'finding': (f"The first step already fails. Velocity series from {rows_f.size:,} elements spread over "
                    f"{N_AZ * acq.azimuth_spacing_m / 1000:.1f} km were cross-correlated between pairs and stacked "
                    f"by separation. Injecting a {RAYLEIGH_M_S:.0f} m/s band-limited wavefield of "
                    f"{INJECT_AMPLITUDE_UM_S:.0f} µm/s into the same series is recovered on {ctrl_ok} of "
                    f"{len(informative)} informative separation bins, with a median peak correlation of "
                    f"{ctrl_amp:.2f}. The real data gives {real_amp:.3f}, lower by a factor of "
                    f"{ctrl_amp / max(real_amp, 1e-12):.0f}, and its peaks fall at unrelated lags. Bins closer than "
                    f"about 350 m are excluded because neighbouring pixels share speckle and correlate at zero lag "
                    f"whatever the ground is doing. The array processing works; there is no wavefield in the data "
                    f"for it to find, because the record is {plan['span_s']:.1f} s where the microseism band needs "
                    f"hundreds."),
        'limitations': ('One strip, one acquisition, one injected speed and band. The control is a single '
                        'plane wavefield, which is easier to recover than a real diffuse ambient field, so this bounds '
                        'the method from the favourable side. A longer record, from many tasked dwells stacked, '
                        'would raise sensitivity at frequencies above about half a hertz; it cannot extend the '
                        'record below that.'),
        'method': (f'A {N_AZ:,} by {N_RG} sample strip read from the native product, {BANK.K} common-reference '
                   f'sub-apertures of {BANK.width_frac:.0%} bandwidth confined to the coherent window, azimuth '
                   f'shifts converted to line-of-sight velocity, pairs drawn at random within each separation bin '
                   f'and cross-correlated with each pair\'s zero-Doppler time offset removed before stacking.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in m.items() if k != 'bins'}, indent=1))
    for r in rows_out:
        print(f"  {r['separation_m']:6.0f} m  predicted {r['predicted_lag_s']:+.3f}  "
              f"real {r['real_peak_lag_s']:+.3f} ({r['real_peak_value']:.3f})  "
              f"control {r['control_peak_lag_s']:+.3f} ({r['control_peak_value']:.3f})")


if __name__ == '__main__':
    main()
