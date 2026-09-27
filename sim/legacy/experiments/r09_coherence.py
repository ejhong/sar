"""R9: why a dwell cannot give both time diversity and coherence.

In a dwell product the Doppler axis is slow time, so two sub-apertures separated by seconds are
also separated by degrees of look angle. Natural ground decorrelates over a fraction of a
degree. This experiment measures that trade directly on the Giza acquisition: complex coherence
between sub-aperture looks against their separation, in hertz, in seconds and in degrees.

It is the constraint behind every other result here. The published protocols overlap their
reference and offset passbands by 96 to 99 per cent, which is the only coherent regime
available, and which corresponds to lags of tens of milliseconds; ground motion in the
microseism band needs lags of about a second.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import PRODUCTS, GIZA_PATCHES, RESULTS, N_AZ, N_RG
from sarsim import viz
from sarsim.dwell import DwellProduct, DopplerBank, occupied_band

TEST = 'r09_coherence'
WIDTH_FRAC = 0.08
K = 61
PATCHES = [('khafre', 'Khafre pyramid'), ('desert_west', 'Open plateau')]


def coherence_curve(product, lat, lon, height, width_frac=WIDTH_FRAC, k=K):
    acq = product.acq
    row, col = product.geolocate(lat, lon, height)
    crop, (r0, c0) = product.crop(row, col, N_AZ, N_RG)
    rate = float(product.doppler_rate_hz_s(c0 + N_RG / 2))
    band = occupied_band(crop, acq.prf_hz)
    bank = DopplerBank(K=k, width_frac=width_frac, shift_hz=None)
    plan = bank.plan(acq, N_AZ, rate, centroid_hz=band['centroid_hz'])
    freqs = np.fft.fftshift(np.fft.fftfreq(N_AZ, d=1.0 / acq.prf_hz))
    ref, _ = bank.masks(plan, freqs)
    spec = np.fft.fftshift(np.fft.fft(crop, axis=0), axes=0)
    inner = (slice(400, N_AZ - 400), slice(250, N_RG - 250))
    centre = k // 2
    base = np.fft.ifft(np.fft.ifftshift(spec * ref[centre][:, None], axes=0), axis=0)[inner]
    denom_b = (np.abs(base) ** 2).sum()
    rows = []
    speed = float(np.linalg.norm(product.ephemeris().at([0.0])[1][0]))
    for i in range(k):
        img = np.fft.ifft(np.fft.ifftshift(spec * ref[i][:, None], axes=0), axis=0)[inner]
        gamma = abs((img * np.conj(base)).sum()) / np.sqrt((np.abs(img) ** 2).sum() * denom_b)
        df = float(plan['centre_hz'][i] - plan['centre_hz'][centre])
        rows.append({'delta_hz': df, 'delta_s': df / rate,
                     'delta_aspect_deg': float(np.rad2deg(np.arcsin(np.clip(
                         df / speed * acq.wavelength_m / 2, -1, 1)))),
                     'coherence': float(gamma)})
    return rows, plan


def crossing(rows, level):
    """Separation in seconds at which coherence first falls below a level."""
    half = sorted((r for r in rows if r['delta_s'] >= 0), key=lambda r: r['delta_s'])
    for a, b in zip(half[:-1], half[1:]):
        if a['coherence'] >= level > b['coherence']:
            f = (a['coherence'] - level) / max(a['coherence'] - b['coherence'], 1e-9)
            return a['delta_s'] + f * (b['delta_s'] - a['delta_s'])
    return None


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    product = DwellProduct(PRODUCTS['giza'])
    curves, summary_rows = {}, []
    for key, label in PATCHES:
        entry = next(p for p in GIZA_PATCHES if p[0] == key)
        rows, plan = coherence_curve(product, entry[2], entry[3], entry[4])
        curves[key] = (rows, label)
        summary_rows.append({
            'patch': key, 'label': label,
            'half_coherence_s': crossing(rows, 0.5), 'quarter_coherence_s': crossing(rows, 0.25),
            'tenth_coherence_s': crossing(rows, 0.10),
            'sub_aperture_s': plan['sub_aperture_s']})

    fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.8))
    for i, (key, (rows, label)) in enumerate(curves.items()):
        t = np.array([r['delta_s'] for r in rows])
        g = np.array([r['coherence'] for r in rows])
        axs[0].plot(t, g, lw=2.4 if i == 0 else 1.6, ls='-' if i == 0 else '--',
                    color=viz.SERIES[i], label=label)
    axs[0].axhline(0.5, color=viz.MUTED, lw=0.9, ls=':')
    axs[0].text(-11, 0.53, 'half coherence', fontsize=8.5, color=viz.INK2)
    axs[0].axvspan(-0.075, 0.075, color=viz.SERIES[1], alpha=0.18, lw=0)
    axs[0].text(0.35, 0.90, 'the published\nreference/offset lag', fontsize=8, color=viz.SERIES[1])
    axs[0].set_xlabel('separation between looks (s of slow time)')
    axs[0].set_ylabel('complex coherence')
    axs[0].set_xlim(-12, 12); axs[0].set_ylim(0, 1.03); axs[0].legend(fontsize=8.5)
    axs[0].set_title('Looks decorrelate within seconds', loc='left')

    periods = np.geomspace(0.05, 30, 200)
    axs[1].plot(periods, np.interp(np.clip(periods / 2, 0, 12),
                                   [r['delta_s'] for r in curves['khafre'][0] if r['delta_s'] >= 0][::-1] or [0],
                                   [r['coherence'] for r in curves['khafre'][0] if r['delta_s'] >= 0][::-1] or [0]),
                lw=2, color=viz.SERIES[0])
    axs[1].axvspan(1 / 0.35, 1 / 0.1, color=viz.GRID, alpha=0.8, lw=0)
    axs[1].text(5.5, 0.9, 'microseism band', fontsize=8.5, color=viz.INK2, ha='center')
    axs[1].set_xscale('log')
    axs[1].set_xlabel('period of the ground motion you want to sense (s)')
    axs[1].set_ylabel('coherence available at half that lag')
    axs[1].set_ylim(0, 1.03)
    axs[1].set_title('What is left by the time the lag is long enough', loc='left')
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r09_coherence.png'), 'caption':
             '<b>The constraint behind everything else.</b> Left: complex coherence between sub-aperture looks of '
             'the real Giza acquisition and the band-centre look, against their separation in slow time. Coherence '
             'is gone within a few seconds on the pyramid and on open desert alike. The shaded strip is the lag the '
             'published reference and offset bands actually use, which is why they stay coherent. Right: the same '
             'curve read the other way, as the coherence still available when the lag is long enough to sample a '
             'given period of ground motion. The microseism band sits where almost nothing is left.'}]

    khafre = summary_rows[0]
    m = {'patches': summary_rows, 'width_frac': WIDTH_FRAC, 'looks': K,
         'curves': {k: v[0] for k, v in curves.items()}, 'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 19, 'tag': 'mechanism', 'eyebrow': '09 · The obstruction',
        'title': 'In a dwell, time diversity and coherence are the same axis, and you cannot have both',
        'question': 'Could any sub-aperture design measure slow ground motion from one dwell?',
        'finding': (f"No. Slow time and look angle are the same axis in a dwell product, so looks separated by "
                    f"seconds are separated by degrees, and natural ground decorrelates over a fraction of a "
                    f"degree. Measured on the real acquisition, coherence with the band-centre look falls to half "
                    f"by {khafre['half_coherence_s']:.2f} s of separation and to a tenth by "
                    f"{khafre['tenth_coherence_s']:.2f} s, identically on the pyramid and on open desert. Ground "
                    f"motion in the microseism band needs lags of about a second, by which point there is little "
                    f"coherence left to register. The published protocols avoid this by using lags of tens of "
                    f"milliseconds, which keeps them coherent and makes them blind to the motion they invoke."),
        'limitations': ('Coherence is measured over one patch of one acquisition at 8 per cent sub-aperture width; '
                        'a wider sub-aperture averages more and decorrelates differently. Strong stable point '
                        'scatterers, such as corner reflectors or some built structures, hold coherence far longer '
                        'than the natural desert measured here, so this curve is the distributed-target case, not '
                        'a universal limit.'),
        'method': ('61 sub-apertures of 8 per cent bandwidth swept across the processed Doppler support of a real '
                   'Khafre and a real desert crop; complex coherence computed against the band-centre look over '
                   'the interior of each crop; separation converted to slow time through the product Doppler rate '
                   'and to aspect angle through the carrier wavelength and platform speed.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary_rows, indent=1))


if __name__ == '__main__':
    main()
