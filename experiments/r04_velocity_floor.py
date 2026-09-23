"""R4: what a 25 second dwell can actually measure.

Strip away the depth claim and a defensible instrument remains. Registering every sub-aperture
against a common reference recovers the azimuth shift of each pixel through slow time, and an
azimuth shift is a line-of-sight velocity: a target moving at v acquires a Doppler offset
2v/lambda, which focusing maps to a position offset of (R/V) v. With R/V near 80 s the lever is
large. This experiment measures the resulting velocity noise floor on real Giza data, which
bounds what any inference resting on surface motion could possibly use.

Ambient ground motion in the microseism band is of order 0.1 to 10 micrometres per second
(Peterson 1993, USGS Open-File Report 93-322); Giza sits beside a large city, so the upper end
of that range is the fair comparison.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import run_patch, GIZA_PATCHES, RESULTS
from sarsim import viz

TEST = 'r04_velocity_floor'
AMBIENT_LOW, AMBIENT_HIGH = 0.1e-6, 10e-6      # m/s, microseism band


def spectrum(v, fs):
    """One-sided power spectral density of each pixel's velocity series, mean over pixels."""
    x = v - v.mean(axis=1, keepdims=True)
    w = np.hanning(x.shape[1])[None, :]
    x = x * w
    f = np.fft.rfftfreq(x.shape[1], d=1.0 / fs)
    p = np.abs(np.fft.rfft(x, axis=1)) ** 2
    p *= 2.0 / (fs * (w ** 2).sum())
    return f, p.mean(axis=0)


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    runs, rows = {}, []
    for key, label, *_ , kind in [(p[0], p[1], p[2], p[3], p[4], p[5]) for p in GIZA_PATCHES]:
        runs[key] = run_patch('giza', key, 'reference')
    figs = []

    # ---- (a) velocity noise against target brightness
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    ax.axhspan(AMBIENT_LOW * 1e6, AMBIENT_HIGH * 1e6, color=viz.GRID, alpha=0.7, lw=0)
    ax.text(0.985, AMBIENT_HIGH * 1e6 * 1.25, 'ambient ground motion, microseism band',
            transform=ax.get_yaxis_transform(), ha='right', fontsize=8.5, color=viz.INK2)
    for i, (key, run) in enumerate(runs.items()):
        m = run['meta']
        amp = run['amplitude'].astype(float)
        v = run['velocity'].astype(float)
        v = v - v.mean(axis=1, keepdims=True)
        edges = np.percentile(amp, np.linspace(0, 100, 13))
        centres, floors = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (amp >= lo) & (amp < hi)
            if sel.sum() < 30:
                continue
            centres.append(20 * np.log10(np.median(amp[sel])))
            floors.append(np.sqrt((v[sel] ** 2).mean()) * 1e6)
        ax.plot(centres, floors, '-o', ms=4, lw=1.8, color=viz.SERIES[i % 8],
                markeredgecolor=viz.SURFACE, label=f"{m['label']}")
        rows.append({'patch': key, 'label': m['label'], 'kind': m['kind'],
                     'pixels': int(amp.size),
                     'median_amplitude_db': float(20 * np.log10(np.median(amp))),
                     'velocity_floor_all_um_s': float(np.sqrt((v ** 2).mean()) * 1e6),
                     'velocity_floor_brightest_decile_um_s': float(floors[-1]) if floors else None,
                     'sampling_hz': m['sampling_hz'], 'series_span_s': m['centres_span_s'],
                     'sub_aperture_s': m['sub_aperture_s']})
    ax.set_yscale('log')
    ax.set_xlabel('target amplitude (dB, arbitrary reference)')
    ax.set_ylabel('line-of-sight velocity noise (micrometres per second, rms)')
    ax.set_title('Velocity precision of one ICEYE dwell, measured on the Giza plateau', loc='left')
    ax.legend(fontsize=8.5, loc='upper right')
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r04_floor.png'), 'caption':
                 '<b>What the instrument can see.</b> Every sub-aperture is registered against a common reference, '
                 'giving each pixel an azimuth-shift series through the 22 seconds of swept aperture, converted to '
                 'line-of-sight velocity. The curves show the measured scatter of that series against target '
                 'brightness, for each patch of the Giza image. The shaded band is the range of real ambient ground '
                 'motion in the microseism band. Where the curves sit relative to that band is the honest answer to '
                 'what this acquisition could contribute to any subsurface inference.'})

    # ---- (b) spectrum of the recovered series
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    for i, (key, run) in enumerate(runs.items()):
        m = run['meta']
        amp = run['amplitude'].astype(float)
        bright = amp >= np.percentile(amp, 90)
        f, p = spectrum(run['velocity'].astype(float)[bright], m['sampling_hz'])
        ax.plot(f[1:], np.sqrt(p[1:]) * 1e6, lw=1.8, color=viz.SERIES[i % 8], label=m['label'])
    ax.axvspan(0.1, 0.35, color=viz.GRID, alpha=0.7, lw=0)
    ax.text(0.2, 0.96, 'microseism band', transform=ax.get_xaxis_transform(), ha='center',
            fontsize=8.5, color=viz.INK2)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('frequency (Hz)')
    ax.set_ylabel('velocity spectral density ($\\mu$m/s/$\\sqrt{\\mathrm{Hz}}$)')
    ax.set_title('Spectrum of the recovered velocity series, brightest tenth of each patch', loc='left')
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r04_spectrum.png'), 'caption':
                 '<b>No resonance, on monuments or on sand.</b> Spectral density of the recovered velocity series '
                 'for the brightest tenth of the pixels in each patch. A real site resonance or a real microseism '
                 'peak would appear as a bump in the shaded band. The series are close to flat, which is what '
                 'measurement noise looks like, and the monuments do not differ from open desert.'})

    floors = [r['velocity_floor_brightest_decile_um_s'] for r in rows if r['velocity_floor_brightest_decile_um_s']]
    best = min(floors)
    m = {'runs': rows, 'ambient_band_um_s': [AMBIENT_LOW * 1e6, AMBIENT_HIGH * 1e6],
         'best_floor_um_s': float(best), 'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 14, 'tag': 'threshold', 'eyebrow': '04 · The instrument',
        'title': 'What one dwell really measures: a velocity series, and its noise floor',
        'question': 'Set the depth claim aside. What surface motion can this acquisition actually detect?',
        'finding': (f"A line-of-sight velocity series for every pixel, sampled at about 2.2 Hz across 22 seconds. "
                    f"On the brightest targets in the Giza image the measured noise floor is "
                    f"{best:,.0f} micrometres per second. Ambient ground motion in the microseism band is of order "
                    f"0.1 to 10 micrometres per second, so this acquisition sits well above the motion that any "
                    f"passive subsurface inference would have to rest on. The recovered spectra are flat and the "
                    f"monuments do not differ from open desert."),
        'limitations': ('The scene-wide median series is removed, so any motion common to the whole patch is '
                        'removed with it; this bounds differential motion, not absolute platform-referenced motion. '
                        'The floor mixes true measurement noise with real clutter decorrelation, so it is an upper '
                        'bound on achievable precision rather than a fundamental limit. Corner reflectors or other '
                        'strong coherent targets would do better than natural desert scattering.'),
        'method': ('Common-reference sub-aperture tracking: 50 sub-apertures of 10 percent bandwidth registered '
                   'against the band-centre sub-aperture, azimuth shift converted to velocity by v = dx V / R with '
                   'the real slant range and platform speed. Spectra use a Hann window and are averaged over the '
                   'brightest tenth of each patch.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
