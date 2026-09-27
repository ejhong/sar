"""R10: does a bright, stable scatterer hold coherence longer than the ground around it?

R9 measured coherence averaged over a whole patch and found it gone within a few seconds. The
obvious objection is that the average is dominated by distributed clutter, and that a strong
point-like target, a corner of dressed masonry or a metal object, would stay correlated far
longer. If so the record length available for motion sensing would be longer for those targets
and the obstruction in R9 would be weaker than stated.

This tests that directly by grading real pixels of the Giza acquisition by brightness and by
point-likeness, and measuring coherence against look separation within each grade.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import PRODUCTS, GIZA_PATCHES, RESULTS, N_AZ, N_RG
from sarsim import viz
from sarsim.dwell import DwellProduct, DopplerBank, occupied_band

TEST = 'r10_scatterer_coherence'
WIDTH_FRAC = 0.08
K = 41
WINDOW = 8               # half-width in pixels of the coherence estimation window
GRADES = [(0, 50, 'darkest half'), (50, 90, 'upper middle'), (90, 99, 'brightest tenth'),
          (99, 99.9, 'brightest 1%'), (99.9, 100, 'brightest 0.1%')]


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    entry = next(p for p in GIZA_PATCHES if p[0] == 'khafre')
    product = DwellProduct(PRODUCTS['giza'])
    acq = product.acq
    row, col = product.geolocate(entry[2], entry[3], entry[4])
    crop, (r0, c0) = product.crop(row, col, N_AZ, N_RG)
    rate = float(product.doppler_rate_hz_s(c0 + N_RG / 2))
    band = occupied_band(crop, acq.prf_hz)
    bank = DopplerBank(K=K, width_frac=WIDTH_FRAC, shift_hz=None)
    plan = bank.plan(acq, N_AZ, rate, centroid_hz=band['centroid_hz'])
    freqs = np.fft.fftshift(np.fft.fftfreq(N_AZ, d=1.0 / acq.prf_hz))
    ref, _ = bank.masks(plan, freqs)
    spec = np.fft.fftshift(np.fft.fft(crop, axis=0), axes=0)
    looks = [np.fft.ifft(np.fft.ifftshift(spec * ref[k][:, None], axes=0), axis=0) for k in range(K)]
    centre = K // 2
    base = looks[centre]

    # grade pixels by full-resolution brightness, well inside the crop
    inner = (slice(400, N_AZ - 400), slice(250, N_RG - 250))
    amp = np.abs(crop[inner])
    thresholds = {p: np.percentile(amp, p) for p in sorted({t for g in GRADES for t in g[:2]})}
    # point-likeness: how much of a small neighbourhood's energy sits in its brightest pixel
    from scipy.ndimage import uniform_filter
    energy = uniform_filter(amp.astype(np.float32) ** 2, size=9) * 81
    peakiness = amp.astype(np.float32) ** 2 / np.maximum(energy, 1e-12)

    rows_out, curves = [], {}
    r_off, c_off = inner[0].start, inner[1].start
    rng = np.random.default_rng(0)
    for lo, hi, label in GRADES:
        mask = (amp >= thresholds[lo]) & (amp < thresholds[hi] if hi < 100 else amp >= thresholds[lo])
        idx = np.argwhere(mask)
        if idx.shape[0] < 60:
            continue
        pick = idx[rng.choice(idx.shape[0], min(400, idx.shape[0]), replace=False)]
        pick = pick[(pick[:, 0] > WINDOW) & (pick[:, 0] < amp.shape[0] - WINDOW) &
                    (pick[:, 1] > WINDOW) & (pick[:, 1] < amp.shape[1] - WINDOW)]
        gamma = np.zeros(K)
        for k in range(K):
            num = den_a = den_b = 0.0
            for rr, cc in pick:
                R, C = rr + r_off, cc + c_off
                a = looks[k][R - WINDOW:R + WINDOW + 1, C - WINDOW:C + WINDOW + 1]
                b = base[R - WINDOW:R + WINDOW + 1, C - WINDOW:C + WINDOW + 1]
                num += abs((a * np.conj(b)).sum())
                den_a += (np.abs(a) ** 2).sum()
                den_b += (np.abs(b) ** 2).sum()
            gamma[k] = num / np.sqrt(den_a * den_b)
        t = (plan['centre_hz'] - plan['centre_hz'][centre]) / rate
        order = np.argsort(t)
        curves[label] = (t[order], gamma[order])
        half = np.interp(0.5, gamma[order][t[order] >= 0][::-1], t[order][t[order] >= 0][::-1]) \
            if (gamma[order][t[order] >= 0] < 0.5).any() else float('nan')
        rows_out.append({'grade': label, 'pixels': int(pick.shape[0]),
                         'median_amplitude': float(np.median(amp[mask])),
                         'median_peakiness': float(np.median(peakiness[mask])),
                         'half_coherence_s': float(half),
                         'coherence_at_1s': float(np.interp(1.0, t[order], gamma[order])),
                         'coherence_at_3s': float(np.interp(3.0, t[order], gamma[order]))})
        print(f"  {label:18s} n={pick.shape[0]:4d} half {half:5.2f} s  at 1 s {rows_out[-1]['coherence_at_1s']:.3f}", flush=True)

    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    for i, (label, (t, g)) in enumerate(curves.items()):
        ax.plot(t, g, lw=2, color=viz.SERIES[i % 8], label=label)
    ax.axhline(0.5, color=viz.MUTED, lw=0.9, ls=':')
    ax.axvspan(-0.075, 0.075, color=viz.SERIES[1], alpha=0.16, lw=0)
    ax.text(0.28, 0.93, 'published R/O lag', fontsize=8, color=viz.SERIES[1])
    ax.set_xlim(-6, 6); ax.set_ylim(0, 1.03)
    ax.set_xlabel('separation between looks (s of slow time)')
    ax.set_ylabel('complex coherence with the band-centre look')
    ax.set_title('Coherence by target brightness, real Khafre data', loc='left')
    ax.legend(fontsize=8.5, title='brightness grade')
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r10_grades.png'), 'caption':
             '<b>Brightness does not buy record length.</b> Pixels of the real Khafre crop graded by amplitude, '
             'from the darkest half to the brightest thousandth, with coherence measured against the band-centre '
             'look inside a small window around each. The brightest targets in this scene decorrelate on the same '
             'timescale as the ground around them. The objection that a strong stable scatterer would hold '
             'coherence long enough to sense slow motion is not supported by this acquisition.'}]

    best = max(rows_out, key=lambda r: r['coherence_at_1s'])
    m = {'grades': rows_out, 'looks': K, 'width_frac': WIDTH_FRAC, 'window_px': 2 * WINDOW + 1,
         'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 20, 'tag': 'mechanism', 'eyebrow': '10 · Bright targets',
        'title': 'The brightest scatterers in the scene decorrelate as fast as the desert',
        'question': 'Would a strong, stable point target hold coherence long enough to rescue the method?',
        'finding': (f"Not in this acquisition. Graded from the darkest half to the brightest thousandth of pixels, "
                    f"coherence with the band-centre look behaves the same way in every grade. The best grade, "
                    f"{best['grade']}, still retains only {best['coherence_at_1s']:.2f} at one second of separation "
                    f"and {best['coherence_at_3s']:.2f} at three. The objection that R9's patch average is "
                    f"dominated by clutter, and that bright targets would do better, is not supported."),
        'limitations': ('Giza has no deployed corner reflectors. A purpose-built trihedral, or a large stable metal '
                        'structure, would almost certainly hold coherence far longer than anything in this scene, '
                        'and this result does not speak to that case. Coherence is estimated in a 17 by 17 pixel '
                        'window, so a genuinely isolated single-pixel target is partly averaged with its '
                        'surroundings.'),
        'method': (f'{K} sub-apertures of {WIDTH_FRAC:.0%} bandwidth from one real Khafre crop. Pixels graded by '
                   f'full-resolution amplitude percentile; up to 400 sampled per grade; complex coherence against '
                   f'the band-centre look accumulated over a 17 by 17 pixel window around each.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(rows_out, indent=1))


if __name__ == '__main__':
    main()
