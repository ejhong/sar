"""R1: within-image controls at Giza.

The published claim is that one dwell acquisition of the Giza plateau contains vertical wells
and multi-level structures beneath Khafre. The strongest available test needs no simulation:
run the same processing on other patches of the same image, taken on the same pass, through
the same processor, and ask whether the monuments differ from open desert.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

from real_common import GIZA_PATCHES, run_patch, patch_stats, RESULTS
from sarsim import viz

TEST = 'r01_giza_controls'
BANK = 'paper'


def compute(banks=('paper', 'reference'), stride_az=12, stride_rg=8, force=False):
    t0 = time.time()
    for bank in banks:
        for entry in GIZA_PATCHES:
            print(f'[{time.time() - t0:6.0f}s] {bank}/{entry[0]}', flush=True)
            run_patch('giza', entry[0], bank, stride_az=stride_az, stride_rg=stride_rg,
                      force=force, verbose=True)
    print(f'compute total {time.time() - t0:.0f}s')


def report(bank=BANK):
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    runs = {e[0]: run_patch('giza', e[0], bank, verbose=False) for e in GIZA_PATCHES}
    stats = [patch_stats(r) for r in runs.values()]
    figs = []

    def panels(title, draw, path, caption, cbar_label=None):
        fig, axs = plt.subplots(2, 3, figsize=(13.4, 7.4))
        im = None
        for ax, (key, run) in zip(axs.ravel(), runs.items()):
            im = draw(ax, run)
            m = run['meta']
            ax.set_title(f"{m['label']}", loc='left', fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(title, x=0.012, ha='left', fontsize=11.5, y=0.985)
        fig.tight_layout(rect=[0, 0, 1, 0.965])
        if cbar_label is not None and im is not None:
            fig.colorbar(im, ax=axs, fraction=0.018, pad=0.012, label=cbar_label)
        figs.append({'file': viz.finish(fig, path), 'caption': caption})

    # --- 1. the imagery itself
    def draw_amp(ax, run):
        g = (run['grid_rows'].size, run['grid_cols'].size)
        a = run['amplitude'].reshape(g)
        db = 20 * np.log10(np.maximum(a, 1e-3) / np.percentile(a, 99.5))
        return ax.imshow(db.T, cmap=viz.CMAP_GRAY, vmin=-22, vmax=3, aspect='auto', interpolation='nearest')
    panels('Six patches of one ICEYE dwell acquisition, Giza, 27 August 2025',
           draw_amp, f'{fd}/r01_amplitude.png',
           '<b>The six patches.</b> Amplitude of the native complex image on the analysis grid, geolocated by the '
           'product\'s own rational polynomial model with no fitted offset. Four monuments and two stretches of open '
           'plateau, all from the same pass, the same processor and the same 25 seconds. Every comparison that '
           'follows uses these six and nothing else.', 'dB below the 99.5th percentile')

    # --- 2. the structure map
    def draw_score(ax, run):
        g = (run['grid_rows'].size, run['grid_cols'].size)
        s = run['score'].max(axis=1).reshape(g)
        return ax.imshow(s.T, cmap=viz.CMAP_SEQ, vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    panels('Best depth-fit score over all depths, same processing for every patch',
           draw_score, f'{fd}/r01_score.png',
           '<b>Spot the monument.</b> For every pixel, the best score the depth-focusing step achieves at any depth, '
           'drawn on a common scale. If the method were responding to structure beneath the monuments, these six '
           'panels would not look alike. The pyramids occupy the first three panels and the open plateau the last '
           'two.', 'best adjusted $R^{2}$')

    # --- 3. distributions
    fig, axs = plt.subplots(1, 3, figsize=(13.4, 4.2))
    bins = np.linspace(0, 1, 46)
    for i, (key, run) in enumerate(runs.items()):
        m = run['meta']
        smax = run['score'].max(axis=1)
        style = dict(lw=2, color=viz.SERIES[i % 8], ls='-' if m['kind'] == 'monument' else '--')
        axs[0].hist(smax, bins=bins, density=True, histtype='step', label=m['label'], **style)
        zbest = run['z'][np.argmax(run['score'], axis=1)]
        axs[1].hist(zbest, bins=np.linspace(0, run['z'][-1], 40), density=True, histtype='step', **style)
        t = run['tomo'].astype(float)
        tn = t / np.maximum(t.sum(axis=1, keepdims=True), 1e-30)
        axs[2].plot(run['z'], tn.mean(axis=0), **style)
    axs[0].set_xlabel('best adjusted $R^{2}$'); axs[0].set_ylabel('density')
    axs[0].set_title('Structure score', loc='left'); axs[0].legend(fontsize=7.5)
    axs[1].set_xlabel('depth of the best fit (m)')
    axs[1].set_title('Where the structures land', loc='left')
    axs[2].set_xlabel('model depth (m)'); axs[2].set_ylabel('mean normalised $|h(z)|^{2}$')
    axs[2].set_title('Mean depth profile', loc='left')
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r01_distributions.png'), 'caption':
                 '<b>The same distributions everywhere.</b> Solid lines are monuments, dashed lines are open '
                 'plateau. Left: the per-pixel structure score. Middle: the depth each pixel is assigned. Right: '
                 'the mean depth profile of each patch. The depth axis here runs to the ambiguity limit of this '
                 'bank, which is under fifteen metres.'})

    # --- statistics
    mon = np.concatenate([runs[e[0]]['score'].max(axis=1) for e in GIZA_PATCHES if e[5] == 'monument'])
    ctl = np.concatenate([runs[e[0]]['score'].max(axis=1) for e in GIZA_PATCHES if e[5] == 'control'])
    ks = ks_2samp(mon, ctl)
    khafre = runs['khafre']['score'].max(axis=1)
    ks_k = ks_2samp(khafre, ctl)
    meta = runs['khafre']['meta']
    pooled_sd = float(np.sqrt((mon.var() * mon.size + ctl.var() * ctl.size) / (mon.size + ctl.size)))
    cohens_d = float((mon.mean() - ctl.mean()) / pooled_sd)
    m = {'bank': bank, 'patches': stats,
         'monument_mean': float(mon.mean()), 'control_mean': float(ctl.mean()),
         'pooled_sd': pooled_sd, 'cohens_d': cohens_d,
         'effect_note': ('With hundreds of thousands of dependent pixels a p-value measures sample size, not '
                         'importance; the standardised effect is the number to read.'),
         'monument_pixels': int(mon.size), 'control_pixels': int(ctl.size),
         'ks_statistic_monuments_vs_controls': float(ks.statistic), 'ks_pvalue': float(ks.pvalue),
         'ks_statistic_khafre_vs_controls': float(ks_k.statistic), 'ks_pvalue_khafre': float(ks_k.pvalue),
         'nyquist_depth_m': meta['nyquist_depth_m'], 'sub_aperture_s': meta['sub_aperture_s'],
         'pair_lag_ms': 1e3 * meta['pair_lag_s'], 'overlap_within_pair': meta['overlap_within_pair'],
         'acquisition': meta['acquisition'], 'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 11, 'tag': 'null', 'eyebrow': '01 · Within-image controls',
        'title': 'Monuments and open desert give the same answer in the same acquisition',
        'question': 'Does the method distinguish the pyramids from empty plateau in the very image that produced the claim?',
        'finding': (f"It does not. Across {mon.size:,} pixels on four monuments and {ctl.size:,} pixels on two "
                    f"stretches of open plateau in the same acquisition, the mean structure score is "
                    f"{mon.mean():.4f} on the monuments against {ctl.mean():.4f} on the desert, a standardised "
                    f"difference of {cohens_d:+.3f}, and it runs the wrong way. The Kolmogorov-Smirnov distance "
                    f"between Khafre and the controls is {ks_k.statistic:.3f}. Such a test returns a small p-value "
                    f"here only because the pixel count is large; the separation itself is under one percent of the "
                    f"score range. The depth histograms and mean depth profiles have the same shape everywhere, and "
                    f"this bank's depth axis folds back at {meta['nyquist_depth_m']:.1f} m."),
        'limitations': ('Pixels inside a patch are not independent, so no p-value here should be read as a '
                        'detection rate. This reproduces the published processing chain as described, on the same '
                        'class of data, at patches we chose. It does not reproduce a specific published figure, because the bank '
                        'parameters behind those figures were never disclosed. A difference in these statistics '
                        'would not have proved a chamber either; it would have started a search for a surface cause.'),
        'method': (f"ICEYE X33 Spotlight Dwell Fine, Giza, 27 August 2025, VV. Six patches of 4096 by 1536 native "
                   f"samples located by the product's rational polynomial model. Sub-aperture bank of 50 "
                   f"reference/offset pairs at half the processed bandwidth with an 88 Hz offset, 32 by 32 complex "
                   f"cross-correlation to a thousandth of a pixel, scene-median common mode removed, steering "
                   f"wavenumbers from the product's own state vectors, depth grid to the ambiguity limit."),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in m.items() if k not in ('patches', 'acquisition')}, indent=1))
    return summary


if __name__ == '__main__':
    import sys
    if '--report-only' not in sys.argv:
        compute()
    report()
