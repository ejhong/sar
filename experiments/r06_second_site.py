"""R6: the same processing at a second megalithic site, on another satellite and continent.

Sacsayhuaman is a different sensor (ICEYE X35), a different pass geometry, a different
hemisphere and 3.6 km of Andean altitude. Its image contains cut-stone terrace walls, a bare
rock outcrop, a dense modern street grid and an open hillside. If the processing responds to
buried structure, these four should not behave alike; if it responds to its own design, they
will, in the same way the Giza patches did.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

from real_common import SACSAYHUAMAN_PATCHES, run_patch, patch_stats, RESULTS
from sarsim import viz

TEST = 'r06_second_site'
BANK = 'paper'


def compute(stride_az=12, stride_rg=8, force=False):
    for entry in SACSAYHUAMAN_PATCHES:
        print(f'  {entry[0]}', flush=True)
        run_patch('sacsayhuaman', entry[0], BANK, stride_az=stride_az, stride_rg=stride_rg,
                  force=force, verbose=True)


def report():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    runs = {e[0]: run_patch('sacsayhuaman', e[0], BANK, verbose=False) for e in SACSAYHUAMAN_PATCHES}
    figs = []

    fig, axs = plt.subplots(2, 4, figsize=(14.2, 7.0))
    for col, (key, run) in enumerate(runs.items()):
        g = (run['grid_rows'].size, run['grid_cols'].size)
        a = run['amplitude'].reshape(g)
        db = 20 * np.log10(np.maximum(a, 1e-3) / np.percentile(a, 99.5))
        axs[0, col].imshow(db.T, cmap=viz.CMAP_GRAY, vmin=-22, vmax=3, aspect='auto', interpolation='nearest')
        axs[0, col].set_title(run['meta']['label'], loc='left', fontsize=9.5)
        s = run['score'].max(axis=1).reshape(g)
        im = axs[1, col].imshow(s.T, cmap=viz.CMAP_SEQ, vmin=0, vmax=1, aspect='auto', interpolation='nearest')
        for ax in (axs[0, col], axs[1, col]):
            ax.set_xticks([]); ax.set_yticks([])
    axs[0, 0].set_ylabel('image amplitude', fontsize=9)
    axs[1, 0].set_ylabel('best depth-fit score', fontsize=9)
    fig.suptitle('Sacsayhuaman, ICEYE X35, 22 August 2025: four patches of one acquisition',
                 x=0.012, ha='left', fontsize=11.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.colorbar(im, ax=axs, fraction=0.016, pad=0.012, label='best adjusted $R^{2}$')
    figs.append({'file': viz.finish(fig, f'{fd}/r06_panels.png'), 'caption':
                 '<b>A second site, the same outcome.</b> Top: the native image amplitude, where the zigzag terrace '
                 'walls, the Rodadero outcrop, the Cusco street grid and an open hillside are all plainly '
                 'distinguishable. Bottom: the best depth-fit score the same processing assigns to each pixel. The '
                 'surface tells them apart at a glance; the depth output does not.'})

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2))
    bins = np.linspace(0, 1, 46)
    for i, (key, run) in enumerate(runs.items()):
        m = run['meta']
        ls = '-' if m['kind'] == 'monument' else ('-.' if m['kind'] == 'urban' else '--')
        axs[0].hist(run['score'].max(axis=1), bins=bins, density=True, histtype='step', lw=2,
                    color=viz.SERIES[i % 8], ls=ls, label=m['label'])
        t = run['tomo'].astype(float)
        tn = t / np.maximum(t.sum(axis=1, keepdims=True), 1e-30)
        axs[1].plot(run['z'], tn.mean(axis=0), lw=2, color=viz.SERIES[i % 8], ls=ls)
    axs[0].set_xlabel('best adjusted $R^{2}$'); axs[0].set_ylabel('density'); axs[0].legend(fontsize=8)
    axs[0].set_title('Structure score', loc='left')
    axs[1].set_xlabel('model depth (m)'); axs[1].set_ylabel('mean normalised $|h(z)|^{2}$')
    axs[1].set_title('Mean depth profile', loc='left')
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r06_distributions.png'), 'caption':
                 '<b>Cut stone, city and hillside overlap.</b> Score distributions and mean depth profiles for the '
                 'four Sacsayhuaman patches. A dense modern city and a bare Andean slope are about as different as '
                 'two radar scenes can be, and the depth output separates them no better than it separates either '
                 'from megalithic walls.'})

    stats = [patch_stats(r) for r in runs.values()]
    mon = np.concatenate([runs[e[0]]['score'].max(axis=1) for e in SACSAYHUAMAN_PATCHES if e[5] == 'monument'])
    ctl = runs['slope_north']['score'].max(axis=1)
    urb = runs['cusco']['score'].max(axis=1)
    ks = ks_2samp(mon, ctl)
    meta = runs['walls']['meta']
    m = {'patches': stats, 'monument_mean': float(mon.mean()), 'control_mean': float(ctl.mean()),
         'urban_mean': float(urb.mean()),
         'ks_statistic_walls_vs_slope': float(ks.statistic), 'ks_pvalue': float(ks.pvalue),
         'nyquist_depth_m': meta['nyquist_depth_m'], 'acquisition': meta['acquisition'],
         'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 16, 'tag': 'null', 'eyebrow': '06 · Second site',
        'title': 'Another satellite, another continent, the same undifferentiated output',
        'question': 'Does the result depend on Giza, on the ICEYE X33, or on the method?',
        'finding': (f"On the method. An ICEYE X35 dwell of Sacsayhuaman, a descending left-looking pass over the "
                    f"Andes at 3.6 km altitude, gives cut-stone terrace walls a mean structure score of "
                    f"{mon.mean():.3f}, an open hillside {ctl.mean():.3f} and the dense street grid of Cusco "
                    f"{urb.mean():.3f}. The Kolmogorov-Smirnov distance between walls and hillside is "
                    f"{ks.statistic:.3f}. The depth axis for this acquisition repeats every "
                    f"{meta['nyquist_depth_m']:.1f} m."),
        'limitations': ('Steep Andean terrain makes single-height RPC placement much less precise than at flat '
                        'Giza, so each patch was checked against the image rather than trusted from coordinates. '
                        'No claim is made about what lies under Sacsayhuaman; this is a control on the processing.'),
        'method': ('Identical pipeline, bank and grid to the Giza controls, applied to four patches of one ICEYE '
                   'X35 Spotlight Dwell Fine acquisition of 22 August 2025, with steering wavenumbers from that '
                   "product's own state vectors."),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in m.items() if k not in ('patches', 'acquisition')}, indent=1))


if __name__ == '__main__':
    import sys
    if '--report-only' not in sys.argv:
        compute()
    report()
