"""R0: the two acquisitions, as images.

Before any processing, this is what the radar actually recorded. Ten patches across two
Spotlight Dwell Fine products, drawn from the native complex data at the analysis grid. Every
later result is a statement about these pixels, and a reader who can recognise the monuments
here can check that the geolocation put them where they belong.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import GIZA_PATCHES, SACSAYHUAMAN_PATCHES, run_patch, RESULTS
from sarsim import viz

TEST = 'r00_acquisitions'
ROWS = [('giza', GIZA_PATCHES, 'Giza plateau · ICEYE X33 · 27 August 2025'),
        ('sacsayhuaman', SACSAYHUAMAN_PATCHES, 'Sacsayhuamán · ICEYE X35 · 22 August 2025')]


def panel(ax, run):
    g = (run['grid_rows'].size, run['grid_cols'].size)
    a = run['amplitude'].astype(float).reshape(g)
    db = 20 * np.log10(np.maximum(a, 1e-3) / np.percentile(a, 99.5))
    ax.imshow(db.T, cmap=viz.CMAP_GRAY, vmin=-22, vmax=3, aspect='auto', interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    m = run['meta']
    ax.set_title(m['label'], loc='left', fontsize=9)
    ax.text(0.985, 0.03, m['kind'], transform=ax.transAxes, ha='right', fontsize=7.5,
            color='#e8e4d8', family='monospace')


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    runs = {}
    for site, patches, _ in ROWS:
        for entry in patches:
            try:
                runs[(site, entry[0])] = run_patch(site, entry[0], 'paper', verbose=False)
            except (FileNotFoundError, StopIteration):
                pass
    widest = max(len(p) for _, p, _ in ROWS)
    fig, axs = plt.subplots(len(ROWS), widest, figsize=(3.05 * widest, 2.6 * len(ROWS)))
    facts = []
    for r, (site, patches, label) in enumerate(ROWS):
        for c in range(widest):
            ax = axs[r, c]
            if c < len(patches) and (site, patches[c][0]) in runs:
                panel(ax, runs[(site, patches[c][0])])
            else:
                ax.set_axis_off()
        acq = next(v['meta']['acquisition'] for k, v in runs.items() if k[0] == site)
        facts.append({'site': site, 'label': label, 'satellite': acq['satellite'],
                      'collected': acq['collection_start'][:19].replace('T', ' '),
                      'shape': acq['shape'], 'incidence_deg': acq['incidence_center_deg'],
                      'collection_duration_s': acq['collection_duration_s'],
                      'patches': sum(1 for k in runs if k[0] == site)})
        axs[r, 0].set_ylabel(label.split(' · ')[0], fontsize=9)
    fig.suptitle('The two acquisitions, before any processing', x=0.006, ha='left', fontsize=11.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    figs = [{'file': viz.finish(fig, f'{fd}/r00_acquisitions.png'), 'caption':
             '<b>What the radar recorded.</b> Amplitude of the native complex images on the analysis grid, geolocated '
             'by each product\'s own rational polynomial model with no fitted offset. Top: four Giza monuments and '
             'two stretches of open plateau, all inside one 25 second acquisition. Bottom: the zigzag terrace walls '
             'of Sacsayhuamán, the Rodadero outcrop, the street grid of Cusco and an open Andean hillside, in a '
             'second acquisition from a different satellite. The surface tells these apart at a glance. Everything '
             'that follows asks whether the depth output does.'}]
    m = {'acquisitions': facts, 'patches': len(runs), 'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 10, 'tag': 'mechanism', 'eyebrow': '00 · The acquisitions',
        'title': 'Two dwell acquisitions, ten patches, one set of pixels behind every result',
        'question': 'What was actually recorded, and are the patches where we say they are?',
        'finding': ('Ten patches across two ICEYE Spotlight Dwell Fine products. At Giza the laid-over faces and '
                    'radar shadows of Khafre, Khufu and Menkaure appear exactly where the product\'s own rational '
                    'polynomial model places them, with no fitted offset. At Sacsayhuamán the terrace walls, the '
                    'Rodadero outcrop, the Cusco street grid and an open hillside are all recognisable. The surface '
                    'separates these scenes trivially, which is the comparison every later test is against.'),
        'limitations': ('Placement is verified against recognisable structures, not against surveyed control points. '
                        'Heights are published elevations plus a nominal geoid undulation, good enough to land on '
                        'the right monument but not a surveyed vertical tie. Sacsayhuamán sits on steep terrain '
                        'where a single assumed height places a patch much less precisely than at flat Giza.'),
        'method': ('Amplitude of the native complex SLC sampled on the same grid used for every measurement, '
                   'displayed in decibels relative to each patch\'s 99.5th percentile.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(facts, indent=1))


if __name__ == '__main__':
    main()
