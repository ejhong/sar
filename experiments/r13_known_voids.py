"""R13: the positive control. Can the method see voids that are certainly there?

Giza's Eastern and Western Cemeteries are dense fields of mastabas whose burial shafts were
excavated and surveyed, typically 5 to 30 m deep. Hundreds of them sit within a few hundred
metres of open plateau that has none. That comparison needs placement only to about fifty
metres, which the product's own geolocation comfortably achieves, so it sidesteps the
metre-level registration problem that blocks a single-chamber test.

If this method has any sensitivity to shallow voids, the cemeteries must differ from the
plateau in the depth band where the shafts are. Two things are checked: whether the depth
profiles differ at all, and whether any difference is specific to the shaft band or is simply
the cemeteries being brighter at the surface.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

from real_common import run_patch, RESULTS
from sarsim import viz

TEST = 'r13_known_voids'
TOMBS = [('cemetery_east', 'Eastern Cemetery'), ('cemetery_west', 'Western Cemetery')]
CONTROLS = [('desert_west', 'Open plateau, west'), ('plateau_nw', 'Open plateau, north-west')]
SHAFT_BAND_M = (5.0, 30.0)


def profiles(run):
    t = run['tomo'].astype(float)
    return t / np.maximum(t.sum(axis=1, keepdims=True), 1e-30)


def band_share(run, band=SHAFT_BAND_M):
    """Fraction of each pixel's depth profile that falls in the shaft-depth band."""
    z = run['z']
    sel = (z >= band[0]) & (z <= band[1])
    return profiles(run)[:, sel].sum(axis=1)


def main(stride_az=12, stride_rg=8):
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    runs = {}
    for key, _ in TOMBS + CONTROLS:
        print(f'  {key}', flush=True)
        runs[key] = run_patch('giza', key, 'paper', stride_az=stride_az, stride_rg=stride_rg, verbose=True)

    tomb_share = np.concatenate([band_share(runs[k]) for k, _ in TOMBS])
    ctrl_share = np.concatenate([band_share(runs[k]) for k, _ in CONTROLS])
    tomb_amp = np.concatenate([runs[k]['amplitude'].astype(float) for k, _ in TOMBS])
    ctrl_amp = np.concatenate([runs[k]['amplitude'].astype(float) for k, _ in CONTROLS])
    ks = ks_2samp(tomb_share, ctrl_share)
    pooled = np.sqrt((tomb_share.var() * tomb_share.size + ctrl_share.var() * ctrl_share.size) /
                     (tomb_share.size + ctrl_share.size))
    d_raw = float((tomb_share.mean() - ctrl_share.mean()) / pooled)

    # brightness-matched comparison: the cemeteries are far brighter at the surface, so any
    # difference has to survive matching on amplitude before it can mean anything about depth
    edges = np.percentile(np.concatenate([tomb_amp, ctrl_amp]), np.linspace(2, 98, 25))
    matched_t, matched_c = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        a = tomb_share[(tomb_amp >= lo) & (tomb_amp < hi)]
        b = ctrl_share[(ctrl_amp >= lo) & (ctrl_amp < hi)]
        n = min(a.size, b.size)
        if n >= 40:
            matched_t.append(a[:n]); matched_c.append(b[:n])
    mt = np.concatenate(matched_t) if matched_t else np.array([])
    mc = np.concatenate(matched_c) if matched_c else np.array([])
    d_matched = float((mt.mean() - mc.mean()) /
                      np.sqrt((mt.var() * mt.size + mc.var() * mc.size) / (mt.size + mc.size))) if mt.size else float('nan')

    # --- the patches themselves: surface, then the shaft-band map, on one shared scale
    order = TOMBS + CONTROLS
    shares = {k: band_share(runs[k]) for k, _ in order}
    lo_s, hi_s = np.percentile(np.concatenate(list(shares.values())), [3, 97])
    fig, axs = plt.subplots(2, 4, figsize=(14.2, 7.2))
    im = None
    for c, (key, label) in enumerate(order):
        run = runs[key]
        g = (run['grid_rows'].size, run['grid_cols'].size)
        a = run['amplitude'].astype(float).reshape(g)
        db = 20 * np.log10(np.maximum(a, 1e-3) / np.percentile(a, 99.5))
        axs[0, c].imshow(db.T, cmap=viz.CMAP_GRAY, vmin=-22, vmax=3, aspect='auto', interpolation='nearest')
        axs[0, c].set_title(label, loc='left', fontsize=9.5)
        im = axs[1, c].imshow(shares[key].reshape(g).T, cmap=viz.CMAP_SEQ, vmin=lo_s, vmax=hi_s,
                              aspect='auto', interpolation='nearest')
        for ax in (axs[0, c], axs[1, c]):
            ax.set_xticks([]); ax.set_yticks([])
        axs[0, c].text(0.985, 0.04, run['meta']['kind'], transform=axs[0, c].transAxes, ha='right',
                       fontsize=7.5, color='#e8e4d8', family='monospace')
    axs[0, 0].set_ylabel('surface amplitude', fontsize=9)
    axs[1, 0].set_ylabel(f'share in {SHAFT_BAND_M[0]:.0f}–{SHAFT_BAND_M[1]:.0f} m band', fontsize=9)
    fig.suptitle('Two cemeteries of surveyed burial shafts, two patches of bare plateau',
                 x=0.012, ha='left', fontsize=11.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 0.955, 0.955])
    cax = fig.add_axes([0.965, 0.1, 0.011, 0.36])
    fig.colorbar(im, cax=cax, label='shaft-band share')
    maps_fig = {'file': viz.finish(fig, f'{fd}/r13_maps.png'), 'caption':
                '<b>The mastaba rows are obvious; the shafts beneath them are not.</b> Top: surface amplitude of '
                'the four patches, where the rectangular mastaba blocks of the Eastern and Western Cemeteries are '
                'unmistakable against featureless plateau. Bottom: for the same pixels, how much of each depth '
                'profile falls in the 5 to 30 m band where the surveyed burial shafts actually are, on one shared '
                'scale. The structures on the surface print clearly. The hundreds of voids beneath them do not.'}

    fig, axs = plt.subplots(1, 3, figsize=(13.6, 4.3))
    z = runs[TOMBS[0][0]]['z']
    for i, (key, label) in enumerate(TOMBS + CONTROLS):
        style = dict(lw=2, color=viz.SERIES[i % 8], ls='-' if i < len(TOMBS) else '--')
        axs[0].plot(z, profiles(runs[key]).mean(axis=0), label=label, **style)
        axs[1].hist(band_share(runs[key]), bins=np.linspace(0, 0.6, 50), density=True,
                    histtype='step', **style)
    axs[0].axvspan(*SHAFT_BAND_M, color=viz.GRID, alpha=0.75, lw=0)
    axs[0].text(np.mean(SHAFT_BAND_M), 0.97, 'surveyed shaft depths', transform=axs[0].get_xaxis_transform(),
                ha='center', va='top', fontsize=8.5, color=viz.INK2)
    axs[0].set_xlabel('model depth (m)'); axs[0].set_ylabel('mean normalised $|h(z)|^{2}$')
    axs[0].set_title('Depth profiles', loc='left'); axs[0].legend(fontsize=7.5)
    axs[1].set_xlabel('share of the profile in the shaft band'); axs[1].set_ylabel('density')
    axs[1].set_title('Shaft-band share per pixel', loc='left')
    axs[2].bar([0, 1], [d_raw, d_matched], color=[viz.SERIES[1], viz.SERIES[0]], width=0.5)
    axs[2].axhline(0, color=viz.AXIS, lw=1)
    for x, v in ((0, d_raw), (1, d_matched)):
        axs[2].text(x, v + (0.01 if v >= 0 else -0.03), f'{v:+.3f}', ha='center', fontsize=9.5, color=viz.INK2)
    axs[2].set_xticks([0, 1]); axs[2].set_xticklabels(['as measured', 'brightness-matched'], fontsize=9)
    axs[2].set_ylabel('standardised difference, tombs minus plateau')
    axs[2].set_ylim(min(-0.25, d_raw * 1.4, d_matched * 1.4), max(0.25, d_raw * 1.4, d_matched * 1.4))
    axs[2].set_title('Effect size in the shaft band', loc='left')
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r13_cemeteries.png'), 'caption':
             '<b>Hundreds of real voids, and nothing to see.</b> The Eastern and Western Cemeteries are dense '
             'fields of excavated burial shafts 5 to 30 m deep; the two control patches are bare plateau a few '
             'hundred metres away. Left: mean depth profiles, with the surveyed shaft band shaded. Middle: the '
             'share of each pixel\'s profile falling in that band. Right: the standardised difference between tombs '
             'and plateau, before and after matching on surface brightness. If the method could see shallow voids, '
             'this is where it would show.'}]
    figs.insert(0, maps_fig)

    m = {'shaft_band_m': list(SHAFT_BAND_M),
         'tomb_pixels': int(tomb_share.size), 'control_pixels': int(ctrl_share.size),
         'tomb_band_share_mean': float(tomb_share.mean()), 'control_band_share_mean': float(ctrl_share.mean()),
         'standardised_difference': d_raw, 'standardised_difference_brightness_matched': d_matched,
         'matched_pixels': int(mt.size), 'ks_statistic': float(ks.statistic), 'ks_pvalue': float(ks.pvalue),
         'tomb_median_amplitude': float(np.median(tomb_amp)),
         'control_median_amplitude': float(np.median(ctrl_amp)),
         'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 24, 'tag': 'null', 'eyebrow': '13 · Positive control',
        'title': 'Hundreds of surveyed burial shafts produce no signal at their known depths',
        'question': 'Can the method find voids that are certainly there, at depths that are certainly known?',
        'finding': (f"No. Across {tomb_share.size:,} pixels over the Eastern and Western Cemeteries and "
                    f"{ctrl_share.size:,} pixels of bare plateau, the share of each pixel's depth profile falling "
                    f"in the surveyed shaft band of 5 to 30 m is {tomb_share.mean():.4f} over the tombs against "
                    f"{ctrl_share.mean():.4f} over the plateau, a standardised difference of {d_raw:+.3f}. The "
                    f"cemeteries are far brighter at the surface, with median amplitude "
                    f"{np.median(tomb_amp):.0f} against {np.median(ctrl_amp):.0f}, so the comparison was repeated "
                    f"on brightness-matched pixels, where the difference is {d_matched:+.3f}. This is the easiest "
                    f"target the site offers: shallow, numerous, excavated and surveyed."),
        'limitations': ('Shaft depths vary between tombs and many shafts are filled or collapsed, so the band is a '
                        'range rather than a target depth. Placement is verified by recognising the mastaba rows '
                        'in the image, not by surveyed control points, which is enough at this scale but would not '
                        'be for a single chamber. A null here bounds sensitivity to shallow voids under these '
                        'conditions; it does not prove no void anywhere could ever be detected.'),
        'method': ('Identical pipeline, bank and grid to the within-image controls, applied to two cemetery patches '
                   'and two plateau patches of the same acquisition. The statistic is the fraction of each pixel\'s '
                   'normalised depth profile lying in the 5 to 30 m band, compared directly and again after '
                   'matching the two populations on surface amplitude in 24 bins.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in m.items()}, indent=1))


if __name__ == '__main__':
    main()
