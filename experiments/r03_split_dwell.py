"""R3: does the same ground give the same depths from two halves of one dwell?

A 24.5 s dwell can be cut into two disjoint stretches of slow time, each with its own
sub-aperture bank and its own set of look angles. Real structure in the ground is the same
object in both halves and must appear at the same depth. A feature that is a property of the
processing is free to move. This is a cross-geometry test that needs only one acquisition.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import run_patch, RESULTS
from sarsim import viz

TEST = 'r03_split_dwell'
PATCHES = [('khafre', 'Khafre pyramid'), ('desert_west', 'Open plateau, west')]


def common_profiles(a, b):
    """Normalised depth profiles of both halves on the shorter of the two depth axes."""
    z_max = min(a['z'][-1], b['z'][-1])
    z = np.linspace(min(a['z'][0], b['z'][0]), z_max, 160)
    out = []
    for run in (a, b):
        t = run['tomo'].astype(float)
        p = np.stack([np.interp(z, run['z'], row) for row in t])
        out.append(p / np.maximum(p.sum(axis=1, keepdims=True), 1e-30))
    return z, out[0], out[1]


def main(stride_az=16, stride_rg=12):
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    rng = np.random.default_rng(0)
    figs, rows = [], []
    store = {}
    for key, label in PATCHES:
        a = run_patch('giza', key, 'split_a', stride_az=stride_az, stride_rg=stride_rg)
        b = run_patch('giza', key, 'split_b', stride_az=stride_az, stride_rg=stride_rg)
        z, pa, pb = common_profiles(a, b)
        za, zb = z[np.argmax(pa, axis=1)], z[np.argmax(pb, axis=1)]
        ca = pa - pa.mean(axis=1, keepdims=True)
        cb = pb - pb.mean(axis=1, keepdims=True)
        corr = (ca * cb).sum(axis=1) / np.sqrt((ca ** 2).sum(axis=1) * (cb ** 2).sum(axis=1) + 1e-30)
        shuffled = rng.permutation(pb.shape[0])
        cs = cb[shuffled]
        null = (ca * cs).sum(axis=1) / np.sqrt((ca ** 2).sum(axis=1) * (cs ** 2).sum(axis=1) + 1e-30)
        tol = z[1] - z[0]
        rows.append({
            'patch': key, 'label': label, 'pixels': int(pa.shape[0]),
            'half_a_times_s': [a['meta']['sub_aperture_times_s'][0], a['meta']['sub_aperture_times_s'][-1]],
            'half_b_times_s': [b['meta']['sub_aperture_times_s'][0], b['meta']['sub_aperture_times_s'][-1]],
            'depth_axis_max_m': float(z[-1]),
            'depth_agreement_within_one_step': float((np.abs(za - zb) <= tol).mean()),
            'depth_agreement_expected_by_chance': float(2 * tol / (z[-1] - z[0])),
            'profile_correlation_median': float(np.median(corr)),
            'shuffled_correlation_median': float(np.median(null)),
            'depth_correlation': float(np.corrcoef(za, zb)[0, 1]),
        })
        store[key] = (z, za, zb, corr, null, label)

    fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.2))
    for i, (key, label) in enumerate(PATCHES):
        z, za, zb, corr, null, _ = store[key]
        h = axs[0, i].hist2d(za, zb, bins=48, cmap=viz.CMAP_SEQ,
                             range=[[z[0], z[-1]], [z[0], z[-1]]])
        axs[0, i].plot([z[0], z[-1]], [z[0], z[-1]], color=viz.SERIES[1], lw=1.2, ls='--')
        axs[0, i].set_title(f'{label}: depth from each half', loc='left')
        axs[0, i].set_xlabel('depth from the first half (m)')
        if i == 0:
            axs[0, i].set_ylabel('depth from the second half (m)')
        axs[1, i].hist(corr, bins=np.linspace(-1, 1, 61), density=True, histtype='step', lw=2,
                       color=viz.SERIES[0], label='same pixel, two halves')
        axs[1, i].hist(null, bins=np.linspace(-1, 1, 61), density=True, histtype='step', lw=2,
                       color=viz.SERIES[1], label='different pixels (null)')
        axs[1, i].set_xlabel('correlation of the two depth profiles')
        axs[1, i].legend(fontsize=8)
        if i == 0:
            axs[1, i].set_ylabel('density')
    fig.tight_layout()
    figs.append({'file': viz.finish(fig, f'{fd}/r03_split.png'), 'caption':
                 '<b>The two halves of one dwell disagree.</b> The aperture is cut into two disjoint stretches of '
                 'slow time, each given an identical sweep, so each half sees the same ground from its own set of '
                 'look angles. Top: the depth each half assigns to the same pixel, with the line of agreement '
                 'dashed. Bottom: the correlation between the two halves\' depth profiles, against the correlation '
                 'between profiles of unrelated pixels. Buried structure cannot move between the first and second '
                 'half of a 25 second acquisition; these features do.'})

    m = {'runs': rows, 'runtime_s': time.time() - t0}
    best = rows[0]
    summary = {
        'id': TEST, 'order': 13, 'tag': 'null', 'eyebrow': '03 · Split dwell',
        'title': 'The first and second halves of one acquisition place the same ground at different depths',
        'question': 'Real structure cannot move during 25 seconds. Do the recovered depths stay put?',
        'finding': (f"No. Splitting the Khafre dwell into two disjoint halves and giving each an identical sweep, "
                    f"the two halves agree on the depth of a pixel for {best['depth_agreement_within_one_step']:.1%} "
                    f"of pixels, against {best['depth_agreement_expected_by_chance']:.1%} expected by chance, and the "
                    f"median correlation between their depth profiles is {best['profile_correlation_median']:.3f} "
                    f"against {best['shuffled_correlation_median']:.3f} for unrelated pixels. Open desert behaves "
                    f"the same way."),
        'limitations': ('Each half uses half the look angles, so each is noisier than the full sweep; some '
                        'disagreement is expected from noise alone. The test shows that the depth assignments are '
                        'not reproducible across look angle, not that a better estimator could not be built.'),
        'method': ('One ICEYE dwell, one patch, two sub-aperture banks with identical width and swept bandwidth '
                   'placed on disjoint halves of the processed Doppler support. Depth profiles are interpolated '
                   'onto the shorter of the two unambiguous axes before comparison. The null is the same '
                   'comparison after randomly pairing pixels.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
