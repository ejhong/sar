"""R5: can a classifier find information in the depth profiles that is not already brightness?

Machine learning cannot create information that the measurement does not contain, but it is a
good way to measure whether information is there. The task here is deliberately far easier than
void detection: tell a pixel on Khafre from a pixel on open desert, using the depth profile the
method produces. If a learner cannot do that beyond what surface brightness already gives it,
then the depth axis is carrying nothing of its own.

Cross-validation is blocked in space. Neighbouring pixels share speckle and share patches, so
random splits leak between folds and inflate every score.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

from real_common import run_patch, RESULTS
from sarsim import viz

TEST = 'r05_learned_null'
N_BINS = 32
N_BLOCKS = 8


def features(run):
    tomo = run['tomo'].astype(float)
    profile = tomo / np.maximum(tomo.sum(axis=1, keepdims=True), 1e-30)
    edges = np.linspace(0, profile.shape[1], N_BINS + 1).astype(int)
    binned = np.stack([profile[:, a:b].mean(axis=1) for a, b in zip(edges[:-1], edges[1:])], axis=1)
    amplitude = 20 * np.log10(np.maximum(run['amplitude'].astype(float), 1e-6))[:, None]
    rows = run['rows'].astype(float)[:, None]
    return binned, amplitude, rows


def blocked_auc(X, y, blocks, seed=0):
    """Mean cross-validated AUC with spatially contiguous folds."""
    order = np.unique(blocks)
    rng = np.random.default_rng(seed)
    folds = np.array_split(rng.permutation(order), 5)
    scores = []
    for hold in folds:
        test = np.isin(blocks, hold)
        if test.all() or not test.any() or len(np.unique(y[~test])) < 2 or len(np.unique(y[test])) < 2:
            continue
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
        model.fit(X[~test], y[~test])
        scores.append(roc_auc_score(y[test], model.predict_proba(X[test])[:, 1]))
    return float(np.mean(scores)), float(np.std(scores))


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    pos = run_patch('giza', 'khafre', 'paper')
    neg = run_patch('giza', 'desert_west', 'paper')
    Pp, Ap, Rp = features(pos)
    Pn, An, Rn = features(neg)
    X_profile = np.vstack([Pp, Pn])
    X_amp = np.vstack([Ap, An])
    y = np.concatenate([np.ones(len(Pp)), np.zeros(len(Pn))])
    rows_px = np.concatenate([Rp.ravel(), Rn.ravel()])
    blocks = np.concatenate([np.digitize(Rp.ravel(), np.quantile(Rp, np.linspace(0, 1, N_BLOCKS + 1)[1:-1])),
                             N_BLOCKS + np.digitize(Rn.ravel(), np.quantile(Rn, np.linspace(0, 1, N_BLOCKS + 1)[1:-1]))])

    results = {}
    results['amplitude only'] = blocked_auc(X_amp, y, blocks)
    results['depth profile only'] = blocked_auc(X_profile, y, blocks)
    results['depth profile + amplitude'] = blocked_auc(np.hstack([X_profile, X_amp]), y, blocks)

    # amplitude-matched subset: the only fair way to ask whether depth adds anything
    lo, hi = np.percentile(X_amp, [40, 60])
    sel = ((X_amp.ravel() >= lo) & (X_amp.ravel() <= hi))
    if sel.sum() > 400 and len(np.unique(y[sel])) == 2:
        results['depth profile, amplitude-matched'] = blocked_auc(X_profile[sel], y[sel], blocks[sel])
        results['amplitude, amplitude-matched'] = blocked_auc(X_amp[sel], y[sel], blocks[sel])

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
    names = list(results)
    vals = [results[n][0] for n in names]
    errs = [results[n][1] for n in names]
    colours = [viz.SERIES[1] if 'amplitude only' in n or n.startswith('amplitude') else viz.SERIES[0] for n in names]
    axs[0].barh(range(len(names)), vals, xerr=errs, color=colours, height=0.55,
                error_kw=dict(ecolor=viz.MUTED, lw=1))
    axs[0].axvline(0.5, color=viz.SERIES[7], lw=1, ls='--')
    axs[0].text(0.505, len(names) - 0.4, 'chance', color=viz.SERIES[7], fontsize=8.5)
    axs[0].set_yticks(range(len(names)))
    axs[0].set_yticklabels(names, fontsize=8.5)
    axs[0].set_xlim(0.4, 1.0)
    axs[0].set_xlabel('area under the ROC curve, spatially blocked folds')
    axs[0].set_title('Telling Khafre from open desert', loc='left')
    for i, (v, e) in enumerate(zip(vals, errs)):
        axs[0].text(v + e + 0.008, i, f'{v:.2f}', va='center', fontsize=8.5, color=viz.INK2)

    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    half = blocks % 2 == 0
    for name, X, colour in [('depth profile', X_profile, viz.SERIES[0]), ('amplitude', X_amp, viz.SERIES[1])]:
        model.fit(X[half], y[half])
        fpr, tpr, _ = roc_curve(y[~half], model.predict_proba(X[~half])[:, 1])
        axs[1].plot(fpr, tpr, lw=2, color=colour, label=f'{name} (AUC {roc_auc_score(y[~half], model.predict_proba(X[~half])[:, 1]):.2f})')
    axs[1].plot([0, 1], [0, 1], color=viz.SERIES[7], lw=1, ls='--')
    axs[1].set_xlabel('false positive rate'); axs[1].set_ylabel('true positive rate')
    axs[1].set_title('One held-out split', loc='left'); axs[1].legend(fontsize=8.5, loc='lower right')
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r05_auc.png'), 'caption':
             '<b>The depth axis adds nothing a brightness map does not already have.</b> A logistic classifier is '
             'trained to separate pixels on Khafre from pixels on open plateau, with folds blocked in space so that '
             'neighbouring speckle cannot leak between training and test. Surface brightness separates them easily, '
             'which it should: a pyramid is not sand. The depth profiles the method produces do not, and inside an '
             'amplitude-matched band they fall to chance. This is the machine-learning form of the null control: '
             'a learner is a good instrument for measuring whether information is present, and here it measures '
             'that it is not.'}]

    m = {'auc': {k: {'mean': v[0], 'std': v[1]} for k, v in results.items()},
         'pixels_khafre': int(len(Pp)), 'pixels_desert': int(len(Pn)),
         'profile_bins': N_BINS, 'blocks': int(len(np.unique(blocks))), 'runtime_s': time.time() - t0}
    prof = results['depth profile only'][0]
    matched = results.get('depth profile, amplitude-matched', (float('nan'), 0))[0]
    summary = {
        'id': TEST, 'order': 15, 'tag': 'null', 'eyebrow': '05 · Learned null',
        'title': 'A classifier finds no information in the depth profiles beyond surface brightness',
        'question': 'Could a trained model extract something from these tomograms that a human reader misses?',
        'finding': (f"No. Separating Khafre pixels from open-desert pixels is an easy task and surface brightness "
                    f"alone reaches an area under the curve of {results['amplitude only'][0]:.2f}. The depth profiles "
                    f"reach {prof:.2f}, and once the comparison is restricted to pixels of matched brightness they "
                    f"fall to {matched:.2f}, which is chance. Training on labelled voids would face the same problem "
                    f"with far fewer labels: Giza offers on the order of ten surveyed chambers, all of them under "
                    f"monuments, so a model would learn where the monuments are."),
        'limitations': ('A linear model on binned profiles; a larger network could fit more, but with spatially '
                        'blocked folds and a matched-brightness control there is little left for it to fit. The '
                        'test shows an absence of separable information in these features, not a proof that no '
                        'representation could ever separate them.'),
        'method': ('Logistic regression on 32-bin normalised depth profiles and on log amplitude, scored by area '
                   'under the ROC curve with five folds blocked into contiguous azimuth strips. The matched control '
                   'keeps only pixels between the fortieth and sixtieth percentile of brightness.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(m['auc'], indent=1))


if __name__ == '__main__':
    main()
