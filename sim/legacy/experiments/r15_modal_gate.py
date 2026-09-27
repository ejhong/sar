"""R15: the modal ("spring") branch, and whether it can break the depth ambiguity.

The derivative protocol has a second branch often described as a spring or resonance model.
Section 10 of its mathematical development fits, inside every contiguous 25-pair window, the
ordered elliptical harmonic model

    q_n = c + u cos(2*pi*m*n/W) + v sin(2*pi*m*n/W),   m = 1..6,

and admits a window to the depth inversion only when the fitted ellipse passes a quality gate.
The hope is that this extra structure supplies information the plain depth fit lacks, and so
resolves the repetition of the depth axis.

Four things are tested, using the protocol's own modal.py and focus.py unmodified, fed with real
sub-pixel registration vectors measured from the Giza acquisition:

A. Whether the repetition is a property of the depth model rather than of the data.
B. What the frozen gate accepts on real ground, on open desert, and on structureless input.
C. Whether gating changes the depth axis.
D. What sets the depth scale.

A mechanical resonance is a different idea from this modal branch, and is checked separately
against the frequency band a single dwell can carry.
"""
import hashlib
import json
import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import RESULTS, run_patch
from sarsim import viz

TEST = 'r15_modal_gate'
V17 = os.environ.get('BIONDI_V17', os.path.expanduser('~/tmp/sar/biondi_v17'))
# sha256 of the two stages as published, so a changed package cannot silently alter this result
EXPECTED = {'modal.py': '7735bb58ca3d7b7d34a0b0983f56eae23fc93acd46f359d5f9feade2dc27d20f',
            'focus.py': '58dbce38595bba41021badc187d86f20294a72dd5a2b147dd3f728a19edb735f'}
# the frozen Branch B gate, from section 12 of the mathematical development
FROZEN = {'adj_r2': 0.25, 'axis_ratio': 0.10, 'minor_px': 0.005}
SHAPE = {'adj_r2': 0.25, 'axis_ratio': 0.10, 'minor_px': 0.0}   # same shape rule, no absolute size
W = 25
MODES = (1, 6)
N_POS_MODAL = 600
N_POS_FOCUS = 200


def load_protocol():
    src = os.path.join(V17, 'src')
    if not os.path.isdir(src):
        raise SystemExit(f'set BIONDI_V17 to the unpacked v1.7 package (looked in {V17})')
    for name, want in EXPECTED.items():
        got = hashlib.sha256(open(os.path.join(src, 'biondi_protocol', name), 'rb').read()).hexdigest()
        if got != want:
            raise SystemExit(f'{name} sha256 {got} != published {want}')
    sys.path.insert(0, src)
    from biondi_protocol import modal, focus
    return modal, focus


def adjusted(r2, w=W, params=6):
    """Adjusted R-squared for the joint two-component fit: 2W observations, 6 free parameters."""
    n = 2 * w
    return 1.0 - (1.0 - r2) * (n - 1) / (n - params)


def modal_run(modal, q):
    """The published modal fit, reduced to one best-of-six record per window."""
    spec = {'window_length_w': W, 'mode_min': MODES[0], 'mode_max': MODES[1]}
    tab = modal.fit_modes(np.asarray(q, float), spec)['table']
    adj = adjusted(tab['r2'])
    key = tab['position'].astype(np.int64) * 10_000 + tab['window_start']
    uniq, inv = np.unique(key, return_inverse=True)
    best_adj = np.full(uniq.size, -np.inf)
    np.maximum.at(best_adj, inv, adj)
    return tab, adj, inv, uniq.size, best_adj


def apply_gate(tab, adj, inv, n_win, gate):
    ok = (adj >= gate['adj_r2']) & (tab['axis_ratio'] >= gate['axis_ratio']) & (tab['minor'] >= gate['minor_px'])
    win_ok = np.zeros(n_win, bool)
    np.logical_or.at(win_ok, inv, ok)
    return win_ok


def auc(a, b):
    """Probability a random value from a exceeds one from b; 0.5 means no discrimination."""
    x = np.concatenate([a, b])
    r = np.argsort(np.argsort(x)) + 1.0
    n1 = a.size
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * b.size))


def main():
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST)
    os.makedirs(fd, exist_ok=True)
    viz.style()
    modal, focus = load_protocol()

    # the protocol's own frozen bank design: K=50, half-band masks, 88 Hz separation
    khafre = run_patch('giza', 'khafre', 'paper')
    desert = run_patch('giza', 'desert_west', 'paper')
    kz = np.asarray(khafre['kz'], float)
    spacings = np.diff(np.sort(kz))
    dkz = float(np.median(spacings))
    lattice_rec = 2 * np.pi / abs(dkz)
    lam_s = float(khafre['meta']['lam_s_m'])
    uniformity = float(spacings.std() / abs(spacings.mean()))
    print(f'  kz spacing {dkz:.5f} rad/m (spread {uniformity*100:.1f}%) -> lattice repeat '
          f'{lattice_rec:.2f} m, unique interval {lattice_rec/2:.2f} m', flush=True)

    rng = np.random.default_rng(20250923)
    qk, qd = np.asarray(khafre['q'], float), np.asarray(desert['q'], float)
    pick_k = rng.choice(qk.shape[0], min(N_POS_MODAL, qk.shape[0]), replace=False)
    pick_d = rng.choice(qd.shape[0], min(N_POS_MODAL, qd.shape[0]), replace=False)
    qk_s, qd_s = qk[pick_k], qd[pick_d]
    scale = float(np.std(qk_s))
    q_noise = rng.normal(0.0, scale, qk_s.shape)                     # white, no structure at all
    walk = np.cumsum(rng.normal(0.0, scale, qk_s.shape), axis=1)     # smooth, but still meaningless
    q_walk = walk * (scale / float(np.std(walk)))                    # matched in spread to the real data

    # ---- A: is the repetition a property of the model rather than the data? ----
    # design(z + Z) = design(z) @ R for a fixed 2x2 R when the wavenumbers lie on a lattice, so the
    # column space, and therefore every least-squares residual, is unchanged.
    d0 = np.column_stack([np.cos(kz * 11.0), np.sin(kz * 11.0)])
    d1 = np.column_stack([np.cos(kz * (11.0 + lattice_rec)), np.sin(kz * (11.0 + lattice_rec))])
    R, *_ = np.linalg.lstsq(d0, d1, rcond=None)
    factor_resid = float(np.abs(d1 - d0 @ R).max() / np.abs(d1).max())

    step = lattice_rec / 200
    fk = focus.focus_depth(qk_s[:N_POS_FOCUS], kz, {'z_min_m': 0.0, 'z_max_m': 3 * lattice_rec,
                                                    'z_step_m': step})
    z = np.asarray(fk['z_m'], float)
    sc = np.asarray(fk['score'], float)

    def at_lag(mat, lag):
        n = mat.shape[1] - lag
        a, b = mat[:, :n].ravel(), mat[:, lag:lag + n].ravel()
        return float(np.corrcoef(a, b)[0, 1]), float(np.abs(a - b).max())

    # the protocol says to measure recurrence numerically when the ladder is uneven; do that too
    lags = np.arange(int(0.5 * lattice_rec / step), int(1.6 * lattice_rec / step))
    corrs = np.array([at_lag(sc, int(L))[0] for L in lags])
    best_lag = int(lags[int(np.argmax(corrs))])
    meas_rec = best_lag * step
    rep_corr, rep_max = at_lag(sc, int(round(lattice_rec / step)))
    best_corr, best_max = at_lag(sc, best_lag)
    half_corr, _ = at_lag(sc, int(round(lattice_rec / step)) // 2)
    print(f'  A: column-space residual {factor_resid:.2e}; at the lattice repeat corr {rep_corr:.6f} '
          f'(max diff {rep_max:.1e}); best numerical repeat {meas_rec:.2f} m corr {best_corr:.6f}; '
          f'half-period corr {half_corr:.3f}', flush=True)

    # ---- B: what does the gate accept? ----
    runs = {}
    cases = []
    for q_, lab in ((qk_s, 'Khafre pyramid'), (qd_s, 'Open plateau'),
                    (q_noise, 'White noise'), (q_walk, 'Random walk')):
        tab, adj, inv, n_win, best_adj = modal_run(modal, q_)
        runs[lab] = best_adj
        row = {'label': lab, 'windows': int(n_win),
               'frozen_accept': float(apply_gate(tab, adj, inv, n_win, FROZEN).mean()),
               'shape_accept': float(apply_gate(tab, adj, inv, n_win, SHAPE).mean()),
               'median_best_adj_r2': float(np.median(best_adj)),
               'median_minor_px': float(np.median(tab['minor'])),
               'median_axis_ratio': float(np.median(tab['axis_ratio']))}
        cases.append(row)
        print(f"  B: {lab:16s} frozen {row['frozen_accept']*100:5.1f}%  shape-only "
              f"{row['shape_accept']*100:5.1f}%  median best adjR2 {row['median_best_adj_r2']:+.3f}  "
              f"median minor {row['median_minor_px']:.5f} px", flush=True)

    threshold_ratio = FROZEN['minor_px'] / cases[0]['median_minor_px']
    auc_site = auc(runs['Khafre pyramid'], runs['Open plateau'])
    auc_noise = auc(runs['Khafre pyramid'], runs['White noise'])
    print(f'  B: frozen minor-axis threshold is {threshold_ratio:.0f}x the measured vectors; '
          f'modal score AUC pyramid vs plateau {auc_site:.3f}, vs white noise {auc_noise:.3f}', flush=True)

    # ---- C: does gating change the depth axis? ----
    tab, adj, inv, n_win, _ = modal_run(modal, qk_s[:N_POS_FOCUS])
    win_ok = apply_gate(tab, adj, inv, n_win, SHAPE)
    acc = win_ok.reshape(N_POS_FOCUS, qk_s.shape[1] - W + 1).any(axis=1)
    if acc.sum() >= 10:
        g_corr, g_max = at_lag(sc[acc], int(round(lattice_rec / step)))
    else:
        g_corr, g_max = float('nan'), float('nan')
    print(f'  C: {int(acc.sum())}/{N_POS_FOCUS} positions pass the shape gate; '
          f'their depth response repeats with corr {g_corr:.6f} (max diff {g_max:.1e})', flush=True)

    # ---- D: what sets the depth scale ----
    peak = z[np.argmax(sc, axis=1)]
    lam_rows = [{'lambda_eff_m': float(L), 'repeat_m': float(lattice_rec * L / lam_s),
                 'median_peak_m': float(np.median(peak) * L / lam_s)}
                for L in (0.03107, 0.2405679137, 1.0, 7.6)]
    for r in lam_rows:
        print(f"  D: lambda_eff {r['lambda_eff_m']:.5f} m -> repeat {r['repeat_m']:8.1f} m, "
              f"median peak {r['median_peak_m']:8.1f} m", flush=True)

    # ---- a true mechanical resonance, against the band a dwell can carry ----
    coh_s = 0.98                                   # measured half-coherence lag, R9
    aperture_s = float(khafre['meta']['centres_span_s'])
    vs = 2200.0                                    # shear speed in Giza limestone, m/s
    chambers = [("King's Chamber roof", 43.0), ("Subterranean chamber", 30.0),
                ("Typical mastaba shaft", 12.0), ("Claimed deep structure", 648.0)]
    band_hi = 1.0 / (2 * coh_s)                    # fastest rate independent looks can be formed
    band_lo = 1.0 / aperture_s                     # slowest, set by the record length
    res = []
    for nm, d in chambers:
        f = vs / (4 * d)
        res.append({'name': nm, 'depth_m': d, 'quarter_wave_hz': f,
                    'in_band': bool(band_lo <= f <= band_hi), 'too_fast_by': float(f / band_hi)})
        print(f"    {nm:26s} {d:6.1f} m -> {f:8.2f} Hz  "
              f"{'in band' if res[-1]['in_band'] else str(round(f / band_hi)) + 'x too fast'}", flush=True)
    print(f'  resonance band a dwell can carry: {band_lo:.3f} to {band_hi:.3f} Hz', flush=True)

    # ---------------- figures ----------------
    fig, axs = plt.subplots(1, 2, figsize=(12.2, 4.1))
    axs[0].plot(z, np.median(sc, axis=0), lw=1.6, color=viz.SERIES[0])
    for j in range(1, 4):
        axs[0].axvline(j * lattice_rec, color=viz.SERIES[7], lw=1.1, ls='--')
    axs[0].axvline(lattice_rec / 2, color=viz.INK2, lw=1.0, ls=':')
    axs[0].text(lattice_rec / 2, 0.985, ' unique interval ends', transform=axs[0].get_xaxis_transform(),
                fontsize=8, color=viz.INK2, va='top')
    axs[0].text(lattice_rec, 0.02, f'  repeats every {lattice_rec:.1f} m', transform=axs[0].get_xaxis_transform(),
                fontsize=8, color=viz.SERIES[7])
    axs[0].set_xlabel('model depth (m)'); axs[0].set_ylabel('median fit score')
    axs[0].set_title('The same depth answer, over and over', loc='left')
    s = slice(None, None, 37)
    axs[1].scatter(sc[:, :sc.shape[1] - best_lag].ravel()[s],
                   sc[:, best_lag:].ravel()[s], s=2, color=viz.SERIES[0], alpha=0.35, lw=0)
    lim = [float(sc.min()), float(sc.max())]
    axs[1].plot(lim, lim, color=viz.AXIS, lw=1, ls='--')
    axs[1].set_xlabel('score at depth z'); axs[1].set_ylabel(f'score at depth z + {meas_rec:.1f} m')
    axs[1].set_title(f'Correlation {best_corr:.5f}', loc='left')
    fig.tight_layout()
    repeat_fig = {'file': viz.finish(fig, f'{fd}/r15_repeat.png'), 'caption':
                  f'<b>The depth axis is periodic by construction.</b> Left: the median fit score of '
                  f'{N_POS_FOCUS} real Khafre positions, scored by the published depth stage over three '
                  f'repeat lengths. Dashed lines mark multiples of {lattice_rec:.1f} m. Right: every score '
                  f'plotted against the score one repeat deeper. They agree to a correlation of '
                  f'{best_corr:.5f}, so a peak at any depth has an equally good twin {meas_rec:.1f} m below '
                  f'it, and another below that. Nothing in the data causes this; the two model columns at '
                  f'those depths span the same plane.'}

    fig, axs = plt.subplots(1, 2, figsize=(12.2, 4.0))
    labs = [c['label'] for c in cases]
    x = np.arange(len(labs))
    axs[0].bar(x - 0.19, [c['frozen_accept'] * 100 for c in cases], width=0.36,
               color=viz.SERIES[7], label='frozen gate, as published')
    axs[0].bar(x + 0.19, [c['shape_accept'] * 100 for c in cases], width=0.36,
               color=viz.SERIES[0], label='shape criteria only')
    axs[0].set_xticks(x); axs[0].set_xticklabels(labs, fontsize=8.5)
    axs[0].set_ylabel('windows accepted (%)'); axs[0].legend(fontsize=8)
    axs[0].set_title('What the modal gate lets through', loc='left')
    bins = np.linspace(-0.4, 1.0, 60)
    for i, lab in enumerate(labs):
        axs[1].hist(runs[lab], bins=bins, density=True, histtype='step', lw=1.8,
                    color=viz.SERIES[[0, 1, 2, 3][i]], label=lab)
    axs[1].axvline(FROZEN['adj_r2'], color=viz.INK2, ls=':', lw=1.2)
    axs[1].text(FROZEN['adj_r2'], 0.97, ' gate', transform=axs[1].get_xaxis_transform(),
                fontsize=8, color=viz.INK2, va='top')
    axs[1].set_xlabel('best adjusted $R^2$ over modes 1 to 6')
    axs[1].set_ylabel('density'); axs[1].legend(fontsize=8)
    axs[1].set_title(f'Pyramid against plateau: AUC {auc_site:.3f}', loc='left')
    fig.tight_layout()
    gate_fig = {'file': viz.finish(fig, f'{fd}/r15_gate.png'), 'caption':
                f'<b>The gate does not know where the pyramid is.</b> Left: the share of windows admitted, '
                f'at the published settings and again with the absolute size threshold removed. The frozen '
                f'threshold asks for a minor axis of {FROZEN["minor_px"]} pixels, about '
                f'{threshold_ratio:.0f} times larger than the displacements this acquisition actually '
                f'contains, so it admits nothing. Right: the modal score itself, where the pyramid and open '
                f'plateau distributions sit on top of one another, AUC {auc_site:.3f} against 0.5 for a coin. '
                f'A random walk with no structure whatsoever scores far higher than either, because the fit '
                f'rewards smoothness along the sweep rather than anything underground.'}

    m = {'kz_spacing_rad_m': dkz, 'kz_spacing_spread': uniformity,
         'lattice_repeat_m': lattice_rec, 'unique_interval_m': lattice_rec / 2,
         'measured_repeat_m': meas_rec, 'lambda_eff_m': lam_s,
         'column_space_residual': factor_resid, 'repeat_correlation': rep_corr,
         'repeat_max_diff': rep_max, 'best_repeat_correlation': best_corr,
         'half_period_correlation': half_corr,
         'frozen_gate': FROZEN, 'shape_gate': SHAPE, 'window_w': W, 'modes': list(MODES),
         'cases': cases, 'minor_threshold_ratio': threshold_ratio,
         'auc_pyramid_vs_plateau': auc_site, 'auc_pyramid_vs_noise': auc_noise,
         'gated_positions': int(acc.sum()), 'gated_total': N_POS_FOCUS,
         'gated_repeat_correlation': g_corr, 'gated_repeat_max_diff': g_max,
         'lambda_rows': lam_rows, 'resonance_band_hz': [band_lo, band_hi],
         'resonance': res, 'shear_speed_m_s': vs, 'aperture_span_s': aperture_s,
         'coherence_half_life_s': coh_s,
         'positions_modal': int(qk_s.shape[0]), 'positions_focus': N_POS_FOCUS,
         'package_sha256': EXPECTED, 'runtime_s': time.time() - t0}

    noise = next(c for c in cases if c['label'] == 'White noise')
    pyr = next(c for c in cases if c['label'] == 'Khafre pyramid')
    plat = next(c for c in cases if c['label'] == 'Open plateau')
    walk_row = next(c for c in cases if c['label'] == 'Random walk')
    summary = {
        'id': TEST, 'order': 25, 'tag': 'mechanism', 'eyebrow': '15 · The modal branch',
        'title': 'The modal gate picks which pixels get a depth, not which depth is right',
        'question': 'Does the elliptical modal branch, the spring-like part of the method, resolve the repeating depth axis?',
        'finding': (
            f"No, and it cannot. The depth stage fits each measured vector sequence to cos(Kz z) and "
            f"sin(Kz z). "
            f"The steering wavenumbers sit on a near-regular ladder of {dkz:.5f} radians per metre, so the two "
            f"model columns at any depth z and at z + {lattice_rec:.1f} m span the same plane, and a "
            f"least-squares fit cannot tell them apart. "
            f"Run on real measured vectors the protocol's own focus stage scores one repeat apart at a "
            f"correlation of {best_corr:.5f}, which leaves only {lattice_rec/2:.1f} m of depth ever "
            f"distinguishable, and its own section 11.4 says so. "
            f"The modal branch runs before that stage and only decides which windows are admitted, so it "
            f"changes the population and leaves the axis untouched: windows that pass the gate repeat at "
            f"{g_corr:.5f}. "
            f"The gate is not a structure detector either, since at its frozen settings it admits "
            f"{pyr['frozen_accept']*100:.0f} per cent here, asking for a minor axis {threshold_ratio:.0f} times "
            f"larger than this acquisition's displacements, and with that size threshold removed it admits "
            f"{pyr['shape_accept']*100:.1f} per cent on the pyramid against {plat['shape_accept']*100:.1f} per "
            f"cent on open plateau, separating them at AUC {auc_site:.3f} where 0.5 is a coin toss, while a "
            f"random walk containing nothing at all passes {walk_row['shape_accept']*100:.0f} per cent of the "
            f"time. "
            f"Its harmonic index counts cycles across sub-aperture position rather than vibrations in time, as "
            f"the protocol's own symbol table states, and a real mechanical resonance fails separately: a "
            f"chamber {chambers[1][1]:.0f} m down rings near {res[1]['quarter_wave_hz']:.0f} Hz while a dwell "
            f"that decorrelates in about a second carries only {band_lo:.2f} to {band_hi:.2f} Hz."),
        'limitations': (
            'The exact-repetition argument holds for an evenly spaced ladder of steering wavenumbers. This '
            f'geometry is even to {uniformity*100:.1f} per cent, so the repeat is very close but not perfect, '
            'which is why the correlation is quoted rather than an identity. The gate rates come from our '
            'registration vectors, not theirs, because their own registration returns whole pixels and would '
            'yield nothing to gate at all; using ours is the more favourable choice. The acquisition is an '
            'ICEYE dwell, while their frozen values describe a Capella scene, so the absolute size threshold '
            'may suit their data better than ours, and the shape-only comparison is included for that reason. '
            'The resonance figures use one representative shear speed and the simplest quarter-wave model, so '
            'they bound the band rather than predict a spectrum.'),
        'method': (
            f"Package v1.7 was used unmodified, with the sha256 of its modal and focus stages checked before "
            f"import. Real sub-pixel registration vectors were measured on the Giza acquisition with the "
            f"protocol's frozen bank design of 50 pairs, half-band masks and 88 Hz separation, at "
            f"{qk_s.shape[0]:,} positions on Khafre and the same number on open plateau. Those vectors were "
            f"passed to its own modal fit at window 25 over modes 1 to 6, and to its own depth focus across "
            f"three repeat lengths. Two gates are reported: the frozen Branch B rule of adjusted R-squared at "
            f"least {FROZEN['adj_r2']}, axis ratio at least {FROZEN['axis_ratio']} and minor axis at least "
            f"{FROZEN['minor_px']} pixels, and the same rule with the absolute size threshold removed."),
        'figures': [repeat_fig, gate_fig], 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    with open(os.path.join(fd, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: m[k] for k in ('lattice_repeat_m', 'measured_repeat_m', 'unique_interval_m',
                                        'best_repeat_correlation', 'gated_repeat_correlation',
                                        'auc_pyramid_vs_plateau', 'minor_threshold_ratio')}, indent=1))


if __name__ == '__main__':
    main()
