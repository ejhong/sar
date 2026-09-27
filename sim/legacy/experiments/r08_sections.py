"""R8: the picture, from real data.

Statistics decide the question but they are hard to look at. This is the view the published
slides showed: a vertical section, horizontal axis across the ground, vertical axis labelled as
depth, drawn in the same colour scheme. Three patches of one real acquisition are rendered side
by side on a shared scale, and the depth axis is extended past the ambiguity limit of the
sub-aperture bank so that what happens there is visible rather than asserted.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

from real_common import run_patch, RESULTS
from sarsim import viz
from sarsim.tomo import nyquist_depth

TEST = 'r08_sections'
PANELS = [('khafre', 'Khafre pyramid'), ('khufu', 'Khufu pyramid'), ('desert_west', 'Open plateau')]
PERIODS = 3.0          # how many ambiguity periods of depth to draw
RANGE_BAND = 12        # grid columns averaged to suppress speckle


def section(run, periods=PERIODS, band=RANGE_BAND):
    """Azimuth-by-depth section, averaged over a band of range columns.

    Two depth scales matter and they differ by a factor of two. The steering wavenumbers are
    uniformly spaced, so |h(z)| is exactly periodic with 2 pi / dKz: below that the picture is
    the same picture again. The derivative protocol's real two-component fit cannot tell z from
    -z and so searches only half of it, which is the quantity `nyquist_depth` returns.
    """
    n_az, n_rg = run['grid_rows'].size, run['grid_cols'].size
    q = run['q'].astype(float).reshape(n_az, n_rg, -1, 2)
    kz = run['kz']
    z_half = nyquist_depth(np.sort(kz))
    z_amb = 2.0 * z_half                      # the exact repetition period of |h(z)|
    z = np.linspace(z_amb / 240, periods * z_amb, int(240 * periods))
    amp = run['amplitude'].astype(float).reshape(n_az, n_rg)
    centre = int(np.argmax(gaussian_filter(amp.mean(axis=0), 3)))
    lo, hi = max(centre - band // 2, 0), min(centre + band // 2, n_rg)
    strip = q[:, lo:hi].reshape(-1, q.shape[2], 2)
    A = np.exp(-1j * np.outer(kz, z))
    power = (np.abs((strip[..., 0] + 1j * strip[..., 1]) @ A) ** 2).reshape(n_az, hi - lo, z.size)
    return z, z_amb, power.mean(axis=1), amp[:, lo:hi].mean(axis=1)


def repetition_correlation(z, power, period):
    """Correlation of the log section with itself shifted by one period; 1.0 means an exact copy."""
    step = int(round(period / (z[1] - z[0])))
    if step < 1 or step >= power.shape[1]:
        return float('nan')
    a = np.log10(power[:, :-step] + 1e-20)
    b = np.log10(power[:, step:] + 1e-20)
    a = a - a.mean()
    b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a ** 2).sum() * (b ** 2).sum()))


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)
    runs = {k: run_patch('giza', k, 'paper', verbose=False) for k, _ in PANELS}
    data = {k: section(runs[k]) for k, _ in PANELS}
    spacing = runs['khafre']['meta']['acquisition']['azimuth_spacing_m'] * runs['khafre']['meta']['stride'][0]
    figs = []

    # shared colour scale across every panel, so the panels are comparable by eye
    stacked = np.concatenate([np.log10(d[2] + 1e-20).ravel() for d in data.values()])
    vmin, vmax = np.percentile(stacked, [4, 99.3])

    fig, axs = plt.subplots(2, 3, figsize=(14.2, 7.6), sharey='row',
                            gridspec_kw={'height_ratios': [1, 4.2], 'hspace': 0.06, 'wspace': 0.07})
    for col, (key, label) in enumerate(PANELS):
        z, z_amb, power, amp = data[key]
        x = (np.arange(power.shape[0]) - power.shape[0] / 2) * spacing
        axs[0, col].plot(x, 20 * np.log10(np.maximum(amp, 1e-3)), color=viz.INK2, lw=1)
        axs[0, col].set_xlim(x[0], x[-1]); axs[0, col].set_xticks([])
        axs[0, col].set_title(label, loc='left', fontsize=10.5)
        if col == 0:
            axs[0, col].set_ylabel('surface\namplitude (dB)', fontsize=8.5)
        im = axs[1, col].imshow(gaussian_filter(np.log10(power + 1e-20), (1.2, 1.6)).T,
                                origin='upper', aspect='auto', cmap='jet', vmin=vmin, vmax=vmax,
                                extent=[x[0], x[-1], z[-1], z[0]], interpolation='bicubic')
        for k in range(1, int(PERIODS) + 1):
            axs[1, col].axhline(k * z_amb, color='white', lw=0.9, ls='--', alpha=0.85)
        axs[1, col].axhline(z_amb / 2, color='white', lw=0.7, ls=':', alpha=0.6)
        axs[1, col].set_xlabel('ground distance across the patch (m)')
        if col == 0:
            axs[1, col].set_ylabel('model depth (m), $\\lambda_s$ = 0.48 m')
            axs[1, col].text(x[0] + 4, z_amb * 0.96, 'the picture repeats below here', color='white',
                             fontsize=8.5, va='bottom')
            axs[1, col].text(x[0] + 4, z_amb * 0.46, 'real-fit ambiguity limit', color='white',
                             fontsize=7.5, va='bottom', alpha=0.85)
    fig.suptitle('One ICEYE dwell of Giza, 27 August 2025: what the method draws beneath three patches',
                 x=0.012, ha='left', fontsize=11.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 0.94, 0.955])
    cax = fig.add_axes([0.955, 0.12, 0.012, 0.52])
    fig.colorbar(im, cax=cax, label='$\\log_{10}|h(z)|^{2}$')
    figs.append({'file': viz.finish(fig, f'{fd}/r08_sections.png'), 'caption':
                 '<b>The view the slides showed, from the real acquisition.</b> Top: the surface radar brightness '
                 'along each strip, so the monument is visible where it actually is. Bottom: the depth output '
                 'beneath it, in the same jet colour scheme and smoothing the published figures used, on one shared '
                 'scale. Vertical banding of the kind read as wells and levels is present under all three patches, '
                 'including open plateau. The dashed lines mark where the steering basis comes back to itself: '
                 'below the first one the section is not merely similar but identical, correlating with the band '
                 'above it at 1.000. The fainter dotted line is the shallower limit at which the derivative '
                 "protocol's real two-component fit can no longer tell one depth from its negative."})

    # --- the same picture, relabelled
    key = 'khafre'
    z, z_amb, power, amp = data[key]
    fig, axs = plt.subplots(1, 2, figsize=(12.2, 5.0), sharey=False)
    x = (np.arange(power.shape[0]) - power.shape[0] / 2) * spacing
    img = gaussian_filter(np.log10(power + 1e-20), (1.2, 1.6)).T
    for ax, lam, title in [(axs[0], 0.48, 'as the 2022 paper labels it, $\\lambda_s$ = 0.48 m'),
                           (axs[1], 0.48 * 648 / (PERIODS * z_amb), f'relabelled to reach 648 m')]:
        scale = lam / 0.48
        ax.imshow(img, origin='upper', aspect='auto', cmap='jet', vmin=vmin, vmax=vmax,
                  extent=[x[0], x[-1], z[-1] * scale, z[0] * scale], interpolation='bicubic')
        ax.set_title(title, loc='left', fontsize=10)
        ax.set_xlabel('ground distance across the patch (m)')
        ax.set_ylabel('model depth (m)')
    fig.suptitle('Khafre: identical measurement, two depth labels', x=0.012, ha='left', fontsize=11.5, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    lam_needed = 0.48 * 648 / (PERIODS * z_amb)
    figs.append({'file': viz.finish(fig, f'{fd}/r08_relabel.png'), 'caption':
                 f'<b>Same pixels, same numbers, different axis.</b> The left panel is the Khafre section with the '
                 f'depth axis the 2022 paper implies. The right panel is the identical array, relabelled by '
                 f'declaring a sound wavelength of {lam_needed:.1f} m instead of 0.48 m, which places the bottom of '
                 f'the picture at 648 m. Nothing was reprocessed between the two panels. The depth of a feature in '
                 f'this method is a choice of constant, not a measurement.'})

    m = {'panels': [{'patch': k, 'label': lab,
                     'repetition_period_m': float(data[k][1]),
                     'real_fit_ambiguity_limit_m': float(data[k][1]) / 2,
                     'repetition_correlation': repetition_correlation(data[k][0], data[k][2], data[k][1]),
                     'half_period_correlation': repetition_correlation(data[k][0], data[k][2], data[k][1] / 2),
                     'ground_extent_m': float(data[k][2].shape[0] * spacing),
                     'range_columns_averaged': RANGE_BAND} for k, lab in PANELS],
         'periods_drawn': PERIODS, 'lambda_s_for_648m': float(lam_needed),
         'shared_colour_scale': [float(vmin), float(vmax)], 'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 18, 'tag': 'artifact', 'eyebrow': '08 · The picture',
        'title': 'Rendered like the published slides, the real data shows the same thing everywhere',
        'question': 'What does this actually look like, and can a reader tell whether it is working?',
        'finding': (f"It looks like the published slides, under the pyramids and under empty plateau alike. "
                    f"Vertical banded columns of the kind read as wells appear in all three patches of the same "
                    f"acquisition. The steering wavenumbers are uniformly spaced to within a hundred-thousandth, so "
                    f"the section is exactly periodic: at Khafre it repeats every {data['khafre'][1]:.1f} m with a "
                    f"correlation of 1.000, meaning every structure drawn below that depth is a copy of one above "
                    f"it. Relabelling the axis by declaring a sound wavelength of {lam_needed:.1f} m instead of "
                    f"0.48 m puts the bottom of the same unchanged array at 648 m."),
        'limitations': ('These are our bank parameters, not the undisclosed ones behind the published images, so '
                        'this is a like-for-like rendering of the method rather than a reproduction of a specific '
                        'figure. The strips are 12 range columns wide and span about 165 m of ground, which is '
                        'less than the full base of Khafre. Smoothing and colour scaling follow the published '
                        'presentation, which flatters any such image.'),
        'method': ('Cached trajectories from the within-image control runs, refocused on a depth grid extended to '
                   'three ambiguity periods, averaged over a band of range columns centred on the brightest part of '
                   'each patch, smoothed and drawn in jet on one shared logarithmic scale.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(m, indent=1))


if __name__ == '__main__':
    main()
