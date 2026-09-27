"""R14: the budget. Four independent requirements, each measured against what is available.

Every earlier test answers one question. This puts them on one axis. For the published claim to
be true, four separate things must hold at once. Each is checked here against a measured or
textbook quantity, and each fails on its own, by the factor shown. They are independent: fixing
any one leaves the others untouched.
"""
import json
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from real_common import RESULTS
from sarsim import viz

TEST = 'r14_budget'
C = 299792458.0
F0 = 9.6e9
# Dry limestone at X-band. Permittivity and loss tangent from standard dielectric compilations
# for dry carbonate rock; the range covers the spread usually quoted.
LIMESTONE = dict(eps=(6.0, 9.0), tand=(0.008, 0.03))
SAND = dict(eps=(2.5, 3.2), tand=(0.002, 0.008))
RAYLEIGH_M_S = 1800.0


def penetration_m(eps, tand, f=F0):
    """Amplitude 1/e depth for a low-loss dielectric, delta = 1 / alpha."""
    alpha = (2 * np.pi * f / C) * np.sqrt(eps) * (tand / 2.0)
    return 1.0 / alpha


def main():
    viz.style()
    t0 = time.time()
    fd = os.path.join(RESULTS, TEST, 'figs')
    os.makedirs(fd, exist_ok=True)

    r09 = json.load(open(os.path.join(RESULTS, 'r09_coherence', 'summary.json')))['metrics']
    r04 = json.load(open(os.path.join(RESULTS, 'r04_velocity_floor', 'summary.json')))['metrics']
    record = r04['runs'][0]['series_span_s']
    floor = r04['best_floor_um_s']

    pen_lo = penetration_m(LIMESTONE['eps'][1], LIMESTONE['tand'][1])
    pen_hi = penetration_m(LIMESTONE['eps'][0], LIMESTONE['tand'][0])
    sand_hi = penetration_m(SAND['eps'][0], SAND['tand'][0])

    # to sample depth d, surface waves need wavelength about 3d; lateral resolution is about lam/2
    def resolution_m(depth):
        return 1.5 * depth

    legs = [
        dict(name='Radar reach',
             need='see 648 m of rock', need_v=648.0,
             have='X-band penetrates dry limestone', have_v=pen_hi, unit='m', bigger_is_better=True,
             note=f'Amplitude 1/e depth is {pen_lo * 100:.0f} to {pen_hi * 100:.0f} cm at 9.6 GHz. Nothing below '
                  f'the surface is illuminated, so every claim must come from surface motion. Dry sand reaches '
                  f'{sand_hi:.1f} m, which is where longer-wavelength penetration stories come from.'),
        dict(name='Record length',
             need='resolve the microseism band', need_v=150.0,
             have='one dwell, inside the coherent window', have_v=record, unit='s', bigger_is_better=True,
             note=f'A five-second period needs tens of cycles. Coherence halves at '
                  f'{r09["patches"][0]["half_coherence_s"]:.2f} s, which caps the usable record at '
                  f'{record:.2f} s.'),
        dict(name='Motion sensitivity',
             need='ambient ground motion', need_v=10.0,
             have='measured floor on the brightest targets', have_v=floor, unit='µm/s', bigger_is_better=False,
             note=f'Ambient motion in the microseism band is 0.1 to 10 µm/s. The measured floor is '
                  f'{floor:.0f} µm/s, and it is the same on monuments as on desert.'),
        dict(name='Lateral resolution',
             need='resolve a 10 m well', need_v=10.0,
             have='best possible from waves that reach 648 m', have_v=resolution_m(648.0), unit='m', bigger_is_better=False,
             note='Surface waves sample to about a third of a wavelength and resolve about half of one. Reaching '
                  '648 m needs kilometre wavelengths. This holds for a perfect dense array with unlimited record '
                  'and has nothing to do with radar.'),
    ]
    # Two legs want a large number and two want a small one; the shortfall is always the ratio
    # of the worse quantity to the better, so every bar reads in the same direction.
    for leg in legs:
        leg['shortfall'] = (leg['need_v'] / leg['have_v'] if leg['bigger_is_better']
                            else leg['have_v'] / leg['need_v'])

    fig, ax = plt.subplots(figsize=(11.2, 4.9))
    y = np.arange(len(legs))[::-1]
    for i, leg in enumerate(legs):
        yy = y[i]
        lo, hi = sorted((leg['need_v'], leg['have_v']))
        ax.plot([lo, hi], [yy, yy], color=viz.GRID, lw=7, solid_capstyle='butt', zorder=1)
        ax.scatter([leg['have_v']], [yy], s=95, color=viz.SERIES[0], zorder=3,
                   edgecolor=viz.SURFACE, linewidth=1.4, label='available' if i == 0 else None)
        ax.scatter([leg['need_v']], [yy], s=95, color=viz.SERIES[7], zorder=3, marker='D',
                   edgecolor=viz.SURFACE, linewidth=1.4, label='required' if i == 0 else None)
        ax.text(hi * 1.5, yy, f"short by {leg['shortfall']:,.0f}×", va='center', fontsize=10,
                color=viz.INK, family=viz.plt.rcParams['font.sans-serif'][0] if False else None)
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l['name']}\n{l['unit']}" for l in legs], fontsize=9.5)
    ax.set_xlabel('required against available, log scale, each leg in its own units')
    ax.set_xlim(0.05, 3e5)
    ax.set_title('Four independent requirements, four independent failures', loc='left')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(axis='x', alpha=0.35)
    fig.tight_layout()
    figs = [{'file': viz.finish(fig, f'{fd}/r14_budget.png'), 'caption':
             '<b>Every leg fails on its own.</b> For the published claim to hold, four separate things must be '
             'true at once. The radar would have to reach hundreds of metres into rock, where X-band reaches '
             'centimetres. The record would have to be long enough to resolve the microseism band, where one dwell '
             'gives under two seconds. The measurement would have to reach ambient ground motion, where the '
             'measured floor is nearly ten times above the top of that range. And surface waves able to sample '
             '648 m would have to resolve a ten-metre well, where the best they can do is hundreds of metres. '
             'These are independent: repairing any one changes none of the others.'}]

    m = {'legs': [{k: v for k, v in leg.items()} for leg in legs],
         'penetration_limestone_cm': [pen_lo * 100, pen_hi * 100],
         'penetration_dry_sand_m': sand_hi,
         'carrier_ghz': F0 / 1e9, 'rayleigh_m_s': RAYLEIGH_M_S, 'runtime_s': time.time() - t0}
    summary = {
        'id': TEST, 'order': 22, 'tag': 'artifact', 'eyebrow': '14 · The budget',
        'title': 'Four independent requirements, each missed by a wide margin',
        'question': 'What would have to be true for the published claim to hold, and by how much is each missed?',
        'finding': (f"X-band penetrates {pen_lo * 100:.0f} to {pen_hi * 100:.0f} cm of dry limestone, so nothing "
                    f"below the surface is illuminated and every claim must rest on surface motion. Resolving the "
                    f"microseism band needs a record of order 150 s; the coherent window allows {record:.2f} s. "
                    f"Ambient ground motion is 0.1 to 10 µm/s; the measured floor is {floor:.0f} µm/s. Surface "
                    f"waves long enough to sample 648 m cannot resolve better than about "
                    f"{resolution_m(648.0):.0f} m laterally, against wells said to be 10 m across. Each leg fails "
                    f"alone, and they are independent."),
        'limitations': ('Dielectric constants for dry carbonate rock span a range and the penetration figure moves '
                        'with them, but not by the three orders of magnitude the claim would need. The resolution '
                        'leg assumes surface-wave sensing; a body-wave method at the frequencies needed to resolve '
                        'ten metres cannot reach the depth at all, so that route is worse rather than better. The '
                        'record and sensitivity legs are measured on one sensor and one scene type.'),
        'method': ('Penetration from the low-loss attenuation constant at 9.6 GHz for dry limestone. Record length '
                   'and motion sensitivity measured in R9 and R4 on the real acquisition. Lateral resolution from '
                   'the standard surface-wave depth-sensitivity and half-wavelength resolution relations.'),
        'figures': figs, 'metrics': m, 'date': time.strftime('%Y-%m-%d'),
    }
    d = os.path.join(RESULTS, TEST)
    with open(os.path.join(d, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    for leg in legs:
        print(f"  {leg['name']:20s} need {leg['need_v']:>9,.1f} {leg['unit']:5s} "
              f"have {leg['have_v']:>9,.2f}  short by {leg['shortfall']:,.0f}x")


if __name__ == '__main__':
    main()
