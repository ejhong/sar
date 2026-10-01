"""The signal-to-noise ratio per resolution cell: the measured backscatter over a specified noise floor (Phase 2).

The oracle bounds (P2-25, P2-27, P2-29, P2-30, P2-32, P2-36, P2-38) need the ratio of the ground's backscatter to the
receiver noise per resolution cell. Half of it is measured, half is not:

- the ground: sigma0 over 10 m cells of the Giza site in the image itself (calibration factor x mean |z|^2 x sin of the
  local incidence), its median measured, from the radiometry block the desktop writes into each acquisition record
  (`sim/katabasis/compose/acquisition_radiometry.py`);
- the noise: neither product carries a noise-equivalent sigma-zero (NESZ), so the noise is ICEYE's specification for
  the product's imaging mode, which has changed between versions of the product documentation (SPECIFICATION below).
  The darkest measured cells give only a ceiling on the noise (they hold noise and any signal left), which bounds the
  ratio from below, never from above.

So every SNR here is conditional: measured backscatter combined with an assumed, specification-based noise level, never
a measurement of this acquisition's noise. The scenarios are kept side by side. The headline is the most favourable
specified for the product's mode (the best end of the newest documentation), which is the conservative choice for an
exclusion: a lower noise can only help a detector.

History (audit, sixth and eighth reviews): an earlier value, -26.7 dB, was Dwell's best end in documentation 6.0.8,
applied to a Dwell Fine product; the sixth review, reading only version 6.0.0 (-18 to -15 dB for both modes), wrongly
called it untraceable and withdrew it; the eighth restores the version and mode history below.
"""
from __future__ import annotations

import json
from pathlib import Path

RECORDS = Path(__file__).resolve().parents[2] / 'sites' / 'acquisitions'

# ICEYE Product Documentation, Product Specification, 2. Imaging Modes, Table 2-11 ("Collection performance attributes
# for Dwell ... imaging modes"), each version as published; best (more negative) to worst, at scene centre.
SPECIFICATION = {
    '6.0.0': {
        'url': 'https://sar.iceye.com/6.0.0/productspecification/imagingmodes/',
        'quote': 'Noise Equivalent Sigma-Zero [dBm2/m2]: Dwell -18 to -15; Dwell Fine -18 to -15',
        'nesz_db': {'SpotlightDwell': (-18.0, -15.0), 'SpotlightDwellFine': (-18.0, -15.0)},
    },
    '6.0.8': {
        'url': 'https://sar.iceye.com/6.0.8/productspecification/imagingmodes/',
        'quote': 'Noise Equivalent Sigma-Zero, scene centre [dBm2/m2]: Dwell -26.7 to -15.6; Dwell Fine -23.7 to -12.6; '
                 'Dwell Precise -20.7 to -11.3 (NESZ ... scene-centre values; ranges from the best to the worst value '
                 'across performant incidence angles and the satellites producing each product)',
        'nesz_db': {'SpotlightDwell': (-26.7, -15.6), 'SpotlightDwellFine': (-23.7, -12.6),
                    'SpotlightDwellPrecise': (-20.7, -11.3)},
    },
}
NEWEST = '6.0.8'


def snr_per_cell(name: str, bright_ground_db: float = 0.0) -> dict:
    """The SNR per cell for acquisition `name`, in dB, each a scenario (measured ground over a specified noise):
    `headline_db`, the site's median ground over the best specified noise for the product's mode in the newest
    documentation (the conservative choice for an exclusion); `bright_db`, ground as bright as `bright_ground_db` (about
    the brightest natural ground, assumed) over the same; `scenarios`, every version's best and worst for the mode; and
    `measured_floor_lower_db`, the median ground over the measured noise ceiling (a lower bound on the ratio). With each
    input and its source."""
    rec = json.loads((RECORDS / f'{name}.json').read_text())
    rad = rec.get('radiometry')
    if rad is None:
        raise KeyError(f'{name} has no radiometry block (written on the desktop, where the products are)')
    mode = rec['source']['mode']
    nf = rad['noise_floor']
    meas = nf['measured_upper_bound']
    ground = meas['sigma0_db_median']
    scen = []
    for ver, sp in SPECIFICATION.items():
        best, worst = sp['nesz_db'][mode]
        scen += [{'label': f'documentation {ver}, best', 'nesz_db': best, 'snr_db': ground - best, 'source': sp['url']},
                 {'label': f'documentation {ver}, worst', 'nesz_db': worst, 'snr_db': ground - worst, 'source': sp['url']}]
    best_new = SPECIFICATION[NEWEST]['nesz_db'][mode][0]
    any_mode_best = min(v[0] for v in SPECIFICATION[NEWEST]['nesz_db'].values())
    inc = rad['product']
    return {
        'headline_db': ground - best_new,
        'bright_db': bright_ground_db - best_new,
        'any_dwell_mode_best_db': ground - any_mode_best,
        'scenarios': scen,
        'measured_floor_lower_db': ground - meas['sigma0_db_min'],
        'ground_sigma0_median_db': ground,
        'bright_ground_db_assumed': bright_ground_db,
        'mode': mode,
        'nesz_headline_db': best_new,
        'nesz_headline_source': f"ICEYE Product Documentation {NEWEST}, Table 2-11, {mode}: best end "
                                f"({SPECIFICATION[NEWEST]['url']})",
        'nesz_in_product': bool(nf['in_product']),
        'scene_centre_note': f"the specification's values are at scene centre; the site lies near the centre of a swath "
                             f"whose incidence runs {inc['incidence_near']:.1f} to {inc['incidence_far']:.1f} deg",
        'ground_source': meas['note'],
        'status': 'conditional: measured backscatter over an assumed, specification-based noise level',
        'acquisition': name,
    }
