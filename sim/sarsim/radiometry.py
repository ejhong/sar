"""The signal-to-noise ratio per resolution cell, from an acquisition's own radiometry (Phase 2).

The oracle bounds (P2-27, P2-29, P2-32, P2-36, P2-38) need the ratio of the ground's backscatter to the receiver noise
per resolution cell. It is read here, once, from the radiometry block the desktop writes into each acquisition record
(`sim/katabasis/compose/acquisition_radiometry.py`, from the product and its delivery files):

- the ground: sigma0 over 10 m cells of the Giza site in the image itself (calibration factor x mean |z|^2 x sin of the
  local incidence), its median measured;
- the noise: neither product carries a noise-equivalent sigma-zero, so ICEYE's specified range for Dwell is used
  (Product Documentation 6.0.0, Table 2-11: -18 to -15 dB, scene centre), its best end for a generous bound; the darkest
  measured cells give only a ceiling on the noise (they hold noise and any signal left), which bounds the ratio from
  below, not above.

An earlier value, -26.7 dB for Dwell's best NESZ, had no traceable source and is withdrawn (audit, sixth review).
"""
from __future__ import annotations

import json
from pathlib import Path

RECORDS = Path(__file__).resolve().parents[2] / 'sites' / 'acquisitions'


def snr_per_cell(name: str, bright_ground_db: float = 0.0) -> dict:
    """The SNR per cell for acquisition `name`, in dB: `nominal`, the site's median ground over the specified best noise
    floor; `generous`, ground as bright as `bright_ground_db` (about the brightest natural ground, assumed) over the same;
    `specified_worst`, the median ground over the specified worst floor; and `measured_floor_lower`, the median ground over
    the measured noise ceiling (a lower bound on the ratio). With each input and its source."""
    rec = json.loads((RECORDS / f'{name}.json').read_text())
    rad = rec.get('radiometry')
    if rad is None:
        raise KeyError(f'{name} has no radiometry block (written on the desktop, where the products are)')
    nf = rad['noise_floor']
    best, worst = min(nf['specification']['nesz_db']), max(nf['specification']['nesz_db'])
    meas = nf['measured_upper_bound']
    ground = meas['sigma0_db_median']
    return {
        'nominal_db': ground - best,
        'generous_db': bright_ground_db - best,
        'specified_worst_db': ground - worst,
        'measured_floor_lower_db': ground - meas['sigma0_db_min'],
        'ground_sigma0_median_db': ground,
        'bright_ground_db_assumed': bright_ground_db,
        'nesz_specified_db': [best, worst],
        'nesz_in_product': bool(nf['in_product']),
        'nesz_source': nf['specification']['source'],
        'ground_source': meas['note'],
        'acquisition': name,
    }
