"""P2-00 · Where on the Giza plateau could the ground shake hard enough for the radar?

    uv run python experiments/p2_00_ambient_giza.py

The ambient motion the Phase 2 radar would look at, from sourced values
only (sites/ambient.json):

* microseisms: 0.049 um/s (median, 0.1-0.3 Hz), uniform over the scene;
* traffic: a vehicle on the nearest road (OpenStreetMap, sites/giza/roads.json)
  shakes the ground at the FTA's level 15 m away, a typical bus or truck
  36 um/s, one over a bump 100 um/s, falling off as (15/r)^0.5 e^(-a(r-15))
  (Caltrans Eq. 6) with a between 0 and the marly-chalk 0.1 per metre.

Worst case by construction: every road carries a truck at its point nearest
each spot. Within 15 m of a road, levels are held at the 15 m value (the
reference distance; closer is higher and unsourced). The result: where on the
plateau, if anywhere, ambient motion reaches the 77 um/s floor measured on
the real ICEYE dwell.
"""
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from katabasis.runs import Run

SITES = Path(__file__).resolve().parents[2] / 'sites'
FLOOR = 77e-6              # m/s: the velocity floor measured on the real Giza dwell (legacy R4)
BOX = (-640.0, 960.0, -760.0, 640.0)   # the Giza site's extent, x0 x1 y0 y1
H = 2.0                    # m, map spacing


def densify(pts: np.ndarray, step: float = 1.0) -> np.ndarray:
    out = [pts[:1]]
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(1, int(np.ceil(np.linalg.norm(b - a) / step)))
        out.append(a + (b - a) * (np.arange(1, n + 1) / n)[:, None])
    return np.concatenate(out)


def main():
    lib = json.loads((SITES / 'ambient.json').read_text())
    roads = json.loads((SITES / 'giza' / 'roads.json').read_text())
    site = json.loads((SITES / 'giza' / 'site.json').read_text())
    cul, prop = lib['cultural'], lib['propagation']
    v15 = {'typical': cul['bus_or_truck_typical']['value'], 'bump': cul['bus_or_truck_over_bump']['value']}
    alphas = {'rock': 0.0, 'marly chalk': prop['damping_marly_chalk_per_m']['value']}
    gamma = prop['geometric_exponent_surface_wave']['value']
    micro = lib['regional']['microseism_vertical_0.1_0.3_hz']['value']
    params = {'floor_m_s': FLOOR, 'box': BOX, 'spacing_m': H, 'roads': roads['source'], 'osm_timestamp': roads['osm_timestamp'],
              'levels_15m_m_s': v15, 'damping_per_m': alphas, 'geometric_exponent': gamma, 'microseism_m_s': micro,
              'worst_case': 'a vehicle on every road at its point nearest each spot; held at the 15 m level within 15 m'}
    with Run('p2_00_ambient_giza', 'Ambient motion on the Giza plateau against the radar floor', params) as run:
        veh = [np.array(r['points']) for r in roads['roads'] if r['kind'] == 'vehicle']
        pts = np.concatenate([densify(p) for p in veh])
        tree = cKDTree(pts)
        x = np.arange(BOX[0], BOX[1] + H / 2, H)
        y = np.arange(BOX[2], BOX[3] + H / 2, H)
        X, Y = np.meshgrid(x, y, indexing='xy')
        dist, _ = tree.query(np.stack([X.ravel(), Y.ravel()], 1))
        dist = dist.reshape(X.shape)
        r = np.maximum(dist, 15.0)
        cases = {}
        for vk, v in v15.items():
            for ak, a in alphas.items():
                vel = v * (15.0 / r) ** gamma * np.exp(-a * (r - 15.0))
                above = vel >= FLOOR
                # the distance from a road within which the floor is reached (analytic, for the strip width)
                rr = np.linspace(15, 400, 20000)
                vv = v * (15.0 / rr) ** gamma * np.exp(-a * (rr - 15.0))
                reach = float(rr[vv >= FLOOR].max()) if v >= FLOOR else 0.0
                cases[f'{vk} · {ak}'] = {'level_15m_um_s': v * 1e6, 'damping_per_m': a, 'reach_m': reach,
                                         'area_fraction_above_floor': float(above.mean())}
        # the monuments and the claim
        places = {s['id']: s['shape']['centre'][:2] for s in site['structures']}
        claim = json.loads((SITES / 'bench-khafre-claim' / 'site.json').read_text())
        near = {k: float(tree.query(np.array(c))[0]) for k, c in places.items()}
        # the roads, clipped and thinned for the plan view
        plan = []
        for rd in roads['roads']:
            p = np.array(rd['points'])
            inside = (p[:, 0] >= BOX[0]) & (p[:, 0] <= BOX[1]) & (p[:, 1] >= BOX[2]) & (p[:, 1] <= BOX[3])
            if inside.sum() >= 2:
                plan.append({'kind': rd['kind'], 'class': rd['class'], 'points': p[inside].round(1).tolist()})
        worst = cases['bump · rock']
        run.save({
            'cases': cases, 'distance_to_nearest_road_m': near, 'microseism_over_floor': micro / FLOOR,
            'plan': {'box': BOX, 'roads': plan,
                     'pyramids': [{'id': s['id'], 'centre': s['shape']['centre'][:2], 'base': s['shape']['base']}
                                  for s in site['structures'] if s['shape']['type'] == 'pyramid']},
            'claim_site': claim.get('name'),
            'finding': (f"Even with a truck on every road at its nearest point, ambient motion reaches the radar's 77 um/s "
                        f"floor only within {worst['reach_m']:.0f} m of a road, and only for a truck over a bump "
                        f"({100 * worst['area_fraction_above_floor']:.0f}% of the plateau); a typical truck never does. "
                        f"Khafre's centre is {near['khafre']:.0f} m from the nearest mapped road; the microseisms are "
                        f"{FLOOR / micro:,.0f} times below the floor everywhere."),
        })


if __name__ == '__main__':
    main()
