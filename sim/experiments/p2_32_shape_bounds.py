"""P2-32 · Shape discrimination, bounded: can any way of reading one image tell one underground layout from another?

    uv run python experiments/p2_32_shape_bounds.py

The detection bounds (P2-25, P2-29, P2-30) limit telling a cavity from none, without needing its depth. The same
mathematics bounds telling two layouts apart: if the images two layouts produce are nearly the same in distribution, no
algorithm can tell them apart, however it reconstructs. Depth ambiguity does not rescue a reader: a uniformly reliable
shape mapper must distinguish two layouts even after the allowed vertical transformations, and the closest pair of
admissible depth variants is at least as hard to tell apart as the pair at the same depth, so the same-depth bound
computed here holds for every allowance.

Pairs (the shapes the published pictures are read for):
- a branching (T) tunnel against no tunnel;
- an L-shaped tunnel against a straight tunnel;
- a plain column (a shaft 3 m square, 4 to 30 m down) against the same column with a spiral of small rooms round it;
- an L-shaped against a branching tunnel; and an isolated room against none, for reference.

Each layout's imprint on the ground's motion is the lab's elastic solver: statically (settling under a uniform strain,
P2-28) for the ambient levels below 8 Hz; dynamically (P2-26's surface wave, 6-120 Hz) for strong shaking close by. The
difference between two layouts' imprints is what any reader would have to see. Three layers, as for detection:
- tight and certificate, through empty ground: TV(A, B) <= TV(A, none) + TV(none, B), each layout certified against the
  same empty-ground reference (P2-25's static_cases; P2-30's certificate from the peak phase over the whole image, with
  the common background's covariance floor). Subtracting two layouts' displacement fields does not remove the reference
  covariance: in the speckle model the same added phase is more or less distinguishable depending on the motion already
  there (an independent review's check), so the difference kernel's values are kept only as a diagnostic;
- the oracle told everything but which layout is there, on the difference itself: given the reflectivity and the
  background motion the two worlds' echoes differ in mean by exp(i k u_A) - exp(i k u_B), of size 2 |sin(k (u_A - u_B) / 2)|
  whatever the background, so the difference is exact there: TV <= erf(sqrt(d) / 2), d from the difference's phase
  energy at 30 dB per cell (P2-29), any texture, an upper bound (Jensen over the phases).
The directions of the ambient waves are sampled every 15 degrees; the peak over directions is enlarged by sec(15 deg)
(a function a + b cos 2 phi + c sin 2 phi sampled at that spacing can peak at most that factor above its samples).
Reported beside each: how many times the signal would have to grow for the oracle to tell the pair apart at 95% / 5%,
and how many times more signal discrimination needs than detecting the first layout. Only these layouts, this rock and
these sources are covered.
"""
import copy
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from katabasis.runs import RESULTS, Run, load, memo
from sarsim import information as inf
from sarsim import finite as fin
from sarsim.acquisition import DwellGeometry

RID = 'p2_32_shape_bounds'
SITES = Path(__file__).resolve().parents[2] / 'sites'
HELIX = [((4.5 * np.cos(a), 4.5 * np.sin(a), -6.0 - 2.5 * k), (2.0, 2.0, 2.0))
         for k, a in enumerate(np.deg2rad(np.arange(0, 360, 45)))]
EXTRA = {
    'none': [],
    'straight_tunnel': [((0.0, 0.0, -10.0), (28.0, 2.5, 2.5))],
    'column': [((0.0, 0.0, -17.0), (3.0, 3.0, 26.0))],
    'column_spiral': [((0.0, 0.0, -17.0), (3.0, 3.0, 26.0))] + [(tuple(map(float, c)), s) for c, s in HELIX],
}
PAIRS = [('branching_tunnel', 'none'), ('L_tunnel', 'straight_tunnel'), ('column_spiral', 'column'),
         ('L_tunnel', 'branching_tunnel'), ('room', 'none')]
DIRECTIONS = np.deg2rad(np.arange(0, 180, 15.0))
GUARD = 1 / np.cos(DIRECTIONS[1] - DIRECTIONS[0])   # sec(15 deg): the peak between sampled directions
GRID_ALLOWANCE = 0.11                              # P2-26's grid check: at most 10.5% in amplitude on halving the grid
HALF = 60.0
TAPER = (50.0, 59.0)
SNR_DB = 30.0
D_TARGET = 3.29
DYN_TAPER = (32.0, 39.0)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def static_layer(g, kern, bands, m25, snr, floor):
    """All three layers for one (difference) kernel, per ambient band."""
    axes = m25.pixel_axes(g, half=HALF)
    maps = m25.static_maps(kern, g, DIRECTIONS, axes)
    n = len(maps)
    dg = g.dr / np.sin(g.theta)
    k0 = 4 * np.pi / g.lam
    unit, pads = {}, {}
    for f in sorted({b['f_hz'] for b in bands}):
        pads[f] = m25.grid_pad(g, *maps[0].shape, f)
        unit[f] = [inf.fisher_grid(M, g, f, pad=pads[f]) for M in maps]
    fb_unit = float(np.mean([inf.fisher_bound(M, (g.dx, dg), g) for M in maps]))
    peak_unit = GUARD * k0 * max(float(np.abs(M).max()) for M in maps)
    worst_unit = sum(np.abs(M) for M in maps) / np.sqrt(n)
    rows = []
    for b in bands:
        e0 = b['strain_amplitude']
        parts = [(A * e0 ** 2 / n, (A + B) * e0 ** 2 / n) for A, B in unit[b['f_hz']]]
        eps = inf.remainder(e0 * worst_unit, g, pads[b['f_hz']])
        EF, kl, x = m25.ambient_kl(parts, eps)
        tight = m25.summary(EF, kl)
        d = 0.5 * snr * fb_unit * e0 ** 2
        Q = fin.peak_phase_energy(m_bins(g), peak_unit * e0)
        cert = fin.certificate_from_energy(Q, floor)
        rows.append({'case': b['case'], 'strain_amplitude': e0, 'f_hz': b['f_hz'], 'tight_ceiling': tight['ceiling'],
                     'tight_fisher': tight['fisher'], 'remainder_eps': eps,
                     'oracle_d': d, 'oracle_tv': fin.oracle_tv_from_mean_energy(d),
                     'oracle_growth_to_target': float(D_TARGET / np.sqrt(2 * d)) if d > 0 else None,
                     'certificate_tv': cert['tv_upper'], 'certificate_tpr_at_005': min(1.0, 0.05 + cert['tv_upper']),
                     'phase_energy_unit': fb_unit, 'peak_phase_unit_rad': peak_unit})
    return rows, fb_unit


def m_bins(g):
    src = json.loads((SITES / 'acquisitions' / 'giza-20250827.json').read_text())['source']['shape']
    cnt = lambda n, fr: int(np.count_nonzero(np.abs(np.fft.fftfreq(n)) <= fr / 2))
    return cnt(src[0], g.band_frac) * cnt(src[1], g.band_frac_r)


def background_floor(p204, amb, g):
    """P2-30's floor: the three regional Rayleigh harmonics' LOS envelope as a phase over the whole image."""
    reg = amb['regional']
    hv = p204['hv']
    vel = [(reg['microseism_vertical_0.1_0.3_hz']['value'], 0.2), (reg['vertical_1_3_hz']['value'], 2.0),
           (reg['vertical_3_8_hz']['value'], 5.0)]
    peak_up = sum(np.sqrt(2) * v / (2 * np.pi * f) for v, f in vel)
    phase = 4 * np.pi / g.lam * np.sqrt(1 + hv * hv) * peak_up
    return fin.baseline_floor_from_phase_energy(fin.peak_phase_energy(m_bins(g), phase))


def dynamic_layer(g, m26, layouts, pairs, amb, snr):
    """The lorry beside each layout for the whole pass: P2-26's surface wave with each layout, the difference of two
    layouts' imprints (both against the same run without any cavity), the oracle at the worst frequency and side, the
    scattered difference carried unattenuated to the scene's farthest corner."""
    from katabasis.compose import load_site
    from katabasis.compose.site import parse_site
    site = load_site('bench-void')
    base_model = m26.model

    def model(site_, chamber):
        if not (isinstance(chamber, dict) and 'boxes' in chamber):
            return base_model(site_, chamber)
        raw = copy.deepcopy(site_.raw)
        f0 = raw['features'][0]
        feats = []
        for k, (c, s) in enumerate(chamber['boxes']):
            f = copy.deepcopy(f0)
            f['id'] = f'part_{k}'
            f['shape'] = dict(f['shape'], centre=list(c), size=list(s))
            f['fill'] = 'air'
            feats.append(f)
        raw['features'] = feats
        return parse_site(raw, site_.directory)
    m26.model = model
    spec = m26.MODELS['main']
    los = np.asarray(g.los_enu)
    r0, (xs, ys), _ = m26.run_case(site, spec, 'line', None)
    n_pad = 2 ** int(np.ceil(np.log2(len(r0.t) * 2)))
    U0, f = m26.spectra(r0, n_pad)
    nx, ny = len(xs), len(ys)
    U0 = U0.reshape(nx, ny, 3, -1)
    keep = (f >= m26.F_LOW) & (f <= spec['f_top'])
    f = f[keep]
    U0 = U0[..., keep]
    ref = np.sqrt((np.abs(U0[:, :, 2]) ** 2).mean(axis=(0, 1)))
    H = {}
    for name, boxes in layouts.items():
        if not boxes:
            H[name] = 0.0
            continue
        r1, _, _ = m26.run_case(site, spec, 'line', {'boxes': [[list(c), list(s)] for c, s in boxes]})
        U1, _ = m26.spectra(r1, n_pad)
        Us = U1.reshape(nx, ny, 3, -1)[..., keep] - U0
        H[name] = m26.rotations(Us, los) / ref                           # [direction, x, y, f] per unit incident
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    R = np.hypot(X, Y)
    tp = np.clip((DYN_TAPER[1] - R) / (DYN_TAPER[1] - DYN_TAPER[0]), 0, 1)
    dA = float(xs[1] - xs[0]) ** 2
    k0 = 4 * np.pi / g.lam
    cells = inf.cells_per_m2(g)
    src = json.loads((SITES / 'acquisitions' / 'giza-20250827.json').read_text())['source']['shape']
    corner = 0.5 * float(np.hypot(src[0] * g.dx, src[1] * g.dr / np.sin(g.theta)))
    cul = amb['cultural']
    band = cul['band_hz']['value']
    v = cul['bus_or_truck_over_bump']['value']
    sel = (f >= band[0]) & (f <= band[1])
    rows = []
    for a, b in pairs:
        D = (H[a] if not np.isscalar(H[a]) else 0) - (H[b] if not np.isscalar(H[b]) else 0)
        D = np.asarray(D)
        if D.ndim == 0:
            continue
        near = (np.abs(D) ** 2 * (tp ** 2)[None, :, :, None]).sum(axis=(1, 2)) * dA          # [direction, f]
        sig = np.max([(np.abs(D[:, (R >= r) & (R < r + 5), :]) ** 2).sum(axis=1) * dA / 5.0
                      for r in (25.0, 30.0, 35.0)], axis=0)                                 # [direction, f], each its own
        w2 = (2 * np.pi * f) ** 2
        nb = 2 * cells * k0 ** 2 * near / w2[None, :]
        fb = 2 * cells * k0 ** 2 * sig * (corner - 32.0) / w2[None, :]
        tot = (np.sqrt(nb) + np.sqrt(fb)) ** 2
        tot_sel = np.where(sel[None, :], tot, 0)
        j, k = np.unravel_index(np.argmax(tot_sel), tot.shape)
        allow2 = (1 + GRID_ALLOWANCE) ** 2                      # P2-26's grid check, the thin roof's (the larger)
        d = 0.5 * snr * 2 * v * v * float(tot[j, k]) * allow2
        dn = 0.5 * snr * 2 * v * v * float(np.where(sel[None, :], nb, 0).max()) * allow2
        rows.append({'pair': f'{a} vs {b}', 'worst_f_hz': float(f[k]), 'side_deg': 90 * int(j), 'sides_sampled': int(tot.shape[0]),
                     'oracle_d': d,
                     'oracle_tv': fin.oracle_tv_from_mean_energy(d),
                     'oracle_growth_to_target': float(D_TARGET / np.sqrt(2 * d)) if d > 0 else None,
                     'oracle_tv_within_32m': fin.oracle_tv_from_mean_energy(dn),
                     'largest_difference_over_motion': float(np.abs(D[..., sel]).max())})
    m26.model = base_model
    return rows


def main():
    m25 = module('p225', Path(__file__).with_name('p2_25_information_bound.py'))
    m26 = module('p226', Path(__file__).with_name('p2_26_imprint_spectrum.py'))
    m28 = module('p228', Path(__file__).with_name('p2_28_two_depths.py'))
    m31 = module('p231', Path(__file__).with_name('p2_31_shape_test.py'))
    m25.TAPER_STATIC = TAPER
    g = DwellGeometry.from_record('giza-20250827')
    amb = json.loads((SITES / 'ambient.json').read_text())
    p204 = load('p2_04_chamber_imprint')
    p225 = load('p2_25_information_bound')
    mat = json.loads((SITES / 'materials.json').read_text())['materials']
    rows_ = mat if isinstance(mat, list) else [dict(v, id=k) for k, v in mat.items()]
    rock = next(m for m in rows_ if m['id'] == 'limestone-mokattam')
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    host = (val('rho') * val('vp') ** 2 - 2 * val('rho') * val('vs') ** 2, val('rho') * val('vs') ** 2)
    layouts = {**m31.BOXES, **EXTRA}
    snr = 10 ** (SNR_DB / 10) / g.band_frac
    bands = [{'case': r['case'], 'strain_amplitude': r['strain_amplitude'], 'f_hz': r['f_hz']} for r in p225['static']]
    params = {'layouts': {k: [{'centre': list(c), 'size': list(s)} for c, s in v] for k, v in layouts.items()},
              'pairs': PAIRS, 'directions_deg': list(np.rad2deg(DIRECTIONS)), 'half_width_m': HALF, 'taper_m': TAPER,
              'snr_db': SNR_DB, 'target_delta': D_TARGET, 'ambient_bands': bands,
              'depth_allowance': 'the same-depth pair bounds every allowed vertical transformation'}
    with Run(RID, 'Shape discrimination, bounded', params) as run:
        t0 = time.time()
        kern = {}
        for name, boxes in layouts.items():
            if boxes:
                kern[name], _ = m31.layout_kernels(m28, name, boxes, host)
        print(f"  static kernels: {time.time() - t0:.0f} s", flush=True)
        zero = {k: np.zeros_like(next(iter(kern.values()))[k]) for k in ('Kxx', 'Kyy', 'Kxy')}
        get = lambda n: kern[n] if n in kern else {'xy': next(iter(kern.values()))['xy'], **zero}
        floor = background_floor(p204, amb, g)
        # each layout against empty ground (the common reference), all bands; then each pair's difference for the oracle
        single = {}
        for n in sorted({x for pair in PAIRS for x in pair} - {'none'}):
            t1 = time.time()
            single[n] = memo(RID, f'single-{n}', lambda: static_layer(g, get(n), bands, m25, snr, floor), __file__,
                             version=f'static-v2-{GUARD:.6f}')
            print(f"  {n} against none: tight {single[n][0][0]['tight_ceiling']:.2e} ({time.time() - t1:.0f} s)", flush=True)
        static = []
        for a, b in PAIRS:
            t1 = time.time()
            KA, KB = get(a), get(b)
            diff = {'xy': KA['xy'], **{k: KA[k] - KB[k] for k in ('Kxx', 'Kyy', 'Kxy')}}
            rows, fb = memo(RID, f'pair-{a}-{b}', lambda: static_layer(g, diff, bands, m25, snr, floor), __file__,
                            version=f'static-v2-{GUARD:.6f}')
            for i, r in enumerate(rows):
                ra = single[a][0][i]
                rb = single[b][0][i] if b != 'none' else None
                r['tight_ceiling_difference_diagnostic'] = r.pop('tight_ceiling')
                r['certificate_tv_difference_diagnostic'] = r.pop('certificate_tv')
                r.pop('certificate_tpr_at_005')
                r['tight_ceiling'] = ra['tight_ceiling'] + (rb['tight_ceiling'] if rb else 0.0)
                r['certificate_tv'] = min(1.0, ra['certificate_tv'] + (rb['certificate_tv'] if rb else 0.0))
                r['certificate_tpr_at_005'] = min(1.0, 0.05 + r['certificate_tv'])
                r['via'] = 'TV(A, none) + TV(none, B) for tight and certificate; the difference for the oracle'
            fa = single[a][1]
            static.append({'pair': f'{a} vs {b}', 'rows': rows,
                           'discrimination_signal_factor': float(np.sqrt(fa / fb)) if fb > 0 else None})
            print(f"  {a} vs {b}: oracle TV {rows[0]['oracle_tv']:.2e}, tight {rows[0]['tight_ceiling']:.2e} under the "
                  f"microseisms ({time.time() - t1:.0f} s)", flush=True)
        t1 = time.time()
        dynamic = dynamic_layer(g, m26, layouts, PAIRS, amb, snr)
        print(f"  dynamic: {time.time() - t1:.0f} s", flush=True)
        mic = [p['rows'][0] for p in static]
        worst_dyn = max(dynamic, key=lambda r: r['oracle_tv']) if dynamic else None
        spiral = next(p for p in static if p['pair'].startswith('column_spiral'))
        lvs = next(p for p in static if p['pair'].startswith('L_tunnel vs straight'))
        finding = (
            f"Telling one layout from another is bounded as detection is, and depth ambiguity does not rescue a reader: the "
            f"same-depth bound holds for every allowed vertical transformation. Under Giza's regional microseism level the "
            f"oracle, told everything but which layout is there, tells an L-shaped from a straight tunnel with detection "
            f"minus false alarm at most {lvs['rows'][0]['oracle_tv']:.1e}, and a column with a spiral from a plain column at "
            f"most {spiral['rows'][0]['oracle_tv']:.1e}; the tight bound gives {lvs['rows'][0]['tight_ceiling']:.1e} and "
            f"{spiral['rows'][0]['tight_ceiling']:.1e} (each layout against empty ground, summed); the loose certificate with "
            f"the background, the same way, caps discrimination at "
            f"{100 * max(r['certificate_tpr_at_005'] for r in mic):.3f}% at 5% for every pair. Discrimination needs "
            f"{min(p['discrimination_signal_factor'] for p in static if not p['pair'].endswith('none')):.1f} to "
            f"{max(p['discrimination_signal_factor'] for p in static if not p['pair'].endswith('none')):.1f} "
            f"times the signal that detecting the first layout needs."
            + (f" With a lorry bouncing beside the layouts all pass, known exactly, at the worst of four sampled sides, "
               f"the oracle's bound for the hardest-to-exclude pair ({worst_dyn['pair']}) reaches {worst_dyn['oracle_tv']:.2f} at "
               f"{worst_dyn['worst_f_hz']:.0f} Hz"
               + (f", so this argument does not exclude telling them apart at 95% / 5%" if worst_dyn['oracle_growth_to_target'] < 1
                  else f", still excluding 95% / 5% until the signal grows {worst_dyn['oracle_growth_to_target']:.1f} times")
               + ": unresolved, not a demonstrated opportunity, for shape as for detection." if worst_dyn else "")
            + " Only these layouts, this rock and these sources are covered.")
        run.save({'static': static, 'single': {n: v[0] for n, v in single.items()}, 'dynamic': dynamic,
                  'background_floor': floor, 'direction_guard': float(GUARD), 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
