"""P2-26 · The chamber's imprint at every frequency, and how long a room rings: the elastic solver.

    uv run python experiments/p2_26_imprint_spectrum.py

P2-04 finds the chamber's imprint on the ground's motion statically, which holds while the waves are far longer than
the chamber is deep (below a few hertz). A cavity could ring: near a resonance its imprint could be many times the
static answer. Here the lab's validated elastic solver (P1-01) finds the imprint's amplitude and phase at every
frequency from 6 to 240 Hz, in rock with no intrinsic attenuation (which could only damp a resonance), and measures
directly how long each room rings.

A. Transfer maps (records 0.3-0.5 s, a short Gaussian force pulse):
1. A surface wave from the side: a line of vertical forces across the whole model, 75 m west of the chamber's axis,
   with and without the chamber. Waves from the other three sides follow from the bench's symmetry.
2. A P wave from below: a plane of vertical forces across the whole model, 60 m down; 3. an S wave: the same, east.
4. The favourable chamber of P2-22 (a 10 m room whose roof is 5 m down), and 5. the bench room 30 m down (for P2-28),
   under the surface wave.
6. The surface wave on a smaller model at 1 m and at 0.5 m spacing, for both rooms: the grid check, and the band to
   240 Hz, where a 6 m cavity's own modes begin and where the favourable room's 5 m roof would flex.
The imprint is the difference with and without the chamber, which removes everything the two share. It is recorded at
every surface point within 40 m of the axis (25 m on the small model) and set against the incident motion:
H(x, f) = (imprint's slant-range velocity at x) / (incident velocity's rms over the recorded ground). The width of the
largest response's envelope (the maximum of |H| over the map, frequency by frequency) is reported as the ratio of its
peak frequency to its half-power band: a measure of how broad the response is, not the Q of any one mode, and not
measurable where the band is cut by the model's upper limit.

B. Ring-down (records of 4 s, 2 s at 0.5 m), because a short record cannot exclude a long-lived narrow resonance (a
60 Hz mode with Q = 600 decays by 1/e in 3.2 s):
7. The room's own modes, excited directly: a vertical force pulse on each room's ceiling, recorded on the ceiling,
   floor, wall and the ground above, beside the same force in intact rock (the control). Reported per record: the
   spectrum's peaks and their half-power widths (the modal Q, at 0.25 Hz resolution), and the energy's decay after the
   pulse (the time to fall 60 dB and the Q it implies at the record's dominant frequency).
8. The surface wave again for 4 s over each room: how fast the imprint at the ground dies once the wave has passed, and
   its spectrum at 0.25 Hz resolution.

Waves from the side at low frequency are compared with P2-04's static kernels, which they must approach. Only the
modelled rooms, excitations and windows are covered; fractured, layered or coupled structures are not.
"""
import copy
import json
import time
from pathlib import Path

import numpy as np

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run, load
from katabasis.seismic.elastic3d import Medium, Receivers, Result, Simulation, Source, _nearest
from sarsim.acquisition import DwellGeometry

RID = 'p2_26_imprint_spectrum'
SITE = 'bench-void'
BENCH = {'centre': [0, 0, -15], 'size': [6, 6, 6]}
FAVOURABLE = {'centre': [0, 0, -10], 'size': [10, 10, 10]}
DEEP = {'centre': [0, 0, -30], 'size': [6, 6, 6]}
MODELS = {
    # extent (x, y, z), spacing, absorbing layer (m), recorded half-width, line source x, plane source z, record, top
    'main': {'extent': ((-110.0, 110.0), (-80.0, 80.0), (-100.0, 3.0)), 'h': 1.0, 'pml_m': 20.0, 'half': 40.0,
             'line_x': -75.0, 'plane_z': -60.0, 'record_s': 0.5, 'f_top': 120.0},
    'small': {'extent': ((-55.0, 55.0), (-40.0, 40.0), (-40.0, 3.0)), 'h': 1.0, 'pml_m': 10.0, 'half': 25.0,
              'line_x': -38.0, 'plane_z': None, 'record_s': 0.3, 'f_top': 120.0},
    'small_fine': {'extent': ((-55.0, 55.0), (-40.0, 40.0), (-40.0, 3.0)), 'h': 0.5, 'pml_m': 10.0, 'half': 25.0,
                   'line_x': -38.0, 'plane_z': None, 'record_s': 0.3, 'f_top': 240.0},
    # ring-down: long records, receivers only where they are needed
    'long': {'extent': ((-70.0, 70.0), (-50.0, 50.0), (-70.0, 3.0)), 'h': 1.0, 'pml_m': 20.0, 'half': 20.0,
             'step_m': 2.0, 'line_x': -50.0, 'plane_z': None, 'record_s': 4.0, 'f_top': 120.0},
    'ring': {'extent': ((-60.0, 60.0), (-60.0, 60.0), (-60.0, 3.0)), 'h': 1.0, 'pml_m': 20.0, 'half': 0.0,
             'record_s': 4.0, 'f_top': 120.0},
    'ring_fine': {'extent': ((-30.0, 30.0), (-30.0, 30.0), (-32.0, 3.0)), 'h': 0.5, 'pml_m': 10.0, 'half': 0.0,
                  'record_s': 2.0, 'f_top': 240.0},
}
CASES = [
    # name, model, incidence, chamber
    ('surface', 'main', 'line', BENCH),
    ('p_below', 'main', 'plane_up', BENCH),
    ('s_below', 'main', 'plane_east', BENCH),
    ('surface_favourable', 'main', 'line', FAVOURABLE),
    ('surface_deep', 'main', 'line', DEEP),
    ('surface_small', 'small', 'line', BENCH),
    ('surface_small_fine', 'small_fine', 'line', BENCH),
    ('surface_favourable_small', 'small', 'line', FAVOURABLE),
    ('surface_favourable_small_fine', 'small_fine', 'line', FAVOURABLE),
]
LONG_CASES = [('long_bench', BENCH), ('long_favourable', FAVOURABLE)]
RING_CASES = [('ring_bench', 'ring', BENCH), ('ring_favourable', 'ring', FAVOURABLE),
              ('ring_favourable_fine', 'ring_fine', FAVOURABLE)]
F_LOW = 6.0                          # Hz: below this the models are too small for the waves (and P2-04 holds)
RECORD_EVERY = 4
CACHE = RESULTS / 'cache' / 'p2_26'


def cached(key, compute):
    """A solver run kept on disk under a hash of everything that defines it, so a run cut short resumes where it
    stopped (the records are the solver's own; nothing downstream is cached)."""
    import hashlib
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / (hashlib.sha1(json.dumps(key, sort_keys=True, default=str).encode()).hexdigest()[:16] + '.npz')
    if path.is_file():
        z = np.load(path, allow_pickle=False)
        res = Result(z['t'], z['traces'], None, None, float(z['dt']), {})
        return res, json.loads(str(z['extra'])), True
    res, extra = compute()
    np.savez(path, t=res.t, traces=res.traces.astype(np.float32), dt=res.dt, extra=json.dumps(extra, default=float))
    return res, extra, False
PASSED_S = 0.12                      # s: in the long records, the surface wave has crossed the recorded ground by then


def model(site, chamber):
    raw = copy.deepcopy(site.raw)
    if chamber is None:
        raw['features'] = []
    else:
        raw['features'][0]['shape'] = dict(raw['features'][0]['shape'], centre=chamber['centre'], size=chamber['size'])
    return parse_site(raw, site.directory)


def pulse(nt, dt, f_top):
    """A Gaussian force pulse whose spectrum falls to exp(-1/2) at f_top: flat from zero frequency, compact in time."""
    s = 1.0 / (2 * np.pi * f_top)
    t = np.arange(nt) * dt
    return np.exp(-0.5 * ((t - 6 * s) / s) ** 2).astype(np.float32)


def run_case(site, spec, incidence, chamber, verbose=True):
    key = {'what': 'case', 'site': SITE, 'spec': spec, 'incidence': incidence, 'chamber': chamber, 'every': RECORD_EVERY}

    def compute():
        res, (xs, ys), info = _run_case(site, spec, incidence, chamber, verbose)
        return res, {'xs': list(map(float, xs)), 'ys': list(map(float, ys)), 'info': info}
    res, extra, hit = cached(key, compute)
    if hit and verbose:
        print(f"    {incidence} {'chamber' if chamber else 'none'} h={spec['h']}: from the run cache", flush=True)
    return res, (np.asarray(extra['xs']), np.asarray(extra['ys'])), extra['info']


def _run_case(site, spec, incidence, chamber, verbose=True):
    g = Grid.covering(*spec['extent'], spec['h'])
    med = Medium.from_model(voxelise(model(site, chamber), g, heterogeneity=False))
    pml = int(round(spec['pml_m'] / spec['h']))
    sim = Simulation(med, pml_width=pml, f0=30.0, pml_vp=3300.0)
    nt = int(np.ceil(spec['record_s'] / sim.dt))
    w = pulse(nt + 1, sim.dt, spec['f_top'])
    h = spec['h']
    # receivers: the first solid cell's centre under every surface cell centre within the recorded half-width
    step = int(round(spec.get('step_m', spec['h']) / spec['h']))
    pick = lambda a: a[(np.abs(a) <= spec['half']) & (np.round((np.abs(a) - h / 2) / h).astype(int) % step == 0)]
    xs, ys = pick(g.x), pick(g.y)                                              # symmetric about the axis
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    top = np.argmax(med.solid, axis=2)
    ix = np.searchsorted(g.x, xs)
    iy = np.searchsorted(g.y, ys)
    ztop = g.origin[2]
    Z = ztop - top[np.ix_(ix, iy)] * h
    rec = Receivers(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1))
    sources, on_step = [], None
    if incidence == 'line':
        z0 = ztop - top[ix[len(ix) // 2], iy[len(iy) // 2]] * h                  # the first solid cell, as receivers
        sources = [Source((spec['line_x'], y, z0), w, 'force', (0.0, 0.0, -h)) for y in g.y]   # 1 N per metre of line
    else:
        i, j, k = _nearest(g, (0.0, 0.0, spec['plane_z']))
        b = (med.bz if incidence == 'plane_up' else med.bx)
        comp = 2 if incidence == 'plane_up' else 0
        layers = [(k, 0.5), (k - 1, 0.5)] if comp == 2 else [(k, 1.0)]
        dt32 = np.float32(sim.dt)

        def on_step(it, s):
            if it + 1 >= len(w):
                return
            for kk, share in layers:
                s.v[comp][:, :, kk] += dt32 * np.float32(share / h) * w[it + 1] * b[:, :, kk]
    t0 = time.time()
    res = sim.run(sources, rec, nt, record_every=RECORD_EVERY, on_step=on_step)
    info = {'runtime_s': round(time.time() - t0, 1), 'steps': nt, 'dt_s': sim.dt, 'cells': int(np.prod(g.shape)),
            'pml_cells': pml, 'points_per_s_wavelength_at_top': 1830.0 / spec['f_top'] / h}
    if verbose:
        print(f"    {incidence} {'chamber' if chamber else 'none'} h={h}: {info['runtime_s']:.0f} s", flush=True)
    return res, (xs, ys), info


def spectra(res, n_pad):
    """Complex spectra (n, 3, nf) of the recorded velocities, and the frequencies."""
    dt = res.t[1] - res.t[0]
    F = np.fft.rfft(res.traces.astype(np.float64), n=n_pad, axis=2) * dt
    return F, np.fft.rfftfreq(n_pad, dt)


def rotations(U, los):
    """Slant-range velocity maps (receding positive) of the imprint U [nx, ny, 3, nf] for the wave from the west (as run)
    and the same wave turned to arrive from the south, east and north: the bench is symmetric under quarter turns."""
    E, N, Up = U[:, :, 0], U[:, :, 1], U[:, :, 2]
    maps = []
    # a quarter turn anticlockwise takes (x, y) -> (-y, x) and (uE, uN) -> (-uN, uE); on a square symmetric grid that
    # is np.rot90 of the arrays
    for q in range(4):
        e, n, u = E, N, Up
        for _ in range(q):
            e, n = -np.rot90(n, 1, axes=(0, 1)), np.rot90(e, 1, axes=(0, 1))
            u = np.rot90(u, 1, axes=(0, 1))
        maps.append(-(los[0] * e + los[1] * n + los[2] * u))
    return np.stack(maps)                       # [4, nx, ny, nf]


def run_ring(site, spec, chamber, control=False):
    key = {'what': 'ring', 'site': SITE, 'spec': spec, 'chamber': chamber, 'control': control, 'every': RECORD_EVERY}

    def compute():
        res, names, info = _run_ring(site, spec, chamber, control)
        return res, {'names': names, 'info': info}
    res, extra, hit = cached(key, compute)
    if hit:
        print(f"    ring {'control' if control else 'room'}: from the run cache", flush=True)
    return res, extra['names'], extra['info']


def _run_ring(site, spec, chamber, control=False):
    """A vertical force pulse on the room's ceiling (or at the same point in intact rock), recorded on the ceiling,
    the floor, a wall and the ground above: the room's own modes and how long they last."""
    g = Grid.covering(*spec['extent'], spec['h'])
    med = Medium.from_model(voxelise(model(site, None if control else chamber), g, heterogeneity=False))
    pml = int(round(spec['pml_m'] / spec['h']))
    sim = Simulation(med, pml_width=pml, f0=30.0, pml_vp=3300.0)
    nt = int(np.ceil(spec['record_s'] / sim.dt))
    w = pulse(nt + 1, sim.dt, spec['f_top'])
    h = spec['h']
    c, size = np.asarray(chamber['centre'], float), np.asarray(chamber['size'], float)
    top, bottom = c[2] + size[2] / 2, c[2] - size[2] / 2
    x0 = h / 2                                                                  # the cell centre next to the axis
    points = {'ceiling': (x0, x0, top + h / 2), 'floor': (x0, x0, bottom - h / 2),
              'wall': (size[0] / 2 + h / 2, x0, c[2] + h / 2), 'ground_above': (x0, x0, -h / 2),
              'ground_10m': (10 + x0, x0, -h / 2)}
    rec = Receivers(np.array(list(points.values())))
    t0 = time.time()
    res = sim.run([Source(points['ceiling'], w, 'force', (0.0, 0.0, -1.0))], rec, nt, record_every=RECORD_EVERY)
    return res, list(points), {'runtime_s': round(time.time() - t0, 1), 'steps': nt, 'dt_s': sim.dt,
                               'cells': int(np.prod(g.shape))}


def peaks_and_widths(f, P, f_lo, f_hi, floor_db=-30.0):
    """Local maxima of a power spectrum within [f_lo, f_hi] standing within floor_db of its largest value, each with its
    half-power width (linear interpolation) and Q = f / width; a width under two bins is reported as a lower bound."""
    sel = np.where((f >= f_lo) & (f <= f_hi))[0]
    top = P[sel].max()
    if not top > 0:
        return []
    df = f[1] - f[0]
    out = []
    for i in sel[1:-1]:
        if not (P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] >= top * 10 ** (floor_db / 10)):
            continue
        half = P[i] / 2
        lo = i
        while lo > 0 and P[lo] > half:
            lo -= 1
        hi = i
        while hi < len(P) - 1 and P[hi] > half:
            hi += 1
        fl = f[lo] + (half - P[lo]) / (P[lo + 1] - P[lo]) * df if P[lo + 1] != P[lo] else f[lo]
        fh = f[hi - 1] + (P[hi - 1] - half) / (P[hi - 1] - P[hi]) * df if P[hi - 1] != P[hi] else f[hi]
        width = max(fh - fl, 1e-12)
        out.append({'f_hz': float(f[i]), 'power_db': float(10 * np.log10(P[i] / top)), 'width_hz': float(width),
                    'q': float(f[i] / width), 'resolved': bool(width >= 2 * df)})
    return out


def decay(v, dt, t_from):
    """The energy decay after t_from (Schroeder's backward integral of v^2): the time to fall 60 dB, from a line fitted
    between -5 and -35 dB, and the dominant frequency of the record after t_from."""
    e = (v[t_from:] ** 2)[::-1].cumsum()[::-1]
    if not np.isfinite(e[0]) or e[0] <= 0:
        return None, None, None
    db = 10 * np.log10(e / e[0] + 1e-30)
    tt = np.arange(len(db)) * dt
    fit = (db <= -5) & (db >= -35)
    if fit.sum() < 5:
        return None, None, None
    slope = np.polyfit(tt[fit], db[fit], 1)[0]                                  # dB per second (negative)
    t60 = -60.0 / slope
    sp = np.abs(np.fft.rfft(v[t_from:])) ** 2
    fr = np.fft.rfftfreq(len(v) - t_from, dt)
    f_dom = float(fr[1:][np.argmax(sp[1:])])
    return float(t60), f_dom, float(np.pi * f_dom * t60 / (3 * np.log(10)))


def main():
    site = load_site(SITE)
    geo = DwellGeometry.from_record('giza-20250827')
    los = np.asarray(geo.los_enu)
    p204 = load('p2_04_chamber_imprint')
    params = {'site': SITE, 'models': MODELS, 'cases': [c[:3] + (c[3],) for c in CASES], 'long_cases': LONG_CASES,
              'ring_cases': RING_CASES, 'f_low_hz': F_LOW, 'passed_s': PASSED_S,
              'rock': 'limestone-mokattam, no intrinsic attenuation', 'line_of_sight': list(los),
              'pml_vp_m_s': 3300.0, 'record_every': RECORD_EVERY}
    with Run(RID, "The chamber's imprint at every frequency, and how long a room rings", params) as run:
        out, maps, none = {}, {}, {}
        for name, mname, incidence, chamber in CASES:
            spec = MODELS[mname]
            print(f'  {name}', flush=True)
            if (mname, incidence) not in none:                  # the run without a chamber, shared by the cases
                none[(mname, incidence)] = run_case(site, spec, incidence, None)
            r0, (xs, ys), i0 = none[(mname, incidence)]
            r1, _, i1 = run_case(site, spec, incidence, chamber)
            n_pad = 2 ** int(np.ceil(np.log2(len(r0.t) * 2)))
            U0, f = spectra(r0, n_pad)
            U1, _ = spectra(r1, n_pad)
            nx, ny = len(xs), len(ys)
            U0 = U0.reshape(nx, ny, 3, -1)
            Us = U1.reshape(nx, ny, 3, -1) - U0
            keep = (f >= F_LOW) & (f <= spec['f_top'])
            f = f[keep]
            U0, Us = U0[..., keep], Us[..., keep]
            # the incident motion: its rms over the recorded ground (vertical; east for the S wave), which a node in
            # the incident field at one point cannot make small
            comp = 0 if incidence == 'plane_east' else 2
            ref = np.sqrt((np.abs(U0[:, :, comp]) ** 2).mean(axis=(0, 1)))       # [nf]
            inc_los = -np.einsum('xyc...,c->xy...', U0, los)
            H = (rotations(Us, los) if incidence == 'line' else -np.einsum('xyc...,c->xy...', Us, los)[None]) / ref
            dA = spec['h'] ** 2
            A = (np.abs(H) ** 2).sum(axis=(1, 2)).mean(axis=0) * dA                # m^2, averaged over directions
            peak = np.abs(H).max(axis=(0, 1, 2))
            X, Y = np.meshgrid(xs, ys, indexing='ij')
            R = np.hypot(X, Y)
            rel = np.linalg.norm(Us, axis=2) / np.sqrt((np.linalg.norm(U0, axis=2) ** 2).mean(axis=(0, 1)))
            k_pk = int(np.argmax(peak))
            above = peak ** 2 >= 0.5 * peak[k_pk] ** 2                            # the half-power band round the largest
            lo_i, hi_i = k_pk, k_pk
            while lo_i > 0 and above[lo_i - 1]:
                lo_i -= 1
            while hi_i < len(f) - 1 and above[hi_i + 1]:
                hi_i += 1
            width = f[hi_i] - f[lo_i]
            truncated = bool(hi_i == len(f) - 1 or lo_i == 0)                     # the band cut by the model's range
            # local maxima standing 20% above everything within 10 Hz either side: what a narrow resonance would draw
            sharp = [float(f[i]) for i in range(1, len(f) - 1)
                     if peak[i] >= peak[max(0, i - 1)] and peak[i] >= peak[min(len(f) - 1, i + 1)]
                     and peak[i] > 1.2 * np.max(np.r_[peak[(f >= f[i] - 10) & (f <= f[i] - 5)], 0.0])
                     and peak[i] > 1.2 * np.max(np.r_[peak[(f <= f[i] + 10) & (f >= f[i] + 5)], 0.0])]
            case = {'incidence': incidence, 'chamber': chamber, 'model': mname, 'runs': [i0, i1], 'f_hz': f,
                    'imprint_area_m2': A, 'peak_transfer': peak,
                    'peak_over_motion': rel.max(axis=(0, 1)),          # 3-component, over the motion's rms
                    'peak_f_hz': float(f[k_pk]), 'peak_value': float(peak[k_pk]),
                    'half_power_band_hz': [float(f[lo_i]), float(f[hi_i])],
                    'envelope_width_ratio': None if truncated or width <= 0 else float(f[k_pk] / width),
                    'band_truncated': truncated, 'sharp_peaks_hz': sharp,
                    'half_max_radius_m_at_peak': float(R[np.abs(H[..., k_pk]).max(axis=0) >= 0.5 * peak[k_pk]].max()),
                    'incident_ref_spectrum': np.abs(ref), 'incident_los_rms_over_ref': np.sqrt(
                        (np.abs(inc_los) ** 2).mean(axis=(0, 1))) / np.abs(ref)}
            out[name] = case
            maps[name] = {'x': xs, 'y': ys, 'f': f, 'H': H.astype(np.complex64)}
            print(f"    peak |H| {case['peak_value']:.3g} at {case['peak_f_hz']:.0f} Hz", flush=True)

        # B8. the surface wave for 4 s: how the imprint at the ground dies once the wave has passed
        spec = MODELS['long']
        print('  long records', flush=True)
        r0, (xs, ys), i0 = run_case(site, spec, 'line', None)
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        near = (np.hypot(X, Y) <= 10.0).ravel()
        dt = r0.t[1] - r0.t[0]
        k_pass = int(PASSED_S / dt)
        longs = {}
        for name, chamber in LONG_CASES:
            r1, _, i1 = run_case(site, spec, 'line', chamber)
            d = (r1.traces - r0.traces).astype(np.float64)[near]                  # the imprint near the axis [n, 3, nt]
            e = (d ** 2).sum(axis=(0, 1))                                          # its energy against time
            inc = (r0.traces.astype(np.float64)[near] ** 2).sum(axis=(0, 1))
            axis = int(np.argmin(np.hypot(X, Y).ravel()[near]))
            v = -np.einsum('ct,c->t', d[axis], los)                                # slant-range velocity over the axis
            t60, f_dom, q_eq = decay(v, dt, k_pass)
            P = np.abs(np.fft.rfft(v)) ** 2
            fr = np.fft.rfftfreq(len(v), dt)
            late = e[k_pass:].sum() / e.sum()
            longs[name] = {'runs': [i0, i1], 'record_s': float(r1.t[-1]),
                           'energy_after_pass_fraction': float(late),
                           'imprint_energy_db_vs_time': {'t_s': r1.t[::20], 'db': 10 * np.log10(e[::20] / e.max() + 1e-30)},
                           'incident_energy_db_vs_time': {'t_s': r0.t[::20], 'db': 10 * np.log10(inc[::20] / inc.max() + 1e-30)},
                           't60_after_pass_s': t60, 'dominant_f_hz': f_dom, 'q_equivalent': q_eq,
                           'peaks': peaks_and_widths(fr, P, F_LOW, spec['f_top'])}
            print(f"    {name}: energy after the wave {late:.2e}, t60 {t60}", flush=True)

        # B7. the rooms' own modes, excited on the ceiling, beside intact rock
        rings = {}
        for name, mname, chamber in RING_CASES:
            spec = MODELS[mname]
            print(f'  {name}', flush=True)
            res, names, info = run_ring(site, spec, chamber)
            ctl, _, info_c = run_ring(site, spec, chamber, control=True)
            dt = res.t[1] - res.t[0]
            k0 = int(0.05 / dt)                                                     # after the pulse and its direct wave
            recs = {}
            for j, nm in enumerate(names):
                vz = res.traces[j, 2].astype(np.float64)
                vc = ctl.traces[j, 2].astype(np.float64)
                P = np.abs(np.fft.rfft(vz)) ** 2
                fr = np.fft.rfftfreq(len(vz), dt)
                t60, f_dom, q_eq = decay(vz, dt, k0)
                t60c, _, _ = decay(vc, dt, k0)
                pk = peaks_and_widths(fr, P, F_LOW, spec['f_top'])
                recs[nm] = {'t60_s': t60, 'dominant_f_hz': f_dom, 'q_equivalent': q_eq, 'control_t60_s': t60c,
                            'peaks': pk, 'largest_resolved_q': max([x['q'] for x in pk if x['resolved']], default=None),
                            'unresolved_peaks': [x['f_hz'] for x in pk if not x['resolved']],
                            'spectrum_db': {'f_hz': fr[(fr >= F_LOW) & (fr <= spec['f_top'])][::4],
                                            'db': 10 * np.log10(P[(fr >= F_LOW) & (fr <= spec['f_top'])][::4]
                                                                / max(P.max(), 1e-300) + 1e-30)}}
            rings[name] = {'model': mname, 'chamber': chamber, 'runs': [info, info_c], 'records': recs,
                           'resolution_hz': float(fr[1] - fr[0])}
            qs = [r['largest_resolved_q'] for r in recs.values() if r['largest_resolved_q']]
            print(f"    largest resolved modal Q {max(qs) if qs else None}", flush=True)

        # the static limit (P2-04): a Rayleigh wave's strain hv V / c, the kernel's LOS displacement per unit strain,
        # velocity = 2 pi f x displacement; per unit vertical velocity: H_static = 2 pi f hv K / c
        hv, c_r = p204['hv'], p204['rayleigh_speed_m_s']
        k_map = np.asarray(p204['imprint_map_los_per_strain_m'])                  # rms over directions, 2 m grid
        A_static_per_f2 = (2 * np.pi * hv / c_r) ** 2 * (k_map ** 2).sum() * 4.0  # m^2 / Hz^2 (2 m cells)
        s = out['surface']
        f = np.asarray(s['f_hz'])
        low = (f >= F_LOW) & (f <= 12.0)
        static_ratio = np.sqrt(np.asarray(s['imprint_area_m2'])[low] / (A_static_per_f2 * f[low] ** 2))
        grid = {}
        for a, b in (('surface_small', 'surface_small_fine'), ('surface_favourable_small', 'surface_favourable_small_fine')):
            fa, fb = np.asarray(out[a]['f_hz']), np.asarray(out[b]['f_hz'])
            common = fa[(fa >= F_LOW) & (fa <= 120.0)]
            Aa = np.interp(common, fa, out[a]['imprint_area_m2'])
            Ab = np.interp(common, fb, out[b]['imprint_area_m2'])
            grid[a] = {'f_hz': common, 'area_ratio_fine_over_coarse': Ab / Aa,
                       'max_amplitude_change': float(np.max(np.abs(np.sqrt(Ab / Aa) - 1)))}
        np.savez_compressed(Path(run.dir) / 'maps.npz', **{f'{k}_{q}': v[q] for k, v in maps.items() for q in v})
        sb, fav = out['surface'], out['surface_favourable']
        fine = out['surface_small_fine']
        sharp = {k: v['sharp_peaks_hz'] for k, v in out.items() if v['sharp_peaks_hz']}
        w = lambda c: 'its band cut by the model\'s upper limit' if c['band_truncated'] else \
            f"envelope width ratio {c['envelope_width_ratio']:.1f}"
        ring_q = {k: max([r['largest_resolved_q'] or 0 for r in v['records'].values()]) for k, v in rings.items()}
        # decay times against intact rock's: the same pulse with no room decays through the model's own residual floor
        # (absorbing edges, grid dispersion above the band), so a room's decay time means something only beside it
        t60_ratio = [r['t60_s'] / r['control_t60_s'] for v in rings.values() for r in v['records'].values()
                     if r['t60_s'] and r['control_t60_s']]
        ctl_t60 = [r['control_t60_s'] for v in rings.values() for r in v['records'].values() if r['control_t60_s']]
        unresolved = sorted({f for v in rings.values() for r in v['records'].values() for f in r['unresolved_peaks']})
        finding = (
            (f"In these models no response peaks narrowly: no case shows a peak standing 20% above its surroundings "
             f"between {F_LOW:.0f} and {fine['f_hz'][-1]:.0f} Hz. " if not sharp else
             f"Peaks standing 20% above their surroundings: {json.dumps(sharp)}. ")
            + f"Under a surface wave the imprint's largest value is {sb['peak_value']:.3f} of the motion, at "
            f"{sb['peak_f_hz']:.0f} Hz ({w(sb)}); {fav['peak_value']:.3f} for the room under a 5 m roof ({w(fav)}); "
            f"{out['surface_deep']['peak_value']:.3f} for the bench room 30 m down. Waves from below: "
            f"{out['p_below']['peak_value']:.3f} (P, {w(out['p_below'])}) and {out['s_below']['peak_value']:.3f} "
            f"(S, {w(out['s_below'])}). Rung on its ceiling and recorded for {MODELS['ring']['record_s']:.0f} s, the "
            f"largest modal Q resolved from the spectra (0.25 Hz bins) is {max(ring_q.values()):.1f}, "
            + (f"with peaks narrower than a bin at {unresolved} Hz; " if unresolved else
               "with no peak within 30 dB of the largest narrower than a bin; ")
            + f"the energy's decay time is set by the model's residual floor, not by the rooms (intact rock with no room "
            f"{min(ctl_t60):.1f} to {max(ctl_t60):.1f} s; a room's at most {max(t60_ratio):.2f} of intact rock's), so a "
            f"long-lived mode weaker than that floor is not excluded; over {MODELS['long']['record_s']:.0f} s of a "
            f"passing surface wave, "
            f"{100 * max(v['energy_after_pass_fraction'] for v in longs.values()):.2f}% of the imprint's energy "
            f"at the ground comes after the wave has passed. At the low end the imprint approaches P2-04's static answer "
            f"({static_ratio.min():.2f} to {static_ratio.max():.2f} of it between {F_LOW:.0f} and 12 Hz). Halving the "
            f"grid changes the imprint by at most {100 * grid['surface_small']['max_amplitude_change']:.0f}% below "
            f"120 Hz ({100 * grid['surface_favourable_small']['max_amplitude_change']:.0f}% for the thin roof).")
        run.save({'cases': out, 'long': longs, 'rings': rings,
                  'static': {'imprint_area_per_hz2_m2': A_static_per_f2, 'f_hz': f[low],
                             'dynamic_over_static': static_ratio}, 'grid_check': grid, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
