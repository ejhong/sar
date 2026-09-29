"""P2-26 · The chamber's imprint at every frequency: waves from the side and from below, by the elastic solver.

    uv run python experiments/p2_26_imprint_spectrum.py

P2-04 finds the chamber's imprint on the ground's motion statically, which holds while the waves are far longer than
the chamber is deep (below a few hertz). A cavity could ring. If it had a resonance, its imprint near that frequency
could be many times the static answer, and that is the one physical route by which the imprint could grow. Here the
lab's validated elastic solver (P1-01) finds the imprint at every frequency from 6 to 240 Hz, in rock with no
intrinsic attenuation, which could only damp a resonance:

1. A surface wave from the side: a line of vertical forces across the whole model, 75 m west of the chamber's axis,
   driven by a short Gaussian pulse (its spectrum flat to within a factor of 3 up to 120 Hz), with and without the
   chamber. Waves from the other three sides follow from the bench's symmetry.
2. A P wave from below: a plane of vertical forces across the whole model, 60 m down.
3. An S wave from below: the same plane pushing east.
4. The favourable chamber of P2-22 (a 10 m room whose roof is 5 m down), under the surface wave.
5. The surface wave on a smaller model at 1 m and at 0.5 m spacing, for both rooms: the grid check, and the band to
   240 Hz, where a 6 m cavity's own modes begin (a shear wavelength of twice its size near 150 Hz) and where the
   favourable room's 5 m roof would flex.

The imprint is the difference between the runs with and without the chamber, which removes everything the two share:
the source, the boundaries' small reflections, the grid's dispersion. It is recorded at every surface point within
40 m of the axis (25 m on the small model) and set against the incident motion:
H(x, f) = (imprint's slant-range velocity at x) / (incident velocity's rms over the recorded ground), per frequency. P2-25 needs its
integral over the ground, A(f) = the integral of |H|^2 dA, and its maps.

Waves from the side at low frequency are compared with P2-04's static kernels, which they must approach.
"""
import copy
import json
import time
from pathlib import Path

import numpy as np

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run, load
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, _nearest
from sarsim.acquisition import DwellGeometry

RID = 'p2_26_imprint_spectrum'
SITE = 'bench-void'
BENCH = {'centre': [0, 0, -15], 'size': [6, 6, 6]}
FAVOURABLE = {'centre': [0, 0, -10], 'size': [10, 10, 10]}
MODELS = {
    # extent (x, y, z), spacing, absorbing layer (m), recorded half-width, line source x, plane source z, record, top
    'main': {'extent': ((-110.0, 110.0), (-80.0, 80.0), (-100.0, 3.0)), 'h': 1.0, 'pml_m': 20.0, 'half': 40.0,
             'line_x': -75.0, 'plane_z': -60.0, 'record_s': 0.5, 'f_top': 120.0},
    'small': {'extent': ((-55.0, 55.0), (-40.0, 40.0), (-40.0, 3.0)), 'h': 1.0, 'pml_m': 10.0, 'half': 25.0,
              'line_x': -38.0, 'plane_z': None, 'record_s': 0.3, 'f_top': 120.0},
    'small_fine': {'extent': ((-55.0, 55.0), (-40.0, 40.0), (-40.0, 3.0)), 'h': 0.5, 'pml_m': 10.0, 'half': 25.0,
                   'line_x': -38.0, 'plane_z': None, 'record_s': 0.3, 'f_top': 240.0},
}
CASES = [
    # name, model, incidence, chamber
    ('surface', 'main', 'line', BENCH),
    ('p_below', 'main', 'plane_up', BENCH),
    ('s_below', 'main', 'plane_east', BENCH),
    ('surface_favourable', 'main', 'line', FAVOURABLE),
    ('surface_small', 'small', 'line', BENCH),
    ('surface_small_fine', 'small_fine', 'line', BENCH),
    ('surface_favourable_small', 'small', 'line', FAVOURABLE),
    ('surface_favourable_small_fine', 'small_fine', 'line', FAVOURABLE),
]
F_LOW = 6.0                          # Hz: below this the models are too small for the waves (and P2-04 holds)
RECORD_EVERY = 4


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
    g = Grid.covering(*spec['extent'], spec['h'])
    med = Medium.from_model(voxelise(model(site, chamber), g, heterogeneity=False))
    pml = int(round(spec['pml_m'] / spec['h']))
    sim = Simulation(med, pml_width=pml, f0=30.0, pml_vp=3300.0)
    nt = int(np.ceil(spec['record_s'] / sim.dt))
    w = pulse(nt + 1, sim.dt, spec['f_top'])
    h = spec['h']
    # receivers: the first solid cell's centre under every surface cell centre within the recorded half-width
    xs = g.x[np.abs(g.x) <= spec['half']]
    ys = g.y[np.abs(g.y) <= spec['half']]
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


def main():
    site = load_site(SITE)
    geo = DwellGeometry.from_record('giza-20250827')
    los = np.asarray(geo.los_enu)
    p204 = load('p2_04_chamber_imprint')
    params = {'site': SITE, 'models': MODELS, 'cases': [c[:3] + (c[3],) for c in CASES], 'f_low_hz': F_LOW,
              'rock': 'limestone-mokattam, no intrinsic attenuation', 'line_of_sight': list(los),
              'pml_vp_m_s': 3300.0, 'record_every': RECORD_EVERY}
    with Run(RID, "The chamber's imprint at every frequency", params) as run:
        out, maps = {}, {}
        for name, mname, incidence, chamber in CASES:
            spec = MODELS[mname]
            print(f'  {name}', flush=True)
            r0, (xs, ys), i0 = run_case(site, spec, incidence, None)
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
            # local maxima standing 20% above everything within 10 Hz either side: what a resonance would draw
            sharp = [float(f[i]) for i in range(1, len(f) - 1)
                     if peak[i] >= peak[max(0, i - 1)] and peak[i] >= peak[min(len(f) - 1, i + 1)]
                     and peak[i] > 1.2 * np.max(np.r_[peak[(f >= f[i] - 10) & (f <= f[i] - 5)], 0.0])
                     and peak[i] > 1.2 * np.max(np.r_[peak[(f <= f[i] + 10) & (f >= f[i] + 5)], 0.0])]
            case = {'incidence': incidence, 'chamber': chamber, 'model': mname, 'runs': [i0, i1], 'f_hz': f,
                    'imprint_area_m2': A, 'peak_transfer': peak,
                    'peak_over_motion': rel.max(axis=(0, 1)),          # 3-component, over the motion's rms
                    'peak_f_hz': float(f[k_pk]), 'peak_value': float(peak[k_pk]),
                    'half_power_band_hz': [float(f[lo_i]), float(f[hi_i])],
                    'q_of_largest': float(f[k_pk] / width) if width > 0 else None,
                    'band_edge_largest': bool(k_pk == len(f) - 1 or k_pk == 0), 'sharp_peaks_hz': sharp,
                    'half_max_radius_m_at_peak': float(R[np.abs(H[..., k_pk]).max(axis=0) >= 0.5 * peak[k_pk]].max()),
                    'incident_ref_spectrum': np.abs(ref), 'incident_los_rms_over_ref': np.sqrt(
                        (np.abs(inc_los) ** 2).mean(axis=(0, 1))) / np.abs(ref)}
            out[name] = case
            maps[name] = {'x': xs, 'y': ys, 'f': f, 'H': H.astype(np.complex64)}
            print(f"    peak |H| {case['peak_value']:.3g} at {case['peak_f_hz']:.0f} Hz", flush=True)

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
        ab = np.asarray(sb['peak_transfer'])
        fine = out['surface_small_fine']
        sharp = {k: v['sharp_peaks_hz'] for k, v in out.items() if v['sharp_peaks_hz']}
        q = lambda c: 'at the band edge' if c['band_edge_largest'] else f"Q {c['q_of_largest']:.1f}"
        finding = (
            (f"The chamber does not ring: no case shows a peak standing 20% above its surroundings between {F_LOW:.0f} "
             f"and {fine['f_hz'][-1]:.0f} Hz. " if not sharp else
             f"Peaks standing 20% above their surroundings: {json.dumps(sharp)}. ")
            + f"Under a surface wave the imprint's largest value is {sb['peak_value']:.3f} of the motion over the axis, "
            f"at {sb['peak_f_hz']:.0f} Hz ({q(sb)}); {fav['peak_value']:.3f} for the favourable room ({q(fav)}). "
            f"Waves from below: at most {out['p_below']['peak_value']:.3f} (P, {q(out['p_below'])}) and "
            f"{out['s_below']['peak_value']:.3f} (S, {q(out['s_below'])}). At the low end the imprint approaches "
            f"P2-04's static answer ({static_ratio.min():.2f} to {static_ratio.max():.2f} of it between {F_LOW:.0f} and "
            f"12 Hz). Halving the grid changes the imprint by at most "
            f"{100 * grid['surface_small']['max_amplitude_change']:.0f}% below 120 Hz "
            f"({100 * grid['surface_favourable_small']['max_amplitude_change']:.0f}% for the favourable room).")
        run.save({'cases': out, 'static': {'imprint_area_per_hz2_m2': A_static_per_f2, 'f_hz': f[low],
                                          'dynamic_over_static': static_ratio}, 'grid_check': grid,
                  'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
