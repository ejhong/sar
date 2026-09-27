"""P1-01 · The elastic solver against exact solutions, at publication size.

    uv run python experiments/p1_01_forward_validation.py

1. Stokes: a point force in a whole space; three receivers, all components.
2. Rayleigh: the surface wave's speed against the Rayleigh equation, at 8-20
   points per minimum S wavelength.
3. Reciprocity across a void.
4. Absorbing boundaries: energy left in the grid after the waves have gone.
"""
import time

import numpy as np
from scipy.signal import correlate

from katabasis.compose.grid import Grid
from katabasis.runs import Run
from katabasis.seismic.analytic import rayleigh_speed, stokes_displacement
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source, ricker

VP, VS, RHO = 3300.0, 1830.0, 2400.0          # Mokattam limestone (sites/materials.json)


def medium(g, air=None):
    air = np.zeros(g.shape, bool) if air is None else air
    return Medium.from_arrays(g, np.where(air, 0, VP), np.where(air, 0, VS), np.where(air, 0, RHO), air=air)


def ncc(a, b):
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def stokes():
    h, f0 = 2.0, 30.0
    g = Grid.covering((-90, 90), (-90, 90), (-90, 90), h)
    sim = Simulation(medium(g), pml_width=12, faces=(True,) * 6, f0=f0)
    T = 0.13
    nt = int(T / sim.dt)
    t = np.arange(nt) * sim.dt
    w = ricker(t, f0) * 1e9
    src = np.array([1.0, 1.0, 1.0])
    recs = np.array([[1.0, 1.0, 61.0], [61.0, 1.0, 1.0], [41.0, 35.0, 29.0]])
    labels = ['above the force', 'broadside', 'oblique']
    t0 = time.time()
    res = sim.run([Source(tuple(src), w, 'force', (0, 0, 1))], Receivers(recs), nt)
    wall = time.time() - t0
    panels = []
    worst = {'ncc': 1.0, 'amp': 0.0}
    for n, r in enumerate(recs):
        u = stokes_displacement(t + sim.dt / 2, w, r - src, 2, VP, VS, RHO)
        v = np.gradient(u, t, axis=1)
        for c, comp in enumerate(('east', 'north', 'up')):
            if np.abs(v[c]).max() < 1e-2 * np.abs(v).max():
                continue
            cc = ncc(v[c], res.traces[n, c])
            amp = float(np.abs(res.traces[n, c]).max() / np.abs(v[c]).max())
            worst['ncc'] = min(worst['ncc'], cc)
            worst['amp'] = max(worst['amp'], abs(amp - 1))
            if comp in ('up', 'east') and len(panels) < 4 and not (labels[n] == 'broadside' and comp == 'east'):
                step = max(1, nt // 220)
                scale = float(np.abs(v[c]).max())
                panels.append({'receiver': labels[n], 'component': comp, 'offset_m': float(np.linalg.norm(r - src)),
                               't_ms': (t[::step] * 1e3).round(2), 'exact': (v[c][::step] / scale).round(4),
                               'simulated': (res.traces[n, c][::step] / scale).round(4), 'ncc': cc, 'amplitude_ratio': amp})
    ppw = VS / (2.5 * f0) / h
    return {'grid': list(g.shape), 'spacing_m': h, 'f0_hz': f0, 'points_per_wavelength': ppw, 'steps': nt,
            'wall_s': wall, 'mcell_steps_per_s': np.prod(g.shape) * nt / wall / 1e6,
            'worst_ncc': worst['ncc'], 'worst_amplitude_error': worst['amp'], 'panels': panels}


def rayleigh():
    rows = []
    cr = rayleigh_speed(VP, VS)
    h = 2.0
    for ppw in (8, 10, 12, 15, 20):
        fmax = VS / (ppw * h)
        f0 = fmax / 2.5
        lam = cr / f0
        g = Grid.covering((-30, 9 * lam + 30), (-40, 40), (-1.5 * lam, 6), h)
        air = np.broadcast_to(g.z[None, None, :] > 0, g.shape).copy()
        sim = Simulation(medium(g, air), pml_width=14, f0=f0)
        nt = int((9 * lam / cr + 3 / f0) / sim.dt)
        t = np.arange(nt) * sim.dt
        k0 = int(np.argmax(sim.m.solid[0, 0]))
        zs = g.z[k0]
        xs = np.round(np.arange(3, 9) * lam / h) * h + g.origin[0] % h
        res = sim.run([Source((0.0, 0.0, zs - h / 2), ricker(t, f0) * 1e9, 'force', (0, 0, -1))],
                      Receivers(np.array([[x, 0.0, zs] for x in xs])), nt)
        tr = res.traces[:, 2].astype(float)
        up = 8
        fine = lambda s: np.interp(np.arange(0, nt, 1 / up), np.arange(nt), s)
        lags = [(np.argmax(correlate(fine(tr[a]), fine(tr[0]))) - (nt * up - 1)) / up * sim.dt for a in range(1, len(xs))]
        c = 1 / np.polyfit(xs[1:] - xs[0], lags, 1)[0]
        rows.append({'points_per_wavelength': ppw, 'measured_m_s': c, 'exact_m_s': cr, 'error_percent': 100 * (c / cr - 1)})
        print(f'  rayleigh ppw {ppw}: {100 * (c / cr - 1):+.2f}%')
    return {'rows': rows, 'exact_m_s': cr}


def reciprocity():
    h = 2.0
    g = Grid.covering((-60, 60), (-60, 60), (-60, 60), h)
    air = np.zeros(g.shape, bool)
    X, Y, Z = np.meshgrid(g.x, g.y, g.z, indexing='ij')
    air[(np.abs(X) < 8) & (np.abs(Y) < 8) & (np.abs(Z) < 8)] = True
    sim = Simulation(medium(g, air), pml_width=12, faces=(True,) * 6, f0=30)
    nt = int(0.12 / sim.dt)
    t = np.arange(nt) * sim.dt
    w = ricker(t, 30.0) * 1e9
    A, B = (-25.0, -9.0, -13.0), (27.0, 11.0, 9.0)
    r1 = sim.run([Source(A, w, 'force', (1, 0, 0))], Receivers(np.array([B])), nt)
    sim.reset()
    r2 = sim.run([Source(B, w, 'force', (0, 0, 1))], Receivers(np.array([A])), nt)
    a, b = r1.traces[0, 2], r2.traces[0, 0]
    step = max(1, nt // 220)
    s = float(np.abs(a).max())
    return {'ncc': ncc(a, b), 'amplitude_ratio': float(np.abs(a).max() / np.abs(b).max()),
            't_ms': (t[::step] * 1e3).round(2), 'forward': (a[::step] / s).round(4), 'reverse': (b[::step] / s).round(4),
            'void_m': 16.0}


def absorbing():
    h, f0 = 2.0, 30.0
    g = Grid.covering((-60, 60), (-60, 60), (-60, 60), h)
    sim = Simulation(medium(g), pml_width=12, faces=(True,) * 6, f0=f0)
    t = np.arange(4000) * sim.dt
    w = ricker(t, f0) * 1e9
    energy, times = [], []
    chunk = 20
    for n in range(0, 700, chunk):
        sim.run([Source((1.0, 1.0, 1.0), w[n:], 'force', (0, 0, 1))], Receivers(np.array([[1.0, 1.0, 1.0]])), chunk)
        energy.append(sim.energy())
        times.append((n + chunk) * sim.dt * 1e3)
    e = np.array(energy)
    return {'t_ms': np.array(times).round(2), 'energy_relative': (e / e.max()).round(8),
            'final_relative': float(e[-1] / e.max()), 'pml_cells': 12}


def main():
    with Run('p1_01_forward_validation', 'The elastic solver against exact solutions',
             {'vp_m_s': VP, 'vs_m_s': VS, 'rho_kg_m3': RHO, 'material': 'limestone-mokattam'}) as run:
        print('stokes'); st = stokes()
        print('rayleigh'); ra = rayleigh()
        print('reciprocity'); rc = reciprocity()
        print('absorbing'); ab = absorbing()
        run.save({'stokes': st, 'rayleigh': ra, 'reciprocity': rc, 'absorbing': ab,
                  'finding': (f"Every component of the point-force response matches the exact solution with correlation "
                              f"≥ {st['worst_ncc']:.3f} and amplitude within {100 * st['worst_amplitude_error']:.1f}%; the Rayleigh "
                              f"wave is within {abs(ra['rows'][1]['error_percent']):.1f}% of its exact speed at ten points per wavelength; "
                              f"reciprocity across a void holds to correlation {rc['ncc']:.4f}.")})


if __name__ == '__main__':
    main()
