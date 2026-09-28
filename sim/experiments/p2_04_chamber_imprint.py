"""P2-04 · The chamber's imprint on the ground's own motion.

    uv run python experiments/p2_04_chamber_imprint.py

Whatever a satellite could use to find the chamber passively, it would have to be the difference the
chamber makes to how the ground moves. Ambient waves are far longer than the chamber and its depth
(microseisms are 15 km long, waves at 3 Hz half a kilometre), so around the chamber the rock is strained
almost uniformly and almost statically, and the imprint is a static elastic problem: open a void in rock
under a uniform horizontal stress and see how the surface moves.

1. Static kernels, by the lab's own solver (katabasis.seismic.static): the one-chamber bench (a 6 m air
   cube, its centre 15 m down, in bare Mokattam limestone, sites/bench-void) loaded in turn with a unit
   s_xx, s_yy and s_xy, run to rest. The same loading without the chamber must not move at all (it does
   not). Checked against grid spacing (1 m, 0.5 m) and domain size, and against Eshelby's void with
   Okada's point sources (an equal-volume sphere at the chamber's centre).
2. The imprint of a Rayleigh wave. A plane Rayleigh wave with vertical velocity amplitude V and phase
   speed c strains the surface by hv V / c along its path (hv its horizontal-to-vertical ratio): the
   strain does not depend on frequency. The imprint's line-of-sight velocity is omega times the kernel's
   displacement for that strain, averaged over the waves' directions (microseisms arrive from all sides).
3. The ambient levels, from sites/ambient.json and M1-01: Giza's microseisms (measured median), the
   noisiest stations on Earth (Peterson's high model), 1-3 Hz (measured), and a truck over a bump 15 m
   away (FTA), the last only roughly since at 15 Hz the waves are not long compared with the depth.
4. The claimed structures by the scaling the solver confirms (volume over depth squared).
"""
import copy
import json
import time
from pathlib import Path

import numpy as np

from katabasis.ambient.field import rayleigh_hv
from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import Run, load
from katabasis.seismic.analytic import moment_surface_displacement, rayleigh_speed, void_moment
from katabasis.seismic.arrays import snap_to_ground, surface_grid
from katabasis.seismic.elastic3d import Medium
from katabasis.seismic.static import relax
from sarsim.acquisition import DwellGeometry

SITES = Path(__file__).resolve().parents[2] / 'sites'
SITE = 'bench-void'
STRESS = 1e6                     # Pa; the problem is linear, so any value will do
DURATION = 0.5                   # s of damped relaxation
HALF = 40.0                      # m: the receiver grid's half-width
SPACING = 2.0                    # m between receivers
AZIMUTHS = np.deg2rad(np.arange(0, 180, 7.5))


def media(H, half):
    site = load_site(SITE)
    g = Grid.covering((-half, half), (-half, half), (-60.0, 4.0), H)
    full = Medium.from_model(voxelise(site, g, heterogeneity=False))
    raw = copy.deepcopy(site.raw)
    raw['features'] = []
    bg = Medium.from_model(voxelise(parse_site(raw, site.directory), g, heterogeneity=False))
    return site, g, full, bg


def kernel(med, rec, stress, H):
    t0 = time.time()
    u, hist, dt, nt = relax(med, stress, rec, DURATION, pml_width=int(round(12 / H)))
    drift = float(np.abs(hist['u'][:, :, -1] - hist['u'][:, :, -2]).max() / max(np.abs(u).max(), 1e-30))
    return u, {'runtime_s': time.time() - t0, 'steps': nt, 'dt_s': dt, 'last_step_change': drift}


def main():
    site, g, full, bg = media(1.0, 60.0)
    mat = json.loads((SITES / 'materials.json').read_text())
    rock = next(m for m in (mat['materials'] if isinstance(mat.get('materials'), list) else
                            [dict(v, id=k) for k, v in mat['materials'].items()]) if m['id'] == 'limestone-mokattam')
    val = lambda q: rock[q]['value'] if isinstance(rock[q], dict) else rock[q]
    vp, vs, rho = val('vp'), val('vs'), val('rho')
    mu, lam = rho * vs ** 2, rho * vp ** 2 - 2 * rho * vs ** 2
    nu, E = lam / (2 * (lam + mu)), mu * (3 * lam + 2 * mu) / (lam + mu)
    chamber = next(f for f in site.raw['features'] if f['id'] == 'chamber')['shape']
    depth, volume = -chamber['centre'][2], float(np.prod(chamber['size']))
    amb = json.loads((SITES / 'ambient.json').read_text())
    m1 = load('m1_01_ambient_levels')
    geo = DwellGeometry.from_record('giza-20250827')
    params = {'site': SITE, 'rock': 'limestone-mokattam', 'vp': vp, 'vs': vs, 'rho': rho, 'poisson': nu,
              'chamber': chamber, 'stress_pa': STRESS, 'relaxation_s': DURATION, 'damping_per_s': 30.0,
              'receivers': f'{2 * HALF:.0f} m square at {SPACING:.0f} m', 'line_of_sight': list(geo.los_enu)}
    with Run('p2_04_chamber_imprint', 'The chamber\'s imprint on the ground\'s own motion', params) as run:
        rec = surface_grid(g, full.solid, HALF, SPACING)
        xy = rec[:, :2]
        n = int(round(2 * HALF / SPACING)) + 1
        centre = int(np.argmin(np.hypot(xy[:, 0], xy[:, 1])))
        runs = {}
        u0, info0 = kernel(bg, rec, {'xx': STRESS, 'yy': STRESS}, 1.0)
        runs['control_no_chamber'] = {**info0, 'max_displacement_m': float(np.abs(u0).max())}
        K = {}
        for comp in ('xx', 'yy', 'xy'):
            K[comp], runs[comp] = kernel(full, rec, {comp: STRESS}, 1.0)
            print(f'  kernel s_{comp}: {runs[comp]["runtime_s"]:.0f} s', flush=True)
        iso = K['xx'] + K['yy']
        # convergence: the isotropic loading against grid spacing and domain size
        conv = [{'spacing_m': 1.0, 'half_width_m': 60.0, 'epicentre_uz_m': float(iso[centre, 2]),
                 'ring_max_uz_m': float(iso[:, 2].max())}]
        inner = np.all(np.abs(xy) <= 30.0, axis=1)                 # receivers well clear of any absorbing layer
        c_in = int(np.argmin(np.hypot(xy[inner, 0], xy[inner, 1])))
        conv[0]['ring_max_uz_m'] = float(iso[inner, 2].max())
        for H, half in ((1.0, 50.0), (0.5, 50.0), (1.0, 90.0)):
            _, gg, ff, _ = media(H, half)
            rr = snap_to_ground(gg, ff.solid, xy[inner])
            uu, info = kernel(ff, rr, {'xx': STRESS, 'yy': STRESS}, H)
            conv.append({'spacing_m': H, 'half_width_m': half, 'epicentre_uz_m': float(uu[c_in, 2]),
                         'ring_max_uz_m': float(uu[:, 2].max()), 'runtime_s': info['runtime_s']})
            print(f'  convergence H={H} half={half}: {uu[c_in, 2] * 1e6:.2f} um ({info["runtime_s"]:.0f} s)', flush=True)
        # the analytic comparison: Eshelby's sphere of the same volume, Okada's point sources
        M = void_moment(np.diag([STRESS, STRESS, 0.0]), volume, nu)
        ua = moment_surface_displacement(M, depth, xy[:, 0], xy[:, 1], lam, mu).T
        on_axis = np.abs(xy[:, 1] - xy[centre, 1]) < 1e-6
        profile = {'x_m': xy[on_axis, 0], 'solver_uz_m': iso[on_axis, 2], 'analytic_uz_m': ua[on_axis, 2],
                   'solver_ux_m': iso[on_axis, 0], 'analytic_ux_m': ua[on_axis, 0]}
        # per unit strain (plane stress at the surface)
        Kxx = E / (1 - nu ** 2) * (K['xx'] + nu * K['yy']) / STRESS
        Kyy = E / (1 - nu ** 2) * (K['yy'] + nu * K['xx']) / STRESS
        Kxy = 2 * mu * K['xy'] / STRESS
        los = np.asarray(geo.los_enu)
        def rayleigh_los(phi):
            """LOS displacement (m) per unit strain along a Rayleigh wave travelling towards azimuth phi."""
            s, c = np.sin(phi), np.cos(phi)
            u = Kxx * s * s + Kyy * c * c + Kxy * s * c
            return -u @ los
        maps = np.stack([rayleigh_los(p) for p in AZIMUTHS])        # [azimuth, receiver]
        rms_map = np.sqrt((maps ** 2).mean(axis=0))                  # the imprint averaged over directions
        peak_per_strain = float(rms_map.max())
        hv = rayleigh_hv(vp, vs)
        c_rock = rayleigh_speed(vp, vs)
        reg = amb['regional']
        cases = [
            ('Giza microseisms, measured', reg['microseism_vertical_0.1_0.3_hz']['value'], 0.2,
             reg['microseism_phase_speed']['value'], 'M1-01 median, 0.1-0.3 Hz; phase speed sites/ambient.json'),
            ("The world's noisiest stations", m1['peterson']['band_rms_um_s']['0.1-0.3 Hz']['high'] * 1e-6, 0.2,
             reg['microseism_phase_speed']['value'], 'Peterson (1993) high-noise model, 0.1-0.3 Hz'),
            ('Giza 1-3 Hz, measured', reg['vertical_1_3_hz']['value'], 2.0, c_rock,
             'M1-01 median; the rock\'s own Rayleigh speed, slower than the real waves, so the strain is generous'),
            ('A truck over a bump 15 m away', amb['cultural']['bus_or_truck_over_bump']['value'], 15.0, c_rock,
             'FTA level at 15 m; at 15 Hz the waves are 110 m long, so the static kernel is only a rough guide'),
        ]
        rows = []
        for label, V, f, c, src in cases:
            strain = hv * V / c
            v_imprint = 2 * np.pi * f * peak_per_strain * strain
            rows.append({'case': label, 'vertical_velocity_m_s': V, 'frequency_hz': f, 'phase_speed_m_s': c,
                         'strain': strain, 'imprint_los_velocity_m_s': v_imprint,
                         'imprint_over_motion': v_imprint / V, 'source': src})
        # the claimed structures, by the volume / depth^2 scaling (the solver's chamber as the reference)
        claim = json.loads((SITES / 'bench-khafre-claim' / 'site.json').read_text())
        deep = next(f for f in claim['features'] if f['id'] == 'deep-structure')['shape']
        scale = lambda V, d: (V / d ** 2) / (volume / depth ** 2)
        claimed = [{'structure': 'claimed deep structure (80 m cube at 1,220 m)',
                    'volume_m3': float(np.prod(deep['size'])), 'depth_m': -deep['centre'][2],
                    'relative_to_bench': scale(np.prod(deep['size']), -deep['centre'][2])}]
        finding = (
            f"The chamber barely touches the ground's own motion. Opening it in rock under a uniform horizontal stress "
            f"and letting the rock settle, the lab's solver finds the surface above it sinks by "
            f"{abs(iso[centre, 2]) * 1e6:.1f} um per MPa (Eshelby's void with Okada's point sources: "
            f"{abs(ua[centre, 2]) * 1e6:.1f} um; a finer grid: {abs(conv[2]['epicentre_uz_m']) * 1e6:.1f} um). A passing "
            f"Rayleigh wave strains the ground by hv V / c whatever its frequency, so the chamber's imprint on the "
            f"line-of-sight velocity is {rows[0]['imprint_over_motion']:.1e} of Giza's microseism motion: "
            f"{rows[0]['imprint_los_velocity_m_s'] * 1e6:.1e} um/s. At the noisiest stations on Earth it would be "
            f"{rows[1]['imprint_los_velocity_m_s'] * 1e6:.1e} um/s, at 1-3 Hz {rows[2]['imprint_los_velocity_m_s'] * 1e6:.1e} um/s, "
            f"and under a truck 15 m away roughly {rows[3]['imprint_los_velocity_m_s'] * 1e6:.1e} um/s.")
        run.save({'runs': runs, 'convergence': conv, 'profile': profile, 'hv': hv, 'rayleigh_speed_m_s': c_rock,
                  'kernel_grid': {'n': n, 'spacing_m': SPACING, 'half_width_m': HALF},
                  'imprint_map_los_per_strain_m': rms_map.reshape(n, n),
                  'imprint_peak_los_per_strain_m': peak_per_strain, 'cases': rows, 'claimed': claimed,
                  'finding': finding})
        np.savez_compressed(run.dir / 'kernels.npz', xy=xy, Kxx=Kxx, Kyy=Kyy, Kxy=Kxy, los=los)
        print(finding)


if __name__ == '__main__':
    main()
