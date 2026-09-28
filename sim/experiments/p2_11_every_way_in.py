"""P2-11 · Every way a buried room could reach a radar.

    uv run python experiments/p2_11_every_way_in.py

Phase 2 tested one route, the one the published method claims: vibrations in a single image. A viewer that
looks underground from any radar data needs every route closed or open, so this experiment asks the general
question. A radar records echoes, and an echo depends on two things only: how strongly each scatterer
returns the wave, and where each scatterer is while it is watched. A buried room can therefore reach the
record in four ways, and no others:

A. the wave itself reaches the room and comes back (penetration);
B. the room changes how the ground above it moves while one image is taken (vibration, Phase 2's subject);
C. the room changes how the ground deforms between images, under loads that vary: tides, air pressure,
   the year's heat, the water table (interferometry across passes);
D. the room changes the surface's own material: its temperature or moisture.

Each is bounded here generously to the room, for three rooms in open Giza limestone: the bench chamber (6 m,
centre 15 m down, the shallowest and most favourable case in this lab), a large hall (20 m, centre 30 m),
and the claimed deep structure (an 80 m cube at 1,220 m). Where a number comes from outside the lab its
source is listed in SOURCES and recorded with the result.

A: measured penetration into the most transparent ground on Earth, hyperarid sand, is taken as an upper
   bound for limestone, which conducts 50 to 200 times more (Davis and Annan). Two-way loss to depth d is
   8.686 d / delta_p dB for a one-way power penetration depth delta_p. A lower bound for limestone at any
   band comes from its conductivity alone: sigma sqrt(mu0 / eps) nepers of power per metre.
B: the room's imprint on the ground's velocity (P2-04), under the ground's own trembling and under the
   strongest shaking there is, against the smallest motions ever measured from orbit (the 2026 papers) and
   the lab's detection floors (P2-05).
C: the room's static imprint under each load, from the analytic Eshelby void with Okada point sources (P2-04
   found this 1.4 times the solver's, so it favours the room), against interferometric precision, and the
   number of passes a detection would need.
D: the surface temperature anomaly the room makes by blocking the year's heat wave.
"""
import numpy as np

from katabasis.runs import Run, load
from katabasis.seismic.analytic import moment_surface_displacement, void_moment
from sarsim.acquisition import DwellGeometry

RID = 'p2_11_every_way_in'

ROOMS = [
    {'id': 'bench', 'name': 'the bench chamber', 'size_m': 6.0, 'centre_m': 15.0},
    {'id': 'hall', 'name': 'a large hall', 'size_m': 20.0, 'centre_m': 30.0},
    {'id': 'claim', 'name': 'the claimed deep structure', 'size_m': 80.0, 'centre_m': 1220.0},
]

# A: penetration. One-way power penetration depth in hyperarid sand, the most transparent ground measured.
BANDS = [
    {'id': 'X', 'name': 'X band (ICEYE, TerraSAR-X)', 'wavelength_m': 0.031, 'sand_penetration_m': 0.3,
     'source': 'martone2014', 'note': 'about 20 to 30 cm into dry desert soil; the larger value is used'},
    {'id': 'L', 'name': 'L band (ALOS-2, NISAR)', 'wavelength_m': 0.24, 'sand_penetration_m': 2.0,
     'source': 'mccauley1982', 'note': 'buried river valleys seen under 1 to 2 m of sand; the larger value is used'},
    {'id': 'P', 'name': 'P band (Biomass)', 'wavelength_m': 0.69, 'sand_penetration_m': 5.0,
     'source': 'esa2025', 'note': 'up to 5 m below the desert surface; used as the penetration depth itself, '
                                   'which is more than the observation implies'},
]
LOSS_BUDGET_DB = 30.0                  # a generous dynamic range: bright surface to the image's noise floor
LIMESTONE = {'sigma_s_m': 0.5e-3, 'eps_r': 8.0}     # the least conductive, most permittive end of Davis and Annan
BIOMASS_BANDWIDTH_HZ = 6e6             # the 432-438 MHz allocation: slant resolution c / 2B

# B: vibration in one image
EARTHQUAKE_F_HZ = 2.0
STRONGEST_PGV_M_S = 3.0                # among the largest ground velocities ever recorded, near a fault
SMALLEST_MEASURED_M_S = 0.66e-3        # Vattulainen 2026, test 7: RMS radial velocity of the smallest tested motion
SMALLEST_MEASURED_ERROR_M_S = 0.91e-3  # its RMS error

# C: loads between passes (generous amplitudes)
TIDAL_STRAIN = 5e-8                    # horizontal, both directions: the solid-earth tide's largest
AIR_PRESSURE_PA = 3000.0               # a 30 hPa swing
WATER_TABLE_PA = 1.0e4                 # a 1 m swing of the water table, as if the room sat below it
YEAR_T_K = 12.0                        # annual amplitude of the ground surface's temperature
ALPHA_PER_K = 8e-6                     # limestone's linear thermal expansion (upper end)
KAPPA_M2_S = 1.3e-6                    # limestone's thermal diffusivity (upper end: deeper heat wave)
YEAR_S = 3.156e7
INSAR_MM = 0.5                         # one measurement's precision, a well-behaved persistent scatterer
SIGMA_K = 3.0

SOURCES = {
    'martone2014': 'Martone, Brautigam, Rizzoli, Yague-Martinez and Krieger (2014), Enhancing interferometric SAR '
                   'performance over sandy areas: experience from the TanDEM-X mission, IEEE JSTARS',
    'mccauley1982': 'McCauley et al. (1982), Subsurface valleys and geoarcheology of the eastern Sahara revealed by '
                    'Shuttle radar, Science 218, 1004-1020; Paillou et al., ALOS/PALSAR over the Sahara',
    'esa2025': 'ESA (23 June 2025), Biomass satellite returns striking first images of forests and more',
    'davis1989': 'Davis and Annan (1989), Ground-penetrating radar for high-resolution mapping of soil and rock '
                 'stratigraphy, Geophysical Prospecting 37, 531-551 (table of typical properties)',
    'vattulainen2026': 'Vattulainen et al. (2026), Assessment of spaceborne SAR micro-motion measurement for '
                       'vibration-based SHM, IEEE Access 14, 6043-6064, doi:10.1109/ACCESS.2026.3652346',
    'lotti2026': 'Lotti et al. (2026), Monitoring bridge vibrations via spaceborne SAR micro-Doppler, Structural '
                 'Control and Health Monitoring, doi:10.1155/stc/3858095',
}


def rock():
    p = load('p2_04_chamber_imprint')['manifest']['params']
    rho, vp, vs = p['rho'], p['vp'], p['vs']
    mu = rho * vs ** 2
    lam = rho * vp ** 2 - 2 * mu
    nu = lam / (2 * (lam + mu))
    return {'lam': lam, 'mu': mu, 'nu': nu, 'E': 2 * mu * (1 + nu), 'M': lam + 2 * mu}


def imprint_los(stress, room, r, los):
    """Peak line-of-sight surface displacement (m) of a void of the room's volume in a remote stress (Pa)."""
    V = room['size_m'] ** 3
    M = void_moment(stress, V, r['nu'])
    half = max(60.0, 2.5 * room['centre_m'])
    x = np.linspace(-half, half, 161)
    X, Y = np.meshgrid(x, x, indexing='ij')
    u = moment_surface_displacement(M, room['centre_m'], X, Y, r['lam'], r['mu'])
    d = np.tensordot(np.asarray(los), u, axes=1)
    return float(np.abs(d).max())


def loads(r):
    """The loads between passes as remote stress tensors (Pa) at depth z (m), in the (x, y, up) frame."""
    E, nu = r['E'], r['nu']
    delta = np.sqrt(2 * KAPPA_M2_S / (2 * np.pi / YEAR_S))
    tide = lambda z: np.diag([1, 1, 0]) * E * TIDAL_STRAIN / (1 - nu)
    air = lambda z: -AIR_PRESSURE_PA * np.diag([nu / (1 - nu), nu / (1 - nu), 1.0])
    water = lambda z: -WATER_TABLE_PA * np.eye(3)
    heat = lambda z: -np.diag([1, 1, 0]) * E * ALPHA_PER_K * YEAR_T_K * np.exp(-z / delta) / (1 - nu)
    return delta, {'tide': (tide, 'solid-earth tide, horizontal strain 5e-8'),
                   'air': (air, 'air pressure, a 30 hPa swing'),
                   'water': (water, 'the water table, a 1 m swing around the room'),
                   'heat': (heat, "the year's heat wave, 12 K at the surface, taken at the room's roof")}


def main():
    g = DwellGeometry.from_record('giza-20250827')
    r = rock()
    imp = load('p2_04_chamber_imprint')
    bud = load('p2_05_budget')
    vs = load('p2_09_virtual_sensors')
    params = {'rooms': ROOMS, 'bands': BANDS, 'loss_budget_db': LOSS_BUDGET_DB, 'limestone': LIMESTONE,
              'earthquake_f_hz': EARTHQUAKE_F_HZ, 'strongest_pgv_m_s': STRONGEST_PGV_M_S,
              'smallest_measured_m_s': SMALLEST_MEASURED_M_S,
              'loads': {'tidal_strain': TIDAL_STRAIN, 'air_pa': AIR_PRESSURE_PA, 'water_pa': WATER_TABLE_PA,
                        'year_t_k': YEAR_T_K, 'alpha_per_k': ALPHA_PER_K, 'kappa_m2_s': KAPPA_M2_S},
              'insar_mm': INSAR_MM, 'sigma_k': SIGMA_K, 'geometry': 'giza-20250827'}
    with Run(RID, 'Every way a buried room could reach a radar', params) as run:
        los = g.los_enu
        # ---- A: the wave reaches the room
        eta = 376.73 / np.sqrt(LIMESTONE['eps_r'])
        cond_db_per_m = 2 * 10 * np.log10(np.e) * LIMESTONE['sigma_s_m'] * eta      # two-way, power
        A = []
        for room in ROOMS:
            roof = room['centre_m'] - room['size_m'] / 2
            row = {'room': room['id'], 'roof_m': roof, 'limestone_floor_db': cond_db_per_m * roof, 'bands': {}}
            for b in BANDS:
                loss = 10 * np.log10(np.e) * 2 * roof / b['sand_penetration_m']
                row['bands'][b['id']] = {'two_way_loss_db_if_sand': loss,
                                         'reach_m': LOSS_BUDGET_DB / (10 * np.log10(np.e) * 2) * b['sand_penetration_m'],
                                         'short_by_db': loss - LOSS_BUDGET_DB}
            A.append(row)
        biomass_res = 2.998e8 / (2 * BIOMASS_BANDWIDTH_HZ)
        # ---- B: vibration within one image
        micro = imp['cases'][0]
        hf = imp['cases'][2]                                     # 1-3 Hz: the ratio at 2 Hz on the rock's Rayleigh speed
        claimed_rel = imp['claimed'][0]['relative_to_bench']
        lab = [{'target': row['target'], 'floor_m_s': row['a_min_m_s'], 'shortfall': row['shortfall'],
                'shake_needed_m_s': row['shake_needed_m_s']}
               for row in bud['rows'] if row['case'] == micro['case']]
        B = {
            'trembling_m_s': micro['vertical_velocity_m_s'], 'imprint_m_s': micro['imprint_los_velocity_m_s'],
            'imprint_over_motion_microseisms': micro['imprint_over_motion'],
            'imprint_over_motion_2hz': hf['imprint_over_motion'],
            'smallest_measured_m_s': SMALLEST_MEASURED_M_S, 'smallest_measured_error_m_s': SMALLEST_MEASURED_ERROR_M_S,
            'measured_over_imprint': SMALLEST_MEASURED_M_S / micro['imprint_los_velocity_m_s'],
            'measured_over_trembling': SMALLEST_MEASURED_M_S / micro['vertical_velocity_m_s'],
            'shaking_needed_bench_m_s': SMALLEST_MEASURED_M_S / hf['imprint_over_motion'],
            'shaking_needed_claim_m_s': SMALLEST_MEASURED_M_S / (hf['imprint_over_motion'] * claimed_rel),
            'strongest_pgv_m_s': STRONGEST_PGV_M_S,
            'imprint_under_strongest_m_s': STRONGEST_PGV_M_S * hf['imprint_over_motion'],
            'lab_floors': lab,
            'virtual_sensor_gain_complex': vs['gains']['complex'],
        }
        # ---- C: deformation between passes
        delta, L = loads(r)
        C = {'heat_wave_depth_m': delta, 'insar_mm': INSAR_MM, 'rooms': {}}
        # the ground's own seasonal breathing: uniform uplift of the free surface under the year's heat wave
        breathing = ALPHA_PER_K * (1 + r['nu']) / (1 - r['nu']) * YEAR_T_K * delta / np.sqrt(2)
        C['surface_breathing_m'] = breathing
        for room in ROOMS:
            roof = max(room['centre_m'] - room['size_m'] / 2, 0.5)
            rows = {}
            for key, (fn, label) in L.items():
                z = roof if key == 'heat' else room['centre_m']
                u = imprint_los(fn(z), room, r, los)
                n_needed = 2 * (SIGMA_K * INSAR_MM * 1e-3 / max(u, 1e-30)) ** 2
                rows[key] = {'label': label, 'imprint_m': u, 'passes_needed': n_needed,
                             'years_daily': n_needed / 365.25, 'breathing_over_imprint': breathing / max(u, 1e-30)}
            C['rooms'][room['id']] = rows
        # ---- D: the surface's temperature
        D = {}
        for room in ROOMS:
            a = room['size_m'] / 2
            dT = YEAR_T_K * np.exp(-2 * room['centre_m'] / delta) * (a / room['centre_m']) ** 2
            D[room['id']] = {'surface_anomaly_k': float(dT)}
        # ---- the answer, route by route, for the bench (the most favourable room)
        bench_C = C['rooms']['bench']
        best_C = max(bench_C.values(), key=lambda v: v['imprint_m'])
        finding = (
            f"No route carries a Giza room to an X-band radar, and none carries the claimed structure to any radar. "
            f"A: the wave gets about {BANDS[0]['sand_penetration_m'] * 100:.0f} cm into the driest sand at X band, "
            f"{BANDS[1]['sand_penetration_m']:.0f} m at L band and {BANDS[2]['sand_penetration_m']:.0f} m at P band; "
            f"limestone is less transparent, so the bench chamber's roof at {A[0]['roof_m']:.0f} m lies at least "
            f"{A[0]['bands']['X']['two_way_loss_db_if_sand']:.0f} dB down at X band and the claimed structure "
            f"{A[2]['bands']['X']['two_way_loss_db_if_sand']:,.0f} dB, against about {LOSS_BUDGET_DB:.0f} dB an image can span; "
            f"only P band could reach a shallow room, blurred into its {biomass_res:.0f} m resolution cell with the surface. "
            f"B: the smallest motion ever measured from orbit, on a corner reflector, is "
            f"{B['measured_over_imprint']:.0e} times the bench chamber's imprint and {B['measured_over_trembling']:,.0f} times "
            f"the ground's whole trembling; lifting the imprint to it needs shaking of {B['shaking_needed_bench_m_s']:.1f} m/s "
            f"during the 25 s pass, with reflectors installed over the spot. C: between passes the largest load, "
            f"{best_C['label'].split(',')[0]}, moves the ground over the bench chamber by {best_C['imprint_m'] * 1e6:.2f} um; "
            f"at {INSAR_MM} mm a measurement, a detection needs {best_C['passes_needed']:.0e} passes, and the ground's own "
            f"seasonal breathing ({breathing * 1e3:.2f} mm) is {best_C['breathing_over_imprint']:,.0f} times larger. D: the "
            f"chamber changes the surface's temperature by {D['bench']['surface_anomaly_k'] * 1e3:.2f} mK. Any processing, "
            f"gates or learning included, is a function of the record and cannot add what these routes do not carry.")
        run.save({'sources': SOURCES, 'rooms': ROOMS, 'A_penetration': {'rows': A, 'bands': BANDS,
                  'limestone_two_way_db_per_m_lower_bound': cond_db_per_m, 'biomass_slant_resolution_m': biomass_res},
                  'B_vibration': B, 'C_between_passes': C, 'D_surface': D, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
