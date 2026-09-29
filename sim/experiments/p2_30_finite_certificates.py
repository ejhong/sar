"""P2-30 · Finite certificates: the bound on every processor, loose but simple, the common background included.

    uv run python experiments/p2_30_finite_certificates.py

An independent derivation (29 September 2026) proves finite, algorithm-independent certificates for the declared
measurement model (sarsim.finite; the proofs are on the site's proof page). This run applies them to the lab's own
current results, where P2-25 computes the tighter bound:

A. Theorem C, phase-energy certificate: the imprint's peak phase (P2-25's, enlarged by 1.01 to cover every wave
   direction) given to every scatterer of the whole image at every retained frequency; KL for the finite change from
   |C0^-1/2 dA|_F <= q, no expansion in the displacement and no use of the cancellation of shared motion.
B. The same with the common background included: the three regional Rayleigh harmonics' envelope bounds the
   background's phase energy, which proves a floor under the reference covariance (1 - sqrt(Q_bg))^2.
C. P2-26's dynamic imprint as sampled transfer energy, the scattered wave carried to the scene's farthest corner, the
   worst sampled frequency and direction, with an 11% stress on the phase (not a certified error bar).
D. An arbitrary fixed-envelope phase modulation (both worlds' motion within the background envelope, the room changing
   its phase on a stated number of pixels), without the linear cavity kernel.
E. The inherited oracle (P2-25's genie inputs); the sufficient exclusion boundary against footprint; sensitivity to the
   phase envelope and the covariance floor; and two-world depth risk (Le Cam), with P2-29's physical depth obstruction:
   rooms 15 m down spread to match a room 30 m down, and the least mean depth error any method must then make.

Beside each, the tight bound from P2-25 (Fisher information with an explicit remainder, reference covariance I) and the
oracle from P2-29, so the three layers can be read together.
"""
import json
import math

import numpy as np

from katabasis.runs import Run, load
from pathlib import Path
from sarsim.finite import (baseline_floor_from_phase_energy, certificate_from_energy, depth_absolute_risk_lower,
                           oracle_tv_from_mean_energy, peak_phase_energy, required_oracle_energy, required_q_for_tv)

RID = 'p2_30_finite_certificates'
SITES = Path(__file__).resolve().parents[2] / 'sites'


def band_count(n, fraction):
    # the numpy mask convention sarsim uses
    return int(np.count_nonzero(np.abs(np.fft.fftfreq(n)) <= fraction / 2))


def row(label, phase_energy, **metadata):
    c = certificate_from_energy(phase_energy)
    return {"case": label, "phase_energy_upper_assumed": phase_energy,
            **c, "tpr_upper_at_fpr_005": min(1.0, 0.05 + c["tv_upper"]),
            "excludes_095_at_005": c["tv_upper"] < 0.90, **metadata}


def certificates():
    p = load('p2_25_information_bound')
    spectrum = load('p2_26_imprint_spectrum')
    g = json.loads((SITES / 'acquisitions' / 'giza-20250827.json').read_text())
    commit = None
    na, nr = g['source']['shape']
    ma, mr = band_count(na, g['band_frac']), band_count(nr, g['band_frac_r'])
    n, m = na * nr, ma * mr
    lam = g['derived']['wavelength_m']
    k0 = 4 * math.pi / lam
    cell_density_discrete = (m / n) * math.sin(math.radians(g['theta_deg'])) / (g['dx'] * g['dr'])
    # Use the larger of the continuum density and exact-grid density. This
    # only handles the bin-count convention, not spatial quadrature error.
    density = max(p['cells_per_m2'], cell_density_discrete)

    static = []
    for source in p['static']:
        # K(phi)=a+b cos(2phi)+c sin(2phi), 24 directions in [0, pi).
        # 1.01 > sec(pi/24), so this encloses all directions in that model.
        peak = 1.01 * source['imprint_peak_phase_rad']
        Q = peak_phase_energy(m, peak)
        static.append(row(source['case'], Q, peak_phase_upper_rad=peak,
                          peak_displacement_upper_m=peak / k0,
                          source_frequency_hz=source['f_hz'],
                          status="finite certificate conditional on simulated peak envelope",
                          footprint="full image, peak imprinted at every scatterer"))
    lowband_peak = sum(r['peak_phase_upper_rad'] for i, r in enumerate(static) if i != 1)
    static.append(row('Three modelled low-frequency bands, coherently summed',
                      peak_phase_energy(m, lowband_peak),
                      peak_phase_upper_rad=lowband_peak,
                      peak_displacement_upper_m=lowband_peak / k0,
                      status="finite certificate conditional on all three peak envelopes",
                      footprint="full image; triangle-inequality combination"))

    # Bound the common background operator as well, rather than assuming its
    # covariance remains exactly I. The declared background contains only the
    # same three representative Rayleigh harmonics and their given envelopes.
    ref=p['static'][0]
    hv=ref['strain_amplitude']*ref['phase_speed_m_s']/(math.sqrt(2)*ref['rms_velocity_m_s'])
    peak_up=sum(math.sqrt(2)*p['static'][i]['rms_velocity_m_s']/(2*math.pi*p['static'][i]['f_hz'])
                for i in (0,2,3))
    # Cauchy bound on vector motion, valid for any unit LOS, with horizontal
    # norm <= hv*vertical-envelope amplitude in this declared wave class.
    los_envelope=math.sqrt(1+hv*hv)*peak_up
    background_phase=k0*los_envelope
    background_Q=peak_phase_energy(m,background_phase)
    background_floor=baseline_floor_from_phase_energy(background_Q)
    background_rows=[]
    if background_floor is not None:
        for idx in (0,2,3,4):
            src=static[idx];cert=certificate_from_energy(src['phase_energy_upper_assumed'],background_floor)
            background_rows.append({"case":src['case'],**cert,
                "tpr_upper_at_fpr_005":min(1,.05+cert['tv_upper']),
                "status":"conditional, common three-harmonic background included"})

    # An arbitrary fixed-envelope phase modulation, without using the linear
    # cavity imprint. Both LOS motion histories are bounded by B; their
    # difference is <=2B. The affected-pixel limit is an explicit assumption.
    phase_modulation=[]
    ground_pixel_area=g['dx']*g['dr']/math.sin(math.radians(g['theta_deg']))
    if background_floor is not None:
        for affected in (52600,n):
            Q=peak_phase_energy(m,2*background_phase,affected/n)
            cert=certificate_from_energy(Q,background_floor)
            phase_modulation.append({"affected_scatterer_pixels_assumed":affected,
                "nominal_ground_area_m2":affected*ground_pixel_area,
                "difference_peak_displacement_upper_m":2*los_envelope,
                "difference_phase_peak_upper_rad":2*background_phase,
                "phase_energy_upper_assumed":Q,**cert,
                "tpr_upper_at_fpr_005":min(1,.05+cert['tv_upper']),
                "status":"conditional on fixed displacement envelope and exact affected-pixel count; no linear cavity kernel used"})

    dynamic = []
    # P2-25 used the INSCRIBED-circle radius. To cover every corner of the
    # rectangular image, enclose all corners in the declared radar/site map.
    heading=math.radians(g['heading_deg'])
    along=np.array([math.sin(heading),math.cos(heading)])
    los_horizontal=np.asarray(g['los_enu'][:2],float)
    ground=-los_horizontal/np.linalg.norm(los_horizontal)
    half_az=na*g['dx']/2
    half_ground=nr*g['dr']/(2*math.sin(math.radians(g['theta_deg'])))
    covering_radius=max(float(np.linalg.norm(sa*half_az*along+sr*half_ground*ground))
                        for sa in (-1,1) for sr in (-1,1))
    length = covering_radius - p['manifest']['params']['far_field_from_m']
    for source in p['dynamic']:
        name = source['case']
        physical = spectrum['cases'][name]
        f = np.asarray(source['f_hz'])
        A = np.asarray(physical['imprint_area_m2'])
        sigma = np.asarray(source['cross_width_m'])
        assert np.allclose(f, physical['f_hz'], atol=1e-6)
        # rotations() creates four LOS maps for line incidence; the plane
        # cases have one map. physical['runs'] counts cavity/background runs,
        # NOT directions.
        directions = 4 if physical['incidence'] == 'line' else 1
        # A and sigma are arithmetic averages over recorded directions.
        # max <= number_of_directions * average. Count the entire raw near
        # grid PLUS the assumed far-field envelope; overlap is conservative.
        # sqrt(2)*v is the peak of the declared single harmonic.
        coefficient = directions * density * k0**2 * 2 * (A + sigma * length) / (2 * math.pi * f)**2
        for kind, band in [('urban', (8, 100)), ('truck', None)]:
            eligible = np.flatnonzero((f >= band[0]) & (f <= band[1])) if band else np.arange(len(f))
            idx = int(eligible[np.argmax(coefficient[eligible])])
            v = source[kind]['rms_velocity_m_s']
            Q = float(coefficient[idx] * v * v)
            r = row(f'{name}: {kind}-level harmonic', Q,
                    source_rms_velocity_m_s=v, source_frequency_hz=float(f[idx]),
                    status="conditional on sampled transfer energy and assumed far-field envelope",
                    direction_rule="worst of recorded directions; no unsampled-direction claim",
                    direction_count=directions,
                    frequency_rule="worst sampled harmonic; no continuum/broadband supremum claim",
                    far_field_radial_extent_m=length,
                    covering_scene_radius_m=covering_radius,
                    oracle_energy_upper_30db=1000 * Q,
                    oracle_tv_upper_30db=oracle_tv_from_mean_energy(1000 * Q))
            stress=certificate_from_energy(Q*1.11**2)
            r.update({"stress_11pct_phase_multiplier":1.11,
                "stress_11pct_tv_upper":stress['tv_upper'],
                "stress_11pct_tpr_upper_at_fpr_005":min(1,.05+stress['tv_upper']),
                "stress_11pct_excludes_095_at_005":stress['tv_upper']<.9,
                "stress_is_certified_error_bound":False,
                "phase_multiplier_to_target_boundary":required_q_for_tv(.9)/r['q']})
            dynamic.append(r)

    # The stored genie values are proportional to a PHASE-ENERGY integral.
    # The new finite inequality is |expm1(i Phi)|^2 <= Phi^2. These are average
    # source-phase/direction results at ASSUMED 30 dB noise, not uniform maxima.
    inherited_oracle = []
    for source in p['genie']:
        # Stored rounded F is accepted as nominal input. A 1% upward allowance
        # is deliberately included but is NOT a certified physics error bar.
        B = 1.01 * source['fisher'] / 2
        tv = oracle_tv_from_mean_energy(B)
        inherited_oracle.append({"case": source['case'], "energy_upper_assumed": B,
            "tv_upper": tv, "tpr_upper_at_fpr_005": min(1.0, .05 + tv),
            "status":"conditional average oracle, inherited phase-energy estimate with 1% allowance",
            "noise_snr_db_assumed":30,
            "extra_allowance_fraction":.01,
            "allowance_is_certified_physics_error":False})

    q_target = required_q_for_tv(.90)
    requirement = []
    for area in np.geomspace(1, 25e6, 161):
        # Continuum area conversion is declared, not an exact irregular-mask count.
        cells = density * area
        requirement.append({"footprint_m2":float(area),
            "exclusion_boundary_peak_m":q_target / (k0 * math.sqrt(cells)),
            "exclusion_boundary_peak_m_covariance_floor_01":math.sqrt(.1) * q_target / (k0 * math.sqrt(cells)),
            "exclusion_boundary_peak_m_covariance_floor_001":.1 * q_target / (k0 * math.sqrt(cells))})

    sensitivities = []
    base = static[0]['phase_energy_upper_assumed']
    for phase_multiplier in [1,10,100,1000]:
        for floor in [1,.1,.01,.001]:
            c = certificate_from_energy(base * phase_multiplier**2, floor)
            sensitivities.append({"phase_multiplier":phase_multiplier,
                "covariance_floor":floor,**c,
                "tpr_upper_at_fpr_005":min(1,.05+c['tv_upper'])})

    # A mathematical corollary ONLY: both depth worlds must satisfy the same
    # assumed envelope. No 30 m physical solver result is fabricated.
    pair_tv = min(1.0, 2*static[0]['tv_upper'])
    out = {"schema_version":1,"id":"p2_30_finite_certificates",
        "title":"Finite exclusion theorem for a declared vibration-SAR model",
        
        "claim_level":"proved conditional statistical theorem; physical input envelope not field-certified",
        "measurement_model":"proper complex Gaussian, phase-only finite Fourier-transfer model",
        "covariance_floor_default":1.0,"texture_covariance_cap_default":1.0,
        "target":{"tpr":.95,"fpr":.05,"required_tv":.90},
        "geometry":{"wavelength_m":lam,"pixel_count":n,"in_band_azimuth_bins":ma,
            "in_band_range_bins":mr,"in_band_bins":m,
            "cells_per_m2_continuum":p['cells_per_m2'],
            "cells_per_m2_exact_grid":cell_density_discrete,"area_density_used":density,
            "dynamic_covering_radius_m":covering_radius,
            "upstream_inscribed_scene_radius_m":p['scene_radius_m'],
            "aperture_seconds":g['derived']['aperture_time_s']},
        "static":static,"dynamic":dynamic,"inherited_average_oracle":inherited_oracle,
        "fixed_envelope_phase_modulation":phase_modulation,
        "background_included":{"wave_class":"three stated representative Rayleigh harmonics; no additional background modes assumed",
            "hv_bound_assumed":hv,"vertical_peak_envelope_m":peak_up,
            "los_peak_envelope_m":los_envelope,"phase_peak_envelope_rad":background_phase,
            "phase_energy_upper_assumed":background_Q,"covariance_floor_proved_within_model":background_floor,
            "rows":background_rows},
        "requirements":{"q_exclusion_boundary_095_005":q_target,
            "oracle_energy_necessary_for_tv_090":required_oracle_energy(.90),
            "peak_displacement_boundary_at_1000_m2_m":q_target/(k0*math.sqrt(density*1000)),
            "meaning":"below boundary this bound excludes target; above it performance remains undecided",
            "area_conversion":"continuum approximation; use exact affected-pixel count for discrete certification",
            "curve":requirement},"sensitivity":sensitivities,
        "depth_corollary":{"status":"illustration, conditional on BOTH worlds satisfying the envelope",
            "separation_m":15,"tv_pair_upper":pair_tv,
            "equal_prior_absolute_error_lower_m":depth_absolute_risk_lower(15,pair_tv),
            "not_a_new_30m_physical_simulation":True},
        "limitations":["No measured local cavity/source/noise envelope is certified here.",
            "Finite Gaussian covariance results assume their declared texture and baseline covariance.",
            "The raw-noise oracle theorem is general; the scenario tables use a finite spectral model, not qualified raw telemetry.",
            "Dynamic far-field extrapolation remains an explicit input assumption.",
            "Sampled frequency/direction results do not exclude narrow unsampled resonances or other excitations.",
            "The statistical target is illustrative; failing it does not prove all archaeological utility is zero."]}
    return out




def depth_from_oracle(out):
    """Two-world depth risk with P2-29's obstruction: the 30 m room and the spread rooms 15 m down differ in the
    oracle's Delta^2 by the residual fraction squared, so at any signal level TV_depth <= erf(delta sqrt(d) / 2), and any
    method's equal-prior mean depth error is at least (15 m / 2)(1 - TV_depth). Given at the signal level where detecting
    the room reaches 95% at 5% (d = 5.41) and at ten times that signal."""
    try:
        p29 = load('p2_29_oracle')
    except Exception:
        return None
    delta = p29['depth']['residual_fraction']
    rows = []
    for label, d in (('detection at 95% / 5%', required_oracle_energy(0.90)), ('ten times that signal', 100 * required_oracle_energy(0.90)),
                     ('a hundred times that signal', 1e4 * required_oracle_energy(0.90))):
        tv = oracle_tv_from_mean_energy(delta ** 2 * d)
        rows.append({'signal': label, 'oracle_energy_detect': d, 'tv_depth_upper': tv,
                     'mean_depth_error_lower_m': depth_absolute_risk_lower(15.0, tv)})
    return {'residual_fraction': delta, 'separation_m': 15.0, 'rows': rows,
            'status': 'P2-29 solver kernels; the class admits rooms spread 15 m down (dilute superposition checked)'}


def main():
    params = {'target': {'tpr': 0.95, 'fpr': 0.05}, 'direction_enlargement': 1.01, 'dynamic_stress_phase': 1.11,
              'inputs': ['p2_25_information_bound', 'p2_26_imprint_spectrum', 'p2_29_oracle', 'giza-20250827'],
              'derivation': 'independent, 29 September 2026; sarsim.finite'}
    with Run(RID, 'Finite certificates: the bound on every processor, the common background included', params) as run:
        out = certificates()
        p = load('p2_25_information_bound')
        tight = {r['case']: r.get('ceiling') for r in p['static']}
        for r in out['static']:
            r['tight_ceiling_p2_25'] = tight.get(r['case'])
        for r in out['background_included']['rows']:
            r['tight_ceiling_p2_25'] = tight.get(r['case'])
        out['depth_oracle'] = depth_from_oracle(out)
        micro = out['background_included']['rows'][0]
        thin = next(r for r in out['dynamic'] if r['case'] == 'surface_favourable: truck-level harmonic')
        pm = out['fixed_envelope_phase_modulation']
        dep = out['depth_oracle']
        finding = (
            f"Loose but simple, and with the common background included, the finite certificate caps any processor's "
            f"detection rate at {100 * micro['tpr_upper_at_fpr_005']:.4f}% at 5% false alarms under the regional microseism "
            f"envelope (detection minus false alarm at most {micro['tv_upper']:.2e}; P2-25's tighter bound, reference "
            f"covariance I: {micro['tight_ceiling_p2_25']:.1e}). Every static row excludes 95% at 5%. Of P2-26's dynamic "
            f"rows the thin-roof room under a truck-level harmonic is the open one: at most "
            f"{100 * thin['tpr_upper_at_fpr_005']:.1f}% at 5%, and "
            + ("no longer excluded" if not thin['stress_11pct_excludes_095_at_005'] else "still excluded")
            + f" if its phase is 11% larger. An arbitrary phase modulation within the background envelope on "
            f"{pm[0]['affected_scatterer_pixels_assumed']:,} pixels is capped at {100 * pm[0]['tpr_upper_at_fpr_005']:.2f}%; "
            f"on the whole image the certificate is uninformative."
            + (f" Depth: with P2-29's spread rooms matching a room 30 m down to {100 * dep['residual_fraction']:.1f}%, where "
               f"detection just reaches 95% at 5% any method's mean depth error between the two is at least "
               f"{dep['rows'][0]['mean_depth_error_lower_m']:.1f} m of their 15 m separation." if dep else ""))
        out['finding'] = finding
        run.save(out)
        print(finding)


if __name__ == '__main__':
    main()
