"""Fixed-design, coupled physical benchmarks for motion and subsurface inference.

No real satellite observations are used. Read results/validation/design.json for
the assumptions and seeds; it is written before simulations begin. Evaluators
receive observations and a frozen candidate library, never the evaluation depth.
"""
from pathlib import Path
from dataclasses import asdict
import csv
import hashlib
import json
import platform
import sys
import time
import numpy as np
import scipy
from scipy.signal import periodogram

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sarsim import Geometry, SubapBank, viz
from sarsim.waves import WaveModel, surface_wave, sample_wave
from sarsim.measurement import ReflectorAcquisition, WaveDictionary, centred, wilson
from sarsim.pipeline import trajectories
from sarsim.tomo import default_depths, focus_windows, focus_branch_b

OUT = ROOT/'results'/'validation'
CACHE = ROOT/'results'/'cache'/'physical_validation'
DESIGN = {
    "version": 2, "seed": 20260920, "status": "synthetic controlled benchmark",
    "dictionary_depths_m": [0]+list(range(40, 261, 20)),
    "evaluation_depths_m": [70, 110, 150, 190, 230],
    "dictionary_speeds_m_s": [1800, 2000, 2200],
    "dictionary_spacing_m": 5., "evaluation_spacing_m": 2.5,
    "reference_frequency_hz": 25., "radius_m": 20.,
    "surface_amplitudes_um": [30., 100., 1000.],
    "amplitude_definition": "peak cavity-free SH displacement at the sampled receivers; shared force scaling, never target normalization",
    "receiver_snr_db": 20., "calibration_null_trials": 300, "test_null_trials": 100,
    "trials_per_depth": 20, "nominal_false_alarm_rate": .01,
    "depth_success_tolerance_m": 20.,
    "receiver_indices": [0, 3, 7, 10, 13, 17, 20],
    "source": {"type": "Ricker body force", "x_m": -150., "z_m": 10., "start_s": -.45},
    "assumptions": [
        "2D homogeneous SH wave equation; horizontal infinite cylindrical cavity, not vertical columns",
        "known active source timing and location; no intrinsic attenuation or heterogeneous geology",
        "seven known, separated coherent point reflectors; noiseless stationary-phase focusing",
        "independent 20 dB complex noise per reflector and aperture sample; no distributed clutter, atmosphere or autofocus",
        "phase readout assumes known range supports; it is an optimistic measurement, not the Biondi processor",
        "physical template inverse shares a wave-equation family with the forward model but uses a different mesh and disjoint depths",
        "noise trials measure conditional receiver uncertainty; they are not independent sites or geological realizations"
    ],
    "stress_cases": ["wrong wave speed", "wrong cavity radius", "wrong source position"],
    "frequency_depth_matrix": {"frequencies_hz": [8., 16., 25.], "depths_m": [50., 90., 130.],
                               "reason": "all test depths lie inside the fixed steering search interval"},
    "snr_extension": {"levels_db": [20., 40., 60.], "reference_amplitude_um": 1000.,
                      "reason": "Documented follow-up to low recovery in the initial 20 dB run; separate seeds and null calibration at each level. No formally preregistered field experiment is claimed."},
}


def write_json(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def write_csv(name, rows):
    with (OUT/name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def waveform(acq, depth=0, spacing=5., speed=2000., frequency=25., radius=20., source_x=-150.):
    model = WaveModel(depth=depth or None, spacing=spacing, speed=speed,
                      frequency=frequency, radius=radius, source_x=source_x)
    t, x, u = surface_wave(model, cache=CACHE)
    indices = DESIGN['receiver_indices']
    return x[indices], sample_wave(t, u[:, indices], acq.time)


def observed(acq, wave, count, seed, snr_db=20.):
    rng = np.random.default_rng(seed)
    return np.stack([acq.read_phase(acq.histories(wave, snr_db, rng)) for _ in range(count)])


def build_dictionary(acq, speeds):
    waves, depths, velocity = [], [], []
    for speed in speeds:
        for depth in DESIGN['dictionary_depths_m']:
            _, u = waveform(acq, depth, speed=speed)
            waves.append(u)
            depths.append(depth)
            velocity.append(speed)
    return WaveDictionary(waves, depths, velocity)


def motion_benchmark(acq):
    """Known motion through focusing/readout, with separate threshold noise."""
    rng = np.random.default_rng(DESIGN['seed']+1)
    t = acq.time
    freqs = np.fft.rfftfreq(len(t), np.median(np.diff(t)))
    band = (freqs >= .5) & (freqs <= 8)
    window = np.hanning(len(t))

    def measure(amplitude_um, frequency, trials):
        rows = []
        for _ in range(trials):
            phase = rng.uniform(-np.pi, np.pi)
            truth = amplitude_um*1e-6*np.sin(2*np.pi*frequency*t+phase)
            signal = acq.histories(truth[:, None], 20., rng)
            recovered = acq.read_phase(signal)[:, 0]
            power = np.abs(np.fft.rfft((recovered-recovered.mean())*window))[band]**2
            intensity = np.abs(signal[:, 0])**2
            ip = np.abs(np.fft.rfft((intensity-intensity.mean())*window))[band]**2
            rows.append((power.max(), ip.max(), freqs[band][power.argmax()]))
        return np.array(rows)

    null = measure(0, 2, 500)
    thresholds = np.quantile(null[:, :2], .99, axis=0, method='higher')
    rows = []
    for frequency in [1., 2., 4.]:
        for amplitude in [0., 30., 100., 300., 1000.]:
            values = measure(amplitude, frequency, 120)
            hits = int((values[:, 0] > thresholds[0]).sum())
            rows.append({"frequency_hz": frequency, "ground_amplitude_um": amplitude,
                         "los_amplitude_um": amplitude*acq.projection,
                         "trials": len(values), "phase_detection_rate": hits/len(values),
                         "phase_ci_low": wilson(hits, len(values))[0],
                         "phase_ci_high": wilson(hits, len(values))[1],
                         "intensity_detection_rate": float((values[:, 1] > thresholds[1]).mean()),
                         "detected_and_frequency_correct_rate": float(((values[:, 0] > thresholds[0]) & (np.abs(values[:, 2]-frequency) < .3)).mean())})
    write_csv('motion_trials.csv', rows)
    # An actual focused SLC round trip, tested separately from the many noise trials.
    ground = .001*np.sin(2*np.pi*2*t)[:, None]
    channels = acq.histories(ground, np.inf, rng)
    slc = acq.focus(channels, [0.])
    recovered = acq.read_phase(acq.separate(slc, [0.]))
    error = float(np.max(np.abs(recovered-centred(ground))))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for i, frequency in enumerate([1., 2., 4.]):
        group = [r for r in rows if r['frequency_hz'] == frequency and r['ground_amplitude_um'] > 0]
        axes[0].semilogx([r['los_amplitude_um'] for r in group], [r['phase_detection_rate'] for r in group], 'o-', color=viz.SERIES[i], label=f'{frequency:g} Hz')
    axes[0].axhline(.01, color=viz.SERIES[7], ls=':', label='Nominal 1% false alarm')
    axes[0].set(xlabel='Line-of-sight amplitude (µm)', ylabel='Detection fraction', ylim=(-.04, 1.04), title='An optimistic coherent phase receiver')
    axes[0].set_xticks([10, 30, 100, 300], labels=['10', '30', '100', '300'])
    axes[0].minorticks_off()
    axes[0].legend(fontsize=8)
    axes[1].plot(t, ground[:, 0]*1e6, color=viz.SERIES[0], lw=2, label='Physical input')
    axes[1].plot(t[::8], recovered[::8, 0]*1e6, '.', color=viz.SERIES[1], ms=3, label='Recovered through focused SLC')
    axes[1].set(xlabel='Aperture time (s)', ylabel='Ground displacement (µm)', title='Known-motion positive control')
    axes[1].legend(fontsize=8, loc='lower left', frameon=True, facecolor='white', framealpha=.95)
    fig.savefig(OUT/'motion_validation.png', dpi=160)
    plt.close(fig)
    return {"rows": rows, "roundtrip_max_error_m": error, "duration_s": float(t[-1]-t[0]),
            "samples": len(t), "projection_factor": acq.projection,
            "note": "Noise trials use the algebraically identical separated-channel readout; the focused-SLC round trip is explicitly verified. Favorable known reflectors, not a demonstrated ground-scatterer extraction method."}


def depth_benchmark(acq):
    nominal = build_dictionary(acq, [2000.])
    flexible = build_dictionary(acq, DESIGN['dictionary_speeds_m_s'])
    positions, coarse_null = waveform(acq)
    _, fine_null = waveform(acq, spacing=2.5)
    reference_peak = np.max(np.abs(coarse_null))
    dictionaries = {'known_speed': nominal, 'speed_search': flexible}
    all_rows, groups, examples = [], [], {}
    families = [
        ('matched_assumptions', {}),
        ('speed_1800', {'speed': 1800.}),
        ('speed_2200', {'speed': 2200.}),
        ('radius_30', {'radius': 30.}),
        ('source_shift_30m', {'source_x': -120.}),
    ]
    for amp_index, amplitude in enumerate(DESIGN['surface_amplitudes_um']):
        scale = amplitude*1e-6/reference_peak
        null_calibration = observed(acq, coarse_null*scale, DESIGN['calibration_null_trials'], 1000+amp_index)
        thresholds = {name: float(np.quantile(est.predict(null_calibration)['improvement'], .99, method='higher'))
                      for name, est in dictionaries.items()}
        for family_index, (family, params) in enumerate(families):
            # Stress cases at the strongest prescribed displacement only.
            if family_index and amplitude != DESIGN['surface_amplitudes_um'][-1]:
                continue
            for depth_index, depth in enumerate([0]+DESIGN['evaluation_depths_m']):
                _, physical = waveform(acq, depth, spacing=2.5, **params)
                physical *= scale
                trials = DESIGN['test_null_trials'] if depth == 0 else DESIGN['trials_per_depth']
                seed = 10000+amp_index*1000+family_index*100+depth_index
                measurement = observed(acq, physical, trials, seed)
                for estimator_name, estimator in dictionaries.items():
                    prediction = estimator.predict(measurement)
                    for trial in range(trials):
                        accepted = bool(prediction['improvement'][trial] > thresholds[estimator_name])
                        estimate = float(prediction['depth'][trial])
                        all_rows.append({"family": family, "amplitude_um": amplitude, "estimator": estimator_name,
                                         "true_depth_m": depth, "trial": trial, "seed": seed,
                                         "estimated_depth_m": estimate, "accepted": accepted,
                                         "estimated_speed_m_s": float(prediction['speed'][trial]),
                                         "score": float(prediction['improvement'][trial]),
                                         "calibrated_threshold": thresholds[estimator_name],
                                         "within_20m": bool(depth > 0 and accepted and abs(estimate-depth) <= 20)})
                if family == 'matched_assumptions' and amplitude == 1000.:
                    pred = nominal.predict(centred(physical))
                    examples[str(depth)] = {"oracle_estimate_m": float(pred['depth'][0]),
                                           "phase_estimates_m": nominal.predict(measurement)['depth'].tolist()}
                    if depth == 150:
                        # Full image measurement: carries the independently generated cavity motion.
                        channels = acq.histories(physical, 20., np.random.default_rng(seed))
                        slc = acq.focus(channels, positions)
                        readback = acq.read_phase(acq.separate(slc, positions))
                        direct = acq.read_phase(channels)
                        examples['focused_slc_roundtrip_error_m'] = float(np.max(np.abs(readback-direct)))
                        fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
                        axes[0].imshow(20*np.log10(np.abs(slc)/np.abs(slc).max()+1e-6), vmin=-60, vmax=0, cmap='gray', origin='lower', aspect='auto')
                        axes[0].set(title='Focused complex image: seven resolved reflectors', xlabel='Range pixel', ylabel='Azimuth pixel')
                        centre = 3
                        axes[1].plot(acq.time, physical[:, centre]*1e6, color=viz.SERIES[0], lw=2, label='Cavity-induced surface motion')
                        axes[1].plot(acq.time[::3], readback[::3, centre]*1e6, '.', alpha=.4, ms=2, color=viz.SERIES[1], label='Noisy SLC phase readout')
                        axes[1].set(xlabel='Aperture time (s)', ylabel='Ground displacement (µm)', title='The coupled measurement at one reflector')
                        axes[1].legend(fontsize=8)
                        fig.savefig(OUT/'coupled_measurement.png', dpi=160)
                        plt.close(fig)
            print(f'  depth benchmark: {family}, {amplitude:g} µm', flush=True)
    write_csv('depth_trials.csv', all_rows)
    for family, amplitude, estimator in sorted(set((r['family'], r['amplitude_um'], r['estimator']) for r in all_rows)):
        selected = [r for r in all_rows if (r['family'], r['amplitude_um'], r['estimator']) == (family, amplitude, estimator)]
        positive = [r for r in selected if r['true_depth_m'] > 0]
        negative = [r for r in selected if r['true_depth_m'] == 0]
        accepted = [r for r in positive if r['accepted']]
        false_alarms = sum(r['accepted'] for r in negative)
        groups.append({"family": family, "amplitude_um": amplitude, "estimator": estimator,
                       "positive_trials": len(positive), "null_trials": len(negative),
                       "detection_rate": len(accepted)/len(positive),
                       "detection_ci": wilson(len(accepted), len(positive)),
                       "false_alarm_rate": false_alarms/len(negative),
                       "false_alarm_ci": wilson(false_alarms, len(negative)),
                       "detected_within_20m_rate": sum(r['within_20m'] for r in positive)/len(positive),
                       "median_depth_error_detected_m": float(np.median([abs(r['estimated_depth_m']-r['true_depth_m']) for r in accepted])) if accepted else None,
                       "median_depth_error_all_m": float(np.median([abs(r['estimated_depth_m']-r['true_depth_m']) for r in positive]))})
    write_json('depth_summary.json', groups)
    _, coarse = waveform(acq, 150)
    _, fine = waveform(acq, 150, spacing=2.5)
    convergence = {"total_waveform_correlation": float(np.corrcoef(coarse.ravel(), fine.ravel())[0, 1]),
                   "cavity_perturbation_correlation": float(np.corrcoef((coarse-coarse_null).ravel(), (fine-fine_null).ravel())[0, 1]),
                   "relative_perturbation_norm_coarse": float(np.linalg.norm(coarse-coarse_null)/np.linalg.norm(coarse_null)),
                   "relative_perturbation_norm_fine": float(np.linalg.norm(fine-fine_null)/np.linalg.norm(fine_null))}
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    colors = [viz.SERIES[0], viz.SERIES[1], viz.SERIES[2]]
    for ax, amplitude, color in zip(axes, DESIGN['surface_amplitudes_um'], colors):
        selected = [r for r in all_rows if r['family'] == 'matched_assumptions' and r['amplitude_um'] == amplitude and r['estimator'] == 'known_speed' and r['true_depth_m'] > 0]
        for depth in DESIGN['evaluation_depths_m']:
            group = [r for r in selected if r['true_depth_m'] == depth]
            for accepted in [False, True]:
                data = [r for r in group if r['accepted'] == accepted]
                ax.scatter([depth+(r['trial']-9.5)*.65 for r in data], [r['estimated_depth_m'] for r in data], marker='o' if accepted else 'x', s=16, alpha=.55, color=color if accepted else '#9a9a9a')
        ax.plot([0, 270], [0, 270], '--', color='#777', lw=1)
        ax.set(xlim=(35, 270), ylim=(-8, 280), xlabel='Held-out true depth (m)', title=f'{amplitude:g} µm reference motion')
    axes[0].set_ylabel('Physical-template estimate (m)')
    fig.suptitle('Off-grid targets, finer forward mesh · coloured = accepted; grey × = rejected')
    fig.savefig(OUT/'depth_recovery.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    labels = ['Matched', 'Speed 1800', 'Speed 2200', 'Radius 30 m', 'Source shifted']
    for j, estimator in enumerate(dictionaries):
        selected = [next(r for r in groups if r['family'] == family and r['amplitude_um'] == 1000. and r['estimator'] == estimator) for family, _ in families]
        xx = np.arange(len(selected))+(j-.5)*.34
        axes[0].bar(xx, [r['detected_within_20m_rate'] for r in selected], width=.32, color=viz.SERIES[j], label='Fixed speed' if j == 0 else 'Speed searched')
        axes[1].bar(xx, [r['false_alarm_rate'] for r in selected], width=.32, color=viz.SERIES[j])
    for ax in axes:
        ax.set_xticks(np.arange(5), labels, rotation=25, ha='right')
        ax.set_ylim(0, 1.05)
    axes[0].set(title='Detected and localized within 20 m', ylabel='Fraction of cavity trials')
    axes[1].set(title='Cavity-free false alarms', ylabel='Fraction of null trials')
    axes[1].axhline(.01, color='#777', ls=':')
    axes[0].legend(fontsize=8)
    fig.suptitle('Assumptions matter · strongest prescribed motion, independent receiver noise')
    fig.savefig(OUT/'assumption_stress.png', dpi=160)
    plt.close(fig)
    return {"groups": groups, "convergence": convergence, "oracle_and_roundtrip": examples,
            "reference_peak_arbitrary_units": float(reference_peak),
            "threshold_note": "Thresholds use coarse-mesh cavity-free calibration only; evaluation nulls use the fine mesh. A failed held-out false-alarm check is reported, not retuned away."}


def frequency_depth_benchmark(acq):
    """Change actual depth and actual excitation independently in the wave solver."""
    bank = SubapBank(K=50, sub_frac=.5, delta_hz=88.)
    kz = acq.geometry.kz(bank.bands(acq.geometry, acq.shape[0])['nu_c'], .48)
    z = default_depths(kz)
    rows = []
    for frequency in [8., 16., 25.]:
        _, reference = waveform(acq, frequency=frequency)
        scale = .001/np.abs(reference).max()
        for depth in DESIGN['frequency_depth_matrix']['depths_m']:
            positions, wave = waveform(acq, depth, spacing=2.5, frequency=frequency)
            signal = acq.histories(wave*scale, 20., np.random.default_rng(int(30000+frequency*100+depth)))
            slc = acq.focus(signal, positions)
            rr, cc = acq.pixel_locations(positions)
            # No common motion is subtracted: known ideal reflectors, favorable zero nuisance.
            q = trajectories(slc, acq.geometry, bank, np.rint(rr).astype(int), np.rint(cc).astype(int), upsample=100, verbose=False)
            score, _ = focus_windows(q, kz, z)
            selected, _, gates = focus_branch_b(q, kz, z)
            for target in range(len(positions)):
                passed = bool(gates.full[target].any())
                rows.append({"frequency_hz": frequency, "true_depth_m": depth, "reflector": target,
                             "ungated_model_peak_m": float(z[score[target].argmax()]),
                             "ungated_score": float(score[target].max()), "accepted": passed,
                             "accepted_window_count": int(gates.full[target].sum()),
                             "gated_model_peak_m": float(z[np.nanargmax(selected[target])]) if passed else None})
            print(f'  depth/frequency: {frequency:g} Hz, {depth:g} m', flush=True)
    write_csv('frequency_depth_trials.csv', rows)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
    for i, frequency in enumerate([8., 16., 25.]):
        group = [r for r in rows if r['frequency_hz'] == frequency]
        for r in group:
            axes[0].scatter(r['true_depth_m']+(r['reflector']-3)*2.5, r['ungated_model_peak_m'], marker='o' if r['accepted'] else 'x', color=viz.SERIES[i], s=28, alpha=.65)
        axes[0].plot([], [], 'o', color=viz.SERIES[i], label=f'{frequency:g} Hz source')
    axes[0].plot([0, 260], [0, 260], '--', color='#777', lw=1)
    axes[0].set(xlabel='True cavity depth (m)', ylabel='Ungated model-axis peak (m)', title='Same physical depth; different source frequencies', ylim=(0, max(260, float(z[-1]))))
    axes[0].legend(fontsize=8)
    for i, depth in enumerate(DESIGN['frequency_depth_matrix']['depths_m']):
        counts = [sum(r['accepted'] for r in rows if r['frequency_hz'] == f and r['true_depth_m'] == depth) for f in [8., 16., 25.]]
        axes[1].plot([8, 16, 25], counts, 'o-', color=viz.SERIES[i], label=f'{depth:g} m cavity')
    axes[1].set(xlabel='Source frequency (Hz)', ylabel='Reflectors with ≥1 accepted window', title='The complete ellipse filter, applied separately', ylim=(-.3, 7.3))
    axes[1].legend(fontsize=8)
    fig.savefig(OUT/'frequency_depth.png', dpi=160)
    plt.close(fig)
    return {"trials": len(rows), "accepted_reflectors": sum(r['accepted'] for r in rows),
            "scene_count": 9, "rows": rows, "wavelength_m": .48,
            "model_depth_max_m": float(z[-1]),
            "scope": "Physical cavities feed the same ideal reflector acquisition. Depth-axis outputs use our broad-overlap patch registration with known zero nuisance, not a claim to reproduce the original method. Crosses are rejected diagnostic peaks, not detections."}


def noise_sweep(acq):
    nominal = build_dictionary(acq, [2000.])
    flexible = build_dictionary(acq, DESIGN['dictionary_speeds_m_s'])
    _, reference = waveform(acq)
    scale = .001/np.abs(reference).max()
    rows, trial_rows = [], []
    for snr in DESIGN['snr_extension']['levels_db']:
        calibration = observed(acq, reference*scale, 300, 60000+int(snr), snr)
        estimators = {'known_speed': nominal, 'speed_search': flexible}
        thresholds = {key: float(np.quantile(est.predict(calibration)['improvement'], .99, method='higher')) for key, est in estimators.items()}
        for family, params in [('matched_assumptions', {}), ('speed_1800', {'speed': 1800.}), ('source_shift_30m', {'source_x': -120.})]:
            for d in [0]+DESIGN['evaluation_depths_m']:
                _, wave = waveform(acq, d, spacing=2.5, **params)
                seed = 70000+int(snr)*100+int(d)+({'matched_assumptions': 0, 'speed_1800': 20000, 'source_shift_30m': 40000}[family])
                measurement = observed(acq, wave*scale, 100 if d == 0 else 20, seed, snr)
                for name, estimator in estimators.items():
                    prediction = estimator.predict(measurement)
                    for i, estimated in enumerate(prediction['depth']):
                        accepted = bool(prediction['improvement'][i] > thresholds[name])
                        trial_rows.append({"family": family, "snr_db": snr, "estimator": name, "true_depth_m": d,
                                           "estimated_depth_m": float(estimated), "trial": i, "seed": seed,
                                           "accepted": accepted, "within_20m": bool(d > 0 and accepted and abs(estimated-d) <= 20),
                                           "score": float(prediction['improvement'][i]), "threshold": thresholds[name]})
            for name in estimators:
                subset = [r for r in trial_rows if r['family'] == family and r['snr_db'] == snr and r['estimator'] == name]
                null = [r for r in subset if r['true_depth_m'] == 0]
                pos = [r for r in subset if r['true_depth_m'] > 0]
                hits = sum(r['accepted'] for r in pos)
                fp = sum(r['accepted'] for r in null)
                rows.append({"family": family, "snr_db": snr, "estimator": name,
                             "positive_trials": len(pos), "null_trials": len(null),
                             "detection_rate": hits/len(pos), "false_alarm_rate": fp/len(null),
                             "detection_ci": wilson(hits, len(pos)), "false_alarm_ci": wilson(fp, len(null)),
                             "detected_within_20m_rate": sum(r['within_20m'] for r in pos)/len(pos),
                             "median_depth_error_all_m": float(np.median([abs(r['estimated_depth_m']-r['true_depth_m']) for r in pos]))})
        print(f'  noise sweep: {snr:g} dB', flush=True)
    write_csv('snr_trials.csv', trial_rows)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
    for i, (family, label) in enumerate([('matched_assumptions', 'Matched assumptions'), ('speed_1800', 'Wrong wave speed'), ('source_shift_30m', 'Wrong source position')]):
        selected = [r for r in rows if r['family'] == family and r['estimator'] == 'known_speed']
        axes[0].plot([r['snr_db'] for r in selected], [r['detected_within_20m_rate'] for r in selected], 'o-', color=viz.SERIES[i], label=label)
        axes[1].plot([r['snr_db'] for r in selected], [r['false_alarm_rate'] for r in selected], 'o-', color=viz.SERIES[i])
    for ax in axes:
        ax.set(xlabel='Per-sample coherent-reflector SNR (dB)', xticks=[20, 40, 60], ylim=(-.04, 1.04))
    axes[0].set(ylabel='Cavity trials detected within 20 m', title='Where depth information survives measurement')
    axes[1].set(ylabel='Cavity-free trials falsely accepted', title='Stronger signals also expose model errors')
    axes[0].legend(fontsize=8)
    axes[1].axhline(.01, ls=':', color='#777')
    fig.savefig(OUT/'noise_transition.png', dpi=160)
    plt.close(fig)
    return {"rows": rows, "reference_ground_amplitude_um": 1000.,
            "scope": "Follow-up SNR sweep; each noise threshold uses separate null trials before evaluating the finer-grid targets. These SNRs are supplied assumptions, not measured satellite performance."}


def ambiguity_benchmark():
    """Can source freedom imitate a cavity at a single frequency?

    Fit only seven receiver positions. Report errors at fourteen additional
    positions as well as at the fitted points. Complex source weights are fully
    unconstrained, so this is a conditional counterexample, not a plausible
    explanation of the Khafre data or a universal impossibility result.
    """
    fit_indices = np.array(DESIGN['receiver_indices'])
    holdout = np.array([i for i in range(21) if i not in fit_indices])
    sources = [-150., -100., -50., 0., 50., 100., 150.]
    rows, example = [], None
    for frequency in [8., 25.]:
        def transfer(depth=None, source=-150.):
            t, x, wave = surface_wave(WaveModel(depth=depth, frequency=frequency, source_x=source), CACHE)
            harmonic = np.exp(-2j*np.pi*frequency*t)
            return x, harmonic @ wave*np.median(np.diff(t))
        x, target = transfer(150.)
        library = np.column_stack([transfer(source=source)[1] for source in sources])
        for number, indices in [(1, [0]), (3, [0, 3, 6]), (7, list(range(7)))]:
            H = library[:, indices]
            weights = np.linalg.lstsq(H[fit_indices], target[fit_indices], rcond=1e-10)[0]
            prediction = H @ weights
            for split, selected in [('fitted_receivers', fit_indices), ('held_out_receivers', holdout)]:
                error = np.linalg.norm(prediction[selected]-target[selected])/np.linalg.norm(target[selected])
                rows.append({"frequency_hz": frequency, "source_count": number, "split": split,
                             "relative_complex_misfit": float(error),
                             "maximum_source_weight": float(np.abs(weights).max()),
                             "source_weight_l2": float(np.linalg.norm(weights)),
                             "condition_number": float(np.linalg.cond(H[fit_indices]))})
            if frequency == 25. and number == 7:
                example = (x, target, prediction, fit_indices)
    write_csv('ambiguity_trials.csv', rows)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout='constrained')
    for i, frequency in enumerate([8., 25.]):
        for split, linestyle in [('fitted_receivers', '-'), ('held_out_receivers', '--')]:
            group = [r for r in rows if r['frequency_hz'] == frequency and r['split'] == split]
            axes[0].semilogy([r['source_count'] for r in group], [max(1e-15, r['relative_complex_misfit']) for r in group], 'o'+linestyle, color=viz.SERIES[i], label=f'{frequency:g} Hz, '+('fitted' if split.startswith('fitted') else 'held out'))
    axes[0].set(xlabel='Independent complex surface-source weights', ylabel='Relative complex-wavefield mismatch', title='No-cavity alternatives with flexible forcing', xticks=[1, 3, 7])
    axes[0].legend(fontsize=8)
    x, target, prediction, fitted = example
    scale = np.abs(target).max()
    axes[1].plot(x, target.real/scale, color=viz.SERIES[0], label='Cavity + one source')
    axes[1].plot(x, prediction.real/scale, '--', color=viz.SERIES[1], label='No cavity + seven fitted sources')
    axes[1].plot(x[fitted], target.real[fitted]/scale, 'o', color=viz.SERIES[0], label='Fitted receiver positions')
    axes[1].set(xlabel='Surface receiver position (m)', ylabel='Real component / target peak', title='Agreement at sensors need not generalize')
    axes[1].legend(fontsize=8)
    fig.savefig(OUT/'source_ambiguity.png', dpi=160)
    plt.close(fig)
    return {"rows": rows, "cavity_depth_m": 150., "fitted_receivers": 7, "held_out_receivers": 14,
            "scope": "Single-frequency displacement before radar noise. Unrestricted complex source weights, not an ambient-source probability model. Additional receivers, bandwidth and source constraints may distinguish these alternatives."}


def independent_null_benchmark():
    """Independent image realizations preserve the full correlated filter bank."""
    from sarsim.synth import synthesize
    from sarsim.scene import point_targets
    from sarsim.gates import two_arc_correction, ellipse_windows
    geom, bank = Geometry(), SubapBank()
    shape = (512, 256)
    target_rows = np.arange(64, 449, 64)
    target_cols = np.full(len(target_rows), 128)
    arc_rows = np.rint(np.linspace(40, 472, 10)).astype(int)
    rows = np.concatenate((target_rows, arc_rows, arc_rows))
    cols = np.concatenate((target_cols, np.full(10, 176), np.full(10, 80)))
    kz = geom.kz(bank.bands(geom, shape[0])['nu_c'], .48)
    z = default_depths(kz, n=100)
    results = []
    for scene in range(32):
        rng = np.random.default_rng(40000+scene)
        scatterers = None
        family = 'speckle' if scene < 16 else 'speckle_and_points'
        if scene >= 16:
            scatterers = point_targets(rng.uniform(-45, 45, 10), rng.uniform(-40, 40, 10), np.zeros(10), rng.uniform(2, 8, 10), rng=50000+scene)
        slc = synthesize(scatterers, geom, shape, ground_intensity=.1, seed=40000+scene)
        Y = trajectories(slc, geom, bank, rows, cols, upsample=100, verbose=False)
        q, qc = two_arc_correction(Y[:7], Y[7:17], Y[17:])
        gates = ellipse_windows(q)
        score, _ = focus_windows(q, kz, z)
        results.append({"scene": scene, "seed": 40000+scene, "family": family, "targets": 7,
                        "full_gate_pixels": int(gates.full.any(axis=1).sum()),
                        "shape_gate_pixels": int(gates.shape.any(axis=1).sum()),
                        "any_full_gate": bool(gates.full.any()), "mean_ungated_best_score": float(score.max(axis=1).mean()),
                        "median_arc_disagreement_px": float(np.median(qc['disagreement']))})
        if scene % 8 == 7:
            print(f'  independent scene controls: {scene+1}/32', flush=True)
    write_csv('independent_scene_controls.csv', results)
    groups = []
    for family in ['speckle', 'speckle_and_points']:
        subset = [r for r in results if r['family'] == family]
        hits = sum(r['any_full_gate'] for r in subset)
        groups.append({"family": family, "scenes": len(subset), "scenes_with_any_gate_pass": hits,
                       "scene_pass_rate": hits/len(subset), "scene_pass_ci": wilson(hits, len(subset)),
                       "pixel_passes": sum(r['full_gate_pixels'] for r in subset)})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for i, family in enumerate(['speckle', 'speckle_and_points']):
        subset = [r for r in results if r['family'] == family]
        xx = np.arange(len(subset))+i*18
        axes[0].bar(xx, [r['mean_ungated_best_score'] for r in subset], color=viz.SERIES[i], label=family.replace('_', ' '))
        axes[1].scatter(xx, [r['full_gate_pixels'] for r in subset], color=viz.SERIES[i], s=25)
    axes[0].set(xlabel='Independent scene realization', ylabel='Mean best ungated fit score', title='Structured fits on motionless surfaces')
    axes[1].set(xlabel='Independent scene realization', ylabel='Targets passing the full filter / 7', title='Apply the entire selection rule', ylim=(-.3, 7.3))
    axes[0].legend(fontsize=8)
    fig.savefig(OUT/'independent_controls.png', dpi=160)
    plt.close(fig)
    return {"scene_count": 32, "groups": groups, "rows": results,
            "method": "Independent band-limited speckle, with/without ideal static points; 7 target patches per image. Separate 10-sample flanks 48 px from the target strip, per-side geometric median, W25 per-window gating. No permutation of correlated looks. The scene is the independent unit; gate acceptance is not a validated cavity detector."}


def main():
    start = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    viz.style()
    write_json('design.json', DESIGN)
    acq = ReflectorAcquisition(Geometry())
    metrics = {"date": time.strftime('%Y-%m-%d'), "design": DESIGN,
               "acquisition": {"geometry": acq.geometry.summary(), "shape": list(acq.shape),
                               "ground_to_los_projection": acq.projection,
                               "aperture_sample_count": len(acq.time)},
               "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}}
    for name, run in [('motion', lambda: motion_benchmark(acq)), ('depth', lambda: depth_benchmark(acq)),
                      ('noise_sweep', lambda: noise_sweep(acq)),
                      ('frequency_depth', lambda: frequency_depth_benchmark(acq)),
                      ('ambiguity', ambiguity_benchmark), ('independent_controls', independent_null_benchmark)]:
        print(f'Running {name}', flush=True)
        metrics[name] = run()
        write_json('metrics.partial.json', metrics)
    metrics['elapsed_seconds'] = time.monotonic()-start
    metrics['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path(__file__)]+sorted((ROOT/'sarsim').glob('*.py'))}
    write_json('metrics.json', metrics)
    print(f'Completed physical validation in {metrics["elapsed_seconds"]:.1f}s', flush=True)


if __name__ == '__main__':
    main()
