"""P2-06 · The published depth is the frequency of the ground's motion.

    uv run python experiments/p2_06_depth_is_frequency.py

The published focusing (Biondi & Malanga 2022, eqs. 21-24; sarsim.tomo) projects each patch's trajectory
across the sub-aperture pairs on exp(j Kz_k z), with Kz_k = 4 pi B_k / (lambda_s R sin theta) and B_k the
platform's offset at pair k's slow time t_k. On a dwell B grows linearly with t_k (the platform moves at a
steady speed), so Kz_k = Kz' t_k and the focusing is a Fourier transform over slow time: motion at frequency
f can only appear at

    z = 2 pi f / Kz' = f lambda_s R sin theta / (2 V_s),

whatever is under the ground. This experiment gives the pipeline the best trajectory it could ever have,
the patch's exact line-of-sight velocity at every pair, with no noise and none of P2-03's blindness, and
reads where it puts the motion:

1. on the real Giza geometry (P2-01) with the 2022 paper's bank (K = 50, half-band sub-apertures, 88 Hz
   offset) and its declared lambda_s = 0.48 m: depth against frequency, and where faster motion folds;
2. over two chambers, 15 m and 30 m down, and over open ground, under the same microseism: the imprint's
   amplitude differs, its place on the depth axis does not;
3. the same with the 2026 replication protocol's windowed fit (branch A).
"""
import numpy as np

from katabasis.runs import Run
from katabasis.seismic.analytic import moment_surface_displacement, void_moment
from sarsim.acquisition import DwellGeometry
from sarsim.subap import SubapBank
from sarsim.tomo import default_depths, focus_paper, focus_windows, kz_for_bank, nyquist_depth

LAM_S = 0.48                    # m, the 2022 paper's declared "sound" wavelength (legacy real_common.LAM_S)
N_ROWS = 113692                 # the Giza product's azimuth samples, for the bank's bin snapping
FREQS = [0.1, 0.2, 0.3, 0.45, 0.6, 1.0, 1.5, 1.9, 2.3, 3.0]
SITE_TONES = {'giza_ground': (0.6, 'Giza ground resonance (ELGabry et al. 2026)'),
              'khufu': (2.3, "Khufu's own resonance (ELGabry et al. 2026)"),
              'microseism': (0.2, 'microseism peak')}


def trajectory(t, f, phase=0.3):
    """The patch's exact line-of-sight velocity at each pair (azimuth component), nothing in range."""
    q = np.zeros((1, len(t), 2))
    q[0, :, 0] = np.cos(2 * np.pi * f * t + phase)
    return q


def main():
    g = DwellGeometry.from_record('giza-20250827')
    bank = SubapBank()
    t = g.nu_to_time(bank.bands(g, N_ROWS)['nu_c'])
    kz = kz_for_bank(bank, g, LAM_S, N_ROWS)
    slope = float(np.polyfit(t, kz, 1)[0])
    m_per_hz = 2 * np.pi / abs(slope)
    closed = LAM_S * g.R0 * np.sin(g.theta) / (2 * g.V_platform)
    fs = (bank.K - 1) / (t.max() - t.min())
    z = np.linspace(0.02, 0.97 * nyquist_depth(kz), 800)
    params = {'geometry': 'giza-20250827', 'lam_s_m': LAM_S, 'bank': {'K': bank.K, 'sub_frac': bank.sub_frac,
              'delta_hz': bank.delta_hz}, 'frequencies_hz': FREQS}
    with Run('p2_06_depth_is_frequency', 'The published depth is the frequency of the ground\'s motion', params) as run:
        rows = []
        for f in FREQS:
            tomo = focus_paper(trajectory(t, f)[..., 0] + 1j * trajectory(t, f)[..., 1], kz, z)[0]
            alias = abs(f - fs * round(f / fs))
            rows.append({'f_hz': f, 'depth_m': float(z[np.argmax(tomo)]), 'predicted_m': float(m_per_hz * alias),
                         'folded': bool(alias != f)})
        # two chambers and open ground under the same 0.2 Hz microseism: only the amplitude differs
        lam, mu = 1.006e10, 8.04e9
        nu = lam / (2 * (lam + mu))
        M = void_moment(np.diag([1.0, 1.0, 0.0]), 216.0, nu)
        profiles = {}
        for label, depth in (('chamber 15 m down', 15.0), ('chamber 30 m down', 30.0)):
            amp = abs(moment_surface_displacement(M, depth, np.zeros(1), np.zeros(1), lam, mu)[2, 0])
            q = amp * trajectory(t, 0.2)
            tomo = focus_paper(q[..., 0] + 1j * q[..., 1], kz, z)[0]
            score, _ = focus_windows(q, kz, z, W=25)
            profiles[label] = {'epicentre_amplitude_per_pa_m': float(amp), 'depth_m': float(z[np.argmax(tomo)]),
                               'windows_depth_m': float(z[np.argmax(score[0])]), 'tomogram': tomo / tomo.max()}
        q = trajectory(t, 0.2)
        tomo = focus_paper(q[..., 0] + 1j * q[..., 1], kz, z)[0]
        profiles['open ground, the wave itself'] = {'depth_m': float(z[np.argmax(tomo)]), 'tomogram': tomo / tomo.max()}
        same = max(np.abs(p['tomogram'] - profiles['chamber 15 m down']['tomogram']).max() for p in profiles.values())
        tones = {k: {'label': lab, 'f_hz': f, 'depth_m': float(m_per_hz * abs(f - fs * round(f / fs)))}
                 for k, (f, lab) in SITE_TONES.items()}
        finding = (
            f"The published depth is the frequency of the ground's motion, relabelled. On the real Giza geometry, with "
            f"the 2022 paper's bank and its lambda_s of {LAM_S} m, the steering grows by {abs(slope):.3f} rad/m every "
            f"second of slow time, so the focusing is a Fourier transform over slow time and puts motion at "
            f"{m_per_hz:.2f} m of depth per hertz (lambda_s R sin theta / 2 V_s = {closed:.2f} m). Given the exact "
            f"motion at every pair, it draws a 0.2 Hz microseism at {rows[1]['depth_m']:.2f} m over a chamber 15 m "
            f"down, over one 30 m down and over open ground alike (profiles identical to {same:.0e}); only the "
            f"brightness changes. Its pairs sample slow time at {fs:.1f} Hz, so its whole depth axis, 0 to "
            f"{nyquist_depth(kz):.1f} m, is the band 0 to {fs / 2:.1f} Hz and anything faster folds back: Giza's own "
            f"0.6 Hz resonance would be drawn at {tones['giza_ground']['depth_m']:.1f} m and "
            f"Khufu's 2.3 Hz at {tones['khufu']['depth_m']:.1f} m.")
        run.save({'kz_slope_rad_m_s': slope, 'metres_per_hz': m_per_hz, 'closed_form_metres_per_hz': closed,
                  'pair_sampling_hz': fs, 'nyquist_depth_m': float(nyquist_depth(kz)), 'pair_times_s': t,
                  'depth_axis_m': z, 'rows': rows, 'profiles': profiles, 'profile_max_difference': same,
                  'site_tones_m': tones, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
