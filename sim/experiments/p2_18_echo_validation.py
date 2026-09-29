"""P2-18 · The radar simulator, checked against pulses: echoes pulse by pulse, focused by back-projection.

    uv run --with scikit-image python experiments/p2_18_echo_validation.py

Every motion test in this lab so far wrote focused images directly in the azimuth-frequency domain (the fast
synthesizer, and injections on real images), mapping each Doppler frequency to a moment of the pass and applying
motion there as a phase. An independent review asked for that to be checked against a calculation of what a radar
records. sarsim.echo computes, for each pulse, the range-compressed echo of every scatterer from its exact distance at
that moment, the scatterer displaced by its own motion, and focuses by time-domain back-projection on the product's
own grid (the 2022 X13 geometry: 24.56 s, Ka -5,401 Hz/s, V_eff 7,354 m/s, 300 MHz, 1,000 pulses a second).

A. Resolution: a still point's response on a grid eight times finer than the product's, against theory.
B. A Doppler slice is a stretch of the pass: the frequency slice of the whole image against the image formed from only
   the pulses of that stretch, for speckle and a bright point.
C. The fast synthesizer against pulses: the same scatterers, still, imaged both ways.
D. A constant line-of-sight velocity shifts the image in azimuth by R v / V_s: pulses against theory.
E. A positive control, as the 2026 measurements did it: a corner reflector 50 dB above clutter at full resolution (a
   trihedral's usual margin; about 23 dB in a 50 ms look) vibrating at 1 and
   5 mm/s, at 1, 2 and 4 Hz, read from one image by tracking its magnitude through short looks.
F. The gated reconstruction's six-second pairs (its own masks and registration) on a vibrating bright point, from
   pulses and from the injection P2-15 used.
"""
import dataclasses
import json

import numpy as np

import p2_13_gates as m
import p2_15_pair_response as r15
from katabasis.runs import Run
from sarsim.acquisition import DwellGeometry
from sarsim.echo import EchoSetup, PointScene, simulate
from sarsim.looks import inject_region_motion, look_masks, looks, velocity_series

RID = 'p2_18_echo_validation'
PRF = 1000.0


def corr(a, b):
    a, b = a.ravel(), b.ravel()
    return float(abs(np.vdot(a, b)) / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real))


def width(p, spacing):
    p = p / p.max()
    i = int(p.argmax())
    lvl = 1 / np.sqrt(2)
    lo, hi = i, i
    while p[lo] > lvl:
        lo -= 1
    while p[hi] > lvl:
        hi += 1
    fl = lo + (lvl - p[lo]) / (p[lo + 1] - p[lo])
    fr = hi - 1 + (p[hi - 1] - lvl) / (p[hi - 1] - p[hi])
    return float((fr - fl) * spacing)


def patch(st, rng, n, half_u_m, half_r_m, bright=None):
    """Random scatterers over a patch (azimuth half-length in ground metres, range half-width), and optionally a bright
    one at the centre; u is converted to the equivalent slant-plane coordinate."""
    u = rng.uniform(-half_u_m, half_u_m, n) * st.V_eff / st.V_g
    r = rng.uniform(-half_r_m, half_r_m, n)
    a = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    if bright is not None:
        u, r, a = np.r_[u, 0.0], np.r_[r, 0.0], np.r_[a, bright]
    return PointScene.still(u, r, a)


def main():
    g = DwellGeometry.from_record(m.ACQ)
    st = EchoSetup.from_geometry(g, prf=PRF)
    params = {'acquisition': m.ACQ, 'prf_simulated': PRF, 'V_eff_m_s': st.V_eff, 'bandwidth_hz': st.bandwidth_hz,
              'dwell_s': st.dwell_s, 'gated_module_sha256': m.sha256(m.MODULE)}
    with Run(RID, 'The radar simulator, checked against pulses', params) as run:
        rng = np.random.default_rng(18)
        out = {}

        # A. resolution
        fine = dataclasses.replace(st, row_dt=st.row_dt / 8, dr=st.dr / 8)
        img = np.abs(simulate(fine, PointScene.still([0.0], [0.0], [1.0]), rows=256, cols=128))
        iy, ix = np.unravel_index(img.argmax(), img.shape)
        prof = img[:, ix] / img[:, ix].max()
        j = iy + 1
        while prof[j + 1] < prof[j]:
            j += 1
        k = j
        while prof[k + 1] > prof[k]:
            k += 1
        out['resolution'] = {'azimuth_m': width(img[:, ix], g.dx / 8), 'azimuth_theory_m': 0.886 * g.V / (g.Ka * st.dwell_s),
                             'range_m': width(img[iy, :], g.dr / 8), 'range_theory_m': 0.886 * st.rho,
                             'azimuth_first_sidelobe_db': float(20 * np.log10(prof[k]))}
        print('A', out['resolution'], flush=True)

        # B. a Doppler slice is a stretch of the pass
        sc = patch(st, rng, 3000, 8.0, 3.0, bright=30.0)
        rows, cols = 1024, 32
        full = simulate(st, sc, rows, cols)
        slices = []
        for centre in (-9.0, -6.0, -3.0, 0.0, 3.0, 6.0, 9.0):
            W = 4.0
            part = simulate(st, sc, rows, cols, pulses=(centre - W / 2, centre + W / 2))
            sl = looks(full, look_masks(g, rows, [centre], W))[0]
            mirror = looks(full, look_masks(g, rows, [-centre], W))[0] if centre else None
            slices.append({'centre_s': centre, 'corr_same': corr(part, sl),
                           'corr_mirror': corr(part, mirror) if mirror is not None else None})
        out['slice_is_stretch'] = slices
        print('B', [(s['centre_s'], round(s['corr_same'], 4)) for s in slices], flush=True)

        # C. the fast synthesizer against pulses, the same still scatterers
        from sarsim import synthesize
        from sarsim.scene import Scatterers
        rows_c, cols_c = 2048, 32
        scc = patch(st, rng, 2000, 30.0, 3.0, bright=30.0)
        pulse_img = simulate(st, scc, rows_c, cols_c)
        x_ground = scc.u * st.V_g / st.V_eff
        z = np.zeros_like(x_ground)
        # in the SLC convention each scatterer carries the phase of its own range, -4 pi r / lambda (the synthesizer's
        # scatterer phase is its reflectivity alone), so the same scene is given that phase to compare pixel by pixel
        scat = Scatterers(x=x_ground, y=(scc.r + 0.0) / np.sin(g.theta), z=z, amp=np.abs(scc.a),
                          phase=np.angle(scc.a) - 4 * np.pi * scc.r / g.lam,
                          iso=np.ones_like(z), flash=z, nu0=z, sig_nu=np.ones_like(z), vib_amp=z, vib_freq=z, vib_phase=z,
                          label=np.zeros(len(z), int))
        fast_img = synthesize(scat, g, (rows_c, cols_c))
        out['fast_vs_pulses'] = {'corr_complex': corr(pulse_img, fast_img), 'corr_magnitude':
                                 float(np.corrcoef(np.abs(pulse_img).ravel(), np.abs(fast_img).ravel())[0, 1])}
        print('C', out['fast_vs_pulses'], flush=True)

        # D. constant velocity -> azimuth shift R v / V_s
        v = 5e-3
        rows_d, cols_d = 1024, 16
        still = simulate(st, PointScene.still([0.0], [0.0], [1.0]), rows_d, cols_d)
        # a steady recession v: dr(t) = v t, written as a very slow sinusoid of the same slope at t = 0
        f_slow = 1e-4
        moving = PointScene.still([0.0], [0.0], [1.0]).moving(np.array([True]), ar=v / (2 * np.pi * f_slow), f=f_slow)
        mv = simulate(st, moving, rows_d, cols_d)
        up = 64
        def peak(im):
            prof_ = np.abs(im[:, cols_d // 2])
            F = np.fft.fft(prof_)
            Fz = np.zeros(len(F) * up, complex)
            h = len(F) // 2
            Fz[:h], Fz[-h:] = F[:h], F[-h:]
            fineprof = np.abs(np.fft.ifft(Fz))
            return fineprof.argmax() / up
        shift_px = peak(mv) - peak(still)
        out['steady_velocity'] = {'v_m_s': v, 'shift_m': float(shift_px * g.dx),
                                  'theory_m': float(v * 2 * g.V / (g.lam * g.Ka_signed)),
                                  'velocity_read_back_m_s': float(g.shift_to_velocity(shift_px * g.dx))}
        print('D', out['steady_velocity'], flush=True)

        # E. positive control: a vibrating corner reflector, read from one image by magnitudes through short looks
        rows_e, cols_e = 4096, 16
        base = patch(st, rng, 1500, 90.0, 2.0)
        amp_cr = np.sqrt(10 ** 5.0 * 1500 / (180.0 / g.resolution * 4.0 / 0.44))    # ~50 dB above the clutter's mean
        refl = PointScene.still(np.r_[base.u, 0.0], np.r_[base.r, 0.0], np.r_[base.a, amp_cr])
        sel = np.zeros(len(refl.u), bool)
        sel[-1] = True
        W_look = 0.05
        centres = np.arange(-2.0, 2.0 + 1e-9, 0.025)
        masks_e = look_masks(g, rows_e, centres, W_look)
        rr = slice(rows_e // 2 - 400, rows_e // 2 + 400)
        cc = slice(cols_e // 2 - 2, cols_e // 2 + 3)
        control = []
        for vamp in (1e-3, 5e-3):
            for f in (1.0, 2.0, 4.0):
                img_m = simulate(st, refl.moving(sel, ar=vamp / (2 * np.pi * f), f=f), rows_e, cols_e)
                vel = velocity_series(looks(img_m, masks_e), g, rr, cc, envelope=True)
                A = np.stack([np.cos(2 * np.pi * f * centres), np.sin(2 * np.pi * f * centres), np.ones_like(centres)], 1)
                coef, *_ = np.linalg.lstsq(A, vel, rcond=None)
                got = float(np.hypot(coef[0], coef[1]))
                truth = vamp * abs(np.sinc(f * W_look))          # the look's own average of the velocity
                resid = vel - A @ coef
                sigma = float(np.std(resid) * np.sqrt(2.0 / len(centres)))    # the fitted amplitude's standard error
                control.append({'v_m_s': vamp, 'f_hz': f, 'recovered_m_s': got, 'truth_look_average_m_s': float(truth),
                                'gain': got / truth, 'gain_se': sigma / truth, 'residual_rms_m_s': float(np.std(resid)),
                                'resolved': got > 3 * sigma})
                print('E', control[-1], flush=True)
        out['positive_control'] = {'look_s': W_look, 'looks': len(centres), 'reflector_scr_db': 50.0, 'rows': control}

        # F. the six-second pairs on a vibrating bright point: pulses against the injection P2-15 used
        profile = json.loads(m.PROFILE.read_text())
        gated = m.load_gated()
        masks_f, records = gated.frequency_masks(r15.SHAPE[0], 1.0 / g.prf, profile)
        t_pair = np.array([q['pair_center_hz'] for q in records]) / g.Ka_signed
        reg = profile['processing']['registration']
        rows_f, cols_f = r15.SHAPE
        pt_still = simulate(st, PointScene.still([0.0], [0.0], [1.0]), rows_f, cols_f)
        pt_still /= np.abs(pt_still).max()
        Y0p, _ = r15.register_series(pt_still, masks_f, reg, gated)
        inj_still = r15.point(g, r15.SHAPE)
        Y0i, _ = r15.register_series(inj_still, masks_f, reg, gated)
        pairs = []
        for vamp in (2e-3, 2e-2):
            for f in (0.26, 1.0691, 3.66):
                sc_m = PointScene.still([0.0], [0.0], [1.0]).moving(np.array([True]), ar=vamp / (2 * np.pi * f), f=f)
                img_p = simulate(st, sc_m, rows_f, cols_f)
                img_p /= np.abs(img_p).max()
                Y1p, _ = r15.register_series(img_p, masks_f, reg, gated)
                d = lambda t, f=f, vamp=vamp: (vamp / (2 * np.pi * f)) * np.sin(2 * np.pi * f * t)
                img_i = inject_region_motion(inj_still, g, np.ones(rows_f), d)
                Y1i, _ = r15.register_series(img_i, masks_f, reg, gated)
                ap, _ = r15.fit(Y1p - Y0p, t_pair, f)
                ai, _ = r15.fit(Y1i - Y0i, t_pair, f)
                pairs.append({'v_m_s': vamp, 'f_hz': f, 'pulses_px': float(ap[1]), 'injection_px': float(ai[1])})
                print('F', pairs[-1], flush=True)
        out['pairs_pulses_vs_injection'] = pairs

        res = out['resolution']
        sl = [s['corr_same'] for s in out['slice_is_stretch']]
        mir = [s['corr_mirror'] for s in out['slice_is_stretch'] if s['corr_mirror'] is not None]
        clean = [c for c in control if c['resolved']]
        noisy = [c for c in control if not c['resolved']]
        noisy_mm = ', '.join(sorted({f"{c['v_m_s'] * 1e3:.0f}" for c in noisy}))
        finding = (
            f"Pulse by pulse, on the 2022 geometry, the simulator focuses a point to {res['azimuth_m'] * 100:.2f} cm in azimuth "
            f"(theory {res['azimuth_theory_m'] * 100:.2f}) with a first sidelobe of {res['azimuth_first_sidelobe_db']:.1f} dB, and to "
            f"{res['range_m']:.3f} m in range against {res['range_theory_m']:.3f} m for a flat 300 MHz band (the dwell's 16 degrees "
            f"widen the range spectrum). A frequency slice of the whole image matches the image formed from only that stretch of "
            f"pulses to {min(sl):.3f}-{max(sl):.3f} in complex correlation, and the mirrored slice to at most {max(mir):.3f}. The fast "
            f"synthesizer's image of the same still scatterers matches the pulses to {out['fast_vs_pulses']['corr_complex']:.3f}. A "
            f"steady {out['steady_velocity']['v_m_s'] * 1e3:.0f} mm/s recession shifts the image by "
            f"{out['steady_velocity']['shift_m'] * 100:.2f} cm against {out['steady_velocity']['theory_m'] * 100:.2f} cm from R v / V_s. "
            f"As a positive control, a corner reflector 50 dB above clutter is read from one image by magnitudes through "
            f"{W_look * 1000:.0f} ms looks: vibrating at "
            + ', '.join(f"{c['v_m_s'] * 1e3:.0f} mm/s and {c['f_hz']:.0f} Hz at {c['gain']:.2f} ± {c['gain_se']:.2f}" for c in clean)
            + " of its look-averaged velocity"
            + (f"; at {noisy_mm} mm/s it lies "
               f"within the per-look noise ({np.mean([c['residual_rms_m_s'] for c in noisy]) * 1e3:.0f} mm/s), as it would need a "
               f"brighter reflector or a longer series" if noisy else "")
            + ". Through the six-second pairs, pulses and the injection give the bright point's pair shifts of "
            + '; '.join(f"{q['pulses_px']:.4f} and {q['injection_px']:.4f} px ({q['v_m_s'] * 1e3:.0f} mm/s, {q['f_hz']:.2f} Hz)" for q in pairs)
            + ".")
        run.save({**out, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
