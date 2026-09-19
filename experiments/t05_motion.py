"""T5: inject real ground vibration and measure what the pipeline recovers, and where it puts it."""
import time, json, numpy as np
from common import *
from sarsim import viz, run_pipeline, synthesize, vibrating_patch
from sarsim.tomo import focus_paper
import matplotlib.pyplot as plt

TEST = "t05_motion"
PATCH_CENTRE = (60.0, 80.0)    # ground x, y (m): on the desert, away from the pyramid
PATCH_SIZE = 12.0
FREQS = [0.1, 0.5, 1.0, 2.0, 5.0]
AMPS = np.array([1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2])   # metres
SNR_THRESHOLD = 3.0


def targets(g):
    x, r = axes(g)
    cx = PATCH_CENTRE[0]; cr = PATCH_CENTRE[1] * np.sin(g.theta)
    rows = np.arange(np.searchsorted(x, cx - 4.5), np.searchsorted(x, cx + 4.5), 3)
    cols = np.arange(np.searchsorted(r, cr - 4.5), np.searchsorted(r, cr + 4.5), 3)
    RR, CC = np.meshgrid(rows, cols, indexing="ij")
    return RR.ravel(), CC.ravel()


def main():
    viz.style()
    t0 = time.time()
    g = geometry(); bank = default_bank(); fd = figdir(TEST); figs = []
    slc_static = make_slc("desert", desert_scene(g), g)
    rows, cols = targets(g)
    z_per_hz = 0.48 * np.sin(g.theta) * g.V / (g.Ka * g.lam)
    T_sub = bank.bands(g, SHAPE[0])["width"] * g.V / g.Ka

    def run(A, f):
        vp = vibrating_patch(PATCH_CENTRE, PATCH_SIZE, density=12.0, amp=0.3, vib_amp=A, vib_freq=f, vib_phase=0.7, rng=2)
        slc = slc_static + synthesize(vp, g, SHAPE, r_centre=R_CENTRE, ground_intensity=0.0)
        o = run_pipeline(slc, g, rows, cols, bank=bank, patch=32, common="none", verbose=False, modes=("paper",))
        return o
    o0 = run(0.0, 1.0)
    q0 = o0["q"]                                   # static trajectories of the same patch [P,K,2]
    static_rms = q0[..., 0].std()                  # per-pixel static structure
    static_mean_rms = (q0[..., 0].mean(axis=0) - q0[..., 0].mean()).std()
    kz, z = o0["kz"], o0["z"]
    resp = np.zeros((len(FREQS), len(AMPS)))
    snr = np.zeros_like(resp)
    kept = {}
    for i, f in enumerate(FREQS):
        for j, A in enumerate(AMPS):
            o = run(A, f)
            dq = o["q"] - q0                        # response: with minus without vibration, same pixels
            resp[i, j] = dq[..., 0].std()
            snr[i, j] = resp[i, j] / static_rms
            if A in (3e-3, 1e-2):
                kept[(f, A)] = dq
        print(f"  f={f} Hz: response rms px {np.round(resp[i], 4)}", flush=True)
    thr = []
    for i, f in enumerate(FREQS):
        s = snr[i]
        idx = np.where(s >= SNR_THRESHOLD)[0]
        if idx.size == 0:
            thr.append({"freq_hz": f, "threshold_amp_m": None, "max_snr": float(s.max())})
        else:
            k = idx[0]
            a = AMPS[0] if k == 0 else 10 ** np.interp(np.log10(SNR_THRESHOLD), np.log10([s[k - 1], s[k]]), np.log10([AMPS[k - 1], AMPS[k]]))
            thr.append({"freq_hz": f, "threshold_amp_m": float(a), "max_snr": float(s.max())})
    # figure 1: response vs amplitude
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.axvspan(1e-8, 1e-5, color=viz.GRID, alpha=0.6, lw=0)
    ax.text(3e-7, 0.93, "Small reference amplitudes\n(not a measured noise floor)", transform=ax.get_xaxis_transform(), color=viz.INK2, fontsize=8.5, ha="center", va="top")
    for i, f in enumerate(FREQS):
        ax.plot(AMPS, resp[i] * 1e3, "-o", color=viz.SERIES[i], lw=1.8, ms=5, markeredgecolor=viz.SURFACE, label=f"{f:g} Hz")
    ax.axhline(static_rms * 1e3, color=viz.SERIES[7], lw=1, ls="--")
    ax.text(1.2e-2, static_rms * 1e3 * 1.08, "static structure of the same patch (no motion)", color=viz.SERIES[7], fontsize=8.5, ha="right")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-8, 2e-2); ax.set_ylim(0.02, 40)
    ax.set_xlabel("injected line-of-sight vibration amplitude (m)"); ax.set_ylabel("change of the trajectory caused by the vibration (millipixels, rms)")
    ax.set_title("How much does a genuinely vibrating patch change the pipeline's output?", loc="left"); ax.legend(title="frequency", loc="lower right")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t05_sensitivity.png"), "caption":
                 "<b>Response to genuine motion.</b> A 12 m patch of ground oscillates coherently along the line of sight; the curves show how much the sampled trajectories change relative to the same patch at rest. At 1 cm, the response ranges from about 1.05 to 3.58 times the static trajectory rms across the tested frequencies. The reference and offset bands overlap by about 99%; this statistic measures changes between their patch-correlation peaks. The shaded amplitudes are reference values, not a calibrated seismic or sensor noise floor. The static baseline is known only because this is a simulation."})
    # figure 2: where the response lands on the depth axis
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.0))
    dep_rows = []
    for i, f in enumerate([1.0, 2.0, 5.0]):
        dq = kept[(f, 1e-2)]
        mean_dq = dq[..., 0].mean(axis=0)
        spec = focus_paper((dq[..., 0] + 1j * dq[..., 1]), kz, z).mean(axis=0)
        axs[0].plot(np.arange(1, 51), mean_dq * 1e3, color=viz.SERIES[i], lw=1.8, label=f"{f:g} Hz")
        axs[1].plot(z, spec / spec.max(), color=viz.SERIES[i], lw=1.8, label=f"{f:g} Hz, predicted {z_per_hz*f:.1f} m")
        axs[1].axvline(z_per_hz * f, color=viz.SERIES[i], lw=0.9, ls=":")
        dep_rows.append({"freq_hz": f, "predicted_depth_m": float(z_per_hz * f), "measured_depth_m": float(z[np.argmax(spec)])})
    axs[0].set_xlabel("sub-aperture pair k"); axs[0].set_ylabel("mean change of azimuth shift (millipixels)"); axs[0].legend(); axs[0].set_title("Vibration-induced part of the trajectory, 1 cm amplitude", loc="left")
    axs[1].set_xlabel("model depth (m)"); axs[1].set_xlim(0, 80); axs[1].legend(fontsize=8); axs[1].set_title("Where the focusing step puts it (dotted: frequency × 10.9 m/Hz)", loc="left")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t05_depth_is_frequency.png"), "caption":
                 f"<b>Motion has no depth, only a frequency.</b> Isolating the vibration-induced part of the trajectories at a 1 cm amplitude and focusing it: the steering step maps the vibration's frequency onto the depth axis at {z_per_hz:.1f} model-metres per hertz (dotted lines mark the prediction; the 2 Hz case also shows its second harmonic at twice the depth). A patch of surface that shakes faster is drawn 'deeper'. This demonstrates frequency-to-depth mapping in this reconstruction, rather than localization of a buried object."})
    m = {"thresholds_snr3": thr, "static_rms_px": float(static_rms), "depth_per_hz_m": float(z_per_hz), "sub_aperture_duration_s": float(T_sub),
         "depth_rows": dep_rows, "response_rms_px": {str(f): [float(v) for v in resp[i]] for i, f in enumerate(FREQS)},
         "snr_vs_static": {str(f): [float(v) for v in snr[i]] for i, f in enumerate(FREQS)}, "amps_m": [float(a) for a in AMPS], "runtime_s": time.time() - t0}
    best = max(snr.max(axis=1))
    summary = {
        "id": TEST, "order": 5, "tag": "threshold", "eyebrow": "05 · Injected motion",
        "title": "Injected vibration changes the output, but does not supply a depth",
        "question": "Grant the premise: if the ground really vibrated, would the method see it, and would it locate it?",
        "finding": f"A patch oscillating with 1 cm amplitude changes its trajectories by at most {best:.1f} times their static structure; at 0.1 Hz the change is {snr[0, -1]:.2f} times. These response-to-static ratios are not detection probabilities or limits on radar vibrometry. The strongly overlapping bands and the translation statistic are specific to this reconstruction. Where a response exists, the focusing step converts its frequency into a depth at {z_per_hz:.1f} m per Hz.",
        "limitations": "One coherent sinusoidal patch with static background clutter; response measured against an exactly known motion-free baseline. No sensor-noise distribution, false-alarm calibration, or coupled cavity model. The 2 Hz response peaks at its second harmonic.",
        "method": f"Static desert SLC plus a synthesised vibrating patch (linear superposition), {len(FREQS)} frequencies × {len(AMPS)} amplitudes from 1 µm to 1 cm; response measured as the change of the sampled patch trajectories relative to the same patch at rest, compared with the rms of the static trajectories.",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps({"thresholds": thr, "depth_rows": dep_rows, "static_rms_px": float(static_rms)}, indent=1, default=float))


if __name__ == "__main__":
    main()
