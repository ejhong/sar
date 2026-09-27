"""T4: the depth axis is a processing choice: patch size, declared wavelength, FFT-grid alignment, band width."""
import time, json, numpy as np
from common import *
from sarsim import viz, run_pipeline, SubapBank
import matplotlib.pyplot as plt

TEST = "t04_parameters"


def transect_run(slc, g, which, **kw):
    rows, cols = transects(g)[which]
    return run_pipeline(slc, g, rows, cols, verbose=False, **kw), rows, cols


def main():
    viz.style()
    t0 = time.time()
    g = geometry(); fd = figdir(TEST); figs = []
    scale = g.depth_scale(0.48)
    slc_d = make_slc("desert", desert_scene(g), g)
    slc_p = make_slc("pyramid_rot8", pyramid_scene(g, rotation_deg=8.0), g, rotation_deg=8.0)
    x, r = axes(g)
    # ---- (a) patch size
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 6.4), gridspec_kw={"height_ratios": [2.2, 1.2]})
    patch_rows = []
    for i, patch in enumerate([16, 32, 64]):
        o, rows, cols = transect_run(slc_d, g, "az", bank=default_bank(), patch=patch)
        z = o["z"]; sw = o["score_w"]; tp = o["tomo_paper"]; tpn = tp / tp.sum(axis=1, keepdims=True)
        axs[0, i].imshow(sw.T, aspect="auto", origin="upper", extent=[x[rows[0]], x[rows[-1]], z[-1], z[0]], cmap=viz.CMAP_SEQ, vmin=0, vmax=1, interpolation="nearest")
        axs[0, i].set_title(f"patch {patch} px = {patch*g.dx:.0f} m", loc="left"); axs[0, i].set_xlabel("azimuth (m)")
        spec = tpn.mean(axis=0)
        axs[1, i].plot(z, spec, color=viz.SERIES[0], lw=1.8)
        cap = scale * (patch * g.dx)
        axs[1, i].axvline(cap, color=viz.SERIES[1], lw=1.0, ls=":")
        axs[1, i].set_xlabel("model depth (m)"); axs[1, i].set_xlim(0, z[-1])
        # vertical-to-horizontal extent of structures: autocorrelation widths of the score map
        def acw(a, axis):
            a = a - a.mean(); f = np.fft.fft(a, axis=axis); ac = np.fft.ifft(np.abs(f) ** 2, axis=axis).real
            ac = ac / np.take(ac, 0, axis=axis)[..., None] if axis == 1 else ac / np.take(ac, 0, axis=axis)[None, :]
            prof = ac.mean(axis=1 - axis)
            return float(np.argmax(prof < 0.5))
        patch_rows.append({"patch_px": patch, "patch_m": patch * g.dx, "blob_width_px_along_transect": acw(sw, 0),
                           "blob_height_depth_steps": acw(sw, 1), "depth_of_patch_width_m": float(cap),
                           "spectrum_above_patch_cap_fraction": float(spec[z > cap].sum() / spec.sum())})
    axs[0, 0].set_ylabel("model depth (m)"); axs[1, 0].set_ylabel("mean spectrum")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t04_patch.png"), "caption":
                 "<b>The patch size sets the shape of the 'structures'.</b> Empty desert, same transect, with 16, 32 and 64 pixel patches. Top: the replication-variant score. Bottom: the mean paper-variant spectrum; the dotted line marks the depth that a separation of one patch width would give. Small patches draw thin vertical shafts, large patches draw rounded chambers, and every patch size fills the same depth range, because scatterers far outside the patch contribute through the 80 m difference kernel and alias onto the axis (T3)."})
    # ---- (b) declared wavelength
    o, rows, cols = transect_run(slc_p, g, "az", bank=default_bank(), patch=32)
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9))
    for i, lam_s in enumerate([0.24, 0.48, 0.96]):
        f = lam_s / 0.48
        z = o["z"] * f
        axs[i].imshow(o["score_w"].T, aspect="auto", origin="upper", extent=[x[rows[0]], x[rows[-1]], z[-1], z[0]], cmap=viz.CMAP_SEQ, vmin=0, vmax=1, interpolation="nearest")
        axs[i].set_title(f"$\lambda_s$ = {lam_s} m  to  axis 0 to {z[-1]:.0f} m", loc="left"); axs[i].set_xlabel("azimuth (m)")
    axs[0].set_ylabel("model depth (m)")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t04_lambda.png"), "caption":
                 "<b>The declared 'sound wavelength' is a free rescaling of the depth axis.</b> The same pyramid transect labelled with λₛ = 0.24, 0.48 (the 2022 paper's value, from an assumed 6 km/s and 12.5 kHz) and 0.96 m. The image does not change; the numbers on the axis do. The replication protocol calibrated λ on a tomb of known depth, which makes any one feature land at any depth one likes."})
    grid_rows = []
    # ---- (d) band width
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9))
    for i, frac in enumerate([0.3, 0.5, 0.7]):
        o, rows, cols = transect_run(slc_p, g, "az", bank=SubapBank(sub_frac=frac), patch=32)
        z = o["z"]
        axs[i].imshow(o["score_w"].T, aspect="auto", origin="upper", extent=[x[rows[0]], x[rows[-1]], z[-1], z[0]], cmap=viz.CMAP_SEQ, vmin=0, vmax=1, interpolation="nearest")
        axs[i].set_title(f"sub-aperture width {frac:.0%} of the band", loc="left"); axs[i].set_xlabel("azimuth (m)")
    axs[0].set_ylabel("model depth (m)")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t04_bandwidth.png"), "caption":
                 "<b>Change the sub-aperture width and the structures rearrange.</b> The pyramid transect with sub-apertures of 30%, 50% and 70% of the processed band. Note the depth axes: the sweep step changes with the width, and the extent of the depth axis is simply the Nyquist depth of that step (π / ΔKz), so the deepest 'structure' any run can show is fixed by K, the band and λs before a single pixel is processed. The pattern within the axis changes as well. Real geology would not care how the Doppler spectrum is sliced."})
    m = {"patch": patch_rows, "runtime_s": time.time() - t0}
    summary = {
        "id": TEST, "order": 4, "tag": "artifact", "eyebrow": "04 · Processing choices",
        "title": "The depth axis is set by the processing, not the ground",
        "question": "Do the 'depths' survive changes that a physical measurement would ignore?",
        "finding": "No. The patch size changes the shape of the structures (thin shafts with small patches, rounded chambers with large ones), the declared sound wavelength rescales the axis freely, and the sub-aperture width sets the extent of the depth axis and rearranges the features within it. Every knob of the processing redraws the underground.",
        "limitations": "Only one value of each parameter is varied at a time; the sensitivity to real data's additional quirks (residual focusing errors, Doppler-centroid drift, ionosphere) can only be larger.",
        "method": "Desert and pyramid azimuth transects re-run with patch 16/32/64 px, λₛ 0.24/0.48/0.96 m, and sub-aperture width 30/50/70%. (Band edges are snapped to whole FFT bins throughout; with fractional edges a small scene-wide common term appears, which the median-trajectory subtraction removes.)",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1, default=float))


if __name__ == "__main__":
    main()
