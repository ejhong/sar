"""Method illustration: the Doppler spectrum with reference/offset bands, and sub-aperture images of the pyramid."""
import json, os, numpy as np
from common import *
from sarsim import viz
import matplotlib.pyplot as plt


def main():
    viz.style()
    g = geometry(); bank = default_bank()
    slc = make_slc("pyramid_rot8", pyramid_scene(g, rotation_deg=8.0), g, rotation_deg=8.0)
    x, r = axes(g)
    spec = np.fft.fft(slc, axis=0)
    nu = np.fft.fftfreq(slc.shape[0], d=g.dx)
    order = np.argsort(nu)
    pw = (np.abs(spec) ** 2).mean(axis=1)[order]
    b = bank.bands(g, slc.shape[0])
    fd = os.path.join(RESULTS, "method", "figs"); os.makedirs(fd, exist_ok=True)
    fig = plt.figure(figsize=(13.5, 7.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.9])
    ax = fig.add_subplot(gs[0, :])
    ax.plot(nu[order], 10 * np.log10(pw / pw.max() + 1e-12), color=viz.INK2, lw=1)
    for i, k in enumerate([0, 24, 49]):
        lo = b["nu_lo_ref"][k]; w = b["width"]
        ax.axvspan(lo, lo + w, color=viz.SERIES[i], alpha=0.18, lw=0)
        ax.axvline(lo + b["dnu"], color=viz.SERIES[i], lw=1, ls="--")
        ax.axvline(lo + w + b["dnu"], color=viz.SERIES[i], lw=1, ls="--")
        ax.text(lo + w / 2, -2.5, f"pair k = {k+1}", ha="center", color=viz.SERIES[i], fontsize=9)
    ax.set_xlim(-2.1, 2.1); ax.set_ylim(-40, 1)
    ax.set_xlabel("azimuth spatial frequency ν (cycles/m)   =   Doppler = V·ν   =   slow time   =   look angle")
    ax.set_ylabel("power (dB)")
    ax.set_title("Azimuth spectrum of the SLC with three of the 50 reference bands (shaded) and their offset bands (dashed edges, +88 Hz)", loc="left")
    masks = bank.masks(g, slc.shape[0])
    rs = slice(np.searchsorted(x, -135), np.searchsorted(x, 135)); cs = slice(np.searchsorted(r, -125), np.searchsorted(r, 60))
    for i, k in enumerate([0, 24, 49]):
        R, _ = bank.pair(spec, k, masks)
        axi = fig.add_subplot(gs[1, i])
        viz.imshow_slc(axi, R[rs, cs], x[rs], r[cs], floor=-35)
        axi.set_title(f"sub-aperture image, pair k = {k+1} (look angle {g.nu_to_aspect_deg(b['nu_c'][k]):+.1f}°)", loc="left")
        if i:
            axi.set_ylabel("")
    fig.tight_layout()
    path = viz.finish(fig, f"{fd}/m00_subapertures.png")
    with open(os.path.join(RESULTS, "method", "summary.json"), "w") as f:
        json.dump({"id": "method", "figures": [{"file": path, "caption":
            "<b>Sub-apertures.</b> Top: the image's azimuth spectrum, whose axis is simultaneously Doppler frequency, slow time along the 3.9 s aperture and look angle across 2.8°. Three of the fifty reference/offset band pairs are drawn. Bottom: the pyramid seen through those three bands. Each is a lower-resolution image from a slightly different angle; the speckle and the bright block facets change between them. The pipeline's 'micro-motion' is the sub-pixel shift between the reference and offset versions of each 32 × 32 patch, tracked across the fifty pairs."}]}, f, indent=2)
    print("method figure written")


if __name__ == "__main__":
    main()
