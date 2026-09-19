"""Generate the opening figure directly from the simulated SLC and analysis arrays."""
from pathlib import Path
import numpy as np
from common import geometry, pyramid_scene, make_slc, analyze_scene, axes, ROOT
from sarsim import viz
import matplotlib.pyplot as plt


def main():
    g = geometry()
    scene = pyramid_scene(g, rotation_deg=8)
    slc = make_slc("pyramid_rot8", scene, g, rotation_deg=8.0)
    analysis = analyze_scene("pyramid_rot8", slc, g)
    x, r = axes(g)
    z = analysis["z"]
    spectrum = analysis["az_tomo_paper"]
    spectrum = spectrum / np.maximum(spectrum.sum(axis=1, keepdims=True), 1e-30)
    viz.style()
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.25), layout="constrained", gridspec_kw={"width_ratios": [1, 1.4]})
    viz.imshow_slc(axs[0], slc, x, r)
    axs[0].axhline(-20, color="#ad573b", lw=1, ls="--")
    axs[0].set(title="01  /  A stationary surface", xlabel="Azimuth (m)", ylabel="Slant range (m)")
    rows = analysis["az_rows"]
    im = axs[1].imshow(spectrum.T, aspect="auto", origin="upper", interpolation="nearest",
                       extent=[x[rows[0]], x[rows[-1]], z[-1], z[0]], cmap=viz.CMAP_SEQ,
                       vmin=0, vmax=np.percentile(spectrum, 99.5))
    axs[1].set(title="02  /  Its reconstructed depth spectrum", xlabel="Azimuth along the dashed transect (m)", ylabel="Model depth (m)")
    fig.colorbar(im, ax=axs[1], shrink=.75, label="Relative spectral power")
    destination = Path(ROOT) / "results" / "overview"
    destination.mkdir(exist_ok=True)
    viz.finish(fig, destination / "overview.png")


if __name__ == "__main__":
    main()
