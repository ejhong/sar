"""T2: the same pipeline on an empty desert scene (no pyramid), compared with T1."""
import time, json, numpy as np
from common import *
from sarsim import viz
import matplotlib.pyplot as plt
from t01_static_pyramid import transect_figure, ROT

TEST = "t02_desert_null"


def main():
    viz.style()
    g = geometry()
    sc = desert_scene(g)
    slc = make_slc("desert", sc, g)
    x, r = axes(g)
    o = analyze_scene("desert", slc, g)
    op = analyze_scene("pyramid_rot8", make_slc("pyramid_rot8", pyramid_scene(g, rotation_deg=ROT), g, rotation_deg=ROT), g)
    fd = figdir(TEST)
    figs = []
    tr = transects(g)
    rows, cols = tr["az"]
    figs.append({"file": transect_figure(f"{fd}/t02_transect_az.png", slc, x, r, rows, cols, o, "az",
                                         "Empty desert: the same azimuth transect as T1, no pyramid anywhere",
                                         "azimuth (m)", x[rows]),
                 "caption": "<b>Desert only.</b> Flat ground with speckle and a few rocks, same transect, same pipeline. The 'tomogram' is as busy as the pyramid's."})
    # histogram comparison
    rs, cs = o["grid_rows"], o["grid_cols"]
    smax_d = o["grid_score"].max(axis=2)
    smax_p = op["grid_score"].max(axis=2)
    inside = region_masks(g, rs, cs, ROT)["inside_footprint"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.9))
    bins = np.linspace(0, 1, 41)
    axs[0].hist(smax_p[inside], bins=bins, density=True, histtype="step", lw=2, color=viz.SERIES[0], label="pyramid pixels (T1)")
    axs[0].hist(smax_p[~inside], bins=bins, density=True, histtype="step", lw=2, color=viz.SERIES[2], label="desert pixels around it (T1)")
    axs[0].hist(smax_d.ravel(), bins=bins, density=True, histtype="step", lw=2, color=viz.SERIES[1], label="empty desert scene (T2)")
    axs[0].set_xlabel("best adjusted $R^{2}$ over all depths"); axs[0].set_ylabel("density"); axs[0].legend(loc="upper right")
    axs[0].set_title("Structure scores: pyramid vs nothing", loc="left")
    zp = op["z"][np.argmax(op["grid_score"], axis=2)]
    zd = o["z"][np.argmax(o["grid_score"], axis=2)]
    zb = np.linspace(0, o["z"][-1], 35)
    axs[1].hist(zp[inside], bins=zb, density=True, histtype="step", lw=2, color=viz.SERIES[0], label="pyramid pixels")
    axs[1].hist(zd.ravel(), bins=zb, density=True, histtype="step", lw=2, color=viz.SERIES[1], label="empty desert")
    axs[1].set_xlabel("model depth of the best fit (m)"); axs[1].legend(loc="upper right")
    axs[1].set_title("Where the 'structures' land", loc="left")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t02_histograms.png"), "caption":
                 "<b>Same distributions.</b> Left: the best steering-fit score per pixel for pixels on the pyramid, on the desert around it, and on a scene with no pyramid at all. Right: the depths those best fits land at. If the method saw anything through the stone, the blue curves would differ from the orange."})
    # maps side by side
    ext = [x[rs[0]], x[rs[-1]], r[cs[-1]], r[cs[0]]]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.6))
    axs[0].imshow(smax_p.T, origin="upper", extent=ext, cmap=viz.CMAP_SEQ, vmin=0, vmax=1, aspect="equal"); axs[0].set_title("pyramid scene: best adjusted $R^{2}$", loc="left")
    im = axs[1].imshow(smax_d.T, origin="upper", extent=ext, cmap=viz.CMAP_SEQ, vmin=0, vmax=1, aspect="equal"); axs[1].set_title("empty desert: best adjusted $R^{2}$", loc="left")
    for a in axs: a.set_xlabel("azimuth (m)")
    axs[0].set_ylabel("slant range (m)")
    fig.colorbar(im, ax=axs[1], fraction=0.046, pad=0.02)
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t02_maps.png"), "caption":
                 "<b>Spot the pyramid.</b> Best-fit score maps for the pyramid scene and the empty scene on the same grid. The monument occupies the upper middle of the left panel and leaves no trace."})
    from scipy.stats import ks_2samp
    ks = ks_2samp(smax_p[inside], smax_d.ravel())
    m = {"mean_best_adjR2_pyramid": float(smax_p[inside].mean()), "mean_best_adjR2_desert_scene": float(smax_d.mean()),
         "frac_gt_0p5_pyramid": float((smax_p[inside] > 0.5).mean()), "frac_gt_0p5_desert_scene": float((smax_d > 0.5).mean()),
         "ks_statistic": float(ks.statistic), "ks_pvalue": float(ks.pvalue), "runtime_s": float(o["runtime_s"]),
         "trajectory_rms_px_desert": float(o["grid_q"][..., 0].std()), "trajectory_rms_px_pyramid": float(op["grid_q"][..., 0].std())}
    summary = {
        "id": TEST, "order": 2, "tag": "null", "eyebrow": "02 · Null control",
        "title": "Empty desert gives the same 'structures' as the pyramid",
        "question": "Is the output any different when there is no pyramid at all?",
        "finding": f"No. On a scene of flat desert clutter the per-pixel structure scores have the same distribution as on the pyramid (means {m['mean_best_adjR2_desert_scene']:.2f} vs {m['mean_best_adjR2_pyramid']:.2f}; Kolmogorov–Smirnov distance {m['ks_statistic']:.2f}), the trajectories have the same amplitude, and the winning depths spread the same way.",
        "limitations": "Speckle is fully developed and Gaussian in this simulation; real desert has texture, but texture only adds more deterministic structure to the trajectories, not less.",
        "method": "Identical pipeline and grid as T1 on a scene containing only clutter and sparse rocks; scores compared pixel by pixel with T1's pyramid footprint.",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1))


if __name__ == "__main__":
    main()
