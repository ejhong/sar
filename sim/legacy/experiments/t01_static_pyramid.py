"""T1: a static, motionless stepped pyramid on flat ground put through the reconstructed pipeline."""
import time, json, numpy as np
from common import *
from sarsim import viz
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

TEST = "t01_static_pyramid"
ROT = 8.0


def transect_figure(path, slc, x, r, rows, cols, o, prefix, title, axis_label, coord):
    fig, axs = plt.subplots(3, 1, figsize=(10, 8.2), gridspec_kw={"height_ratios": [1, 2.2, 2.2]})
    amp = np.abs(slc[rows, cols])
    axs[0].plot(coord, amp, color=viz.SERIES[0], lw=0.9)
    axs[0].set_ylabel("|SLC|")
    axs[0].set_title(title, loc="left")
    axs[0].set_xlim(coord[0], coord[-1])
    tp = o[f"{prefix}_tomo_paper"]
    tpn = tp / tp.sum(axis=1, keepdims=True)
    z = o["z"]
    im = axs[1].imshow(tpn.T, aspect="auto", origin="upper", extent=[coord[0], coord[-1], z[-1], z[0]],
                       cmap=viz.CMAP_SEQ, vmin=0, vmax=np.percentile(tpn, 99.5), interpolation="nearest")
    axs[1].set_ylabel("model depth (m)")
    axs[1].set_title("paper variant: normalised $|a(z)^{H}Y|^{2}$ per pixel", loc="left")
    sw = o[f"{prefix}_score_w"]
    axs[2].imshow(sw.T, aspect="auto", origin="upper", extent=[coord[0], coord[-1], z[-1], z[0]],
                  cmap=viz.CMAP_SEQ, vmin=0, vmax=1, interpolation="nearest")
    axs[2].set_ylabel("model depth (m)")
    axs[2].set_xlabel(axis_label)
    axs[2].set_title("replication variant: max-over-windows adjusted $R^{2}$ (0 to 1)", loc="left")
    for a in axs[1:]:
        a.set_xlim(coord[0], coord[-1])
    fig.tight_layout(h_pad=1.2)
    return viz.finish(fig, path)


def main():
    viz.style()
    t0 = time.time()
    g = geometry()
    sc = pyramid_scene(g, rotation_deg=ROT)
    slc = make_slc("pyramid_rot8", sc, g, rotation_deg=ROT)
    x, r = axes(g)
    o = analyze_scene("pyramid_rot8", slc, g)
    fd = figdir(TEST)
    figs = []
    # 1. SLC with transects and grid box
    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    viz.imshow_slc(ax, slc, x, r)
    tr = transects(g)
    for name, (rows, cols), col in [("az", tr["az"], viz.SERIES[1]), ("rg", tr["rg"], viz.SERIES[2])]:
        ax.plot(x[rows], r[cols], color=col, lw=1.6)
    ax.add_patch(Rectangle((GRID_X[0], GRID_R[0]), GRID_X[1] - GRID_X[0], GRID_R[1] - GRID_R[0],
                           fill=False, ec=viz.SERIES[3], lw=1.2, ls="--"))
    ax.text(x[tr["az"][0][-1]] - 2, r[tr["az"][1][0]] - 3, "azimuth transect", color=viz.SERIES[1], ha="right", fontsize=8.5)
    ax.text(x[tr["rg"][0][0]] + 3, r[tr["rg"][1][0]] + 6, "range transect", color=viz.SERIES[2], fontsize=8.5)
    ax.set_title("Simulated SLC amplitude (dB): stepped pyramid, 136 m high, 1 m courses, heading offset 8°", loc="left")
    figs.append({"file": viz.finish(fig, f"{fd}/t01_slc.png"), "caption":
                 "<b>The input.</b> A simulated X-band single-look complex image at 0.25 m pixels of a motionless stepped pyramid on flat desert clutter. Layover folds the near face and apex toward the radar (up). Nothing exists below the ground plane. The orange and green lines are the transects analysed below; the dashed box is the 2-D grid."})
    # 2/3. transects
    rows, cols = tr["az"]
    figs.append({"file": transect_figure(f"{fd}/t01_transect_az.png", slc, x, r, rows, cols, o, "az",
                                         f"Azimuth transect at slant range {TRANSECT_AZ_R:.0f} m (crosses both side faces and the desert between)",
                                         "azimuth (m)", x[rows]),
                 "caption": "<b>Azimuth transect.</b> Top: image amplitude along the line (the two bright stretches are the side faces). Middle and bottom: the two published depth-focusing variants applied to the same 50-step trajectories. Both fill the whole depth axis with blobs, streaks and layered bands. The desert between the faces is as 'structured' as the pyramid."})
    rows, cols = tr["rg"]
    figs.append({"file": transect_figure(f"{fd}/t01_transect_rg.png", slc, x, r, rows, cols, o, "rg",
                                         f"Range transect at azimuth {TRANSECT_RG_X:.0f} m (apex layover, near face, ground)",
                                         "slant range (m)", r[cols]),
                 "caption": "<b>Range transect.</b> Same treatment down the middle of the pyramid. The laid-over apex and near face (left of about −55 m) and the desert (right) produce indistinguishable 'tomograms'."})
    # 4. maps
    rs, cs = o["grid_rows"], o["grid_cols"]
    S = o["grid_score"]
    smax = S.max(axis=2)
    zbest = o["z"][np.argmax(S, axis=2)]
    ext = [x[rs[0]], x[rs[-1]], r[cs[-1]], r[cs[0]]]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.9))
    sub = slc[rs[0]:rs[-1] + 1, cs[0]:cs[-1] + 1]
    viz.imshow_slc(axs[0], sub, x[rs[0]:rs[-1] + 1], r[cs[0]:cs[-1] + 1])
    axs[0].set_title("SLC amplitude (dB)", loc="left")
    im1 = axs[1].imshow(smax.T, origin="upper", extent=ext, cmap=viz.CMAP_SEQ, vmin=0, vmax=1, aspect="equal")
    axs[1].set_title("best adjusted $R^{2}$ over all depths", loc="left")
    im2 = axs[2].imshow(zbest.T, origin="upper", extent=ext, cmap=viz.CMAP_SEQ, vmin=0, vmax=o["z"][-1], aspect="equal")
    axs[2].set_title("model depth of the best fit (m)", loc="left")
    for a in axs[1:]:
        a.set_xlabel("azimuth (m)")
        a.set_yticklabels([])
    fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.02)
    fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.02)
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t01_maps.png"), "caption":
                 "<b>Maps over the whole pyramid at 1 m spacing.</b> Left: the image. Middle: for every pixel, the best steering-fit score at any depth. Right: the depth that wins. Neither map knows where the pyramid is: the score is as high on empty desert as on the monument, and the winning depths are scattered across the full axis."})
    # 5. example trajectories
    q = o["grid_q"].reshape(rs.size, cs.size, 50, 2)
    picks = {"near face (layover)": (0.0, -80.0), "side face": (-84.0, -10.0), "desert": (0.0, 45.0)}
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.7))
    for i, (name, (px, pr)) in enumerate(picks.items()):
        ii = int(np.argmin(np.abs(x[rs] - px))); jj = int(np.argmin(np.abs(r[cs] - pr)))
        t = q[ii, jj]
        axs[i].plot(t[:, 0], t[:, 1], "-", color=viz.SERIES[i], lw=1.2)
        axs[i].scatter(t[:, 0], t[:, 1], s=14, color=viz.SERIES[i], edgecolor=viz.SURFACE, linewidth=1, zorder=3)
        axs[i].set_title(f"{name}: x={px:.0f} m, r={pr:.0f} m", loc="left")
        axs[i].set_xlabel("azimuth shift (px)")
        axs[i].axhline(0, color=viz.GRID, lw=0.8); axs[i].axvline(0, color=viz.GRID, lw=0.8)
        lim = np.abs(t).max() * 1.15
        axs[i].set_xlim(-lim, lim); axs[i].set_ylim(-lim, lim); axs[i].set_aspect("equal")
    axs[0].set_ylabel("range shift (px)")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t01_trajectories.png"), "caption":
                 "<b>What a 'micro-motion trajectory' looks like.</b> The 50 reference/offset displacement vectors for one pixel on the near face, one on a side face and one on bare desert, in pixels (1 px = 0.25 m). Amplitudes are a few thousandths of a pixel, about a millimetre, and the loops look the same whether or not there is a pyramid under the patch. The scene never moved: these are deterministic properties of static scattering."})
    # metrics
    masks = region_masks(g, rs, cs, ROT)
    inside = masks["inside_footprint"]
    m = {"grid_pixels": int(smax.size), "runtime_s": float(o["runtime_s"]),
         "mean_best_adjR2_pyramid": float(smax[inside].mean()), "mean_best_adjR2_desert": float(smax[~inside].mean()),
         "frac_adjR2_gt_0p5_pyramid": float((smax[inside] > 0.5).mean()), "frac_adjR2_gt_0p5_desert": float((smax[~inside] > 0.5).mean()),
         "trajectory_rms_px": float(o["grid_q"][..., 0].std()), "trajectory_rms_mm": float(o["grid_q"][..., 0].std() * g.dx * 1e3),
         "scatterers": int(sc.n)}
    summary = {
        "id": TEST, "order": 1, "tag": "artifact", "eyebrow": "01 · Static scene",
        "title": "A motionless pyramid on solid ground produces a full 'tomogram'",
        "question": "If nothing moves and nothing is buried, does the method still draw structures at depth?",
        "finding": f"Yes. With every scatterer frozen and no subsurface at all, both published focusing variants fill the depth axis with blobs, streaks and layered bands. The best-fit score averages {m['mean_best_adjR2_pyramid']:.2f} on the pyramid and {m['mean_best_adjR2_desert']:.2f} on bare desert, and the winning depths scatter across the whole axis.",
        "limitations": "A simulation with an idealised, straight-line spotlight geometry and a stylised pyramid. The exact pattern of blobs depends on the random scene; the fact that a pattern appears does not.",
        "method": f"Scene: {sc.n:,} point scatterers (course edges with aspect-dependent 'flash', diffuse face scatter, rocks) plus speckle clutter, synthesised into a 1280 × 1024 SLC. Pipeline: K = 50 reference/offset sub-aperture pairs (half-band, 88 Hz offset snapped to 4 FFT bins), 32 × 32 complex cross-correlation to 1/1000 px, global common-mode removal, λₛ = 0.48 m. Grid: {smax.size:,} pixels at 1 m, run time {o['runtime_s']/60:.1f} min.",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1))


if __name__ == "__main__":
    main()
