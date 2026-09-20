"""T6: are static trajectories 'ordered', and does the replication protocol's ellipse gate pass motionless ground?"""
import time, json, numpy as np
from common import *
from sarsim import viz
from sarsim.tomo import focus_windows
from sarsim.gates import ellipse_gate, R2_MIN, BA_MIN, B_MIN
import matplotlib.pyplot as plt

TEST = "t06_ordering"
N_PIX = 1500
N_PERM = 60
W = 25
MODES = range(1, 7)
def main():
    viz.style()
    t0 = time.time()
    g = geometry(); fd = figdir(TEST); figs = []
    rng = np.random.default_rng(0)
    res = {}
    for tag, name in [("pyramid_rot8", "pyramid pixels"), ("desert", "empty desert")]:
        scene = pyramid_scene(g, rotation_deg=8.0) if tag.startswith("pyr") else desert_scene(g)
        slc = make_slc(tag, scene, g, **({"rotation_deg": 8.0} if tag.startswith("pyr") else {}))
        o = analyze_scene(tag, slc, g)
        q = o["grid_q"]; z = o["z"]; kz = o["kz"]
        if tag.startswith("pyr"):
            inside = region_masks(g, o["grid_rows"], o["grid_cols"], 8.0)["inside_footprint"].ravel()
            pool = np.where(inside)[0]
        else:
            pool = np.arange(q.shape[0])
        pick = rng.choice(pool, N_PIX, replace=False)
        qs = q[pick]
        true = focus_windows(qs, kz, z)[0].max(axis=1)
        null = np.zeros((N_PERM, N_PIX))
        for p in range(N_PERM):
            null[p] = focus_windows(qs[:, rng.permutation(q.shape[1]), :], kz, z)[0].max(axis=1)
        pval = (1 + (null >= true[None, :]).sum(axis=0)) / (N_PERM + 1)
        r2, ba, b, shape, full = ellipse_gate(qs)
        r2s, bas, bs, shapes, fulls = ellipse_gate(qs[:, rng.permutation(q.shape[1]), :])
        res[name] = {"true": true, "null": null, "pval": pval, "r2": r2, "ba": ba, "b": b, "shape": shape, "full": full,
                     "shape_shuffled": shapes, "r2_shuffled": r2s}
        print(f"  {name}: {time.time()-t0:.0f}s  perm-significant {np.mean(pval < 0.05):.3f}  gate pass {shape.mean():.3f} (shuffled {shapes.mean():.3f})", flush=True)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.0))
    bins = np.linspace(0, 1, 41)
    for i, (name, d) in enumerate(res.items()):
        axs[0].hist(d["true"], bins=bins, density=True, histtype="step", lw=2, color=viz.SERIES[i], label=f"{name}, true order")
        axs[0].hist(d["null"].ravel(), bins=bins, density=True, histtype="step", lw=2, ls="--", color=viz.SERIES[i], label=f"{name}, shuffled")
    axs[0].set_xlabel(r"best adjusted $R^{2}$ over all depths"); axs[0].set_ylabel("density"); axs[0].legend(fontsize=8)
    axs[0].set_title("True and shuffled orders give similar score distributions", loc="left")
    fr = [float((d["pval"] < 0.05).mean()) for d in res.values()]
    axs[1].bar(range(2), fr, color=[viz.SERIES[0], viz.SERIES[1]], width=0.55)
    axs[1].set_xticks(range(2)); axs[1].set_xticklabels(list(res.keys()))
    axs[1].axhline(0.05, color=viz.SERIES[7], lw=1, ls="--"); axs[1].text(1.3, 0.09, "5% reference", color=viz.SERIES[7], fontsize=8.5, ha="right")
    axs[1].set_ylim(0, 0.5); axs[1].set_ylabel("fraction below heuristic shuffle threshold")
    axs[1].set_title("Historical shuffle diagnostic only", loc="left")
    for i, v in enumerate(fr):
        axs[1].text(i, v + 0.012, f"{v:.1%}", ha="center", fontsize=9, color=viz.INK2)
    gp = [float(d["shape"].mean()) for d in res.values()]; gps = [float(d["shape_shuffled"].mean()) for d in res.values()]
    xx = np.arange(2)
    axs[2].bar(xx - 0.17, gp, width=0.3, color=[viz.SERIES[0], viz.SERIES[1]], label="true order")
    axs[2].bar(xx + 0.17, gps, width=0.3, color=[viz.SERIES[0], viz.SERIES[1]], alpha=0.45, label="shuffled")
    axs[2].set_xticks(xx); axs[2].set_xticklabels(list(res.keys())); axs[2].set_ylim(0, 0.5)
    axs[2].set_ylabel("fraction passing the shape-only gate"); axs[2].set_title("Per-window shape selection, before amplitude floor", loc="left")
    for i in range(2):
        axs[2].text(i - 0.17, gp[i] + 0.012, f"{gp[i]:.0%}", ha="center", fontsize=9, color=viz.INK2)
        axs[2].text(i + 0.17, gps[i] + 0.012, f"{gps[i]:.0%}", ha="center", fontsize=9, color=viz.INK2)
    axs[2].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t06_permutation.png"), "caption":
                 "<b>Shape and amplitude are separate checks.</b> The shape-only panel omits the 0.005 px minor-axis floor; no sampled pixel passes the complete rule. Each window is selected and gated separately. Shuffling overlapping looks changes their covariance, so the two shuffle panels are historical diagnostics, not a calibrated physical false-alarm test."})
    m = {"frac_perm_significant": dict(zip(res.keys(), fr)), "ellipse_shape_gate_pass": dict(zip(res.keys(), gp)),
         "ellipse_shape_gate_pass_shuffled": dict(zip(res.keys(), gps)),
         "ellipse_full_gate_pass": {k: float(d["full"].mean()) for k, d in res.items()},
         "median_minor_axis_px": {k: float(np.median(d["b"])) for k, d in res.items()},
         "median_best_ellipse_adjR2": {k: float(np.median(d["r2"])) for k, d in res.items()},
         "n_pixels": N_PIX, "n_permutations": N_PERM, "runtime_s": time.time() - t0}
    summary = {
        "id": TEST, "order": 6, "tag": "null", "eyebrow": "06 · Ordering & gates",
        "title": "The shape gate passes some static pixels; the full gate rejects them",
        "question": "Does motionless ground pass the complete per-window ellipse filter, including its amplitude floor?",
        "finding": f"Shape-only acceptance is {gp[0]:.1%} and {gp[1]:.1%}; full acceptance is {m['ellipse_full_gate_pass']['pyramid pixels']:.1%} and {m['ellipse_full_gate_pass']['empty desert']:.1%}.",
        "limitations": "One realization per scene and dependent spatial samples. Look shuffles do not preserve the physical null covariance. Independent scene controls and physical positives are reported in the new validation suite.",
        "method": "Best mode 1–6 per W25 window, then shape and amplitude gating; any accepted window admits a pixel. A retained historical diagnostic uses 60 look permutations on 1,500 pixels per scene.",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1, default=float))


if __name__ == "__main__":
    main()
