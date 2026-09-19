"""T6: are static trajectories 'ordered', and does the replication protocol's ellipse gate pass motionless ground?"""
import time, json, numpy as np
from common import *
from sarsim import viz
from sarsim.tomo import focus_windows
import matplotlib.pyplot as plt

TEST = "t06_ordering"
N_PIX = 1500
N_PERM = 60
W = 25
MODES = range(1, 7)
R2_MIN, BA_MIN, B_MIN = 0.25, 0.1, 0.005     # protocol v1.5 branch-B gate


def ellipse_gate(q):
    """Best-of-(windows × modes) ellipse fit q_n = c + u cos + v sin, as in the replication protocol."""
    P, K, _ = q.shape
    nW = K - W + 1
    n = np.arange(W)
    best_r2 = np.full(P, -np.inf); best_ba = np.zeros(P); best_b = np.zeros(P)
    for w in range(nW):
        Qw = q[:, w:w + W, :]
        for m in MODES:
            X = np.stack([np.ones(W), np.cos(2 * np.pi * m * n / W), np.sin(2 * np.pi * m * n / W)], 1)
            Pinv = np.linalg.pinv(X)
            B = np.einsum('ik,pkj->pij', Pinv, Qw)
            fit = np.einsum('ki,pij->pkj', X, B)
            ssres = ((Qw - fit) ** 2).sum(axis=(1, 2))
            sstot = ((Qw - Qw.mean(axis=1, keepdims=True)) ** 2).sum(axis=(1, 2))
            r2 = 1 - ssres / np.maximum(sstot, 1e-30)
            adj = 1 - (1 - r2) * (2 * W - 1) / (2 * W - 6 - 1)
            sv = np.linalg.svd(np.transpose(B[:, 1:, :], (0, 2, 1)), compute_uv=False)
            better = adj > best_r2
            best_r2 = np.where(better, adj, best_r2)
            best_ba = np.where(better, sv[:, 1] / np.maximum(sv[:, 0], 1e-30), best_ba)
            best_b = np.where(better, sv[:, 1], best_b)
    shape = (best_r2 >= R2_MIN) & (best_ba >= BA_MIN)
    return best_r2, best_ba, best_b, shape, shape & (best_b >= B_MIN)


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
    axs[1].axhline(0.05, color=viz.SERIES[7], lw=1, ls="--"); axs[1].text(1.3, 0.09, "5% = chance", color=viz.SERIES[7], fontsize=8.5, ha="right")
    axs[1].set_ylim(0, 0.5); axs[1].set_ylabel("fraction of pixels 'ordered' (p < 0.05)")
    axs[1].set_title("Rejections near the nominal 5% rate", loc="left")
    for i, v in enumerate(fr):
        axs[1].text(i, v + 0.012, f"{v:.1%}", ha="center", fontsize=9, color=viz.INK2)
    gp = [float(d["shape"].mean()) for d in res.values()]; gps = [float(d["shape_shuffled"].mean()) for d in res.values()]
    xx = np.arange(2)
    axs[2].bar(xx - 0.17, gp, width=0.3, color=[viz.SERIES[0], viz.SERIES[1]], label="true order")
    axs[2].bar(xx + 0.17, gps, width=0.3, color=[viz.SERIES[0], viz.SERIES[1]], alpha=0.45, label="shuffled")
    axs[2].set_xticks(xx); axs[2].set_xticklabels(list(res.keys())); axs[2].set_ylim(0, 0.5)
    axs[2].set_ylabel("fraction passing the ellipse gate"); axs[2].set_title("Shape-only gate: about 1 in 10 passes", loc="left")
    for i in range(2):
        axs[2].text(i - 0.17, gp[i] + 0.012, f"{gp[i]:.0%}", ha="center", fontsize=9, color=viz.INK2)
        axs[2].text(i + 0.17, gps[i] + 0.012, f"{gps[i]:.0%}", ha="center", fontsize=9, color=viz.INK2)
    axs[2].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t06_permutation.png"), "caption":
                 "<b>Shape and amplitude are separate checks.</b> Left and middle: a permutation test asks whether the structure score depends on the true sub-aperture order. On both simulated scenes, rejection rates are close to the nominal level for this statistic. Right: the replication protocol's branch-B gate (best ellipse fit over 26 windows × 6 modes with adjusted R² ≥ 0.25 and axis ratio ≥ 0.1) still passes about a tenth of motionless pixels, shuffled or not, after searching 156 window/mode combinations. The protocol's absolute minor-axis floor of 0.005 px is met by none of these pixels; therefore the full gate rejects all sampled pixels in this run."})
    m = {"frac_perm_significant": dict(zip(res.keys(), fr)), "ellipse_shape_gate_pass": dict(zip(res.keys(), gp)),
         "ellipse_shape_gate_pass_shuffled": dict(zip(res.keys(), gps)),
         "ellipse_full_gate_pass": {k: float(d["full"].mean()) for k, d in res.items()},
         "median_minor_axis_px": {k: float(np.median(d["b"])) for k, d in res.items()},
         "median_best_ellipse_adjR2": {k: float(np.median(d["r2"])) for k, d in res.items()},
         "n_pixels": N_PIX, "n_permutations": N_PERM, "runtime_s": time.time() - t0}
    summary = {
        "id": TEST, "order": 6, "tag": "null", "eyebrow": "06 · Ordering & gates",
        "title": "The shape gate passes some static pixels; the full gate rejects them",
        "question": "The replication protocol keeps only pixels whose trajectories trace ordered elliptical loops. Does motionless ground pass, and is there any ordering to find?",
        "finding": f"This statistic finds little evidence of ordering: shuffling the 50 sub-apertures leaves the structure scores unchanged, with {fr[0]:.1%} of pyramid pixels and {fr[1]:.1%} of desert pixels 'significant' at the 5% level, i.e. chance. The protocol's shape gate nevertheless passes {gp[0]:.0%} of pyramid and {gp[1]:.0%} of desert pixels, and just as many after shuffling. The implemented full gate, including the 0.005 px minor-axis floor, accepts none of the sampled pixels. Shape-only acceptance must not be reported as acceptance by the full protocol. Real Giza trajectories, if they turn out to be ordered, would be showing something this simulation does not contain, and the first question would be whether empty plateau shows it too.",
        "limitations": "Sixty permutations on 1,500 pixels from one realization per scene; spatially dependent samples share each permutation. A near-chance rejection rate is not proof that every form of ordering is absent. The full gate accepts zero pixels at the tested amplitude scale.",
        "method": "Per-pixel best score in true order versus 60 random orderings, plus the protocol's branch-B ellipse fit (26 windows × modes 1 to 6, adjusted R² with 6 parameters, axes from the fitted cos/sin vectors) on 1,500 pixels on the pyramid footprint (T1) and 1,500 on the empty desert (T2).",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1, default=float))


if __name__ == "__main__":
    main()
