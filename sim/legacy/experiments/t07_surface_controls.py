"""T7: paired surface-target controls, with fixed processing and unsmoothed metrics.

The design is saved before any output is inspected. This is a sensitivity study
of ideal surface targets, not a reconstruction of Khafre's observed scene.
"""
from pathlib import Path
import csv
import json
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c
from sarsim import concat, point_targets, run_pipeline, viz
from sarsim.scene import rot2d
from t06_ordering import ellipse_gate, R2_MIN, BA_MIN, B_MIN
import matplotlib.pyplot as plt

TEST = "t07_surface_controls"
ROTATION = 8.0
STRENGTHS = (1.0, 4.0, 8.0)
SEEDS = (101, 102, 103)
STRIDE = 10  # 2.5 m; every scene uses the same full spatial grid
PATCH = 32
WAVELENGTH = 0.48
HALF_PROBE = 4.0  # fixed 8 m square in image coordinates
X_LIMITS = (-135, 135)
R_LIMITS = (-125, 85)


def layouts(g):
    """Original arrangement and three seeded, separated perimeter layouts.

    Random relative amplitudes are uniform on [0.5, 1.5]. The same positions,
    relative amplitudes and phases are reused at all three strength settings.
    Added targets are ideal isotropic ground points, without terrain shadowing.
    """
    half = c.PYRAMID["base"] / 2
    x = np.array([-half, half, half, -half, 0, half, 0, -half])
    y = np.array([-half, -half, half, half, -half, 0, half, 0])
    x, y = rot2d(x, y, ROTATION)
    phase = np.random.default_rng(3).uniform(0, 2 * np.pi, 8)
    result = [{"id": "original", "label": "Original layout", "seed": None,
               "x": x, "y": y, "relative_amp": np.ones(8), "phase": phase}]
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        positions = []
        for _ in range(10000):
            side = int(rng.integers(4))
            along = rng.uniform(-half, half)
            px, py = ((along, -half), (half, along), (along, half), (-half, along))[side]
            px, py = rot2d(px, py, ROTATION)
            if all(np.hypot(px - qx, (py - qy) * np.sin(g.theta)) >= 12 for qx, qy in positions):
                positions.append((float(px), float(py)))
            if len(positions) == 8:
                break
        if len(positions) != 8:
            raise RuntimeError("Could not place separated perimeter targets")
        x, y = np.array(positions).T
        result.append({"id": f"random_{seed}", "label": f"Random layout {seed}", "seed": seed,
                       "x": x, "y": y, "relative_amp": rng.uniform(0.5, 1.5, 8), "phase": phase.copy()})
    return result


def paired_region_metrics(power, baseline, mask):
    """Compare regional mean raw depth spectra at exactly the same pixels.

    Equality is not suppression. Neither spectrum is normalized by its own
    maximum, and there is no Gaussian smoothing, thresholding or window search.
    """
    if not np.any(mask):
        raise ValueError("The paired region contains no grid targets")
    current = power[mask].mean(axis=0)
    reference = baseline[mask].mean(axis=0)
    if reference.sum() <= 0:
        raise ValueError("The baseline has no power in the paired region")
    return {
        "sampled_pixels": int(mask.sum()),
        "mean_raw_power": float(current.mean()),
        "baseline_mean_raw_power": float(reference.mean()),
        "power_ratio_to_baseline": float(current.sum() / reference.sum()),
        "fraction_depth_bins_below_baseline": float(np.mean(current < reference)),
    }


def run_scene(slc, g, bank, rows, cols):
    path = Path(c.cache_path("t07_surface_control", slc=c.array_fingerprint(slc),
                            geometry=g.__dict__, bank=bank.__dict__, rows=c.array_fingerprint(rows),
                            cols=c.array_fingerprint(cols), patch=PATCH, wavelength=WAVELENGTH)).with_suffix(".npz")
    if path.exists():
        with np.load(path, allow_pickle=False) as saved:
            return {key: saved[key] for key in saved.files}
    rr, cc = np.meshgrid(rows, cols, indexing="ij")
    output = run_pipeline(slc, g, rr.ravel(), cc.ravel(), bank=bank, patch=PATCH,
                          lam_s=WAVELENGTH, common="global", modes=("paper",), verbose=True)
    result = {"power": output["tomo_paper"].reshape(rows.size, cols.size, -1),
              "q": output["q"], "z": output["z"]}
    np.savez(path, **result)
    return result


def gate_counts(q):
    shape = full = 0
    for start in range(0, len(q), 1024):
        *_, accepted_shape, accepted_full = ellipse_gate(q[start:start + 1024])
        shape += int(accepted_shape.sum())
        full += int(accepted_full.sum())
    return {"shape_gate_count": shape, "full_gate_count": full,
            "shape_gate_fraction": shape / len(q), "full_gate_fraction": full / len(q)}


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def figures(out, g, xs, rs, z, samples, runs, all_layouts, reference_slc):
    """Use the first seeded random layout, never a result-selected example."""
    fd = out / "figs"
    fd.mkdir(exist_ok=True)
    base = samples["baseline"]["power"]
    plan_reference = float(np.median(base.sum(axis=2)))
    depth_reference = float(np.median(base))
    amplitude_reference = float(np.abs(reference_slc).max())
    section_r = float(all_layouts[0]["y"][4] * np.sin(g.theta))
    section_index = int(np.argmin(np.abs(rs - section_r)))
    extent = [xs[0], xs[-1], rs[-1], rs[0]]
    chosen = ("baseline", "original_8", "random_101_8")
    names = ("Unchanged pyramid", "Original layout · strength 8", "First random layout · strength 8")
    fig, axes = plt.subplots(3, 3, figsize=(13, 11), layout="constrained")
    images = []
    for row, (key, name) in enumerate(zip(chosen, names)):
        sample = samples[key]
        power = sample["power"]
        amplitude = 20 * np.log10(np.maximum(np.abs(sample["slc"]) / amplitude_reference, 1e-12))
        total = 10 * np.log10(np.maximum(power.sum(axis=2) / plan_reference, 1e-12))
        section = 10 * np.log10(np.maximum(power[:, section_index, :] / depth_reference, 1e-12))
        images = [axes[row, 0].imshow(amplitude.T, extent=extent, origin="upper", aspect="equal",
                                     cmap=viz.CMAP_GRAY, vmin=-50, vmax=5, interpolation="nearest"),
                  axes[row, 1].imshow(total.T, extent=extent, origin="upper", aspect="equal",
                                     cmap="viridis", vmin=-15, vmax=15, interpolation="nearest"),
                  axes[row, 2].imshow(section.T, extent=[xs[0], xs[-1], z[-1], z[0]], origin="upper",
                                     aspect="auto", cmap="viridis", vmin=-25, vmax=25, interpolation="nearest")]
        for column, title in enumerate(("Input amplitude", "Total unsmoothed power", f"Fixed section at range {rs[section_index]:.1f} m")):
            axes[row, column].set_title(f"{name}\n{title}", loc="left", fontsize=9)
            axes[row, column].set_xlabel("azimuth (m)")
            axes[row, column].set_ylabel("slant range (m)" if column < 2 else "model depth (m), λs = 0.48 m")
        axes[row, 1].axhline(rs[section_index], color="white", lw=0.7, ls="--")
        if sample["layout"] is not None:
            layout = sample["layout"]
            for ax in axes[row, :2]:
                ax.scatter(layout["x"], layout["y"] * np.sin(g.theta), s=25, facecolors="none", edgecolors="white", lw=0.8)
    for column, label in enumerate(("dB relative to unchanged image maximum", "dB relative to baseline median total power", "dB relative to baseline median voxel power")):
        fig.colorbar(images[column], ax=axes[:, column], orientation="horizontal", shrink=0.9, pad=0.04, label=label)
    comparison = viz.finish(fig, str(fd / "t07_controlled_maps.png"))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    for index, layout in enumerate(all_layouts):
        selected = [row for row in runs if row["layout"] == layout["id"]]
        x = [row["strength"] for row in selected]
        axes[0].plot(x, [row["median_point_power_ratio"] for row in selected], "o-", color=viz.SERIES[index], label=layout["label"])
        axes[1].plot(x, [row["median_point_depth_suppression_fraction"] for row in selected], "o-", color=viz.SERIES[index], label=layout["label"])
    axes[0].axhline(1, color=viz.MUTED, lw=1, ls="--")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("median raw power / matched baseline")
    axes[0].set_title("Below 1: less power than the same locations without added points", loc="left", fontsize=9)
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_ylabel("median fraction of depth bins below baseline")
    axes[1].set_title("Suppression across depth; this is not a column detector", loc="left", fontsize=9)
    for ax in axes:
        ax.set_xticks(STRENGTHS)
        ax.set_xlabel("prescribed point-amplitude multiplier")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.35)
    metrics = viz.finish(fig, str(fd / "t07_paired_metrics.png"))
    return [
        {"file": metrics, "alt": "Paired power ratios and depth-bin suppression at three strengths for four prescribed point layouts", "caption": "<b>Paired measurements before smoothing.</b> Each added-point region is compared with exactly the same pixels in the unchanged pyramid. Lines connect three prescribed strengths for one layout; each value is a median over eight ideal targets. The three random layouts use seeds 101–103 and fixed random relative amplitudes. Shared pixels and one background do not provide independent field trials. All per-point measurements are downloadable."},
        {"file": comparison, "alt": "Unchanged pyramid, original point layout and first random point layout, with shared scales and unsmoothed depth sections", "caption": "<b>Same scene, same settings.</b> Unchanged pyramid, original eight-point layout, and the first seeded random layout. Each column shares one baseline reference and fixed colour limits. All panels use nearest-neighbour display with no smoothing; the vertical section is fixed from the original front-midpoint geometry before examining output. Circled points are deliberately added synthetic targets. Brightness is reconstruction power, not acceptance by a detection gate."},
    ], {"amplitude_reference": amplitude_reference, "plan_power_reference": plan_reference,
        "depth_power_reference": depth_reference, "fixed_section_slant_range_m": float(rs[section_index]),
        "amplitude_limits_db": [-50, 5], "plan_limits_db": [-15, 15], "depth_limits_db": [-25, 25]}


def main():
    viz.style()
    started = time.time()
    out = Path(c.RESULTS) / TEST
    out.mkdir(exist_ok=True)
    g, bank = c.geometry(), c.default_bank()
    all_layouts = layouts(g)
    x, r = c.axes(g)
    rows = np.arange(np.searchsorted(x, X_LIMITS[0]), np.searchsorted(x, X_LIMITS[1]), STRIDE)
    cols = np.arange(np.searchsorted(r, R_LIMITS[0]), np.searchsorted(r, R_LIMITS[1]), STRIDE)
    xs, rs = x[rows], r[cols]
    X, R = np.meshgrid(xs, rs, indexing="ij")
    far = np.ones(X.shape, bool)
    targets = []
    for layout in all_layouts:
        for i, (px, py, amp, phase) in enumerate(zip(layout["x"], layout["y"], layout["relative_amp"], layout["phase"])):
            image_r = py * np.sin(g.theta)
            far &= np.hypot(X - px, R - image_r) > 12
            targets.append({"layout": layout["id"], "seed": layout["seed"], "point": i,
                            "x_m": px, "y_m": py, "z_m": 0.0, "slant_range_m": image_r,
                            "relative_amplitude": amp, "phase_rad": phase})
    write_csv(out / "point_design.csv", targets)
    design = {"scene_seed": 0, "rotation_deg": ROTATION, "strengths": list(STRENGTHS),
              "random_seeds": list(SEEDS), "n_points_per_layout": 8, "grid_spacing_m": STRIDE * g.dx,
              "random_layout_rule": "Uniform choice of base side and uniform position along it; reject image separations below 12 m",
              "minimum_random_point_image_separation_m": 12.0, "point_phase_seed": 3,
              "sampled_pixels_per_scene": int(X.size), "probe_width_m": 2 * HALF_PROBE,
              "remote_region": "More than 12 image-metres from every target in every layout",
              "remote_pixels": int(far.sum()), "wavelength_m": WAVELENGTH, "patch_pixels": PATCH,
              "common_mode": "Same rule in every scene: subtract the median trajectory of the full fixed grid",
              "smoothing": "none", "relative_amplitude_distribution": "uniform [0.5, 1.5] for random layouts; 1 for original",
              "added_target_model": "Ideal isotropic points at z=0; no shadowing or object-specific scattering model",
              "gate": {"adjusted_r2_min": R2_MIN, "axis_ratio_min": BA_MIN, "minor_axis_min_px": B_MIN,
                       "window_samples": 25, "modes": [1, 2, 3, 4, 5, 6]},
              "representative_panels": ["baseline", "original_8", "random_101_8"]}
    (out / "design.json").write_text(json.dumps(design, indent=2, allow_nan=False))
    scene = c.pyramid_scene(g, rotation_deg=ROTATION, seed=0)
    background_hash = {name: c.array_fingerprint(array) for name, array in vars(scene).items()}
    slc = c.make_slc("pyramid_rot8", scene, g, rotation_deg=ROTATION)
    print("[1/13] Unchanged pyramid", flush=True)
    base = run_scene(slc, g, bank, rows, cols)
    runs = [{"run": "baseline", "layout": "none", "strength": 0.0, "sampled_pixels": int(X.size),
             "median_point_power_ratio": None, "minimum_point_power_ratio": None, "maximum_point_power_ratio": None,
             "median_point_depth_suppression_fraction": None, "remote_power_ratio": 1.0, **gate_counts(base["q"])}]
    image_slice = np.s_[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
    samples = {"baseline": {"power": base["power"], "slc": slc[image_slice], "layout": None}}
    point_rows = []
    counter = 1
    for layout in all_layouts:
        for strength in STRENGTHS:
            counter += 1
            name = f"{layout['id']}_{strength:g}"
            print(f"[{counter}/13] {name}", flush=True)
            points = point_targets(layout["x"], layout["y"], np.zeros(8), strength * layout["relative_amp"],
                                   phase=layout["phase"], label=50)
            image = c.make_slc("surface_control", concat(scene, points), g, rotation_deg=ROTATION)
            result = run_scene(image, g, bank, rows, cols)
            np.testing.assert_array_equal(result["z"], base["z"])
            per_point = []
            for i, (px, py) in enumerate(zip(layout["x"], layout["y"])):
                mask = (np.abs(X - px) <= HALF_PROBE) & (np.abs(R - py * np.sin(g.theta)) <= HALF_PROBE)
                measurement = paired_region_metrics(result["power"], base["power"], mask)
                record = {"run": name, "layout": layout["id"], "strength": strength, "point": i,
                          "point_amplitude": float(strength * layout["relative_amp"][i]), **measurement}
                per_point.append(record)
                point_rows.append(record)
            ratios = [row["power_ratio_to_baseline"] for row in per_point]
            run = {"run": name, "layout": layout["id"], "strength": strength, "sampled_pixels": int(X.size),
                   "median_point_power_ratio": float(np.median(ratios)),
                   "minimum_point_power_ratio": min(ratios), "maximum_point_power_ratio": max(ratios),
                   "median_point_depth_suppression_fraction": float(np.median([row["fraction_depth_bins_below_baseline"] for row in per_point])),
                   "remote_power_ratio": paired_region_metrics(result["power"], base["power"], far)["power_ratio_to_baseline"],
                   **gate_counts(result["q"])}
            runs.append(run)
            if name in design["representative_panels"]:
                samples[name] = {"power": result["power"], "slc": image[image_slice], "layout": layout}
            print(f"  paired power {run['median_point_power_ratio']:.3f}; full gate {run['full_gate_count']}/{X.size}; elapsed {time.time()-started:.0f}s", flush=True)
            write_csv(out / "scene_runs.csv", runs)
            write_csv(out / "point_measurements.csv", point_rows)
    assert background_hash == {name: c.array_fingerprint(array) for name, array in vars(scene).items()}
    figs, display = figures(out, g, xs, rs, base["z"], samples, runs, all_layouts, slc)
    strong = [run for run in runs if run["strength"] == 8]
    original = next(run for run in strong if run["layout"] == "original")
    random = [run for run in strong if run["layout"] != "original"]
    total_full = sum(run["full_gate_count"] for run in runs)
    gate_finding = ("The implemented full gate accepts no sampled pixels in any of the 13 scenes."
                    if total_full == 0 else f"The implemented full gate accepts {total_full} pixel-scene samples across the 13 scenes; these are dependent samples, not independent detections.")
    summary = {"id": TEST, "order": 7, "tag": "controlled", "eyebrow": "07 · Paired surface controls",
               "title": "Add surface targets; measure what changes before rendering",
               "question": "How much of the column-like output depends on deliberately inserted surface points, and does it pass the implemented acceptance gate?",
               "finding": f"At strength 8, median unsmoothed power near the original eight points is {original['median_point_power_ratio']:.2f} times the same locations in the unchanged pyramid. Across the three random layouts it ranges from {min(row['median_point_power_ratio'] for row in random):.2f} to {max(row['median_point_power_ratio'] for row in random):.2f}. {gate_finding}",
               "limitations": "One background realization, three random perimeter layouts and ideal point targets with prescribed amplitudes. Added points have no terrain-shadowing or object-specific scattering model. Spatial samples overlap, and strength variants reuse the same layouts. Power changes and depth-bin suppression are descriptive responses, not validated column detections. Neither the point count nor the output depth limit predicts a Khafre structure. This tests our reconstruction and T6 gate implementation, not the original data or complete Khafre workflow.",
               "method": "Thirteen paired scenes: unchanged pyramid; original corner/midpoint layout at multipliers 1, 4 and 8; three random perimeter layouts at those same multipliers, with seeded relative amplitudes uniform on [0.5, 1.5]. All reuse the same pyramid, clutter seed, point phases, radar geometry, 50 sub-aperture pairs, 32-pixel patches, 2.5 m target grid, common-mode rule, 0.48 m wavelength and 160 model-depth bins. Regional raw power is compared with the identical baseline pixels in fixed 8 m squares. The T6 shape and full gates are applied to every grid trajectory. No smoothing, isosurfaces, per-scene normalization or depth rescaling is used. The displayed random example is always seed 101.",
               "figures": figs, "metrics": {"design": design, "runs": runs, "display": display,
                                            "total_full_gate_pixel_scene_samples": total_full,
                                            "depth_axis_max_m": float(base["z"][-1]),
                                            "background_unchanged": True, "runtime_s": time.time() - started},
               "downloads": ["scene_runs.csv", "point_measurements.csv", "point_design.csv", "design.json"],
               "date": time.strftime("%Y-%m-%d")}
    c.write_summary(TEST, summary)
    print(summary["finding"], flush=True)


if __name__ == "__main__":
    main()
