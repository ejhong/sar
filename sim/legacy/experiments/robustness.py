"""Additional controls: phase variation, wavelength invariance, and ordering-test power."""
from pathlib import Path
import csv
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sarsim import Geometry, SubapBank, point_targets, synthesize, run_pipeline, focus_windows, viz
import matplotlib.pyplot as plt


def main():
    start = time.time()
    out = ROOT / "results" / "robustness"
    out.mkdir(parents=True, exist_ok=True)
    g = Geometry()
    bank = SubapBank()
    records = []
    reference = None
    for seed in range(10):
        phase = np.random.default_rng(seed).uniform(-np.pi, np.pi, 2)
        for separation in (1.0, 2.0, 3.0, 4.0):
            scene = point_targets([-separation / 2, separation / 2], [0, 0], [0, 0], [1, .8], phase=phase)
            slc = synthesize(scene, g, (1280, 128), ground_intensity=1e-8, seed=seed)
            result = run_pipeline(slc, g, np.array([640]), np.array([64]), bank=bank,
                                  common="none", verbose=False, modes=("paper",))
            measured = float(result["z"][np.argmax(result["tomo_paper"][0])])
            predicted = separation * .48 * np.sin(np.deg2rad(35)) / (299792458 / 9.65e9)
            records.append(dict(seed=seed, separation_m=separation, predicted_depth_m=float(predicted), peak_depth_m=measured))
            if seed == 0 and separation == 3:
                reference = (slc, result)
    with (out / "two_scatterer_seeds.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    slc, ref = reference
    invariance = []
    for wavelength in (.24, .48, .96):
        result = run_pipeline(slc, g, np.array([640]), np.array([64]), bank=bank, lam_s=wavelength,
                              common="none", verbose=False, modes=("paper",))
        scale_error = float(np.max(np.abs(result["z"] - ref["z"] * wavelength / .48)))
        power_error = float(np.max(np.abs(result["tomo_paper"] - ref["tomo_paper"])))
        assert scale_error < 1e-9 and power_error < 1e-12
        invariance.append(dict(wavelength_m=wavelength, max_axis_error_m=scale_error, max_power_error=power_error))
    rng = np.random.default_rng(832)
    kz = np.linspace(-.3, .3, 50)
    z = np.linspace(5, 180, 80)
    z0 = z[20]
    n = 96
    phase = rng.uniform(-np.pi, np.pi, (n, 1))
    wave = kz[None, :] * z0 + phase
    ordered = np.stack([.02 * np.cos(wave), .01 * np.sin(wave)], axis=2)
    ordered += rng.normal(0, .002, ordered.shape)
    noise = rng.normal(0, .02, ordered.shape)
    q = np.concatenate([ordered, noise])
    true = focus_windows(q, kz, z)[0].max(axis=1)
    exceed = np.zeros(2 * n, int)
    for _ in range(99):
        shuffled = focus_windows(q[:, rng.permutation(50)], kz, z)[0].max(axis=1)
        exceed += shuffled >= true
    p = (1 + exceed) / 100
    power = {"ordered_rejection_fraction": float((p[:n] < .05).mean()),
             "noise_rejection_fraction": float((p[n:] < .05).mean()),
             "trajectories_per_class": n, "permutations": 99, "nominal_level": .05,
             "scope": "Synthetic noisy steering sinusoids and independent Gaussian samples; not radar-derived trajectories."}
    assert power["ordered_rejection_fraction"] > .9
    viz.style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    for seed in range(10):
        rows = [r for r in records if r["seed"] == seed]
        axes[0].scatter([r["separation_m"] + (seed - 4.5) * .015 for r in rows],
                        [r["peak_depth_m"] for r in rows], s=18, color="#26756b", alpha=.55)
    axes[0].plot([1, 4], [records[0]["predicted_depth_m"], records[3]["predicted_depth_m"]], color="#9a5b33", label="Spacing prediction")
    axes[0].set(xlabel="Azimuth separation (m)", ylabel="Peak model depth (m)", title="Ten phase realizations per separation")
    axes[0].legend()
    axes[1].bar(["Planted ordering", "Independent noise"], [power["ordered_rejection_fraction"], power["noise_rejection_fraction"]], color=["#26756b", "#9a5b33"], width=.5)
    axes[1].axhline(.05, color="#777", ls=":", label="Nominal 5%")
    axes[1].set(ylabel="Fraction with permutation p < 0.05", ylim=(0, 1.08), title="Can the ordering test detect a positive control?")
    axes[1].legend(loc="upper right")
    fig.savefig(out / "robustness.png"); plt.close(fig)
    errors = np.array([abs(r["peak_depth_m"] - r["predicted_depth_m"]) for r in records])
    report = {"date": time.strftime("%Y-%m-%d"), "seed_runs": records, "wavelength_invariance": invariance,
              "ordering_control": power, "median_peak_error_m": float(np.median(errors)),
              "max_peak_error_m": float(errors.max()), "runs": len(records),
              "limitations": "Ten phase realizations of the same two-point scene, with negligible clutter. Not a multi-scene estimate of cavity-detection performance. The ordering positive control is generated in trajectory space.",
              "runtime_s": time.time() - start}
    (out / "metrics.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    print(json.dumps({k: v for k, v in report.items() if k != "seed_runs"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
