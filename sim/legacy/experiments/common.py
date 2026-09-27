"""Shared configuration, caching and reporting for the experiments."""
import json, os, sys, time, hashlib
import numpy as np
from pathlib import Path
from dataclasses import asdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from sarsim import Geometry, stepped_pyramid, desert_rocks, concat, synthesize, SubapBank  # noqa: E402
from sarsim.synth import pixel_axes  # noqa: E402

RESULTS = os.path.join(ROOT, "results")
CACHE = os.path.join(RESULTS, "cache")
os.makedirs(CACHE, exist_ok=True)

# scene grid: azimuth x range pixels at 0.25 m
SHAPE = (1280, 1024)
R_CENTRE = -18.0          # slant-range grid centre (m): pyramid lays over toward the radar
GROUND_I = 0.02           # desert clutter intensity relative to a unit course scatterer
PYRAMID = dict(base=215.0, height=136.0, course=1.0)


def geometry():
    return Geometry()


def cache_path(name, **params):
    # Include implementation and shared settings, so edits cannot reuse old results.
    source = hashlib.sha256()
    for path in sorted(Path(ROOT, "sarsim").glob("*.py")) + [Path(__file__)]:
        source.update(path.read_bytes())
    params = {"source": source.hexdigest(), "params": params}
    key = hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()[:16]
    return os.path.join(CACHE, f"{name}_{key}.npy")


def array_fingerprint(array):
    array = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str((array.dtype.str, array.shape)).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def pyramid_scene(geom, rotation_deg=8.0, seed=0, **kw):
    ext = (-160, 160, -170, 170)
    p = stepped_pyramid(geom, rotation_deg=rotation_deg, rng=seed, **PYRAMID, **kw)
    rocks = desert_rocks(geom, ext, rng=seed + 100)
    # drop rocks under the pyramid footprint
    from sarsim.scene import PyramidShape
    inside = PyramidShape(PYRAMID["base"], PYRAMID["height"], rotation_deg)(rocks.x, rocks.y) > 0
    return concat(p, rocks.subset(~inside))


def desert_scene(geom, seed=0):
    return desert_rocks(geom, (-160, 160, -170, 170), rng=seed + 100)


def make_slc(name, scat, geom, seed=0, ground=GROUND_I, shape=SHAPE, **params):
    scatterers = None if scat is None else {k: array_fingerprint(v) for k, v in vars(scat).items()}
    path = cache_path(name, seed=seed, ground=ground, shape=shape, geom=asdict(geom),
                      scatterers=scatterers, r_centre=R_CENTRE, **params)
    if os.path.exists(path):
        return np.load(path)
    t0 = time.time()
    slc = synthesize(scat, geom, shape, r_centre=R_CENTRE, ground_intensity=ground, seed=seed, verbose=False)
    np.save(path, slc)
    print(f"  synthesized {name}: {scat.n if scat is not None else 0} scatterers in {time.time() - t0:.0f}s", flush=True)
    return slc


def axes(geom, shape=SHAPE):
    return pixel_axes(geom, shape, 0.0, R_CENTRE)


def write_summary(test_id, summary):
    from experiments.reporting import refine_summary
    summary = refine_summary(summary)
    d = os.path.join(RESULTS, test_id)
    os.makedirs(os.path.join(d, "figs"), exist_ok=True)
    for figure in summary.get("figures", []):
        path = Path(figure["file"])
        if path.is_absolute():
            figure["file"] = str(path.relative_to(ROOT))
    summary["provenance"] = {
        "python": sys.version.split()[0], "numpy": np.__version__,
        "source_sha256": hashlib.sha256(Path(sys.argv[0]).read_bytes()).hexdigest()
            if Path(sys.argv[0]).is_file() else None,
    }
    with open(os.path.join(d, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float, allow_nan=False)
    print(f"wrote {d}/summary.json")


def figdir(test_id):
    d = os.path.join(RESULTS, test_id, "figs")
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------- scene analysis suite
from sarsim import run_pipeline, line_targets, grid_targets  # noqa: E402

GRID_X = (-135.0, 135.0)     # azimuth extent of the 2-D target grid (m)
GRID_R = (-125.0, 75.0)      # slant-range extent (m)
GRID_STRIDE = 4              # pixels (1 m)
TRANSECT_AZ_R = -20.0        # azimuth transect at this slant range (crosses both side faces)
TRANSECT_RG_X = 0.0          # range transect at this azimuth (through the near face)


def default_bank():
    return SubapBank(K=50, sub_frac=0.5, delta_hz=88.0, snap_bins=True, taper_bins=0)


def transects(geom, shape=SHAPE):
    x, r = axes(geom, shape)
    c = int(np.argmin(np.abs(r - TRANSECT_AZ_R)))
    rows_a, cols_a = line_targets(int(np.searchsorted(x, -150)), c, int(np.searchsorted(x, 150)), c)
    rr = int(np.argmin(np.abs(x - TRANSECT_RG_X)))
    rows_r, cols_r = line_targets(rr, int(np.searchsorted(r, -130)), rr, int(np.searchsorted(r, 90)))
    return {"az": (rows_a, cols_a), "rg": (rows_r, cols_r)}


def grid(geom, shape=SHAPE, stride=GRID_STRIDE):
    x, r = axes(geom, shape)
    return grid_targets(int(np.searchsorted(x, GRID_X[0])), int(np.searchsorted(x, GRID_X[1])),
                        int(np.searchsorted(r, GRID_R[0])), int(np.searchsorted(r, GRID_R[1])), stride)


def analyze_scene(tag, slc, geom, bank=None, patch=32, lam_s=0.48, force=False, **kw):
    """Run transects + the 2-D grid through the pipeline; cache the arrays."""
    bank = default_bank() if bank is None else bank
    path = cache_path("analysis_" + tag, bank=bank.__dict__, patch=patch, lam_s=lam_s,
                      slc=array_fingerprint(slc), geom=asdict(geom), **kw)
    path = path.replace(".npy", ".npz")
    if os.path.exists(path) and not force:
        d = np.load(path, allow_pickle=False)
        return {k: d[k] for k in d.files}
    t0 = time.time()
    out = {}
    for name, (rows, cols) in transects(geom).items():
        o = run_pipeline(slc, geom, rows, cols, bank=bank, patch=patch, lam_s=lam_s, verbose=False, **kw)
        for k in ("q", "tomo_paper", "score_w", "best_w", "rows", "cols"):
            out[f"{name}_{k}"] = o[k]
        out["z"] = o["z"]
        out["kz"] = o["kz"]
    rows, cols, (rs, cs) = grid(geom)
    o = run_pipeline(slc, geom, rows, cols, bank=bank, patch=patch, lam_s=lam_s, verbose=True, **kw)
    out["grid_q"] = o["q"]
    out["grid_score"] = o["score_w"].reshape(rs.size, cs.size, -1).astype(np.float32)
    out["grid_paper"] = o["tomo_paper"].reshape(rs.size, cs.size, -1).astype(np.float32)
    out["grid_rows"], out["grid_cols"] = rs, cs
    out["runtime_s"] = time.time() - t0
    np.savez(path, **out)
    print(f"  analysis {tag}: {time.time() - t0:.0f}s", flush=True)
    return out


def region_masks(geom, rs, cs, rotation_deg=8.0, shape=SHAPE):
    """Classify grid pixels by what their patch mostly contains: near face / side faces / desert / apex layover."""
    x, r = axes(geom, shape)
    X, R = np.meshgrid(x[rs], r[cs], indexing="ij")
    from sarsim.scene import PyramidShape
    # ground footprint (unrotated square, rotated by rotation_deg): point on the ground plane
    # a pixel at slant range R with the ground at z=0 lies at y = R / sin(theta)
    y_ground = R / np.sin(geom.theta)
    shape_fn = PyramidShape(PYRAMID["base"], PYRAMID["height"], rotation_deg)
    inside = shape_fn(X, y_ground) > 0
    return {"inside_footprint": inside, "x": X, "r": R}
