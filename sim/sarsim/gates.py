"""Per-window selection in the public derivative protocol, version 1.5.

These are image-trajectory criteria, not probabilities of underground structure.
Keep window identity: selecting a global winner before gating is not equivalent.
"""
from dataclasses import dataclass
import numpy as np

R2_MIN, BA_MIN, B_MIN = 0.25, 0.1, 0.005


@dataclass
class WindowGates:
    adjusted_r2: np.ndarray
    axis_ratio: np.ndarray
    minor_axis: np.ndarray
    mode: np.ndarray
    shape: np.ndarray
    full: np.ndarray


def ellipse_windows(q, W=25, modes=range(1, 7), r2_min=R2_MIN,
                    ratio_min=BA_MIN, minor_min=B_MIN):
    """Choose the best harmonic *within each window*, then apply all gates.

    Arrays have shape [pixel, window]. Degenerate/constant input is rejected.
    """
    q = np.asarray(q, float)
    if q.ndim != 3 or q.shape[2] != 2 or not np.isfinite(q).all():
        raise ValueError("expected finite trajectories [pixel, look, 2]")
    if not 4 <= W <= q.shape[1]:
        raise ValueError("window length must be between 4 and the look count")
    modes = tuple(modes)
    if not modes or any(m <= 0 or m >= W / 2 for m in modes):
        raise ValueError("harmonics must be positive and below Nyquist")
    shape = (q.shape[0], q.shape[1] - W + 1)
    best = np.full(shape, -np.inf)
    ratio = np.zeros(shape)
    minor = np.zeros(shape)
    chosen = np.zeros(shape, int)
    n = np.arange(W)
    for w in range(shape[1]):
        y = q[:, w:w + W]
        total = ((y - y.mean(axis=1, keepdims=True)) ** 2).sum(axis=(1, 2))
        usable = total > np.maximum(1e-30, (y ** 2).sum(axis=(1, 2)) * 1e-14)
        for m in modes:
            X = np.column_stack((np.ones(W), np.cos(2*np.pi*m*n/W), np.sin(2*np.pi*m*n/W)))
            coeff = np.einsum('ik,pkj->pij', np.linalg.pinv(X), y)
            residual = ((y - np.einsum('ki,pij->pkj', X, coeff)) ** 2).sum(axis=(1, 2))
            adj = 1 - residual / np.maximum(total, 1e-30) * (2*W-1)/(2*W-7)
            sv = np.linalg.svd(coeff[:, 1:].transpose(0, 2, 1), compute_uv=False)
            better = usable & (adj > best[:, w])
            best[better, w] = adj[better]
            ratio[better, w] = sv[better, 1] / np.maximum(sv[better, 0], 1e-30)
            minor[better, w] = sv[better, 1]
            chosen[better, w] = m
    shape_pass = (best >= r2_min) & (ratio >= ratio_min)
    return WindowGates(best, ratio, minor, chosen, shape_pass,
                       shape_pass & (minor >= minor_min))


def ellipse_gate(q):
    """Compatibility summary; acceptance means *any* eligible window passes.

    First three outputs describe the largest-R² window for diagnostic use only.
    Its metrics must not be used to reconstruct the final acceptance decision.
    """
    gates = ellipse_windows(q)
    p = np.arange(len(q))
    w = gates.adjusted_r2.argmax(axis=1)
    return (gates.adjusted_r2[p, w], gates.axis_ratio[p, w], gates.minor_axis[p, w],
            gates.shape.any(axis=1), gates.full.any(axis=1))


def geometric_median(vectors, tolerance=1e-10, iterations=200):
    """Robust two-component centre, with a regularized Weiszfeld iteration."""
    vectors = np.asarray(vectors, float)
    if vectors.ndim != 3 or vectors.shape[2] != 2 or not len(vectors) or not np.isfinite(vectors).all():
        raise ValueError("expected finite flank vectors [sample, look, 2]")
    centre = np.median(vectors, axis=0)
    for _ in range(iterations):
        distance = np.linalg.norm(vectors-centre, axis=2)
        weight = 1 / np.maximum(distance, tolerance)
        update = (vectors*weight[..., None]).sum(axis=0)/weight.sum(axis=0)[:, None]
        if np.max(np.abs(update-centre)) < tolerance:
            return update
        centre = update
    return centre


def two_arc_correction(target, plus, minus):
    """Apply separately supplied flank samples; geometry is caller-validated.

    Return the residual and QC quantities without interpreting residual as motion.
    """
    gp, gm = geometric_median(plus), geometric_median(minus)
    if target.shape[1:] != gp.shape or gp.shape != gm.shape:
        raise ValueError("target and both arcs must share their look axes")
    common = (gp+gm)/2
    return target-common, {
        "common": common, "plus": gp, "minus": gm,
        "disagreement": np.linalg.norm(gp-gm, axis=1),
        "plus_scatter": np.median(np.linalg.norm(plus-gp, axis=2), axis=0),
        "minus_scatter": np.median(np.linalg.norm(minus-gm, axis=2), axis=0),
    }
