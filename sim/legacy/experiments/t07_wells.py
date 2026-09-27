"""Illustrative appendix (formerly T7): hand-placed surface points and chosen rendering settings.

This file retains its original name for existing links. It does not explain the
Khafre observations. The controlled replacement is t07_surface_controls.py.
"""
import os, time, json, numpy as np
from common import *
from sarsim import viz, run_pipeline, point_targets, concat
from sarsim.scene import rot2d
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

TEST = "t07_wells"
ROT = 8.0
STRIDE = 5                     # pixels = 1.25 m
LAM_S_SLIDES = 2.2             # chosen to illustrate a ~630 m display extent; not an inferred calibration
BRIGHT_AMP = 8.0


def build_scene(g):
    L2 = PYRAMID["base"] / 2
    cx = np.array([-L2, L2, L2, -L2, 0, L2, 0, -L2]); cy = np.array([-L2, -L2, L2, L2, -L2, 0, L2, 0])
    xr, yr = rot2d(cx, cy, ROT)
    bright = point_targets(xr, yr, np.zeros(8), np.full(8, BRIGHT_AMP), label=50)
    return concat(pyramid_scene(g, rotation_deg=ROT), bright), xr, yr


def volume(g, bank, slc):
    path = cache_path("t07_volume", stride=STRIDE, bank=bank.__dict__,
                      slc=array_fingerprint(slc), geom=g.__dict__).replace(".npy", ".npz")
    x, r = axes(g)
    rows = np.arange(np.searchsorted(x, -135), np.searchsorted(x, 135), STRIDE)
    cols = np.arange(np.searchsorted(r, -125), np.searchsorted(r, 85), STRIDE)
    if os.path.exists(path):
        d = np.load(path)
        return d["H"], d["z"], d["kz"], rows, cols
    RR, CC = np.meshgrid(rows, cols, indexing="ij")
    o = run_pipeline(slc, g, RR.ravel(), CC.ravel(), bank=bank, patch=32, common="global", verbose=True, modes=("paper",))
    z, kz, q = o["z"], o["kz"], o["q"]
    A = np.exp(-1j * np.outer(kz, z))
    H = ((q[..., 0] + 1j * q[..., 1]) @ A).reshape(rows.size, cols.size, -1).astype(np.complex64)
    np.savez(path, H=H, z=z, kz=kz)
    return H, z, kz, rows, cols


def main():
    viz.style()
    t0 = time.time()
    g = geometry(); bank = default_bank(); fd = figdir(TEST); figs = []
    scene, bx, by = build_scene(g)
    slc = make_slc("pyramid_rot8_corners", scene, g, rotation_deg=ROT, corners=8, amp=BRIGHT_AMP)
    x, r = axes(g)
    H, z, kz, rows, cols = volume(g, bank, slc)
    xs, rs = x[rows], r[cols]
    br = by * np.sin(g.theta)                                   # slant range of the bright scatterers
    P = np.log10(np.abs(H) ** 2 + 1e-14)
    V = gaussian_filter(P, sigma=(1.0, 1.0, 1.5))
    zs = z * LAM_S_SLIDES / 0.48                                # deliberately rescaled illustration axis
    ext_map = [xs[0], xs[-1], rs[-1], rs[0]]
    # arris geometry: the near face's two edges run from the near base corners to the laid-over apex
    L2 = PYRAMID["base"] / 2
    sq = np.array([[-L2, -L2], [L2, -L2], [L2, L2], [-L2, L2], [-L2, -L2]])
    qx, qy = rot2d(sq[:, 0], sq[:, 1], ROT)
    apex = (0.0, -PYRAMID["height"] * np.cos(g.theta))
    R_SLICE = -95.0
    tot_s = gaussian_filter(np.log10((np.abs(H) ** 2).sum(axis=2)), 1.5)
    c_sl = int(np.argmin(np.abs(rs - R_SLICE)))
    prof = tot_s[:, c_sl]
    band_x = [xs[np.argmax(np.where(xs < 0, prof, -np.inf))], xs[np.argmax(np.where(xs > 0, prof, -np.inf))]]   # the two bright bands at this range
    # ---------- figure 1: plan view, total power
    tot = np.log10((np.abs(H) ** 2).sum(axis=2))
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4))
    sub = slc[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
    viz.imshow_slc(axs[0], sub, x[rows[0]:rows[-1] + 1], r[cols[0]:cols[-1] + 1])
    axs[0].scatter(bx, br, s=60, facecolor="none", edgecolor=viz.SERIES[1], lw=1.4)
    axs[0].set_title("SLC: stepped pyramid with eight bright corner-like scatterers (circled)", loc="left")
    axs[1].imshow(gaussian_filter(tot, 0.8).T, origin="upper", extent=ext_map, cmap="jet", aspect="equal", interpolation="bicubic")
    axs[1].plot(qx, qy * np.sin(g.theta), color="white", lw=1.2, alpha=0.8)
    axs[1].axhline(R_SLICE, color="white", lw=0.8, ls="--", alpha=0.8)
    axs[1].set_title("plan view of the volume: total trajectory power, jet colormap", loc="left"); axs[1].set_xlabel("azimuth (m)"); axs[1].set_yticklabels([])
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t07_plan.png"), "caption":
                 "<b>Constructed input and rendered power.</b> Eight ideal point scatterers are deliberately added at corners and edge midpoints. Their realism at Khafre is not established. The output includes dark regions near those points and bright face-edge bands. The dashed line marks the selected section."})
    # ---------- figure 2: vertical slices, published style
    fig, axs = plt.subplots(3, 1, figsize=(12.5, 12.5))
    def show(ax, sl, extent, title):
        ax.imshow(sl.T, origin="upper", aspect="auto", extent=extent, cmap="jet", interpolation="bicubic",
                  vmin=np.percentile(V, 3), vmax=np.percentile(V, 99.5))
        ax.set_title(title, loc="left"); ax.set_ylabel("illustrative rescaled depth (m)")
        a2 = ax.twinx(); a2.set_ylim(z[-1], z[0]); a2.set_ylabel("model depth, λs = 0.48 m", color=viz.MUTED); a2.tick_params(colors=viz.MUTED)
    c0 = int(np.argmin(np.abs(rs - R_SLICE)))
    show(axs[0], V[:, c0, :], [xs[0], xs[-1], zs[-1], zs[0]], f"azimuth slice at slant range {rs[c0]:.0f} m, crossing the two bright bands along the near face's upper edges (dotted)")
    for xb in band_x:
        axs[0].axvline(xb, color="white", lw=0.8, ls=":", alpha=0.9)
    j0 = int(np.argmin(np.abs(xs - band_x[0])))
    show(axs[1], V[j0, :, :], [rs[0], rs[-1], zs[-1], zs[0]], f"range slice at azimuth {xs[j0]:.0f} m, through the left band (bright column at {R_SLICE:.0f} m)")
    axs[1].axvline(R_SLICE, color="white", lw=0.8, ls=":", alpha=0.9)
    i_pt = int(np.argmax(bx))                                        # a mid-edge / corner point away from the arrises
    jp = int(np.argmin(np.abs(xs - bx[i_pt])))
    show(axs[2], V[jp, :, :], [rs[0], rs[-1], zs[-1], zs[0]], f"range slice at azimuth {xs[jp]:.0f} m, through a bright point at {br[i_pt]:.0f} m (dark column)")
    axs[2].axvline(br[i_pt], color="white", lw=0.8, ls=":", alpha=0.9); axs[2].set_xlabel("slant range (m)")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t07_slices.png"), "caption":
                 "<b>Sections through an illustrative volume.</b> These views use smoothing and a chosen wavelength of 2.2 m to stretch the displayed depth extent. The right axes retain the 0.48 m model scale. This is not a measured depth calibration for the published images."})
    # ---------- figure 2b: side by side with the published slide
    PUB = os.path.join(RESULTS, TEST, "published")
    import matplotlib.image as mpimg
    pub = mpimg.imread(os.path.join(PUB, "khafre_columns_slide.jpg"))
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [pub.shape[1] / pub.shape[0], 1.75]})
    axs[0].imshow(pub); axs[0].set_axis_off()
    axs[0].set_title("published: Khafre Project tomogram slide, March 2025", loc="left")
    j1 = int(np.argmin(np.abs(xs - band_x[1])))
    sl = V[j1, :, :]
    axs[1].imshow(sl.T, origin="upper", aspect="auto", extent=[rs[0], rs[-1], zs[-1], zs[0]], cmap="jet", interpolation="bicubic",
                  vmin=np.percentile(V, 3), vmax=np.percentile(V, 99.5))
    axs[1].set_title(f"simulation: range slice at azimuth {xs[j1]:.0f} m through the right band, static pyramid", loc="left")
    axs[1].set_xlabel("slant range (m)"); axs[1].set_ylabel("illustrative rescaled depth (m)")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t07_compare.png"), "caption":
                 "<b>A visual comparison only.</b> Left: a Khafre Research Project slide, copyright its authors. Right: a constructed surface-only simulation with chosen display settings. Resemblance does not establish the origin of the published features or reproduce their processing."})
    # ---------- published images: what was shown, and what was drawn from it
    figs.append({"file": os.path.join(PUB, "khafre_composite_slide.jpg"), "caption":
                 "<b>Published composite slide.</b> Plan view of the pyramid footprint with circled 'wells', vertical sections with columns, and hand annotations placing Khafre's known chambers (Belzoni's chamber, the lower chamber) on blobs. Reproduced at reduced size for critique; copyright the Khafre Research Project."})
    figs.append({"file": os.path.join(PUB, "khafre_cad_ramps.jpg"), "caption":
                 "<b>Published CAD illustration.</b> Cylinders and spiral ramps depict the project's interpretation. Their geometry is not independently validated here. Reproduced at reduced size for critique; copyright the Khafre Research Project."})
    figs.append({"file": os.path.join(PUB, "khafre_cad_coils.jpg"), "caption":
                 "<b>Published CAD illustration.</b> Cylinders, coils and basal blocks depict the project's interpretation. Their depths and shapes are not measured by this illustration. Reproduced at reduced size for critique; copyright the Khafre Research Project."})
    # ---------- figure 3: 3-D isosurfaces, bright walls and dark wells
    from skimage import measure
    dx_ = xs[1] - xs[0]; dr_ = rs[1] - rs[0]; dz_ = zs[1] - zs[0]
    def mesh(level):
        verts, faces, _, _ = measure.marching_cubes(V, level=level, step_size=2, allow_degenerate=False)
        pts = np.column_stack([xs[0] + verts[:, 0] * dx_, rs[0] + verts[:, 1] * dr_, zs[0] + verts[:, 2] * dz_])
        return pts[faces]
    light = np.array([-0.5, 0.4, -0.75]); light /= np.linalg.norm(light)
    def shaded(tri, rgb):
        fn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12
        sh = 0.55 + 0.45 * np.clip(np.abs(fn @ light), 0, 1)
        return np.clip(np.array(rgb)[None, :] * sh[:, None], 0, 1)
    Vz = gaussian_filter(P, sigma=(1.0, 1.0, 8.0))                 # smoothed along depth for the 3-D view
    hi = np.percentile(Vz, 98.0); lo = np.percentile(Vz, 2.0)
    def mesh(level):
        verts, faces, _, _ = measure.marching_cubes(Vz, level=level, step_size=2, allow_degenerate=False)
        pts = np.column_stack([xs[0] + verts[:, 0] * dx_, rs[0] + verts[:, 1] * dr_, zs[0] + verts[:, 2] * dz_])
        return pts[faces]
    tri_hi = mesh(hi); tri_lo = mesh(lo)
    fig = plt.figure(figsize=(12, 9.5))
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Poly3DCollection(tri_lo, facecolors=shaded(tri_lo, (0.16, 0.36, 0.62)), edgecolors="none"))
    ax.add_collection3d(Poly3DCollection(tri_hi, facecolors=shaded(tri_hi, (0.92, 0.41, 0.20)), edgecolors="none", alpha=0.85))
    ax.plot(qx, qy * np.sin(g.theta), 0, color=viz.INK2, lw=1.2)
    ax.scatter(bx, br, np.zeros(8), s=28, color=viz.SERIES[1], depthshade=False)
    ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(rs[0], rs[-1]); ax.set_zlim(zs[-1], 0)
    ax.set_box_aspect((xs[-1] - xs[0], rs[-1] - rs[0], zs[-1] * 0.55))
    ax.view_init(elev=24, azim=-58)
    ax.set_xlabel("azimuth (m)"); ax.set_ylabel("slant range (m)"); ax.set_zlabel("depth as labelled (m)")
    ax.set_title("isosurfaces of the volume smoothed along depth: lowest 2% (blue) and highest 2% (orange) of voxels", loc="left")
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t07_isosurface.png"), "caption":
                 "<b>Chosen isosurfaces after depth smoothing.</b> Blue marks the lowest-power regions and orange the highest-power regions. The thresholds and smoothing determine the rendered shapes; all share the finite computed volume. No spiral geometry or published feature is reproduced here."})
    # ---------- figure 4: rungs
    spacing = 0.79 / np.cos(np.deg2rad(ROT))
    scale = g.depth_scale(0.48)
    from t03_mechanism import comb_period
    def region_spec(mask):
        s = (np.abs(H[mask]) ** 2); s = s / s.sum(axis=1, keepdims=True); return s.mean(axis=0)
    X, Rg = np.meshgrid(xs, rs, indexing="ij")
    arris_mask = (tot_s > np.percentile(tot_s, 97)) & (Rg < -60)          # the bright bands
    point_mask = np.zeros(X.shape, bool)
    for i in range(8):
        point_mask |= (np.abs(X - bx[i]) < 3.5) & (np.abs(Rg - br[i]) < 3.5)
    desert_mask = (Rg > 60)
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    rung_rows = []
    for name, mask, col in [("bright bands (walls)", arris_mask, viz.SERIES[1]), ("wells under bright points (dark)", point_mask, viz.SERIES[0]), ("open desert", desert_mask, viz.SERIES[2])]:
        sp = region_spec(mask); d, sc = comb_period(z, sp, dmin=4.0, dmax=12.0, H=10)
        rung_rows.append({"region": name, "pixels": int(mask.sum()), "comb_spacing_m": d, "comb_score": sc})
        ax.plot(z, sp, color=col, lw=1.8, label=f"{name}: comb spacing {d:.1f} m")
    for h in range(1, 14):
        ax.axvline(scale * spacing * h, color=viz.MUTED, lw=0.8, ls=":")
    ax.set_xlim(0, 100); ax.set_xlabel("model depth (m), λs = 0.48 m"); ax.set_ylabel("mean normalised |h(z)|²"); ax.legend(fontsize=8.5)
    ax.set_title(f"Rungs: mean depth spectra by region, dotted at multiples of the projected course spacing ({scale*spacing:.1f} m)", loc="left")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t07_rungs.png"), "caption":
                 f"<b>Spectra in selected regions.</b> Mean depth spectra under the face edges, added points and open desert. Dotted lines are {scale*spacing:.1f} model-metres apart. The fitted periods vary by region and do not cleanly follow one course spacing; this does not measure or explain a published spiral."})
    pt_ratio = float(np.median((np.abs(H[point_mask]) ** 2).sum(axis=1)) / np.median((np.abs(H[desert_mask]) ** 2).sum(axis=1)))
    ar_ratio = float(np.median((np.abs(H[arris_mask]) ** 2).sum(axis=1)) / np.median((np.abs(H[desert_mask]) ** 2).sum(axis=1)))
    d = rung_rows[0]["comb_spacing_m"]; score = rung_rows[0]["comb_score"]
    m = {"n_bright": 8, "bright_amplitude": BRIGHT_AMP, "grid_pixels": int(rows.size * cols.size), "depth_axis_max_m_slides": float(zs[-1]),
         "depth_axis_max_m_paper": float(z[-1]), "power_ratio_point_wells_vs_desert": pt_ratio, "power_ratio_bright_bands_vs_desert": ar_ratio,
         "rungs": rung_rows, "predicted_course_spacing_depth_m": float(scale * spacing), "isosurface_percentiles": [2.0, 98.0], "runtime_s": time.time() - t0}
    summary = {
        "id": TEST, "order": 7, "role": "illustration", "tag": "illustration", "eyebrow": "Appendix · Constructed rendering",
        "title": "Illustration: surface points and rendering choices",
        "question": "What does a deliberately constructed scene look like after chosen rendering settings?",
        "finding": f"Near the hand-placed points, power is {pt_ratio:.2f} times the desert reference; the selected bright bands have {ar_ratio:.2f} times that reference. Smoothing and thresholding make column-like forms in this example.",
        "limitations": "No field data, independently justified surface targets, full acceptance gate or validated depth calibration. Point count and the displayed depth extent were deliberately chosen. This is an illustration, not evidence that the Khafre features share this origin.",
        "method": "T1's pyramid plus eight point scatterers of amplitude 8 at the base corners and mid-edges; paper-variant volume on a 1.25 m grid (36,000 pixels); log-magnitude smoothed with a Gaussian of 1 sample in position and 1.5 in depth; jet colormap and bicubic interpolation to match the published rendering; marching-cubes isosurfaces at the 2nd and 98th percentiles of the volume smoothed along depth.",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1, default=float))


if __name__ == "__main__":
    main()
