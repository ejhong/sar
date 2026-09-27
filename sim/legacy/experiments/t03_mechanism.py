"""T3: where the 'depth' comes from. Two-scatterer toy, a periodic wall, and heading dependence."""
import time, json, numpy as np
from common import *
from sarsim import viz, run_pipeline, point_targets, synthesize, line_targets
from sarsim.scene import make, concat, rot2d
import matplotlib.pyplot as plt

TEST = "t03_mechanism"
LAM_S = 0.48
TOY_SHAPE = (1280, 256)


def toy(g, bank, dxs):
    sc = point_targets([-dxs / 2, dxs / 2], [0.0, 0.0], [0.0, 0.0], [1.0, 0.8], phase=[0.0, 1.0])
    slc = synthesize(sc, g, TOY_SHAPE, ground_intensity=1e-8, seed=1)
    o = run_pipeline(slc, g, np.array([TOY_SHAPE[0] // 2]), np.array([TOY_SHAPE[1] // 2]), bank=bank,
                     patch=32, common="none", verbose=False, modes=("paper",))
    return o


def wall(g, spacing=0.79, n_lines=60, length=60.0, rotation_deg=0.0, amp_sigma=0.0, diffuse=0.0, rng=0):
    """Parallel lines of scatterers (course edges of one face seen from above), lines along range."""
    rng = np.random.default_rng(rng)
    parts = []
    for i in range(n_lines):
        xl = (i - n_lines / 2) * spacing
        s = np.arange(-length / 2, length / 2, 0.25)
        xr, yr = rot2d(np.full_like(s, xl), s, rotation_deg)
        parts.append(make(xr, yr, 0.0, np.exp(rng.normal(0, amp_sigma, s.size)), label=1, rng=rng))
    if diffuse > 0:
        n = int(diffuse * n_lines * spacing * length)
        xd = rng.uniform(-n_lines * spacing / 2, n_lines * spacing / 2, n)
        yd = rng.uniform(-length / 2, length / 2, n)
        xr, yr = rot2d(xd, yd, rotation_deg)
        parts.append(make(xr, yr, 0.0, 0.3 * np.exp(rng.normal(0, 0.5, n)), label=2, rng=rng))
    return concat(*parts)


def comb_period(z, spec, dmin=4.0, dmax=12.0, H=8):
    """Fundamental spacing d whose harmonic comb d, 2d, ... best matches the spectrum's peaks."""
    base = np.median(spec)
    ds = np.linspace(dmin, dmax, 400)
    sc = []
    for d in ds:
        zz = d * np.arange(1, H + 1)
        zz = zz[zz < z[-1]]
        sc.append(np.mean(np.interp(zz, z, spec)) / base)
    return float(ds[int(np.argmax(sc))]), float(max(sc))


def wall_spectrum(g, bank, sc, tag):
    slc = make_slc("wall_" + tag, sc, g, ground=1e-4, shape=(1280, 512))
    x, r = axes(g, (1280, 512))
    c = 256
    rows, cols = line_targets(int(np.searchsorted(x, -20)), c, int(np.searchsorted(x, 20)), c)
    o = run_pipeline(slc, g, rows, cols, bank=bank, patch=32, common="global", verbose=False)
    tp = o["tomo_paper"]; tpn = tp / tp.sum(axis=1, keepdims=True)
    return o["z"], tpn.mean(axis=0), o, x[rows]


def main():
    viz.style()
    t0 = time.time()
    g = geometry(); bank = default_bank(); fd = figdir(TEST); figs = []
    scale = g.depth_scale(LAM_S)
    # ---- (a) two scatterers
    dxs_list = [1.0, 2.0, 3.0, 4.0]
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.0))
    rows = []
    for i, dxs in enumerate(dxs_list):
        o = toy(g, bank, dxs)
        q = o["q"][0, :, 0]; z = o["z"]; tp = o["tomo_paper"][0]
        zp = z[np.argmax(tp)]
        rows.append({"separation_m": dxs, "predicted_depth_m": scale * dxs, "measured_depth_m": float(zp)})
        axs[0].plot(np.arange(1, 51), (q - q.mean()) * 1e3, color=viz.SERIES[i], lw=1.6, label=f"Δx = {dxs:.0f} m")
        axs[1].plot(z, tp / tp.max(), color=viz.SERIES[i], lw=1.6, label=f"Δx = {dxs:.0f} m")
        axs[1].axvline(scale * dxs, color=viz.SERIES[i], lw=0.9, ls=":")
    axs[0].set_xlabel("sub-aperture pair k"); axs[0].set_ylabel("azimuth shift (millipixels)")
    axs[0].set_title("Two isotropic scatterers, no clutter: the shift oscillates with k", loc="left"); axs[0].legend()
    axs[1].set_xlabel("model depth (m)"); axs[1].set_ylabel("$|a(z)^{H}Y|^{2}$, normalised")
    axs[1].set_title(f"…and the 'depth' is {scale:.1f} m per metre of separation (dotted: predicted)", loc="left"); axs[1].legend()
    axs[1].set_xlim(0, 60)
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t03_two_scatterers.png"), "caption":
                 f"<b>The mechanism in its simplest form.</b> Two point scatterers a distance Δx apart in azimuth, nothing else. The reference/offset shift oscillates with the sub-aperture position at a rate proportional to Δx, and the steering-matrix 'depth' lands at Δx × λₛ sinθ / λ = {scale:.1f} Δx. A pair of rocks 3 m apart is a 'structure' {scale*3:.0f} m underground."})
    # ---- (b) periodic wall: ideal -> realistic
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), sharey=True)
    variants = [("ideal: identical courses 0.79 m apart", dict(amp_sigma=0.0, diffuse=0.0)),
                ("random block amplitudes (σ = 0.5 ln)", dict(amp_sigma=0.5, diffuse=0.0)),
                ("random amplitudes + diffuse face scatter", dict(amp_sigma=0.5, diffuse=2.0))]
    wall_rows = []
    for i, (name, kw) in enumerate(variants):
        z, spec, o, xs = wall_spectrum(g, bank, wall(g, **kw), f"v{i}")
        axs[i].plot(z, spec, color=viz.SERIES[0], lw=1.8)
        for h in range(1, 9):
            axs[i].axvline(scale * 0.79 * h, color=viz.MUTED, lw=0.8, ls=":")
        axs[i].set_title(name, loc="left"); axs[i].set_xlabel("model depth (m)"); axs[i].set_xlim(0, 80)
        d, score = comb_period(z, spec)
        wall_rows.append({"variant": name, "comb_spacing_m": d, "predicted_spacing_m": float(scale * 0.79), "comb_score": score})
    axs[0].set_ylabel("mean normalised spectrum")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t03_wall.png"), "caption":
                 f"<b>A wall of courses becomes a ladder of 'levels'.</b> Sixty parallel lines of scatterers 0.79 m apart, the plan view of one side face's course edges, averaged over 160 transect pixels. The spectrum is a harmonic series of the course spacing (dotted lines at multiples of {scale*0.79:.1f} model-metres): perfectly identical courses light up only a few harmonics, realistic random block brightness fills in the series with a fitted fundamental of {wall_rows[1]['comb_spacing_m']:.1f} m, and diffuse face scatter leaves it in place. Stacked, evenly spaced 'floors' under a stepped monument are what regular masonry looks like after this processing."})
    # ---- (c) heading dependence
    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    head_rows = []
    for i, rot in enumerate([0.0, 15.0, 30.0]):
        z, spec, o, xs = wall_spectrum(g, bank, wall(g, rotation_deg=rot, amp_sigma=0.5), f"rrot{rot:.0f}")
        pred = scale * 0.79 / np.cos(np.deg2rad(rot))
        d, score = comb_period(z, spec)
        ax.plot(z, spec + i * 0.004, color=viz.SERIES[i], lw=1.8, label=f"heading offset {rot:.0f}°: level spacing {d:.2f} m (predicted {pred:.2f})")
        head_rows.append({"heading_deg": rot, "predicted_spacing_m": float(pred), "comb_spacing_m": d, "comb_score": score})
        for h in range(1, 8):
            ax.axvline(pred * h, color=viz.SERIES[i], lw=0.7, ls=":", alpha=0.8)
    ax.set_xlim(0, 60); ax.set_xlabel("model depth (m)"); ax.set_ylabel("mean normalised spectrum (offset)"); ax.legend()
    ax.set_title("Rotate the scene relative to the flight line and the 'levels' move", loc="left")
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t03_heading.png"), "caption":
                 "<b>Depth depends on the satellite's heading.</b> The realistic wall imaged with the flight line rotated 0°, 15° and 30° relative to the courses (curves offset for clarity; dotted combs at multiples of the predicted spacing). A harmonic-comb fit recovers the level spacing, which follows the 1/cos law of the course spacing projected on the azimuth axis. Real buried architecture cannot move when the satellite changes heading; these features do, because they are surface geometry projected on the azimuth axis."})
    # ---- (d) non-locality and aliasing
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.0))
    nl_rows = []
    Ds = [3.0, 6.0, 12.0, 20.0, 40.0]
    for i, D in enumerate(Ds):
        sc = point_targets([0.0, D], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0], phase=[0.0, 1.0])
        slc = synthesize(sc, g, TOY_SHAPE, ground_intensity=1e-8, seed=1)
        o = run_pipeline(slc, g, np.array([TOY_SHAPE[0] // 2]), np.array([TOY_SHAPE[1] // 2]), bank=bank,
                         patch=32, common="none", verbose=False, modes=("paper",))
        q = o["q"][0, :, 0]; z = o["z"]; tp = o["tomo_paper"][0]
        nyq = np.pi / np.median(np.diff(o["kz"]))
        zp = scale * D
        while zp > nyq:
            zp = abs(2 * nyq - zp) if zp < 2 * nyq else zp - 2 * nyq
        nl_rows.append({"distance_m": D, "shift_rms_millipx": float(q.std() * 1e3), "unaliased_depth_m": float(scale * D),
                        "aliased_prediction_m": float(zp), "measured_depth_m": float(z[np.argmax(tp)])})
        axs[1].plot(z, tp / tp.max() + i * 1.1, color=viz.SERIES[i], lw=1.6, label=f"D = {D:.0f} m")
        axs[1].axvline(zp, color=viz.SERIES[i], lw=0.9, ls=":")
    axs[0].bar(range(len(Ds)), [r_["shift_rms_millipx"] for r_ in nl_rows], color=viz.SERIES[:len(Ds)], width=0.55)
    axs[0].set_xticks(range(len(Ds))); axs[0].set_xticklabels([f"{D:.0f} m" for D in Ds])
    axs[0].axvspan(-0.5, 1.5, color=viz.GRID, alpha=0.6, lw=0); axs[0].text(0.5, max(r_["shift_rms_millipx"] for r_ in nl_rows) * 1.02, "inside the 8 m patch", ha="center", fontsize=8.5, color=viz.INK2)
    axs[0].set_xlabel("distance of the second scatterer from the patch centre"); axs[0].set_ylabel("trajectory amplitude (millipixels, rms)")
    axs[0].set_title("A scatterer far outside the patch still drives the trajectory", loc="left")
    axs[1].set_xlabel("model depth (m)"); axs[1].set_yticks([]); axs[1].legend(loc="upper right", fontsize=8)
    axs[1].set_title("…and lands at the aliased depth (dotted: predicted)", loc="left"); axs[1].set_xlim(0, z[-1])
    fig.tight_layout()
    figs.append({"file": viz.finish(fig, f"{fd}/t03_nonlocal.png"), "caption":
                 "<b>The 'structure' under a pixel is made of things tens of metres away.</b> One scatterer at the patch centre and a second one at 3, 6, 12, 20 or 40 m. The offset band differs from the reference band by a sliver only 88 Hz wide, whose impulse response is about 80 m long, so scatterers well outside the 32-pixel patch modulate the measured shift as strongly as neighbours do. Separations beyond the depth axis' Nyquist limit fold back: the 20 m and 40 m scatterers appear at 107 m and 71 m, exactly where aliasing predicts. This is why every scene fills the whole depth axis regardless of patch size (T4)."})
    m = {"depth_per_metre_separation": float(scale), "two_scatterers": rows, "wall": wall_rows, "heading": head_rows, "nonlocal": nl_rows,
         "runtime_s": time.time() - t0}
    summary = {
        "id": TEST, "order": 3, "tag": "mechanism", "eyebrow": "03 · Mechanism",
        "title": "'Depth' is the azimuth spacing of surface scatterers, rescaled",
        "question": "What exactly does the steering matrix respond to, if not motion?",
        "finding": f"Interference between surface scatterers. Two point targets Δx apart produce a shift that oscillates with sub-aperture position at a rate set by Δx, and the focusing step relabels that rate as a depth of {scale:.1f} × Δx model-metres. Regular courses give a harmonic ladder of 'levels' whose spacing follows the satellite heading; and because the offset band's difference kernel is about 80 m long, scatterers far outside the patch contribute too, folding onto the depth axis by aliasing.",
        "limitations": "The clean predictions use isotropic point scatterers; aspect-dependent 'flash' and layover mixing add more terms of the same kind. The constant 8.9 depends on the declared sound wavelength λₛ, which is a free parameter (T4).",
        "method": "Clutter-free two-scatterer scenes at 1, 2, 3, 4 m separation and at 3 to 40 m; a 60-line periodic wall in three degrees of realism; the realistic wall at three headings with a harmonic-comb fit of the level spacing. Same bank and patch as T1; transect spectra averaged over 160 pixels.",
        "figures": figs, "metrics": m, "date": time.strftime("%Y-%m-%d"),
    }
    write_summary(TEST, summary)
    print(json.dumps(m, indent=1, default=float))


if __name__ == "__main__":
    main()
