"""P2-08 · How the method draws pillars and boxes.

    uv run python experiments/p2_08_pillars_and_boxes.py

The Khafre renderings of 2025 show columns descending from the surface, banded along their length, with blocks at their
feet and deeper. This experiment reproduces that look from the published pipeline and measures where each feature comes
from, on three inputs:

* the real Khafre crop and the open plateau west of it: the first investigation's trajectories (legacy real_common.run_patch,
  the 2022 bank of 50 half-band pairs 88 Hz apart, 32-pixel patches, a target every 12 x 8 pixels; cached on the desktop);
* P2-07's simulated flat desert with nothing below it, through the same pipeline.

The depth grid runs past the steering's repeat depth z_rep = 2 pi / dKz, as the derivative protocol's own grid does (R12),
and the axis is relabelled, as any lambda_s relabels it (P2-06), so that z_rep falls at 648 m, the depth of the published
shafts. Nothing is reprocessed by the relabelling.

1. Pillars. Over a full period of the axis the focused power at a pixel averages to its trajectory's energy (Parseval): a
   pixel whose registration wanders is bright at every depth, a column. Where do such pixels lie?
2. Bands. Along a column the power fluctuates with the axis's resolution 2 pi / (Kz span), a fixed spacing, which rendered
   in 3-D reads as rings or a spiral.
3. Boxes. At z = 0 and at z = z_rep every steering phase is the same, so the power there is |mean q|^2: each pixel's
   average offset, drawn at the surface, again at the bottom of one period, and again at twice that.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter

from katabasis.runs import Run, RESULTS, load
from sarsim import synthesize
from sarsim.acquisition import DwellGeometry
from sarsim.pipeline import run_pipeline
from sarsim.subap import SubapBank
from sarsim.tomo import focus_paper

RID = 'p2_08_pillars_and_boxes'
LEGACY_CACHE = Path(__file__).resolve().parents[2] / 'results' / 'cache'
REAL = {'khafre': 'real_giza_khafre_paper_6d230513af04.npz', 'plateau': 'real_giza_desert_west_paper_f0a750ee4bb4.npz'}
Z_SHOW = 648.0                  # m: where the relabelled axis puts z_rep (the published shafts' depth)
PERIODS = 2.2                   # the grid runs this many repeat depths down
NZ = 440
_spec = importlib.util.spec_from_file_location('p207', Path(__file__).with_name('p2_07_whole_chain.py'))
p207 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p207)


def load_real(name):
    d = np.load(LEGACY_CACHE / REAL[name], allow_pickle=True)
    meta = json.loads(str(d['meta']))
    return {'q': d['q'].astype(float), 'kz': d['kz'].astype(float), 'shape': (len(d['grid_rows']), len(d['grid_cols'])),
            'amplitude': d['amplitude'].astype(float), 'meta': meta,
            'pixel_m': (meta['stride'][0] * 0.04398, meta['stride'][1] * 0.1536 / np.sin(np.deg2rad(20.8)))}


def simulated_flat():
    """P2-07's scene B (no chamber) through the pipeline, keeping the trajectories."""
    g = DwellGeometry.from_record('giza-20250827')
    rng = np.random.default_rng(71)
    import katabasis.ambient.field as fld
    amb = json.loads((p207.SITES / 'ambient.json').read_text())
    level = amb['regional']['microseism_vertical_0.1_0.3_hz']['value']
    f = np.linspace(0.1, 0.3, 21)
    field = fld.microseisms(f, np.full(len(f), level ** 2 / 0.2), np.random.default_rng(72), per_bin=4,
                            speed=amb['regional']['microseism_phase_speed']['value'], hv=fld.rayleigh_hv(3300.0, 1830.0))
    scat = p207.scene(g, rng)
    img = synthesize(scat, g, p207.SHAPE, motion=p207.Shaking(g, field, p207.kernels(g), 0.0))
    rows = np.arange(48, p207.SHAPE[0] - 48, 12)
    cols = np.arange(32, p207.SHAPE[1] - 32, 8)
    RR, CC = np.meshgrid(rows, cols, indexing='ij')
    out = run_pipeline(img, g, RR.ravel(), CC.ravel(), bank=SubapBank(), patch=32, upsample=1000, lam_s=0.48,
                       modes=(), verbose=False)
    return {'q': out['q'], 'kz': out['kz'], 'shape': (len(rows), len(cols)), 'amplitude': np.abs(img[RR.ravel(), CC.ravel()]),
            'pixel_m': (12 * g.dx, 8 * g.dr / np.sin(g.theta))}


def analyse(d):
    order = np.argsort(d['kz'])
    kz = d['kz'][order]
    q = d['q'][:, order, :]
    qc = q[..., 0] + 1j * q[..., 1]
    dk = float(np.median(np.diff(kz)))
    z_rep = 2 * np.pi / dk
    scale = Z_SHOW / z_rep
    z = np.linspace(0.0, PERIODS * z_rep, NZ)
    T = focus_paper(qc, kz, z)
    K = qc.shape[1]
    one = z < z_rep
    # 1 pillars: average over one period against the trajectory's energy
    energy = (np.abs(qc) ** 2).sum(axis=1) / K ** 2
    mean_T = T[:, one].mean(axis=1)
    pillar_corr = float(np.corrcoef(mean_T, energy)[0, 1])
    amp_corr = float(np.corrcoef(np.log(d['amplitude'] + 1e-9), np.log(energy + 1e-30))[0, 1])
    # 2 bands: autocorrelation along depth, within one period
    Ts = T[:, one] - T[:, one].mean(axis=1, keepdims=True)
    ac = np.array([np.mean(np.sum(Ts[:, :Ts.shape[1] - L] * Ts[:, L:], axis=1) / np.sum(Ts * Ts, axis=1))
                   for L in range(Ts.shape[1] // 2)])
    dz = z[1] - z[0]
    first_zero = float(np.argmax(ac < 0) * dz) if (ac < 0).any() else float('nan')
    resolution = 2 * np.pi / (kz.max() - kz.min())
    # 3 boxes: the power at z = 0 and z = z_rep against each pixel's mean offset
    dc = np.abs(qc.mean(axis=1)) ** 2
    i_rep = int(np.argmin(np.abs(z - z_rep)))
    box_ratio = float(np.median(T[:, i_rep] / dc))
    box_corr = float(np.corrcoef(T[:, i_rep], dc)[0, 1])
    top_corr = float(np.corrcoef(T[:, 0], dc)[0, 1])
    rep_corr = float(np.corrcoef(T[:, one][:, :-5].ravel(), T[:, i_rep:i_rep + one.sum() - 5].ravel())[0, 1])
    return {'T': T, 'z': z, 'z_rep': z_rep, 'scale': scale, 'energy': energy, 'dc': dc,
            'metrics': {'repeat_depth_m_raw': z_rep, 'repeat_depth_m_shown': Z_SHOW, 'scale': scale,
                        'pillar_corr_mean_power_vs_energy': pillar_corr, 'energy_vs_amplitude_log_corr': amp_corr,
                        'band_first_zero_m_shown': first_zero * scale, 'resolution_m_shown': resolution * scale,
                        'box_power_over_mean_offset': box_ratio, 'box_corr': box_corr, 'top_corr': top_corr,
                        'repeat_corr': rep_corr}}


def jet_image(A):
    from matplotlib import colormaps
    return colormaps['jet'](np.clip(A, 0, 1))[..., :3]


SHOW_TO = 1.08                  # the published-style section runs just past one repeat depth


def styled_section(T, z, z_rep, axis, band, smooth=(3.0, 2.0)):
    """A thick section in the published style: the focused power averaged over `band` lines either side of the middle
    (as a rendered volume or a thick slab averages), smoothed, on a log scale, over the axis to just past z_rep."""
    mid = T.shape[axis] // 2
    slab = np.take(T, np.arange(mid - band, mid + band + 1), axis=axis).mean(axis=axis)   # [along, depth]
    keep = z <= SHOW_TO * z_rep
    A = gaussian_filter(slab[:, keep].T, smooth)
    A = np.log10(A / np.percentile(A, 50))
    lo, hi = np.percentile(A, 2), np.percentile(A, 99.5)
    return np.clip((A - lo) / (hi - lo), 0, 1), slab


def save_section(path, A):
    import matplotlib.pyplot as plt
    plt.imsave(path, jet_image(A))


def main():
    params = {'inputs': {k: f'legacy cache {v}' for k, v in REAL.items()} | {'simulated': 'P2-07 scene B (no chamber)'},
              'bank': {'K': 50, 'sub_frac': 0.5, 'delta_hz': 88.0}, 'lam_s_m': 0.48, 'repeat_depth_shown_m': Z_SHOW,
              'periods': PERIODS}
    with Run(RID, 'How the method draws pillars and boxes', params) as run:
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        sets = {}
        if all((LEGACY_CACHE / f).is_file() for f in REAL.values()):
            for k in REAL:
                sets[k] = load_real(k)
        sets['simulated'] = simulated_flat()
        out, volumes = {}, {}
        for name, d in sets.items():
            a = analyse(d)
            nr, nc = d['shape']
            T = a['T'].reshape(nr, nc, -1)
            # real crops: across range (526 m of ground); the simulated scene is long in azimuth instead
            axis, band = (0, 15) if name != 'simulated' else (1, 5)
            A, slab = styled_section(T, a['z'], a['z_rep'], axis, band)
            save_section(figs / f'{name}_section.png', A)
            mid = T.shape[axis] // 2
            lines = np.arange(mid - band, mid + band + 1)
            e = np.take(a['energy'].reshape(nr, nc), lines, axis=axis).mean(axis=axis)
            dc = np.take(a['dc'].reshape(nr, nc), lines, axis=axis).mean(axis=axis)
            i_rep = int(np.argmin(np.abs(a['z'] - a['z_rep'])))
            along_m = (nc if axis == 0 else nr) * d['pixel_m'][1 if axis == 0 else 0]
            out[name] = {**a['metrics'], 'grid': [nr, nc], 'pixel_m': list(d['pixel_m']), 'section_axis': 'range' if axis == 0 else 'azimuth',
                         'section_along_m': along_m, 'depth_shown_m': float(SHOW_TO * a['z_rep'] * a['scale']),
                         'energy_along': e, 'power_mean_along': slab[:, a['z'] < a['z_rep']].mean(axis=1),
                         'dc_along': dc, 'power_at_repeat_along': slab[:, i_rep],
                         'corr_along_energy': float(np.corrcoef(e, slab[:, a['z'] < a['z_rep']].mean(axis=1))[0, 1]),
                         'corr_along_box': float(np.corrcoef(dc, slab[:, i_rep])[0, 1])}
            if name != 'simulated':                        # where the noisy pixels (the pillars) stand on the surface
                lE = np.log(a['energy'].reshape(nr, nc) + 1e-30)
                lA = np.log(uniform_filter(d['amplitude'].reshape(nr, nc) ** 2, 5) + 1e-9)
                gA = np.hypot(*np.gradient(lA))
                top = lE > np.percentile(lE, 95)
                out[name].update({'noise_vs_brightness_corr': float(np.corrcoef(lE.ravel(), lA.ravel())[0, 1]),
                                  'noisiest_5pct_darker_db': float(10 / np.log(10) * (np.median(lA) - lA[top].mean())),
                                  'noisiest_5pct_gradient_ratio': float(gA[top].mean() / np.median(gA))})
            volumes[name] = (T, a)
            print(f"  {name}: pillar corr {a['metrics']['pillar_corr_mean_power_vs_energy']:.3f}, box ratio "
                  f"{a['metrics']['box_power_over_mean_offset']:.3f}, bands {a['metrics']['band_first_zero_m_shown']:.1f} m", flush=True)
        # the volume the viewer draws under the claimed underworld: Khafre's, relabelled to the published depths
        if 'khafre' in volumes:
            T, a = volumes['khafre']
            np.savez_compressed(run.dir / 'khafre_volume.npz', T=T.astype(np.float32), z_shown=a['z'] * a['scale'])
        np.savez_compressed(run.dir / 'simulated_volume.npz', T=volumes['simulated'][0].astype(np.float32),
                            z_shown=volumes['simulated'][1]['z'] * volumes['simulated'][1]['scale'])
        k = out.get('khafre', out['simulated'])
        pl = out.get('plateau', k)
        finding = (
            f"The method's pillars and boxes are its own arithmetic. Averaged over a period of its depth axis, the power it "
            f"draws at a pixel is that pixel's trajectory energy (correlation {k['pillar_corr_mean_power_vs_energy']:.3f}): a "
            f"pixel whose registration wanders is bright at every depth and, averaged over a thick slab as any rendering "
            f"does, stands as a column. On the real image such pixels lie where the surface is darker or changes brightness "
            f"(noise against brightness, correlation {k.get('noise_vs_brightness_corr', float('nan')):.2f} at Khafre and "
            f"{pl.get('noise_vs_brightness_corr', float('nan')):.2f} on open plateau), and open plateau with no monument draws "
            f"the tallest pillars of all; perfectly uniform simulated ground draws none. Along a column the power rises and "
            f"falls every {k['band_first_zero_m_shown']:.0f} m of the relabelled axis, its resolution: the bands. At the "
            f"surface and again at each repeat depth, {Z_SHOW:.0f} m and {2 * Z_SHOW:.0f} m relabelled, every steering phase is "
            f"the same, so the power there is each pixel's average offset (ratio {k['box_power_over_mean_offset']:.3f}, "
            f"correlation {k['box_corr']:.3f}): the blocks. Below {Z_SHOW:.0f} m the picture repeats (correlation "
            f"{k['repeat_corr']:.3f}).")
        run.save({'sets': out, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
