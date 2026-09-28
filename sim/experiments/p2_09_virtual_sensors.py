"""P2-09 · A million virtual sensors: can the image act as a dense seismic array?

    uv run python experiments/p2_09_virtual_sensors.py

The strongest form of the idea: every patch of a radar image is a virtual seismometer, far denser than any real array,
and stacking millions of them should reach motions no single sensor can. The first investigation tested stacking (R11)
but added its test wave to the readings after tracking, so it never asked whether the sensors hear. Here the motion goes
into the image itself.

1. The sensors: P2-07's simulated desert on the real Giza geometry, a sensor every 16 x 8 pixels, each a 32-pixel patch
   read the way the first investigation read a velocity series (R4): sixteen looks of 1.96 s against the centre look, by
   complex correlation (the published trackers) and by correlation of magnitudes.
2. Hearing: a test wave with spatial structure (60 m long, 0.15 Hz, 2 mm/s, unphysically slow so that the sensors see
   different motion) is put into the image; each sensor's reading is regressed on the true motion under it.
3. Stacking: with the ground still, the readings of n sensors are averaged; their scatter falls as n^-1/2 if the sensors
   are independent.
4. The count: how many sensors the chamber's imprint would need (P2-04), against how many independent ones exist over it,
   in the whole image, and on all of Earth's land.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

from katabasis.runs import Run, load
from sarsim import synthesize
from sarsim.acquisition import DwellGeometry
from sarsim.looks import look_masks, looks
from sarsim.track import patch_shifts

RID = 'p2_09_virtual_sensors'
_spec = importlib.util.spec_from_file_location('p207', Path(__file__).with_name('p2_07_whole_chain.py'))
p207 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p207)

STRIDE = (16, 8)
PATCH = 32
CENTRES = np.linspace(-0.93, 0.93, 16)
W = 1.96
WAVE = {'wavelength_m': 60.0, 'f_hz': 0.15, 'v_m_s': 2e-3, 'azimuth_deg': 30.0}
EARTH_LAND_M2 = 1.49e14


class TestWave:
    """A plane wave of vertical velocity with a short wavelength, so neighbouring sensors move differently."""

    def __init__(self, g):
        self.g = g
        k = 2 * np.pi / WAVE['wavelength_m']
        a = np.deg2rad(WAVE['azimuth_deg'])
        self.k = k * np.array([np.sin(a), np.cos(a)])
        self.w = 2 * np.pi * WAVE['f_hz']

    def up(self, east, north, t):
        """Vertical displacement (m) [len(t), len(east)]."""
        ph = np.outer(np.ones_like(t), east * self.k[0] + north * self.k[1]) - self.w * t[:, None]
        return (WAVE['v_m_s'] / self.w) * np.sin(-ph)

    def motion(self):
        los_up = self.g.los_enu[2]
        def f(east, north, t_abs):
            d = np.zeros(t_abs.shape + (3,))
            d[..., 2] = (WAVE['v_m_s'] / self.w) * np.sin(self.w * t_abs - (east * self.k[0] + north * self.k[1])[None, :])
            return d
        return self.g.absolute_motion(f)

    def los_velocity(self, east, north, t):
        return self.g.los_enu[2] * WAVE['v_m_s'] * np.cos(self.w * t[:, None] - (east * self.k[0] + north * self.k[1])[None, :])


def read(img, g, rows, cols, envelope):
    L = looks(img, look_masks(g, img.shape[0], CENTRES, W))
    if envelope:
        L = np.abs(L).astype(np.complex64)
    c = len(CENTRES) // 2
    out = np.zeros((rows.size, len(CENTRES)))
    for k in range(len(CENTRES)):
        dr, _ = patch_shifts(L[c], L[k], rows, cols, patch=PATCH, upsample=100)
        out[:, k] = g.shift_to_velocity(dr * g.dx)
    return out


def main():
    g = DwellGeometry.from_record('giza-20250827')
    rng = np.random.default_rng(71)
    scat = p207.scene(g, rng)
    rows = np.arange(64, p207.SHAPE[0] - 64, STRIDE[0])
    cols = np.arange(32, p207.SHAPE[1] - 32, STRIDE[1])
    RR, CC = np.meshgrid(rows, cols, indexing='ij')
    R, Cc = RR.ravel(), CC.ravel()
    xr = (R - p207.SHAPE[0] // 2) * g.dx
    yr = (Cc - p207.SHAPE[1] // 2) * g.dr / np.sin(g.theta)
    east = xr * g.along_track_en[0] + yr * g.ground_range_en[0]
    north = xr * g.along_track_en[1] + yr * g.ground_range_en[1]
    params = {'geometry': 'giza-20250827', 'sensor_stride_px': STRIDE, 'patch_px': PATCH, 'looks': len(CENTRES), 'look_s': W,
              'test_wave': WAVE, 'sensors': int(R.size)}
    with Run(RID, 'A million virtual sensors: can the image act as a dense seismic array?', params) as run:
        figs = run.dir / 'figs'
        figs.mkdir(exist_ok=True)
        still = synthesize(scat, g, p207.SHAPE)
        wave = TestWave(g)
        moving = synthesize(scat, g, p207.SHAPE, motion=wave.motion())
        truth = wave.los_velocity(east, north, CENTRES).T * np.sinc(WAVE['f_hz'] * W)      # [sensor, look], look-averaged
        truth_rel = truth - truth[:, [len(CENTRES) // 2]]
        res = {}
        for name, env in (('complex', False), ('magnitude', True)):
            s0 = read(still, g, R, Cc, env)
            s1 = read(moving, g, R, Cc, env)
            sig = s1 - s0                                                  # the reading's response to the wave alone
            gain = float(np.sum(sig * truth_rel) / np.sum(truth_rel ** 2))
            res[name] = {'still': s0, 'gain': gain, 'moving_minus_still': sig, 'moving': s1,
                         'noise_um_s': float(np.std(s0) * 1e6)}
            print(f'  {name}: gain {gain:.3f}, noise {np.std(s0) * 1e6:.0f} um/s', flush=True)
        # stacking: sensors whose patches do not overlap (every 2nd row and 4th column of the grid)
        nr, nc = len(rows), len(cols)
        idx = np.arange(nr * nc).reshape(nr, nc)[::2, ::4].ravel()
        stack = {}
        for name in res:
            s = res[name]['still'][idx]
            s = s - s.mean(axis=1, keepdims=True)
            ns = [1, 2, 4, 8, 16, 32, 64, len(idx)]
            ns = sorted(set(n for n in ns if n <= len(idx)))
            srng = np.random.default_rng(5)
            sd = []
            for n in ns:
                draws = [np.std(s[srng.choice(len(idx), n, replace=False)].mean(axis=0)) for _ in range(60)]
                sd.append(float(np.mean(draws)))
            slope = float(np.polyfit(np.log(ns), np.log(sd), 1)[0])
            stack[name] = {'n': ns, 'scatter_m_s': sd, 'slope': slope}
        # the count
        imp = load('p2_04_chamber_imprint')
        sees = load('p2_03_what_a_dwell_sees')
        signal = imp['cases'][0]['imprint_los_velocity_m_s']
        trembling = imp['cases'][0]['vertical_velocity_m_s']
        cell = g.look_resolution(W) * (0.886 / g.kr_band / np.sin(g.theta))                    # m^2, one independent look cell
        footprint = np.pi * 20.0 ** 2
        image = (g.source['shape'][0] * g.dx) * (g.source['shape'][1] * g.dr / np.sin(g.theta))
        looks_in_dwell = g.aperture_time / W
        r_real = sees['real']['block']['gain_envelope']
        sigma1 = res['magnitude']['noise_um_s'] * 1e-6 * np.sqrt(PATCH * g.dx * PATCH * g.dr / np.sin(g.theta) / cell)  # per cell
        need = lambda target, r: (sigma1 / (r * target)) ** 2 / looks_in_dwell
        count = {'independent_cell_m2': cell, 'in_footprint': footprint / cell, 'in_image': image / cell,
                 'on_earth_land': EARTH_LAND_M2 / cell, 'looks_in_dwell': looks_in_dwell,
                 'per_cell_noise_m_s': sigma1, 'magnitude_gain_real': r_real,
                 'needed_for_chamber': need(signal, r_real), 'needed_for_trembling': need(trembling, r_real),
                 'complex_gain': res['complex']['gain']}
        # every sensor's position, the true motion under it and what each way of reading reports, for the viewer
        np.savez_compressed(run.dir / 'sensors.npz', east=east, north=north, rows=R, cols=Cc, looks_s=CENTRES,
                            truth=truth_rel, complex=res['complex']['moving_minus_still'],
                            magnitude=res['magnitude']['moving_minus_still'], grid=np.array([len(rows), len(cols)]))
        # maps for the site: the true motion under each sensor at one look, and what each way of reading reports
        k = 2
        grid = lambda v: v.reshape(nr, nc)
        vmax = np.abs(truth_rel[:, k]).max()
        maps = {}
        for name, arr in (('truth', truth_rel[:, k]), ('complex', res['complex']['moving_minus_still'][:, k]),
                          ('magnitude', res['magnitude']['moving_minus_still'][:, k])):
            fig = p207.save(figs / f'sensors_{name}.png', grid(arr), p207.IMPRINT, -vmax, vmax)
            fig['aspect_m'] = [nc * STRIDE[1] * g.dr / np.sin(g.theta), nr * STRIDE[0] * g.dx]    # width x height of ground
            maps[name] = fig
        finding = (
            f"A radar image is not a dense seismic array. Made into {R.size:,} virtual sensors and shaken by a wave put into the "
            f"image itself, the sensors read by complex correlation, as the published method reads them, recover "
            f"{res['complex']['gain'] * 100:.0f}% of the motion under them, and read by magnitudes "
            f"{res['magnitude']['gain'] * 100:.0f}% (P2-03 found magnitudes catch up to about a fifth over large stretches of "
            f"real, point-like ground). Stacking still sensors does shrink their scatter as n^{stack['magnitude']['slope']:.2f}, "
            f"close to the n^-0.5 of independent sensors, but stacking a deaf sensor gives a precise zero. Where they do hear, "
            f"the chamber's imprint ({signal * 1e6:.1e} um/s) would need about {count['needed_for_chamber']:.0e} independent "
            f"sensors over it; there are about {count['in_footprint']:,.0f} over the chamber, {count['in_image']:.0e} in the whole "
            f"image and {count['on_earth_land']:.0e} on all of Earth's land. Even the ground's own trembling, with no chamber in "
            f"it, would need {count['needed_for_trembling']:.0e}, {count['needed_for_trembling'] / count['in_image']:,.0f} times "
            f"more than the whole image holds.")
        run.save({'sensors': {'n': int(R.size), 'rows': int(nr), 'cols': int(nc)},
                  'gains': {k2: v['gain'] for k2, v in res.items()}, 'noise_um_s': {k2: v['noise_um_s'] for k2, v in res.items()},
                  'stack': stack, 'count': count, 'maps': maps, 'map_look_s': float(CENTRES[k]), 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
