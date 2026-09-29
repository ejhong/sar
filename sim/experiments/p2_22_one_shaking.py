"""P2-22 · One shaking: a vibrator beside the bench, with the chamber and without it.

    uv run python experiments/p2_22_one_shaking.py

The decisive benchmark starts from one physical event, which geophones and the satellite will both read (P2-23). A
vertical vibrator on the surface, 30 m west of the chamber's axis, drives the bench at the gated reconstruction's own
tenth mode (10.69 Hz, from P2-20's record) until the ground's motion is steady. The lab's validated elastic solver
(P1-01) runs it three times on one grid: without a chamber, with the bench's chamber (a 6 m room 15 m down), and with a
more favourable one (a 10 m room whose roof is 5 m down). Recorded: the steady complex velocity of every surface point
within 60 m of the chamber's axis, east, north and up, per newton of force. The solver is linear, so any force scales it.

Checks written with the run: the motion is steady (two successive three-cycle windows agree), and the chamber's
imprint (with minus without) is set against the motion itself. The absorbing boundaries are 20 m thick, an eighth of a
surface wavelength, so the steady field carries some reflection from them; it is the same with and without a chamber.
"""
import copy
import json
import time
from pathlib import Path

import numpy as np

from katabasis.compose import Grid, load_site, voxelise
from katabasis.compose.site import parse_site
from katabasis.runs import RESULTS, Run
from katabasis.seismic.elastic3d import Medium, Receivers, Simulation, Source

RID = 'p2_22_one_shaking'
SITE = 'bench-void'
MODE = 10
SOURCE = (-30.0, 0.0, 0.0)
EXTENT = ((-150.0, 150.0), (-150.0, 150.0), (-100.0, 3.0))
H = 1.0
PML = 20
HALF = 60.0                          # the recorded surface: within this of the chamber's axis, every metre
CYCLES_ON, CYCLES = 3, 12            # a raised-cosine start over three cycles; steady windows are cycles 6-9 and 9-12
RUN_FORCE_N = 1e6                    # run at 1 MN for single precision's sake; fields are reported per newton
RECORD_EVERY = 8
CASES = {
    'none': None,
    'bench': {'centre': [0, 0, -15], 'size': [6, 6, 6]},
    'favourable': {'centre': [0, 0, -10], 'size': [10, 10, 10]},
}


def model(site, chamber):
    raw = copy.deepcopy(site.raw)
    if chamber is None:
        raw['features'] = []
    else:
        raw['features'][0]['shape'] = dict(raw['features'][0]['shape'], centre=chamber['centre'], size=chamber['size'])
    return parse_site(raw, site.directory)


def window_amplitude(traces, t, f, t0, t1):
    """The complex amplitude A of each trace's steady sinusoid, v = Re(A exp(i 2 pi f t)), over [t0, t1)."""
    sel = (t >= t0) & (t < t1)
    e = np.exp(-2j * np.pi * f * t[sel])
    return 2 * (traces[..., sel] * e).mean(axis=-1)


def main():
    p20 = json.loads((RESULTS / 'p2_20_bench_gated' / 'summary.json').read_text())
    f = MODE * p20['mode_hz'] / 2                           # P2-20 records the second mode
    site = load_site(SITE)
    g = Grid.covering(*EXTENT, H)
    xs = np.arange(-HALF, HALF + 1e-9, H)
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    rec = Receivers(np.stack([X.ravel(), Y.ravel(), np.zeros(X.size)], 1))
    params = {'site': SITE, 'frequency_hz': f, 'mode': MODE, 'source_m': SOURCE, 'extent_m': EXTENT, 'spacing_m': H,
              'pml_cells': PML, 'cycles': [CYCLES_ON, CYCLES], 'cases': CASES, 'recorded_half_m': HALF}
    with Run(RID, 'One shaking: a vibrator beside the bench, with the chamber and without it', params) as run:
        fields, checks = {}, {}
        for name, chamber in CASES.items():
            t0 = time.time()
            med = Medium.from_model(voxelise(model(site, chamber), g))
            sim = Simulation(med, pml_width=PML, f0=f)
            nt = int(CYCLES / f / sim.dt)
            tt = np.arange(nt) * sim.dt
            ramp = np.where(tt < CYCLES_ON / f, 0.5 * (1 - np.cos(np.pi * tt * f / CYCLES_ON)), 1.0)
            w = (RUN_FORCE_N * ramp * np.sin(2 * np.pi * f * tt)).astype(np.float32)
            res = sim.run([Source(SOURCE, w, 'force', (0, 0, -1))], rec, nt, record_every=RECORD_EVERY)
            a = window_amplitude(res.traces, res.t, f, 6 / f, 9 / f) / RUN_FORCE_N
            b = window_amplitude(res.traces, res.t, f, 9 / f, 12 / f) / RUN_FORCE_N
            fields[name] = b.reshape(len(xs), len(xs), 3)
            up = np.abs(b[:, 2])
            checks[name] = {'steady_change_max': float(np.max(np.abs(b - a)) / np.max(np.abs(b))),
                            'steady_change_median': float(np.median(np.abs(b[:, 2] - a[:, 2]) / np.maximum(up, 1e-30))),
                            'up_at_axis_m_s_per_n': float(np.abs(fields[name][len(xs) // 2, len(xs) // 2, 2])),
                            'runtime_s': round(time.time() - t0, 1)}
            print(f"  {name}: {checks[name]}", flush=True)

        # the chambers' imprints, against the motion without a chamber
        U0 = fields['none']
        R = np.hypot(X, Y)
        imprint = {}
        for name in ('bench', 'favourable'):
            d = fields[name] - U0
            mag = np.linalg.norm(d, axis=2)
            ref = np.linalg.norm(U0, axis=2)
            imprint[name] = {'max_m_s_per_n': float(mag.max()), 'max_relative': float((mag / ref).max()),
                             'relative_at_axis': float(mag[len(xs) // 2, len(xs) // 2] / ref[len(xs) // 2, len(xs) // 2]),
                             'half_max_radius_m': float(R[mag >= 0.5 * mag.max()].max()),
                             'up_max_m_s_per_n': float(np.abs(d[..., 2]).max())}
        np.savez_compressed(Path(run.dir) / 'fields.npz', x=xs, y=xs, f_hz=f, source=np.array(SOURCE),
                            **{f'U_{k}': v.astype(np.complex64) for k, v in fields.items()})
        ib, iff = imprint['bench'], imprint['favourable']
        finding = (
            f"A vertical vibrator 30 m west of the chamber's axis, at {f:.2f} Hz, moves the ground over the axis by "
            f"{checks['none']['up_at_axis_m_s_per_n'] * 1e9:.2f} nm/s per newton of force. The bench's chamber (6 m, 15 m "
            f"down) changes the surface motion by at most {ib['max_relative'] * 100:.2f}% of it ({ib['max_m_s_per_n'] * 1e9:.3f} "
            f"nm/s per newton), within {ib['half_max_radius_m']:.0f} m of its axis at half that; the favourable chamber (10 m, "
            f"roof 5 m down) by {iff['max_relative'] * 100:.1f}% ({iff['max_m_s_per_n'] * 1e9:.2f} nm/s per newton). The motion "
            f"is steady: successive windows differ by at most {max(c['steady_change_max'] for c in checks.values()) * 100:.1f}% "
            f"of its largest value.")
        run.save({'checks': checks, 'imprint': imprint, 'frequency_hz': f, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
