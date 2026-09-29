"""P2-23 · One shaking, read by geophones and by the satellite: where each finds the chamber.

    uv run --with scikit-image python experiments/p2_23_every_reader.py

P2-22's vibrator shakes the one-chamber bench with its chamber, without it, and with a more favourable one. Here that
one steady motion is read over the same 24.6 s, the 2022 dwell's aperture, by

  geophones      49 vertical geophones 8 m apart over the chamber, in the ambient noise measured at Kottamya (M1-01:
                 its spectrum at 7.7 Hz, the highest frequency it resolves, held flat to the vibrator's 10.69 Hz) and,
                 for a noisy site, Peterson's high-noise model;
  the satellite  the 2022 dwell's geometry, as synthetic products: corner reflectors 50 dB over the ground, laid every
                 3 m of ground range across the chamber (a field test's best case: a look short enough for 10.69 Hz,
                 20 ms, blurs tens of metres along the track, so one row across it is what the looks can separate), and
                 natural ground carrying the bright points the real 2022 image holds over open plateau.

The standard readings are the amplitude of each geophone's record at the vibrator's frequency, and each bright point's
line-of-sight velocity tracked through 20 ms looks of the whole dwell (P2-19's tracker), calibrated here on reflectors
of known motion. The favourable detector knows the vibrator and both candidate answers, the chamber and no chamber, and
asks which the readings fit (the likelihood ratio of two known signals in Gaussian noise): the best any method could do
with these observations. The reflectors are simulated moving with the ground; the speckle around them is held still,
since it carries no usable motion (P2-03). P2-24 runs the gated reconstruction on shaken natural ground.

Stated before the run: the chamber counts as found when the favourable detector tells it from no chamber with 90%
power at a 1% false-alarm rate; each reader's boundary is the smallest vibrator force at which it does. The same force
gives both instruments the same ground motion, so the ratio of their boundaries is what the instrument costs.
"""
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import maximum_filter, median_filter
from scipy.stats import norm

import p2_13_gates as m
from katabasis.ambient.peterson import P_MIN, acceleration_db
from katabasis.lab import PRODUCT
from katabasis.runs import RESULTS, Run
from sarsim.acquisition import DwellGeometry
from sarsim.looks import frequencies, look_masks, region_shift

RID = 'p2_23_every_reader'
ALPHA, POWER = 0.01, 0.90
Z = norm.ppf(1 - ALPHA) + norm.ppf(POWER)          # a deflection this large is found at that power and false-alarm rate
GEOPHONES = np.array([(x, y) for x in range(-24, 25, 8) for y in range(-24, 25, 8)], float)
ROWS, COLS = 8192, 640
LOOK_S = 0.020
SUB_ROWS = 4096                                     # rows the tracker reads about each point
CAL_SCR_DB = [30.0, 40.0, 50.0]
CAL_V = 0.05                                        # m/s: calibration reflectors vibrate this fast along the line of sight
CAL_SEEDS = [230, 231, 232]
REFLECTOR_SCR_DB = 50.0
REFLECTOR_STEP_M = 3.0
REFLECTOR_HALF_M = 30.0
PLATEAU = (-470.0, -600.0)                          # P2-17's open plateau, where the real image's bright points are counted
POP_HALF_M = 60.0
POP_THRESH_DB = [10.0, 15.0, 20.0, 25.0, 30.0, 40.0]
FORCES = np.logspace(2, 12, 41)


def fields():
    d = np.load(RESULTS / 'p2_22_one_shaking' / 'fields.npz')
    x, y = d['x'], d['y']
    interp = {}
    for k in ('none', 'bench', 'favourable'):
        U = d[f'U_{k}']
        interp[k] = [RegularGridInterpolator((x, y), U[..., c]) for c in range(3)]
    return float(d['f_hz']), interp


def sample(interp, pts):
    return np.stack([f(pts) for f in interp], axis=-1)            # complex (n, 3): east, north, up per newton


def noise_asd(f):
    """Velocity amplitude spectral density (m/s/rtHz) of the ambient noise: Kottamya's median (M1-01) at the shortest
    period it resolves, 0.13 s (20 samples a second), and Peterson's high-noise model at its shortest period, 0.1 s,
    each with its acceleration held flat to f."""
    s = json.loads((RESULTS / 'm1_01_ambient_levels' / 'summary.json').read_text())
    c = s['components']['BHZ']
    T = np.array(c['period_s'])
    i = int(np.flatnonzero(T >= 0.125)[0])
    acc = 10 ** (np.array(c['acc_db']['all_p50'], float)[i] / 10)
    return {'kottamya_median': float(np.sqrt(acc) / (2 * np.pi * f)), 'kottamya_period_s': float(T[i]),
            'peterson_high': float(np.sqrt(10 ** (acceleration_db(np.array([P_MIN]), 'high')[0] / 10)) / (2 * np.pi * f)),
            'peterson_period_s': float(P_MIN)}


# ------------------------------------------------------------------ the satellite's products and its tracker --

def band_mask(g):
    fa = np.fft.fftfreq(ROWS, d=1.0 / g.prf)
    kr = np.fft.fftfreq(COLS, d=g.dr)
    return (np.abs(fa) <= g.doppler_band_hz / 2)[:, None] & (np.abs(kr) <= g.kr_band / 2)[None, :]


def speckle(g, rng, mask):
    z = (rng.standard_normal((ROWS, COLS)) + 1j * rng.standard_normal((ROWS, COLS))) / np.sqrt(2)
    out = np.fft.ifft2(np.fft.fft2(z) * mask)
    return out / np.sqrt(mask.mean())                             # mean intensity 1


_KERNEL = {}


def point_image(mask, r, c, scr_db):
    """A point scatterer at pixel (r, c), band-limited like the speckle, its peak intensity scr over the speckle's mean:
    one kernel, shifted."""
    if 'k' not in _KERNEL:
        z = np.zeros((ROWS, COLS), np.complex128)
        z[0, 0] = 1.0
        k = np.fft.ifft2(np.fft.fft2(z) * mask)
        _KERNEL['k'] = k / np.abs(k).max()
    return np.roll(_KERNEL['k'], (r, c), axis=(0, 1)) * np.sqrt(10 ** (scr_db / 10))


def vibrate(img, g, A_los, f):
    """The whole of `img` (one point's image) displaced along the line of sight by Re(D exp(i 2 pi f t)), D = A / (i 2 pi f)
    from the complex velocity amplitude A: the exact phase history, through the Doppler-time mapping t = f_D / Ka."""
    fd = np.fft.fftfreq(ROWS, d=1.0 / g.prf)
    t = fd / g.Ka_signed
    D = A_los / (2j * np.pi * f)
    d = np.real(D * np.exp(2j * np.pi * f * t))
    return np.fft.ifft(np.fft.fft(img, axis=0) * np.exp(-4j * np.pi * d / g.lam)[:, None], axis=0)


def smear_rows(g):
    """Rows over which a 20 ms look spreads a point along the track: the image's sampling over the look's bandwidth."""
    return int(np.ceil(g.prf / (g.Ka * LOOK_S)))


def track(img, g, f, r, c, centres):
    """The line-of-sight velocity of the point at (r, c) through 20 ms looks, and its complex amplitude at f (with a mean
    and a drift fitted beside it) and that amplitude's standard error per quadrature."""
    rs = slice(r - SUB_ROWS // 2, r + SUB_ROWS // 2)
    cs = slice(c - 4, c + 5)
    sub = img[rs, cs]
    S = np.fft.fftshift(np.fft.fft(sub, axis=0), axes=0)
    masks = look_masks(g, SUB_ROWS, centres, LOOK_S)
    look = lambda k: np.fft.ifft(np.fft.ifftshift(S * masks[k][:, None], axes=0), axis=0)
    ref = look(len(centres) // 2)
    half = min(int(0.6 * smear_rows(g)), SUB_ROWS // 2 - 8)
    rows = slice(SUB_ROWS // 2 - half, SUB_ROWS // 2 + half)
    cols = slice(3, 6)
    v = np.array([g.shift_to_velocity(region_shift(ref, look(k), rows, cols, envelope=True) * g.dx) for k in range(len(centres))])
    w = 2 * np.pi * f
    M = np.stack([np.cos(w * centres), -np.sin(w * centres), np.ones_like(centres), centres], 1)
    coef, *_ = np.linalg.lstsq(M, v, rcond=None)
    resid = v - M @ coef
    se = float(np.std(resid) * np.sqrt(2.0 / len(centres)))
    return complex(coef[0], coef[1]), se


def look_centres(g):
    half = g.doppler_band_hz / 2 / g.Ka - LOOK_S
    return np.arange(-half, half + 1e-9, LOOK_S)


# ------------------------------------------------------------------ the real image's bright points --

def bright_population(product_path):
    """Bright points in the real 2022 image over open plateau: local peaks whose intensity stands over the local speckle
    mean (the local median over ln 2) by each threshold, per square metre of ground."""
    from katabasis.compose import load_site
    from katabasis.export.sites import scene
    from sarsim.dwell import DwellProduct
    from sarsim.ortho import frame_to_lla, surface
    sc = scene(load_site('giza'))
    g = DwellGeometry.from_record(m.ACQ)
    p = DwellProduct(product_path)
    x, y = PLATEAU
    z = float(surface(sc, np.array([x]), np.array([y]))[0])
    lat, lon = frame_to_lla(sc['frame']['origin'], np.array([x]), np.array([y]))
    r, c = p.geolocate(float(lat[0]), float(lon[0]), z + 15.5)
    n_rows = int(2 * POP_HALF_M / g.dx)
    n_cols = int(2 * POP_HALF_M * np.sin(g.theta) / g.dr)
    crop, _ = p.crop(r, c, n_rows, n_cols)
    I = np.abs(crop.astype(np.complex64)) ** 2
    mean = median_filter(I, size=(201, 15)) / np.log(2)
    peak = (I == maximum_filter(I, size=(41, 5)))
    ratio_db = 10 * np.log10(np.maximum(I / np.maximum(mean, 1e-30), 1e-30))
    area = (2 * POP_HALF_M) ** 2
    scr = np.sort(ratio_db[peak & (ratio_db >= min(POP_THRESH_DB))])[::-1]
    # what speckle alone gives: pixels whose exponential intensity exceeds each threshold times its mean
    by_speckle = {f'{t:g}': float(I.size * np.exp(-10 ** (t / 10))) for t in POP_THRESH_DB}
    return {'per_m2': {f'{t:g}': float((scr >= t).sum() / area) for t in POP_THRESH_DB},
            'counts': {f'{t:g}': int((scr >= t).sum()) for t in POP_THRESH_DB}, 'speckle_alone': by_speckle,
            'brightest_db': float(scr[0]) if len(scr) else None, 'scr_db': scr.tolist(), 'area_m2': area,
            'crop_px': [n_rows, n_cols]}


def main():
    f, U = fields()
    g = DwellGeometry.from_record(m.ACQ)
    T = float(g.aperture_time)
    asd = noise_asd(f)
    los = np.asarray(g.los_enu, float)
    params = {'acquisition': m.ACQ, 'frequency_hz': f, 'record_s': T, 'alpha': ALPHA, 'power': POWER,
              'geophones': GEOPHONES.tolist(), 'noise_asd_m_s_rthz': asd, 'look_s': LOOK_S, 'cal_scr_db': CAL_SCR_DB,
              'cal_v_m_s': CAL_V, 'cal_seeds': CAL_SEEDS, 'reflector_scr_db': REFLECTOR_SCR_DB,
              'reflector_step_m': REFLECTOR_STEP_M, 'reflector_half_m': REFLECTOR_HALF_M, 'plateau_m': PLATEAU,
              'p2_22_commit': json.loads((RESULTS / 'p2_22_one_shaking' / 'summary.json').read_text())['manifest']['commit']}
    with Run(RID, 'One shaking, read by geophones and by the satellite', params) as run:
        rng = np.random.default_rng(23)
        out = {'noise': asd}
        axis_up = float(np.abs(sample(U['none'], np.array([[0.0, 0.0]]))[0, 2]))

        # --- geophones: the amplitude at f of each record; the favourable detector's deflection per newton
        geo = {}
        for name in ('kottamya_median', 'peterson_high'):
            if not np.isfinite(asd[name]):
                raise ValueError(f'no noise level for {name}')
            sigma = asd[name] / np.sqrt(T)                                       # per quadrature, from a T-second record
            geo[name] = {'sigma_m_s': float(sigma)}
            for ch in ('bench', 'favourable'):
                delta = (sample(U[ch], GEOPHONES) - sample(U['none'], GEOPHONES))[:, 2]
                per_n = float(np.linalg.norm(delta) / sigma)
                f90 = Z / per_n
                # the standard reading at that force: where the residual against no chamber is largest
                far = float(np.hypot(*GEOPHONES[np.argmax(np.abs(delta))]))
                geo[name][ch] = {'deflection_per_n': per_n, 'force_n': f90, 'ground_up_m_s': f90 * axis_up,
                                 'imprint_max_m_s': f90 * float(np.abs(delta).max()), 'largest_residual_from_axis_m': far}
            # a check of the detector's power by drawing the noise (the bench's chamber, at its boundary force)
            delta = (sample(U['bench'], GEOPHONES) - sample(U['none'], GEOPHONES))[:, 2] * geo[name]['bench']['force_n']
            s0 = np.zeros_like(delta)
            draws = []
            for truth in (0, 1):
                n = rng.standard_normal((4000, len(delta))) * sigma + 1j * rng.standard_normal((4000, len(delta))) * sigma
                y = s0 + truth * delta + n
                draws.append(np.real((y - delta / 2) @ np.conj(delta)) / (sigma * np.linalg.norm(delta)))
            thr = np.quantile(draws[0], 1 - ALPHA)
            geo[name]['check_power'] = float(np.mean(draws[1] > thr))
            print(f"  geophones, {name}: {geo[name]}", flush=True)
        out['geophones'] = geo

        # --- the satellite: calibrate the tracker on reflectors of known motion
        mask = band_mask(g)
        centres = look_centres(g)
        gain_true = float(np.sinc(f * LOOK_S))
        cal = {s: [] for s in CAL_SCR_DB}
        cols = np.arange(40, COLS - 40, 24)
        for seed in CAL_SEEDS:
            r2 = np.random.default_rng(seed)
            img = speckle(g, r2, mask)
            truth = []
            for j, c in enumerate(cols):
                scr = CAL_SCR_DB[j % len(CAL_SCR_DB)]
                A = CAL_V * np.exp(2j * np.pi * r2.random())
                img = img + vibrate(point_image(mask, ROWS // 2, int(c), scr), g, A, f)
                truth.append((int(c), scr, A))
            for c, scr, A in truth:
                est, se = track(img, g, f, ROWS // 2, c, centres)
                cal[scr].append({'ratio': est / A, 'se': se})
        # each brightness: the tracker's gain (the brightest fix it best) and its noise per quadrature, from each fit's
        # own residuals and, as a check, from the scatter of the recovered amplitudes about that gain
        calib = {}
        for scr, rows in cal.items():
            ratio = np.array([x['ratio'] for x in rows])
            gain = complex(np.median(ratio.real), np.median(ratio.imag))
            calib[f'{scr:g}'] = {'gain': [gain.real, gain.imag], 'sigma_m_s': float(np.median([x['se'] for x in rows])),
                                 'scatter_m_s': float(np.sqrt(np.mean(np.abs(ratio - gain) ** 2) / 2) * CAL_V), 'n': len(rows)}
        out['tracker'] = {'look_gain_sinc': gain_true, 'calibration': calib, 'smear_rows': smear_rows(g),
                          'smear_m': smear_rows(g) * g.dx}
        print(f"  tracker: {out['tracker']}", flush=True)
        # the tracker's noise against brightness, a straight line in log noise, extended beyond the calibrated range
        fit = np.polyfit(CAL_SCR_DB, np.log10([calib[f'{s:g}']['sigma_m_s'] for s in CAL_SCR_DB]), 1)
        out['tracker']['noise_db_slope'] = float(fit[0] * 20)
        sig = lambda scr_db: float(10 ** np.polyval(fit, scr_db))
        G = complex(*calib[f'{REFLECTOR_SCR_DB:g}']['gain'])

        # --- the field test: reflectors across the chamber in ground range, on its own row
        s_off = np.arange(-REFLECTOR_HALF_M, REFLECTOR_HALF_M + 1e-9, REFLECTOR_STEP_M)
        pts = np.stack([s_off * g.ground_range_en[0], s_off * g.ground_range_en[1]], 1)
        cols_r = (COLS // 2 + np.round(s_off * np.sin(g.theta) / g.dr)).astype(int)
        v_los = {k: sample(U[k], pts) @ los for k in ('none', 'bench', 'favourable')}
        sat = {'reflectors': {'count': len(pts), 'sigma_m_s': sig(REFLECTOR_SCR_DB)}}
        for ch in ('bench', 'favourable'):
            delta = G * (v_los[ch] - v_los['none'])
            per_n = float(np.linalg.norm(delta) / sig(REFLECTOR_SCR_DB))
            f90 = Z / per_n
            sat['reflectors'][ch] = {'deflection_per_n': per_n, 'force_n': f90, 'ground_up_m_s': f90 * axis_up,
                                     'imprint_max_los_m_s': f90 * float(np.abs(v_los[ch] - v_los['none']).max()),
                                     'direct_los_max_m_s': f90 * float(np.abs(v_los['none']).max())}

        # the whole chain, twice: at a force whose motion lies within the tracker's calibration (the measurement control:
        # it must recover the vibrator's own motion), and at the idealised boundary above, where the ground moves far faster
        F_bound = sat['reflectors']['bench']['force_n']
        F_track = CAL_V / float(np.abs(v_los['none']).max())
        chain = {}
        for label, F in (('within_tracking', F_track), ('at_boundary', F_bound)):
            d = float(Z * F / F_bound)                                    # the detector's deflection at this force
            chain[label] = {'force_n': F, 'direct_los_max_m_s': F * float(np.abs(v_los['none']).max()), 'deflection': d}
            for truth in ('none', 'bench'):
                r3 = np.random.default_rng(2323 if truth == 'none' else 2324)
                img = speckle(g, r3, mask)
                for c, A in zip(cols_r, v_los[truth] * F):
                    img = img + vibrate(point_image(mask, ROWS // 2, int(c), REFLECTOR_SCR_DB), g, A, f)
                est = np.array([track(img, g, f, ROWS // 2, int(c), centres)[0] for c in cols_r])
                m0, m1 = G * v_los['none'] * F, G * v_los['bench'] * F
                stat = float(np.real(np.vdot(m1 - m0, est - (m0 + m1) / 2)) / (sig(REFLECTOR_SCR_DB) * np.linalg.norm(m1 - m0)))
                control = float(np.real(np.vdot(m0, est)) / np.vdot(m0, m0).real)
                chain[label][truth] = {'statistic': stat, 'threshold': float(norm.ppf(1 - ALPHA) - d / 2),
                                       'direct_motion_recovered': control}
                print(f"  chain, {label}, {truth}: {chain[label][truth]}", flush=True)
        sat['chain'] = chain

        # --- natural ground: the real image's bright points, the tracker's noise for each, blended along the track
        pop = bright_population(PRODUCT)
        sat['natural'] = {'per_m2': pop['per_m2'], 'counts': pop['counts'], 'speckle_alone': pop['speckle_alone'],
                          'brightest_db': pop['brightest_db'], 'area_m2': pop['area_m2'],
                          'scr_db': [round(x, 1) for x in pop['scr_db'] if x >= min(CAL_SCR_DB)]}
        density = pop['per_m2'][f'{min(CAL_SCR_DB):g}']
        n_pts = rng.poisson(density * (2 * 60.0) ** 2)
        scr_pool = np.array(pop['scr_db'])
        scr_pool = scr_pool[scr_pool >= min(CAL_SCR_DB)]
        xy = rng.uniform(-60.0, 60.0, (n_pts, 2))
        scr = rng.choice(scr_pool, n_pts) if len(scr_pool) else np.zeros(0)
        # points that share a range cell and lie within one look's smear along the track blend into one reading
        a_m = xy @ np.asarray(g.along_track_en)
        col = np.round((xy @ np.asarray(g.ground_range_en)) * np.sin(g.theta) / g.dr).astype(int)
        smear_m = g.dx * smear_rows(g)
        groups = {}
        for i in np.argsort(-scr):
            key = next((k for k in groups if k[0] == col[i] and abs(a_m[k[1]] - a_m[i]) < smear_m), None)
            groups.setdefault(key or (col[i], i), []).append(i)
        for ch in ('bench', 'favourable'):
            d_all = (sample(U[ch], xy) - sample(U['none'], xy)) @ los
            dvec, svec = [], []
            for members in groups.values():
                wts = 10 ** (scr[members] / 10)
                dvec.append(G * np.sum(wts * d_all[members]) / wts.sum())
                svec.append(sig(float(scr[members].max())))
            per_n = float(np.sqrt(np.sum(np.abs(np.array(dvec)) ** 2 / np.array(svec) ** 2))) if dvec else 0.0
            f90 = Z / per_n if per_n > 0 else None                        # None: nothing to track, at any force
            sat['natural'][ch] = {'points': int(n_pts), 'readings': len(groups), 'deflection_per_n': per_n, 'force_n': f90,
                                  'ground_up_m_s': f90 * axis_up if f90 else None}
        out['satellite'] = sat
        print(f"  satellite: {sat}", flush=True)

        # the boundaries on one axis of force, and the power curves along it
        curves = {}
        for label, per_n in (('geophones, quiet site', geo['kottamya_median']['bench']['deflection_per_n']),
                             ('geophones, noisy site', geo['peterson_high']['bench']['deflection_per_n']),
                             ('satellite, reflectors', sat['reflectors']['bench']['deflection_per_n']),
                             ('satellite, natural ground', sat['natural']['bench']['deflection_per_n'])):
            curves[label] = norm.sf(norm.ppf(1 - ALPHA) - FORCES * per_n).tolist()
        np.savez_compressed(Path(run.dir) / 'curves.npz', forces=FORCES, **{k.replace(', ', '_').replace(' ', '_'): np.array(v) for k, v in curves.items()})

        gq, gn = geo['kottamya_median']['bench'], geo['peterson_high']['bench']
        rb, nb = sat['reflectors']['bench'], sat['natural']['bench']
        ct, cb = chain['within_tracking'], chain['at_boundary']
        finding = (
            f"One vibrator at {f:.2f} Hz, 30 m from the bench's chamber, read over the same {T:.1f} s. Geophones over the "
            f"chamber find it with {gq['force_n'] / 1e3:,.1f} kN of force at a quiet site and {gn['force_n'] / 1e3:,.0f} kN at a "
            f"noisy one, when the ground over it moves at {gq['ground_up_m_s'] * 1e9:,.0f} and {gn['ground_up_m_s'] * 1e6:,.1f} "
            f"um/s. The satellite, with {len(pts)} corner reflectors laid across the chamber, needs {rb['force_n']:.2e} N, "
            f"the ground moving at {rb['ground_up_m_s'] * 1e3:,.0f} mm/s: {rb['force_n'] / gq['force_n']:.1e} times the "
            f"geophones' force. "
            + (f"On natural ground, with the bright points the real image holds, {nb['force_n']:.1e} N. " if nb['force_n']
               else f"On natural ground there is nothing to track: over open plateau the real image's brightest point stands "
                    f"{sat['natural']['brightest_db']:.0f} dB over the speckle, and its {sat['natural']['counts']['10']} points of 10 dB "
                    f"or more are what speckle alone gives ({sat['natural']['speckle_alone']['10']:.0f}). ")
            + f"Run whole, the chain recovers {ct['none']['direct_motion_recovered'] * 100:.0f}% of the vibrator's own motion at "
            f"{ct['force_n']:.1e} N, where the reflectors move at up to {ct['direct_los_max_m_s'] * 1e3:.0f} mm/s (the measurement "
            f"control), while the chamber stays out of reach there; at the boundary the reflectors move at up to "
            f"{cb['direct_los_max_m_s']:.2f} m/s and the tracker recovers {cb['none']['direct_motion_recovered'] * 100:.0f}% of it, "
            f"so that boundary is a floor.")
        run.save(dict(out, forces=FORCES.tolist(), curves=curves, axis_up_m_s_per_n=axis_up, finding=finding))
        print(finding)


if __name__ == '__main__':
    main()
