"""P2-27 · What one image would need: the motion a room must add to the ground, under stated conditions.

    uv run python experiments/p2_27_requirement.py

P2-25 bounds what one image holds of the bench chamber under the ambient motion measured or estimated for Giza. This
turns the bound around: for a stated acquisition, scattering model, footprint and performance target, how large must the
motion the room itself adds to the ground be before any way of reading one image could find it? A requirement curve,
conditional on everything stated here. It does not exclude mechanisms that are not specified, nor routes outside the
model.

Conditions.
- Acquisition: the Giza dwell of 27 August 2025 (P2-01); the Giza acquisition of 15 July 2022 for comparison.
- Scattering: fully developed speckle, scatterers white on the image's own pixels; no receiver noise and no reference
  image (both could only raise what the image holds); the model's Doppler-to-time relation (P2-25 checks it against
  pulse-by-pulse physics: where the two differ, the model holds more).
- The room-specific motion: the part of the surface's motion that is there because the room is, along the line of sight,
  D exp(-r^2 / 2 L^2) cos(2 pi f (t + x_a / V_g) + psi), coherent for the whole pass, at the worst phase psi. Its
  footprint is the area within half its peak, 2 pi ln 2 L^2. These curves hold for this family of smooth patterns only:
  what an image holds depends on the motion's spatial pattern, not only its amplitude, frequency and area, and a pattern
  with more structure inside the same envelope can hold far more (computed below: a Gaussian 4 m wide at 0.2 Hz
  multiplied by a cosine of 4 m period holds 39 times the information, though nowhere larger; an independent review's
  example). The requirement for every pattern is on the motion's energy (5). Outside: changes of the surface's
  reflectivity during the pass, the echo from beneath (P2-11), scatterers that are not speckle, more than one image.
- The target: detection rate 95% at false-alarm rate 5% (detection minus false alarm, TV, 0.9), and TV 0.5.

Computed.
1. Necessary, exact within the model: below this motion no test reaches the target. The image's range band is widened to
   its whole pixel grid, which can only add information (the narrower image is a function of the wider), and which splits
   the image into independent lines along track. Each line's KL for the finite change is computed exactly, from the
   determinant of its covariance (no linearisation, no remainder), tabulated against the line's amplitude, and summed over
   the lines. Footprints up to 8 m wide (16 m at 5 Hz and below).
2. Linearised: KL = F / 2 at the worst phase, F exact on the image's own pixels (fisher_grid) or by the continuum
   (fisher_sinusoid, with its bound on the phase swing) for wide footprints, to 300 m; and what the best test needs by
   weak-signal theory (deflection 3.29 for TV 0.9), an estimate, not a bound.
3. Helped: the method handed the ground's reflectivity, every pixel at the calibrated SNR over receiver noise (sarsim.radiometry), KL <= SNR sum_cells
   <Phi^2>; and a corner reflector (50 dB, P2-23) on the footprint's peak, KL <= SCR (k0 D)^2 / 2.
5. For every pattern and time history, on the displacement integral E = integral max_t |du_los|^2 dA (m^4: an integrated
   squared displacement envelope, not a mechanical energy). The oracle of P2-29, told the reflectivity, the excitation
   and both models, with receiver noise at the calibrated SNR per cell: its mean echoes differ by exp(i k u_A) - exp(i k u_0), so
   motion common to both worlds cancels whatever excites it, and E Delta^2 <= 2 snr cells k0^2 E; it bounds every method
   told less, under any shaking. Theorem C (P2-30), fully developed speckle without receiver noise, holds only where the
   no-cavity world's motion is the background whose covariance floor is proved: quiet ground. Both at three targets:
   near chance (TV 0.05), even odds (0.5), 95% found at 5% false alarms (0.9); at every operating point found <= false
   alarms + TV. A change of the ambient motion's phase alone moves the ground by at most twice its amplitude a; to escape
   this exclusion a phase-only mechanism must change more than E_min / (2 a)^2 of ground, under the stated amplitude and
   noise.
6. How much the pattern matters: the Gaussian's F as a share of 4 N <Phi^2>, the general upper bound for its envelope
   (sarsim.information.fisher_bound; not a demonstrated optimum); and the review's modulated Gaussians, exactly on the
   image's pixels.
4. Beside the curve: the ambient motion itself (Giza's regional microseism level, measured 67 km east; the noisiest
   stations on Earth; the measured 1-3 and 3-8 Hz levels; the FTA's urban background and a truck over a bump 15 m away),
   none of them an upper limit at the pyramids; the largest room-specific motion per unit ambient motion in the rooms
   P2-26 models, a marker for mechanisms that only redistribute the ambient motion, within those rooms; and the bench
   room's own modelled imprints (P2-04, P2-26), each also by its displacement integral against 5: the static imprint
   bounded over every strain direction (the largest eigenvalue of its components' Gram matrix) and beyond the solver's
   grid (Eshelby's void as Okada point sources, P2-04's analytic model, scaled to the solver's largest envelope ratio
   times 1.25), with allowances for the grid (10%) and for a random field's peak (3 / sqrt 2); the dynamic imprint
   over the whole 39 m disk, untapered, at each of four sampled incidence sides, and beyond by energy-flux conservation for
   a surface wave without attenuation to the scene's farthest corner, with P2-26's grid allowance.
"""
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri

from katabasis.runs import RESULTS, Run, load, memo
from sarsim import finite as fin
from sarsim.radiometry import snr_per_cell
from sarsim import information as inf
from sarsim.acquisition import DwellGeometry

RID = 'p2_27_requirement'
SITES = Path(__file__).resolve().parents[2] / 'sites'
ACQUISITIONS = ('giza-20250827', 'giza-20220715')
FREQS = (0.2, 2.0, 5.0, 20.0, 80.0)                   # Hz: the microseisms, 1-3 Hz, 3-8 Hz, traffic, the imprint's peak
L_EXACT = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0)   # m
L_EXACT_SLOW = (16.0,)                                # m, at 5 Hz and below, where the paired echoes stay close
SLOW = 5.0
L_LINEAR = tuple(float(v) for v in np.geomspace(0.05, 300.0, 22))
L_GRID_MAX = 4.0                                      # m: fisher_grid on the pixels up to here, the continuum beyond
TARGETS = (0.9, 0.5)
PSI_ALL = (0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4)
PSI_FAST = (0.0, np.pi / 2)                           # above 1 Hz, where the swing over psi is below 0.5%
SNR_GENIE_DB = snr_per_cell('giza-20250827')['nominal_db']   # calibrated per cell (sarsim.radiometry); was 30, untraced
SNR_ORACLE_DB = snr_per_cell('giza-20250827')['nominal_db']   # per cell, calibrated: the site's median sigma0 over ICEYE's best specified NESZ (-18 dB); was 30 dB on an untraced -26.7 dB
TARGETS_ANY = (0.05, 0.5, 0.9)                         # near chance, even odds, 95% found at 5% false alarms
EXTENT = 5.3                                          # the footprint is cut at this many L (exp(-14) of its peak)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def kl_needed(tv):
    """The smallest KL at which the ceiling min(sqrt(KL / 2), sqrt(1 - exp(-KL))) reaches tv."""
    return max(2 * tv * tv, -np.log1p(-tv * tv))


def deflection_needed(tv):
    """The deflection d at which the best test reaches tv for a weak signal: 2 Phi(d / 2) - 1 = tv."""
    return float(2 * ndtri(0.5 + tv / 2))


def area_half(L):
    return 2 * np.pi * np.log(2) * L * L


# ------------------------------------------------------------------ 1. exact, line by line

def line_kl(g, f, L, b, psi):
    """KL, exactly, for one line of the image along track (full range band: independent of its neighbours) whose pixels
    move by b exp(-x^2 / 2 L^2) cos(2 pi f (t + x / V_g) + psi) in slant range: -log det of the in-band covariance
    (unit white speckle, so the covariance without motion is I). Where the moving pixels are few beside the line's
    length, by Sylvester's identity on them alone: C = I - E_S E_S^H + B_S B_S^H."""
    k0 = 4 * np.pi / g.lam
    margin = 1.25 * f * g.V / g.Ka + 3.0                     # room for the paired echoes, which the FFT would wrap
    n = int(np.ceil((2 * EXTENT * L + 2 * margin) / g.dx)) + 64
    x = (np.arange(n) - n / 2) * g.dx
    prof = b * np.exp(-x ** 2 / (2 * L * L))
    S = np.abs(x) < EXTENT * L
    S[n // 2] = True
    in_a, _, nu = inf.band_masks(g, (n, 1))
    t = g.nu_to_time(nu)[in_a]
    fb = np.fft.fftfreq(n)[in_a]
    al = 2 * np.pi * f * x / g.V + psi
    s, nb = int(S.sum()), int(in_a.sum())
    ct, st = np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t)
    if 4 * s * s < nb * n:
        xs = np.nonzero(S)[0]
        ES = np.exp(-2j * np.pi * np.outer(fb, xs)) / np.sqrt(n)
        Phi = k0 * (np.outer(ct, (prof * np.cos(al))[S]) - np.outer(st, (prof * np.sin(al))[S]))
        U = np.hstack([ES, ES * np.exp(-1j * Phi)])
        G = U.conj().T @ U
        G[:s] *= -1
        _, ld = np.linalg.slogdet(np.eye(2 * s) + G)
        return float(-ld)
    Phi = k0 * (np.outer(ct, prof * np.cos(al)) - np.outer(st, prof * np.sin(al)))
    B = np.exp(-2j * np.pi * np.outer(fb, np.arange(n)) - 1j * Phi) / np.sqrt(n)
    Lc = np.linalg.cholesky(B @ B.conj().T)
    return float(-2 * np.sum(np.log(np.real(np.diag(Lc)))))


def columns(g, L):
    """The image's range columns across the footprint (ground spacing dr / sin theta), centred on its peak: each
    column's share of the peak, and how many columns have it."""
    dg = g.dr / np.sin(g.theta)
    j = np.arange(0, int(EXTENT * L / dg) + 1)
    return np.exp(-(j * dg) ** 2 / (2 * L * L)), np.where(j == 0, 1, 2)


def exact_requirement(g, f, L, targets):
    """The smallest D, over the phases psi, at which the summed line KL (an upper bound on the image's KL) reaches each
    target's KL: tabulated against the line's amplitude and interpolated in log-log; below the table, the table's
    smallest KL / b^2 (checked flat); the table extended upward until it brackets the answer."""
    k0 = 4 * np.pi / g.lam
    prof, mult = columns(g, L)
    need = {tv: kl_needed(tv) for tv in targets}
    psis = PSI_ALL if f <= 1.0 else PSI_FAST
    rows = []
    for psi in psis:
        b0 = 1e-3 / k0
        c2 = line_kl(g, f, L, b0, psi) / b0 ** 2
        D1 = np.sqrt(max(need.values()) / (c2 * float(np.sum(mult * prof ** 2))))
        table = {b: line_kl(g, f, L, b, psi) for b in D1 * np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.15])}

        def kl2d(D):
            bs = np.array(sorted(table))
            ks = np.array([table[b] for b in bs])
            q = D * prof
            lo = q < bs[0]
            out = np.empty_like(q)
            out[lo] = ks[0] / bs[0] ** 2 * q[lo] ** 2
            out[~lo] = np.exp(np.interp(np.log(q[~lo]), np.log(bs), np.log(ks)))
            return float(np.sum(mult * out))

        sol = {}
        for tv, kn in need.items():
            while kl2d(max(table)) < kn:                 # extend the table upward until it brackets the target
                b = 1.6 * max(table)
                table[b] = line_kl(g, f, L, b, psi)
            lo, hi = np.log(min(table) * 1e-3), np.log(max(table))
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if kl2d(np.exp(mid)) < kn else (lo, mid)
            sol[tv] = float(np.exp(hi))
        bs = sorted(table)
        rows.append({'psi': psi, 'kl_over_b2_small': c2, 'D': sol,
                     'nonlinearity_at_D': {tv: table_ratio(table, sol[tv], c2) for tv in targets},
                     'table_b_m': bs, 'table_kl': [table[b] for b in bs]})
    worst = {tv: min(r['D'][tv] for r in rows) for tv in targets}
    swing = max(r['D'][TARGETS[0]] for r in rows) / worst[TARGETS[0]] - 1
    return {'D_necessary_m': worst, 'peak_phase_rad': {tv: k0 * d for tv, d in worst.items()},
            'kl_lines_per_D2_small': float(min(r['kl_over_b2_small'] for r in rows) * np.sum(mult * prof ** 2)),
            'columns': int(np.sum(mult)), 'psi_swing': swing, 'phases': rows}


def table_ratio(table, D, c2):
    """KL / (c2 b^2) at the line amplitude D: how far the exact KL has left its quadratic start (1 = not at all)."""
    bs = np.array(sorted(table))
    ks = np.array([table[b] for b in bs])
    return float(np.exp(np.interp(np.log(D), np.log(bs), np.log(ks))) / (c2 * D * D))


# ------------------------------------------------------------------ 2. linearised; 3. helped

def gaussian(g, L, spacing):
    da, dg = spacing
    na = max(int(2 * EXTENT * L / da), 8)
    ng = max(int(2 * EXTENT * L / dg), 4)
    xa = (np.arange(na) - na / 2) * da
    xg = (np.arange(ng) - ng / 2) * dg
    return np.exp(-(xa[:, None] ** 2 + xg[None, :] ** 2) / (2 * L * L))


def linear_requirement(g, f, L, targets):
    """F per unit D^2 (the mean over psi, A, and the largest, A + |B|): exact on the image's pixels for narrow
    footprints, the continuum beyond; the linearised D for each target and the best test's weak-signal D."""
    dg = g.dr / np.sin(g.theta)
    if L <= L_GRID_MAX:
        K = gaussian(g, L, (g.dx, dg))
        A, B = inf.fisher_grid(K, g, f)
        how = 'pixels'
    else:
        step = L / 8
        K = gaussian(g, L, (step, step))
        A, B = inf.fisher_sinusoid(K, (step, step), f, g, swing=True)
        how = 'continuum'
    Fmax = A + B
    return {'fisher_per_m2': A, 'fisher_largest_per_m2': Fmax, 'swing': B / A, 'method': how,
            'D_linearised_m': {tv: float(np.sqrt(2 * kl_needed(tv) / Fmax)) for tv in targets},
            'D_best_test_weak_signal_m': {tv: float(deflection_needed(tv) / np.sqrt(A)) for tv in targets}}


def helped(g, L, targets, snr, scr):
    k0 = 4 * np.pi / g.lam
    cells = inf.cells_per_m2(g)
    return {'D_genie_m': {tv: float(np.sqrt(2 * kl_needed(tv) / (snr * cells * k0 ** 2 * np.pi * L * L))) for tv in targets},
            'D_reflector_m': {tv: float(np.sqrt(2 * kl_needed(tv) / scr) / k0) for tv in targets}}


# ------------------------------------------------------------------ 5. every pattern: the energy requirement

def any_pattern(g, p230, targets):
    """The least displacement integral E = integral max_t |du_los|^2 dA (m^4: an integrated squared displacement envelope,
    not a mechanical energy) at which each bound stops excluding each target, whatever the pattern and time history.
    - oracle: told the reflectivity, the excitation and both models, with receiver noise at SNR_ORACLE_DB per cell; the
      mean echoes differ by exp(i k u_A) - exp(i k u_0), so motion common to both worlds cancels, whatever excites it:
      E Delta^2 <= 8 sum SNR <sin^2(k0 du / 2)> <= 2 snr cells k0^2 E (P2-29's convention), valid for strong shaking;
    - certificate (Theorem C), fully developed speckle, no receiver noise: valid only where the no-cavity world's motion is
      the regional background whose floor P2-30 proves (quiet ground); with the reference taken as the identity, a model."""
    k0 = 4 * np.pi / g.lam
    density = p230['geometry']['area_density_used']
    floor = p230['background_included']['covariance_floor_proved_within_model']
    snr = 10 ** (SNR_ORACLE_DB / 10) / g.band_frac
    cells = inf.cells_per_m2(g)
    out = {'density_per_m2': density, 'covariance_floor': floor, 'oracle_snr_db': SNR_ORACLE_DB,
           'E_min_oracle_m4': {}, 'E_min_m4': {}, 'E_min_identity_m4': {}}
    for tv in targets:
        delta = float(2 * ndtri(0.5 + tv / 2))
        out['E_min_oracle_m4'][tv] = delta * delta / (2 * snr * cells * k0 * k0)
        qb = fin.required_q_for_tv(tv)
        out['E_min_identity_m4'][tv] = qb * qb / (k0 * k0 * density)
        out['E_min_m4'][tv] = out['E_min_identity_m4'][tv] * floor
    return out


def pattern_check(g):
    """The review's example and its neighbours: a Gaussian envelope (L = 4 m) alone and multiplied by cosines along track
    and across, F exactly on the pixels at each frequency; and the Gaussian's F as a share of 4 N <Phi^2>, the general
    upper bound for its envelope (not a demonstrated optimum)."""
    dg = g.dr / np.sin(g.theta)
    L = 4.0
    G = gaussian(g, L, (g.dx, dg))
    xa = (np.arange(G.shape[0]) - G.shape[0] / 2) * g.dx
    xg = (np.arange(G.shape[1]) - G.shape[1] / 2) * dg
    cap = inf.fisher_bound(G, (g.dx, dg), g)
    rows = []
    for f in FREQS:
        A0, _ = inf.fisher_grid(G, g, f)
        row = {'f_hz': f, 'L_m': L, 'gaussian_share_of_envelope_bound': A0 / cap, 'ratio_to_gaussian': {}}
        for label, K in (('cos along track, 4 m period', G * np.cos(2 * np.pi * xa / 4.0)[:, None]),
                         ('cos along track, 1 m period', G * np.cos(2 * np.pi * xa / 1.0)[:, None]),
                         ('cos across track, 4 m period', G * np.cos(2 * np.pi * xg / 4.0)[None, :])):
            A, _ = inf.fisher_grid(K, g, f)
            row['ratio_to_gaussian'][label] = A / A0
        rows.append(row)
    return rows


def static_integral(m25, g, p204):
    """integral max |u_los|^2 dA per unit strain amplitude^2 of the bench room's static imprint, bounded over every strain
    direction and beyond P2-04's grid.
    Directions: u(phi) = sum_c w_c(phi) K_c with w = (sin^2, cos^2, sin cos), |w| <= 1, so integral u_los^2 <= lambda_max(M),
    M the 3x3 Gram matrix of the three components' line-of-sight maps over the disk inside the grid's largest circle.
    Beyond it: the field of Eshelby's void of the same volume as Okada point sources in a half-space (P2-04's analytic
    model), which carries the far field's angular pattern and its r^-2 decay, scaled for each component by the largest
    ratio of the solver's envelope to the analytic one over 8 sectors and four annuli inside the grid (1.3 to 1.6 there),
    times 1.25; integrated on a polar grid to 30 km, each sector's largest value over each shell; lambda_max(M_tail) <=
    tr(M_tail) <= the sum over components. Over 16 to 40 m both fields still decay slower than r^-1 (the intermediate
    field of a source 15 m down), which is why a power law fitted inside the grid cannot bound the tail."""
    from katabasis.seismic.analytic import moment_surface_displacement, void_moment
    kern, _ = m25.static_kernels()
    xy = kern['xy']
    sp = float(np.min(np.diff(np.unique(xy[:, 0]))))
    los = np.asarray(g.los_enu)
    Lc = np.stack([kern[k] @ los for k in ('Kxx', 'Kyy', 'Kxy')])        # [3, points]
    r = np.hypot(xy[:, 0], xy[:, 1])
    th = np.arctan2(xy[:, 1], xy[:, 0])
    R = float(min(np.abs(xy[:, 0]).max(), np.abs(xy[:, 1]).max()) - sp / 2)   # the largest circle inside the grid
    inside = r <= R
    M = Lc[:, inside] @ Lc[:, inside].T * sp * sp
    lam_in = float(np.linalg.eigvalsh(M).max())
    pr = p204['manifest']['params']
    mu, lam = pr['rho'] * pr['vs'] ** 2, pr['rho'] * pr['vp'] ** 2 - 2 * pr['rho'] * pr['vs'] ** 2
    nu = pr['poisson']
    E_y = mu * (3 * lam + 2 * mu) / (lam + mu)
    a = E_y / (1 - nu * nu)
    depth, vol = -pr['chamber']['centre'][2], float(np.prod(pr['chamber']['size']))
    moments = [np.diag([a, nu * a, 0.0]), np.diag([nu * a, a, 0.0]), np.array([[0, 2 * mu, 0], [2 * mu, 0, 0], [0, 0, 0.0]])]
    an = lambda Mc, x, y: moment_surface_displacement(void_moment(Mc, vol, nu), depth, x, y, lam, mu).T @ los
    nsec = 8
    sec = np.floor((th + np.pi) / (2 * np.pi / nsec)).astype(int) % nsec
    radii = np.geomspace(R, 30000.0, 600)
    thetas = np.linspace(-np.pi, np.pi, 257)[:-1]
    RR, TT = np.meshgrid(radii, thetas, indexing='ij')
    sec_p = np.floor((TT + np.pi) / (2 * np.pi / 16)).astype(int) % 16
    shell = np.diff(np.r_[radii, radii[-1] * radii[-1] / radii[-2]])
    ratios, tail = [], 0.0
    for c, Mc in enumerate(moments):
        La = an(Mc, xy[:, 0], xy[:, 1])
        rs = []
        for lo in (8.0, 16.0, 24.0, R - 8.0):
            m = (r >= lo) & (r < lo + 8.0)
            rs += [np.abs(Lc[c][m & (sec == q)]).max() / np.abs(La[m & (sec == q)]).max() for q in range(nsec)]
        s_c = 1.25 * max(rs)
        ratios.append(max(rs))
        Lp = np.abs(an(Mc, (RR * np.cos(TT)).ravel(), (RR * np.sin(TT)).ravel())).reshape(RR.shape)
        env = np.stack([Lp[:, sec_p[0] == q].max(axis=1) for q in range(16)], axis=1)      # [radius, sector]
        integral = float(np.sum(env ** 2 * (radii * shell)[:, None] * (2 * np.pi / 16)))
        edge = float(np.max(env[-1])) * radii[-1] ** 2                          # beyond 30 km, r^-2
        integral += 2 * np.pi * edge ** 2 / (2 * radii[-1] ** 2)
        tail += s_c * s_c * integral
    return {'in_grid_lambda_max_m2': lam_in, 'tail_bound_m2': tail, 'grid_radius_m': R,
            'solver_over_analytic_largest': ratios, 'envelope_margin': 1.25, 'unit_bound_m2': lam_in + tail}


def displacement_integrals(m25, p204, p225, p226, maps, g, amb):
    """integral max_t |u_los|^2 dA (m^4) of the bench room's modelled imprints, with numerical allowances from the grid
    checks (P2-04: 10% in amplitude; P2-26: its own largest change on halving the grid).
    Static, under each ambient band at its peak strain (sqrt 2 times the band's rms, a sinusoid; a random field can peak
    higher, by 3 / sqrt 2 at three standard deviations, kept as an allowance), every strain direction, the tail beyond.
    Dynamic, under a truck over a bump 15 m away, at each of the four sampled incidence sides (0, 90, 180, 270 degrees;
    other sides are not computed): the whole disk within 39 m, untapered, and beyond it the scattered wave by energy-flux
    conservation for a surface wave without attenuation (integral |u|^2 per metre of radius constant, its value the
    largest over the outer annuli), to the scene's farthest corner; the two regions disjoint, so they add."""
    st = static_integral(m25, g, p204)
    out = []
    peak_allow = 3 / np.sqrt(2)
    for c, row in zip(p204['cases'][:3], p225['static'][:3]):
        eps0 = row['strain_amplitude']
        out.append({'label': f"bench room, static, {c['case']}", 'f_hz': row['f_hz'], 'kind': 'ambient',
                    'displacement_integral_m4': st['unit_bound_m2'] * (eps0 * 1.10 * peak_allow) ** 2,
                    'allowances': {'numerical_amplitude': 1.10, 'random_peak_amplitude': peak_allow}})
    v = amb['cultural']['bus_or_truck_over_bump']['value']
    corner = p225['manifest']['params'].get('scene_radius_m')
    if corner is None:
        src = json.loads((SITES / 'acquisitions' / 'giza-20250827.json').read_text())['source']['shape']
        corner = 0.5 * float(np.hypot(src[0] * g.dx, src[1] * g.dr / np.sin(g.theta)))
    grid = {'surface': p226['grid_check'].get('surface_small', {}).get('max_amplitude_change', 0.1),
            'surface_favourable': p226['grid_check'].get('surface_favourable_small', {}).get('max_amplitude_change', 0.11)}
    for key, label in (('surface', 'bench room'), ('surface_favourable', 'room under a 5 m roof')):
        xs, ys, fs, H = maps[f'{key}_x'], maps[f'{key}_y'], maps[f'{key}_f'], maps[f'{key}_H']
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        Rr = np.hypot(X, Y)
        dA = float(xs[1] - xs[0]) ** 2
        disk = Rr <= 39.0
        allow = 1 + grid[key]
        for f in (20.0, 80.0):
            k = int(np.argmin(np.abs(fs - f)))
            w2 = (2 * np.pi * fs[k]) ** 2
            near = (np.abs(H[:, disk, k]) ** 2).sum(axis=1) * dA / w2                                   # [side]
            sig = np.max([(np.abs(H[:, (Rr >= r0) & (Rr < r0 + 5), k]) ** 2).sum(axis=1) * dA / 5.0
                          for r0 in (25.0, 30.0, 34.5)], axis=0) / w2                                    # [side]
            far = sig * (corner - 39.0)
            tot = near + far
            j = int(np.argmax(tot))
            peak2 = 2 * v * v * allow ** 2
            out.append({'label': f'{label}, dynamic, truck over a bump 15 m away', 'f_hz': float(fs[k]), 'kind': 'cultural',
                        'displacement_integral_m4': float(tot[j]) * peak2, 'near_m4': float(near[j]) * peak2,
                        'side_deg': 90 * j, 'sides_sampled': int(len(tot)),
                        'allowances': {'numerical_amplitude': allow}})
    return out, st


# ------------------------------------------------------------------ 4. beside the curve

def markers(amb, p204, p226, maps):
    """The ambient motion's own amplitude (vertical displacement, sqrt 2 x rms / 2 pi f), the amplitude-constrained
    marker, and the bench room's modelled imprints, each with its footprint where it has one."""
    reg, cul = amb['regional'], amb['cultural']
    ambient = [
        {'label': "Giza's microseisms, regional", 'f_hz': 0.2, 'rms_velocity_m_s': reg['microseism_vertical_0.1_0.3_hz']['value'],
         'note': 'measured at Kottamya, 67 km east (M1-01); not an upper limit at the pyramids'},
        {'label': 'the noisiest stations on Earth', 'f_hz': 0.2, 'rms_velocity_m_s': p204['cases'][1]['vertical_velocity_m_s'],
         'note': "Peterson (1993) high-noise model, 0.1-0.3 Hz"},
        {'label': 'Giza 1-3 Hz, regional', 'f_hz': 2.0, 'rms_velocity_m_s': reg['vertical_1_3_hz']['value'],
         'note': 'M1-01; a quiet hill, a lower bound for Giza above 1 Hz'},
        {'label': 'Giza 3-8 Hz, regional', 'f_hz': 5.0, 'rms_velocity_m_s': reg['vertical_3_8_hz']['value'], 'note': 'M1-01'},
        {'label': 'urban background, all of 8-100 Hz at one frequency', 'f_hz': 20.0,
         'rms_velocity_m_s': cul['urban_background']['value'], 'note': 'FTA (1995), Figure 7-3'},
        {'label': 'a truck over a bump, 15 m away', 'f_hz': 20.0, 'rms_velocity_m_s': cul['bus_or_truck_over_bump']['value'],
         'note': 'FTA (1995), Figure 7-3, 50 ft'},
        {'label': 'a truck over a bump, 15 m away', 'f_hz': 80.0, 'rms_velocity_m_s': cul['bus_or_truck_over_bump']['value'],
         'note': 'FTA (1995), Figure 7-3, 50 ft'},
    ]
    for a in ambient:
        a['displacement_amplitude_m'] = np.sqrt(2) * a['rms_velocity_m_s'] / (2 * np.pi * a['f_hz'])
    cases = {k: v for k, v in p226['cases'].items() if 'peak_value' in v} if isinstance(p226['cases'], dict) \
        else {c['case']: c for c in p226['cases'] if 'peak_value' in c}
    name, top = max(((k, v['peak_value']) for k, v in cases.items()), key=lambda kv: kv[1])
    # a change of phase alone, of any size and by any mechanism (linear, nonlinear, parametric), moves the ground by at
    # most twice the ambient amplitude: |A cos(wt + p + dp) - A cos(wt + p)| <= 2A
    constrained = {'ratio': float(max(2.0, top)), 'phase_only_ratio': 2.0, 'modelled_ratio': float(top), 'case': name,
                   'f_hz': float(cases[name]['peak_f_hz']),
                   'note': "the larger of 2 (a change of the ambient motion's phase alone, of any size, by any mechanism) "
                           "and the largest room-specific motion per unit ambient motion in P2-26's rooms, applied at every "
                           "frequency (generous); covers mechanisms that change the ambient motion's phase, or scatter it "
                           "as those rooms do, not mechanisms with gain"}
    # the bench room's imprints: P2-04's static map (LOS per unit strain, rms over directions, 2 m grid)
    m = np.asarray(p204['imprint_map_los_per_strain_m'])
    sp = p204['kernel_grid']['spacing_m']
    a_static = float(np.sum(m >= m.max() / 2) * sp * sp)
    imprints = []
    for c, f in zip(p204['cases'][:3], (0.2, 0.2, 2.0)):
        imprints.append({'label': f"bench room, static, {c['case']}", 'f_hz': f, 'area_half_m2': a_static,
                         'peak_m': float(c['strain'] * np.sqrt(2) * m.max())})
    for key, label in (('surface', 'bench room'), ('surface_favourable', 'room under a 5 m roof')):
        fs, H = maps[f'{key}_f'], np.abs(maps[f'{key}_H'])
        dA = float(maps[f'{key}_x'][1] - maps[f'{key}_x'][0]) ** 2
        for f in (20.0, 80.0):
            k = int(np.argmin(np.abs(fs - f)))
            d = int(np.argmax(H[..., k].reshape(H.shape[0], -1).max(axis=1)))
            Hk = H[d, ..., k]
            v = amb['cultural']['bus_or_truck_over_bump']['value']
            imprints.append({'label': f'{label}, dynamic, truck over a bump 15 m away', 'f_hz': float(fs[k]),
                             'area_half_m2': float(np.sum(Hk >= Hk.max() / 2) * dA),
                             'peak_m': float(Hk.max() * np.sqrt(2) * v / (2 * np.pi * fs[k]))})
    return ambient, constrained, imprints


def sci(v):
    """1.2e6 as '1.2e6', for the finding; plain below ten thousand."""
    if 0.01 <= abs(v) < 1e4:
        return f"{v:.2g}" if abs(v) < 10 else f"{v:,.0f}"
    m, e = f'{v:.1e}'.split('e')
    return f"{m}e{int(e)}"


def on_curve(curve, area, key, tv):
    """D on a curve (a list of rows with area_half_m2 and D dicts) at an area, log-log interpolation."""
    xs = np.array([r['area_half_m2'] for r in curve])
    ys = np.array([r[key][tv] for r in curve])
    return float(np.exp(np.interp(np.log(area), np.log(xs), np.log(ys))))


def main():
    amb = json.loads((SITES / 'ambient.json').read_text())
    p204, p226, p223 = (load(r) for r in ('p2_04_chamber_imprint', 'p2_26_imprint_spectrum', 'p2_23_every_reader'))
    p225, p230 = load('p2_25_information_bound'), load('p2_30_finite_certificates')
    m25 = module('p225', Path(__file__).with_name('p2_25_information_bound.py'))
    maps = np.load(RESULTS / 'p2_26_imprint_spectrum' / 'maps.npz')
    snr = 10 ** (SNR_GENIE_DB / 10)
    scr = 10 ** (p223['manifest']['params']['reflector_scr_db'] / 10)
    params = {'acquisitions': ACQUISITIONS, 'frequencies_hz': FREQS, 'targets_tv': TARGETS,
              'footprint': 'D exp(-r^2 / 2 L^2), area within half the peak', 'widths_exact_m': L_EXACT,
              'widths_exact_slow_m': L_EXACT_SLOW, 'widths_linearised_m': L_LINEAR, 'genie_snr_db': SNR_GENIE_DB,
              'reflector_scr_db': 10 * np.log10(scr), 'scattering': 'fully developed speckle, no reference image',
              'noise_free': True, 'motion': 'coherent for the whole pass, worst phase'}
    with Run(RID, 'What one image would need: the motion a room must add to the ground', params) as run:
        out = {}
        for acq in ACQUISITIONS:
            g = DwellGeometry.from_record(acq)
            per_f = []
            for f in FREQS:
                t0 = time.time()
                widths = L_EXACT + (L_EXACT_SLOW if f <= SLOW else ())
                exact = []
                for L in widths:
                    ex = memo(RID, f'exact-{acq}-{f}-{L}', lambda: exact_requirement(g, f, L, TARGETS), __file__,
                              version='exact-v1')
                    lin = linear_requirement(g, f, L, TARGETS)
                    exact.append({'L_m': L, 'area_half_m2': area_half(L), **ex,
                                  'lines_over_image_small': ex['kl_lines_per_D2_small'] / (lin['fisher_largest_per_m2'] / 2),
                                  'necessary_over_linearised': {tv: ex['D_necessary_m'][tv] / lin['D_linearised_m'][tv]
                                                                for tv in TARGETS}})
                    print(f"  {acq} f={f} L={L}: necessary {ex['D_necessary_m'][0.9]:.3g} m "
                          f"(phase {ex['peak_phase_rad'][0.9]:.2f} rad), linearised {lin['D_linearised_m'][0.9]:.3g}, "
                          f"best test {lin['D_best_test_weak_signal_m'][0.9]:.3g}", flush=True)
                linear = [{'L_m': L, 'area_half_m2': area_half(L), **linear_requirement(g, f, L, TARGETS),
                           **helped(g, L, TARGETS, snr, scr)} for L in L_LINEAR]
                per_f.append({'f_hz': f, 'exact': exact, 'linear': linear, 'runtime_s': round(time.time() - t0, 1)})
            out[acq] = {'look_angle_deg': float(np.degrees(g.theta)), 'cells_per_m2': inf.cells_per_m2(g),
                        'aperture_s': g.aperture_time, 'stretch_s': g.R0 / g.V_platform, 'frequencies': per_f}

        ambient, constrained, imprints = markers(amb, p204, p226, maps)
        main_acq = out[ACQUISITIONS[0]]['frequencies']
        byf = {r['f_hz']: r for r in main_acq}
        # how far each marker falls short of the necessary curve at its own footprint and frequency
        for im in imprints:
            f = min(FREQS, key=lambda q: abs(q - im['f_hz']))
            need = on_curve(byf[f]['exact'], im['area_half_m2'], 'D_necessary_m', 0.9)
            im['necessary_m'] = need
            im['short_by'] = need / im['peak_m']
        for a in ambient:
            ex = byf[a['f_hz']]['exact']
            best = min(ex, key=lambda r: r['D_necessary_m'][0.9])
            a['least_necessary_m'] = best['D_necessary_m'][0.9]
            a['least_necessary_area_m2'] = best['area_half_m2']
            a['constrained_short_by'] = best['D_necessary_m'][0.9] / (constrained['ratio'] * a['displacement_amplitude_m'])
        checks = {
            'lines_over_image': {acq: [float(r['lines_over_image_small']) for fr in out[acq]['frequencies']
                                       for r in fr['exact']] for acq in out},
            'necessary_over_linearised': [float(r['necessary_over_linearised'][0.9]) for fr in main_acq for r in fr['exact']],
            'largest_psi_swing': max(float(r['psi_swing']) for fr in main_acq for r in fr['exact']),
            'largest_nonlinearity': max(float(ph['nonlinearity_at_D'][tv]) for acq in out for fr in out[acq]['frequencies']
                                        for r in fr['exact'] for ph in r['phases'] for tv in TARGETS),
            'grid_vs_continuum': []}
        g = DwellGeometry.from_record(ACQUISITIONS[0])
        dg = g.dr / np.sin(g.theta)
        for f in FREQS:
            for L in (1.0, 2.0, 4.0):
                A, _ = inf.fisher_grid(gaussian(g, L, (g.dx, dg)), g, f)
                Ac = inf.fisher_sinusoid(gaussian(g, L, (L / 8, L / 8)), (L / 8, L / 8), f, g)
                checks['grid_vs_continuum'].append({'f_hz': f, 'L_m': L, 'continuum_over_pixels': Ac / A})

        # 5. every pattern: the displacement-integral requirement, the imprints against it, what each bound implies at
        # every operating point (found <= false alarms + TV), and the area a phase-only mechanism would have to change
        anyp = any_pattern(g, p230, TARGETS_ANY)
        E_or, E_ct = anyp['E_min_oracle_m4'], anyp['E_min_m4']
        integrals, stat = displacement_integrals(m25, p204, p225, p226, maps, g, amb)
        k0 = 4 * np.pi / g.lam
        snr_or = 10 ** (SNR_ORACLE_DB / 10) / g.band_frac
        tv_of = lambda E: float(2 * ndtr(np.sqrt(2 * snr_or * inf.cells_per_m2(g) * k0 * k0 * E) / 2) - 1)
        outcome = lambda tv: 'near chance' if tv < 0.05 else ('95% at 5% excluded' if tv < 0.9 else 'unresolved')
        for im in imprints:
            e = min((x for x in integrals if x['label'] == im['label']), key=lambda x: abs(x['f_hz'] - im['f_hz']))
            E = e['displacement_integral_m4']
            im.update({k: e[k] for k in ('displacement_integral_m4', 'allowances', 'kind') if k in e})
            im['oracle_tv_upper'] = tv_of(E)
            im['outcome'] = outcome(im['oracle_tv_upper'])
            im['short_by_oracle'] = float(np.sqrt(E_or[0.9] / E))
            if e['kind'] == 'ambient':          # the certificate holds where the no-cavity world moves with the background
                im['short_by_certificate'] = float(np.sqrt(E_ct[0.9] / E))
            else:
                im['near_m4'] = e['near_m4']
                im['side_deg'] = e['side_deg']
                im['oracle_tv_upper_near_only'] = tv_of(e['near_m4'])
                im['short_by_oracle_near_only'] = float(np.sqrt(E_or[0.9] / e['near_m4']))
        hv = p204['hv']
        for a in ambient:
            a_los = np.sqrt(1 + hv * hv) * a['displacement_amplitude_m']            # a Rayleigh wave's LOS amplitude, at most
            a['phase_only_area_oracle_m2'] = float(E_or[0.9] / (2 * a_los) ** 2)
            if a['f_hz'] <= 5.0:
                a['phase_only_area_certificate_m2'] = float(E_ct[0.9] / (2 * a_los) ** 2)
        patterns = pattern_check(g)
        for fr in main_acq:
            for r in fr['linear']:
                K = gaussian(g, r['L_m'], (r['L_m'] / 8, r['L_m'] / 8)) if r['method'] == 'continuum' else gaussian(g, r['L_m'], (g.dx, dg))
                sp_ = (r['L_m'] / 8, r['L_m'] / 8) if r['method'] == 'continuum' else (g.dx, dg)
                r['share_of_envelope_bound'] = float(r['fisher_largest_per_m2'] / inf.fisher_bound(K, sp_, g))
        sb = [im for im in imprints if 'static' in im['label']][0]
        dyn = [im for im in imprints if 'dynamic' in im['label']]
        e02 = byf[0.2]['exact']
        e80 = byf[80.0]['exact']
        at = lambda rows, L: next(r for r in rows if r['L_m'] == L)
        p02 = next(r for r in patterns if r['f_hz'] == 0.2)
        rt = p02['ratio_to_gaussian']
        amb_reg = next(a for a in ambient if 'microseisms, regional' in a['label'])
        amb_t80 = next(a for a in ambient if 'truck' in a['label'] and a['f_hz'] == 80.0)
        unres = [im for im in dyn if im['outcome'] == 'unresolved']
        finding = (
            f"Under the stated conditions (the Giza dwell, fully developed speckle, no reference image, motion coherent for "
            f"the whole pass at its worst phase), a detection at 95% with 5% false alarms needs, for one smooth (Gaussian) "
            f"pattern of room-specific line-of-sight motion, at least {at(e02, 8.0)['D_necessary_m'][0.9] * 1e3:.2f} mm at "
            f"0.2 Hz and {at(e80, 8.0)['D_necessary_m'][0.9] * 1e6:.1f} um at 80 Hz over {area_half(8.0):.0f} m^2, exactly "
            f"within the model. That binds only that family: a Gaussian 4 m wide at 0.2 Hz times a cosine of 4 m period along "
            f"track holds {rt['cos along track, 4 m period']:.0f} times its information ({rt['cos along track, 1 m period']:.0f} "
            f"times for 1 m) though nowhere larger, and the Gaussian's is {sci(p02['gaussian_share_of_envelope_bound'])} of the "
            f"general upper bound for its envelope. For every pattern and time history the requirement is on the displacement "
            f"integral E, the integral of max_t |du_los|^2 over the ground (m^4, an integrated squared displacement envelope, "
            f"not a mechanical energy): an oracle told the reflectivity, the excitation and both models, with receiver noise at "
            f"{SNR_ORACLE_DB:.0f} dB per cell, cannot reach 95% at 5% below {sci(E_or[0.9])} m^4 (a spatial rms of temporal peak "
            f"displacement of {np.sqrt(E_or[0.9] / 1000.0) * 1e6:.1f} um over 1,000 m^2), whatever the motion common to both "
            f"worlds, since that cancels from its mean echoes; under quiet ground the certificate, with Giza's background in "
            f"the reference, gives {sci(E_ct[0.9])} m^4 without receiver noise. At every operating point found <= false alarms "
            f"+ TV: below {sci(E_or[0.05])} m^4 the oracle stays within 5 points of chance. The bench room's static imprint "
            f"under the regional microseisms, bounded over every strain direction and beyond the solver's grid (the tail "
            f"by the analytic point-source field scaled {max(stat['solver_over_analytic_largest']) * 1.25:.1f} times) with numerical and peak allowances, falls short of the oracle's requirement "
            f"by {sci(sb['short_by_oracle'])} in amplitude and of the certificate's by {sci(sb['short_by_certificate'])}: "
            f"{sb['outcome']} (oracle TV at most {sb['oracle_tv_upper']:.1e}). Under a truck over a bump 15 m away all pass, "
            f"at the worst of four sampled sides, with the scattered wave carried unattenuated to the scene's corner, the "
            f"dynamic imprints fall short by {sci(min(im['short_by_oracle'] for im in dyn))} to "
            f"{sci(max(im['short_by_oracle'] for im in dyn))}"
            + (f" ({', '.join(im['label'].split(',')[0] + ' at ' + str(round(im['f_hz'])) + ' Hz' for im in unres)} unresolved)" if unres else "")
            + f", and counted within 39 m only by {sci(min(im['short_by_oracle_near_only'] for im in dyn))} to "
            f"{sci(max(im['short_by_oracle_near_only'] for im in dyn))}: under strong shaking the answer rests on how far the "
            f"scattered wave carries. Crossing a threshold demonstrates neither detectability nor usefulness. A mechanism that "
            f"only changes the ambient motion's phase moves the ground by at most twice its amplitude; to escape this exclusion "
            f"under its stated amplitude and noise it must change the phase over more than "
            f"{sci(amb_reg['phase_only_area_oracle_m2'])} m^2 under the regional microseisms (a disk of radius "
            f"{np.sqrt(amb_reg['phase_only_area_oracle_m2'] / np.pi):.0f} m) or {sci(amb_t80['phase_only_area_oracle_m2'])} m^2 "
            f"under the truck at 80 Hz. These are requirements conditional on the model; they say nothing of routes outside it.")
        run.save({'acquisitions': out, 'ambient': ambient, 'constrained': constrained, 'imprints': imprints,
                  'any_pattern': anyp, 'pattern_check': patterns, 'static_integral': stat,
                  'checks': checks, 'kl_needed': {tv: kl_needed(tv) for tv in TARGETS},
                  'deflection_needed': {tv: deflection_needed(tv) for tv in TARGETS}, 'finding': finding})
        print(finding)


if __name__ == '__main__':
    main()
