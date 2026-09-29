"""Finite, algorithm-independent statistical certificates for a declared SAR model (P2-30).

A complement to sarsim.information, from an independent derivation (29 September 2026; its proofs are set out on the
site's proof page and in PROOF.md). sarsim.information computes the Fisher information one image holds and bounds the
finite change with an explicit remainder, tightly, taking the reference covariance as the identity; this module bounds
the finite change more loosely but more simply, straight from a phase-energy envelope, and can include the common
background motion through a proved covariance floor (Theorem C), and bounds the oracle told everything but the answer
(Theorem D, as P2-29 computes it pulse by pulse), and two-world depth risk (Theorem E).

No local Fisher approximation is used. Complex covariance convention:
E[(Y-m)(Y-m).conj().T] = C. Natural logarithms.
This module proves no upper bound on real-world cavity motion or instrument noise.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import ndtr, erfinv


def stable_exp_difference(phase):
    """exp(i*phase)-1, preserving phase changes far below one radian."""
    return np.expm1(1j * np.asarray(phase, float))


def stable_range_difference(satellite_minus_ground, ground_displacement):
    """|b-u|-|b| without subtracting two ~600 km ranges.

    b and u broadcast with last axis length 3; positive means receding.
    Rationalizing the difference preserves tiny differential displacement.
    """
    b, u = np.broadcast_arrays(np.asarray(satellite_minus_ground, float),
                               np.asarray(ground_displacement, float))
    numerator = -2 * np.sum(b * u, axis=-1) + np.sum(u * u, axis=-1)
    denominator = np.linalg.norm(b - u, axis=-1) + np.linalg.norm(b, axis=-1)
    if np.any(denominator == 0):
        raise ValueError("A nonzero sensor/target separation is required")
    return numerator / denominator


def perturbation_certificate(q, mean_energy=0.0):
    """Finite Gaussian covariance certificate from q >= ||C0^-1/2 delta A||F.

    A includes latent-Gaussian signal columns and common noise columns, such
    that C0=A0 A0*. delta A has zero entries in unchanged noise columns.
    Mean energy is ||C0^-1/2 (m1-m0)||^2. Returned TV is an upper bound,
    not the optimal detector's performance. q must already be a valid bound.
    """
    if not math.isfinite(q) or not math.isfinite(mean_energy) or q < 0 or mean_energy < 0:
        raise ValueError("q and mean_energy must be nonnegative and finite")
    rho = 2 * q + q * q
    if rho >= 1:
        return {"q": q, "rho": rho, "kl_upper": None, "tv_upper": 1.0,
                "informative": False}
    kl = rho * rho / (2 * (1 - rho) ** 2) + mean_energy
    tv = min(1.0, math.sqrt(kl / 2))
    return {"q": q, "rho": rho, "kl_upper": kl, "tv_upper": tv,
            "informative": tv < 1}


def certificate_from_energy(phase_energy, covariance_floor=1.0,
                            texture_cap=1.0, mean_energy=0.0):
    """q^2 <= texture_cap * phase_energy / covariance_floor.

    phase_energy bounds sum_{k,j}|T0[k,j]|^2 Phi[k,j]^2. covariance_floor
    is a lower eigenvalue bound for C0; texture_cap is an upper eigenvalue
    bound on the zero-mean Gaussian scatterer covariance in matching units.
    """
    if phase_energy < 0 or covariance_floor <= 0 or texture_cap < 0:
        raise ValueError("Invalid energy, covariance floor, or texture cap")
    q = math.sqrt(texture_cap * phase_energy / covariance_floor)
    return perturbation_certificate(q, mean_energy)


def _g_scalar(e):
    """e-log(1+e), avoiding cancellation when |e| is small."""
    if e <= -1:
        raise ValueError("The changed covariance must be positive definite")
    if abs(e) >= 1e-3:
        return float(e - np.log1p(e))
    # Alternating analytic power series, not a local KL approximation.
    total = 0.0
    power = e * e
    for k in range(2, 80):
        term = (1 if k % 2 == 0 else -1) * power / k
        total += term
        if abs(term) <= abs(total) * 2e-17:
            break
        power *= e
    return max(0.0, float(total))


def whiten_delta(c0, delta_c):
    """L^-1 delta C L^-*, for C0=LL*, using solves instead of an inverse."""
    L = np.linalg.cholesky(c0)
    left = solve_triangular(L, delta_c, lower=True)
    E = solve_triangular(L, left.conj().T, lower=True).conj().T
    return (E + E.conj().T) / 2, L


def finite_complex_kl(c0, delta_c, delta_mean=None):
    """Exact finite proper-complex Gaussian KL, supplied a stable delta C.

    Mathematical formula is exact. This float implementation is numerical;
    it is not a formal interval-arithmetic certificate.
    """
    E, L = whiten_delta(np.asarray(c0, complex), np.asarray(delta_c, complex))
    e = np.linalg.eigvalsh(E)
    covariance_term = sum(_g_scalar(float(x)) for x in e)
    mean_term = 0.0
    if delta_mean is not None:
        dm = solve_triangular(L, np.asarray(delta_mean, complex), lower=True)
        mean_term = float(np.vdot(dm, dm).real)
    rho = float(np.max(np.abs(e), initial=0))
    frob2 = float(np.vdot(E, E).real)
    upper = frob2 / (2 * (1 - rho) ** 2) + mean_term if rho < 1 else None
    return {"kl": covariance_term + mean_term, "mean_energy": mean_term,
            "rho_actual": rho, "frob2_actual": frob2,
            "finite_norm_kl_upper": upper}


def stable_covariance_difference(a0, delta_a):
    """C1-C0 from delta A, retaining tiny perturbations."""
    return (delta_a @ a0.conj().T + a0 @ delta_a.conj().T
            + delta_a @ delta_a.conj().T)


def oracle_tv_from_mean_energy(energy):
    """Finite oracle TV upper bound erf(sqrt(E d)/2), common Gaussian noise.

    If d is constant (e.g. fixed known scene), the conditional TV is exact.
    When latent scene/source variables vary, this applies Jensen to their
    expected d and is generally an upper bound on joint-oracle TV.
    """
    if energy < 0 or not math.isfinite(energy):
        raise ValueError("Energy must be finite and nonnegative")
    return math.erf(math.sqrt(energy) / 2)


def oracle_mean_energy(delta_t, second_moment, noise_covariance):
    """E d = tr(N^-1 delta T E[s s*] delta T*); no Gaussian s needed."""
    delta_t = np.asarray(delta_t, complex)
    M = np.asarray(second_moment, complex)
    N = np.asarray(noise_covariance, complex)
    if M.shape != (delta_t.shape[1],) * 2:
        raise ValueError("Scene second moment does not match transfer columns")
    np.linalg.cholesky(N)  # require strictly positive Gaussian noise
    if np.linalg.eigvalsh((M + M.conj().T) / 2).min() < -1e-12:
        raise ValueError("Scene second moment must be positive semidefinite")
    return max(0.0, float(np.trace(np.linalg.solve(
        N, delta_t @ M @ delta_t.conj().T)).real))


def conditional_gaussian_roc(energy, fpr=0.05):
    """Exact optimal ROC for two fixed means with common proper-CN noise.

    This cannot be substituted for a variable-latent-scene ROC using E d.
    """
    from scipy.special import ndtri
    if not 0 < fpr < 1 or energy < 0:
        raise ValueError("Invalid ROC arguments")
    return float(ndtr(ndtri(fpr) + math.sqrt(2 * energy)))


def peak_phase_energy(in_band_bins, peak_phase, affected_fraction=1.0):
    """Q <= m*(J/n)*peak_phase^2 for a unit-magnitude normalized DFT model."""
    if in_band_bins < 0 or peak_phase < 0 or not 0 <= affected_fraction <= 1:
        raise ValueError("Invalid phase-envelope arguments")
    return float(in_band_bins) * affected_fraction * peak_phase * peak_phase


def baseline_floor_from_phase_energy(background_energy):
    """Lower bound (1-sqrt(Qbg))^2 when a DFT baseline acquires bounded phase.

    If sqrt(Qbg)>=1 the singular-value argument supplies no positive floor.
    This does not imply the actual covariance is singular.
    """
    if background_energy < 0 or not math.isfinite(background_energy):
        raise ValueError("Invalid background phase energy")
    r=math.sqrt(background_energy)
    return (1-r)**2 if r < 1 else None


def required_q_for_tv(tv_target):
    """Where this sufficient exclusion ceases to certify TV < tv_target.

    In the mean-zero covariance theorem: rho/(2(1-rho))=tv_target.
    A larger q does not establish detectability.
    """
    if not 0 < tv_target < 1:
        raise ValueError("Target must be in (0,1)")
    rho = 2 * tv_target / (1 + 2 * tv_target)
    return math.sqrt(1 + rho) - 1


def required_oracle_energy(tv_target):
    if not 0 < tv_target < 1:
        raise ValueError("Target must be in (0,1)")
    return 4 * float(erfinv(tv_target)) ** 2


def depth_absolute_risk_lower(separation, tv_upper):
    """Equal-prior mean absolute depth error >= Delta/2*(1-TV upper)."""
    if separation < 0 or not 0 <= tv_upper <= 1:
        raise ValueError("Invalid depth-risk arguments")
    return separation / 2 * (1 - tv_upper)


def pair_tv_via_common_reference(tv0, tv1):
    return min(1.0, max(0.0, tv0) + max(0.0, tv1))
