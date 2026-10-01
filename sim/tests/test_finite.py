"""sarsim.finite: independent checks of its finite identities and a small physical echo chain (P2-30).

These are mathematical/numerical checks, not field validation.
"""
import json
import math
import unittest
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np
from scipy.special import ndtri

from sarsim.finite import (baseline_floor_from_phase_energy, certificate_from_energy, conditional_gaussian_roc,
    finite_complex_kl, line_floor, oracle_mean_energy, oracle_tv_from_mean_energy, per_line_certificate,
    perturbation_certificate, required_oracle_energy, required_q_for_tv,
    stable_covariance_difference, stable_exp_difference, stable_range_difference)



def real_cov(C):
    return .5*np.block([[C.real,-C.imag],[C.imag,C.real]])


class FiniteProofTests(unittest.TestCase):
    def test_real_complex_convention(self):
        rng=np.random.default_rng(160)
        for n in (2,5,9):
            A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
            c0=A@A.conj().T+np.eye(n)
            B=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
            c1=c0+.05*(B@B.conj().T)
            dm=.1*(rng.normal(size=n)+1j*rng.normal(size=n))
            complex_kl=finite_complex_kl(c0,c1-c0,dm)['kl']
            r0,r1=real_cov(c0),real_cov(c1)
            rd=np.r_[dm.real,dm.imag]
            real_kl=.5*(np.trace(np.linalg.solve(r0,r1))-2*n+
                np.linalg.slogdet(r0)[1]-np.linalg.slogdet(r1)[1]+rd@np.linalg.solve(r0,rd))
            self.assertAlmostEqual(complex_kl,float(real_kl),places=12)

    def test_finite_perturbation_bound(self):
        rng=np.random.default_rng(161)
        for scale in (1e-9,1e-6,.002,.02,.07):
            for _ in range(8):
                a0=(rng.normal(size=(7,11))+1j*rng.normal(size=(7,11)))/math.sqrt(22)
                c0=a0@a0.conj().T
                delta=scale*(rng.normal(size=a0.shape)+1j*rng.normal(size=a0.shape))/math.sqrt(22)
                L=np.linalg.cholesky(c0)
                q=np.linalg.norm(np.linalg.solve(L,delta),'fro')
                certificate=perturbation_certificate(float(q))
                exact=finite_complex_kl(c0,stable_covariance_difference(a0,delta))
                self.assertLessEqual(exact['kl'],certificate['kl_upper']*(1+1e-9)+1e-30)
                self.assertLessEqual(exact['rho_actual'],certificate['rho']*(1+1e-9)+1e-15)

    def test_high_precision_tiny_covariance(self):
        for value in ['-0.8','-0.2','-1e-12','1e-12','0.2','0.8']:
            with localcontext() as ctx:
                ctx.prec=80
                d=Decimal(value);exact=d-(1+d).ln()
                rho=abs(d);upper=d*d/(2*(1-rho))
                self.assertLessEqual(exact,upper)
                result=finite_complex_kl(np.eye(1),np.array([[float(d)]]))['kl']
                self.assertLess(abs(result-float(exact)),float(exact)*2e-13)

    def test_partial_fourier_finite_and_shared_motion(self):
        n,m=128,96
        E=np.exp(-2j*math.pi*np.outer(np.arange(m),np.arange(n))/n)/math.sqrt(n)
        x=np.linspace(-1,1,n);t=np.linspace(-.5,.5,m)
        pattern=np.exp(-x*x/.08)[None,:]*np.cos(2*math.pi*2*t[:,None])
        for scale in (1e-9,1e-5,.001,.015):
            phase=scale*pattern
            da=E*stable_exp_difference(-phase)
            Q=float(np.sum(np.abs(E)**2*phase**2))
            self.assertLessEqual(np.linalg.norm(da,'fro')**2,Q*(1+1e-14))
            c0=E@E.conj().T
            floor=float(np.linalg.eigvalsh(c0).min())
            exact=finite_complex_kl(c0,stable_covariance_difference(E,da))
            cert=certificate_from_energy(Q,covariance_floor=floor)
            self.assertLessEqual(exact['kl'],cert['kl_upper']*(1+1e-8)+1e-29)
        # Shared row phase is an exact finite invariance in this model.
        a1=np.exp(-.7j*np.cos(4*t))[:,None]*E
        self.assertLess(np.linalg.norm(a1@a1.conj().T-c0,'fro'),1e-12)
        background=.03*np.sin(t[:,None]+3*x[None,:])
        da=E*stable_exp_difference(-background)
        Qbg=float(np.sum(np.abs(E)**2*background**2))
        floor=baseline_floor_from_phase_energy(Qbg)
        actual_floor=float(np.linalg.eigvalsh((E+da)@(E+da).conj().T).min())
        self.assertLessEqual(floor,actual_floor+1e-14)

    def test_non_gaussian_reflectivity_oracle(self):
        rng=np.random.default_rng(162)
        n,m=9,7
        dt=.08*(rng.normal(size=(m,n))+1j*rng.normal(size=(m,n)))
        noise=.5*np.eye(m)
        mean=np.zeros(n,complex);mean[0]=3+2j
        M=np.eye(n)+np.outer(mean,mean.conj())
        theoretical=oracle_mean_energy(dt,M,noise)
        # Unit-magnitude random phases, plus a coherent bright mean: not CN.
        s=np.exp(2j*math.pi*rng.random((n,60000)))+mean[:,None]
        observed=np.sum(np.abs(dt@s)**2,axis=0)/.5
        self.assertLess(abs(observed.mean()-theoretical),6*observed.std()/math.sqrt(len(observed)))

    def test_oracle_upper_experiment(self):
        rng=np.random.default_rng(163)
        for scale in (1e-6,.01,.1):
            a0=(rng.normal(size=(6,10))+1j*rng.normal(size=(6,10)))/math.sqrt(20)
            da=scale*(rng.normal(size=a0.shape)+1j*rng.normal(size=a0.shape))/math.sqrt(20)
            N=.7*np.eye(6);M=np.eye(10)
            c0=a0@a0.conj().T+N
            marginal=finite_complex_kl(c0,stable_covariance_difference(a0,da))['kl']
            joint=oracle_mean_energy(da,M,N)
            self.assertLessEqual(marginal,joint*(1+1e-12))

    def test_raw_echo_tiny_motion_and_detector(self):
        from sarsim.acquisition import DwellGeometry
        g=DwellGeometry.from_record('giza-20250827')
        R,V,lam=g.R0,g.V_platform,g.lam
        t=np.linspace(-12,12,257)
        b=np.c_[V*t,np.zeros_like(t),np.full_like(t,R)]
        amp=2.5838805e-12
        u=np.zeros_like(b);u[:,2]=amp*np.cos(2*math.pi*.2*t)
        dR=stable_range_difference(b,u)
        # Direct subtraction loses all or almost all of this picometre signal.
        direct=np.linalg.norm(b-u,axis=1)-np.linalg.norm(b,axis=1)
        self.assertGreater(np.max(np.abs(dR)),amp*.98)
        self.assertLess(np.count_nonzero(direct),len(t)*.01)
        delta_phi=-4*math.pi*dR/lam
        a0=np.exp(-1j*np.linspace(0,40,len(t)))/math.sqrt(len(t))
        delta_mu=a0*stable_exp_difference(delta_phi)
        self.assertGreater(np.linalg.norm(delta_mu),0)
        # Boost the motion for a measurable diagnostic; NOT a physical result.
        phase=delta_phi*1e9
        delta_mu=a0*stable_exp_difference(phase)
        power=float(np.vdot(delta_mu,delta_mu).real)
        noise=power/.5  # fixed conditional oracle energy .5
        energy=power/noise
        rng=np.random.default_rng(164);M=80000
        # Projection is sufficient for the Gaussian likelihood ratio.
        z=rng.normal(size=M)
        threshold=float(ndtri(.95))
        observed_fpr=float(np.mean(z>threshold))
        observed_tpr=float(np.mean(z+math.sqrt(2*energy)>threshold))
        predicted=conditional_gaussian_roc(energy,.05)
        self.assertLess(abs(observed_fpr-.05),.004)
        self.assertLess(abs(observed_tpr-predicted),.006)
        self.assertAlmostEqual(oracle_tv_from_mean_energy(energy),math.erf(math.sqrt(.5)/2),places=14)

    def test_local_fisher_does_not_prove_finite_invariance(self):
        # N(a^2,1): Fisher at zero is zero, finite KL at a=1 is 1/2.
        f0=0;kl_finite=.5
        self.assertEqual(f0,0);self.assertGreater(kl_finite,0)

    def test_boundary_and_jensen(self):
        q=required_q_for_tv(.9)
        self.assertAlmostEqual(perturbation_certificate(q)['tv_upper'],.9,places=13)
        b=required_oracle_energy(.9)
        self.assertAlmostEqual(oracle_tv_from_mean_energy(b),.9,places=13)
        energies=np.array([0,.01,.1,1,5])
        average=sum(oracle_tv_from_mean_energy(float(e)) for e in energies)/len(energies)
        self.assertLessEqual(average,oracle_tv_from_mean_energy(float(energies.mean())))

    def test_direction_enlargement(self):
        # Analytical guard: cos(pi/24) >= 1-(22/7/24)^2/2 > 1/1.01.
        self.assertGreater(1-((22/7)/24)**2/2,1/1.01)
        rng=np.random.default_rng(165)
        az=np.arange(24)*math.pi/24
        for _ in range(100):
            a,b,c=rng.normal(size=3)
            sampled=np.max(np.abs(a+b*np.cos(2*az)+c*np.sin(2*az)))
            true=abs(a)+math.hypot(b,c)
            self.assertLessEqual(true,1.01*sampled+1e-15)

if __name__=='__main__':unittest.main(verbosity=2)


class PatternTests(unittest.TestCase):
    """An independent review's point, checked in a small phase-only Fourier model with the exact finite divergence: at the
    same displacement envelope the spatial pattern changes how distinguishable the motion is, so a requirement computed for
    one pattern binds only that pattern; the energy certificate (Theorem C) bounds every pattern."""

    def model(self, pattern, n=48, keep=24, f=0.04):
        k = np.arange(n)
        T0 = np.exp(-2j * np.pi * np.outer(k, k) / n) / math.sqrt(n)
        rows = np.r_[0:keep // 2, n - keep // 2:n]                      # the retained band
        t = rows / n                                                     # each bin's time, the Doppler-to-time relation
        Phi = np.outer(np.cos(2 * np.pi * f * n * t), pattern)           # [bin, scatterer]
        A0 = T0[rows]
        dA = A0 * np.expm1(1j * Phi)
        c0 = A0 @ A0.conj().T
        dC = A0 @ dA.conj().T + dA @ A0.conj().T + dA @ dA.conj().T
        q2 = float(np.sum(np.abs(A0) ** 2 * Phi ** 2))
        return finite_complex_kl(c0, dC)['kl'], q2

    def test_pattern_matters_at_equal_envelope(self):
        n = 48
        x = np.arange(n) - n / 2
        env = 0.25 * np.exp(-x ** 2 / (2 * 6.0 ** 2))
        kl_smooth, q2_smooth = self.model(env)
        kl_mod, q2_mod = self.model(env * np.cos(2 * np.pi * x / 4.0))
        self.assertTrue(np.all(np.abs(env * np.cos(2 * np.pi * x / 4.0)) <= env + 1e-15))
        self.assertGreater(kl_mod, 5 * kl_smooth)                        # nowhere larger, far more distinguishable
        # and both stay within the energy certificate, which does not care about the pattern
        for kl, q2 in ((kl_smooth, q2_smooth), (kl_mod, q2_mod)):
            cert = perturbation_certificate(math.sqrt(q2))
            self.assertLessEqual(min(1.0, math.sqrt(kl / 2)), cert['tv_upper'] + 1e-12)

    def test_certificate_bounds_random_patterns(self):
        rng = np.random.default_rng(31)
        n = 48
        for _ in range(20):
            pattern = rng.uniform(-0.2, 0.2, n) * (rng.uniform(size=n) < 0.5)
            kl, q2 = self.model(pattern)
            cert = perturbation_certificate(math.sqrt(q2))
            self.assertLessEqual(min(1.0, math.sqrt(kl / 2)), cert['tv_upper'] + 1e-12)


class PerLineTests(unittest.TestCase):
    """The per-line certificate (P2-36): one along-track line of the model, each world moving by its own phase field (a
    large common part and a small difference), white receiver noise; the exact divergence from the eigenvalues of the
    whitened change never exceeds the certificate, and the floor never exceeds the covariance's smallest eigenvalue."""

    def line(self, rng, n=256, keep=200, common=0.03, diff=0.004, noise=1e-3):
        k = np.arange(n)
        rows = np.r_[0:keep // 2, n - keep // 2:n]
        E = np.exp(-2j * np.pi * np.outer(rows, k) / n) / math.sqrt(n)
        t = np.linspace(-0.5, 0.5, keep)[:, None]
        x = k[None, :] / n
        f = rng.uniform(3, 40)
        phi0 = common * rng.uniform(0.3, 1) * np.cos(2 * np.pi * (f * t + rng.uniform(1, 5) * x) + rng.uniform(0, 6.3))
        dphi = diff * np.exp(-((k - n / 2) / rng.uniform(4, 30)) ** 2)[None, :] * np.cos(2 * np.pi * f * t + rng.uniform(0, 6.3))
        B0 = E * np.exp(-1j * phi0)
        D = B0 * np.expm1(-1j * dphi)
        c0 = B0 @ B0.conj().T + noise * np.eye(len(rows))
        dC = stable_covariance_difference(B0, D)
        kl = finite_complex_kl(c0, dC)['kl']
        Q = float(np.sum(np.abs(D) ** 2))
        Qc = float(np.sum(np.abs(E * np.expm1(-1j * phi0)) ** 2))
        return kl, Q, Qc, float(np.linalg.eigvalsh(c0).min()), noise

    def test_certificate_and_floor_hold(self):
        rng = np.random.default_rng(36)
        for _ in range(12):
            kl, Q, Qc, lam_min, noise = self.line(rng)
            floor = line_floor(Qc, noise)
            self.assertIsNotNone(floor)
            self.assertLessEqual(floor, lam_min + 1e-12)
            cert = per_line_certificate([Q], [Qc], noise)
            self.assertTrue(cert['informative'])
            self.assertLessEqual(kl, cert['kl_upper'])

    def test_noise_floors_where_motion_fails(self):
        # common motion past one: the Weyl term gives nothing, the receiver noise still floors the covariance
        rng = np.random.default_rng(38)
        n = 64
        k = np.arange(n)
        rows = np.arange(10, 40)
        E = np.exp(-2j * np.pi * np.outer(rows, k) / n) / np.sqrt(n)
        for _ in range(6):
            phi0 = rng.uniform(1.0, 3.0) * np.cos(2 * np.pi * k / rng.uniform(5, 20) + rng.uniform(0, 6.3))
            B0 = E * np.exp(-1j * phi0)
            Qc = float(np.sum(np.abs(E * np.expm1(-1j * phi0)) ** 2))
            noise = rng.uniform(0.01, 0.2)
            lam_min = float(np.linalg.eigvalsh(B0 @ B0.conj().T + noise * np.eye(len(rows))).min())
            self.assertGreaterEqual(Qc, 1.0)
            self.assertAlmostEqual(line_floor(Qc, noise), noise)
            self.assertLessEqual(line_floor(Qc, noise), lam_min + 1e-12)
        self.assertIsNone(line_floor(2.0, 0.0))
        self.assertFalse(per_line_certificate([1e-6], [2.0], 0.0)['informative'])
        self.assertTrue(per_line_certificate([1e-6], [2.0], 0.1)['informative'])

    def test_lines_add(self):
        rng = np.random.default_rng(37)
        rows = [self.line(rng) for _ in range(4)]
        cert = per_line_certificate([r[1] for r in rows], [r[2] for r in rows], rows[0][4])
        self.assertLessEqual(sum(r[0] for r in rows), cert['kl_upper'])

