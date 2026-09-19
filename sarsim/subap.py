"""Doppler sub-aperture filter bank: K reference/offset pairs swept across the processed band.

Each reference band spans `sub_frac` of the processed azimuth bandwidth; the offset band is the
same width shifted by `delta_hz` (converted to cycles/m with the platform velocity). The pair is
swept rigidly from the low edge to the high edge of the band in K steps, as in the 2022 paper's
Figure 4 and the 2026 replication protocol (K = 50, B_shift = 88 Hz).

Band edges live on the discrete FFT grid of the image. With `snap_bins=True` the offset and the
sweep step are rounded to whole bins so that every pair differs from its neighbours by exactly
the same number of bins; with hard edges and no snapping the bands alternately gain and lose a
bin as the fractional edge advances, which injects a spurious periodic term into every
trajectory (see test T4). `taper_bins` > 0 replaces the hard edge by a raised-cosine ramp.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class SubapBank:
    K: int = 50
    sub_frac: float = 0.5
    delta_hz: float = 88.0
    snap_bins: bool = True
    taper_bins: int = 0

    def bands(self, geom, Nx=None):
        B = geom.nu_band
        Bs = self.sub_frac * B
        dnu = self.delta_hz / geom.V
        if self.snap_bins and Nx is not None:
            dbin = 1.0 / (Nx * geom.dx)
            dnu = max(1, round(dnu / dbin)) * dbin
            Bs = round(Bs / dbin) * dbin
            step = np.floor((B - Bs - dnu) / (self.K - 1) / dbin) * dbin
            start = np.ceil(-B / 2 / dbin) * dbin
            lo = start + np.arange(self.K) * step
        else:
            lo = np.linspace(-B / 2, B / 2 - Bs - dnu, self.K)
        return {"nu_lo_ref": lo, "nu_lo_off": lo + dnu, "width": Bs, "dnu": dnu,
                "nu_c": lo + (Bs + dnu) / 2}

    def _mask(self, nu, lo, width):
        m = ((nu[None, :] >= lo[:, None] - 1e-12) & (nu[None, :] < lo[:, None] + width - 1e-12)).astype(float)
        if self.taper_bins > 0:
            dbin = nu[1] - nu[0]
            t = self.taper_bins * dbin
            ramp_lo = np.clip((nu[None, :] - lo[:, None]) / t, 0, 1)
            ramp_hi = np.clip((lo[:, None] + width - nu[None, :]) / t, 0, 1)
            m = 0.5 * (1 - np.cos(np.pi * ramp_lo)) * 0.5 * (1 - np.cos(np.pi * ramp_hi))
            m[(nu[None, :] < lo[:, None]) | (nu[None, :] >= lo[:, None] + width)] = 0
        return m

    def masks(self, geom, Nx):
        b = self.bands(geom, Nx)
        nu = np.fft.fftfreq(Nx, d=geom.dx)
        order = np.argsort(nu)
        nus = nu[order]
        mR = np.empty((self.K, Nx))
        mO = np.empty((self.K, Nx))
        mR[:, order] = self._mask(nus, b["nu_lo_ref"], b["width"])
        mO[:, order] = self._mask(nus, b["nu_lo_off"], b["width"])
        return mR, mO

    def pair(self, spec_az, k, masks):
        mR, mO = masks
        R = np.fft.ifft(spec_az * mR[k][:, None], axis=0)
        O = np.fft.ifft(spec_az * mO[k][:, None], axis=0)
        return R, O

    def sub_image(self, spec_az, nu_lo, width, geom):
        """A single sub-aperture image for an arbitrary band (used for illustrations)."""
        nu = np.fft.fftfreq(spec_az.shape[0], d=geom.dx)
        m = (nu >= nu_lo) & (nu < nu_lo + width)
        return np.fft.ifft(spec_az * m[:, None], axis=0)
