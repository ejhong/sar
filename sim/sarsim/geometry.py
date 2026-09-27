from dataclasses import dataclass
import numpy as np

C = 299792458.0


@dataclass
class Geometry:
    """Side-looking spotlight geometry with a straight-line, flat-earth approximation.

    Coordinates: x = along-track (azimuth), y = ground range (positive away from the
    radar), z = height. The radar looks in +y from height, at incidence `theta_deg`.
    Slant range of a scatterer relative to the scene centre: r = y sin(theta) - z cos(theta).
    Azimuth spatial frequency nu (cycles/m) <-> Doppler f = V nu <-> slow time t = f/Ka
    <-> platform along-track offset B = V t <-> aspect angle psi = nu lambda / 2.
    """
    f0: float = 9.65e9        # carrier (Hz), X-band
    R0: float = 600e3         # slant range to scene centre (m)
    theta_deg: float = 35.0   # incidence angle
    V: float = 7600.0         # effective platform velocity (m/s)
    dx: float = 0.25          # azimuth pixel spacing (m)
    dr: float = 0.25          # slant-range pixel spacing (m)
    band_frac: float = 0.8    # processed bandwidth as a fraction of Nyquist, both axes

    @property
    def lam(self):
        return C / self.f0

    @property
    def theta(self):
        return np.deg2rad(self.theta_deg)

    @property
    def Ka(self):
        """Magnitude of the Doppler rate (Hz/s)."""
        return 2 * self.V ** 2 / (self.lam * self.R0)

    @property
    def nu_band(self):
        """Processed azimuth spatial bandwidth (cycles/m)."""
        return self.band_frac / self.dx

    @property
    def kr_band(self):
        return self.band_frac / self.dr

    @property
    def doppler_band_hz(self):
        return self.V * self.nu_band

    @property
    def aperture_time(self):
        return self.doppler_band_hz / self.Ka

    @property
    def aspect_span_deg(self):
        return np.rad2deg(self.nu_band * self.lam / 2)

    @property
    def resolution(self):
        """Approximate 3 dB resolution of a rectangular band (m)."""
        return 0.886 / self.nu_band

    def nu_to_time(self, nu):
        return self.V * np.asarray(nu, float) / self.Ka

    def nu_to_baseline(self, nu):
        """Along-track platform offset (m) of the sub-aperture centred at nu."""
        return self.V * self.nu_to_time(nu)

    def nu_to_aspect_deg(self, nu):
        return np.rad2deg(np.asarray(nu, float) * self.lam / 2)

    def kz(self, nu, lam_s):
        """Biondi/Seyfzadeh steering wavenumber for the sub-aperture centred at nu.

        Kz = 4 pi B_perp / (lam_s R0 sin theta), with B_perp the along-track platform
        offset and lam_s the declared "sound" wavelength (0.48 m in the 2022 paper).
        """
        return 4 * np.pi * self.nu_to_baseline(nu) / (lam_s * self.R0 * np.sin(self.theta))

    def depth_scale(self, lam_s):
        """Model-metres of 'depth' per metre of azimuth separation between interfering scatterers."""
        return lam_s * np.sin(self.theta) / self.lam

    def summary(self):
        return {
            "carrier_GHz": self.f0 / 1e9, "wavelength_m": self.lam, "slant_range_km": self.R0 / 1e3,
            "incidence_deg": self.theta_deg, "velocity_m_s": self.V, "pixel_az_m": self.dx,
            "pixel_rg_m": self.dr, "doppler_band_kHz": self.doppler_band_hz / 1e3,
            "aperture_time_s": self.aperture_time, "aspect_span_deg": self.aspect_span_deg,
            "resolution_m": self.resolution,
        }
