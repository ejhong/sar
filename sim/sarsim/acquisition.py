"""The real dwell as the simulator's geometry (Phase 2, P2-01).

The first investigation's simulator used a flat-earth spotlight geometry with one velocity for
everything. A real ICEYE dwell has three:

* V_s, the platform's speed along its orbit (7,661 m/s at Giza);
* V_g, the speed at which the image's zero-Doppler rows sweep the ground: row spacing over row
  time (4.40 cm / 6.27 us = 7,014 m/s), which sets how Doppler maps to image position;
* V_eff = sqrt(V_s V_g), which sets the Doppler rate Ka = -2 V_eff^2 / (lambda R).

The flat-earth model forces V_s = V_g = V_eff and misplaces Ka by nine per cent either way.
`DwellGeometry` keeps all three and the product's own sign convention for slow time,
f = f_dc + Ka (t - t_zd) with Ka < 0, so positive Doppler is earlier. It is built from a product
(`from_product`, desktop only) or from the acquisition record the product leaves in
sites/acquisitions (`from_record`), so everything downstream runs without the imagery.

Two consequences of the real geometry are used throughout Phase 2:

* A look W seconds long spans |Ka| W hertz of Doppler, so its azimuth resolution is
  0.886 V_g / (|Ka| W): resolution times duration is fixed by the orbit (1.24 m s at Giza).
* A line-of-sight velocity v shifts a look in azimuth by 2 V_g v / (lambda |Ka|) = R v / V_s.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

from .geometry import Geometry, C

RECORDS = Path(__file__).resolve().parents[2] / 'sites' / 'acquisitions'


@dataclass
class DwellGeometry(Geometry):
    """Geometry of a real dwell. `V` (inherited) is the ground-sweep speed V_g, so f = V nu holds for the
    image's own sampling; `Ka_signed` is the product's Doppler rate (negative); `band_frac` (inherited) is
    the azimuth processed fraction and `band_frac_r` the range one."""
    V_platform: float = 7600.0
    Ka_signed: float = -5000.0
    band_frac_r: float = 0.8
    prf: float = 0.0
    row_dt: float = 0.0
    heading_deg: float = 0.0
    los_enu: tuple = (0.0, 0.0, 1.0)             # unit vector from the ground towards the satellite
    baseline_coeffs: tuple = ()                  # B_perp(t) = sum c_i t^i, metres, from the state vectors
    name: str = ''
    source: dict = field(default_factory=dict)

    # ---- Doppler and slow time
    @property
    def Ka(self):
        """Magnitude of the Doppler rate (Hz/s), as the flat-earth class defines it."""
        return abs(self.Ka_signed)

    @property
    def kr_band(self):
        return self.band_frac_r / self.dr

    @property
    def V_eff(self):
        return float(np.sqrt(self.Ka * self.lam * self.R0 / 2))

    def nu_to_time(self, nu):
        """Slow time relative to each scatterer's zero-Doppler crossing, product convention (Ka < 0)."""
        return self.V * np.asarray(nu, float) / self.Ka_signed

    def nu_to_baseline(self, nu):
        """Virtual baseline of the look centred at nu: the platform's offset from the aperture centre,
        along track and perpendicular to the line of sight (from the state vectors when recorded)."""
        t = self.nu_to_time(nu)
        if len(self.baseline_coeffs):
            return np.polynomial.polynomial.polyval(t, np.asarray(self.baseline_coeffs, float))
        return self.V_platform * t

    def zero_doppler_time(self, x):
        """Zero-Doppler time (s) of azimuth position x (m) relative to x = 0: rows advance with time."""
        return np.asarray(x, float) / self.V

    # ---- looks
    @property
    def time_space_m_s(self):
        """V_g / |Ka|: a look's azimuth resolution times its duration, over 0.886."""
        return self.V / self.Ka

    def look_resolution(self, duration_s):
        return 0.886 * self.time_space_m_s / np.asarray(duration_s, float)

    def look_duration(self, resolution_m):
        return 0.886 * self.time_space_m_s / np.asarray(resolution_m, float)

    def shift_to_velocity(self, shift_m):
        """Slant-range rate (m/s, receding positive) that shifts a look's image by shift_m in azimuth: the
        motion exp(-4 pi j v t / lambda) with t = V nu / Ka is a phase ramp in nu, a shift of 2 V v / (lambda Ka)."""
        return np.asarray(shift_m, float) * self.lam * self.Ka_signed / (2 * self.V)

    # ---- site frame
    @property
    def along_track_en(self):
        h = np.deg2rad(self.heading_deg)
        return np.array([np.sin(h), np.cos(h)])

    @property
    def ground_range_en(self):
        """Horizontal direction away from the radar (east, north)."""
        los = np.asarray(self.los_enu, float)
        g = -los[:2]
        return g / np.linalg.norm(g)

    def site_to_radar(self, east, north):
        """Site (east, north) metres to the simulator's (azimuth x, ground range y) metres."""
        p = np.stack([np.asarray(east, float), np.asarray(north, float)], axis=-1)
        return p @ self.along_track_en, p @ self.ground_range_en

    def slant_range_increase(self, u_enu):
        """Increase in slant range (m) for ground displacements u_enu [..., 3] (east, north, up)."""
        return -np.asarray(u_enu, float) @ np.asarray(self.los_enu, float)

    def absolute_motion(self, motion_enu):
        """Wrap `motion_enu(x, y, t_abs) -> [len(t), len(x), 3]` displacements (site frame, metres), where
        x, y are site east and north and t_abs is time from the aperture centre, into the synthesizer's
        `motion(x, y, z, t)` callback: its x, y are radar azimuth and ground range and t is slow time relative
        to each scatterer's own zero-Doppler crossing."""
        a, g = self.along_track_en, self.ground_range_en

        def motion(xr, yr, z, t):
            east = xr * a[0] + yr * g[0]
            north = xr * a[1] + yr * g[1]
            t_abs = np.asarray(t, float)[:, None] + self.zero_doppler_time(xr)[None, :]
            return self.slant_range_increase(motion_enu(east, north, t_abs))
        return motion

    # ---- records
    def record(self):
        d = asdict(self)
        d['los_enu'] = list(self.los_enu)
        d['baseline_coeffs'] = list(self.baseline_coeffs)
        d['derived'] = {
            'wavelength_m': self.lam, 'ground_sweep_speed_m_s': self.V, 'effective_speed_m_s': self.V_eff,
            'aperture_time_s': self.aperture_time, 'doppler_band_hz': self.doppler_band_hz,
            'azimuth_resolution_m': self.resolution, 'slant_resolution_m': 0.886 / self.kr_band,
            'time_space_m_s': self.time_space_m_s, 'aspect_span_deg': self.aspect_span_deg,
        }
        return d

    @classmethod
    def from_record(cls, rec):
        if isinstance(rec, (str, Path)):
            p = Path(rec)
            if not p.suffix:
                p = RECORDS / f'{rec}.json'
            rec = json.loads(p.read_text())
        keys = {f for f in cls.__dataclass_fields__}
        d = {k: v for k, v in rec.items() if k in keys}
        d['los_enu'] = tuple(d.get('los_enu', (0, 0, 1)))
        d['baseline_coeffs'] = tuple(d.get('baseline_coeffs', ()))
        return cls(**d)

    @classmethod
    def from_product(cls, path, name=''):
        """Read a dwell product's acquisition (metadata only; no imagery is loaded)."""
        import h5py
        from .dwell import DwellProduct
        from .orbit import lla_to_ecef, perpendicular_baselines
        p = DwellProduct(path)
        a = p.acq
        with h5py.File(path, 'r') as f:
            chirp = float(f['chirp_bandwidth'][()])
            fs = float(f['range_sampling_rate'][()])
            col_c, row_c, lat, lon = (float(v) for v in np.ravel(f['coord_center'][()]))
            heading = float(f['heading'][()])
            wr = f['window_function_range'][()]
            wa = f['window_function_azimuth'][()]
        col = col_c - 1                                          # the product's coordinates count from 1
        ka = float(p.doppler_rate_hz_s(col))
        R0 = float(a.slant_range_m(col))
        eph = p.ephemeris()
        pos, vel = eph.at([0.0])
        target = lla_to_ecef(lat, lon, a.scene_height_m)
        T = a.bandwidth_hz / abs(ka)
        times = np.linspace(-T / 2, T / 2, 49)
        B, R_ref, inc = perpendicular_baselines(eph, times, target, 0.0)
        coeffs = np.polynomial.polynomial.polyfit(times, B, 3)
        resid = B - np.polynomial.polynomial.polyval(times, coeffs)
        # line of sight in the local east-north-up frame at the scene centre
        la, lo = np.deg2rad(lat), np.deg2rad(lon)
        E = np.array([-np.sin(lo), np.cos(lo), 0.0])
        N = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
        U = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
        d = pos[0] - target
        d /= np.linalg.norm(d)
        los = (float(d @ E), float(d @ N), float(d @ U))
        dx = a.azimuth_spacing_m
        g = cls(f0=C / a.wavelength_m, R0=R0, theta_deg=a.incidence_center_deg,
                V=dx / a.azimuth_time_interval_s, dx=dx, dr=a.range_spacing_m,
                band_frac=a.bandwidth_hz / a.prf_hz, V_platform=float(np.linalg.norm(vel[0])), Ka_signed=ka,
                band_frac_r=chirp / fs, prf=a.prf_hz, row_dt=a.azimuth_time_interval_s, heading_deg=heading,
                los_enu=los, baseline_coeffs=tuple(float(c) for c in coeffs), name=name or Path(path).stem,
                source={'product': Path(path).name, 'satellite': a.satellite, 'mode': a.product_type,
                        'polarization': a.polarization, 'look_side': a.look_side,
                        'orbit_direction': a.orbit_direction, 'collection_start': a.collection_start,
                        'collection_end': a.collection_end, 'collection_duration_s': a.collection_duration_s,
                        'image_span_s': a.image_span_s, 'shape': list(a.shape),
                        'scene_centre': {'row': row_c, 'column': col_c, 'latitude': lat, 'longitude': lon,
                                         'height_m': a.scene_height_m},
                        'chirp_bandwidth_hz': chirp, 'range_sampling_rate_hz': fs,
                        'window_range': wr.decode() if isinstance(wr, bytes) else str(wr),
                        'window_azimuth': wa.decode() if isinstance(wa, bytes) else str(wa),
                        'incidence_from_state_vectors_deg': inc, 'slant_range_from_state_vectors_m': R_ref,
                        'baseline_span_m': float(B.max() - B.min()),
                        'baseline_fit_residual_m': float(np.abs(resid).max()),
                        'fields': ['processing_prf', 'azimuth_time_interval', 'total_processed_bandwidth_azimuth',
                                   'carrier_frequency', 'azimuth_ground_spacing', 'slant_range_spacing',
                                   'slant_range_to_first_pixel', 'incidence_center', 'doppler_rate_coeffs',
                                   'chirp_bandwidth', 'range_sampling_rate', 'coord_center', 'heading',
                                   'state vectors (posX..velZ, state_vector_time_utc)']})
        return g
