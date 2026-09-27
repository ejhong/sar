"""Scatterer scenes: a stepped pyramid, desert rocks, vibrating patches, explicit point targets.

Every scatterer carries: position (x, y, z), complex reflectivity (amp, phase), an aspect
response g(nu) = iso + flash * exp(-(nu - nu0)^2 / (2 sig_nu^2)) in azimuth spatial frequency
nu (cycles/m), an optional line-of-sight vibration d(t) = vib_amp sin(2 pi vib_freq t + vib_phase),
and an integer label.
"""
from dataclasses import dataclass
import numpy as np

FIELDS = ("x", "y", "z", "amp", "phase", "iso", "flash", "nu0", "sig_nu",
          "vib_amp", "vib_freq", "vib_phase", "label")


@dataclass
class Scatterers:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    amp: np.ndarray
    phase: np.ndarray
    iso: np.ndarray
    flash: np.ndarray
    nu0: np.ndarray
    sig_nu: np.ndarray
    vib_amp: np.ndarray
    vib_freq: np.ndarray
    vib_phase: np.ndarray
    label: np.ndarray

    @property
    def n(self):
        return int(self.x.size)

    def subset(self, mask):
        return Scatterers(**{f: getattr(self, f)[mask] for f in FIELDS})

    def copy(self):
        return Scatterers(**{f: getattr(self, f).copy() for f in FIELDS})


def make(x, y, z, amp, phase=None, iso=1.0, flash=0.0, nu0=0.0, sig_nu=1.0,
         vib_amp=0.0, vib_freq=0.0, vib_phase=0.0, label=0, rng=None):
    x = np.asarray(x, float).ravel()
    n = x.size
    rng = np.random.default_rng(rng)
    if phase is None:
        phase = rng.uniform(0, 2 * np.pi, n)

    def b(v, dtype=float):
        v = np.asarray(v, dtype)
        return np.broadcast_to(v, (n,)).copy() if v.ndim == 0 else v.ravel().astype(dtype)

    return Scatterers(x, b(y), b(z), b(amp), b(phase), b(iso), b(flash), b(nu0), b(sig_nu),
                      b(vib_amp), b(vib_freq), b(vib_phase), b(label, int))


def concat(*scats):
    scats = [s for s in scats if s is not None and s.n > 0]
    return Scatterers(**{f: np.concatenate([getattr(s, f) for s in scats]) for f in FIELDS})


def rot2d(x, y, deg):
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return x * c - y * s, x * s + y * c


class PyramidShape:
    """Height field of a square pyramid (used for occlusion tests)."""

    def __init__(self, base, height, rotation_deg=0.0, centre=(0.0, 0.0)):
        self.base, self.height, self.rot, self.centre = base, height, rotation_deg, centre

    def __call__(self, x, y):
        xp, yp = rot2d(np.asarray(x) - self.centre[0], np.asarray(y) - self.centre[1], -self.rot)
        h = self.height * (1 - np.maximum(np.abs(xp), np.abs(yp)) / (self.base / 2))
        return np.clip(h, 0, None)


def occluded(scat, height_fn, theta_rad, s_min=0.3, s_max=400.0, n_samples=96):
    """True where the ray from the scatterer toward the radar (-y, rising at 1/tan theta) hits terrain."""
    s = np.linspace(s_min, s_max, n_samples)
    slope = 1.0 / np.tan(theta_rad)
    out = np.zeros(scat.n, bool)
    for i in range(0, scat.n, 20000):
        sl = slice(i, i + 20000)
        xs = scat.x[sl][:, None]
        ys = scat.y[sl][:, None] - s[None, :]
        zr = scat.z[sl][:, None] + slope * s[None, :]
        out[sl] = np.any(height_fn(xs, ys) > zr + 1e-3, axis=1)
    return out


def _flash_params(face_angle_deg, n, geom, tilt_sigma_deg, block_len, rng):
    """Flash centre (cycles/m) and width for vertical block facets whose horizontal normal points
    `face_angle_deg` away from the radar's ground-projected line of sight, with random tilts."""
    psi = np.deg2rad(face_angle_deg + rng.normal(0, tilt_sigma_deg, n))
    nu0 = 2 * np.sin(psi) / geom.lam
    D = rng.uniform(block_len[0], block_len[1], n)
    sig = 0.4 / D
    return nu0, sig


def stepped_pyramid(geom, base=215.0, height=136.0, course=1.0, rotation_deg=0.0, spacing=0.25,
                    centre=(0.0, 0.0), amp_face=(1.0, 0.04, 0.35, 0.35), amp_sigma=0.5,
                    flash=0.8, iso=0.35, tilt_sigma_deg=3.0, block_len=(0.6, 2.0),
                    diffuse_density=1.5, amp_diffuse=0.10, occlude=True, rng=0):
    """A stepped square pyramid of `height/course` courses.

    Faces (before rotation): 0 = near (normal toward the radar), 1 = far, 2 = -x side, 3 = +x side.
    Course-edge scatterers sit on the dihedral corner at the foot of each riser; they carry an
    isotropic term plus a specular 'flash' whose aspect position follows the face orientation.
    Diffuse scatterers are scattered over the faces with isotropic response.
    """
    rng = np.random.default_rng(rng)
    L2 = base / 2
    normals = {0: (0.0, -1.0), 1: (0.0, 1.0), 2: (-1.0, 0.0), 3: (1.0, 0.0)}
    parts = []
    n_courses = int(np.floor(height / course + 1e-9))
    for c in range(n_courses):
        zc = c * course
        w = L2 * (1 - zc / height)
        if w < spacing:
            break
        n_edge = max(int(round(2 * w / spacing)), 2)
        s = (np.arange(n_edge) + 0.5) / n_edge * 2 * w - w
        ones = np.ones_like(s)
        edges = {0: (s, -w * ones), 1: (s, w * ones), 2: (-w * ones, s), 3: (w * ones, s)}
        for face, (ex, ey) in edges.items():
            nx, ny = rot2d(*normals[face], rotation_deg)
            face_angle = np.rad2deg(np.arctan2(nx, -ny))  # 0 = facing the radar
            n = ex.size
            nu0, sig = _flash_params(face_angle, n, geom, tilt_sigma_deg, block_len, rng)
            xr, yr = rot2d(ex, ey, rotation_deg)
            amp = amp_face[face] * np.exp(rng.normal(0, amp_sigma, n))
            parts.append(make(xr + centre[0], yr + centre[1], zc, amp, iso=iso, flash=flash,
                              nu0=nu0, sig_nu=sig, label=face, rng=rng))
    # diffuse face scatterers
    slant = np.hypot(height, L2)
    area = base * slant / 2
    n_diff = int(diffuse_density * area)
    for face in range(4):
        u = rng.uniform(0, 1, n_diff)
        zf = height * (1 - np.sqrt(u))
        w = L2 * (1 - zf / height)
        t = rng.uniform(-1, 1, n_diff) * w
        if face == 0:
            ex, ey = t, -w
        elif face == 1:
            ex, ey = t, w
        elif face == 2:
            ex, ey = -w, t
        else:
            ex, ey = w, t
        xr, yr = rot2d(ex, ey, rotation_deg)
        amp = amp_diffuse * (amp_face[face] / amp_face[0]) * np.exp(rng.normal(0, amp_sigma, n_diff))
        parts.append(make(xr + centre[0], yr + centre[1], zf, amp, iso=1.0, flash=0.0,
                          label=10 + face, rng=rng))
    scat = concat(*parts)
    if occlude:
        shape = PyramidShape(base, height, rotation_deg, centre)
        occ = occluded(scat, shape, geom.theta)
        scat = scat.subset(~occ)
    return scat


def desert_rocks(geom, extent, density=0.02, amp=0.25, amp_sigma=0.7, rng=1):
    """Sparse isotropic point-like rocks on flat ground; extent = (x0, x1, y0, y1)."""
    rng = np.random.default_rng(rng)
    x0, x1, y0, y1 = extent
    n = int(density * (x1 - x0) * (y1 - y0))
    x = rng.uniform(x0, x1, n)
    y = rng.uniform(y0, y1, n)
    a = amp * np.exp(rng.normal(0, amp_sigma, n))
    return make(x, y, 0.0, a, label=20, rng=rng)


def vibrating_patch(centre, size, density=6.0, amp=0.3, vib_amp=1e-3, vib_freq=1.0,
                    vib_phase=0.0, amp_sigma=0.5, label=30, rng=2):
    """A square ground patch whose scatterers all move coherently along the line of sight."""
    rng = np.random.default_rng(rng)
    n = int(density * size * size)
    x = centre[0] + rng.uniform(-size / 2, size / 2, n)
    y = centre[1] + rng.uniform(-size / 2, size / 2, n)
    a = amp * np.exp(rng.normal(0, amp_sigma, n))
    return make(x, y, 0.0, a, vib_amp=vib_amp, vib_freq=vib_freq, vib_phase=vib_phase,
                label=label, rng=rng)


def point_targets(x, y, z, amp, phase=None, label=40, rng=3, **kw):
    return make(x, y, z, amp, phase=phase, label=label, rng=rng, **kw)
