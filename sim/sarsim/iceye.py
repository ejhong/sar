"""Groundwork for the real-data phase: read an ICEYE SLC (HDF5) into the toolkit's conventions.

ICEYE SLC products store the complex image as two datasets `s_i` and `s_q` (azimuth lines ×
range samples) plus scalar metadata as root-level datasets/attributes. Key names have varied
between product versions, so `inspect()` lists everything and `load()` maps the usual names with
fallbacks. Anything not found is reported rather than guessed; verify against the product's
JSON/XML metadata before trusting a geometry.

Typical usage (after the 9.8 GB file is available):
    from sarsim.iceye import inspect, load
    inspect("ICEYE_X33_SLC_SLEDF_951562307_....h5")          # prints keys, shapes, scalars
    slc, geom, meta = load(path, rows=(r0, r1), cols=(c0, c1))  # cropped complex64 [az, rg]
"""
import json
import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None

from .geometry import Geometry

CANDIDATES = {
    "carrier_frequency": ["carrier_frequency", "center_frequency", "radar_center_frequency"],
    "azimuth_spacing": ["azimuth_ground_spacing", "azimuth_spacing", "azimuth_pixel_spacing"],
    "range_spacing": ["slant_range_spacing", "range_spacing", "range_pixel_spacing"],
    "incidence_center": ["incidence_center", "incidence_angle_center", "incidence_angle"],
    "prf": ["processing_prf", "acquisition_prf", "prf"],
    "doppler_bandwidth": ["total_processed_bandwidth_azimuth", "processed_azimuth_bandwidth", "azimuth_bandwidth"],
    "chirp_bandwidth": ["chirp_bandwidth", "range_bandwidth"],
    "look_side": ["look_side", "antenna_pointing"],
    "orbit_direction": ["orbit_direction", "orbit_pass"],
    "start": ["acquisition_start_utc", "zerodoppler_start_utc", "acquisition_start"],
    "end": ["acquisition_end_utc", "zerodoppler_end_utc", "acquisition_end"],
    "slant_range_first": ["slant_range_to_first_pixel", "first_pixel_slant_range", "near_range"],
    "velocity": ["velocity", "platform_velocity", "satellite_velocity"],
    "product_type": ["product_type"],
    "mode": ["acquisition_mode", "imaging_mode", "mode"],
    "coord_center": ["coord_center", "center_coordinates"],
    "first_pixel_time": ["first_pixel_time", "first_pixel_time_utc"],
    "doppler_rate": ["doppler_rate_coeffs", "doppler_rate", "azimuth_fm_rate"],
    "doppler_centroid": ["dc_estimate_coeffs", "doppler_centroid", "dc_coeffs"],
}


def _scalar(v):
    v = np.asarray(v)
    if v.dtype.kind in "SU" or v.dtype == object:
        try:
            return v.item().decode() if isinstance(v.item(), bytes) else v.item()
        except Exception:
            return str(v)
    if v.size == 1:
        return v.item()
    return v.tolist() if v.size <= 12 else f"array{v.shape}"


def inspect(path, max_items=200):
    """Print every dataset and attribute (shape, dtype, small values) so the mapping can be checked."""
    out = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                if obj.size <= 12 or obj.dtype.kind in "SU" or obj.dtype == object:
                    try:
                        out[name] = _scalar(obj[()])
                    except Exception:
                        out[name] = f"dataset {obj.shape} {obj.dtype}"
                else:
                    out[name] = f"dataset {obj.shape} {obj.dtype}"
        f.visititems(visit)
        for k, v in f.attrs.items():
            out[f"@{k}"] = _scalar(v)
    for i, (k, v) in enumerate(out.items()):
        if i >= max_items:
            print(f"... {len(out) - max_items} more")
            break
        print(f"{k}: {v}")
    return out


def _find(f, keys):
    for k in keys:
        if k in f:
            return _scalar(f[k][()])
        if k in f.attrs:
            return _scalar(f.attrs[k])
    return None


def load(path, rows=None, cols=None, sidecar_json=None):
    """Return (slc complex64 [azimuth, range], Geometry, metadata dict) for an optional crop."""
    meta = {}
    with h5py.File(path, "r") as f:
        for name, keys in CANDIDATES.items():
            meta[name] = _find(f, keys)
        si, sq = f["s_i"], f["s_q"]
        r0, r1 = (0, si.shape[0]) if rows is None else rows
        c0, c1 = (0, si.shape[1]) if cols is None else cols
        slc = (si[r0:r1, c0:c1].astype(np.float32) + 1j * sq[r0:r1, c0:c1].astype(np.float32)).astype(np.complex64)
        meta["shape_full"] = tuple(si.shape)
        meta["crop"] = {"rows": (r0, r1), "cols": (c0, c1)}
    if sidecar_json:
        with open(sidecar_json) as fj:
            meta["sidecar"] = json.load(fj)
    missing = [k for k in ("carrier_frequency", "azimuth_spacing", "range_spacing", "incidence_center") if meta.get(k) is None]
    if missing:
        raise KeyError(f"could not find {missing} in {path}; run inspect() and extend CANDIDATES")
    geom = Geometry(f0=float(meta["carrier_frequency"]), theta_deg=float(meta["incidence_center"]),
                    dx=float(meta["azimuth_spacing"]), dr=float(meta["range_spacing"]))
    if meta.get("velocity") is not None and np.isscalar(meta["velocity"]):
        geom.V = float(meta["velocity"])
    if meta.get("slant_range_first") is not None and np.isscalar(meta["slant_range_first"]):
        geom.R0 = float(meta["slant_range_first"])
    # processed bandwidth fraction from the Doppler bandwidth and the azimuth sampling, if present
    if meta.get("doppler_bandwidth") is not None and np.isscalar(meta["doppler_bandwidth"]):
        nyq = geom.V / geom.dx
        geom.band_frac = float(min(0.99, meta["doppler_bandwidth"] / nyq))
    return slc, geom, meta


def estimate_band_fraction(slc, axis=0, thresh_db=-10.0):
    """Estimate the processed bandwidth fraction from the spectrum itself (fallback when metadata is missing)."""
    p = (np.abs(np.fft.fft(slc, axis=axis)) ** 2).mean(axis=1 - axis)
    p_db = 10 * np.log10(p / p.max() + 1e-12)
    return float((p_db > thresh_db).mean())
