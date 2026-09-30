"""Radiometry and a Khafre power map for each Giza pass, read from the ICEYE products on the desktop.

    uv run python -m katabasis.compose.acquisition_radiometry

For giza-20250827 and giza-20220715:

1. A `radiometry` block appended to sites/acquisitions/<name>.json: every calibration field the product carries, verbatim
   (the HDF5's calibration_factor and compensation flags, and for the 2025 pass the SICD XML's radiometric scale
   factors), the incidence across the swath, and what is known of the noise floor. Neither product carries a noise-floor
   field (no NESZ in the HDF5 or the SLC XML; no Radiometric/NoiseLevel in the 2025 SICD XML; no SICD in the 2022
   delivery), so the block says so, gives ICEYE's specified range, and adds an upper bound measured on the image: the
   darkest 10 m cells over the Giza site, whose power holds the noise floor and whatever ground signal is left there.
2. sites/acquisitions/<name>-khafre-power.npz: |z|^2 averaged over cells of about 1.5 m (whole pixels, in the image's
   own geometry: lines along track, samples in slant range) over Khafre's footprint, a 320 m square round the base centre
   (the Giza site's origin), with its grid, the calibration to beta0 and sigma0, and a georeference: the product's RPC,
   the geoid undulation and the residual offset fitted for that pass (P2-12 for 2025, P2-13 for 2022), and tie points
   every 10 m over the footprint giving each ground point's site-frame, geographic and image coordinates.

Only these two passes, and only on the desktop where the products are; P2-01 carries the block forward when it rewrites
a record.
"""
from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RECORDS = ROOT / 'sites' / 'acquisitions'
RESULTS = ROOT / 'sim' / 'results'
PASSES = {
    'giza-20250827': {'product': Path.home() / 'tmp/sar/giza/ICEYE_X33_SLC_SLEDF_951562307_20250827T202654.h5',
                      'delivery': Path.home() / 'tmp/sar/giza/drive-download-20260920T143751Z-1-001.zip',
                      'offset_from': 'p2_12_real_pass'},
    'giza-20220715': {'product': Path.home() / 'tmp/sar/giza2/ICEYE_X13_SLC_SLED_868226_20220715T235744.h5',
                      'delivery': Path.home() / 'tmp/sar/giza2/drive-download-20260928T115627Z-1-001.zip',
                      'offset_from': 'p2_13_gates'},
}
HALF_M = 160.0              # Khafre's base is 215 m; a 320 m square round its centre, the site's origin
CELL_M = 1.5                # the power map's cells, about 1.5 m on the ground each way (Khafre's top lays over 360 m at 21 deg)
TIE_M = 10.0
DARK_CELL_M = 10.0          # cells for the noise-floor bound
SPEC = {'nesz_db': [-18.0, -15.0], 'source': 'ICEYE Product Documentation 6.0.0, Product Specification, 2. Imaging Modes, '
        'Table 2-11 (Dwell and Dwell Fine): "Noise Equivalent Sigma-Zero [dBm2/m2]: -18 to -15" (scene-centre values)'}
H5_FIELDS = ['calibration_factor', 'range_spread_comp_flag', 'ant_elev_corr_flag', 'antenna_pattern_compensation',
             'sample_precision', 'window_function_azimuth', 'window_function_range', 'incidence_near', 'incidence_center',
             'incidence_far', 'avg_scene_height', 'processor_version', 'spec_version', 'product_name']


def plain(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray):
        return v.tolist() if v.size <= 16 else None
    if isinstance(v, np.generic):
        return v.item()
    return v


def product_fields(path: Path) -> dict:
    with h5py.File(path, 'r') as f:
        out = {k: plain(f[k][()]) for k in H5_FIELDS if k in f}
        fsl = np.asarray(f['fsl_compensation'][()])
        out['fsl_compensation'] = {'samples': int(fsl.size), 'all_zero': bool(np.all(fsl == 0))}
        lia = np.asarray(f['local_incidence_angle'][()], float)
        out['local_incidence_angle_deg'] = {'first': float(lia[0]), 'last': float(lia[-1]), 'samples': int(lia.size)}
    return out


def delivery_fields(path: Path) -> dict:
    """The SLC XML's calibration factor and, where the delivery has one, the SICD XML's radiometric block."""
    out = {'files': []}
    if not path.exists():
        return out
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            if not name.endswith('.xml') or name.endswith('.aux.xml'):
                continue
            text = z.read(name).decode('utf-8', 'replace')
            out['files'].append(name)
            if '_SLC_' in name:
                m = re.search(r'<calibration_factor>([^<]+)</calibration_factor>', text)
                out['slc_xml_calibration_factor'] = float(m.group(1)) if m else None
                out['slc_xml_noise_fields'] = sorted(set(re.findall(r'<([A-Za-z_]*(?:noise|nesz|NESZ|Noise)[A-Za-z_]*)>', text)))
            if '_SICD_' in name:
                rad = {}
                for tag in ('RCSSFPoly', 'SigmaZeroSFPoly', 'BetaZeroSFPoly', 'GammaZeroSFPoly'):
                    m = re.search(rf'<{tag}[^>]*>\s*<Coef[^>]*>([^<]+)</Coef>', text)
                    rad[tag] = float(m.group(1)) if m else None
                rad['NoiseLevel_present'] = '<NoiseLevel' in text
                out['sicd_radiometric'] = rad
    return out


def main() -> None:
    from sarsim.dwell import DwellProduct
    from sarsim.ortho import frame_to_lla, multilook, project, surface
    from ..compose import load_site
    from ..export.sites import scene as site_scene
    sc = site_scene(load_site('giza'))
    geoid = float(sc['frame'].get('geoid_m', 15.5))
    for name, cfg in PASSES.items():
        rec_path = RECORDS / f'{name}.json'
        rec = json.loads(rec_path.read_text())
        p = DwellProduct(cfg['product'])
        prod = product_fields(cfg['product'])
        deliv = delivery_fields(cfg['delivery'])
        shift = tuple(json.loads((RESULTS / cfg['offset_from'] / 'summary.json').read_text())['registration']['shift_px'])
        k = float(prod['calibration_factor'])
        with h5py.File(cfg['product'], 'r') as f:
            dx = float(f['azimuth_ground_spacing'][()]) if 'azimuth_ground_spacing' in f else float(rec['dx'])
            dr = float(f['slant_range_spacing'][()])
            lia = np.asarray(f['local_incidence_angle'][()], float)
        theta = np.deg2rad(float(prod['incidence_center']))
        # ---- an upper bound on the noise floor: the darkest 10 m cells over the Giza site
        x0, x1 = sc['extent']['x']
        y0, y1 = sc['extent']['y']
        gx, gy = np.meshgrid(np.linspace(x0, x1, 41), np.linspace(y0, y1, 41))
        gr, gc = project(p, sc, gx, gy, surface(sc, gx, gy), geoid, shift)
        la10, lr10 = max(1, round(DARK_CELL_M / dx)), max(1, round(DARK_CELL_M * np.sin(theta) / dr))
        obs, (r0, c0) = multilook(p, int(gr.min()), int(gr.max()) + 1, int(gc.min()), int(gc.max()) + 1, la10, lr10)
        cols = c0 + lr10 * (np.arange(obs.shape[1]) + 0.5)
        sig = k * obs * np.sin(np.deg2rad(np.interp(cols, np.arange(lia.size), lia)))[None, :]
        good = sig[sig > 0]
        dark = {'cell_m': DARK_CELL_M, 'cells': int(good.size),
                'sigma0_db_min': float(10 * np.log10(good.min())),
                'sigma0_db_p0_1': float(10 * np.log10(np.percentile(good, 0.1))),
                'sigma0_db_p1': float(10 * np.log10(np.percentile(good, 1))),
                'sigma0_db_median': float(10 * np.log10(np.median(good))),
                'note': 'sigma0 = calibration_factor x mean |z|^2 x sin(local incidence), over 10 m cells of the Giza site in '
                        'the image; the darkest cells hold the noise floor plus any signal left, so the noise floor is at '
                        'most their level'}
        rec['radiometry'] = {
            'about': 'Written by sim/katabasis/compose/acquisition_radiometry.py from the product and its delivery XMLs.',
            'product': prod,
            'delivery': deliv,
            'beta0': 'beta0 = calibration_factor x (I^2 + Q^2); sigma0 = beta0 x sin(incidence) (ICEYE); the product sets '
                     'range_spread_comp_flag and leaves antenna_pattern_compensation and fsl_compensation at zero',
            'noise_floor': {'in_product': False,
                            'note': 'No noise-floor field in the product: no NESZ in the HDF5 or the SLC XML, and no '
                                    'Radiometric/NoiseLevel in the SICD XML where the delivery has one.',
                            'specification': SPEC, 'measured_upper_bound': dark},
        }
        rec_path.write_text(json.dumps(rec, indent=1))
        # ---- the power map over Khafre's footprint, in the image's own geometry
        tx = np.arange(-HALF_M, HALF_M + 1e-6, TIE_M)
        TX, TY = np.meshgrid(tx, tx)
        TZ = surface(sc, TX, TY)
        tr, tc = project(p, sc, TX, TY, TZ, geoid, shift)
        la, lr = max(1, round(CELL_M / dx)), max(1, round(CELL_M * np.sin(theta) / dr))
        m = 40
        pw, (pr0, pc0) = multilook(p, int(tr.min()) - m * la // 4, int(tr.max()) + m * la // 4, int(tc.min()) - m,
                                   int(tc.max()) + m, la, lr)
        ccols = pc0 + lr * (np.arange(pw.shape[1]) + 0.5) - 0.5
        inc = np.interp(ccols, np.arange(lia.size), lia)
        lat, lon = frame_to_lla(sc['frame']['origin'], TX, TY)
        with h5py.File(cfg['product'], 'r') as f:
            rpc = {f'rpc_{kk.lower()}': np.asarray(f['RPC/' + kk][()], float) for kk in f['RPC'].keys()}
        out = RECORDS / f'{name}-khafre-power.npz'
        np.savez_compressed(
            out, power_dn2=pw.astype(np.float32), line0=np.int64(pr0), sample0=np.int64(pc0),
            lines_per_cell=np.int64(la), samples_per_cell=np.int64(lr), line_spacing_m=np.float64(dx),
            slant_range_spacing_m=np.float64(dr), incidence_deg=inc.astype(np.float32),
            calibration_factor=np.float64(k),
            tie_x_m=TX, tie_y_m=TY, tie_z_m=TZ, tie_latitude=lat, tie_longitude=lon, tie_height_ellipsoid_m=TZ + geoid,
            tie_line=tr, tie_sample=tc, geoid_m=np.float64(geoid), shift_px=np.asarray(shift, float),
            frame_origin=np.array([sc['frame']['origin']['latitude'], sc['frame']['origin']['longitude']]), **rpc,
            about=np.array('power_dn2[i, j]: the mean of I^2 + Q^2 over product lines line0 + i*lines_per_cell .. +'
                           'lines_per_cell-1 and samples sample0 + j*samples_per_cell .. +samples_per_cell-1. '
                           'beta0 = calibration_factor x power; sigma0 = beta0 x sin(incidence_deg[j]). Tie points: '
                           'site-frame x, y (east, north from Khafre\'s base centre) and z (orthometric, terrain and '
                           'pyramid), their latitude, longitude and ellipsoidal height, and their product line and sample '
                           'through the RPC with shift_px added (the residual offset fitted for this pass).'))
        print(f"{name}: calibration {k:.4g}; darkest 10 m cells {dark['sigma0_db_min']:.1f} dB (0.1% {dark['sigma0_db_p0_1']:.1f}, "
              f"median {dark['sigma0_db_median']:.1f}); power map {pw.shape} cells of {la} x {lr} px -> {out.stat().st_size / 1e3:.0f} kB")


if __name__ == '__main__':
    main()
