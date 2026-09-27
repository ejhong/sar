"""Audit native ICEYE geolocation using small metadata reads only.

Outputs contain selected acquisition metadata and diagnostics, not SAR pixels.
The report is not a validated target registration or a motion/depth result.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sarsim.geolocation import audit_iceye_geometry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', type=Path)
    parser.add_argument('--sicd-xml', type=Path)
    parser.add_argument('--out', type=Path, required=True, help='Small JSON report destination')
    args = parser.parse_args()
    if args.out.resolve() in {args.path.resolve(), args.sicd_xml.resolve() if args.sicd_xml else None}:
        parser.error('Output must not overwrite an input')
    report = audit_iceye_geometry(args.path, args.sicd_xml)
    report['source_sha256'] = {
        name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        for name in ('sarsim/geolocation.py', 'scripts/audit_iceye_geometry.py')
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(f'Report: {args.out}')
    print('Native registration and depth analysis remain unvalidated.')
    print('Centre sensitivity to +10 m height:', report['control_point_checks'][0]['height_plus_10m_shift_px'], 'pixels (line, sample)')


if __name__ == '__main__':
    main()
