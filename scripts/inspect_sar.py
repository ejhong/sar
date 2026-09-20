"""Inspect a huge local SAR file without loading it; optionally read a bounded chip.

python scripts/inspect_sar.py /path/to/scene.h5 --out data/inspection
python scripts/inspect_sar.py /path/to/scene.h5 --rows 1000 1512 --cols 2000 2512 --out data/inspection
Outputs stay in ignored data/ by default and are never picked up by the site build.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sarsim.realdata import inspect_product, read_crop


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', type=Path)
    parser.add_argument('--rows', type=int, nargs=2)
    parser.add_argument('--cols', type=int, nargs=2)
    parser.add_argument('--out', type=Path, default=ROOT/'data'/'inspection')
    args = parser.parse_args()
    if bool(args.rows) != bool(args.cols):
        parser.error('provide both --rows and --cols, or neither')
    args.out.mkdir(parents=True, exist_ok=True)
    if args.rows:
        crop, report = read_crop(args.path, args.rows, args.cols)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
        reference = max(float(np.abs(crop).max()), 1e-30)
        axes[0].imshow(20*np.log10(np.abs(crop)/reference+1e-6), cmap='gray', vmin=-50, vmax=0, origin='lower')
        axes[0].set(title='Bounded native-image preview', xlabel='Crop column', ylabel='Crop row')
        spectral = np.fft.fftshift(np.abs(np.fft.fft2(crop))**2)
        axes[1].imshow(10*np.log10(spectral/max(float(spectral.max()), 1e-30)+1e-10), cmap='magma', vmin=-50, vmax=0, origin='lower')
        axes[1].set(title='Crop spectrum — affected by crop window', xlabel='Column frequency bin', ylabel='Row frequency bin')
        fig.savefig(args.out/'preview.png', dpi=150)
        plt.close(fig)
        np.save(args.out/'complex_crop.npy', crop)
        report['preview_scope'] = 'Image/spectrum QC only; neither axis is a validated aperture-time or depth coordinate.'
    else:
        report = inspect_product(args.path)
    (args.out/'inspection.json').write_text(json.dumps(report, indent=2, allow_nan=False, default=str)+'\n')
    print(json.dumps({k: report[k] for k in ('filename', 'file_bytes', 'format', 'complex_samples', 'depth_analysis_ready')}, indent=2))
    print(f'Report: {args.out / "inspection.json"}')


if __name__ == '__main__':
    main()
