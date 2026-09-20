"""Rerun the numbered experiments, supporting controls, and the generated page."""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    scripts = sorted(path for path in (ROOT / "experiments").glob("t[0-9][0-9]_*.py")
                     if path.name != "t07_wells.py")
    # The staged-column illustration remains archived and is no longer on the page.
    scripts += [ROOT / "experiments" / name for name in ("feasibility.py", "robustness.py", "validation.py")]
    for script in scripts:
        print(f"Running {script.name}", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    subprocess.run([sys.executable, str(ROOT / "build_site.py")], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
