"""Rerun the experiments and the generated page.

By default only the simulation chapter runs, because it is self-contained. The real chapter
needs the two ICEYE products named in `real_common.PRODUCTS` and takes roughly an hour on a
laptop; pass --real to include it, or --real-only to rebuild just that chapter.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"


def run(script):
    print(f"Running {script.name}", flush=True)
    subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)


def simulation():
    # The staged-column illustration remains archived and is no longer on the page.
    scripts = sorted(p for p in EXPERIMENTS.glob("t[0-9][0-9]_*.py") if p.name != "t07_wells.py")
    return scripts + [EXPERIMENTS / n for n in ("feasibility.py", "robustness.py", "validation.py")]


def real():
    return sorted(EXPERIMENTS.glob("r[0-9][0-9]_*.py"))


def main(argv):
    scripts = []
    if "--real-only" not in argv:
        scripts += simulation()
    if "--real" in argv or "--real-only" in argv:
        missing = [name for name in ("r01_giza_controls",) if not (EXPERIMENTS / f"{name}.py").is_file()]
        if missing:
            raise SystemExit(f"missing real-chapter scripts: {missing}")
        scripts += real()
    for script in scripts:
        run(script)
    subprocess.run([sys.executable, str(ROOT / "build_site.py")], cwd=ROOT, check=True)


if __name__ == "__main__":
    main(sys.argv[1:])
