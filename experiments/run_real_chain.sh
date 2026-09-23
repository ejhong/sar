#!/bin/bash
# Wait for the R1 compute pass to finish, then report and run R3-R5.
cd "$(dirname "$0")"
PY=../.venv/bin/python
while pgrep -f "r01_giza_controls.py" > /dev/null; do sleep 15; done
echo "=== R1 report ==="
$PY r01_giza_controls.py --report-only 2>&1 | grep -vE "Warning|savefig"
echo "=== R3 split dwell ==="
$PY r03_split_dwell.py 2>&1 | grep -vE "Warning|savefig"
echo "=== R4 velocity floor ==="
$PY r04_velocity_floor.py 2>&1 | grep -vE "Warning|savefig"
echo "=== R5 learned null ==="
$PY r05_learned_null.py 2>&1 | grep -vE "Warning|savefig"
echo "=== CHAIN DONE ==="
