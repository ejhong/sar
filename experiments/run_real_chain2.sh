#!/bin/bash
# Second stage: after the first chain finishes, run the second site, the positive control,
# and rebuild the page.
cd "$(dirname "$0")"
PY=../.venv/bin/python
while pgrep -f "run_real_chain.sh" > /dev/null; do sleep 20; done
echo "=== R6 second site ==="
$PY r06_second_site.py 2>&1 | grep -vE "Warning|savefig"
echo "=== R7 injected motion ==="
$PY r07_injected_motion.py 2>&1 | grep -vE "Warning|savefig"
echo "=== rebuild site ==="
cd .. && .venv/bin/python build_site.py 2>&1 | tail -3
echo "=== CHAIN2 DONE ==="
