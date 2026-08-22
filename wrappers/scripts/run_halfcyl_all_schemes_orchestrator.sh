#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/ramesh.kolluru/CFD_Solver_withCUDA
export LD_LIBRARY_PATH=/home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
LOG=$ROOT/solutions/results/scheme_sweep_orchestrator.log
exec > >(tee -a "$LOG") 2>&1
echo "==== orchestrator start $(date -Is) ===="

echo "--- GPU-friendly schemes (1O/WENO MOVERS/RICCA/NWSC/RICCA_LLF) ---"
OMP_NUM_THREADS=4 python3 "$ROOT/wrappers/scripts/run_halfcyl_scheme_sweep.py" gpu

echo "--- waiting for any existing RICCA_2O CFD_solver_gpu to finish ---"
while pgrep -f 'CFD_solver_gpu .*RICCA_2O_converge' >/dev/null 2>&1; do
  echo "$(date -Is) RICCA_2O still running..."
  sleep 120
done
# also wait if generic 2O converge still up
while pgrep -f 'Run_HalfCylinder_M6_RICCA_2O_converge' >/dev/null 2>&1; do
  echo "$(date -Is) RICCA_2O converge still running..."
  sleep 120
done

# archive RICCA_2O if present
OUT=$ROOT/solutions/2D_Euler_Solutions/Half_Cylinder/Explicit/Flux_2/GridSize_7
DEST=$ROOT/solutions/results/scheme_sweep/RICCA_2O
mkdir -p "$DEST"
cp -f "$OUT"/Solution_RICCA_2O_MinMod.txt "$OUT"/Error_RICCA_2O_MinMod.txt \
  "$OUT"/results_RICCA_2O_MinMod.txt "$DEST"/ 2>/dev/null || true
cp -f "$OUT"/Half_Cylinder_Grid_Size_7_FluxScheme_2O_MinMod_Dissipation_RICCA.vtk "$DEST"/ 2>/dev/null || true
tail -40 "$ROOT/results/halfcyl_RICCA_2O_converge.log" > "$DEST/tail.log" 2>/dev/null || true
python3 - <<'PY'
import json
from pathlib import Path
st=Path('/home/ramesh.kolluru/CFD_Solver_withCUDA/results/scheme_sweep/status.json')
d=json.loads(st.read_text()) if st.exists() else {}
d['RICCA_2O']={'state':'ok_external','log':'results/halfcyl_RICCA_2O_converge.log'}
st.write_text(json.dumps(d, indent=2)+'\n')
PY

echo "--- CPU / host schemes (LLF/ROE + all 2O + remaining) ---"
OMP_NUM_THREADS=12 python3 "$ROOT/scripts/run_halfcyl_scheme_sweep.py" all

echo "==== orchestrator done $(date -Is) ===="
python3 - <<'PY'
import json
from pathlib import Path
st=json.loads(Path('/home/ramesh.kolluru/CFD_Solver_withCUDA/results/scheme_sweep/status.json').read_text())
ok=sum(1 for v in st.values() if str(v.get('state','')).startswith('ok') or v.get('state')=='skipped_have_solution')
fail=sum(1 for v in st.values() if v.get('state')=='fail')
print(f'SUMMARY: {ok} ok/skipped, {fail} fail, {len(st)} tracked')
for k,v in sorted(st.items()):
    print(f"  {k:20s} {v.get('state')}")
PY
