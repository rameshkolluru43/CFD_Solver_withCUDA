#!/usr/bin/env bash
export LD_LIBRARY_PATH=/home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
exec python3 /home/ramesh.kolluru/CFD_Solver_withCUDA/scripts/run_halfcyl_scheme_sweep.py "$@"
