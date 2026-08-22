#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/ramesh.kolluru/CFD_Solver_withCUDA
LOG=$ROOT/results/scheme_sweep_retries.log
# wait for orchestrator
while pgrep -f 'run_halfcyl_all_schemes_orchestrator' >/dev/null 2>&1; do
  echo "$(date -Is) waiting for orchestrator..." | tee -a "$LOG"
  sleep 60
done
export LD_LIBRARY_PATH=/home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=12
echo "==== retries start $(date -Is) ====" | tee -a "$LOG"
# Force re-run anything not ok / skipped_have_solution
python3 - <<'PY' 2>&1 | tee -a "$LOG"
import json
from pathlib import Path
st=Path('/home/ramesh.kolluru/CFD_Solver_withCUDA/results/scheme_sweep/status.json')
d=json.loads(st.read_text()) if st.exists() else {}
for k,v in list(d.items()):
    if v.get('state') not in ('ok','skipped_have_solution'):
        # delete so runner will execute
        if v.get('state') in ('fail','retry_pending','ok_external','fail_retry'):
            d[k]={'state':'pending_force'}
st.write_text(json.dumps(d, indent=2)+'\n')
print('pending force:', [k for k,v in d.items() if v.get('state')=='pending_force'])
PY
# Temporarily treat pending_force as not-ok (already), run with --force for failed tags only via custom
python3 - <<'PY' 2>&1 | tee -a "$LOG"
import json, os, subprocess, time
from pathlib import Path
import importlib.util
ROOT=Path('/home/ramesh.kolluru/CFD_Solver_withCUDA')
# monkey: call runner main pieces
spec=importlib.util.spec_from_file_location('sweep', ROOT/'scripts'/'run_halfcyl_scheme_sweep.py')
sweep=importlib.util.module_from_spec(spec); spec.loader.exec_module(sweep)
manifest=json.loads(sweep.MANIFEST.read_text())
status=sweep.load_status()
env=os.environ.copy()
env['LD_LIBRARY_PATH']='/home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu:'+env.get('LD_LIBRARY_PATH','')
env['OMP_NUM_THREADS']='12'
for m in manifest:
    tag=m['tag']
    st=status.get(tag,{}).get('state')
    if st in ('ok','skipped_have_solution'):
        continue
    print(f'[RETRY] {tag} (was {st})', flush=True)
    sweep.warm_start(m['dname'], m['order'])
    log=sweep.LOGDIR/f'{tag}.log'
    t0=time.time()
    with open(log,'w') as lf:
        p=subprocess.run([str(sweep.BIN), m['file']], cwd=str(ROOT/'build-cuda'), stdout=lf, stderr=subprocess.STDOUT, env=env)
    dt=time.time()-t0
    text=log.read_text(errors='replace')
    ok=p.returncode==0 and 'completed successfully' in text
    sweep.archive(tag, m['dname'], m['order'], log)
    status[tag]={'state':'ok' if ok else 'fail', 'returncode':p.returncode, 'seconds':round(dt,1), 'log':str(log), 'retry':True}
    sweep.save_status(status)
    print(f'[{"OK  " if ok else "FAIL"}] {tag}  {dt/60:.1f} min', flush=True)
print('Retries done')
for k,v in sorted(status.items()):
    print(f"  {k:20s} {v.get('state')}")
PY
echo "==== retries done $(date -Is) ====" | tee -a "$LOG"
