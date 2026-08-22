#!/usr/bin/env python3
"""Run half-cylinder M6 with all dissipation × reconstruction schemes."""
from __future__ import annotations
import json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path('/home/ramesh.kolluru/CFD_Solver_withCUDA')
BIN = ROOT / 'build-cuda' / 'CFD_solver_gpu'
OUT = ROOT / 'solutions/2D_Euler_Solutions/Half_Cylinder/Explicit/Flux_2/GridSize_7'
ARCH = ROOT / 'solutions' / 'results' / 'scheme_sweep'
LOGDIR = ROOT / 'solutions' / 'results' / 'scheme_sweep_logs'
MANIFEST = ROOT / 'input' / 'json_Files' / 'scheme_sweep' / 'manifest.json'
STATUS = ARCH / 'status.json'
WARM1O = OUT / 'Solution_RICCA_1O_MinMod.txt'

def load_status():
    if STATUS.exists():
        return json.loads(STATUS.read_text())
    return {}

def save_status(st):
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(json.dumps(st, indent=2) + '\n')

def warm_start(dname: str, order: str):
    if order == 'WENO':
        wdir = OUT / 'WENO'
        wdir.mkdir(exist_ok=True)
        # Non-RICCA WENO is unstable if warm-started from a RICCA-WENO field;
        # always seed from robust RICCA 1O. RICCA/RICCA_LLF may use prior WENO.
        if dname in ('RICCA', 'RICCA_LLF'):
            warmw = wdir / 'Solution_WENO_Grid_Size_7.txt'
            src = warmw if warmw.exists() and warmw.stat().st_size > 1000 else WARM1O
        else:
            src = WARM1O
        (wdir / 'Solution_WENO_Grid_Size_7.txt').write_bytes(Path(src).read_bytes())
    else:
        dest = OUT / f'Solution_{dname}_{order}_MinMod.txt'
        dest.write_bytes(WARM1O.read_bytes())

def archive(tag: str, dname: str, order: str, log: Path):
    dest = ARCH / tag
    dest.mkdir(parents=True, exist_ok=True)
    if order == 'WENO':
        for name in ['Solution_WENO_Grid_Size_7.txt', 'Error_WENO_Grid_Size_7.txt', 'results_WENO_Grid_Size_7.txt']:
            f = OUT / 'WENO' / name
            if f.exists():
                (dest / name).write_bytes(f.read_bytes())
        for f in (OUT / 'WENO').glob(f'*Dissipation_{dname}*.vtk'):
            (dest / f.name).write_bytes(f.read_bytes())
    else:
        for name in [f'Solution_{dname}_{order}_MinMod.txt', f'Error_{dname}_{order}_MinMod.txt', f'results_{dname}_{order}_MinMod.txt']:
            f = OUT / name
            if f.exists():
                (dest / name).write_bytes(f.read_bytes())
        for f in OUT.glob(f'*FluxScheme_{order}_*Dissipation_{dname}.vtk'):
            (dest / f.name).write_bytes(f.read_bytes())
    dest.joinpath('tail.log').write_text('\n'.join(log.read_text(errors='replace').splitlines()[-40:]) + '\n')

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'  # all|gpu|cpu
    force = '--force' in sys.argv
    manifest = json.loads(MANIFEST.read_text())
    status = load_status()
    LOGDIR.mkdir(parents=True, exist_ok=True)
    ARCH.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = '/home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu:' + env.get('LD_LIBRARY_PATH', '')
    env.setdefault('OMP_NUM_THREADS', '8')

    # GPU-resident jobs first so they can overlap with any host 2O run
    items = list(manifest)
    if mode == 'all':
        items.sort(key=lambda m: (0 if m.get('gpu_ish') else 1, m['tag']))
    for m in items:
        tag = m['tag']
        if mode == 'gpu' and not m['gpu_ish']:
            continue
        if mode == 'cpu' and m['gpu_ish']:
            continue
        if m.get('skip_reason') == 'already_running' and not force:
            print(f'[SKIP] {tag}: already_running', flush=True)
            status[tag] = {**status.get(tag, {}), 'state': 'skipped_running'}
            save_status(status)
            continue
        if m.get('skip_reason') == 'already_have_solution' and not force:
            print(f'[SKIP] {tag}: already_have_solution — archiving', flush=True)
            # ensure archive copy of existing artifacts
            fake_log = LOGDIR / f'{tag}.log'
            if not fake_log.exists():
                fake_log.write_text('preexisting solution; not re-run\n')
            archive(tag, m['dname'], m['order'], fake_log)
            status[tag] = {**status.get(tag, {}), 'state': 'skipped_have_solution'}
            save_status(status)
            continue
        if status.get(tag, {}).get('state') == 'ok' and not force:
            print(f'[SKIP] {tag}: already ok', flush=True)
            continue
        if status.get(tag, {}).get('state') == 'fail' and not force and '--no-retry' in sys.argv:
            print(f'[SKIP] {tag}: fail (no-retry)', flush=True)
            continue

        print(f'[RUN ] {tag}', flush=True)
        warm_start(m['dname'], m['order'])
        log = LOGDIR / f'{tag}.log'
        t0 = time.time()
        with open(log, 'w') as lf:
            p = subprocess.run([str(BIN), m['file']], cwd=str(ROOT / 'build-cuda'),
                               stdout=lf, stderr=subprocess.STDOUT, env=env)
        dt = time.time() - t0
        text = log.read_text(errors='replace')
        ok = p.returncode == 0 and 'completed successfully' in text
        archive(tag, m['dname'], m['order'], log)
        status[tag] = {'state': 'ok' if ok else 'fail', 'returncode': p.returncode,
                       'seconds': round(dt, 1), 'log': str(log)}
        save_status(status)
        print(f'[{"OK  " if ok else "FAIL"}] {tag}  {dt/60:.1f} min  rc={p.returncode}', flush=True)

    print('Sweep finished.', flush=True)
    # summary
    for k, v in sorted(status.items()):
        print(f"  {k:20s} {v.get('state')}  {v.get('seconds','-')}s", flush=True)

if __name__ == '__main__':
    main()
