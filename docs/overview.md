# CFD Solver with CUDA — Overview

High-level map of the compressible CFD solver, current GPU/CPU status, and the Mach-6 half-cylinder validation path (inviscid and viscous).

## What this code does

Finite-volume solver for 2D compressible **Euler** and **Navier–Stokes** on structured/unstructured meshes (quad-focused for WENO and viscous co-volumes). Time marching is explicit (local time stepping or RK). Flux options include LLF, MOVERS, RICCA, Roe, Van Leer, AUSM, with optional WENO reconstruction.

Executables (CMake):

| Target | Role |
|--------|------|
| `CFD_solver` | 2D CPU (VTK library optional) |
| `CFD_solver_mpi` | 2D MPI |
| `CFD_solver_gpu` | 2D CUDA |
| `CFD_solver_3d` / `CFD_solver_3d_gpu` | 3D host / CUDA |

Build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target CFD_solver_gpu CFD_solver_3d_gpu -j$(nproc)
./build/CFD_solver_3d_gpu input/json_Files/Run_3D_Sod_LLF_smoke_cuda.json
```

2D runs: `cd build` then `./CFD_solver_gpu ../input/json_Files/<run>.json`. Map: [LAYOUT.md](../LAYOUT.md).

Details: [docs/BUILD_AND_RUN.md](docs/BUILD_AND_RUN.md), [docs/CONFIGURATION.md](docs/CONFIGURATION.md).

## Architecture (short)

```
JSON config → mesh (TXT/VTK/GMSH) → Form_Cells / co-volumes
    → IC / restart → BC ghosts each step
    → convective flux (host WENO+RICCA or CUDA 1O schemes)
    → viscous flux (NS; co-volume Green–Gauss face gradients)
    → update / residuals → VTK + solution files
```

Face indexing for structured text grids is **LBRT**: neighbours `[left, bottom, right, top]` with verts `[o,a,b,c]`. Boundary lists store `(cell, face, ghost)`.

## CUDA vs host today

| Path | Status |
|------|--------|
| 1O inviscid net flux (LLF / MOVERS / RICCA / MOVERS_NWSC) | CUDA main-loop |
| WENO reconstruction + RICCA | **Host** (GPU WENO kept off until bit-identical to host) |
| Viscous fluxes | Host (CUDA viscous integration exists; NS runs validated on host) |
| RK / implicit / AMR topology | Host |

See [docs/MAIN_LOOP_STATUS.md](docs/MAIN_LOOP_STATUS.md).

## Half-cylinder M∞ = 6 (validation focus)

- Grid: `input/Grid_Files/Half_Cylinder/Half_Cylinder_481_161.txt` (480×160 cells)
- Inviscid reference: WENO + RICCA, **Pmax ≈ 48.3**, **Mmax ≈ 6**  
  Fields: `solutions/2D_Euler_Solutions/Half_Cylinder/Explicit/Flux_2/GridSize_7/`
- Viscous: isothermal no-slip wall, Re = 1e5, Tw = 2  
  Case JSON: `input/json_Files/Half_Cylinder_481_161_M6_viscous.json`  
  Continue dump: `solutions/results/halfcyl_P3_viscous_20k/`
- Plots: `solutions/plots/half_cylinder_inviscid_viscous/`  
  Regenerate: `python3 wrappers/scripts/plot_halfcyl_inviscid_viscous.py`

### Corner / viscous face-normal fix (2026-07)

**Symptom:** Viscous continues from P3 blew up (Pmax ~137, Mmax ~64, tiny dt).

**Cause:** `Sort_Points_AntiClockWise` plus a geometric face remap misaligned `Face_Normals` with Neighbours/BC indices on ~47 wall cells near the side corners. Wall BC face `1` no longer pointed into the cylinder.

**Fix:**

1. Keep file vertex order `[o,a,b,c]` when reading structured TXT grids.
2. Build face areas/normals in classic LBRT order; flip if not outward.
3. Ghost face midpoints use the same LBRT vertex pairs.
4. Corner secondary neighbours for co-volume gradients corrected; safer face-gradient stencil.
5. Optional dump: `viscous_corner_diagnostics.txt` (written under the solution directory when viscous).

Related robustness fixes in the same campaign:

- Supersonic inlet: prescribe freestream (do not mutate inlet condition state).
- Host WENO wall faces: fall back to 1O when the stencil is degenerate.
- Update(): reject non-positive updates instead of flooring to tiny pressure (avoids fake M ≫ 6).
- Remove double first-order flux evaluation before viscous explicit step.

## Useful run configs

| JSON | Purpose |
|------|---------|
| `json_Files/Half_Cylinder_481_161_M6.json` | Inviscid freestream case params |
| `json_Files/Half_Cylinder_481_161_M6_viscous.json` | Viscous wall / Re / Tw |
| `json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json` | P3 → viscous 500-iter smoke |
| `json_Files/Run_HalfCylinder_M6_WENO_RICCA_P3_viscous_20k.json` | Longer viscous continue |

## Plots and diagnostics

| Artifact | Path |
|----------|------|
| Inviscid / viscous Mach & pressure | `plots/half_cylinder_inviscid_viscous/*.png` |
| Side-by-side compare | `plots/half_cylinder_inviscid_viscous/compare_mach_pressure.png` |
| Corner face/normal dump | `plots/half_cylinder_inviscid_viscous/viscous_corner_diagnostics.txt` |
| WENO order-of-accuracy script | `scripts/test_weno_order_accuracy.py` |

## Documentation index

- [README.md](README.md) — features, build, usage
- [docs/HALF_CYLINDER_VALIDATION.md](docs/HALF_CYLINDER_VALIDATION.md) — M=6 validation recipe
- [docs/README.md](docs/README.md) — full doc table of contents
- [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md) — changelog
- API: `./scripts/update_docs.sh` → `docs/doxygen/html/index.html`

## Contact

Ramesh Kolluru — rameshkolluru43@gmail.com  
Repository: https://github.com/rameshkolluru43/CFD_Solver_withCUDA
