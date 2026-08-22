# 3D Compressible Euler / NS (`solvers/src3d/`)

Independent of the incomplete `solvers/src_3D/` tree. Conservatives `[ρ, ρu, ρv, ρw, ρE]` on a Cartesian hex mesh.

## Build

```bash
cmake -S . -B build && cmake --build build --target CFD_solver_3d CFD_solver_3d_gpu -j
```

| Binary | Role |
|--------|------|
| `CFD_solver_3d` | Host Euler/NS/LES |
| `CFD_solver_3d_gpu` | CUDA first-order resident; MUSCL/WENO/viscous hybrid on host |

Run **from the repository root** so `VTK_Out` under `solutions/` is correct.

## Capabilities

| Feature | Status |
|---------|--------|
| Schemes LLF / RICCA / Roe / Van Leer | Yes |
| RICCA legacy pressure / RICCA-RH sensors | Yes |
| WENO-Z, characteristic, hybrid combo | Yes |
| Order 1O / MUSCL / WENO5 | Yes |
| Laminar NS, WALE LES, trip (`LES_Perturb`) | Yes |
| Inactive/solid cells (`Mesh.active`) | Yes (cavity block mesh) |
| BC: slip / noslip / extrapolate / freestream | Yes |
| Cases: Sod-x, freestream, 15° ramp, SWBLI plate, Mach-6 cavity | Yes |
| VTK output | Yes (legacy VTK text) |

## Example configs (`input/json_Files/`)

Smoke: `Run_3D_Sod_LLF_smoke.json`, `Run_3D_Sod_LLF_smoke_cuda.json`, `Run_3D_Cavity_M6_LES_WALE_smoke.json`, `Run_3D_SWBLI_LES_WALE_smoke.json`

Production-length: `Run_3D_Cavity_M6_LES_WALE.json`, `Run_3D_SWBLI_LES_WALE.json`, ramp `Run_3D_Ramp15_finest_Lx4p5_*`

```bash
./build/CFD_solver_3d input/json_Files/Run_3D_Sod_WENO_smoke.json
./build/CFD_solver_3d_gpu input/json_Files/Run_3D_Sod_LLF_smoke_cuda.json
```

## Out of scope (vs 2D)

- Unstructured / GMSH hex meshes
- Implicit time marching
- Reuse of legacy `src_3D/` sources
