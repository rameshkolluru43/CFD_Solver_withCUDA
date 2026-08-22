# 3D Compressible Euler / NS (`src3d/`)

Full **3D compressible** path (host + CUDA), independent of the incomplete `src_3D/` tree.

## Build

```bash
cmake -S . -B build && cmake --build build --target CFD_solver_3d CFD_solver_3d_gpu -j
```

| Binary | Role |
|--------|------|
| `CFD_solver_3d` | Host Euler/NS |
| `CFD_solver_3d_gpu` | CUDA FO resident; hybrid host flux for MUSCL/WENO/viscous |

## Capabilities

| Feature | Status |
|---------|--------|
| Conservatives `[ρ,ρu,ρv,ρw,ρE]` | Yes |
| Cartesian hex mesh | Yes |
| Schemes LLF / RICCA / Roe / Van Leer | Yes |
| RICCA legacy pressure / RICCA-RH flux sensors | Yes |
| WENO hybrid reconstruction (combo / pressure / off) | Yes |
| Order 1O / MUSCL / WENO5 | Yes |
| Laminar viscous NS | Yes |
| BC: slip / noslip / extrapolate / freestream | Yes |
| Cases: Sod-x, freestream, 15° ramp, **flat-plate SWBLI** | Yes |
| CUDA FO multi-scheme | Yes |
| CUDA hybrid (WENO/MUSCL/viscous) | Yes |
| VTK output | Yes |

## Example configs

- `json_Files/Run_3D_Sod_LLF_smoke.json`
- `json_Files/Run_3D_Sod_ROE_smoke.json`
- `json_Files/Run_3D_Sod_WENO_smoke.json`
- `json_Files/Run_3D_Sod_viscous_smoke.json`
- `json_Files/Run_3D_Freestream_smoke.json`
- `json_Files/Run_3D_Sod_LLF_smoke_cuda.json`

```bash
./build/CFD_solver_3d json_Files/Run_3D_Sod_WENO_smoke.json
./build/CFD_solver_3d_gpu json_Files/Run_3D_Sod_LLF_smoke_cuda.json
```

## Still out of scope (vs production 2D)

- Unstructured / GMSH hex meshes
- Characteristic WENO / LES / implicit
- Drop-in reuse of legacy `src_3D/` sources
