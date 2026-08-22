# Build and Run Guide

## Requirements

- **CMake** ≥ 3.16
- **C++17** compiler (GCC, Clang, or MSVC)
- **CUDA Toolkit** ≥ 11.0 (for `CFD_solver_gpu`)
- **MPI** (optional; for `CFD_solver_mpi`)
- **Boost** (regex)
- **JsonCpp**
- **VTK** ≥ 9.4 (**required** for the default `CFD_solver` CMake target)

## Build

```bash
cd CFD_Solver_withCUDA
cmake -S . -B build
cmake --build build -j$(nproc)

# GPU-focused out-of-tree build
cmake -S . -B build-cuda
cmake --build build-cuda --target CFD_solver_gpu -j$(nproc)
```

Targets when dependencies are found:

| Target | Needs |
|--------|--------|
| `CFD_solver` | VTK + JsonCpp + Boost |
| `CFD_solver_gpu` | CUDA (`nvcc`) + above |
| `CFD_solver_mpi` | VTK + MPI |

Set `CUDA_ARCHITECTURES` in `CMakeLists.txt` to match your GPU (e.g. `75` Turing, `80`/`86` Ampere).

### GPU flux schemes (MOVERS / RICCA / LLF)

On **`CFD_solver_gpu`**, first-order net flux uses:

- `Dissipation_Type` **2 / 4 / 5** → `CUDA_KERNELS/MOVERS_RICCA_Flux_Cuda_Kernels.cu`
- `Dissipation_Type` **1** (LLF) → generic 1O CUDA main-loop path
- Roe / MUSCL 2O → host unless separately wired

**WENO + RICCA** for the validated half-cylinder path runs on the **host** (GPU WENO kept off until bit-identical). See [MAIN_LOOP_STATUS.md](MAIN_LOOP_STATUS.md) and [overview.md](../overview.md).

When resident GPU init succeeds, `Solver.cpp` can take `Resident_GPU_Explicit_Step` for inviscid/NS; otherwise it falls back to the host loop.

### GPU viscous flux and viscous wall BC

When **`Is_Viscous`** / **`Is_Viscous_Wall`** is true:

- CUDA viscous integration: `CUDA_KERNELS/Viscous_Flux_Cuda_Integration.cu`
- CUDA no-slip wall ghosts: `CUDA_KERNELS/Boundary_Conditions_Cuda_Integration.cu`

Validated Mach-6 NS continues after the corner fix are documented as **host WENO+RICCA + host viscous**. Set Re / Tw in the test-case JSON (e.g. Re = 1e5, Tw = 2).

## Run

```bash
./build/CFD_solver input/json_Files/Solver_Config.json
mpirun -np 4 ./build/CFD_solver_mpi input/json_Files/Solver_Config.json
./build-cuda/CFD_solver_gpu input/json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json
```

### Example configs

| Config | Purpose |
|--------|---------|
| `input/json_Files/Test_Config_WENO.json` | Generic WENO |
| `input/json_Files/Test_Config_AMR.json` | Gradient AMR |
| `input/json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json` | M=6 viscous 500-iter smoke (P3 continue) |
| `input/json_Files/Run_HalfCylinder_M6_WENO_RICCA_P3_viscous_20k.json` | Longer viscous continue |
| `input/json_Files/Half_Cylinder_481_161_M6_viscous.json` | Case: Re, Tw, freestream (via `Test_Case_Json`) |

Full half-cylinder recipe: [HALF_CYLINDER_VALIDATION.md](HALF_CYLINDER_VALIDATION.md).

## Output

- Solution / error files under the case output directory
- VTK for ParaView (periodic and/or final)
- Viscous corner dump: `viscous_corner_diagnostics.txt` in the solution dir when viscous
- Plots: `python3 scripts/plot_halfcyl_inviscid_viscous.py` → `plots/half_cylinder_inviscid_viscous/`

## Documentation (Doxygen)

```bash
./scripts/update_docs.sh
# or:
doxygen docs/Doxyfile_Cleaned
```

Open `docs/doxygen/html/index.html` (or `./scripts/view_docs.sh`).

## Troubleshooting

| Symptom | Likely cause |
|---------|----------------|
| Min_dt is zero | CFL / IC / BC; try smaller CFL |
| Ghost memory not created | `Total_No_Cells` / `Cells` resize before ghost construct |
| Viscous blow-up near corners (Pmax ≫ 50, M ≫ 6) | Face normals not LBRT-aligned — use current `Construct_Cell` / no atan2 reorder on structured TXT |
| CUDA / nvcc | Match `CUDA_ARCHITECTURES` to the GPU |
| `doxygen: command not found` | Install Doxygen, or place a binary on `PATH` |
