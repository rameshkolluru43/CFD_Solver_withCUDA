# Build and Run Guide

## Requirements

- **CMake** ≥ 3.16  
- **C++17** compiler (GCC, Clang, or MSVC)  
- **CUDA Toolkit** ≥ 11.0 (for GPU build)  
- **MPI** implementation such as OpenMPI or MPICH (for MPI CPU build)  
- **Boost** (regex)  
- **JsonCpp** (JSON config parsing)  
- **VTK** ≥ 9.4 (optional; for some visualization/IO)

## Build

```bash
cd CFD_Solver_withCUDA
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

When dependencies are found, CMake may produce:

- **CFD_solver** — CPU solver (requires VTK in CMake)  
- **CFD_solver_gpu** — Solver built with `USE_CUDA=1` (requires CUDA toolkit)  
- **CFD_solver_mpi** — CPU solver with MPI rank-partitioned cell loops (requires VTK and MPI)

### GPU flux schemes (MOVERS / RICCA / LLF)

On **`CFD_solver_gpu`**, first-order net flux tries the MOVERS/RICCA CUDA kernels for **`Dissipation_Type`** **2** (MOVERS), **4** (RICCA), and **5** (MOVERS_NWSC) via `CUDA_KERNELS/MOVERS_RICCA_Flux_Cuda_Kernels.cu`, then falls back to the generic 1st-order inviscid CUDA path for **1** (LLF). Second-order runs (`Is_Second_Order` with MUSCL) still use the CPU path in `src/Net_Flux.cpp`. Other dissipation types (Roe, etc.) remain CPU unless separate CUDA kernels are wired in.

### GPU viscous flux and viscous wall BC

When **`Is_Viscous_Wall`** is true (viscous Navier–Stokes):

- **`Evaluate_Viscous_Fluxes()`** uses `CUDA_KERNELS/Viscous_Flux_Cuda_Integration.cu` (Green–Gauss cell gradients + face stress; up to 8 faces/cell).
- **`Apply_Boundary_Conditions()`** applies no-slip wall state on ghost cells via `CUDA_KERNELS/Boundary_Conditions_Cuda_Integration.cu`, then refreshes ghost primitives on the host.

Inlet/exit/symmetry boundaries still run on the CPU. Set **`Is_Viscous_Wall": true`** and a finite **`Re`** in your test JSON to exercise the viscous CUDA path.

## Run

Pass the path to a **JSON configuration file**:

```bash
# CPU
./CFD_solver ../json_Files/Solver_Config.json

# MPI CPU
mpirun -np 4 ./CFD_solver_mpi ../json_Files/Solver_Config.json

# GPU
./CFD_solver_gpu ../json_Files/Solver_Config.json
```

### Example configs

- **Default / WENO**: `../json_Files/Test_Config_WENO.json`  
- **AMR (gradient-based tagging)**: `../json_Files/Test_Config_AMR.json`  

The JSON file sets the test case, grid, solver options, and (optionally) AMR parameters. The test case JSON (e.g. `Half_Cylinder.json`) defines the grid path and boundary data.

## Output

- **Solution / error files**: Written under the test-case output directory (e.g. `../2D_Euler_Solutions/...`).  
- **VTK**: When the run writes solution VTK (e.g. at intervals or at end), open in ParaView.  
- **Console**: Iteration, dt, errors, and (if AMR is on) periodic “AMR: N cells tagged for refinement”.

## Troubleshooting

- **“Min_dt is zero”**: Usually CFL or initial state; try smaller CFL or check initial/boundary conditions.  
- **“Memory for Ghost Cell is not created”**: Ghost list and `Cells` size mismatch; ensure grid pipeline (e.g. `Form_Cells` or VTK path) sets `Total_No_Cells` and resizes `Cells` before constructing ghost cells.  
- **CUDA / nvcc**: Ensure `CUDA_ARCHITECTURES` in CMake matches your GPU (e.g. 75 for Turing, 86 for Ampere).
