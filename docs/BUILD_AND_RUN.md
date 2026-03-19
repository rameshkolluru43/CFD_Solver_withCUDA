# Build and Run Guide

## Requirements

- **CMake** ≥ 3.16  
- **C++17** compiler (GCC, Clang, or MSVC)  
- **CUDA Toolkit** ≥ 11.0 (for GPU build)  
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

Two executables are produced:

- **CFD_solver** — CPU-only solver  
- **CFD_solver_gpu** — Solver with CUDA kernels  

## Run

Pass the path to a **JSON configuration file**:

```bash
# CPU
./CFD_solver ../json_Files/Solver_Config.json

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
