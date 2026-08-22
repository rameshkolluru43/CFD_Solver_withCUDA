# Solvers

| Path | Role |
|------|------|
| `src/` | 2D compressible Euler/NS. Entry: `Main.cpp`. OpenMP; MPI via `CFD_solver_mpi`. |
| `src3d/` | 3D compressible Euler/NS/LES. Binaries `CFD_solver_3d`, `CFD_solver_3d_gpu`. |
| `CUDA_KERNELS/` | CUDA kernels and host wrappers (2D + `Euler3D_Cuda.cu`). |
| `include/` | Shared 2D headers (`definitions.h`, `Solver.h`, …). |
| `src_3D/` | Incomplete legacy 3D tree; do not use for new work. |
| `Metal_Kernels/` | Optional Apple Metal path. |
| `MPI_Src/` | Optional MPI numerical-method bits. |
| `Incompressible_Solver/` | Separate incompressible / SIMPLE sources. |
| `Laplace_Solver/` | Laplace demo. |
| `Parallel_Method_Kernels/` | Host parallel flux experiments. |
| `Basic_Function_Files/` | Older cell/flux helpers. |
| `Main/`, `Flow_Over_Cylinder/` | Legacy makefiles / cylinder drivers. |

CMake reads `solvers/src`, `solvers/src3d`, `solvers/include`, and `solvers/CUDA_KERNELS`. See [LAYOUT.md](../LAYOUT.md).
