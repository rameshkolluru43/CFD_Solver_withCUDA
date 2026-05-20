# Standalone and Low-End CUDA Run Guide

This guide explains how to build and run the solver on a standalone Linux machine, including machines with older or lower-end NVIDIA GPUs.

## 1. Install System Dependencies

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install -y build-essential cmake git pkg-config \
    libboost-regex-dev libjsoncpp-dev libvtk9-dev \
    libopenmpi-dev openmpi-bin
```

For CUDA runs, install an NVIDIA driver and CUDA Toolkit version supported by your GPU:

```bash
nvidia-smi
nvcc --version
```

If `nvidia-smi` works but `nvcc` is missing, the driver is installed but the CUDA Toolkit is not.

## 2. Clone and Configure

```bash
git clone <repo-url> CFD_Solver_withCUDA
cd CFD_Solver_withCUDA
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

This builds the CPU executable when VTK/JsonCpp/Boost are available. If CUDA is available, it also configures `CFD_solver_gpu`. If MPI is available, it configures `CFD_solver_mpi`.

## 3. Build Targets

Build everything:

```bash
cmake --build build -j$(nproc)
```

Or build specific targets:

```bash
cmake --build build --target CFD_solver -j$(nproc)
cmake --build build --target CFD_solver_mpi -j$(nproc)
cmake --build build --target CFD_solver_gpu -j$(nproc)
```

## 4. Run on CPU

From the repository root:

```bash
./build/CFD_solver json_Files/Test_Config.json
```

For a short smoke test:

```bash
./build/CFD_solver json_Files/MPI_Smoke_Test.json
```

Outputs are written under `2D_Euler_Solutions/<test-case>/...`.

## 5. Run with MPI on CPU

```bash
mpirun -np 4 ./build/CFD_solver_mpi json_Files/Test_Config.json
```

On small machines or containers, OpenMPI may report not enough slots. Use fewer ranks or oversubscribe:

```bash
mpirun -np 2 --oversubscribe ./build/CFD_solver_mpi json_Files/MPI_Smoke_Test.json
```

The current MPI implementation uses a replicated mesh with rank-owned cell ranges. It parallelizes solver loops but keeps full mesh/state data on every rank.

## 6. Configure for a Low-End CUDA GPU

The GPU target uses the CMake cache variable `CFD_CUDA_ARCHITECTURES`. Default is `75` for Turing-class GPUs. For older cards, configure with the architecture matching your GPU:

```bash
cmake -S . -B build-cuda -DCMAKE_BUILD_TYPE=Release -DCFD_CUDA_ARCHITECTURES=61
cmake --build build-cuda --target CFD_solver_gpu -j$(nproc)
```

Common NVIDIA architecture values:

| GPU family | Example cards | CMake value |
|------------|---------------|-------------|
| Maxwell | GTX 750 Ti, GTX 9xx | `50` or `52` |
| Pascal | GTX 10xx, Quadro P-series | `61` |
| Volta | Tesla V100 | `70` |
| Turing | GTX 16xx, RTX 20xx, T4 | `75` |
| Ampere | RTX 30xx, A10, A100 | `86` or `80` |
| Ada | RTX 40xx, L4 | `89` |

If your CMake/CUDA version supports it, you can also try:

```bash
cmake -S . -B build-cuda -DCFD_CUDA_ARCHITECTURES=native
```

## 7. Run on CUDA

```bash
./build-cuda/CFD_solver_gpu json_Files/Test_Config.json
```

For lower-end GPUs, start with a small or short case:

```bash
./build-cuda/CFD_solver_gpu json_Files/MPI_Smoke_Test.json
```

Recommended settings for low-end GPUs:

- Use small grids first, such as the `61 x 21` Half Cylinder setup.
- Keep `Total_Iterations` low for smoke tests, for example `50` or `100`.
- Prefer first-order or MUSCL before WENO when debugging performance.
- Use `Dissipation_Type=1` for LLF or `Dissipation_Type=4` for RICCA. The hybrid WENO RICCA/LLF mode is `Dissipation_Type=6`.
- Avoid running many MPI ranks and CUDA at the same time on a small machine unless you know the memory headroom.

## 8. What CUDA Accelerates Today

`CFD_solver_gpu` is not a full device-resident solver yet. It uses the normal solver loop and currently offloads the 1st-order inviscid net flux path for:

- `Dissipation_Type=1` LLF
- `Dissipation_Type=2` MOVERS
- `Dissipation_Type=4` RICCA
- `Dissipation_Type=5` MOVERS_NWSC

Unsupported GPU dispatch cases fall back to CPU behavior where implemented. WENO, 2nd-order MUSCL, viscous terms, AMR, and the RK update are still CPU-side in the main loop.

## 9. Troubleshooting

### CUDA architecture mismatch

If you see errors like `no kernel image is available for execution on the device`, reconfigure with the correct architecture:

```bash
cmake -S . -B build-cuda -DCFD_CUDA_ARCHITECTURES=61
cmake --build build-cuda --target CFD_solver_gpu -j$(nproc)
```

### CUDA out of memory

Use a smaller grid, fewer iterations, and avoid WENO until the basic case runs.

### Missing VTK/JsonCpp/Boost

Install the development packages, then reconfigure from a clean build directory:

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

### MPI refuses to launch ranks

Use fewer ranks or:

```bash
mpirun -np 2 --oversubscribe ./build/CFD_solver_mpi json_Files/MPI_Smoke_Test.json
```

### Verify binaries

```bash
ls -lh build/CFD_solver build/CFD_solver_mpi build-cuda/CFD_solver_gpu
```

