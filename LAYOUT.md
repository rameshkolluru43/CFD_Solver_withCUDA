# Repository layout

All paths are physical folders (no compatibility symlinks). Configure and run from the **repository root** unless noted.

```
.
├── CMakeLists.txt          Top-level build (2D CPU/GPU/MPI + 3D CPU/GPU)
├── LAYOUT.md               This file
├── README.md               Project overview and usage
├── solvers/                Source code
│   ├── src/                2D compressible Euler/NS (`Main.cpp`)
│   ├── include/            2D headers
│   ├── src3d/              3D compressible Euler/NS (use this, not `src_3D/`)
│   ├── CUDA_KERNELS/       CUDA (2D kernels + `Euler3D_Cuda.cu`)
│   ├── Metal_Kernels/      Optional Apple Metal
│   └── …                   Legacy: `src_3D/`, incompressible, Laplace, MPI_Src
├── input/                  Cases and meshes
│   ├── json_Files/         Run JSON (`Run_3D_*.json`, half-cylinder, ramp)
│   ├── Grid_Files/         2D structured TXT/VTK
│   ├── Gmsh_Grids/
│   └── 2D_Half_Cylinder_Generate/
├── wrappers/
│   ├── python/             Python / pybind helpers
│   └── scripts/            Docs, Colab, scheme sweeps, plots
├── tests/                  `Test_Cases/`, unit tests, examples
├── solutions/              VTK + plots (large VTK gitignored)
│   └── 2D_Euler_Solutions/ Cavity, SWBLI, ramp, half-cylinder
├── docs/                   Guides (`overview.md`, BUILD_AND_RUN, …)
├── cmake/                  Extra CMake (Colab)
├── tools/                  Profiling leftovers, local `.deb` extract
└── build/                  CMake output (gitignored)
```

## Binaries (`build/`)

| Target | Needs | Role |
|--------|--------|------|
| `CFD_solver` | Boost, JsonCpp, OpenMP | 2D CPU (VTK optional; TXT meshes work without it) |
| `CFD_solver_mpi` | MPI + above | 2D CPU, replicated mesh, rank-owned cells |
| `CFD_solver_gpu` | CUDA (`nvcc`) | 2D GPU + host fallback |
| `CFD_solver_3d` | OpenMP | 3D host Euler/NS/LES |
| `CFD_solver_3d_gpu` | CUDA | 3D CUDA (resident FO; WENO/viscous hybrid) |

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

If Boost/JsonCpp live in a user prefix:

```bash
cmake -S . -B build -DBOOST_ROOT="$HOME/deps/usr"
```

MPI headers can live under `$HOME/deps/.../openmpi/include` when `libopenmpi-dev` is not installed system-wide. OpenMPI runtime (`libmpi.so.40`, `mpirun`) is still required.

## Working directory

| Solver | Typical cwd | Why |
|--------|-------------|-----|
| 3D `CFD_solver_3d[_gpu]` | **repo root** | `VTK_Out` is `solutions/2D_Euler_Solutions/...` |
| 2D `CFD_solver[_gpu|_mpi]` | **`build/`** | Grid paths `../input/Grid_Files/...`; output under `../solutions/...` |

## Example runs

```bash
# 3D from repo root
./build/CFD_solver_3d input/json_Files/Run_3D_Sod_LLF_smoke.json
./build/CFD_solver_3d_gpu input/json_Files/Run_3D_Cavity_M6_LES_WALE.json

# 2D from build/
cd build
OMP_NUM_THREADS=4 ./CFD_solver ../input/json_Files/MPI_Smoke_Test.json
mpirun --oversubscribe -np 2 ./CFD_solver_mpi ../input/json_Files/MPI_Smoke_Test.json
./CFD_solver_gpu ../input/json_Files/Run_HalfCylinder_M6_WENO_wallfix_smoke50.json
```

2D `Test_Case_Json` values like `"../json_Files/Half_Cylinder.json"` are resolved **relative to the run JSON file** (under `input/json_Files/`), not the cwd.

## Solutions (kept production)

See [solutions/README.md](solutions/README.md). Smoke VTK and duplicate sweep copies were removed.

| Folder | Contents |
|--------|----------|
| `solutions/2D_Euler_Solutions/Cavity_3D/` | Mach-6 open-cavity LES (90° block mesh) |
| `solutions/2D_Euler_Solutions/SWBLI_3D/` | Laminar plate + WALE LES |
| `solutions/2D_Euler_Solutions/Ramp_15_Degree_3D/` | Finest 3D ramp |
| `solutions/2D_Euler_Solutions/Half_Cylinder/` | M=6 GridSize_7 |
| `solutions/plots/` | Half-cylinder figures |

## Docs and scripts

- Guides: [docs/README.md](docs/README.md), [docs/overview.md](docs/overview.md)
- Doxygen: `wrappers/scripts/update_docs.sh` (from repo root)
- 3D code notes: [solvers/src3d/README.md](solvers/src3d/README.md)
