# Documentation Index

## Getting Started

| Document | Description |
|---------|-------------|
| [BUILD_AND_RUN.md](BUILD_AND_RUN.md) | Build requirements, CMake, running CPU/GPU executables |
| [CONFIGURATION.md](CONFIGURATION.md) | JSON configuration reference (simulation, solver, AMR) |
| API (Doxygen) | Run `doxygen Doxyfile_Cleaned` at repo root; open **`docs/doxygen/html/index.html`** |

**Note:** Older runs may have placed HTML under `docs/html/`; the configured output is now **`docs/doxygen/html/`** (see `OUTPUT_DIRECTORY` in `Doxyfile_Cleaned`).

## Mesh and Geometry

| Document | Description |
|---------|-------------|
| [MESH_AND_GRID.md](MESH_AND_GRID.md) | Mesh formats (VTK, TXT, CSV), mixed tri/quad, face-ordered connectivity |

## Features

| Document | Description |
|---------|-------------|
| [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md) | Gradient-based AMR: indicator, tagging, JSON options |

## Flux and Numerical Methods

| Document | Description |
|---------|-------------|
| [Van_Leer_Flux_Implementation.md](Van_Leer_Flux_Implementation.md) | Van Leer flux vector splitting |
| [ROE_2O_Implementation.md](ROE_2O_Implementation.md) | Second-order Roe scheme |
| [Enhanced_ROE_First_Order.md](Enhanced_ROE_First_Order.md) | First-order Roe with entropy fix |
| [AUSM_Flux_Implementation.md](AUSM_Flux_Implementation.md) | AUSM flux scheme |

## Implementation Summaries

| Document | Description |
|---------|-------------|
| [RELEASE_NOTES.md](RELEASE_NOTES.md) | Changelog and recent updates (mixed mesh, AMR, fixes) |
| Metal_Kernels/ (repository) | Optional macOS Metal bridge and `.metal` shaders (not part of default CMake target) |
| [Grid_CUDA_Kernels_Documentation.md](Grid_CUDA_Kernels_Documentation.md) | CUDA grid kernels |
| [Implicit_Solver_Implementation.md](Implicit_Solver_Implementation.md) | Implicit solver |
| [Turbulence_Models_User_Guide.md](Turbulence_Models_User_Guide.md) | Turbulence models |

## Phase and Completion Reports

See `*_Completion_Summary.md` and `Phase*_*.md` in this directory for historical implementation phases and completion reports.
