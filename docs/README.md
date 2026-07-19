# Documentation Index

## Getting Started

| Document | Description |
|---------|-------------|
| [../overview.md](../overview.md) | Project overview, GPU/host status, half-cylinder M=6 validation |
| [BUILD_AND_RUN.md](BUILD_AND_RUN.md) | Build requirements, CMake, running CPU/GPU executables |
| [STANDALONE_LOW_END_CUDA_RUN.md](STANDALONE_LOW_END_CUDA_RUN.md) | Standalone Linux setup and low-end CUDA GPU build/run guide |
| [CONFIGURATION.md](CONFIGURATION.md) | JSON configuration reference (simulation, solver, AMR) |
| [START_HERE_COLAB_TESTING.md](START_HERE_COLAB_TESTING.md) | Colab testing quick start |
| [CONNECT_TO_COLAB.md](CONNECT_TO_COLAB.md) | Connect VS Code/Cursor workflow to Colab |
| [COLAB_COMPILATION_GUIDE.md](COLAB_COMPILATION_GUIDE.md) | Colab CUDA build guide |
| API (Doxygen) | Run `doxygen Doxyfile_Cleaned` at repo root; open **`docs/doxygen/html/index.html`** |

**Note:** Older runs may have placed HTML under `docs/html/`; the configured output is now **`docs/doxygen/html/`** (see `OUTPUT_DIRECTORY` in `Doxyfile_Cleaned`).

## Mesh and Geometry

| Document | Description |
|---------|-------------|
| [MESH_AND_GRID.md](MESH_AND_GRID.md) | Mesh formats (VTK, TXT, CSV), mixed tri/quad, face-ordered connectivity |

## Features

| Document | Description |
|---------|-------------|
| [MAIN_LOOP_STATUS.md](MAIN_LOOP_STATUS.md) | Main solver loop: implemented vs pending (flux, viscous, AMR, GPU) |
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
| [PROJECT_REPORT.md](PROJECT_REPORT.md) | High-level project report |
| [PROJECT_STATUS_UPDATE.md](PROJECT_STATUS_UPDATE.md) | Project status summary |
| [QUICK_REFERENCE.md](QUICK_REFERENCE.md) | CUDA kernel quick reference |
| [TEST_SUMMARY.md](TEST_SUMMARY.md) | Test and validation summary |
| [TESTING_AND_VALIDATION_PLAN.md](TESTING_AND_VALIDATION_PLAN.md) | Validation plan |
| Metal_Kernels/ (repository) | Optional macOS Metal bridge and `.metal` shaders (not part of default CMake target) |
| [Grid_CUDA_Kernels_Documentation.md](Grid_CUDA_Kernels_Documentation.md) | CUDA grid kernels |
| [Implicit_Solver_Implementation.md](Implicit_Solver_Implementation.md) | Implicit solver |
| [Turbulence_Models_User_Guide.md](Turbulence_Models_User_Guide.md) | Turbulence models |

## Phase and Completion Reports

See `*_Completion_Summary.md`, `Phase*_*.md`, `*_GUIDE.md`, and `*_REPORT.md` in this directory for historical implementation phases, setup guides, and completion reports.

## Helper Scripts

Root-level shell helpers were moved into `../scripts/`:

| Script | Purpose |
|--------|---------|
| `../scripts/setup_colab.sh` | Install/setup dependencies in Colab |
| `../scripts/sync_to_colab.sh` | Push/pull files between local workspace and Colab |
| `../scripts/update_docs.sh` | Documentation maintenance helper |
| `../scripts/validate_cuda_syntax.sh` | CUDA syntax validation helper |
| `../scripts/view_docs.sh` | Open generated docs |
