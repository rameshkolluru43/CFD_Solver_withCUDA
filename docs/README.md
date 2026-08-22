# Documentation Index

## Getting Started

| Document | Description |
|---------|-------------|
| [overview.md](overview.md) | Project overview, GPU/host status, half-cylinder M=6 validation |
| [../LAYOUT.md](../LAYOUT.md) | Folder map (solvers, input, wrappers, solutions) |
| [HALF_CYLINDER_VALIDATION.md](HALF_CYLINDER_VALIDATION.md) | Canonical M=6 configs, expected Pmax/Mmax, plots |
| [BUILD_AND_RUN.md](BUILD_AND_RUN.md) | Build requirements, CMake, running CPU/GPU |
| [CONFIGURATION.md](CONFIGURATION.md) | JSON configuration reference (solver, viscous, AMR, Dissipation_Type) |
| [STANDALONE_LOW_END_CUDA_RUN.md](STANDALONE_LOW_END_CUDA_RUN.md) | Low-end CUDA GPU build/run |
| API (Doxygen) | `./wrappers/scripts/update_docs.sh` or `doxygen docs/Doxyfile_Cleaned` → **`docs/doxygen/html/index.html`** |

**Note:** Older trees under `docs/html/` are obsolete. Configured output is **`docs/doxygen/html/`** (`OUTPUT_DIRECTORY` in `docs/Doxyfile_Cleaned`).

### Colab / remote (optional)

| Document | Description |
|---------|-------------|
| [START_HERE_COLAB_TESTING.md](START_HERE_COLAB_TESTING.md) | Colab testing quick start |
| [CONNECT_TO_COLAB.md](CONNECT_TO_COLAB.md) | VS Code/Cursor ↔ Colab |
| [COLAB_COMPILATION_GUIDE.md](COLAB_COMPILATION_GUIDE.md) | Colab CUDA build |

## Mesh, Loop, and Features

| Document | Description |
|---------|-------------|
| [MESH_AND_GRID.md](MESH_AND_GRID.md) | Mesh formats, **LBRT** structured TXT, face normals, mixed tri/quad |
| [MAIN_LOOP_STATUS.md](MAIN_LOOP_STATUS.md) | Main solver loop: host vs CUDA, WENO, viscous, AMR |
| [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md) | Gradient AMR indicator / tagging / split |
| [RELEASE_NOTES.md](RELEASE_NOTES.md) | Changelog (incl. July 2026 corner fix) |

## Flux and Numerical Methods

| Document | Description |
|---------|-------------|
| [Van_Leer_Flux_Implementation.md](Van_Leer_Flux_Implementation.md) | Van Leer flux vector splitting |
| [ROE_2O_Implementation.md](ROE_2O_Implementation.md) | Second-order Roe |
| [Enhanced_ROE_First_Order.md](Enhanced_ROE_First_Order.md) | First-order Roe with entropy fix |
| [AUSM_Flux_Implementation.md](AUSM_Flux_Implementation.md) | AUSM flux |
| [WENO2D_Validation_Report.md](WENO2D_Validation_Report.md) | WENO2D validation + OA script note |

## Implementation / historical

Phase reports, Jan 2026 CUDA roadmaps, and `*_Completion_Summary.md` files in this directory are **historical**. Prefer `overview.md`, `MAIN_LOOP_STATUS.md`, and `RELEASE_NOTES.md` for current status.

| Document | Description |
|---------|-------------|
| [PROJECT_REPORT.md](PROJECT_REPORT.md) | High-level project report (archive) |
| [QUICK_REFERENCE.md](QUICK_REFERENCE.md) | Older CUDA kernel quick reference (archive) |
| [CUDA_IMPLEMENTATION_STATUS_AND_ROADMAP.md](CUDA_IMPLEMENTATION_STATUS_AND_ROADMAP.md) | Jan 2026 roadmap (archive) |
| [Turbulence_Models_User_Guide.md](Turbulence_Models_User_Guide.md) | RANS turbulence models |
| [Implicit_Solver_Implementation.md](Implicit_Solver_Implementation.md) | Implicit solver notes |
| Metal_Kernels/ (repo) | Optional macOS Metal bridge (not default CMake target) |

## Helper Scripts

| Script | Purpose |
|--------|---------|
| `../wrappers/scripts/update_docs.sh` | Regenerate Doxygen → `docs/doxygen/html` |
| `../wrappers/scripts/view_docs.sh` | Open generated HTML docs |
| `../wrappers/scripts/plot_halfcyl_inviscid_viscous.py` | Half-cylinder inviscid/viscous plots |
| `../wrappers/scripts/test_weno_order_accuracy.py` | WENO order-of-accuracy helper |
| `../wrappers/scripts/validate_cuda_syntax.sh` | CUDA syntax check |
| `../wrappers/scripts/setup_colab.sh` | Colab dependency setup |
| `../wrappers/scripts/sync_to_colab.sh` | Sync workspace ↔ Colab |
