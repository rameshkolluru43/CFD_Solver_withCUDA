# CFD Solver with CUDA — Project Report

**Date**: March 2025  
**Repository**: CFD_Solver_withCUDA

## Executive Summary

This project is a 2D compressible CFD solver (Euler and Navier–Stokes) with optional CUDA acceleration. It supports unstructured and mixed (triangle/quadrilateral) meshes, multiple flux and reconstruction schemes, and gradient-based adaptive mesh refinement (AMR) tagging. Documentation has been updated and new reference documents added for mesh support, AMR, configuration, and build/run.

## Capabilities

- **Equations**: 2D Euler and Navier–Stokes (compressible).
- **Meshes**: VTK legacy, structured TXT, CSV; **mixed triangles and quadrilaterals** with per-face connectivity.
- **Flux schemes**: LLF, Roe (1st/2nd), Van Leer, AUSM, RICCA, MOVERS; WENO on quad meshes.
- **Time**: Explicit (Euler, Runge–Kutta), optional implicit; CFL-based time step.
- **AMR**: Gradient-based refinement indicator (density + pressure); cell tagging every N iterations; configurable threshold and max fraction.
- **GPU**: CUDA kernels for flux, reconstruction, and grid operations; CPU and GPU executables built from one tree.

## Documentation Delivered

| Item | Location | Purpose |
|------|----------|---------|
| README update | [README.md](README.md) | Mixed mesh and AMR features; doc links |
| Mesh and grid | [docs/MESH_AND_GRID.md](docs/MESH_AND_GRID.md) | Formats, tri/quad, connectivity |
| AMR | [docs/ADAPTIVE_MESH_REFINEMENT.md](docs/ADAPTIVE_MESH_REFINEMENT.md) | Indicator, tagging, config |
| Configuration | [docs/CONFIGURATION.md](docs/CONFIGURATION.md) | JSON reference |
| Build and run | [docs/BUILD_AND_RUN.md](docs/BUILD_AND_RUN.md) | Build, run, troubleshooting |
| Release notes | [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md) | Changelog and recent changes |
| Doc index | [docs/README.md](docs/README.md) | Index of all docs |
| Project report | [PROJECT_REPORT.md](PROJECT_REPORT.md) | This report |

## Build and Run

```bash
mkdir -p build && cd build && cmake .. && make -j$(nproc)
./CFD_solver ../json_Files/Test_Config_AMR.json      # CPU with AMR
./CFD_solver_gpu ../json_Files/Solver_Config.json     # GPU
```

## Suggested Next Steps

- **AMR**: Implement actual mesh refinement (split quads/tris, new points/cells, solution transfer, re-connectivity).
- **VTK output**: Add `Gradient_Refinement_Indicator` and `RefineMarker` as cell data in solution VTK for visualization.
- **Validation**: Systematic tests on mixed meshes and with AMR tagging enabled.

---

*For technical details, see the docs in [docs/](docs/) and the main [README.md](README.md).*
