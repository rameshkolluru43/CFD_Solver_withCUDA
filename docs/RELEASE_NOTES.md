# Release Notes and Changelog

Summary of recent features and fixes.

---

## Mixed Mesh Support (Triangles and Quadrilaterals)

- **Per-face connectivity**: `Neighbours[f]` is the neighbour across face `f`; face order matches the cell’s cyclic vertex order (VTK path).
- **Generic geometry**: Cell centroid, area (shoelace for 2D polygons), face normals and lengths, and ghost cell centres are computed for any polygon (tri/quad).
- **Boundary lists**: Inlet/exit/wall lists store `(cell, face_index, ghost_id)`; `Cells_Face_Boundary_Type[cell][face]` is used for boundary flags.
- **Flux and time step**: Flux loops and inviscid time step use `numFaces` and face-indexed neighbours; viscous time step and co-volume/implicit assembly are quad-only or skip non-quad cells.
- **WENO**: Enabled only on all-quad meshes; otherwise the solver falls back to first- or second-order flux.

**Files**: `Read_Gmsh_File.cpp`, `Grid_Computations.cpp`, `Net_Flux.cpp`, `Time_Step.cpp`, `Initialize.cpp`, `Limiters.cpp`, `Van_Leer.cpp`, `Ausm_Flux.cpp`, `output_files.cpp`, `Co_Volume_Grid_Computations.cpp`, `Assemble_Matrix.cpp`, and related headers.

---

## Gradient-Based Adaptive Mesh Refinement (AMR)

- **Indicator**: Green–Gauss gradient at cell centre (any polygon); refinement indicator combines \(|\nabla\rho|\) and \(|\nabla P|/P\) scaled by \(\sqrt{A}\).
- **Tagging**: Cells above `AMR_Gradient_Threshold` are marked `Is_Splittable`; optional cap via `AMR_Max_Fraction` (refine top fraction by indicator).
- **Solver hook**: Every `AMR_Period` iterations, indicator is computed and cells are tagged; count is printed. No mesh topology change yet (tagging only).
- **Config**: JSON keys `Enable_AMR`, `AMR_Period`, `AMR_Gradient_Threshold`, `AMR_Max_Fraction` under `Solver`.

**Files**: `Grid_Refine_Functions.cpp`, `Solver.cpp`, `Configuration_Read.cpp`, `Globals.h`, `Initialize.cpp`, `Grid.h`; `json_Files/Test_Config_AMR.json`.

---

## Bug Fixes and Consistency

- **WENO (CPU/CUDA)**: Wave speeds \(S_0\) and \(S_3\) use \(fabs(Vdotn - C)\); right-state polynomial and smoothness indicators corrected in CUDA WENO kernels.
- **LLF wave speeds**: Same \(fabs(Vdotn - C)\) fix in CPU LLF and WENO flux.
- **Ghost cells**: Structured (TXT) path resizes `Cells` to `Total_No_Cells` before constructing ghost cells; `Check_Cells` accepts physical+ghost count.
- **Viscous time step**: Neighbour indices for quads use `Neighbours[0..3]` (was 1..4); non-quad cells use inviscid time step.

---

## Documentation Added/Updated

- **README.md**: Mixed mesh and AMR features; doc table linking to new docs.
- **docs/MESH_AND_GRID.md**: Mesh formats, tri/quad, face-ordered connectivity, boundaries.
- **docs/ADAPTIVE_MESH_REFINEMENT.md**: AMR indicator, tagging, configuration.
- **docs/CONFIGURATION.md**: JSON reference including AMR.
- **docs/BUILD_AND_RUN.md**: Build, run, and troubleshooting.
- **docs/RELEASE_NOTES.md**: This changelog.

---

*For older implementation details (flux schemes, CUDA, phases), see the other documents in `docs/` and root-level `*_SUMMARY.md` / `*_GUIDE.md` files.*
