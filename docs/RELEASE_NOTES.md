# Release Notes and Changelog

Summary of recent features and fixes.

---

## Half-cylinder M∞=6 viscous corner / BC fixes (July 2026)

- **Face normals vs BC faces**: Structured TXT grids keep vertex order `[o,a,b,c]`; `Construct_Cell` builds LBRT face normals (outward-corrected). Ghost face midpoints use the same LBRT edges. Fixes wall normals at side corners (~47 cells) that previously pointed at the y=0 exit plane.
- **Secondary neighbours**: Corner diagonal stencil for co-volume face gradients no longer collapses to cell 0 via invalid ghost walk.
- **Supersonic inlet**: Hard-prescribe freestream; do not mutate `inletCond` velocities during BC fill.
- **Host WENO walls**: Degenerate near-wall stencil falls back to 1O cell/ghost states.
- **Update positivity**: Reject non-positive cell updates (keep prior U) instead of flooring pressure (which produced M ≫ 6).
- **Viscous driver**: Remove double 1O flux evaluation before the explicit NS step.
- **Diagnostics**: `Dump_Viscous_Corner_Diagnostics` writes multi-BC corner face/normal info when viscous.
- **Plots**: `scripts/plot_halfcyl_inviscid_viscous.py` → `plots/half_cylinder_inviscid_viscous/`.
- **Docs**: Root [overview.md](../overview.md); [HALF_CYLINDER_VALIDATION.md](HALF_CYLINDER_VALIDATION.md); MESH/MAIN_LOOP/CONFIGURATION/BUILD refreshed; Doxygen main page + regenerate path fixed.

**Key files**: `Grid_Computations.cpp`, `Grid_Refine_Functions.cpp`, `Initialize.cpp`, `Inlet_Boundary_Conditions.cpp`, `WENO2D.cpp`, `Error_Estimate_Update.cpp`, `Solver.cpp`, `Viscous_Functions.h`.

---

## Documentation refresh (July 2026)

- Added **[HALF_CYLINDER_VALIDATION.md](HALF_CYLINDER_VALIDATION.md)** (configs, expected Pmax/Mmax, plots).
- Updated **MESH_AND_GRID.md** (LBRT / no atan2 reorder), **MAIN_LOOP_STATUS.md** (resident GPU, host WENO, corner fix), **CONFIGURATION.md**, **BUILD_AND_RUN.md**, **docs/README.md**, root **README.md** Doxygen commands.
- **Doxyfile_Cleaned** / **DoxygenMainPage.cpp**: project number and main page for v3.2 July 2026; regenerate with `./scripts/update_docs.sh` → `docs/doxygen/html/`.

---

## Documentation (April 2026)

- **README.md**: Dependencies aligned with CMake (Boost, OpenMP, VTK required for CPU build); Metal_Kernels called out; fixed corrupted section headings; Doxygen HTML path documented as `docs/doxygen/html/`.
- **Doxyfile_Cleaned**: `PROJECT_NUMBER` / `PROJECT_BRIEF` updated; `OUTPUT_DIRECTORY=docs/doxygen`; `INPUT` now includes `Metal_Kernels` and only `docs/DoxygenMainPage.cpp` (avoids accidentally scanning old `docs/html` trees); excludes `docs/html`, generated `docs/doxygen`, and common virtualenv paths; `GENERATE_LATEX=NO` by default for faster runs; `.mm`/`.metal` file patterns added.
- **docs/DoxygenMainPage.cpp**, **docs/README.md**: Main page and index refreshed for mixed mesh, AMR, CUDA vs Metal, and build prerequisites.

---

## AMR (status note)

Gradient AMR tagging is wired in the inviscid loop; quad split/merge via `AMR_Adaptive_Step` is available when enabled. See [MAIN_LOOP_STATUS.md](MAIN_LOOP_STATUS.md) and [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md). Older “tagging only” wording in some archives is outdated.

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
