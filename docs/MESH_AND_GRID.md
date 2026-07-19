# Mesh and Grid Support

This document describes the mesh formats, cell types, and connectivity conventions used by the CFD solver.

## Supported Mesh Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| **VTK Legacy** | `.vtk` | Unstructured 2D mesh; points + CELLS + CELL_TYPES |
| **Structured TXT** | `.txt` | Cartesian-style grid (e.g. Half_Cylinder_61_21.txt) |
| **CSV** | `.csv` | Node lists (x, y) for structured-like quad grids |

The main entry for loading a mesh is via the test case JSON: the `Test_Case_Json` points to a case file that specifies `Grid_File` or mesh paths. The solver also supports `Load_Mesh(configOrMeshPath)` for direct mesh or JSON config paths.

## Cell Types (2D)

- **Triangle** (VTK type 5): 3 nodes per cell, 3 faces.
- **Quadrilateral** (VTK type 8): 4 nodes per cell, 4 faces.

The solver supports **mixed meshes**: triangles and quadrilaterals can coexist. All geometric and flux logic uses the per-cell number of faces (`numFaces`) and face-ordered neighbour lists.

## Connectivity: Face-Ordered Neighbours

For each cell, the following are **aligned by face index**:

- **`Neighbours[f]`**: Cell index on the other side of face `f`, or a **ghost cell index** if the face is on a boundary.
- **Face `f`**: Edge from vertex `f` to vertex `(f+1) % numNodes`, in the order defined by `nodeIndices` / `Cell_Vertices` (anticlockwise in 2D).
- **`Face_Areas[f]`**: Length of face `f` (in 2D).
- **`Face_Normals[2*f]`, `Face_Normals[2*f+1]`**: Outward normal (nx, ny) for face `f`.

So for any face index `f`, `Neighbours[f]` is the neighbour across that face. Boundary faces have `Neighbours[f] >= No_Physical_Cells` (ghost).

## Boundary Lists

Boundary conditions use triplets **(cell index, face index, ghost cell index)**:

- **Inlet_Cells_List**, **Exit_Cells_List**, **Wall_Cells_List**: each entry is `[cell, face_index, ghost_id]`.
- **Cells_Face_Boundary_Type[cell][face]** is `true` when that face is a boundary face (wall/inlet/exit).

Face indices are 0-based and match the cell’s cyclic vertex order.

## Geometry

- **Cell centroid**: From vertices; triangles use centroid, quads/polygons use the same centroid formula as in `Compute_Centroid`.
- **Cell area**: For 2D polygons, area is computed with the **shoelace formula** (signed area in the xy-plane).
- **Ghost cell centres**: Reflected from the cell centre across the **midpoint of the boundary face** (generic for any polygon).

## File Conventions

### VTK (Legacy)

- **POINTS**: `N x y z` then N lines of `x y z`.
- **CELLS**: `N_cells N_total_entries` then for each cell: `num_nodes n0 n1 ...`.
- **CELL_TYPES**: One type per cell (e.g. 5 = triangle, 8 = quadrilateral).

Node indices in VTK are 0-based. The reader builds `Cells` with `nodeIndices` and then fills `Cell_Vertices` from the global point list; vertices are ordered anticlockwise for 2D.

### Structured TXT (LBRT convention)

Used by half-cylinder and similar cases (e.g. `Grid_Files/Half_Cylinder/Half_Cylinder_481_161.txt`). The grid is a structured set of quadrilaterals; ghost count and boundary lists come from the block dimensions. After loading, `Cells` is resized to `Total_No_Cells` (physical + ghost) before constructing ghost cells.

**Critical indexing (July 2026 fix):**

| Index | Neighbour | Face edge (verts `[o,a,b,c]`) |
|-------|-----------|-------------------------------|
| 0 | Left | `c → o` |
| 1 | Bottom | `o → a` |
| 2 | Right | `a → b` |
| 3 | Top | `b → c` |

Rules enforced in `Read_Grid` / `Construct_Cell`:

1. **Keep file vertex order** `[o,a,b,c]`. Do **not** call `Sort_Points_AntiClockWise` / atan2 reorder on structured TXT — that misaligns `Face_Normals` with Neighbours/BC faces near side corners (~47 wall cells on the 481×161 half cylinder).
2. Build face areas/normals in classic **LBRT** order; flip if the normal is not outward from the cell centre.
3. Ghost face midpoints use the same LBRT vertex pairs (`Construct_Cell(cell, face, ghost)`).
4. Corner secondary neighbours for co-volume face gradients use corrected SW/SE/NE/NW logic (`Identify_Neighbours_For_Second_Gradients`); optional dump via `Dump_Viscous_Corner_Diagnostics`.

Without (1)–(3), viscous continues from a good inviscid P3 field can blow up (Pmax ~137, Mmax ~64). After the fix, smoke stays near the inviscid reference (Pmax ≈ 48.9, Mmax ≈ 6.05). See [HALF_CYLINDER_VALIDATION.md](HALF_CYLINDER_VALIDATION.md).

## Limits and Notes

- **WENO** reconstruction is currently implemented for **quad-only** meshes (4 faces). On mixed or triangular meshes, the solver falls back to first- or second-order flux (MUSCL/LLF, etc.) when WENO is selected. Near-wall host WENO faces with a degenerate stencil fall back to first-order cell/ghost states.
- **Viscous** time-step and **co-volume** / **implicit assembly** paths assume quadrilaterals; for non-quad cells they are skipped or use an inviscid-style time step.
- **AMR**: gradient tagging works on all cell types; the inviscid loop can also run quad 1→4 split / 4→1 merge when enabled (`AMR_Adaptive_Step`). See [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md) and [MAIN_LOOP_STATUS.md](MAIN_LOOP_STATUS.md).
