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

### Structured TXT

Used by some test cases (e.g. Half_Cylinder). The grid is built as a structured set of quadrilaterals; ghost count and boundary lists are derived from the block dimensions. After loading, `Cells` is resized to `Total_No_Cells` (physical + ghost) before constructing ghost cells.

## Limits and Notes

- **WENO** reconstruction is currently implemented for **quad-only** meshes (4 faces). On mixed or triangular meshes, the solver falls back to first- or second-order flux (MUSCL/LLF, etc.) when WENO is selected.
- **Viscous** time-step and **co-volume** / **implicit assembly** paths assume quadrilaterals; for non-quad cells they are skipped or use an inviscid-style time step.
- **Refinement (AMR)** tagging is gradient-based and works on all cell types; actual mesh refinement (splitting) is currently implemented only as tagging (no mesh change yet).
