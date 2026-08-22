# Input (cases, grids, generators)

| Path | Role |
|------|------|
| `json_Files/` | Solver run JSON and test-case JSON (`Test_Case_Json`, CFL, scheme, VTK paths). |
| `Grid_Files/` | Structured TXT/VTK meshes referenced by 2D test cases. |
| `Gmsh_Grids/` | GMSH meshes. |
| `2D_Half_Cylinder_Generate/` | Half-cylinder mesh generators. |
| `VTK_Grid_Tests/` | VTK I/O experiments. |

3D Cartesian cases (SWBLI, cavity, ramp) generate the mesh in the solver. 2D cases use `Grid_Files/` (run the 2D binary from `build/`). See [LAYOUT.md](../LAYOUT.md).
