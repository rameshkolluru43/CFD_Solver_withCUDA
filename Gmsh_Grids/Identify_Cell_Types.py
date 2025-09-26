import vtk

# Load the VTK file
file_path = "cylinder_mesh.vtk"  # Update with your file path

# Read the VTK file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(file_path)
reader.Update()
mesh = reader.GetOutput()

# Extract cell types
cell_types = {}
num_cells = mesh.GetNumberOfCells()

for i in range(num_cells):
    cell = mesh.GetCell(i)
    cell_type = cell.GetCellType()
    
    # Count cell types
    if cell_type in cell_types:
        cell_types[cell_type] += 1
    else:
        cell_types[cell_type] = 1

# Mapping of VTK cell types to human-readable names
vtk_cell_types = {
    3: "Line Element",
    5: "Triangle",
    9: "Quad",
    10: "Tetra",
    12: "Hexahedron",
    13: "Wedge",
    14: "Pyramid",
}

# Convert cell types to readable format
cell_types_readable = {vtk_cell_types.get(k, f"Unknown({k})"): v for k, v in cell_types.items()}

# Display results
print("\nIdentified Cell Types:")
for cell_name, count in cell_types_readable.items():
    print(f"{cell_name}: {count} cells")