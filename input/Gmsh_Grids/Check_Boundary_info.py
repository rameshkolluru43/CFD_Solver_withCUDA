import vtk

# Load the VTK file
file_path = "cylinder_mesh.vtk"  # Replace with your file path
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(file_path)
reader.Update()

# Get the unstructured grid
unstructured_grid = reader.GetOutput()

# Access the 'CellEntityIds' array
cell_data = unstructured_grid.GetCellData()
cell_entity_ids = cell_data.GetArray("CellEntityIds")

if cell_entity_ids:
    print("CellEntityIds:")
    for i in range(cell_entity_ids.GetNumberOfTuples()):
        print(f"Cell {i} -> Entity ID: {cell_entity_ids.GetValue(i)}")
else:
    print("No 'CellEntityIds' array found.")