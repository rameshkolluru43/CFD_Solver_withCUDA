import vtk

# Load the VTK file
vtk_file_path = "cylinder_mesh.vtk"  # Replace with the path to your .vtk file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(vtk_file_path)
reader.Update()

# Get the unstructured grid
unstructured_grid = reader.GetOutput()

# Extract available point and cell data arrays
point_data = unstructured_grid.GetPointData()
cell_data = unstructured_grid.GetCellData()

# Print available data arrays
if point_data:
    print("Point Data Arrays:")
    for i in range(point_data.GetNumberOfArrays()):
        print(f" - {point_data.GetArrayName(i)}")

if cell_data:
    print("Cell Data Arrays:")
    for i in range(cell_data.GetNumberOfArrays()):
        print(f" - {cell_data.GetArrayName(i)}")