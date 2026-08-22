import vtk

# Load the VTK file
file_path = "cylinder_mesh.vtk"

# Read the VTK file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(file_path)
reader.Update()

# Get the unstructured grid
mesh = reader.GetOutput()

def get_neighbors(cell_id):
    """
    Get the neighboring cells that share a face with the given cell_id.
    """
    neighbors = vtk.vtkIdList()
    
    # Get the points that form the cell
    cell_point_ids = vtk.vtkIdList()
    mesh.GetCellPoints(cell_id, cell_point_ids)

    # Find neighbors using these points
    mesh.GetCellNeighbors(cell_id, cell_point_ids, neighbors)

    neighbor_ids = [neighbors.GetId(i) for i in range(neighbors.GetNumberOfIds())]
    return neighbor_ids

# Example: Get neighbors for cell ID 0
cell_id = 0  # Change this to the desired cell ID
neighbors = get_neighbors(cell_id)

print(f"Neighbors of cell {cell_id}: {neighbors}")