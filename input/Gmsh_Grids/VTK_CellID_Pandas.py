import vtk
import pandas as pd

# Load the VTK file
file_path = "cylinder_mesh.vtk"

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(file_path)
reader.Update()

# Get the unstructured grid
mesh = reader.GetOutput()

# Dictionary to store neighbors
cell_neighbors = {}

# Iterate through all cells and extract connectivity
for cell_id in range(mesh.GetNumberOfCells()):
    cell = mesh.GetCell(cell_id)
    num_points = cell.GetNumberOfPoints()
    point_ids = [cell.GetPointId(i) for i in range(num_points)]

    # Store point IDs for this cell
    cell_neighbors[cell_id] = {
        "Point IDs": point_ids,
        "Neighbors": set()  # Using a set to prevent duplicate neighbors
    }

# Find neighbors by checking shared edges
for cell_id in range(mesh.GetNumberOfCells()):
    cell_points = set(cell_neighbors[cell_id]["Point IDs"])  # Set of points for the current cell

    for other_cell_id in range(mesh.GetNumberOfCells()):
        if cell_id == other_cell_id:
            continue  # Skip self

        other_cell_points = set(cell_neighbors[other_cell_id]["Point IDs"])

        # Find shared points (common edge)
        shared_points = cell_points.intersection(other_cell_points)

        # For quads, a neighbor must share exactly 2 points (1 edge)
        if len(shared_points) == 2:
            cell_neighbors[cell_id]["Neighbors"].add(other_cell_id)

# Convert to Pandas DataFrame
df_data = []
for cell_id, info in cell_neighbors.items():
    df_data.append({
        "Cell ID": cell_id,
        "Point IDs": info["Point IDs"],
        "Neighbors": list(info["Neighbors"])  # Convert set to list
    })

df = pd.DataFrame(df_data)

# Print all rows properly
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)

# Print the DataFrame
print(df.to_string(index=False))