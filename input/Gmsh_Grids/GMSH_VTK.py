import meshio
import numpy as np

# Read the gmsh mesh file (.msh)
mesh = meshio.read("Test_Domain_Structured.msh")
# Prepare the element numbering for each cell block
element_numbers = []  # this will hold a NumPy array for each block

# Prepare the element numbering: for each cell block, create a numbering array
cell_data = {}
new_data = []  # will collect numbering for each cell block
for cell_block in mesh.cells:
    num_cells = cell_block.data.shape[0]
    # Create a list of element numbers (0-based indexing; add 1 if you prefer 1-based)
    elem_numbers = np.arange(num_cells)
    element_numbers.append(elem_numbers)

# Add the element numbering to the cell data dictionary.
# Note: mesh.cell_data is a dict mapping data names to lists (one per cell block).
mesh.cell_data["ElementNumber"] = element_numbers

# Write out the mesh as a VTK file (legacy or XML format)
meshio.write("Test_Domain_Structured_With_Elements.vtk", mesh,binary=False)