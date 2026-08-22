import meshio

# Input SU2 file
input_file = "cylinder_mesh.su2"  # Replace with your SU2 file
# Output VTK file
output_file = "cylinder_mesh.vtk"

# Read SU2 mesh and write to VTK format
mesh = meshio.read(input_file)
meshio.write(output_file, mesh)

print(f"Converted {input_file} to {output_file}")