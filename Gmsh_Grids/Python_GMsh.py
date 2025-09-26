import gmsh
import vtk
import numpy as np

gmsh.initialize()

# Define geometry parameters
length = 1.0
width = 0.5
center_x = length * 0.25
center_y = width * 0.5
radius = width * 0.2

# Create points
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(length, 0, 0)
p3 = gmsh.model.geo.addPoint(length, width, 0)
p4 = gmsh.model.geo.addPoint(0, width, 0)
center_point = gmsh.model.geo.addPoint(center_x, center_y, 0)

# Create lines
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Create the rectangle surface
curve_loop_rect = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
rectangle_surface = gmsh.model.geo.addPlaneSurface([curve_loop_rect])

# Create the circle arc
c1 = gmsh.model.geo.addCircleArc(p4, center_point, p1)

# Create the circle surface
curve_loop_circ = gmsh.model.geo.addCurveLoop([c1])
circle_surface = gmsh.model.geo.addPlaneSurface([curve_loop_circ])

# Subtract the circle from the rectangle
domain_tag = gmsh.model.geo.cut([(2, rectangle_surface)], [(2, circle_surface)])

# Define physical groups with names for boundary conditions
left_boundary_tag = 101
right_boundary_tag = 102
top_boundary_tag = 103
bottom_boundary_tag = 104
circle_boundary_tag = 105

gmsh.model.addPhysicalGroup(1, [l1], left_boundary_tag)
gmsh.model.setPhysicalGroupName(1, left_boundary_tag, "Inlet")

gmsh.model.addPhysicalGroup(1, [l2], right_boundary_tag)
gmsh.model.setPhysicalGroupName(1, right_boundary_tag, "Exit")

gmsh.model.addPhysicalGroup(1, [l3], top_boundary_tag)
gmsh.model.setPhysicalGroupName(1, top_boundary_tag, "Wall_Top")

gmsh.model.addPhysicalGroup(1, [l4], bottom_boundary_tag)
gmsh.model.setPhysicalGroupName(1, bottom_boundary_tag, "Wall_Bottom")

gmsh.model.addPhysicalGroup(1, [c1], circle_boundary_tag)  # Use the arc 'c1'
gmsh.model.setPhysicalGroupName(1, circle_boundary_tag, "Wall_Circle")

gmsh.model.addPhysicalGroup(2, [domain_tag[0][1]], 3)
gmsh.model.setPhysicalGroupName(2, 3, "Domain")

gmsh.model.geo.synchronize() # Important: synchronize the geometry before meshing

gmsh.model.mesh.generate(2)



gmsh.finalize()