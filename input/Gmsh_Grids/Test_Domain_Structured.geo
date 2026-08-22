// Define corner points
Point(1) = {0 , 0 , 0};
Point(2) = {10, 0, 0};
Point(3) = {10,10,0};
Point(4) = {0, 10, 0};

// Define lines connecting the points
Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 1}; 

// Set number of divisions along each line (for structured meshing)
Transfinite Line {1 , 2 , 3, 4} = 10;

// Define the surface enclosed by the lines
Line Loop(1) ={1, 2, 3, 4};
Plane Surface(1) = {1};

// Enforce structured meshing by specifying the four corner points
Transfinite Surface {1} = {1, 2, 3, 4};

// Ensure the mesh uses quadrilateral (recombined) elements
Recombine Surface {1};