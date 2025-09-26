// Define geometry parameters
length = 1.0;
width = 0.5;
center_x = length * 0.25;
center_y = width * 0.5;
radius = width * 0.2;

// Define rectangle corners
Point(1) = {0, 0, 0};
Point(2) = {length, 0, 0};
Point(3) = {length, width, 0};
Point(4) = {0, width, 0};

// Define circle points
Point(5) = {center_x + radius, center_y, 0};  // Right
Point(6) = {center_x, center_y + radius, 0};  // Top
Point(7) = {center_x - radius, center_y, 0};  // Left
Point(8) = {center_x, center_y - radius, 0};  // Bottom

// Define rectangle edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the circular hole using arcs
Circle(5) = {5, center_x, 6};  
Circle(6) = {6, center_x, 7};  
Circle(7) = {7, center_x, 8};  
Circle(8) = {8, center_x, 5};  

// Define structured quadrilateral mesh around the circle
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

// Create subregions with exact 4-corner quads
Curve Loop(13) = {5, -9, -6, 10};  // Top-Right
Curve Loop(14) = {6, -10, -7, 11};  // Top-Left
Curve Loop(15) = {7, -11, -8, 12};  // Bottom-Left
Curve Loop(16) = {8, -12, -5, 9};  // Bottom-Right

// Define structured quadrilateral mesh regions
Plane Surface(17) = {13};  
Plane Surface(18) = {14};  
Plane Surface(19) = {15};  
Plane Surface(20) = {16};  

// Apply transfinite meshing only on surfaces with 4 corners
Transfinite Surface {17, 18, 19, 20};  
Recombine Surface {17, 18, 19, 20};  // Ensures quadrilateral elements

// Define physical groups for boundary conditions
Physical Curve("Inlet") = {1}; 
Physical Curve("Exit") = {2}; 
Physical Curve("Wall_Top") = {3}; 
Physical Curve("Wall_Bottom") = {4}; 
Physical Curve("Wall_Circle") = {5, 6, 7, 8}; 
Physical Surface("Domain") = {17, 18, 19, 20};  

// Mesh size control
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.05;

// Generate mesh
Mesh 2;