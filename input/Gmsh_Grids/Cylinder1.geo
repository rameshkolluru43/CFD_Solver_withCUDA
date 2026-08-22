// ---- 2D Circular Cylinder Gmsh Mesh ----
// Converted from 3D to 2D (No Extrusion)

cl1 = 6;
cl2 = .03;
cl3 = 10;

radius = 1;
outer = 10;
numinner = 30; 

// Exterior bounding box of mesh
Point(1) = {-30, -30, 0, cl1};
Point(2) = { 50, -30, 0, cl3};
Point(3) = { 50,  30, 0, cl3};
Point(4) = {-30,  30, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cylinder and structured quadrants
Point(5) = {0,   0, 0, cl2};
Point(6) = {radius,  0, 0, cl2};
Point(7) = {0,  radius, 0, cl2};
Point(8) = {-radius, 0, 0, cl2};
Point(9) = {0, -radius, 0, cl2};
Point(10) = {outer,  0, 0, cl2};
Point(11) = {0,  outer, 0, cl2};
Point(12) = {-outer, 0, 0, cl2};
Point(13) = {0, -outer, 0, cl2};

Circle(5) = {8, 5, 9};
Circle(6) = {6, 5, 7};
Circle(7) = {10, 5, 11};
Circle(8) = {12, 5, 13};

Line(9)  = {6, 10};
Line(10) = {8, 12};
Line(11) = {9, 13};
Line(12) = {7, 11};

Circle(13) = {9, 5, 6};
Circle(14) = {8, 5, 7};
Circle(15) = {11, 5, 12};
Circle(16) = {13, 5, 10};

// Transfinite meshing for structured quadrants
Transfinite Line {5,6,7,8,13,14,15,16} = 30;  
Transfinite Line {9,10,11,12} = numinner Using Progression 1.05;    

// Line loops for surfaces
Line Loop(1) = {1, 2, 3, 4, 7, 16, 8, 15}; // Exterior region
Line Loop(2) = {10, 8, -11, -5}; // Bottom-right quad
Line Loop(3) = {7, -12, -6, 9}; // Top-right quad
Line Loop(4) = {-10, 14, 12, 15}; // **Fixed bottom-left quad**
Line Loop(5) = {16, -9, -13, 11}; // Top-left quad

// Define 2D Surfaces
Plane Surface(1) = {1}; // Outer unstructured region
Plane Surface(2) = {2}; // Bottom-right quad
Plane Surface(3) = {3}; // Top-right quad
Plane Surface(4) = {4}; // Bottom-left quad
Plane Surface(5) = {5}; // Top-left quad

// Structured meshing
Transfinite Surface {2,3,4,5};
Recombine Surface {2,3,4,5}; // Convert to quads
Recombine Surface {1};  // Convert outer region to quads (optional)

// Physical Boundary Conditions
Physical Surface("wall") = {1,2,3,4,5};  
Physical Surface("periodic_0_r") = {1,2,3,4,5};