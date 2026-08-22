// ---- 2D Circular Cylinder Gmsh Grid ----
// Modified to be strictly 2D (no extrusion)
// --------------------------------------------

cl1 = 6;
cl2 = .01;
cl3 = 10;

radius = 1;
outer = 10;
numinner = 30;

// Exterior (bounding box) of mesh
Point(1) = {-20, -20, 0, cl1};
Point(2) = { 60, -20, 0, cl3};
Point(3) = { 60,  20, 0, cl3};
Point(4) = {-20,  20, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Circle & surrounding structured-quad region
Point(5) = {0,   0, 0, cl2};
Point(6) = {0,  radius, 0, cl2};
Point(7) = {0, -radius, 0, cl2};
Point(8) = {0,  outer, 0, cl2};
Point(9) = {0, -outer, 0, cl2};
Point(10) = {radius,  0, 0, cl2};
Point(11) = {-radius, 0, 0, cl2};
Point(12) = {outer,  0, 0, cl2};
Point(13) = {-outer, 0, 0, cl2};

Circle(5) = {7, 5, 10};
Circle(6) = {6, 5, 11};
Circle(7) = {8, 5, 13};
Circle(8) = {9, 5, 12};

Line(9)  = {6, 8};
Line(10) = {7, 9};
Line(11)  = {10, 12};
Line(12) = {11, 13};

Circle(13) = {10, 5, 6};
Circle(14) = {11, 5, 7};
Circle(15) = {13, 5, 9};
Circle(16) = {12, 5, 8};

// Transfinite constraints to enforce structured quad meshing
Transfinite Line {5,6,7,8,13,14,15,16} = 20;
Transfinite Line {9,10,11,12} = numinner Using Progression 1.1;
Transfinite Line {1,2,3,4} = 40 Using Progression 1.0;
Transfinite Line {2,3} = 42 Using Progression 1.0; // Ensure exit boundary is structured

// Define structured quad regions
Line Loop(1) = {1, 2, 3, 4, 7, 16, 8, 15}; // Exterior
Line Loop(2) = {10, 8, -11, -5}; // Right-hand structured region
Line Loop(3) = {7, -12, -6, 9}; // Left-hand structured region
Line Loop(4) = {-10, -14, 12, 15}; // Right-hand structured region
Line Loop(5) = {16, -9, -13, 11}; // Left-hand structured region

Plane Surface(1) = {1}; // Outer unstructured region
Plane Surface(2) = {2}; // RH inner structured region
Plane Surface(3) = {3}; // LH inner structured region
Plane Surface(4) = {4}; // RH inner structured region
Plane Surface(5) = {5}; // LH inner structured region

// Enforce structured quad mesh
Transfinite Surface{2,3,4,5};
Recombine Surface {2,3,4,5}; // Ensure quads only
Recombine Surface {1,2,3,4,5}; // Outer unstructured quad region

// Force quadrilateral-only meshing

// Apply boundary conditions
Physical Line("Bottom") = {1};
Physical Line("Right")  = {2};
Physical Line("Top")    = {3};
Physical Line("Left")   = {4};
Physical Line("Circle") = {5,6};

// Define fluid domain
Physical Surface("FLUID") = {1,2,3,4,5};
