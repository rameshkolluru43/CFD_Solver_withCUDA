
// Enable OpenCASCADE kernel for Boolean operations
SetFactory("OpenCASCADE");

// Define geometry parameters
length = 10.0;
width = 5;
center_x = length * 0.25;
center_y = width * 0.5;
radius = width * 0.02;

numCirclePts = 60;

// Create rectangle (outer domain)
Rectangle(1) = {0, 0, 0, length, width};

// Create circle (inner void)
Disk(2) = {center_x, center_y, 0, radius, radius};

// Perform Boolean Difference to subtract the circle
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// Extract boundary lines of the rectangle
Curve Loop(4) = {1, 2, 3, 4}; // Outer rectangle boundary

// Apply transfinite meshing with 60x30 grid resolution
Transfinite Curve{1, 3} = 60 Using Progression 1; // X-direction (horizontal)
Transfinite Curve{2, 4} = 30 Using Progression 1; // Y-direction (vertical)
Transfinite Curve{5} = numCirclePts Using Progression 1; // Circular void

// Ensure quadrilateral elements
Transfinite Surface{3}; // Explicitly set transfinite surface
Recombine Surface{3};
Mesh.RecombineAll = 1;


// Define physical groups (boundary conditions)
Physical Surface("FluidDomain") = {3}; // Final meshed domain
Physical Curve("Left") = {1};
Physical Curve("Right") = {3};
Physical Curve("Bottom") = {2};
Physical Curve("Top") = {4};


// Mesh generation
Mesh 2;
