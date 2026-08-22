SetFactory("OpenCASCADE");

lc1 = 1e-2;
numSegments = 20; // Discretize the cylinder into segments

// Points defining the cylinder and domain
Point(1) = {   0,   0 , 0, lc1}; // Cylinder centre
Point(2) = {   0, -0.05, 0, lc1}; // Cylinder bottom point
Point(3) = {   0,  0.05, 0, lc1}; // Cylinder top point
Point(4) = {-1.1, -1.1, 0, 10*lc1}; // Bottom-left
Point(5) = {-1.1,  1.1, 0, 5*lc1};  // Top-left
Point(6) = { 5.1,  1.1, 0, 5*lc1};  // Top-right
Point(7) = { 5.1, -1.1, 0, 5*lc1};  // Bottom-right

// Discretize the cylinder with line segments instead of a continuous circle
Spline(1) = {2,3};
Spline(2) = {3,2};

// Defining external boundaries
Line(3)   = {4,5}; // Inlet
Line(4)   = {5,6}; // Top boundary
Line(5)   = {6,7}; // Outlet
Line(6)   = {7,4}; // Bottom boundary

// Creating structured quadrilateral mesh by splitting into smaller transfinite patches
Line Loop(1) = {3,4,5,6};
Line Loop(2) = {1,2};
Plane Surface(1) = {1,2};

// Assigning boundary conditions
Physical Line("inlet", 1) = {3};
Physical Line("outlet", 2) = {5};
Physical Line("bottom", 3) = {6};
Physical Line("top", 4) = {4};
Physical Line("cylinder", 5) = {1,2};
Physical Surface("flow") = {1};

// Enforce structured quad mesh using transfinite
Transfinite Line {3,4,5,6,1,2} = numSegments Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1}; // Ensure quads only

// Boundary layer around the cylinder
Field[1] = BoundaryLayer;
Field[1].EdgesList = {1,2};
Field[1].hfar = 0.01;
Field[1].hwall_n = 0.001;
Field[1].thickness = 0.01;
Field[1].ratio = 1.1;
Field[1].IntersectMetrics = 1;
Field[1].Quads = 1;
BoundaryLayer Field = 1;

// Adding boundary layer for inlet and outlet
Field[2] = BoundaryLayer;
Field[2].EdgesList = {3,5}; // Apply to inlet and outlet
Field[2].hfar = 0.01;
Field[2].hwall_n = 0.001;
Field[2].thickness = 0.005;
Field[2].ratio = 1.1;
Field[2].IntersectMetrics = 1;
Field[2].Quads = 1;
BoundaryLayer Field = 2;

// Extrusion from top boundary
// Extrude{0,0.05,0} {
//     Line{4}; Layers{5}; Recombine;}

// Extrusion from bottom boundary
// Extrude{0,-0.05,0} {
//     Line{6}; Layers{5}; Recombine;}

Mesh.Algorithm = 6;
Recombine Surface {1};
