// script to generate gmsh grids

Point(1) = {0 , 0 , 0};
Point(2) = {10, 0, 0};
Point(3) = {10,10,0};
Point(4) = {0, 10, 0};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = { 4, 1}; 

Transfinite Line {1 , 2 , 3, 4} = 10;

Line Loop(1) ={ 1, 2, 3, 4};

Plane Surface(1) = {1};


//Recombine Surface {1} 0.4;