#include "Incompressible_Solver.h"

// Global variables implementation for standalone build
// These are normally defined in Grid_Computations.cpp

#include "Incompressible_Solver_Standalone.h"

// Global variables implementation for standalone build
// These are normally defined in Grid_Computations.cpp

// Global variables definitions (normally from Grid_Computations.cpp)
std::vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;
std::vector<int> Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
std::vector<double> Vertices;
int Total_No_Cells = 0, No_Physical_Cells = 0;
bool Is_2D_Flow = true;