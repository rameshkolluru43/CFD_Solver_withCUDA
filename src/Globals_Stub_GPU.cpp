#include "Globals.h"

// Define global containers required by GPU main stub
vector<V_D> Cells_DelU, U_Cells, Primitive_Cells, Cells_Net_Flux, R_Cell, Cells_Viscous_Flux;
vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;
vector<V_D> U_Cells_RK_1, U_Cells_RK_2, A_x, A_y, A, A_x_L, A_x_R, A_y_T, A_y_B;
V_D Error, b, Global_Primitive, Global_U, Average_Convective_Flux, Dissipative_Flux, Vertices;
V_D d_U, d_F, Mod_Alpha;
vector<vector<bool>> Cells_Face_Boundary_Type;

double Limiter_Zeta = 0.0, Limiter_Zeta1 = 0.0;
int nx_1 = 0, nx_2 = 0, ny_1 = 0, ny_2 = 0;
vector<string> gridFiles;
string gridDir, Test_Case_Name, GridVTKFile, Flow_Type, Test_Case_JSON_File, Test_Case_Config_File;
string Solver_Name, Description, Author, GeometryType;
int Solver_Type = 0;
bool Is_Conservative = false, Is_Viscous = false;
bool Is_Viscous_Wall = false, Is_2D_Flow = true, Is_Inlet_SubSonic = true, Is_Exit_SubSonic = true, Enable_Far_Field = false, Is_Second_Order = false;
bool Is_Implicit_Method = false, Is_MOVERS_1 = false, Enable_Entropy_Fix = false, Is_Time_Dependent = false, has_Symmetry_BC = false, Time_Accurate = false;
bool Local_Time_Stepping = false, Non_Dimensional_Form = false, Is_WENO = false, Is_Char = false;
std::string Solution_File = "gpu_solution.dat";
std::string Grid_File, Initial_Solution_File, Error_File, Limiter_File, Final_Solution_File, Grid_Vtk_File, CF_File, BCFileName, InitCondFileName;
int Total_No_Cells = 0, No_Cartesian_Cells = 0, No_Polar_Cells = 0, No_Physical_Cells = 10, Grid_Type = 0, Case_Type = 0, Exit_Type = 0, Inlet_Type = 0, Initialize_Type = 0, Procedure_Type = 0, Numerical_Method = 0, Method_Type = 0, Test_Case = 1;
int No_Ghost_Cells = 0, Cells_in_Plane = 0, nx_p = 0, ny_p = 0, nz_p = 0, nx_c = 0, ny_c = 0, nz_c = 0, iterations = 0, Area_Weighted_Average = 0, Total_Iterations = 5, Limiter_Case = 0;
int Dissipation_Type = 0, Flux_Type = 0, Grid_Size = 0, Viscous_Time_Case = 0, NUM_FLUX_COMPONENTS = 4;

double global_temp = 0.0, Pressure_Static_Exit = 0.0, Pressure_Total_Inlet = 0.0, Temperature_Total_Inlet = 0.0, Cell_Minimum_Length = 1.0, CFL = 0.5, Re = 0.0, Pr = 0.0, Inv_Pr = 0.0, Inv_Re = 0.0, Qx = 0.0, Qy = 0.0;
double Max_dt = 0.0, Min_dt = 0.0;
double u_ref = 1.0, t_ref = 1.0, L_ref = 1.0, M_ref = 1.0, mu_ref = 1.0, P_ref = 1.0, Rho_ref = 1.0, K1 = 0.0, M = 0.0, cp_ref = 0.0, R_ref = 0.0;
double Pressure_Static_Inlet = 0.0, Rho_Static_Inlet = 1.0, Temperature_Static_Inlet = 300.0, Inlet_Mach_No = 0.1, u = 0.0, v = 0.0, V_Magnitude_Inlet = 0.0, V_Magnitude = 0.0, C_Acoustic = 0.0, Terminating_Time = 0.0, Total_Time = 0.0;
double Mod_Alpha0 = 0.0, Mod_Alpha1 = 0.0, Mod_Alpha2 = 0.0, Mod_Alpha3 = 0.0, d_F_0 = 0.0, d_F_1 = 0.0, d_F_2 = 0.0, d_F_3 = 0.0, d_U_0 = 0.0, d_U_1 = 0.0, d_U_2 = 0.0, d_U_3 = 0.0;
double Lambda_Max = 0.0, Lambda_Min = 0.0, Max1 = 0.0, Max2 = 0.0, Min1 = 0.0, Min2 = 0.0, d_Var = 0.0;
double P_inf = 101325.0, Rho_inf = 1.0, R_inf = 0.0, M_inf = 0.0, u_inf = 0.0, T_inf = 300.0, mu_star = 0.0, q_inf = 0.0;
int Neighbour_1 = -1, Neighbour_2 = -1, Neighbour_3 = -1, Neighbour_4 = -1, Neighbour_5 = -1, Neighbour_6 = -1, Neighbour_7 = -1, Neighbour_8 = -1;
V_D u_Gradient, v_Gradient, T_Gradient, Rho_Gradient, P_Gradient, Viscous_Flux, Grad_Var;
double mu = 0.0, K = 0.0;

// Variables for face normals and areas
double nx = 0.0, ny = 0.0, dl = 0.0, dA = 0.0;

// Left and Right State Variables
double Rho_L = 0.0, T_L = 0.0, P_L = 0.0, u_L = 0.0, v_L = 0.0, Vdotn_L = 0.0, Vmag_L = 0.0, M_L = 0.0, C_L = 0.0, H_L = 0.0;
V_D Flux_L, U_L, CF;
double Rho_R = 0.0, T_R = 0.0, P_R = 0.0, u_R = 0.0, v_R = 0.0, Vdotn_R = 0.0, Vmag_R = 0.0, M_R = 0.0, C_R = 0.0, H_R = 0.0, Epsilon = 0.0;
V_D Flux_R, U_R, S;

// Roe average variables
double Roe_u = 0.0, Roe_v = 0.0, Roe_Rho = 0.0, Roe_P = 0.0, Roe_e = 0.0, Roe_Vmag = 0.0, U_avg = 0.0, Roe_a = 0.0, Roe_H = 0.0, Un = 0.0, Ut = 0.0;

// Eigenvalue variables
double mev_L = 0.0, mev_R = 0.0, max_eigen_value = 0.0;

// Cell lists for boundary conditions
V_I Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List, Far_Field_Out_Flow, Far_Field_InFlow;
vector<bool> v_bool;

// Sparse matrix data
vector<int> row_indices, col_indices;
vector<double> values;

// Grid reader variables
int numNodes = 0, nodeIndex = 0, PointCellType = 0, LineCellType = 0, TriangleCellType = 0, QuadrilateralCellType = 0, HexahedronCellType = 0, TetrahedronCellType = 0, WedgeCellType = 0;

// Boundary info structure
BoundaryInfo gBoundaryInfo;

// Boundary and geometry condition structures
InletCondition inletCond;
ExitCondition exitCond;
InitialCondition initCond;
WallCondition wallCond;
GeometryParams geomParams;
MeshParams meshParams;

// Simple stubs to satisfy linker for GPU-only prototype build
void readJSON(const std::string &) {}
void readTestCaseJSON(const std::string &) {}
void Write_Solution(std::string &, int) {}
