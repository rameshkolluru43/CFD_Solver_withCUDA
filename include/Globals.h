#ifndef GLOBALS_H
#define GLOBALS_H

#include "definitions.h"

/**
 * @file definitions.h
 * @brief Header file containing global constants, structures, and function declarations for the solver.
 *
 * This file includes all necessary libraries, defines global constants, declares global variables,
 * and provides function prototypes required for the solver. It also contains detailed documentation
 * for the structures and functions used in the computational domain and boundary condition handling.
 *
 * @section Libraries
 * Includes standard C++ libraries, VTK libraries, Boost, and JSONCPP for various functionalities.
 *
 * @section Global Constants
 * Defines mathematical constants, physical constants, and other macros used throughout the code.
 *
 * @section Structures
 * - `BoundaryInfo`: Stores boundary information for different regions.
 * - `InletCondition`: Represents inlet boundary conditions.
 * - `ExitCondition`: Represents exit boundary conditions.
 * - `InitialCondition`: Represents initial conditions for the simulation.
 * - `Cell`: Represents a computational cell with properties, connectivity, and geometric attributes.
 *
 * @section Global Variables
 * Declares numerous global variables for storing simulation data, grid information, solver parameters,
 * and boundary conditions.
 *
 * @section Functions
 * - Grid and Cell Construction:
 *   - Functions for reading grid files, constructing cells, identifying neighbors, and computing centroids.
 * - Solver Functions:
 *   - Functions for solving inviscid and viscous flows, updating variables, and evaluating fluxes.
 * - Boundary Conditions:
 *   - Functions for applying various boundary conditions (e.g., wall, inlet, exit, symmetry).
 * - Utility Functions:
 *   - Functions for vector operations, sorting points, and evaluating gradients.
 * - Output Functions:
 *   - Functions for writing solution data, error files, and VTK files.
 * - Test Case Functions:
 *   - Functions for creating test case directories and handling specific test cases.
 *
 * @section Notes
 * - The file uses extensive use of `extern` to declare global variables and structures.
 * - The file is designed to be included in multiple source files for modularity.
 * - Ensure proper initialization of global variables before using them in the solver.
 */
/*@brief Contains the definition of the Cell structure used in the solver.
        * /

    /*
     * @struct Cell
     * @brief Represents a computational cell in the solver.
     *
     * The Cell structure contains information about the cell's properties,
     * connectivity, and geometric attributes. It is used to store and manage
     * data related to the computational domain.
     *
     * @var Cell::cellType
     * Type of the cell (e.g., triangular, quadrilateral, etc.).
     *
     * @var Cell::cellID
     * Unique identifier for the cell.
     *
     * @var Cell::Dimension
     * Dimensionality of the cell (e.g., 2D, 3D).
     *
     * @var Cell::ParentCellID
     * Identifier of the parent cell, if applicable.
     *
     * @var Cell::nodeIndices
     * Indices of the nodes that define the cell.
     *
     * @var Cell::Neighbours
     * Indices of neighboring cells.
     *
     * @var Cell::faceID
     * Indices of the faces associated with the cell.
     *
     * @var Cell::Secondary_Neighbours
     * Indices of secondary neighboring cells.
     *
     * @var Cell::hasBoundaryface
     * Indicates whether the cell has a boundary face.
     *
     * @var Cell::Is_Splittable
     * Indicates whether the cell can be split.
     *
     * @var Cell::Area
     * Area of the cell.
     *
     * @var Cell::Inv_Area
     * Inverse of the cell area.
     *
     * @var Cell::Face_Areas
     * Vector of areas of the cell's faces.
     *
     * @var Cell::Face_Normals
     * Vector of normals of the cell's faces.
     *
     * @var Cell::Cell_Center
     * Coordinates of the cell's center.
     *
     * @var Cell::Cell_Center_Distances
     * Distances from the cell center to other relevant points.
     *
     * @var Cell::Cell_Vertices
     * Coordinates of the vertices of the cell.
     *
     * @var Cell::Cell_Face_Distances
     * Distances from the cell center to the faces.
     *
     * @var Cell::Cell_Areas
     * Vector of areas associated with the cell.
     *
     * @brief Constructor to initialize a Cell object.
     *
     * @fn Cell::Cell()
     * Default constructor that initializes an empty Cell object.
     *
     * @fn Cell::Cell(int numNodes, int id, const vector<int> &indices)
     * Parameterized constructor to initialize a Cell object with specific
     * attributes.
     *
     * @param numNodes Number of nodes in the cell.
     * @param id Unique identifier for the cell.
     * @param indices Indices of the nodes that define the cell.
     */
/*
 * @struct InitialCondition
 * @brief Represents the initial conditions for a simulation.
 *
 * This structure contains various physical properties and parameters
 * that define the initial state of the system being simulated.
 *
 * @var InitialCondition::P
 * Pressure at the initial condition.
 *
 * @var InitialCondition::Rho
 * Density at the initial condition.
 *
 * @var InitialCondition::M
 * Mach number at the initial condition.
 *
 * @var InitialCondition::u
 * Velocity component in the x-direction at the initial condition.
 *
 * @var InitialCondition::v
 * Velocity component in the y-direction at the initial condition.
 *
 * @var InitialCondition::T
 * Temperature at the initial condition.
 *
 * @var InitialCondition::T_inf
 * Reference temperature (ambient) at the initial condition.
 *
 * @var InitialCondition::P_inf
 * Reference pressure (ambient) at the initial condition.
 *
 * @var InitialCondition::Rho_inf
 * Reference density (ambient) at the initial condition.
 */

/*
 * @struct InletCondition
 * @brief Represents the inlet boundary conditions for a simulation.
 *
 * This structure contains various physical properties and parameters
 * that define the state of the system at the inlet boundary.
 *
 * @var InletCondition::testCase
 * Identifier for the test case being simulated.
 *
 * @var InletCondition::P
 * Pressure at the inlet boundary.
 *
 * @var InletCondition::Rho
 * Density at the inlet boundary.
 *
 * @var InletCondition::M
 * Mach number at the inlet boundary.
 *
 * @var InletCondition::u
 * Velocity component in the x-direction at the inlet boundary.
 *
 * @var InletCondition::v
 * Velocity component in the y-direction at the inlet boundary.
 *
 * @var InletCondition::T
 * Temperature at the inlet boundary.
 *
 * @var InletCondition::T_inf
 * Reference temperature (ambient) at the inlet boundary.
 *
 * @var InletCondition::P_inf
 * Reference pressure (ambient) at the inlet boundary.
 *
 * @var InletCondition::Rho_inf
 * Reference density (ambient) at the inlet boundary.
 */

// vector<V_D> represents a 2D array such as a[i][j]

struct hash_pair
{
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2> &p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1); // Bitwise XOR and shift to avoid collisions
    }
};

/*------------------------------------------Global Variables--------------------------------------------------------*/
// extern string pwd;
typedef vector<double> V_D;
typedef vector<int> V_I;

// Struct to store boundary information
// Fallback for vtkIdType if VTK disabled
#ifndef USE_VTK
using vtkIdType = long long;
#endif

#ifndef __CUDACC__
struct BoundaryInfo
{
    std::unordered_set<vtkIdType> leftBoundary;
    std::unordered_set<vtkIdType> rightBoundary;
    std::unordered_set<vtkIdType> topBoundary;
    std::unordered_set<vtkIdType> bottomBoundary;
    std::unordered_set<vtkIdType> cylinderBoundary;
};
extern BoundaryInfo gBoundaryInfo;
#else
struct BoundaryInfo {};
#endif

// Structures to hold boundary conditions
struct InletCondition
{
    std::string type;
    int testCase;
    double P;
    double Rho;
    double M;
    double u;
    double v;
    double T;
    double T_inf;
    double P_inf;
    double Rho_inf;
};

struct ExitCondition
{
    std::string type;
    double P;
    double Rho;
    double T;
    double M;
    double u;
    double v;
    // Extend with additional fields as needed
};

struct InitialCondition
{
    std::string type;
    double P;
    double Rho;
    double M;
    double u;
    double v;
    double T;
    double T_inf;
    double P_inf;
    double Rho_inf;
};

struct WallCondition
{
    std::string type;
    double M;
    double u;
    double v;
    double T;
    double P;
    double Rho;
};

struct GeometryParams
{
    double radius, length, thickness;
};

struct MeshParams
{
    int gridSize, nx, ny;
    std::string meshType;
};
extern MeshParams meshParams;
struct GeneralParams
{
    double L_ref, M_ref, mu_ref, P_ref, Rho_ref, cp_ref, R_ref;
    string Solver_Name, Description, Author, GeometryType;
    string Test_Case_Name;
    int version;
};
struct Cell
{
    int cellType, cellID, Dimension, ParentCellID, NoBoundaryFaces, numFaces, numNodes;
    vector<int> nodeIndices, Neighbours, faceID, Secondary_Neighbours;
    V_D Diagonal_Vector;
    bool hasBoundaryface = false, Is_Splittable = false, has_Wall_Face = false, has_Inlet_Face = false, has_Exit_Face = false, has_Symmetry_Face = false;
    double Area, Inv_Area, del_t;
    V_D Face_Areas, Face_Normals, Cell_Center, Cell_Center_Distances, Cell_Vertices, Cell_Face_Distances, Cell_Areas, Cell_Center_Vector;
    bool Left_Face = false, Right_Face = false, Top_Face = false, Bottom_Face = false, Interior_Face = false;
    /** Per-face boundary kind for mixed tri/quad: 0=INTERNAL, 1=LEFT, 2=RIGHT, 3=TOP, 4=BOTTOM, 5=WALL. Size = numFaces. */
    vector<int> Face_Boundary_Kind;
    // Constructor to initialize Cell
    Cell() {}
    Cell(int nNodes, int id, const vector<int> &indices)
        : numNodes(nNodes), cellID(id), nodeIndices(indices) {}
};

// Header file to store all the global variables used in the code
//  2D array for storing Cell Face Normals, Area Componets,Conservative variables and Convective Flux(Fc), Viscous Flux (Fv), Primitive Variables and Net Flux is (Fc+Fv)
extern vector<V_D> Cells_DelU, U_Cells, Primitive_Cells, Cells_Net_Flux, R_Cell, Cells_Viscous_Flux;

extern vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;

extern vector<V_D> U_Cells_RK_1, U_Cells_RK_2, A_x, A_y, A, A_x_L, A_x_R, A_y_T, A_y_B;

extern V_D Error, b, Global_Primitive, Global_U, Average_Convective_Flux, Dissipative_Flux, Vertices;
extern V_D d_U, d_F, Mod_Alpha;
extern vector<vector<bool>> Cells_Face_Boundary_Type;
extern double Limiter_Zeta, Limiter_Zeta1;
extern int nx_1, nx_2, ny_1, ny_2;
extern vector<string> gridFiles;
extern string gridDir, Test_Case_Name, GridVTKFile, Flow_Type, Test_Case_JSON_File, Test_Case_Config_File;
extern int Total_No_Cells, No_Cartesian_Cells, No_Polar_Cells, No_Physical_Cells, Grid_Type, Case_Type, Exit_Type, Inlet_Type, Initialize_Type, Procedure_Type, Numerical_Method, Method_Type, Test_Case;

extern int No_Ghost_Cells, Cells_in_Plane, nx_p, ny_p, nz_p, nx_c, ny_c, nz_c, iterations, Area_Weighted_Average, Total_Iterations, Limiter_Case;

extern double global_temp, Pressure_Static_Exit, Pressure_Total_Inlet, Temperature_Total_Inlet, Cell_Minimum_Length, CFL, Re, Pr, Inv_Pr, Inv_Re, Qx, Qy;

// Variables for storing the minimum and maximum values of timesetp
extern double Max_dt, Min_dt;
// reference Variables for non-dimensionalization in NS equations
extern double u_ref, t_ref, L_ref, M_ref, mu_ref, P_ref, Rho_ref, K1, M, cp_ref, R_ref;

extern double Pressure_Static_Inlet, Rho_Static_Inlet, Temperature_Static_Inlet, Inlet_Mach_No, u, v, V_Magnitude_Inlet, V_Magnitude, C_Acoustic, Terminating_Time, Total_Time;

extern double Mod_Alpha0, Mod_Alpha1, Mod_Alpha2, Mod_Alpha3, d_F_0, d_F_1, d_F_2, d_F_3, d_U_0, d_U_1, d_U_2, d_U_3;

extern double Lambda_Max, Lambda_Min, Max1, Max2, Min1, Min2, d_Var;
// Free Stream Varialbles
extern double P_inf, Rho_inf, R_inf, M_inf, u_inf, T_inf, mu_star, q_inf;
// Indicies for Neighbouring grid points of a given cell

extern int Neighbour_1, Neighbour_2, Neighbour_3, Neighbour_4, Neighbour_5, Neighbour_6, Neighbour_7, Neighbour_8;

// Variables used for Viscous Fluxes and Gradients
extern V_D u_Gradient, v_Gradient, T_Gradient, Rho_Gradient, P_Gradient, Viscous_Flux, Grad_Var;

// Variables for viscosity and thermal Conductivity
extern double mu, K, M_inf;

// variables to read the face normals and face lenghts and face area in 2D
extern double nx, ny, dl, dA;

// Left State Variables
extern double Rho_L, T_L, P_L, u_L, v_L, Vdotn_L, Vmag_L, M_L, C_L, H_L;
extern V_D Flux_L, U_L, CF;

// Right State Varialbes
extern double Rho_R, T_R, P_R, u_R, v_R, Vdotn_R, Vmag_R, M_R, C_R, H_R, Epsilon;
extern V_D Flux_R, U_R, S;

// Variables for storing the cell indicies of the wall, inlet, exit and symmetry cells
extern V_I Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List, Far_Field_Out_Flow, Far_Field_InFlow;

extern vector<bool> v_bool;

extern int Dissipation_Type, Flux_Type, Grid_Size, Viscous_Time_Case, NUM_FLUX_COMPONENTS;
extern int Solver_Type;
extern string Solver_Name, Description, Author, GeometryType;
// boolean variables for tagging various conditions
extern bool Is_Viscous_Wall, Is_2D_Flow, Is_Inlet_SubSonic, Is_Exit_SubSonic, Enable_Far_Field, Is_Second_Order;
extern bool Is_Implicit_Method, Is_MOVERS_1, Enable_Entropy_Fix, Is_Time_Dependent, has_Symmetry_BC, Time_Accurate;
extern bool Enable_AMR;
extern int AMR_Period;
extern double AMR_Gradient_Threshold, AMR_Max_Fraction;
extern vector<double> Gradient_Refinement_Indicator;
extern bool Local_Time_Stepping, Non_Dimensional_Form, Is_WENO, Is_Char, Is_Conservative, Is_Viscous;

extern string Grid_File, Initial_Solution_File, Solution_File, Error_File, Limiter_File, Final_Solution_File, Grid_Vtk_File, CF_File, BCFileName, InitCondFileName;

extern double mev_L, mev_R, max_eigen_value;
extern InletCondition inletCond;
extern ExitCondition exitCond;
extern InitialCondition initCond;

//	Variables for Roe Averages
extern double Roe_u, Roe_v, Roe_Rho, Roe_P, Roe_e, Roe_Vmag, U_avg, Roe_a, Roe_H, Un, Ut;

// Define vectors to store sparse matrix data
extern vector<int> row_indices;
extern vector<int> col_indices;
extern vector<double> values;
// Define extern variables used in gmsh grid reader
extern int numNodes, nodeIndex, PointCellType, LineCellType, TriangleCellType, QuadrilateralCellType, HexahedronCellType, TetrahedronCellType, WedgeCellType;

#endif // #ifndef GLOBALS_H