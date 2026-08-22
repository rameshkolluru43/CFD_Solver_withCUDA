#ifndef GLOBALS_H
#define GLOBALS_H

#include "definitions.h"
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

using namespace std;

/**
 * @file Globals.h
 * @brief Global variables and data structures for 3D CFD solver
 *
 * This file contains all global variables, data structures, and type definitions
 * used throughout the 3D CFD solver. It includes cell definitions, boundary
 * conditions, flow variables, and solver parameters extended for 3D computations.
 */

// Forward declarations
struct Cell;
struct BoundaryInfo;
struct InletCondition;
struct ExitCondition;
struct InitialCondition;
struct WallCondition;
struct GeometryParams;
struct MeshParams;
struct GeneralParams;

// Type definitions for convenience
typedef vector<double> V_D;
typedef vector<int> V_I;
typedef vector<vector<double>> VV_D;
typedef vector<vector<int>> VV_I;

// Fallback for vtkIdType if VTK disabled
#ifndef USE_VTK
using vtkIdType = long long;
#endif

/**
 * @struct hash_pair
 * @brief Hash function for pairs used in unordered containers
 */
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

/**
 * @struct BoundaryInfo
 * @brief Structure to store 3D boundary information
 *
 * Contains sets of boundary cell IDs for all six faces of a 3D domain:
 * left/right (x-direction), bottom/top (y-direction), back/front (z-direction),
 * plus special geometries like cylinders or spheres.
 */
struct BoundaryInfo
{
    std::unordered_set<vtkIdType> leftBoundary;     // x-min boundary
    std::unordered_set<vtkIdType> rightBoundary;    // x-max boundary
    std::unordered_set<vtkIdType> topBoundary;      // y-max boundary
    std::unordered_set<vtkIdType> bottomBoundary;   // y-min boundary
    std::unordered_set<vtkIdType> backBoundary;     // z-min boundary
    std::unordered_set<vtkIdType> frontBoundary;    // z-max boundary
    std::unordered_set<vtkIdType> cylinderBoundary; // Cylindrical/curved boundaries
    std::unordered_set<vtkIdType> sphereBoundary;   // Spherical boundaries
};

/**
 * @struct InletCondition
 * @brief Inlet boundary conditions for 3D flow
 *
 * Contains all physical properties at inlet boundary including
 * the additional w-velocity component for 3D flow.
 */
struct InletCondition
{
    std::string type;
    int testCase;
    double P;       // Pressure
    double Rho;     // Density
    double M;       // Mach number
    double u;       // x-velocity component
    double v;       // y-velocity component
    double w;       // z-velocity component (3D extension)
    double T;       // Temperature
    double T_inf;   // Reference temperature
    double P_inf;   // Reference pressure
    double Rho_inf; // Reference density
};

/**
 * @struct ExitCondition
 * @brief Exit boundary conditions for 3D flow
 */
struct ExitCondition
{
    std::string type;
    double P;   // Pressure
    double Rho; // Density
    double T;   // Temperature
    double M;   // Mach number
    double u;   // x-velocity component
    double v;   // y-velocity component
    double w;   // z-velocity component (3D extension)
};

/**
 * @struct InitialCondition
 * @brief Initial conditions for 3D flow field
 */
struct InitialCondition
{
    std::string type;
    double P;       // Pressure
    double Rho;     // Density
    double M;       // Mach number
    double u;       // x-velocity component
    double v;       // y-velocity component
    double w;       // z-velocity component (3D extension)
    double T;       // Temperature
    double T_inf;   // Reference temperature
    double P_inf;   // Reference pressure
    double Rho_inf; // Reference density
};

/**
 * @struct WallCondition
 * @brief Wall boundary conditions for 3D flow
 */
struct WallCondition
{
    std::string type;
    double M;   // Mach number
    double u;   // x-velocity component
    double v;   // y-velocity component
    double w;   // z-velocity component (3D extension)
    double T;   // Temperature
    double P;   // Pressure
    double Rho; // Density
};

/**
 * @struct GeometryParams
 * @brief Geometric parameters for 3D configurations
 */
struct GeometryParams
{
    double radius, length, width, height, thickness;
    double x_center, y_center, z_center; // 3D center coordinates
};

/**
 * @struct MeshParams
 * @brief 3D mesh parameters
 */
struct MeshParams
{
    int gridSize, nx, ny, nz; // Added nz for 3D
    std::string meshType;
    double dx, dy, dz; // Grid spacing in each direction
};

/**
 * @struct GeneralParams
 * @brief General solver parameters
 */
struct GeneralParams
{
    double L_ref, M_ref, mu_ref, P_ref, Rho_ref, cp_ref, R_ref;
    string Solver_Name, Description, Author, GeometryType;
    string Test_Case_Name;
    int version;
};

/**
 * @struct Cell
 * @brief 3D Cell structure for hexahedral elements
 *
 * Extended from 2D to handle 6 faces per cell with corresponding
 * areas, normals, and connectivity information for 3D geometries.
 */
struct Cell
{
    // Basic cell properties
    int cellType, cellID, Dimension, ParentCellID, NoBoundaryFaces;
    int numFaces, numNodes, numEdges;

    // Connectivity information
    vector<int> nodeIndices, Neighbours, faceID, Secondary_Neighbours;
    vector<int> edgeIndices; // Additional for 3D

    // Geometric properties
    V_D Diagonal_Vector;
    double Volume, Inv_Volume, del_t; // Changed from Area to Volume for 3D

    // Face-related arrays (extended to 6 faces for 3D)
    V_D Face_Areas;          // Size 6 for hexahedral cells
    V_D Face_Normals;        // Size 18 (3 components × 6 faces: nx,ny,nz for each face)
    V_D Face_Centers;        // Size 18 (3 coordinates × 6 faces)
    V_D Cell_Face_Distances; // Size 6 (distance from cell center to each face)

    // Cell center and vertices
    V_D Cell_Center;        // Size 3 (x, y, z coordinates)
    V_D Cell_Center_Vector; // Size 3
    V_D Cell_Center_Distances;
    V_D Cell_Vertices; // Size 24 (3 coordinates × 8 vertices for hexahedron)
    V_D Cell_Areas;    // For surface area calculations

    // Boundary face flags (extended for 3D)
    bool hasBoundaryface = false;
    bool Is_Splittable = false;
    bool has_Wall_Face = false;
    bool has_Inlet_Face = false;
    bool has_Exit_Face = false;
    bool has_Symmetry_Face = false;

    // 3D boundary face identification
    bool Left_Face = false;   // x-min face
    bool Right_Face = false;  // x-max face
    bool Top_Face = false;    // y-max face
    bool Bottom_Face = false; // y-min face
    bool Back_Face = false;   // z-min face (3D extension)
    bool Front_Face = false;  // z-max face (3D extension)
    bool Interior_Face = false;

    // Constructors
    Cell()
    {
        numFaces = NUM_FACES_3D; // 6 faces for hexahedral cells
        Face_Areas.resize(NUM_FACES_3D, 0.0);
        Face_Normals.resize(NUM_FACES_3D * 3, 0.0); // 3 components per face
        Face_Centers.resize(NUM_FACES_3D * 3, 0.0);
        Cell_Face_Distances.resize(NUM_FACES_3D, 0.0);
        Cell_Center.resize(3, 0.0);
        Cell_Center_Vector.resize(3, 0.0);
    }

    Cell(int nNodes, int id, const vector<int> &indices)
        : numNodes(nNodes), cellID(id), nodeIndices(indices)
    {
        numFaces = NUM_FACES_3D;
        Face_Areas.resize(NUM_FACES_3D, 0.0);
        Face_Normals.resize(NUM_FACES_3D * 3, 0.0);
        Face_Centers.resize(NUM_FACES_3D * 3, 0.0);
        Cell_Face_Distances.resize(NUM_FACES_3D, 0.0);
        Cell_Center.resize(3, 0.0);
        Cell_Center_Vector.resize(3, 0.0);
    }
};

/*------------------------------------------Global Variables--------------------------------------------------------*/

// Boundary information
extern BoundaryInfo gBoundaryInfo;

// Condition structures
extern InletCondition inletCond;
extern ExitCondition exitCond;
extern InitialCondition initCond;
extern WallCondition wallCond;
extern GeometryParams geomParams;
extern MeshParams meshParams;
extern GeneralParams generalParams;

// Cell arrays
extern vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;

// 3D arrays for storing cell data (extended for 5 conservative variables in 3D)
extern vector<V_D> Cells_DelU, U_Cells, Primitive_Cells, Cells_Net_Flux, R_Cell, Cells_Viscous_Flux;
extern vector<V_D> U_Cells_RK_1, U_Cells_RK_2, U_Cells_RK_3, U_Cells_RK_4; // Extended RK stages

// 3D flux arrays
extern vector<V_D> A_x, A_y, A_z;                            // Added A_z for 3D
extern vector<V_D> A_x_L, A_x_R, A_y_T, A_y_B, A_z_F, A_z_B; // Added z-direction fluxes

// Global vectors
extern V_D Error, b, Global_Primitive, Global_U, Average_Convective_Flux, Dissipative_Flux, Vertices;
extern V_D d_U, d_F, Mod_Alpha;
extern vector<vector<bool>> Cells_Face_Boundary_Type;

// Grid parameters (extended for 3D)
extern int nx_1, nx_2, ny_1, ny_2, nz_1, nz_2; // Added z-direction grid parameters
extern int nx_c, ny_c, nz_c;                   // Cell counts in each direction
extern int nx_p, ny_p, nz_p;                   // Point counts in each direction
extern int Total_No_Cells, No_Cartesian_Cells, No_Polar_Cells, No_Physical_Cells;
extern int No_Ghost_Cells, Cells_in_Volume; // Changed from Cells_in_Plane to Cells_in_Volume

// Solver parameters
extern double Limiter_Zeta, Limiter_Zeta1;
extern vector<string> gridFiles;
extern string gridDir, Test_Case_Name, GridVTKFile, Flow_Type;
extern string Test_Case_JSON_File, Test_Case_Config_File;
extern int Grid_Type, Case_Type, Exit_Type, Inlet_Type, Initialize_Type;
extern int Procedure_Type, Numerical_Method, Method_Type, Test_Case;
extern int Area_Weighted_Average, Total_Iterations, Limiter_Case;
extern int iterations;

// Physical parameters
extern double global_temp, Pressure_Static_Exit, Pressure_Total_Inlet, Temperature_Total_Inlet;
extern double Cell_Minimum_Length, CFL, Re, Pr, Inv_Pr, Inv_Re;
extern double Qx, Qy, Qz; // Added heat flux in z-direction

// Time stepping
extern double Max_dt, Min_dt;
extern double u_ref, t_ref, L_ref, M_ref, mu_ref, P_ref, Rho_ref, K1, M, cp_ref, R_ref;
extern double Pressure_Static_Inlet, Rho_Static_Inlet, Temperature_Static_Inlet;
extern double Inlet_Mach_No, u, v, w; // Added w-velocity component
extern double V_Magnitude_Inlet, V_Magnitude, C_Acoustic, Terminating_Time, Total_Time;

// Modified alpha variables (extended for 3D)
extern double Mod_Alpha0, Mod_Alpha1, Mod_Alpha2, Mod_Alpha3, Mod_Alpha4; // Added 4th component
extern double d_F_0, d_F_1, d_F_2, d_F_3, d_F_4;                          // Added 4th flux component
extern double d_U_0, d_U_1, d_U_2, d_U_3, d_U_4;                          // Added 4th conservative variable

// Eigenvalue and flux variables
extern double Lambda_Max, Lambda_Min, Max1, Max2, Max3, Min1, Min2, Min3, d_Var;        // Added Max3, Min3
extern double P_inf, Rho_inf, R_inf, M_inf, u_inf, v_inf, w_inf, T_inf, mu_star, q_inf; // Added v_inf, w_inf

// Neighbor indices (extended for 3D hexahedral cells - up to 26 neighbors)
extern int Neighbour_1, Neighbour_2, Neighbour_3, Neighbour_4, Neighbour_5, Neighbour_6;
extern int Neighbour_7, Neighbour_8, Neighbour_9, Neighbour_10, Neighbour_11, Neighbour_12;

// Gradient variables (extended for 3D)
extern V_D u_Gradient, v_Gradient, w_Gradient, T_Gradient, Rho_Gradient, P_Gradient; // Added w_Gradient
extern V_D Viscous_Flux, Grad_Var;

// Material properties
extern double mu, K, M_inf;

// 3D face normal and area variables
extern double nx, ny, nz, dl, dA, dV; // Added nz and dV (volume element)

// Left state variables (extended for 3D)
extern double Rho_L, T_L, P_L, u_L, v_L, w_L, Vdotn_L, Vmag_L, M_L, C_L, H_L; // Added w_L
extern V_D Flux_L, U_L, CF;

// Right state variables (extended for 3D)
extern double Rho_R, T_R, P_R, u_R, v_R, w_R, Vdotn_R, Vmag_R, M_R, C_R, H_R, Epsilon; // Added w_R
extern V_D Flux_R, U_R, S;

// Boundary cell lists
extern V_I Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
extern V_I Far_Field_Out_Flow, Far_Field_InFlow;
extern V_I Back_Cells_List, Front_Cells_List; // Added for 3D boundaries

// Boolean and control variables
extern vector<bool> v_bool;
extern int Dissipation_Type, Flux_Type, Grid_Size, Viscous_Time_Case, NUM_FLUX_COMPONENTS;
extern int Solver_Type;
extern string Solver_Name, Description, Author, GeometryType;

// Boolean flags
extern bool Is_Viscous_Wall, Is_3D_Flow, Is_Inlet_SubSonic, Is_Exit_SubSonic; // Changed from Is_2D_Flow to Is_3D_Flow
extern bool Enable_Far_Field, Is_Second_Order, Is_Implicit_Method, Is_MOVERS_1;
extern bool Enable_Entropy_Fix, Is_Time_Dependent, has_Symmetry_BC, Time_Accurate;
extern bool Local_Time_Stepping, Non_Dimensional_Form, Is_WENO, Is_Char;
extern bool Is_Conservative, Is_Viscous;

// File names
extern string Grid_File, Initial_Solution_File, Solution_File, Error_File, Limiter_File;
extern string Final_Solution_File, Grid_Vtk_File, CF_File, BCFileName, InitCondFileName;

// Roe averaging variables (extended for 3D)
extern double mev_L, mev_R, max_eigen_value;
extern double Roe_u, Roe_v, Roe_w, Roe_Rho, Roe_P, Roe_e, Roe_Vmag; // Added Roe_w
extern double U_avg, Roe_a, Roe_H, Un, Ut, Uz;                      // Added Uz component

// Sparse matrix data
extern vector<int> row_indices, col_indices;
extern vector<double> values;

// Mesh reader variables (extended for 3D cell types)
extern int numNodes, nodeIndex;
extern int PointCellType, LineCellType, TriangleCellType, QuadrilateralCellType;
extern int TetrahedronCellType, HexahedronCellType, WedgeCellType, PyramidCellType, PrismCellType;

#endif // GLOBALS_H