/**
 * @file LES_Models.h
 * @brief Header file for Large Eddy Simulation (LES) subgrid-scale models
 * @author CFD Solver Team
 * @date 2026
 *
 * This file contains declarations for LES subgrid-scale (SGS) models:
 * - Smagorinsky model (classic)
 * - Dynamic Smagorinsky model (Germano et al.)
 * - WALE model (Wall-Adapting Local Eddy-viscosity)
 * - Vreman model
 * - Sigma model
 *
 * Supports both CPU and GPU (CUDA) implementations.
 */

#ifndef LES_MODELS_H
#define LES_MODELS_H

#include "definitions.h"
#include "Globals.h"

//=============================================================================
// LES MODEL ENUMERATIONS
//=============================================================================

enum class LES_SGS_Model
{
    NONE = 0,
    SMAGORINSKY = 1,
    DYNAMIC_SMAGORINSKY = 2,
    WALE = 3,
    VREMAN = 4,
    SIGMA = 5
};

//=============================================================================
// LES MODEL CONSTANTS
//=============================================================================

namespace LES_Constants
{
    // Smagorinsky model constants
    constexpr double Cs_default = 0.1;      // Smagorinsky constant (typical: 0.1-0.2)
    constexpr double Cs_channel = 0.1;       // For channel flow
    constexpr double Cs_isotropic = 0.17;    // For isotropic turbulence
    constexpr double Cs_mixing_layer = 0.12; // For mixing layers

    // Dynamic Smagorinsky constants
    constexpr double Cs_min = 0.0;           // Minimum Cs (backscatter clipping)
    constexpr double Cs_max = 0.3;           // Maximum Cs (stability limit)
    constexpr double test_filter_ratio = 2.0; // Test filter to grid filter ratio

    // WALE model constants
    constexpr double Cw_default = 0.5;       // WALE constant (typical: 0.5)
    constexpr double Cw_squared = 0.25;      // Cw^2

    // Vreman model constants
    constexpr double Cv_default = 0.07;      // Vreman constant (typical: 0.07)

    // Sigma model constants
    constexpr double Csigma_default = 1.35;  // Sigma model constant

    // Common constants
    constexpr double Pr_t = 0.9;             // Turbulent Prandtl number
    constexpr double Sc_t = 0.7;             // Turbulent Schmidt number
    constexpr double kappa = 0.41;           // Von Karman constant
    constexpr double A_plus = 26.0;          // Van Driest damping constant
}

//=============================================================================
// LES VARIABLES STRUCTURE
//=============================================================================

struct LES_Variables
{
    double nu_sgs;           // Subgrid-scale kinematic viscosity
    double mu_sgs;           // Subgrid-scale dynamic viscosity
    double k_sgs;            // Subgrid-scale kinetic energy (for some models)
    double delta;            // Filter width (grid scale)
    double S_mag;            // Strain rate magnitude |S|
    double Omega_mag;        // Vorticity magnitude |Omega|
    
    // Dynamic Smagorinsky specific
    double Cs_dynamic;       // Dynamically computed Cs
    double Lij[3][3];        // Leonard stress tensor
    double Mij[3][3];        // Resolved stress tensor
    
    // WALE specific
    double Sd_ij[3][3];      // Traceless symmetric part of g_ij^2
    double Sd_mag;           // |Sd|
    
    // Vreman specific
    double B_beta;           // Vreman B parameter
    double alpha_ij[3][3];   // Velocity gradient tensor
    
    // Sigma specific
    double sigma_1, sigma_2, sigma_3;  // Singular values of velocity gradient
    
    // Effective properties
    double mu_effective;     // Total effective viscosity (laminar + SGS)
    double k_effective;      // Total effective thermal conductivity
    
    // Wall damping
    double y_plus;           // Dimensionless wall distance
    double f_damping;        // Van Driest damping function
    
    // Statistics
    double Q_criterion;      // Q-criterion for vortex identification
    double lambda2;          // Lambda-2 criterion
    
    LES_Variables() : nu_sgs(0.0), mu_sgs(0.0), k_sgs(0.0), delta(0.0),
                      S_mag(0.0), Omega_mag(0.0), Cs_dynamic(0.1),
                      Sd_mag(0.0), B_beta(0.0),
                      sigma_1(0.0), sigma_2(0.0), sigma_3(0.0),
                      mu_effective(0.0), k_effective(0.0),
                      y_plus(0.0), f_damping(1.0),
                      Q_criterion(0.0), lambda2(0.0)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                Lij[i][j] = 0.0;
                Mij[i][j] = 0.0;
                Sd_ij[i][j] = 0.0;
                alpha_ij[i][j] = 0.0;
            }
    }
};

//=============================================================================
// GLOBAL LES VARIABLES (extern declarations)
//=============================================================================

extern LES_SGS_Model current_les_model;
extern vector<LES_Variables> les_vars;
extern bool Is_LES;
extern bool Use_Wall_Damping;
extern bool Use_GPU_LES;
extern double Cs_Smagorinsky;
extern double Cw_WALE;
extern double Cv_Vreman;
extern double Csigma;
extern double Pr_turb_LES;

//=============================================================================
// INITIALIZATION FUNCTIONS
//=============================================================================

void Initialize_LES_Model(LES_SGS_Model model);
void Initialize_LES_Variables();
void Set_LES_Filter_Width();
void Calculate_LES_Filter_Width(int cell_index);

//=============================================================================
// SMAGORINSKY MODEL FUNCTIONS
//=============================================================================

void Calculate_Smagorinsky_Viscosity(int cell_index);
void Calculate_Smagorinsky_Viscosity_All();
double Smagorinsky_Nu_SGS(double delta, double S_mag, double Cs);

//=============================================================================
// DYNAMIC SMAGORINSKY MODEL FUNCTIONS
//=============================================================================

void Calculate_Dynamic_Smagorinsky_Viscosity(int cell_index);
void Calculate_Dynamic_Smagorinsky_All();
void Apply_Test_Filter(int cell_index, const double *field, double *filtered);
void Calculate_Leonard_Stress(int cell_index);
void Calculate_Resolved_Stress(int cell_index);
double Calculate_Dynamic_Cs(int cell_index);
void Average_Dynamic_Cs_Homogeneous();

//=============================================================================
// WALE MODEL FUNCTIONS
//=============================================================================

void Calculate_WALE_Viscosity(int cell_index);
void Calculate_WALE_Viscosity_All();
void Calculate_Sd_Tensor(int cell_index);
double WALE_Nu_SGS(double delta, double S_mag, double Sd_mag, double Cw);

//=============================================================================
// VREMAN MODEL FUNCTIONS
//=============================================================================

void Calculate_Vreman_Viscosity(int cell_index);
void Calculate_Vreman_Viscosity_All();
void Calculate_Vreman_B(int cell_index);
double Vreman_Nu_SGS(double delta, double B_beta, double alpha_sq, double Cv);

//=============================================================================
// SIGMA MODEL FUNCTIONS
//=============================================================================

void Calculate_Sigma_Viscosity(int cell_index);
void Calculate_Sigma_Viscosity_All();
void Calculate_Singular_Values(int cell_index);
double Sigma_Nu_SGS(double delta, double sigma1, double sigma2, double sigma3, double Csigma);

//=============================================================================
// COMMON UTILITY FUNCTIONS
//=============================================================================

// Velocity gradient computations
void Calculate_Velocity_Gradient_Tensor(int cell_index, double G[3][3]);
void Calculate_Strain_Rate_Tensor_LES(int cell_index, double S[3][3]);
void Calculate_Rotation_Rate_Tensor(int cell_index, double Omega[3][3]);
double Calculate_Strain_Rate_Magnitude_LES(const double S[3][3]);
double Calculate_Vorticity_Magnitude_LES(const double Omega[3][3]);

// Filter width calculations
double Calculate_Filter_Width_Cube_Root(int cell_index);
double Calculate_Filter_Width_Max_Edge(int cell_index);
double Calculate_Filter_Width_Face_Area(int cell_index);

// Wall damping functions
void Calculate_Van_Driest_Damping(int cell_index);
double Van_Driest_Function(double y_plus);
void Calculate_Wall_Distance_LES(int cell_index);
void Calculate_Wall_Distances_LES_All();

// Effective viscosity
void Update_LES_Effective_Viscosity();
void Update_LES_Effective_Conductivity();
double Get_Effective_Viscosity_LES(int cell_index);
double Get_Effective_Conductivity_LES(int cell_index);

//=============================================================================
// SGS STRESS TENSOR FUNCTIONS
//=============================================================================

void Calculate_SGS_Stress_Tensor(int cell_index, double tau_sgs[3][3]);
void Calculate_SGS_Heat_Flux(int cell_index, double q_sgs[3]);
void Add_SGS_Stress_To_Viscous_Flux(int cell_index);

//=============================================================================
// FILTERING OPERATIONS
//=============================================================================

void Apply_Box_Filter(int cell_index, const vector<double> &field, double &filtered);
void Apply_Gaussian_Filter(int cell_index, const vector<double> &field, double &filtered);
void Apply_Top_Hat_Filter(int cell_index, const vector<double> &field, double &filtered);
void Filter_Velocity_Field(vector<double> &u_filtered, vector<double> &v_filtered);

//=============================================================================
// VORTEX IDENTIFICATION
//=============================================================================

void Calculate_Q_Criterion(int cell_index);
void Calculate_Lambda2_Criterion(int cell_index);
void Calculate_Delta_Criterion(int cell_index);
void Identify_Vortex_Structures();

//=============================================================================
// STATISTICS AND OUTPUT
//=============================================================================

void Calculate_LES_Statistics();
void Write_LES_Variables(const string &filename);
void Write_Instantaneous_Field(const string &filename, int timestep);
void Accumulate_Time_Averaged_Statistics(double dt);
void Reset_LES_Statistics();

//=============================================================================
// BOUNDARY CONDITIONS
//=============================================================================

void Apply_LES_Wall_BC(int cell_index);
void Apply_LES_Inlet_BC(int cell_index);
void Apply_LES_Outlet_BC(int cell_index);
void Apply_LES_Periodic_BC(int cell_index);
void Apply_All_LES_Boundary_Conditions();

// Inflow turbulence generation
void Generate_Synthetic_Turbulence_Inlet();
void Apply_Recycling_Rescaling_BC();
void Apply_Digital_Filter_Inflow();

//=============================================================================
// MAIN SOLVER INTERFACE
//=============================================================================

void Solve_LES_Step(double dt);
void Update_LES_Variables(double dt);
bool LES_Solver_Available();
void Setup_LES_From_Config(const Json::Value &config);
/** Stash LES block from driver JSON; apply after grid/cells exist. */
void Store_Pending_LES_Config(const Json::Value &root);
void Apply_Pending_LES_Setup();
/** Overwrite Primitive_Cells[*][8/9] with μ_eff / k_eff for viscous fluxes. */
void Inject_LES_Viscosity_Into_Primitives();

//=============================================================================
// GPU (CUDA) INTERFACE FUNCTIONS
//=============================================================================

#ifdef USE_CUDA
bool LES_GPU_Available();
void Initialize_LES_GPU();
void Finalize_LES_GPU();
bool Calculate_SGS_Viscosity_GPU();
bool Calculate_Smagorinsky_GPU();
bool Calculate_WALE_GPU();
bool Calculate_Vreman_GPU();
bool Calculate_Dynamic_Smagorinsky_GPU();
void Copy_LES_Data_To_Device();
void Copy_LES_Data_From_Device();
#endif

//=============================================================================
// VALIDATION AND VERIFICATION
//=============================================================================

void Validate_LES_Implementation();
void Test_Decaying_Isotropic_Turbulence();
void Test_Channel_Flow_LES();
void Test_Backward_Facing_Step_LES();
void Test_Circular_Cylinder_LES();

#endif // LES_MODELS_H
