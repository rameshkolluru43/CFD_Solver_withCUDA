/**
 * @file Turbulence_Models.h
 * @brief Header file for RANS turbulence models implementation
 * @author Ramesh Kolluru
 * @date 2025-01-XX
 *
 * This file contains the declarations for K-epsilon and K-omega turbulence models
 * for 2D compressible flow simulations. The implementation includes:
 * - Standard K-epsilon model
 * - Standard K-omega model
 * - Wilcox K-omega model
 * - SST K-omega model
 */

#ifndef TURBULENCE_MODELS_H
#define TURBULENCE_MODELS_H

#include "definitions.h"
#include "Globals.h"

//=============================================================================
// TURBULENCE MODEL CONSTANTS
//=============================================================================

// K-epsilon model constants (Standard model)
namespace KEpsilon
{
    const double C_mu = 0.09;
    const double C1_eps = 1.44;
    const double C_2 = 1.92;
    const double sigma_k = 1.0;
    const double sigma_e = 1.3;
    const double kappa = 0.41; // Von Karman constant
    const double E = 9.793;    // Wall function constant
}

// K-omega model constants (Standard Wilcox model)
namespace KOmega
{
    const double beta_star = 0.09;
    const double beta = 0.075;
    const double sigma = 0.5;
    const double sigma_star = 0.5;
    const double alpha = 5.0 / 9.0;
    const double kappa = 0.41; // Von Karman constant
    const double E = 9.793;    // Wall function constant
}

// SST K-omega model constants
namespace SST
{
    // Set 1 constants (K-omega)
    const double alpha_1 = 5.0 / 9.0;
    const double beta_1 = 0.075;
    const double beta_star = 0.09;
    const double sigma_k1 = 0.85;
    const double sigma_w1 = 0.5;

    // Set 2 constants (K-epsilon)
    const double alpha_2 = 0.44;
    const double beta_2 = 0.0828;
    const double sigma_k2 = 1.0;
    const double sigma_w2 = 0.856;

    // Common constants
    const double a1 = 0.31;
    const double kappa = 0.41;
    const double E = 9.793;
}

//=============================================================================
// TURBULENCE MODEL ENUMERATIONS
//=============================================================================

enum class TurbulenceModel
{
    LAMINAR = 0,
    K_EPSILON = 1,
    K_OMEGA_WILCOX = 2,
    K_OMEGA_SST = 3
};

//=============================================================================
// TURBULENCE VARIABLES STRUCTURE
//=============================================================================

struct TurbulenceVariables
{
    double k;        // Turbulent kinetic energy
    double epsilon;  // Turbulent dissipation rate
    double omega;    // Specific dissipation rate
    double mut;      // Turbulent viscosity
    double nut;      // Turbulent kinematic viscosity
    double y_plus;   // Dimensionless wall distance
    double y_wall;   // Wall distance
    double tau_wall; // Wall shear stress
    double u_tau;    // Friction velocity

    // Production and destruction terms
    double Pk;       // Production of k
    double Pepsilon; // Production of epsilon
    double Pomega;   // Production of omega

    // Residual storage for implicit/RK4 solvers
    double k_residual;
    double epsilon_residual;
    double omega_residual;

    // Effective transport properties
    double mu_effective;
    double k_effective;

    // Constructor
    TurbulenceVariables() : k(0.0), epsilon(0.0), omega(0.0), mut(0.0),
                            nut(0.0), y_plus(0.0), y_wall(0.0), tau_wall(0.0), u_tau(0.0),
                            Pk(0.0), Pepsilon(0.0), Pomega(0.0),
                            k_residual(0.0), epsilon_residual(0.0), omega_residual(0.0),
                            mu_effective(0.0), k_effective(0.0) {}
};

//=============================================================================
// GLOBAL TURBULENCE VARIABLES
//=============================================================================

extern TurbulenceModel current_turbulence_model;
extern vector<TurbulenceVariables> turbulence_vars;
extern bool use_wall_functions;
extern double y_plus_target;

//=============================================================================
// FUNCTION DECLARATIONS
//=============================================================================

// Initialization functions
void Initialize_Turbulence_Model(TurbulenceModel model);
void Initialize_Turbulence_Variables();
void Set_Initial_Turbulence_Conditions();

// K-epsilon model functions
void Calculate_KEpsilon_Production_Terms(int cell_index);
void Solve_KEpsilon_Transport_Equations(int cell_index, double dt);
void Calculate_KEpsilon_Turbulent_Viscosity(int cell_index);
void Apply_KEpsilon_Wall_Functions(int cell_index);

// K-omega model functions
void Calculate_KOmega_Production_Terms(int cell_index);
void Solve_KOmega_Transport_Equations(int cell_index, double dt);
void Calculate_KOmega_Turbulent_Viscosity(int cell_index);
void Apply_KOmega_Wall_Functions(int cell_index);

// SST K-omega model functions
void Calculate_SST_Production_Terms(int cell_index);
void Solve_SST_Transport_Equations(int cell_index, double dt);
void Calculate_SST_Turbulent_Viscosity(int cell_index);
void Calculate_SST_Blending_Functions(int cell_index, double &F1, double &F2);
void Apply_SST_Wall_Functions(int cell_index);

// Common utility functions
void Calculate_Strain_Rate_Tensor(int cell_index, double S[3][3]);
void Calculate_Vorticity_Tensor(int cell_index, double Omega[3][3]);
void Calculate_Turbulent_Production(int cell_index, const double S[3][3]);
void Calculate_Wall_Distance(int cell_index);
void Calculate_Wall_Distances_All_Cells();
void Calculate_Wall_Shear_Stress(int cell_index);
void Calculate_Y_Plus(int cell_index);
double Calculate_Laminar_Viscosity(int cell_index);
double Calculate_Strain_Rate_Magnitude(int cell_index);
double Calculate_Cell_Distance(int cell1, int cell2);
int Find_Upstream_Neighbor(int cell_index);
double Calculate_Vorticity_Magnitude(int cell_index);
double Calculate_SST_F2_Function(int cell_index);

// Diffusion term functions
double Calculate_KEpsilon_K_Diffusion(int cell_index);
double Calculate_KEpsilon_Epsilon_Diffusion(int cell_index);
double Calculate_KOmega_K_Diffusion(int cell_index);
double Calculate_KOmega_Omega_Diffusion(int cell_index);
double Calculate_SST_K_Diffusion(int cell_index, double sigma_k);
double Calculate_SST_Omega_Diffusion(int cell_index, double sigma_omega);
double Calculate_SST_Cross_Diffusion(int cell_index);
void Calculate_Scalar_Gradient(int cell_index, const string &variable, double grad[3]);

// Initialization functions (model-specific)
void Initialize_KEpsilon_Variables();
void Initialize_KOmega_Variables();

// Model-specific boundary conditions
void Apply_KEpsilon_Inlet_BC(int cell_index);
void Apply_KEpsilon_Outlet_BC(int cell_index);
void Apply_KOmega_Inlet_BC(int cell_index);
void Apply_KOmega_Outlet_BC(int cell_index);

// Helper functions
void Apply_All_Turbulence_Boundary_Conditions();
bool Is_Inlet_Cell(int cell_index);
bool Is_Outlet_Cell(int cell_index);
bool Is_Wall_Cell(int cell_index);
bool Is_Symmetry_Cell(int cell_index);
int Find_Interior_Neighbor(int cell_index);
void Apply_Low_Re_Wall_Treatment(int cell_index);

// Integration functions
bool Viscous_Solver_With_Turbulence(string &Error_Filename, string &Sol_Filename);
void Setup_Turbulence_From_Config(const string &config_file);
double Calculate_Turbulence_Error(const string &variable);
void Apply_Turbulence_Limiters_All_Cells();

// Boundary condition functions
void Apply_Turbulence_Inlet_BC(int cell_index);
void Apply_Turbulence_Outlet_BC(int cell_index);
void Apply_Turbulence_Wall_BC(int cell_index);
void Apply_Turbulence_Symmetry_BC(int cell_index);

// Source term functions
void Add_Turbulence_Source_Terms(int cell_index, vector<double> &source_k,
                                 vector<double> &source_epsilon,
                                 vector<double> &source_omega);

// Limiters and corrections
void Apply_Turbulence_Limiters(int cell_index);
void Apply_Realizability_Constraints(int cell_index);
void Calculate_Turbulent_Time_Scale(int cell_index);

// Update functions
void Update_Turbulence_Variables(double dt);
void Update_Effective_Viscosity();
void Update_Turbulent_Thermal_Conductivity();

// Output and diagnostics
void Write_Turbulence_Variables(const string &filename);
void Calculate_Turbulence_Statistics();
bool Check_Turbulence_Convergence();

// Validation and verification
void Validate_Turbulence_Implementation();
void Test_Flat_Plate_Boundary_Layer();
void Test_Channel_Flow();
void Test_Backward_Facing_Step();

#endif // TURBULENCE_MODELS_H