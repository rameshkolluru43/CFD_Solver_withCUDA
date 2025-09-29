/**
 * @file K_Omega_Model.cpp
 * @brief Implementation of K-omega turbulence model for 2D compressible flows
 * @author Ramesh Kolluru
 * @date 2025-01-XX
 *
 * This file implements the K-omega two-equation turbulence model including:
 * - Standard Wilcox K-omega model
 * - SST (Shear Stress Transport) K-omega model
 *
 * The model solves transport equations for turbulent kinetic energy (k)
 * and specific dissipation rate (omega).
 *
 * Model Equations:
 * ∂k/∂t + ∇·(ρuk) = ∇·[(μ + σ*μt)∇k] + Pk - β*ρkω
 * ∂ω/∂t + ∇·(ρuω) = ∇·[(μ + σ*μt)∇ω] + α(ω/k)Pk - βρω²
 *
 * where μt = ρk/ω (for standard model)
 */

#include "Turbulence_Models.h"
#include "Solver.h"
#include "Boundary_Conditions.h"
#include "Utilities.h"
#include <algorithm>
#include <cmath>

//=============================================================================
// K-OMEGA MODEL INITIALIZATION
//=============================================================================

/**
 * @brief Initialize K-omega model variables for all cells
 */
void Initialize_KOmega_Variables()
{
    cout << "Initializing K-omega turbulence model..." << endl;

    // Resize turbulence variables vector
    turbulence_vars.resize(Total_Cells);

    // Set initial values based on inlet conditions
    double k_init = 0.001 * pow(Inlet_Condition.u, 2); // 0.1% turbulence intensity
    double L_turb = 0.1 * characteristic_length;       // Turbulent length scale
    double omega_init = k_init / (sqrt(KOmega::beta_star) * L_turb);

    for (int i = 0; i < Total_Cells; i++)
    {
        if (Cells[i].cellType != GHOST_CELL)
        {
            turbulence_vars[i].k = max(k_init, 1e-10);
            turbulence_vars[i].omega = max(omega_init, 1e-6);
            turbulence_vars[i].mut = 0.0;
            turbulence_vars[i].nut = 0.0;

            // Calculate initial turbulent viscosity
            Calculate_KOmega_Turbulent_Viscosity(i);
        }
    }

    cout << "K-omega model initialized with:" << endl;
    cout << "  Initial k = " << k_init << " m²/s²" << endl;
    cout << "  Initial ω = " << omega_init << " 1/s" << endl;
}

//=============================================================================
// TURBULENT VISCOSITY CALCULATION
//=============================================================================

/**
 * @brief Calculate turbulent viscosity using K-omega model
 * @param cell_index Index of the current cell
 */
void Calculate_KOmega_Turbulent_Viscosity(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    double k = turbulence_vars[cell_index].k;
    double omega = turbulence_vars[cell_index].omega;
    double rho = Cells[cell_index].Conservative_Variables[0];

    // Ensure positive values
    k = max(k, 1e-10);
    omega = max(omega, 1e-6);

    if (current_turbulence_model == TurbulenceModel::K_OMEGA_SST)
    {
        // SST model: μt = ρa₁k/max(a₁ω, ΩF₂)
        double vorticity_magnitude = Calculate_Vorticity_Magnitude(cell_index);
        double F2 = Calculate_SST_F2_Function(cell_index);

        double denominator = max(SST::a1 * omega, vorticity_magnitude * F2);
        turbulence_vars[cell_index].mut = rho * SST::a1 * k / denominator;
    }
    else
    {
        // Standard K-omega: μt = ρk/ω
        turbulence_vars[cell_index].mut = rho * k / omega;
    }

    turbulence_vars[cell_index].nut = turbulence_vars[cell_index].mut / rho;

    // Apply realizability constraint
    double S_magnitude = Calculate_Strain_Rate_Magnitude(cell_index);
    if (S_magnitude > 0)
    {
        double realizability_limit = 2.0 * k / S_magnitude;
        turbulence_vars[cell_index].nut = min(turbulence_vars[cell_index].nut, realizability_limit);
        turbulence_vars[cell_index].mut = turbulence_vars[cell_index].nut * rho;
    }

    // Limit maximum turbulent viscosity
    double mu_laminar = Calculate_Laminar_Viscosity(cell_index);
    turbulence_vars[cell_index].mut = min(turbulence_vars[cell_index].mut, 1000.0 * mu_laminar);
}

//=============================================================================
// PRODUCTION TERMS CALCULATION
//=============================================================================

/**
 * @brief Calculate production terms for K-omega model
 * @param cell_index Index of the current cell
 */
void Calculate_KOmega_Production_Terms(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    // Calculate strain rate tensor
    double S[3][3];
    Calculate_Strain_Rate_Tensor(cell_index, S);

    // Calculate strain rate magnitude
    double S_magnitude = 0.0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            S_magnitude += S[i][j] * S[i][j];
        }
    }
    S_magnitude = sqrt(2.0 * S_magnitude);

    // Calculate production of turbulent kinetic energy
    double mut = turbulence_vars[cell_index].mut;
    turbulence_vars[cell_index].Pk = mut * S_magnitude * S_magnitude;

    // Limit production to prevent excessive values
    double k = turbulence_vars[cell_index].k;
    double omega = turbulence_vars[cell_index].omega;
    double rho = Cells[cell_index].Conservative_Variables[0];

    if (current_turbulence_model == TurbulenceModel::K_OMEGA_SST)
    {
        // SST production limiter
        double max_production = 10.0 * KOmega::beta_star * rho * k * omega;
        turbulence_vars[cell_index].Pk = min(turbulence_vars[cell_index].Pk, max_production);
    }
    else
    {
        // Standard K-omega production limiter
        double max_production = 10.0 * KOmega::beta * rho * k * omega;
        turbulence_vars[cell_index].Pk = min(turbulence_vars[cell_index].Pk, max_production);
    }
}

//=============================================================================
// TRANSPORT EQUATIONS SOLVER
//=============================================================================

/**
 * @brief Solve K-omega transport equations
 * @param cell_index Index of the current cell
 * @param dt Time step size
 */
void Solve_KOmega_Transport_Equations(int cell_index, double dt)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    double rho = Cells[cell_index].Conservative_Variables[0];
    double k = turbulence_vars[cell_index].k;
    double omega = turbulence_vars[cell_index].omega;
    double Pk = turbulence_vars[cell_index].Pk;

    // Calculate diffusion terms
    double dk_diff = Calculate_KOmega_K_Diffusion(cell_index);
    double domega_diff = Calculate_KOmega_Omega_Diffusion(cell_index);

    if (current_turbulence_model == TurbulenceModel::K_OMEGA_SST)
    {
        Solve_SST_Transport_Equations(cell_index, dt);
        return;
    }

    // Standard K-omega model
    // K-equation source terms
    double Sk = Pk - KOmega::beta_star * rho * k * omega;

    // Omega-equation source terms
    double Somega = KOmega::alpha * (omega / k) * Pk - KOmega::beta * rho * omega * omega;

    // Time integration (explicit Euler)
    double k_new = k + dt * (dk_diff + Sk) / rho;
    double omega_new = omega + dt * (domega_diff + Somega) / rho;

    // Apply positivity constraints
    turbulence_vars[cell_index].k = max(k_new, 1e-10);
    turbulence_vars[cell_index].omega = max(omega_new, 1e-6);

    // Update turbulent viscosity
    Calculate_KOmega_Turbulent_Viscosity(cell_index);
}

//=============================================================================
// SST K-OMEGA MODEL IMPLEMENTATION
//=============================================================================

/**
 * @brief Solve SST K-omega transport equations
 * @param cell_index Index of the current cell
 * @param dt Time step size
 */
void Solve_SST_Transport_Equations(int cell_index, double dt)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    double rho = Cells[cell_index].Conservative_Variables[0];
    double k = turbulence_vars[cell_index].k;
    double omega = turbulence_vars[cell_index].omega;
    double Pk = turbulence_vars[cell_index].Pk;

    // Calculate blending functions
    double F1, F2;
    Calculate_SST_Blending_Functions(cell_index, F1, F2);

    // Blended constants
    double sigma_k = F1 * SST::sigma_k1 + (1.0 - F1) * SST::sigma_k2;
    double sigma_omega = F1 * SST::sigma_w1 + (1.0 - F1) * SST::sigma_w2;
    double alpha = F1 * SST::alpha_1 + (1.0 - F1) * SST::alpha_2;
    double beta = F1 * SST::beta_1 + (1.0 - F1) * SST::beta_2;

    // Calculate diffusion terms with blended constants
    double dk_diff = Calculate_SST_K_Diffusion(cell_index, sigma_k);
    double domega_diff = Calculate_SST_Omega_Diffusion(cell_index, sigma_omega);

    // Cross-diffusion term for omega equation
    double CDkomega = Calculate_SST_Cross_Diffusion(cell_index);

    // Source terms
    double Sk = Pk - SST::beta_star * rho * k * omega;
    double Somega = alpha * (omega / k) * Pk - beta * rho * omega * omega +
                    2.0 * (1.0 - F1) * SST::sigma_w2 * rho * CDkomega;

    // Time integration
    double k_new = k + dt * (dk_diff + Sk) / rho;
    double omega_new = omega + dt * (domega_diff + Somega) / rho;

    // Apply constraints
    turbulence_vars[cell_index].k = max(k_new, 1e-10);
    turbulence_vars[cell_index].omega = max(omega_new, 1e-6);

    // Update turbulent viscosity
    Calculate_KOmega_Turbulent_Viscosity(cell_index);
}

/**
 * @brief Calculate SST blending functions F1 and F2
 * @param cell_index Index of the current cell
 * @param F1 Reference to F1 blending function
 * @param F2 Reference to F2 blending function
 */
void Calculate_SST_Blending_Functions(int cell_index, double &F1, double &F2)
{
    double k = turbulence_vars[cell_index].k;
    double omega = turbulence_vars[cell_index].omega;
    double rho = Cells[cell_index].Conservative_Variables[0];
    double mu = Calculate_Laminar_Viscosity(cell_index);

    double y_wall = turbulence_vars[cell_index].y_wall;

    // Calculate CDkomega
    double CDkomega = Calculate_SST_Cross_Diffusion(cell_index);
    CDkomega = max(CDkomega, 1e-20);

    // Calculate arguments for F1
    double arg1_1 = sqrt(k) / (SST::beta_star * omega * y_wall);
    double arg1_2 = 500.0 * mu / (rho * omega * y_wall * y_wall);
    double arg1_3 = 4.0 * SST::sigma_w2 * rho * k / (CDkomega * y_wall * y_wall);

    double arg1 = min(max(arg1_1, arg1_2), arg1_3);
    F1 = tanh(pow(arg1, 4));

    // Calculate argument for F2
    double arg2_1 = 2.0 * sqrt(k) / (SST::beta_star * omega * y_wall);
    double arg2_2 = 500.0 * mu / (rho * omega * y_wall * y_wall);

    double arg2 = max(arg2_1, arg2_2);
    F2 = tanh(pow(arg2, 2));
}

/**
 * @brief Calculate cross-diffusion term for SST model
 * @param cell_index Index of the current cell
 * @return Cross-diffusion term CDkomega
 */
double Calculate_SST_Cross_Diffusion(int cell_index)
{
    double CDkomega = 0.0;

    // Calculate gradients of k and omega
    double grad_k[3], grad_omega[3];
    Calculate_Scalar_Gradient(cell_index, "k", grad_k);
    Calculate_Scalar_Gradient(cell_index, "omega", grad_omega);

    // CDkomega = (∇k · ∇ω) / ω
    double omega = turbulence_vars[cell_index].omega;
    for (int i = 0; i < 3; i++)
    {
        CDkomega += grad_k[i] * grad_omega[i];
    }
    CDkomega = CDkomega / omega;

    return max(CDkomega, 0.0);
}

//=============================================================================
// DIFFUSION TERMS CALCULATION
//=============================================================================

/**
 * @brief Calculate diffusion term for k-equation in K-omega model
 * @param cell_index Index of the current cell
 * @return Diffusion contribution to k-equation
 */
double Calculate_KOmega_K_Diffusion(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return 0.0;

    double diffusion = 0.0;
    double mu_laminar = Calculate_Laminar_Viscosity(cell_index);
    double mut = turbulence_vars[cell_index].mut;
    double mu_eff = mu_laminar + KOmega::sigma_star * mut;

    // Face-based diffusion calculation
    for (int face = 0; face < Cells[cell_index].No_of_Faces; face++)
    {
        int neighbor = Cells[cell_index].Neighbours[face];
        if (neighbor >= 0 && neighbor < Total_Cells)
        {
            double k_neighbor = turbulence_vars[neighbor].k;
            double k_cell = turbulence_vars[cell_index].k;

            double face_area = Cells[cell_index].Face_Areas[face];
            double distance = Calculate_Cell_Distance(cell_index, neighbor);

            diffusion += mu_eff * (k_neighbor - k_cell) * face_area / distance;
        }
    }

    return diffusion;
}

/**
 * @brief Calculate diffusion term for omega-equation in K-omega model
 * @param cell_index Index of the current cell
 * @return Diffusion contribution to omega-equation
 */
double Calculate_KOmega_Omega_Diffusion(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return 0.0;

    double diffusion = 0.0;
    double mu_laminar = Calculate_Laminar_Viscosity(cell_index);
    double mut = turbulence_vars[cell_index].mut;
    double mu_eff = mu_laminar + KOmega::sigma * mut;

    // Face-based diffusion calculation
    for (int face = 0; face < Cells[cell_index].No_of_Faces; face++)
    {
        int neighbor = Cells[cell_index].Neighbours[face];
        if (neighbor >= 0 && neighbor < Total_Cells)
        {
            double omega_neighbor = turbulence_vars[neighbor].omega;
            double omega_cell = turbulence_vars[cell_index].omega;

            double face_area = Cells[cell_index].Face_Areas[face];
            double distance = Calculate_Cell_Distance(cell_index, neighbor);

            diffusion += mu_eff * (omega_neighbor - omega_cell) * face_area / distance;
        }
    }

    return diffusion;
}

//=============================================================================
// WALL FUNCTIONS IMPLEMENTATION
//=============================================================================

/**
 * @brief Apply wall functions for K-omega model
 * @param cell_index Index of the wall-adjacent cell
 */
void Apply_KOmega_Wall_Functions(int cell_index)
{
    if (!Cells[cell_index].hasBoundaryface)
        return;

    // Calculate wall parameters
    Calculate_Wall_Distance(cell_index);
    Calculate_Wall_Shear_Stress(cell_index);
    Calculate_Y_Plus(cell_index);

    double y_plus = turbulence_vars[cell_index].y_plus;
    double u_tau = turbulence_vars[cell_index].u_tau;
    double y_wall = turbulence_vars[cell_index].y_wall;
    double rho = Cells[cell_index].Conservative_Variables[0];
    double mu = Calculate_Laminar_Viscosity(cell_index);

    if (y_plus > 30.0)
    { // Log-law region
        // K at the wall
        double k_wall = u_tau * u_tau / sqrt(KOmega::beta_star);

        // Omega at the wall using log-law
        double omega_wall = u_tau / (sqrt(KOmega::beta_star) * KOmega::kappa * y_wall);

        turbulence_vars[cell_index].k = k_wall;
        turbulence_vars[cell_index].omega = omega_wall;
    }
    else if (y_plus < 1.0)
    { // Viscous sublayer
        // Near-wall treatment
        turbulence_vars[cell_index].k = 0.0;

        // Omega in viscous sublayer
        double omega_viscous = 60.0 * mu / (rho * KOmega::beta * y_wall * y_wall);
        turbulence_vars[cell_index].omega = omega_viscous;
    }

    // Update turbulent viscosity
    Calculate_KOmega_Turbulent_Viscosity(cell_index);
}

//=============================================================================
// BOUNDARY CONDITIONS
//=============================================================================

/**
 * @brief Apply inlet boundary conditions for K-omega model
 * @param cell_index Index of the inlet cell
 */
void Apply_KOmega_Inlet_BC(int cell_index)
{
    double U_inlet = sqrt(Inlet_Condition.u * Inlet_Condition.u +
                          Inlet_Condition.v * Inlet_Condition.v);
    double Tu = 0.05;                       // 5% turbulence intensity
    double L = 0.1 * characteristic_length; // Turbulent length scale

    double k_inlet = 1.5 * pow(Tu * U_inlet, 2);
    double omega_inlet = k_inlet / (sqrt(KOmega::beta_star) * L);

    turbulence_vars[cell_index].k = k_inlet;
    turbulence_vars[cell_index].omega = omega_inlet;

    Calculate_KOmega_Turbulent_Viscosity(cell_index);
}

/**
 * @brief Apply outlet boundary conditions for K-omega model
 * @param cell_index Index of the outlet cell
 */
void Apply_KOmega_Outlet_BC(int cell_index)
{
    // Zero gradient boundary condition
    int upstream_neighbor = Find_Upstream_Neighbor(cell_index);
    if (upstream_neighbor >= 0)
    {
        turbulence_vars[cell_index].k = turbulence_vars[upstream_neighbor].k;
        turbulence_vars[cell_index].omega = turbulence_vars[upstream_neighbor].omega;
    }

    Calculate_KOmega_Turbulent_Viscosity(cell_index);
}

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Calculate vorticity magnitude for SST model
 * @param cell_index Index of the current cell
 * @return Vorticity magnitude
 */
double Calculate_Vorticity_Magnitude(int cell_index)
{
    double Omega[3][3];
    Calculate_Vorticity_Tensor(cell_index, Omega);

    double vorticity_magnitude = 0.0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vorticity_magnitude += Omega[i][j] * Omega[i][j];
        }
    }

    return sqrt(2.0 * vorticity_magnitude);
}

/**
 * @brief Calculate F2 blending function for SST turbulent viscosity
 * @param cell_index Index of the current cell
 * @return F2 blending function value
 */
double Calculate_SST_F2_Function(int cell_index)
{
    double F1, F2;
    Calculate_SST_Blending_Functions(cell_index, F1, F2);
    return F2;
}

/**
 * @brief Calculate scalar gradient for a given variable
 * @param cell_index Index of the current cell
 * @param variable Variable name ("k" or "omega")
 * @param grad Array to store gradient components
 */
void Calculate_Scalar_Gradient(int cell_index, const string &variable, double grad[3])
{
    // Initialize gradient
    grad[0] = grad[1] = grad[2] = 0.0;

    double value_center;
    if (variable == "k")
    {
        value_center = turbulence_vars[cell_index].k;
    }
    else if (variable == "omega")
    {
        value_center = turbulence_vars[cell_index].omega;
    }
    else
    {
        return;
    }

    // Green-Gauss gradient calculation
    for (int face = 0; face < Cells[cell_index].No_of_Faces; face++)
    {
        int neighbor = Cells[cell_index].Neighbours[face];
        double value_neighbor;

        if (neighbor >= 0 && neighbor < Total_Cells)
        {
            if (variable == "k")
            {
                value_neighbor = turbulence_vars[neighbor].k;
            }
            else
            {
                value_neighbor = turbulence_vars[neighbor].omega;
            }
        }
        else
        {
            value_neighbor = value_center; // Zero gradient at boundary
        }

        double face_value = 0.5 * (value_center + value_neighbor);
        double face_area = Cells[cell_index].Face_Areas[face];

        // Face normal components (assuming they are stored)
        double nx = Cells[cell_index].Face_Normals[face][0];
        double ny = Cells[cell_index].Face_Normals[face][1];
        double nz = Cells[cell_index].Face_Normals[face][2];

        grad[0] += face_value * nx * face_area;
        grad[1] += face_value * ny * face_area;
        grad[2] += face_value * nz * face_area;
    }

    // Divide by cell volume
    double cell_volume = Cells[cell_index].Volume;
    grad[0] /= cell_volume;
    grad[1] /= cell_volume;
    grad[2] /= cell_volume;
}