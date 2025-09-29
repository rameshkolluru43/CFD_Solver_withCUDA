/**
 * @file K_Epsilon_Model.cpp
 * @brief Implementation of K-epsilon turbulence model for 2D compressible flows
 * @author Ramesh Kolluru
 * @date 2025-01-XX
 *
 * This file implements the standard K-epsilon two-equation turbulence model
 * for Reynolds-Averaged Navier-Stokes (RANS) simulations. The model solves
 * transport equations for turbulent kinetic energy (k) and turbulent
 * dissipation rate (epsilon).
 *
 * Model Equations:
 * ∂k/∂t + ∇·(ρuk) = ∇·[(μ + μt/σk)∇k] + Pk - ρε
 * ∂ε/∂t + ∇·(ρuε) = ∇·[(μ + μt/σε)∇ε] + C1(ε/k)Pk - C2ρ(ε²/k)
 *
 * where μt = ρCμk²/ε
 */

#include "Turbulence_Models.h"
#include "Solver.h"
#include "Boundary_Conditions.h"
#include "Utilities.h"
#include <algorithm>
#include <cmath>

//=============================================================================
// K-EPSILON MODEL INITIALIZATION
//=============================================================================

/**
 * @brief Initialize K-epsilon model variables for all cells
 */
void Initialize_KEpsilon_Variables()
{
    cout << "Initializing K-epsilon turbulence model..." << endl;

    // Resize turbulence variables vector
    turbulence_vars.resize(Total_Cells);

    // Set initial values based on inlet conditions or default values
    double k_init = 0.001 * pow(Inlet_Condition.u, 2); // 0.1% turbulence intensity
    double epsilon_init = pow(k_init, 1.5) / (0.1 * characteristic_length);

    for (int i = 0; i < Total_Cells; i++)
    {
        if (Cells[i].cellType != GHOST_CELL)
        {
            turbulence_vars[i].k = max(k_init, 1e-10);
            turbulence_vars[i].epsilon = max(epsilon_init, 1e-10);
            turbulence_vars[i].mut = 0.0;
            turbulence_vars[i].nut = 0.0;

            // Calculate initial turbulent viscosity
            Calculate_KEpsilon_Turbulent_Viscosity(i);
        }
    }

    cout << "K-epsilon model initialized with:" << endl;
    cout << "  Initial k = " << k_init << " m²/s²" << endl;
    cout << "  Initial ε = " << epsilon_init << " m²/s³" << endl;
}

//=============================================================================
// TURBULENT VISCOSITY CALCULATION
//=============================================================================

/**
 * @brief Calculate turbulent viscosity using K-epsilon model
 * @param cell_index Index of the current cell
 */
void Calculate_KEpsilon_Turbulent_Viscosity(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    double k = turbulence_vars[cell_index].k;
    double epsilon = turbulence_vars[cell_index].epsilon;
    double rho = Cells[cell_index].Conservative_Variables[0];

    // Ensure positive values
    k = max(k, 1e-10);
    epsilon = max(epsilon, 1e-10);

    // Calculate turbulent viscosity: μt = ρCμk²/ε
    turbulence_vars[cell_index].mut = rho * KEpsilon::C_mu * k * k / epsilon;
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
 * @brief Calculate production terms for K-epsilon model
 * @param cell_index Index of the current cell
 */
void Calculate_KEpsilon_Production_Terms(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    // Calculate strain rate tensor components
    double S[3][3];
    Calculate_Strain_Rate_Tensor(cell_index, S);

    // Calculate strain rate magnitude: S = sqrt(2SijSij)
    double S_magnitude = 0.0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            S_magnitude += S[i][j] * S[i][j];
        }
    }
    S_magnitude = sqrt(2.0 * S_magnitude);

    // Calculate production of turbulent kinetic energy: Pk = μt * S²
    double mut = turbulence_vars[cell_index].mut;
    turbulence_vars[cell_index].Pk = mut * S_magnitude * S_magnitude;

    // Limit production to prevent excessive values
    double k = turbulence_vars[cell_index].k;
    double epsilon = turbulence_vars[cell_index].epsilon;
    double rho = Cells[cell_index].Conservative_Variables[0];

    // Apply Kato-Launder correction if needed
    double max_production = 10.0 * rho * epsilon;
    turbulence_vars[cell_index].Pk = min(turbulence_vars[cell_index].Pk, max_production);
}

//=============================================================================
// TRANSPORT EQUATIONS SOLVER
//=============================================================================

/**
 * @brief Solve K-epsilon transport equations using explicit time integration
 * @param cell_index Index of the current cell
 * @param dt Time step size
 */
void Solve_KEpsilon_Transport_Equations(int cell_index, double dt)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return;

    double rho = Cells[cell_index].Conservative_Variables[0];
    double k = turbulence_vars[cell_index].k;
    double epsilon = turbulence_vars[cell_index].epsilon;
    double Pk = turbulence_vars[cell_index].Pk;

    // Calculate diffusion terms (simplified for structured grids)
    double dk_diff = Calculate_KEpsilon_K_Diffusion(cell_index);
    double depsilon_diff = Calculate_KEpsilon_Epsilon_Diffusion(cell_index);

    // K-equation source terms
    double Sk = Pk - rho * epsilon;

    // Epsilon-equation source terms
    double Sepsilon = KEpsilon::C_1 * (epsilon / k) * Pk - KEpsilon::C_2 * rho * (epsilon * epsilon / k);

    // Time integration (explicit Euler)
    double k_new = k + dt * (dk_diff + Sk) / rho;
    double epsilon_new = epsilon + dt * (depsilon_diff + Sepsilon) / rho;

    // Apply positivity constraints
    turbulence_vars[cell_index].k = max(k_new, 1e-10);
    turbulence_vars[cell_index].epsilon = max(epsilon_new, 1e-10);

    // Update turbulent viscosity
    Calculate_KEpsilon_Turbulent_Viscosity(cell_index);
}

//=============================================================================
// DIFFUSION TERMS CALCULATION
//=============================================================================

/**
 * @brief Calculate diffusion term for k-equation
 * @param cell_index Index of the current cell
 * @return Diffusion contribution to k-equation
 */
double Calculate_KEpsilon_K_Diffusion(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return 0.0;

    double diffusion = 0.0;
    double mu_laminar = Calculate_Laminar_Viscosity(cell_index);
    double mut = turbulence_vars[cell_index].mut;
    double mu_eff = mu_laminar + mut / KEpsilon::sigma_k;

    // Calculate diffusion using face-based approach
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
 * @brief Calculate diffusion term for epsilon-equation
 * @param cell_index Index of the current cell
 * @return Diffusion contribution to epsilon-equation
 */
double Calculate_KEpsilon_Epsilon_Diffusion(int cell_index)
{
    if (Cells[cell_index].cellType == GHOST_CELL)
        return 0.0;

    double diffusion = 0.0;
    double mu_laminar = Calculate_Laminar_Viscosity(cell_index);
    double mut = turbulence_vars[cell_index].mut;
    double mu_eff = mu_laminar + mut / KEpsilon::sigma_e;

    // Calculate diffusion using face-based approach
    for (int face = 0; face < Cells[cell_index].No_of_Faces; face++)
    {
        int neighbor = Cells[cell_index].Neighbours[face];
        if (neighbor >= 0 && neighbor < Total_Cells)
        {
            double epsilon_neighbor = turbulence_vars[neighbor].epsilon;
            double epsilon_cell = turbulence_vars[cell_index].epsilon;

            double face_area = Cells[cell_index].Face_Areas[face];
            double distance = Calculate_Cell_Distance(cell_index, neighbor);

            diffusion += mu_eff * (epsilon_neighbor - epsilon_cell) * face_area / distance;
        }
    }

    return diffusion;
}

//=============================================================================
// WALL FUNCTIONS IMPLEMENTATION
//=============================================================================

/**
 * @brief Apply wall functions for K-epsilon model
 * @param cell_index Index of the wall-adjacent cell
 */
void Apply_KEpsilon_Wall_Functions(int cell_index)
{
    if (!Cells[cell_index].hasBoundaryface)
        return;

    // Calculate wall distance and friction velocity
    Calculate_Wall_Distance(cell_index);
    Calculate_Wall_Shear_Stress(cell_index);
    Calculate_Y_Plus(cell_index);

    double y_plus = turbulence_vars[cell_index].y_plus;
    double u_tau = turbulence_vars[cell_index].u_tau;
    double rho = Cells[cell_index].Conservative_Variables[0];

    if (y_plus > 30.0)
    { // Use log-law region
        // Standard wall function approach
        double k_wall = u_tau * u_tau / sqrt(KEpsilon::C_mu);
        double epsilon_wall = pow(u_tau, 3) / (KEpsilon::kappa * turbulence_vars[cell_index].y_wall);

        turbulence_vars[cell_index].k = k_wall;
        turbulence_vars[cell_index].epsilon = epsilon_wall;
    }
    else if (y_plus < 1.0)
    { // Viscous sublayer
        // Near-wall treatment for low-Re
        double y_wall = turbulence_vars[cell_index].y_wall;
        double mu_laminar = Calculate_Laminar_Viscosity(cell_index);

        turbulence_vars[cell_index].k = 0.0;
        turbulence_vars[cell_index].epsilon = 2.0 * mu_laminar * turbulence_vars[cell_index].k /
                                              (rho * y_wall * y_wall);
    }
    // For buffer layer (1 < y+ < 30), use blending or enhanced wall treatment

    // Update turbulent viscosity
    Calculate_KEpsilon_Turbulent_Viscosity(cell_index);
}

//=============================================================================
// BOUNDARY CONDITIONS
//=============================================================================

/**
 * @brief Apply inlet boundary conditions for K-epsilon model
 * @param cell_index Index of the inlet cell
 */
void Apply_KEpsilon_Inlet_BC(int cell_index)
{
    // Set inlet turbulence values based on turbulence intensity and length scale
    double U_inlet = sqrt(Inlet_Condition.u * Inlet_Condition.u +
                          Inlet_Condition.v * Inlet_Condition.v);
    double Tu = 0.05;                       // 5% turbulence intensity (should be input parameter)
    double L = 0.1 * characteristic_length; // Turbulent length scale

    double k_inlet = 1.5 * pow(Tu * U_inlet, 2);
    double epsilon_inlet = KEpsilon::C_mu * pow(k_inlet, 1.5) / L;

    turbulence_vars[cell_index].k = k_inlet;
    turbulence_vars[cell_index].epsilon = epsilon_inlet;

    Calculate_KEpsilon_Turbulent_Viscosity(cell_index);
}

/**
 * @brief Apply outlet boundary conditions for K-epsilon model
 * @param cell_index Index of the outlet cell
 */
void Apply_KEpsilon_Outlet_BC(int cell_index)
{
    // Zero gradient boundary condition
    int upstream_neighbor = Find_Upstream_Neighbor(cell_index);
    if (upstream_neighbor >= 0)
    {
        turbulence_vars[cell_index].k = turbulence_vars[upstream_neighbor].k;
        turbulence_vars[cell_index].epsilon = turbulence_vars[upstream_neighbor].epsilon;
    }

    Calculate_KEpsilon_Turbulent_Viscosity(cell_index);
}

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Calculate strain rate magnitude for a given cell
 * @param cell_index Index of the current cell
 * @return Strain rate magnitude
 */
double Calculate_Strain_Rate_Magnitude(int cell_index)
{
    double S[3][3];
    Calculate_Strain_Rate_Tensor(cell_index, S);

    double S_magnitude = 0.0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            S_magnitude += S[i][j] * S[i][j];
        }
    }

    return sqrt(2.0 * S_magnitude);
}

/**
 * @brief Calculate laminar viscosity for a given cell
 * @param cell_index Index of the current cell
 * @return Laminar viscosity
 */
double Calculate_Laminar_Viscosity(int cell_index)
{
    double rho = Cells[cell_index].Conservative_Variables[0];
    double E = Cells[cell_index].Conservative_Variables[4];
    double u = Cells[cell_index].Conservative_Variables[1] / rho;
    double v = Cells[cell_index].Conservative_Variables[2] / rho;
    double w = Cells[cell_index].Conservative_Variables[3] / rho;

    double kinetic_energy = 0.5 * (u * u + v * v + w * w);
    double internal_energy = E / rho - kinetic_energy;
    double T = internal_energy / cv;

    // Sutherland's law
    return C_1 * pow(T, 1.5) / (T + T_S_Mu);
}

/**
 * @brief Calculate distance between two cell centers
 * @param cell1 Index of first cell
 * @param cell2 Index of second cell
 * @return Distance between cell centers
 */
double Calculate_Cell_Distance(int cell1, int cell2)
{
    double dx = Cells[cell1].Cell_Center[0] - Cells[cell2].Cell_Center[0];
    double dy = Cells[cell1].Cell_Center[1] - Cells[cell2].Cell_Center[1];
    double dz = Cells[cell1].Cell_Center[2] - Cells[cell2].Cell_Center[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Find upstream neighbor for outlet boundary condition
 * @param cell_index Index of the outlet cell
 * @return Index of upstream neighbor cell
 */
int Find_Upstream_Neighbor(int cell_index)
{
    // Simple implementation - find neighbor with lowest x-coordinate
    // Should be improved based on actual flow direction
    int upstream = -1;
    double min_x = 1e10;

    for (int i = 0; i < Cells[cell_index].No_of_Faces; i++)
    {
        int neighbor = Cells[cell_index].Neighbours[i];
        if (neighbor >= 0 && neighbor < Total_Cells)
        {
            if (Cells[neighbor].Cell_Center[0] < min_x)
            {
                min_x = Cells[neighbor].Cell_Center[0];
                upstream = neighbor;
            }
        }
    }

    return upstream;
}