#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Utilities.h"

/**
 * @file Ausm_Flux.cpp
 * @brief 3D AUSM (Advection Upstream Splitting Method) flux implementation
 *
 * This file implements the AUSM flux scheme for 3D Euler/Navier-Stokes equations.
 * AUSM splits the flux into convective and pressure parts, providing excellent
 * performance for low Mach number flows and good shock-capturing capabilities.
 */

/**
 * @brief AUSM flux calculation for 3D hexahedral cell
 * @param Cell_No Current cell index
 *
 * Computes AUSM flux for all 6 faces of a 3D hexahedral cell.
 * The method splits convective and pressure fluxes using Mach number
 * splitting functions for robust performance across all speed regimes.
 */
void Ausm_Flux_3D(const int &Cell_No)
{
    // Initialize net flux for this cell to zero
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Get neighbor cell indices for all six faces
    vector<int> neighbors(NUM_FACES_3D);
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        neighbors[face] = Cells[Cell_No].Neighbours[face];
    }

    // Process each face of the hexahedral cell
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor_idx = neighbors[face];

        // Skip boundary faces (handled separately)
        if (neighbor_idx < 0 || neighbor_idx >= No_Physical_Cells)
        {
            continue;
        }

        // Calculate AUSM flux for this face
        V_D face_flux(NUM_CONSERVATIVE_VARS, 0.0);
        Compute_AUSM_Face_Flux_3D(Cell_No, neighbor_idx, face, face_flux);

        // Add contribution to net flux
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            Cells_Net_Flux[Cell_No][i] += face_flux[i];
        }
    }

    // Normalize by cell volume
    double inv_volume = Cells[Cell_No].Inv_Volume;
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        Cells_Net_Flux[Cell_No][i] *= inv_volume;
    }
}

/**
 * @brief Compute AUSM flux across a single face in 3D
 * @param left_cell_idx Left cell index
 * @param right_cell_idx Right cell index
 * @param face_idx Face index
 * @param face_flux Output flux vector
 */
void Compute_AUSM_Face_Flux_3D(const int &left_cell_idx, const int &right_cell_idx,
                               const int &face_idx, V_D &face_flux)
{
    // Get face geometry
    double nx = Cells[left_cell_idx].Face_Normals[3 * face_idx];
    double ny = Cells[left_cell_idx].Face_Normals[3 * face_idx + 1];
    double nz = Cells[left_cell_idx].Face_Normals[3 * face_idx + 2];
    double face_area = Cells[left_cell_idx].Face_Areas[face_idx];

    // Left state (current cell)
    double rho_L = U_Cells[left_cell_idx][RHO_INDEX];
    double rhou_L = U_Cells[left_cell_idx][RHU_INDEX];
    double rhov_L = U_Cells[left_cell_idx][RHV_INDEX];
    double rhow_L = U_Cells[left_cell_idx][RHW_INDEX];
    double rhoE_L = U_Cells[left_cell_idx][ENERGY_INDEX];

    double u_L = rhou_L / rho_L;
    double v_L = rhov_L / rho_L;
    double w_L = rhow_L / rho_L;
    double p_L = (GAMMA - 1.0) * (rhoE_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L + w_L * w_L));
    double a_L = sqrt(GAMMA * p_L / rho_L);
    double H_L = (rhoE_L + p_L) / rho_L;

    // Right state (neighbor cell)
    double rho_R = U_Cells[right_cell_idx][RHO_INDEX];
    double rhou_R = U_Cells[right_cell_idx][RHU_INDEX];
    double rhov_R = U_Cells[right_cell_idx][RHV_INDEX];
    double rhow_R = U_Cells[right_cell_idx][RHW_INDEX];
    double rhoE_R = U_Cells[right_cell_idx][ENERGY_INDEX];

    double u_R = rhou_R / rho_R;
    double v_R = rhov_R / rho_R;
    double w_R = rhow_R / rho_R;
    double p_R = (GAMMA - 1.0) * (rhoE_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R + w_R * w_R));
    double a_R = sqrt(GAMMA * p_R / rho_R);
    double H_R = (rhoE_R + p_R) / rho_R;

    // Face-normal velocities
    double Vn_L = u_L * nx + v_L * ny + w_L * nz;
    double Vn_R = u_R * nx + v_R * ny + w_R * nz;

    // Interface sound speed (Roe average)
    double sqrt_rho_L = sqrt(rho_L);
    double sqrt_rho_R = sqrt(rho_R);
    double sum_sqrt_rho = sqrt_rho_L + sqrt_rho_R;

    double a_interface = (sqrt_rho_L * a_L + sqrt_rho_R * a_R) / sum_sqrt_rho;

    // Mach numbers
    double M_L = Vn_L / a_interface;
    double M_R = Vn_R / a_interface;

    // AUSM+ splitting functions
    double M_plus, M_minus, P_plus, P_minus;

    // Mach number splitting
    if (fabs(M_L) <= 1.0)
    {
        M_plus = 0.25 * (M_L + 1.0) * (M_L + 1.0);
        P_plus = 0.25 * (M_L + 1.0) * (M_L + 1.0) * (2.0 - M_L);
    }
    else
    {
        M_plus = 0.5 * (M_L + fabs(M_L));
        P_plus = 0.5 * (1.0 + (M_L > 0.0 ? 1.0 : 0.0));
    }

    if (fabs(M_R) <= 1.0)
    {
        M_minus = -0.25 * (M_R - 1.0) * (M_R - 1.0);
        P_minus = 0.25 * (M_R - 1.0) * (M_R - 1.0) * (2.0 + M_R);
    }
    else
    {
        M_minus = 0.5 * (M_R - fabs(M_R));
        P_minus = 0.5 * (1.0 - (M_R > 0.0 ? 1.0 : 0.0));
    }

    // Interface Mach number and pressure
    double M_interface = M_plus + M_minus;
    double P_interface = P_plus * p_L + P_minus * p_R;

    // Mass flux
    double mass_flux;
    if (M_interface > 0.0)
    {
        mass_flux = a_interface * M_interface * rho_L;
    }
    else
    {
        mass_flux = a_interface * M_interface * rho_R;
    }

    // Convective flux
    if (mass_flux > 0.0)
    {
        // Upwind from left state
        face_flux[RHO_INDEX] = mass_flux;
        face_flux[RHU_INDEX] = mass_flux * u_L + P_interface * nx;
        face_flux[RHV_INDEX] = mass_flux * v_L + P_interface * ny;
        face_flux[RHW_INDEX] = mass_flux * w_L + P_interface * nz;
        face_flux[ENERGY_INDEX] = mass_flux * H_L;
    }
    else
    {
        // Upwind from right state
        face_flux[RHO_INDEX] = mass_flux;
        face_flux[RHU_INDEX] = mass_flux * u_R + P_interface * nx;
        face_flux[RHV_INDEX] = mass_flux * v_R + P_interface * ny;
        face_flux[RHW_INDEX] = mass_flux * w_R + P_interface * nz;
        face_flux[ENERGY_INDEX] = mass_flux * H_R;
    }

    // Scale by face area
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        face_flux[i] *= face_area;
    }

    // Apply AUSM+ enhancements if needed
    if (Enable_AUSM_Plus_Enhancements)
    {
        Apply_AUSM_Plus_Corrections_3D(left_cell_idx, right_cell_idx, face_idx,
                                       M_L, M_R, p_L, p_R, face_flux);
    }
}

/**
 * @brief Apply AUSM+ corrections for improved accuracy
 * @param left_cell Left cell index
 * @param right_cell Right cell index
 * @param face_idx Face index
 * @param M_L Left Mach number
 * @param M_R Right Mach number
 * @param p_L Left pressure
 * @param p_R Right pressure
 * @param face_flux Face flux to be modified
 */
void Apply_AUSM_Plus_Corrections_3D(const int &left_cell, const int &right_cell,
                                    const int &face_idx, const double &M_L, const double &M_R,
                                    const double &p_L, const double &p_R, V_D &face_flux)
{
    // AUSM+ pressure diffusion term
    const double Ku = 0.75;   // Pressure velocity coupling parameter
    const double sigma = 1.0; // Pressure diffusion parameter

    double face_area = Cells[left_cell].Face_Areas[face_idx];

    // Pressure diffusion correction
    double p_diff = sigma * fabs(p_R - p_L) / min(p_L, p_R);
    double M_p = -Ku * max(1.0 - sigma * max(M_L, M_R), 0.0) * p_diff;

    // Apply correction to momentum components
    double nx = Cells[left_cell].Face_Normals[3 * face_idx];
    double ny = Cells[left_cell].Face_Normals[3 * face_idx + 1];
    double nz = Cells[left_cell].Face_Normals[3 * face_idx + 2];

    face_flux[RHU_INDEX] += M_p * (p_R - p_L) * nx * face_area;
    face_flux[RHV_INDEX] += M_p * (p_R - p_L) * ny * face_area;
    face_flux[RHW_INDEX] += M_p * (p_R - p_L) * nz * face_area;
}

/**
 * @brief AUSM flux with preconditioning for low Mach numbers
 * @param Cell_No Current cell index
 */
void Ausm_Preconditioned_Flux_3D(const int &Cell_No)
{
    // Initialize net flux
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Calculate reference Mach number for preconditioning
    double M_ref = Calculate_Reference_Mach_3D(Cell_No);

    // Preconditioning parameter
    double beta = min(max(M_ref, 0.01), 1.0);
    double beta_sq = beta * beta;

    // Process each face with preconditioning
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor_idx = Cells[Cell_No].Neighbours[face];

        if (neighbor_idx >= 0 && neighbor_idx < No_Physical_Cells)
        {
            V_D face_flux(NUM_CONSERVATIVE_VARS, 0.0);

            // Compute base AUSM flux
            Compute_AUSM_Face_Flux_3D(Cell_No, neighbor_idx, face, face_flux);

            // Apply preconditioning matrix
            Apply_Preconditioning_3D(face_flux, beta_sq);

            // Add to net flux
            for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
            {
                Cells_Net_Flux[Cell_No][i] += face_flux[i];
            }
        }
    }

    // Normalize by cell volume
    double inv_volume = Cells[Cell_No].Inv_Volume;
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        Cells_Net_Flux[Cell_No][i] *= inv_volume;
    }
}

/**
 * @brief Calculate reference Mach number for preconditioning
 * @param Cell_No Cell index
 * @return Reference Mach number
 */
double Calculate_Reference_Mach_3D(const int &Cell_No)
{
    double rho = U_Cells[Cell_No][RHO_INDEX];
    double u = U_Cells[Cell_No][RHU_INDEX] / rho;
    double v = U_Cells[Cell_No][RHV_INDEX] / rho;
    double w = U_Cells[Cell_No][RHW_INDEX] / rho;

    double V_mag = sqrt(u * u + v * v + w * w);
    double p = (GAMMA - 1.0) * (U_Cells[Cell_No][ENERGY_INDEX] -
                                0.5 * rho * V_mag * V_mag);
    double a = sqrt(GAMMA * p / rho);

    return V_mag / a;
}

/**
 * @brief Apply preconditioning matrix to flux
 * @param flux Flux vector to be preconditioned
 * @param beta_sq Preconditioning parameter squared
 */
void Apply_Preconditioning_3D(V_D &flux, const double &beta_sq)
{
    // Simplified preconditioning - scale energy equation
    // More sophisticated preconditioning matrices can be implemented
    double scaling_factor = 1.0 / beta_sq;
    flux[ENERGY_INDEX] *= scaling_factor;
}

/**
 * @brief AUSM-DV (Dissipation Vector) scheme for 3D
 * @param Cell_No Current cell index
 */
void Ausm_DV_Flux_3D(const int &Cell_No)
{
    // Initialize net flux
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Process each face
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor_idx = Cells[Cell_No].Neighbours[face];

        if (neighbor_idx >= 0 && neighbor_idx < No_Physical_Cells)
        {
            V_D face_flux(NUM_CONSERVATIVE_VARS, 0.0);

            // Compute AUSM flux with dissipation vector
            Compute_AUSM_DV_Face_Flux_3D(Cell_No, neighbor_idx, face, face_flux);

            // Add to net flux
            for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
            {
                Cells_Net_Flux[Cell_No][i] += face_flux[i];
            }
        }
    }

    // Normalize by cell volume
    double inv_volume = Cells[Cell_No].Inv_Volume;
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        Cells_Net_Flux[Cell_No][i] *= inv_volume;
    }
}

/**
 * @brief Compute AUSM-DV flux with dissipation vector
 * @param left_cell Left cell index
 * @param right_cell Right cell index
 * @param face_idx Face index
 * @param face_flux Output flux
 */
void Compute_AUSM_DV_Face_Flux_3D(const int &left_cell, const int &right_cell,
                                  const int &face_idx, V_D &face_flux)
{
    // First compute standard AUSM flux
    Compute_AUSM_Face_Flux_3D(left_cell, right_cell, face_idx, face_flux);

    // Add dissipation vector corrections
    V_D dissipation(NUM_CONSERVATIVE_VARS, 0.0);
    Compute_Dissipation_Vector_3D(left_cell, right_cell, face_idx, dissipation);

    // Add dissipation to flux
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        face_flux[i] += dissipation[i];
    }
}

/**
 * @brief Compute dissipation vector for AUSM-DV
 * @param left_cell Left cell index
 * @param right_cell Right cell index
 * @param face_idx Face index
 * @param dissipation Output dissipation vector
 */
void Compute_Dissipation_Vector_3D(const int &left_cell, const int &right_cell,
                                   const int &face_idx, V_D &dissipation)
{
    // Get face area
    double face_area = Cells[left_cell].Face_Areas[face_idx];

    // Calculate pressure and velocity differences
    double p_L = (GAMMA - 1.0) * (U_Cells[left_cell][ENERGY_INDEX] -
                                  0.5 * (U_Cells[left_cell][RHU_INDEX] * U_Cells[left_cell][RHU_INDEX] + U_Cells[left_cell][RHV_INDEX] * U_Cells[left_cell][RHV_INDEX] + U_Cells[left_cell][RHW_INDEX] * U_Cells[left_cell][RHW_INDEX]) /
                                      U_Cells[left_cell][RHO_INDEX]);

    double p_R = (GAMMA - 1.0) * (U_Cells[right_cell][ENERGY_INDEX] -
                                  0.5 * (U_Cells[right_cell][RHU_INDEX] * U_Cells[right_cell][RHU_INDEX] + U_Cells[right_cell][RHV_INDEX] * U_Cells[right_cell][RHV_INDEX] + U_Cells[right_cell][RHW_INDEX] * U_Cells[right_cell][RHW_INDEX]) /
                                      U_Cells[right_cell][RHO_INDEX]);

    // Dissipation coefficients
    const double C_p = 0.25; // Pressure dissipation
    const double C_u = 0.75; // Velocity dissipation

    // Calculate dissipation terms
    double p_diff = C_p * fabs(p_R - p_L);

    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        double u_diff = U_Cells[right_cell][i] - U_Cells[left_cell][i];
        dissipation[i] = -C_u * fabs(u_diff) * face_area;
    }

    // Add pressure dissipation to momentum equations
    double nx = Cells[left_cell].Face_Normals[3 * face_idx];
    double ny = Cells[left_cell].Face_Normals[3 * face_idx + 1];
    double nz = Cells[left_cell].Face_Normals[3 * face_idx + 2];

    dissipation[RHU_INDEX] += p_diff * nx * face_area;
    dissipation[RHV_INDEX] += p_diff * ny * face_area;
    dissipation[RHW_INDEX] += p_diff * nz * face_area;
}