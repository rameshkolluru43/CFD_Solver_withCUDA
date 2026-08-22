#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"
#include "Utilities.h"

/**
 * @file LLF.cpp
 * @brief 3D Local Lax-Friedrichs (LLF) flux implementation
 *
 * This file implements the Local Lax-Friedrichs (Rusanov) flux scheme for 3D
 * Euler/Navier-Stokes equations. LLF is a simple and robust flux scheme that
 * provides good stability and shock-capturing capabilities.
 */

/**
 * @brief Local Lax-Friedrichs flux calculation for 3D hexahedral cell
 * @param Cell_No Current cell index
 *
 * Computes LLF flux for all 6 faces of a 3D hexahedral cell.
 * The method uses maximum wave speeds to compute artificial dissipation.
 */
void LLF_3D(const int &Cell_No)
{
    cout << "Evaluating 3D LLF Flux for Cell " << Cell_No << endl;

    // Initialize net flux for this cell
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Process each face of the hexahedral cell
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor_idx = Cells[Cell_No].Neighbours[face];

        // Skip boundary faces (handled separately)
        if (neighbor_idx < 0 || neighbor_idx >= No_Physical_Cells)
        {
            continue;
        }

        // Calculate LLF flux for this face
        V_D face_flux(NUM_CONSERVATIVE_VARS, 0.0);
        Compute_LLF_Face_Flux_3D(Cell_No, neighbor_idx, face, face_flux);

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
 * @brief Compute LLF flux across a single face in 3D
 * @param left_cell_idx Left cell index
 * @param right_cell_idx Right cell index
 * @param face_idx Face index
 * @param face_flux Output flux vector
 */
void Compute_LLF_Face_Flux_3D(const int &left_cell_idx, const int &right_cell_idx,
                              const int &face_idx, V_D &face_flux)
{
    // Initialize variables
    double max_eigen_value = 0.0;
    double u_L = 0.0, v_L = 0.0, w_L = 0.0, u_R = 0.0, v_R = 0.0, w_R = 0.0;
    double nx = 0.0, ny = 0.0, nz = 0.0;
    double Vdotn_L = 0.0, Vdotn_R = 0.0;
    double C_L = 0.0, C_R = 0.0;

    // Get face geometry
    nx = Cells[left_cell_idx].Face_Normals[3 * face_idx];
    ny = Cells[left_cell_idx].Face_Normals[3 * face_idx + 1];
    nz = Cells[left_cell_idx].Face_Normals[3 * face_idx + 2];
    double face_area = Cells[left_cell_idx].Face_Areas[face_idx];

    // Left state variables
    V_D U_L(NUM_CONSERVATIVE_VARS);
    U_L[RHO_INDEX] = U_Cells[left_cell_idx][RHO_INDEX];
    U_L[RHU_INDEX] = U_Cells[left_cell_idx][RHU_INDEX];
    U_L[RHV_INDEX] = U_Cells[left_cell_idx][RHV_INDEX];
    U_L[RHW_INDEX] = U_Cells[left_cell_idx][RHW_INDEX];
    U_L[ENERGY_INDEX] = U_Cells[left_cell_idx][ENERGY_INDEX];

    u_L = Primitive_Cells[left_cell_idx][PRIM_U];
    v_L = Primitive_Cells[left_cell_idx][PRIM_V];
    w_L = Primitive_Cells[left_cell_idx][PRIM_W];
    C_L = sqrt(GAMMA * Primitive_Cells[left_cell_idx][PRIM_P] /
               Primitive_Cells[left_cell_idx][PRIM_RHO]);

    // Right state variables
    V_D U_R(NUM_CONSERVATIVE_VARS);
    U_R[RHO_INDEX] = U_Cells[right_cell_idx][RHO_INDEX];
    U_R[RHU_INDEX] = U_Cells[right_cell_idx][RHU_INDEX];
    U_R[RHV_INDEX] = U_Cells[right_cell_idx][RHV_INDEX];
    U_R[RHW_INDEX] = U_Cells[right_cell_idx][RHW_INDEX];
    U_R[ENERGY_INDEX] = U_Cells[right_cell_idx][ENERGY_INDEX];

    u_R = Primitive_Cells[right_cell_idx][PRIM_U];
    v_R = Primitive_Cells[right_cell_idx][PRIM_V];
    w_R = Primitive_Cells[right_cell_idx][PRIM_W];
    C_R = sqrt(GAMMA * Primitive_Cells[right_cell_idx][PRIM_P] /
               Primitive_Cells[right_cell_idx][PRIM_RHO]);

    // Normal velocities
    Vdotn_L = u_L * nx + v_L * ny + w_L * nz;
    Vdotn_R = u_R * nx + v_R * ny + w_R * nz;

    // Wave speed evaluation for 3D
    V_D wave_speeds(6, 0.0);
    wave_speeds[0] = fabs(Vdotn_L - C_L); // Left acoustic wave (backward)
    wave_speeds[1] = fabs(Vdotn_L + C_L); // Left acoustic wave (forward)
    wave_speeds[2] = fabs(Vdotn_L);       // Left entropy/shear wave
    wave_speeds[3] = fabs(Vdotn_R - C_R); // Right acoustic wave (backward)
    wave_speeds[4] = fabs(Vdotn_R + C_R); // Right acoustic wave (forward)
    wave_speeds[5] = fabs(Vdotn_R);       // Right entropy/shear wave

    // Find maximum wave speed
    double Max_L = max({wave_speeds[0], wave_speeds[1], wave_speeds[2]});
    double Max_R = max({wave_speeds[3], wave_speeds[4], wave_speeds[5]});
    max_eigen_value = max(Max_L, Max_R);

    // Compute physical fluxes
    V_D Flux_L(NUM_CONSERVATIVE_VARS), Flux_R(NUM_CONSERVATIVE_VARS);
    Compute_Physical_Flux_3D(U_L, Primitive_Cells[left_cell_idx], nx, ny, nz, Flux_L);
    Compute_Physical_Flux_3D(U_R, Primitive_Cells[right_cell_idx], nx, ny, nz, Flux_R);

    // LLF flux: F_LLF = 0.5 * (F_L + F_R) - 0.5 * λ_max * (U_R - U_L)
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        face_flux[i] = 0.5 * (Flux_L[i] + Flux_R[i]) -
                       0.5 * max_eigen_value * (U_R[i] - U_L[i]);
        face_flux[i] *= face_area;
    }
}

/**
 * @brief Second-order LLF flux with MUSCL reconstruction
 * @param Cell_No Current cell index
 */
void LLF_Second_Order_3D(const int &Cell_No)
{
    cout << "Evaluating 2nd Order 3D LLF Flux for Cell " << Cell_No << endl;

    // Initialize net flux
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Process each face
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor_idx = Cells[Cell_No].Neighbours[face];

        if (neighbor_idx >= 0 && neighbor_idx < No_Physical_Cells)
        {
            V_D face_flux(NUM_CONSERVATIVE_VARS, 0.0);
            Compute_LLF_Second_Order_Face_Flux_3D(Cell_No, neighbor_idx, face, face_flux);

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
 * @brief Compute second-order LLF flux with reconstruction
 * @param left_cell Left cell index
 * @param right_cell Right cell index
 * @param face_idx Face index
 * @param face_flux Output flux
 */
void Compute_LLF_Second_Order_Face_Flux_3D(const int &left_cell, const int &right_cell,
                                           const int &face_idx, V_D &face_flux)
{
    // Get reconstructed left and right states
    V_D U_L_recon(NUM_CONSERVATIVE_VARS), U_R_recon(NUM_CONSERVATIVE_VARS);
    V_D Prim_L_recon(NUM_PRIMITIVE_VARS), Prim_R_recon(NUM_PRIMITIVE_VARS);

    // Base states
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        U_L_recon[i] = U_Cells[left_cell][i];
        U_R_recon[i] = U_Cells[right_cell][i];
    }
    for (int i = 0; i < NUM_PRIMITIVE_VARS; i++)
    {
        Prim_L_recon[i] = Primitive_Cells[left_cell][i];
        Prim_R_recon[i] = Primitive_Cells[right_cell][i];
    }

    // Apply MUSCL reconstruction if second-order is enabled
    if (Is_Second_Order)
    {
        Apply_MUSCL_Reconstruction_LLF_3D(left_cell, right_cell, face_idx,
                                          U_L_recon, U_R_recon,
                                          Prim_L_recon, Prim_R_recon);
    }

    // Compute LLF flux with reconstructed states
    Compute_LLF_Flux_With_States_3D(U_L_recon, U_R_recon, Prim_L_recon, Prim_R_recon,
                                    left_cell, face_idx, face_flux);
}

/**
 * @brief Apply MUSCL reconstruction for LLF scheme
 * @param left_cell Left cell index
 * @param right_cell Right cell index
 * @param face_idx Face index
 * @param U_L_recon Reconstructed left conservative variables
 * @param U_R_recon Reconstructed right conservative variables
 * @param Prim_L_recon Reconstructed left primitive variables
 * @param Prim_R_recon Reconstructed right primitive variables
 */
void Apply_MUSCL_Reconstruction_LLF_3D(const int &left_cell, const int &right_cell,
                                       const int &face_idx,
                                       V_D &U_L_recon, V_D &U_R_recon,
                                       V_D &Prim_L_recon, V_D &Prim_R_recon)
{
    // Get face center
    V_D face_center(3);
    face_center[0] = Cells[left_cell].Face_Centers[3 * face_idx];
    face_center[1] = Cells[left_cell].Face_Centers[3 * face_idx + 1];
    face_center[2] = Cells[left_cell].Face_Centers[3 * face_idx + 2];

    // Calculate gradients for left cell
    VV_D gradients_L(NUM_PRIMITIVE_VARS, V_D(3, 0.0));
    Calculate_Cell_Gradients_3D(left_cell, gradients_L);

    // Calculate gradients for right cell
    VV_D gradients_R(NUM_PRIMITIVE_VARS, V_D(3, 0.0));
    Calculate_Cell_Gradients_3D(right_cell, gradients_R);

    // Apply limiters
    VV_D limited_gradients_L(NUM_PRIMITIVE_VARS, V_D(3, 0.0));
    VV_D limited_gradients_R(NUM_PRIMITIVE_VARS, V_D(3, 0.0));

    Apply_Gradient_Limiters_3D(left_cell, gradients_L, limited_gradients_L);
    Apply_Gradient_Limiters_3D(right_cell, gradients_R, limited_gradients_R);

    // Reconstruct left state at face
    V_D dr_L(3);
    dr_L[0] = face_center[0] - Cells[left_cell].Cell_Center[0];
    dr_L[1] = face_center[1] - Cells[left_cell].Cell_Center[1];
    dr_L[2] = face_center[2] - Cells[left_cell].Cell_Center[2];

    for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
    {
        double gradient_contribution = DOT_PRODUCT_3D(limited_gradients_L[var], dr_L);
        Prim_L_recon[var] += gradient_contribution;
    }

    // Reconstruct right state at face
    V_D dr_R(3);
    dr_R[0] = face_center[0] - Cells[right_cell].Cell_Center[0];
    dr_R[1] = face_center[1] - Cells[right_cell].Cell_Center[1];
    dr_R[2] = face_center[2] - Cells[right_cell].Cell_Center[2];

    for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
    {
        double gradient_contribution = DOT_PRODUCT_3D(limited_gradients_R[var], dr_R);
        Prim_R_recon[var] += gradient_contribution;
    }

    // Convert back to conservative variables
    Convert_Primitive_to_Conservative_3D(Prim_L_recon, U_L_recon);
    Convert_Primitive_to_Conservative_3D(Prim_R_recon, U_R_recon);
}

/**
 * @brief Compute LLF flux with given left and right states
 * @param U_L Left conservative variables
 * @param U_R Right conservative variables
 * @param Prim_L Left primitive variables
 * @param Prim_R Right primitive variables
 * @param cell_idx Cell index for face geometry
 * @param face_idx Face index
 * @param face_flux Output flux
 */
void Compute_LLF_Flux_With_States_3D(const V_D &U_L, const V_D &U_R,
                                     const V_D &Prim_L, const V_D &Prim_R,
                                     const int &cell_idx, const int &face_idx,
                                     V_D &face_flux)
{
    // Get face geometry
    double nx = Cells[cell_idx].Face_Normals[3 * face_idx];
    double ny = Cells[cell_idx].Face_Normals[3 * face_idx + 1];
    double nz = Cells[cell_idx].Face_Normals[3 * face_idx + 2];
    double face_area = Cells[cell_idx].Face_Areas[face_idx];

    // Calculate wave speeds
    double u_L = Prim_L[PRIM_U], v_L = Prim_L[PRIM_V], w_L = Prim_L[PRIM_W];
    double u_R = Prim_R[PRIM_U], v_R = Prim_R[PRIM_V], w_R = Prim_R[PRIM_W];

    double C_L = sqrt(GAMMA * Prim_L[PRIM_P] / Prim_L[PRIM_RHO]);
    double C_R = sqrt(GAMMA * Prim_R[PRIM_P] / Prim_R[PRIM_RHO]);

    double Vn_L = u_L * nx + v_L * ny + w_L * nz;
    double Vn_R = u_R * nx + v_R * ny + w_R * nz;

    // Maximum wave speed
    double lambda_max = max({fabs(Vn_L) + C_L, fabs(Vn_R) + C_R});

    // Compute physical fluxes
    V_D Flux_L(NUM_CONSERVATIVE_VARS), Flux_R(NUM_CONSERVATIVE_VARS);
    Compute_Physical_Flux_3D(U_L, Prim_L, nx, ny, nz, Flux_L);
    Compute_Physical_Flux_3D(U_R, Prim_R, nx, ny, nz, Flux_R);

    // LLF flux
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        face_flux[i] = 0.5 * (Flux_L[i] + Flux_R[i]) -
                       0.5 * lambda_max * (U_R[i] - U_L[i]);
        face_flux[i] *= face_area;
    }
}

/**
 * @brief Calculate cell gradients using Green-Gauss method
 * @param cell_idx Cell index
 * @param gradients Output gradients for each primitive variable
 */
void Calculate_Cell_Gradients_3D(const int &cell_idx, VV_D &gradients)
{
    // Initialize gradients to zero
    for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
    {
        gradients[var].assign(3, 0.0);
    }

    double inv_volume = Cells[cell_idx].Inv_Volume;

    // Green-Gauss gradient calculation
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor = Cells[cell_idx].Neighbours[face];

        // Get face normal and area
        double nx = Cells[cell_idx].Face_Normals[3 * face];
        double ny = Cells[cell_idx].Face_Normals[3 * face + 1];
        double nz = Cells[cell_idx].Face_Normals[3 * face + 2];
        double area = Cells[cell_idx].Face_Areas[face];

        // Get face values (simple averaging)
        V_D face_values(NUM_PRIMITIVE_VARS);
        if (neighbor >= 0 && neighbor < No_Physical_Cells)
        {
            for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
            {
                face_values[var] = 0.5 * (Primitive_Cells[cell_idx][var] +
                                          Primitive_Cells[neighbor][var]);
            }
        }
        else
        {
            // Boundary face - use cell values
            for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
            {
                face_values[var] = Primitive_Cells[cell_idx][var];
            }
        }

        // Accumulate gradient contributions
        for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
        {
            gradients[var][0] += face_values[var] * nx * area;
            gradients[var][1] += face_values[var] * ny * area;
            gradients[var][2] += face_values[var] * nz * area;
        }
    }

    // Normalize by volume
    for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
    {
        gradients[var][0] *= inv_volume;
        gradients[var][1] *= inv_volume;
        gradients[var][2] *= inv_volume;
    }
}

/**
 * @brief Apply gradient limiters to prevent spurious oscillations
 * @param cell_idx Cell index
 * @param gradients Input gradients
 * @param limited_gradients Output limited gradients
 */
void Apply_Gradient_Limiters_3D(const int &cell_idx, const VV_D &gradients,
                                VV_D &limited_gradients)
{
    for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
    {
        // Calculate limiter based on neighboring cell values
        double limiter = Calculate_Gradient_Limiter_3D(cell_idx, var, gradients[var]);

        // Apply limiter
        limited_gradients[var][0] = limiter * gradients[var][0];
        limited_gradients[var][1] = limiter * gradients[var][1];
        limited_gradients[var][2] = limiter * gradients[var][2];
    }
}

/**
 * @brief Calculate gradient limiter for a specific variable
 * @param cell_idx Cell index
 * @param var_idx Variable index
 * @param gradient Gradient vector
 * @return Limiter value (0 to 1)
 */
double Calculate_Gradient_Limiter_3D(const int &cell_idx, const int &var_idx,
                                     const V_D &gradient)
{
    // Simplified Barth-Jespersen limiter
    double phi = 1.0;
    double cell_value = Primitive_Cells[cell_idx][var_idx];

    // Find min and max values among neighbors
    double min_val = cell_value;
    double max_val = cell_value;

    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor = Cells[cell_idx].Neighbours[face];
        if (neighbor >= 0 && neighbor < No_Physical_Cells)
        {
            double neighbor_val = Primitive_Cells[neighbor][var_idx];
            min_val = min(min_val, neighbor_val);
            max_val = max(max_val, neighbor_val);
        }
    }

    // Check each face for monotonicity
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Calculate reconstructed value at face
        V_D dr(3);
        dr[0] = Cells[cell_idx].Face_Centers[3 * face] - Cells[cell_idx].Cell_Center[0];
        dr[1] = Cells[cell_idx].Face_Centers[3 * face + 1] - Cells[cell_idx].Cell_Center[1];
        dr[2] = Cells[cell_idx].Face_Centers[3 * face + 2] - Cells[cell_idx].Cell_Center[2];

        double reconstructed_val = cell_value + DOT_PRODUCT_3D(gradient, dr);

        // Apply Barth-Jespersen limiting
        if (reconstructed_val > cell_value)
        {
            if (max_val > cell_value)
            {
                phi = min(phi, (max_val - cell_value) / (reconstructed_val - cell_value));
            }
            else
            {
                phi = 0.0;
            }
        }
        else if (reconstructed_val < cell_value)
        {
            if (min_val < cell_value)
            {
                phi = min(phi, (min_val - cell_value) / (reconstructed_val - cell_value));
            }
            else
            {
                phi = 0.0;
            }
        }
    }

    return max(0.0, phi);
}