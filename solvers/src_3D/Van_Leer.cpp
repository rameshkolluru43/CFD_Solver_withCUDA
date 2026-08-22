#include "definitions.h"
#include "Globals.h"
#include "Numerical_Schemes.h"
#include "Primitive_Computational.h"
#include "Utilities.h"

/**
 * @file Van_Leer.cpp
 * @brief 3D Van Leer flux vector splitting scheme implementation
 *
 * This file implements the Van Leer flux vector splitting method for 3D Euler/Navier-Stokes equations.
 * The scheme splits convective fluxes based on the sign of velocity components and Mach number,
 * providing upwind bias and maintaining positivity of density and internal energy.
 */

/**
 * @brief Van Leer flux vector splitting for 3D Euler equations
 * @param Cell_No Current cell index
 * @param Face_No Current face index (0-5 for hexahedral cells)
 * @param Net_Flux Output flux vector for the face
 *
 * Computes Van Leer flux across a face in 3D using flux vector splitting approach.
 * Handles all three coordinate directions (x, y, z) with appropriate face normal projections.
 */
void Van_Leer_Flux_3D(const int &Cell_No, const int &Face_No, V_D &Net_Flux)
{
    // Reset flux vector
    Net_Flux.assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Get face properties
    double face_area = Cells[Cell_No].Face_Areas[Face_No];
    double nx = Cells[Cell_No].Face_Normals[3 * Face_No];
    double ny = Cells[Cell_No].Face_Normals[3 * Face_No + 1];
    double nz = Cells[Cell_No].Face_Normals[3 * Face_No + 2];

    // Get left and right states
    V_D U_L(NUM_CONSERVATIVE_VARS), U_R(NUM_CONSERVATIVE_VARS);
    V_D Primitive_L(NUM_PRIMITIVE_VARS), Primitive_R(NUM_PRIMITIVE_VARS);

    // Extract states from current cell and neighbor
    Get_Left_Right_States_3D(Cell_No, Face_No, U_L, U_R, Primitive_L, Primitive_R);

    // Left state variables
    double rho_L = Primitive_L[PRIM_RHO];
    double u_L = Primitive_L[PRIM_U];
    double v_L = Primitive_L[PRIM_V];
    double w_L = Primitive_L[PRIM_W];
    double p_L = Primitive_L[PRIM_P];
    double c_L = sqrt(GAMMA * p_L / rho_L);
    double H_L = (U_L[ENERGY_INDEX] + p_L) / rho_L; // Specific enthalpy

    // Right state variables
    double rho_R = Primitive_R[PRIM_RHO];
    double u_R = Primitive_R[PRIM_U];
    double v_R = Primitive_R[PRIM_V];
    double w_R = Primitive_R[PRIM_W];
    double p_R = Primitive_R[PRIM_P];
    double c_R = sqrt(GAMMA * p_R / rho_R);
    double H_R = (U_R[ENERGY_INDEX] + p_R) / rho_R;

    // Compute face-normal velocities
    double Vn_L = u_L * nx + v_L * ny + w_L * nz;
    double Vn_R = u_R * nx + v_R * ny + w_R * nz;

    // Compute Mach numbers in face-normal direction
    double M_L = Vn_L / c_L;
    double M_R = Vn_R / c_R;

    // Van Leer flux splitting
    V_D Flux_L_plus(NUM_CONSERVATIVE_VARS, 0.0);
    V_D Flux_R_minus(NUM_CONSERVATIVE_VARS, 0.0);

    // Left state contribution (positive flux)
    if (M_L >= 1.0)
    {
        // Supersonic outflow - take full flux
        Flux_L_plus[RHO_INDEX] = rho_L * Vn_L;
        Flux_L_plus[RHU_INDEX] = rho_L * u_L * Vn_L + p_L * nx;
        Flux_L_plus[RHV_INDEX] = rho_L * v_L * Vn_L + p_L * ny;
        Flux_L_plus[RHW_INDEX] = rho_L * w_L * Vn_L + p_L * nz;
        Flux_L_plus[ENERGY_INDEX] = rho_L * H_L * Vn_L;
    }
    else if (M_L <= -1.0)
    {
        // Supersonic inflow - no contribution
        // Flux_L_plus remains zero
    }
    else
    {
        // Subsonic case - Van Leer splitting
        double M_plus = 0.25 * (M_L + 1.0) * (M_L + 1.0);
        double factor = rho_L * c_L * M_plus;

        Flux_L_plus[RHO_INDEX] = factor;
        Flux_L_plus[RHU_INDEX] = factor * (u_L + (2.0 - M_L) * c_L * nx / GAMMA);
        Flux_L_plus[RHV_INDEX] = factor * (v_L + (2.0 - M_L) * c_L * ny / GAMMA);
        Flux_L_plus[RHW_INDEX] = factor * (w_L + (2.0 - M_L) * c_L * nz / GAMMA);
        Flux_L_plus[ENERGY_INDEX] = factor * (H_L + (2.0 - M_L) * (2.0 - M_L) * c_L * c_L / (2.0 * (GAMMA - 1.0)));
    }

    // Right state contribution (negative flux)
    if (M_R <= -1.0)
    {
        // Supersonic inflow - take full flux
        Flux_R_minus[RHO_INDEX] = rho_R * Vn_R;
        Flux_R_minus[RHU_INDEX] = rho_R * u_R * Vn_R + p_R * nx;
        Flux_R_minus[RHV_INDEX] = rho_R * v_R * Vn_R + p_R * ny;
        Flux_R_minus[RHW_INDEX] = rho_R * w_R * Vn_R + p_R * nz;
        Flux_R_minus[ENERGY_INDEX] = rho_R * H_R * Vn_R;
    }
    else if (M_R >= 1.0)
    {
        // Supersonic outflow - no contribution
        // Flux_R_minus remains zero
    }
    else
    {
        // Subsonic case - Van Leer splitting
        double M_minus = -0.25 * (M_R - 1.0) * (M_R - 1.0);
        double factor = rho_R * c_R * M_minus;

        Flux_R_minus[RHO_INDEX] = factor;
        Flux_R_minus[RHU_INDEX] = factor * (u_R + (2.0 + M_R) * c_R * nx / GAMMA);
        Flux_R_minus[RHV_INDEX] = factor * (v_R + (2.0 + M_R) * c_R * ny / GAMMA);
        Flux_R_minus[RHW_INDEX] = factor * (w_R + (2.0 + M_R) * c_R * nz / GAMMA);
        Flux_R_minus[ENERGY_INDEX] = factor * (H_R + (2.0 + M_R) * (2.0 + M_R) * c_R * c_R / (2.0 * (GAMMA - 1.0)));
    }

    // Combine left and right contributions
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        Net_Flux[i] = (Flux_L_plus[i] + Flux_R_minus[i]) * face_area;
    }

    // Apply entropy fix if needed
    if (Enable_Entropy_Fix)
    {
        Apply_Entropy_Fix_Van_Leer_3D(M_L, M_R, c_L, c_R, Net_Flux, face_area);
    }
}

/**
 * @brief Van Leer flux computation for all faces of a 3D cell
 * @param Cell_No Current cell index
 */
void Van_Leer_Cell_Flux_3D(const int &Cell_No)
{
    V_D Face_Flux(NUM_CONSERVATIVE_VARS);

    // Initialize net flux for this cell
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Compute flux for each of the 6 faces
    for (int Face_No = 0; Face_No < NUM_FACES_3D; Face_No++)
    {
        Van_Leer_Flux_3D(Cell_No, Face_No, Face_Flux);

        // Add contribution to net flux (considering face orientation)
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            Cells_Net_Flux[Cell_No][i] += Face_Flux[i];
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
 * @brief Get left and right states for flux calculation at a face
 * @param Cell_No Current cell index
 * @param Face_No Current face index
 * @param U_L Left state conservative variables
 * @param U_R Right state conservative variables
 * @param Prim_L Left state primitive variables
 * @param Prim_R Right state primitive variables
 */
void Get_Left_Right_States_3D(const int &Cell_No, const int &Face_No,
                              V_D &U_L, V_D &U_R, V_D &Prim_L, V_D &Prim_R)
{
    // Left state is current cell
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        U_L[i] = U_Cells[Cell_No][i];
    }
    for (int i = 0; i < NUM_PRIMITIVE_VARS; i++)
    {
        Prim_L[i] = Primitive_Cells[Cell_No][i];
    }

    // Right state is neighboring cell
    int neighbor_idx = Cells[Cell_No].Neighbours[Face_No];

    if (neighbor_idx >= 0 && neighbor_idx < Total_No_Cells)
    {
        // Interior face - use neighbor values
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            U_R[i] = U_Cells[neighbor_idx][i];
        }
        for (int i = 0; i < NUM_PRIMITIVE_VARS; i++)
        {
            Prim_R[i] = Primitive_Cells[neighbor_idx][i];
        }
    }
    else
    {
        // Boundary face - apply boundary conditions
        Apply_Boundary_Condition_3D(Cell_No, Face_No, U_R, Prim_R);
    }

    // Apply higher-order reconstruction if enabled
    if (Is_Second_Order)
    {
        Apply_Second_Order_Reconstruction_3D(Cell_No, Face_No, U_L, U_R, Prim_L, Prim_R);
    }
}

/**
 * @brief Apply entropy fix to Van Leer flux
 * @param M_L Left state Mach number
 * @param M_R Right state Mach number
 * @param c_L Left state speed of sound
 * @param c_R Right state speed of sound
 * @param Net_Flux Flux vector to be modified
 * @param face_area Face area
 */
void Apply_Entropy_Fix_Van_Leer_3D(const double &M_L, const double &M_R,
                                   const double &c_L, const double &c_R,
                                   V_D &Net_Flux, const double &face_area)
{
    // Entropy fix parameter
    const double epsilon = 0.1;

    // Check if entropy fix is needed (near sonic conditions)
    if (fabs(M_L) < 1.0 + epsilon && fabs(M_R) < 1.0 + epsilon)
    {
        // Apply Harten's entropy fix
        double lambda_avg = 0.5 * (c_L + c_R);
        double entropy_correction = epsilon * lambda_avg * face_area;

        // Apply correction to mass flux
        if (fabs(Net_Flux[RHO_INDEX]) < entropy_correction)
        {
            double sign = (Net_Flux[RHO_INDEX] >= 0.0) ? 1.0 : -1.0;
            Net_Flux[RHO_INDEX] = sign * entropy_correction;
        }
    }
}

/**
 * @brief Apply second-order reconstruction for Van Leer scheme
 * @param Cell_No Current cell index
 * @param Face_No Current face index
 * @param U_L Left state conservative variables (modified)
 * @param U_R Right state conservative variables (modified)
 * @param Prim_L Left state primitive variables (modified)
 * @param Prim_R Right state primitive variables (modified)
 */
void Apply_Second_Order_Reconstruction_3D(const int &Cell_No, const int &Face_No,
                                          V_D &U_L, V_D &U_R, V_D &Prim_L, V_D &Prim_R)
{
    // Get face center and cell centers for reconstruction
    V_D face_center(3);
    face_center[0] = Cells[Cell_No].Face_Centers[3 * Face_No];
    face_center[1] = Cells[Cell_No].Face_Centers[3 * Face_No + 1];
    face_center[2] = Cells[Cell_No].Face_Centers[3 * Face_No + 2];

    V_D cell_center_L = Cells[Cell_No].Cell_Center;

    // Calculate gradients for left cell
    V_D grad_rho_L(3), grad_u_L(3), grad_v_L(3), grad_w_L(3), grad_p_L(3);
    Calculate_Gradients_3D(Cell_No, grad_rho_L, grad_u_L, grad_v_L, grad_w_L, grad_p_L);

    // Reconstruct left state at face center
    V_D dr_L(3);
    dr_L[0] = face_center[0] - cell_center_L[0];
    dr_L[1] = face_center[1] - cell_center_L[1];
    dr_L[2] = face_center[2] - cell_center_L[2];

    // Apply limiter
    double limiter_L = Calculate_Limiter_3D(Cell_No, Face_No);

    // Update primitive variables with limited gradients
    Prim_L[PRIM_RHO] += limiter_L * DOT_PRODUCT_3D(grad_rho_L, dr_L);
    Prim_L[PRIM_U] += limiter_L * DOT_PRODUCT_3D(grad_u_L, dr_L);
    Prim_L[PRIM_V] += limiter_L * DOT_PRODUCT_3D(grad_v_L, dr_L);
    Prim_L[PRIM_W] += limiter_L * DOT_PRODUCT_3D(grad_w_L, dr_L);
    Prim_L[PRIM_P] += limiter_L * DOT_PRODUCT_3D(grad_p_L, dr_L);

    // Convert back to conservative variables
    Convert_Primitive_to_Conservative_3D(Prim_L, U_L);

    // Similar reconstruction for right state
    int neighbor_idx = Cells[Cell_No].Neighbours[Face_No];
    if (neighbor_idx >= 0 && neighbor_idx < No_Physical_Cells)
    {
        V_D cell_center_R = Cells[neighbor_idx].Cell_Center;
        V_D grad_rho_R(3), grad_u_R(3), grad_v_R(3), grad_w_R(3), grad_p_R(3);
        Calculate_Gradients_3D(neighbor_idx, grad_rho_R, grad_u_R, grad_v_R, grad_w_R, grad_p_R);

        V_D dr_R(3);
        dr_R[0] = face_center[0] - cell_center_R[0];
        dr_R[1] = face_center[1] - cell_center_R[1];
        dr_R[2] = face_center[2] - cell_center_R[2];

        double limiter_R = Calculate_Limiter_3D(neighbor_idx, Face_No);

        Prim_R[PRIM_RHO] += limiter_R * DOT_PRODUCT_3D(grad_rho_R, dr_R);
        Prim_R[PRIM_U] += limiter_R * DOT_PRODUCT_3D(grad_u_R, dr_R);
        Prim_R[PRIM_V] += limiter_R * DOT_PRODUCT_3D(grad_v_R, dr_R);
        Prim_R[PRIM_W] += limiter_R * DOT_PRODUCT_3D(grad_w_R, dr_R);
        Prim_R[PRIM_P] += limiter_R * DOT_PRODUCT_3D(grad_p_R, dr_R);

        Convert_Primitive_to_Conservative_3D(Prim_R, U_R);
    }
}

/**
 * @brief Convert primitive to conservative variables in 3D
 * @param Prim Primitive variables [rho, u, v, w, p]
 * @param U Conservative variables [rho, rho*u, rho*v, rho*w, rho*E]
 */
void Convert_Primitive_to_Conservative_3D(const V_D &Prim, V_D &U)
{
    double rho = Prim[PRIM_RHO];
    double u = Prim[PRIM_U];
    double v = Prim[PRIM_V];
    double w = Prim[PRIM_W];
    double p = Prim[PRIM_P];

    U[RHO_INDEX] = rho;
    U[RHU_INDEX] = rho * u;
    U[RHV_INDEX] = rho * v;
    U[RHW_INDEX] = rho * w;
    U[ENERGY_INDEX] = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

/**
 * @brief Calculate 3D gradients for a cell
 * @param Cell_No Cell index
 * @param grad_rho Density gradient
 * @param grad_u U-velocity gradient
 * @param grad_v V-velocity gradient
 * @param grad_w W-velocity gradient
 * @param grad_p Pressure gradient
 */
void Calculate_Gradients_3D(const int &Cell_No, V_D &grad_rho, V_D &grad_u,
                            V_D &grad_v, V_D &grad_w, V_D &grad_p)
{
    // Green-Gauss gradient calculation using face values
    grad_rho.assign(3, 0.0);
    grad_u.assign(3, 0.0);
    grad_v.assign(3, 0.0);
    grad_w.assign(3, 0.0);
    grad_p.assign(3, 0.0);

    double inv_volume = Cells[Cell_No].Inv_Volume;

    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Get face normal and area
        double nx = Cells[Cell_No].Face_Normals[3 * face];
        double ny = Cells[Cell_No].Face_Normals[3 * face + 1];
        double nz = Cells[Cell_No].Face_Normals[3 * face + 2];
        double area = Cells[Cell_No].Face_Areas[face];

        // Get face values (average of left and right states)
        int neighbor = Cells[Cell_No].Neighbours[face];
        double rho_face, u_face, v_face, w_face, p_face;

        if (neighbor >= 0 && neighbor < No_Physical_Cells)
        {
            rho_face = 0.5 * (Primitive_Cells[Cell_No][PRIM_RHO] + Primitive_Cells[neighbor][PRIM_RHO]);
            u_face = 0.5 * (Primitive_Cells[Cell_No][PRIM_U] + Primitive_Cells[neighbor][PRIM_U]);
            v_face = 0.5 * (Primitive_Cells[Cell_No][PRIM_V] + Primitive_Cells[neighbor][PRIM_V]);
            w_face = 0.5 * (Primitive_Cells[Cell_No][PRIM_W] + Primitive_Cells[neighbor][PRIM_W]);
            p_face = 0.5 * (Primitive_Cells[Cell_No][PRIM_P] + Primitive_Cells[neighbor][PRIM_P]);
        }
        else
        {
            // Boundary face - use cell values
            rho_face = Primitive_Cells[Cell_No][PRIM_RHO];
            u_face = Primitive_Cells[Cell_No][PRIM_U];
            v_face = Primitive_Cells[Cell_No][PRIM_V];
            w_face = Primitive_Cells[Cell_No][PRIM_W];
            p_face = Primitive_Cells[Cell_No][PRIM_P];
        }

        // Accumulate gradient contributions
        grad_rho[0] += rho_face * nx * area;
        grad_rho[1] += rho_face * ny * area;
        grad_rho[2] += rho_face * nz * area;

        grad_u[0] += u_face * nx * area;
        grad_u[1] += u_face * ny * area;
        grad_u[2] += u_face * nz * area;

        grad_v[0] += v_face * nx * area;
        grad_v[1] += v_face * ny * area;
        grad_v[2] += v_face * nz * area;

        grad_w[0] += w_face * nx * area;
        grad_w[1] += w_face * ny * area;
        grad_w[2] += w_face * nz * area;

        grad_p[0] += p_face * nx * area;
        grad_p[1] += p_face * ny * area;
        grad_p[2] += p_face * nz * area;
    }

    // Normalize by volume
    for (int i = 0; i < 3; i++)
    {
        grad_rho[i] *= inv_volume;
        grad_u[i] *= inv_volume;
        grad_v[i] *= inv_volume;
        grad_w[i] *= inv_volume;
        grad_p[i] *= inv_volume;
    }
}

/**
 * @brief Calculate limiter function for 3D reconstruction
 * @param Cell_No Cell index
 * @param Face_No Face index
 * @return Limiter value (0 to 1)
 */
double Calculate_Limiter_3D(const int &Cell_No, const int &Face_No)
{
    // Simplified limiter - can be extended to various limiter types
    switch (Limiter_Case)
    {
    case LIMITER_NONE:
        return 1.0;

    case LIMITER_MINMOD:
        return Calculate_Minmod_Limiter_3D(Cell_No, Face_No);

    case LIMITER_VAN_LEER:
        return Calculate_Van_Leer_Limiter_3D(Cell_No, Face_No);

    default:
        return 1.0;
    }
}

/**
 * @brief Calculate Van Leer limiter in 3D
 * @param Cell_No Cell index
 * @param Face_No Face index
 * @return Van Leer limiter value
 */
double Calculate_Van_Leer_Limiter_3D(const int &Cell_No, const int &Face_No)
{
    // This is a simplified implementation
    // In practice, you'd compute gradients and apply Van Leer limiting

    int neighbor = Cells[Cell_No].Neighbours[Face_No];
    if (neighbor < 0 || neighbor >= No_Physical_Cells)
        return 1.0;

    // Calculate pressure ratio as limiting parameter
    double p_cell = Primitive_Cells[Cell_No][PRIM_P];
    double p_neighbor = Primitive_Cells[neighbor][PRIM_P];

    double r = p_neighbor / (p_cell + SMALL_NUMBER);

    // Van Leer limiter function
    return (r + fabs(r)) / (1.0 + fabs(r));
}

/**
 * @brief Calculate Minmod limiter in 3D
 * @param Cell_No Cell index
 * @param Face_No Face index
 * @return Minmod limiter value
 */
double Calculate_Minmod_Limiter_3D(const int &Cell_No, const int &Face_No)
{
    // Simplified Minmod limiter
    return 0.5; // Conservative choice
}