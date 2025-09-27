#include "definitions.h"
#include "Globals.h"
#include "Numerical_Schemes.h"
#include "Primitive_Computational.h"
#include "Utilities.h"

/**
 * @file Roe_Scheme.cpp
 * @brief 3D Roe approximate Riemann solver implementation
 *
 * This file implements the Roe flux-difference splitting scheme for 3D Euler/Navier-Stokes equations.
 * The scheme uses Roe-averaged states to compute eigenvalues and eigenvectors for upwind flux calculation.
 */

/**
 * @brief Roe flux computation for 3D Euler equations
 * @param Cell_No Current cell index
 * @param Face_No Current face index (0-5 for hexahedral cells)
 * @param Net_Flux Output flux vector for the face
 *
 * Computes Roe flux across a face in 3D using linearized Riemann solver approach.
 * Handles all three coordinate directions with proper eigenvalue decomposition.
 */
void Roe_Flux_3D(const int &Cell_No, const int &Face_No, V_D &Net_Flux)
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

    Get_Left_Right_States_3D(Cell_No, Face_No, U_L, U_R, Primitive_L, Primitive_R);

    // Left state variables
    double rho_L = Primitive_L[PRIM_RHO];
    double u_L = Primitive_L[PRIM_U];
    double v_L = Primitive_L[PRIM_V];
    double w_L = Primitive_L[PRIM_W];
    double p_L = Primitive_L[PRIM_P];
    double c_L = sqrt(GAMMA * p_L / rho_L);
    double H_L = (U_L[ENERGY_INDEX] + p_L) / rho_L;

    // Right state variables
    double rho_R = Primitive_R[PRIM_RHO];
    double u_R = Primitive_R[PRIM_U];
    double v_R = Primitive_R[PRIM_V];
    double w_R = Primitive_R[PRIM_W];
    double p_R = Primitive_R[PRIM_P];
    double c_R = sqrt(GAMMA * p_R / rho_R);
    double H_R = (U_R[ENERGY_INDEX] + p_R) / rho_R;

    // Face-normal velocities
    double Vn_L = u_L * nx + v_L * ny + w_L * nz;
    double Vn_R = u_R * nx + v_R * ny + w_R * nz;

    // Compute Roe averages
    double sqrt_rho_L = sqrt(rho_L);
    double sqrt_rho_R = sqrt(rho_R);
    double sum_sqrt_rho = sqrt_rho_L + sqrt_rho_R;

    double rho_roe = sqrt_rho_L * sqrt_rho_R;
    double u_roe = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / sum_sqrt_rho;
    double v_roe = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) / sum_sqrt_rho;
    double w_roe = (sqrt_rho_L * w_L + sqrt_rho_R * w_R) / sum_sqrt_rho;
    double H_roe = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / sum_sqrt_rho;

    double Vn_roe = u_roe * nx + v_roe * ny + w_roe * nz;
    double V2_roe = u_roe * u_roe + v_roe * v_roe + w_roe * w_roe;
    double c_roe = sqrt((GAMMA - 1.0) * (H_roe - 0.5 * V2_roe));

    // Compute eigenvalues
    double lambda1 = Vn_roe - c_roe;
    double lambda2 = Vn_roe;
    double lambda3 = Vn_roe;
    double lambda4 = Vn_roe;
    double lambda5 = Vn_roe + c_roe;

    // Apply entropy fix
    if (Enable_Entropy_Fix)
    {
        Apply_Roe_Entropy_Fix_3D(lambda1, lambda2, lambda3, lambda4, lambda5,
                                 Vn_L, Vn_R, c_L, c_R);
    }

    // Compute flux differences
    V_D dU(NUM_CONSERVATIVE_VARS);
    dU[RHO_INDEX] = U_R[RHO_INDEX] - U_L[RHO_INDEX];
    dU[RHU_INDEX] = U_R[RHU_INDEX] - U_L[RHU_INDEX];
    dU[RHV_INDEX] = U_R[RHV_INDEX] - U_L[RHV_INDEX];
    dU[RHW_INDEX] = U_R[RHW_INDEX] - U_L[RHW_INDEX];
    dU[ENERGY_INDEX] = U_R[ENERGY_INDEX] - U_L[ENERGY_INDEX];

    // Compute wave strengths (characteristic variables)
    V_D alpha(NUM_CONSERVATIVE_VARS);
    Compute_Wave_Strengths_3D(dU, rho_roe, u_roe, v_roe, w_roe, H_roe, c_roe,
                              nx, ny, nz, alpha);

    // Compute dissipation terms
    V_D dissipation(NUM_CONSERVATIVE_VARS, 0.0);

    // Wave 1: rho - c (acoustic wave moving left)
    if (fabs(lambda1) > EPSILON)
    {
        double coeff = fabs(lambda1) * alpha[0];
        dissipation[RHO_INDEX] += coeff;
        dissipation[RHU_INDEX] += coeff * (u_roe - c_roe * nx);
        dissipation[RHV_INDEX] += coeff * (v_roe - c_roe * ny);
        dissipation[RHW_INDEX] += coeff * (w_roe - c_roe * nz);
        dissipation[ENERGY_INDEX] += coeff * (H_roe - c_roe * Vn_roe);
    }

    // Wave 2: rho (entropy wave)
    if (fabs(lambda2) > EPSILON)
    {
        double coeff = fabs(lambda2) * alpha[1];
        dissipation[RHO_INDEX] += coeff;
        dissipation[RHU_INDEX] += coeff * u_roe;
        dissipation[RHV_INDEX] += coeff * v_roe;
        dissipation[RHW_INDEX] += coeff * w_roe;
        dissipation[ENERGY_INDEX] += coeff * 0.5 * V2_roe;
    }

    // Wave 3: rho*u_t1 (shear wave 1)
    if (fabs(lambda3) > EPSILON)
    {
        double coeff = fabs(lambda3) * alpha[2];
        V_D t1(3);
        Compute_Tangent_Vector_1_3D(nx, ny, nz, t1);

        dissipation[RHU_INDEX] += coeff * t1[0];
        dissipation[RHV_INDEX] += coeff * t1[1];
        dissipation[RHW_INDEX] += coeff * t1[2];
        dissipation[ENERGY_INDEX] += coeff * (u_roe * t1[0] + v_roe * t1[1] + w_roe * t1[2]);
    }

    // Wave 4: rho*u_t2 (shear wave 2)
    if (fabs(lambda4) > EPSILON)
    {
        double coeff = fabs(lambda4) * alpha[3];
        V_D t2(3);
        Compute_Tangent_Vector_2_3D(nx, ny, nz, t2);

        dissipation[RHU_INDEX] += coeff * t2[0];
        dissipation[RHV_INDEX] += coeff * t2[1];
        dissipation[RHW_INDEX] += coeff * t2[2];
        dissipation[ENERGY_INDEX] += coeff * (u_roe * t2[0] + v_roe * t2[1] + w_roe * t2[2]);
    }

    // Wave 5: rho + c (acoustic wave moving right)
    if (fabs(lambda5) > EPSILON)
    {
        double coeff = fabs(lambda5) * alpha[4];
        dissipation[RHO_INDEX] += coeff;
        dissipation[RHU_INDEX] += coeff * (u_roe + c_roe * nx);
        dissipation[RHV_INDEX] += coeff * (v_roe + c_roe * ny);
        dissipation[RHW_INDEX] += coeff * (w_roe + c_roe * nz);
        dissipation[ENERGY_INDEX] += coeff * (H_roe + c_roe * Vn_roe);
    }

    // Compute physical fluxes
    V_D Flux_L(NUM_CONSERVATIVE_VARS), Flux_R(NUM_CONSERVATIVE_VARS);
    Compute_Physical_Flux_3D(U_L, Primitive_L, nx, ny, nz, Flux_L);
    Compute_Physical_Flux_3D(U_R, Primitive_R, nx, ny, nz, Flux_R);

    // Combine fluxes: F_roe = 0.5 * (F_L + F_R) - 0.5 * |A| * (U_R - U_L)
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        Net_Flux[i] = 0.5 * (Flux_L[i] + Flux_R[i]) - 0.5 * dissipation[i];
        Net_Flux[i] *= face_area;
    }
}

/**
 * @brief Compute wave strengths for Roe scheme in 3D
 * @param dU Conservative variable differences
 * @param rho_roe Roe-averaged density
 * @param u_roe Roe-averaged u-velocity
 * @param v_roe Roe-averaged v-velocity
 * @param w_roe Roe-averaged w-velocity
 * @param H_roe Roe-averaged enthalpy
 * @param c_roe Roe-averaged speed of sound
 * @param nx Face normal x-component
 * @param ny Face normal y-component
 * @param nz Face normal z-component
 * @param alpha Output wave strengths
 */
void Compute_Wave_Strengths_3D(const V_D &dU, const double &rho_roe,
                               const double &u_roe, const double &v_roe, const double &w_roe,
                               const double &H_roe, const double &c_roe,
                               const double &nx, const double &ny, const double &nz,
                               V_D &alpha)
{
    double Vn_roe = u_roe * nx + v_roe * ny + w_roe * nz;
    double V2_roe = u_roe * u_roe + v_roe * v_roe + w_roe * w_roe;

    // Compute tangent vectors
    V_D t1(3), t2(3);
    Compute_Tangent_Vector_1_3D(nx, ny, nz, t1);
    Compute_Tangent_Vector_2_3D(nx, ny, nz, t2);

    double dVn = (dU[RHU_INDEX] * nx + dU[RHV_INDEX] * ny + dU[RHW_INDEX] * nz) / rho_roe;
    double dVt1 = (dU[RHU_INDEX] * t1[0] + dU[RHV_INDEX] * t1[1] + dU[RHW_INDEX] * t1[2]) / rho_roe;
    double dVt2 = (dU[RHU_INDEX] * t2[0] + dU[RHV_INDEX] * t2[1] + dU[RHW_INDEX] * t2[2]) / rho_roe;

    double drho = dU[RHO_INDEX];
    double dE = dU[ENERGY_INDEX];

    // Wave strengths
    alpha[0] = 0.5 * (drho * (Vn_roe + c_roe) - (dU[RHU_INDEX] * nx + dU[RHV_INDEX] * ny + dU[RHW_INDEX] * nz) - c_roe * (dE - 0.5 * rho_roe * V2_roe * drho / rho_roe)) / c_roe;
    alpha[1] = drho - (dE - 0.5 * rho_roe * V2_roe * drho / rho_roe) / (c_roe * c_roe);
    alpha[2] = rho_roe * dVt1;
    alpha[3] = rho_roe * dVt2;
    alpha[4] = 0.5 * (drho * (Vn_roe - c_roe) + (dU[RHU_INDEX] * nx + dU[RHV_INDEX] * ny + dU[RHW_INDEX] * nz) - c_roe * (dE - 0.5 * rho_roe * V2_roe * drho / rho_roe)) / c_roe;
}

/**
 * @brief Compute first tangent vector orthogonal to face normal
 * @param nx Face normal x-component
 * @param ny Face normal y-component
 * @param nz Face normal z-component
 * @param t1 Output first tangent vector
 */
void Compute_Tangent_Vector_1_3D(const double &nx, const double &ny, const double &nz, V_D &t1)
{
    // Find a vector not parallel to normal
    if (fabs(nx) < 0.9)
    {
        t1[0] = 1.0;
        t1[1] = 0.0;
        t1[2] = 0.0;
    }
    else
    {
        t1[0] = 0.0;
        t1[1] = 1.0;
        t1[2] = 0.0;
    }

    // Compute cross product with normal to get tangent
    V_D n = {nx, ny, nz};
    V_D temp = t1;
    CROSS_PRODUCT_3D(t1, n, temp);
    NORMALIZE_3D(t1);
}

/**
 * @brief Compute second tangent vector orthogonal to normal and first tangent
 * @param nx Face normal x-component
 * @param ny Face normal y-component
 * @param nz Face normal z-component
 * @param t2 Output second tangent vector
 */
void Compute_Tangent_Vector_2_3D(const double &nx, const double &ny, const double &nz, V_D &t2)
{
    V_D n = {nx, ny, nz};
    V_D t1(3);
    Compute_Tangent_Vector_1_3D(nx, ny, nz, t1);

    // Second tangent is cross product of normal and first tangent
    CROSS_PRODUCT_3D(t2, n, t1);
    NORMALIZE_3D(t2);
}

/**
 * @brief Compute physical flux in 3D
 * @param U Conservative variables
 * @param Prim Primitive variables
 * @param nx Face normal x-component
 * @param ny Face normal y-component
 * @param nz Face normal z-component
 * @param Flux Output flux vector
 */
void Compute_Physical_Flux_3D(const V_D &U, const V_D &Prim,
                              const double &nx, const double &ny, const double &nz,
                              V_D &Flux)
{
    double rho = Prim[PRIM_RHO];
    double u = Prim[PRIM_U];
    double v = Prim[PRIM_V];
    double w = Prim[PRIM_W];
    double p = Prim[PRIM_P];

    double Vn = u * nx + v * ny + w * nz;

    Flux[RHO_INDEX] = rho * Vn;
    Flux[RHU_INDEX] = rho * u * Vn + p * nx;
    Flux[RHV_INDEX] = rho * v * Vn + p * ny;
    Flux[RHW_INDEX] = rho * w * Vn + p * nz;
    Flux[ENERGY_INDEX] = (U[ENERGY_INDEX] + p) * Vn;
}

/**
 * @brief Apply entropy fix to Roe eigenvalues
 * @param lambda1 First eigenvalue (u-c)
 * @param lambda2 Second eigenvalue (u)
 * @param lambda3 Third eigenvalue (u)
 * @param lambda4 Fourth eigenvalue (u)
 * @param lambda5 Fifth eigenvalue (u+c)
 * @param Vn_L Left state normal velocity
 * @param Vn_R Right state normal velocity
 * @param c_L Left state speed of sound
 * @param c_R Right state speed of sound
 */
void Apply_Roe_Entropy_Fix_3D(double &lambda1, double &lambda2, double &lambda3,
                              double &lambda4, double &lambda5,
                              const double &Vn_L, const double &Vn_R,
                              const double &c_L, const double &c_R)
{
    const double epsilon = 0.1;

    // Harten's entropy fix for acoustic waves
    double lambda_L1 = Vn_L - c_L;
    double lambda_R1 = Vn_R - c_R;
    double lambda_L5 = Vn_L + c_L;
    double lambda_R5 = Vn_R + c_R;

    // Fix for first acoustic wave
    if (lambda1 < lambda_R1 && lambda1 > lambda_L1)
    {
        double delta = max(0.0, max(lambda_R1 - lambda1, lambda1 - lambda_L1));
        lambda1 = 0.5 * (lambda1 * lambda1 / delta + delta);
    }
    else
    {
        lambda1 = fabs(lambda1);
    }

    // Fix for second acoustic wave
    if (lambda5 < lambda_R5 && lambda5 > lambda_L5)
    {
        double delta = max(0.0, max(lambda_R5 - lambda5, lambda5 - lambda_L5));
        lambda5 = 0.5 * (lambda5 * lambda5 / delta + delta);
    }
    else
    {
        lambda5 = fabs(lambda5);
    }

    // Absolute values for entropy and shear waves
    lambda2 = fabs(lambda2);
    lambda3 = fabs(lambda3);
    lambda4 = fabs(lambda4);
}

/**
 * @brief Roe flux computation for all faces of a 3D cell
 * @param Cell_No Current cell index
 */
void Roe_Cell_Flux_3D(const int &Cell_No)
{
    V_D Face_Flux(NUM_CONSERVATIVE_VARS);

    // Initialize net flux for this cell
    Cells_Net_Flux[Cell_No].assign(NUM_CONSERVATIVE_VARS, 0.0);

    // Compute flux for each of the 6 faces
    for (int Face_No = 0; Face_No < NUM_FACES_3D; Face_No++)
    {
        Roe_Flux_3D(Cell_No, Face_No, Face_Flux);

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
 * @brief Second-order Roe scheme with MUSCL reconstruction
 * @param Cell_No Current cell index
 * @param Face_No Current face index
 * @param Net_Flux Output flux vector
 */
void Roe_Second_Order_3D(const int &Cell_No, const int &Face_No, V_D &Net_Flux)
{
    // Get reconstructed left and right states
    V_D U_L(NUM_CONSERVATIVE_VARS), U_R(NUM_CONSERVATIVE_VARS);
    V_D Prim_L(NUM_PRIMITIVE_VARS), Prim_R(NUM_PRIMITIVE_VARS);

    Get_Left_Right_States_3D(Cell_No, Face_No, U_L, U_R, Prim_L, Prim_R);

    // Apply MUSCL reconstruction
    if (Is_Second_Order)
    {
        Apply_Second_Order_Reconstruction_3D(Cell_No, Face_No, U_L, U_R, Prim_L, Prim_R);
    }

    // Use the standard Roe flux with reconstructed states
    // The reconstruction is already applied in Get_Left_Right_States_3D
    Roe_Flux_3D(Cell_No, Face_No, Net_Flux);
}

/**
 * @brief MUSCL reconstruction for higher-order accuracy
 * @param Cell_No Current cell index
 * @param gradients Cell gradients for reconstruction
 * @param limited_gradients Output limited gradients
 */
void Apply_MUSCL_Reconstruction_3D(const int &Cell_No, const VV_D &gradients, VV_D &limited_gradients)
{
    limited_gradients.resize(NUM_PRIMITIVE_VARS);
    for (int var = 0; var < NUM_PRIMITIVE_VARS; var++)
    {
        limited_gradients[var].resize(3, 0.0);

        // Apply limiter to each gradient component
        for (int dir = 0; dir < 3; dir++)
        {
            double unlimited_grad = gradients[var][dir];
            double limiter = Calculate_Limiter_3D(Cell_No, 0); // Simplified
            limited_gradients[var][dir] = limiter * unlimited_grad;
        }
    }
}