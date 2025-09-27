#include "definitions.h"
#include "Globals.h"
#include "Weno.h"
#include "Primitive_Computational.h"
#include "Utilities.h"
#include "Timestep.h"
#include <cmath>    // For std::isfinite and other math functions
#include <iostream> // For error messages

/**
 * @file WENO3D.cpp
 * @brief 3D WENO (Weighted Essentially Non-Oscillatory) scheme implementation
 *
 * This file implements the WENO scheme for 3D hexahedral cells, providing
 * high-order accurate reconstruction with essentially non-oscillatory properties.
 * The scheme can operate in both conservative and characteristic variable spaces.
 */

/**
 * @brief Fifth-order WENO reconstruction for 3D applications
 * @param a First stencil point
 * @param b Second stencil point
 * @param c Third stencil point (central)
 * @param d Fourth stencil point
 * @param e Fifth stencil point
 * @param shift Direction flag (0 for left, 1 for right)
 * @param U Output reconstructed value
 */
void WENO_Reconstruction_3D(double &a, double &b, double &c, double &d, double &e, int &shift, double &U)
{
    // Check for NaN or infinite values in input
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d) || !std::isfinite(e))
    {
        std::cout << "Warning: Non-finite values in 3D WENO reconstruction input" << std::endl;
        U = c; // Fall back to central value
        return;
    }

    double d0 = 0.0, d1 = 0.0, d2 = 0.0, b0 = 0.0, b1 = 0.0, b2 = 0.0, a0 = 0.0, a1 = 0.0, a2 = 0.0;
    double v0 = 0.0, v1 = 0.0, v2 = 0.0, w0 = 0.0, w1 = 0.0, w2 = 0.0, sum = 0.0;
    double epsilon = 1e-6; // WENO parameter for avoiding division by zero

    int p = 2; // Order of accuracy parameter

    switch (shift)
    {
    case 0: // for Left Values (upwind-biased)
        d0 = 3.0 / 10.0;
        d1 = 6.0 / 10.0;
        d2 = 1.0 / 10.0;

        v2 = (2.0 * a - 7.0 * b + 11.0 * c) / 6.0;
        v1 = (-1.0 * b + 5.0 * c + 2.0 * d) / 6.0;
        v0 = (2.0 * c + 5.0 * d - 1.0 * e) / 6.0;
        break;

    case 1: // for Right Values (downwind-biased)
        d0 = 1.0 / 10.0;
        d1 = 6.0 / 10.0;
        d2 = 3.0 / 10.0;

        v2 = (-1.0 * a + 5.0 * b + 2.0 * c) / 6.0;
        v1 = (2.0 * b + 5.0 * c - 1.0 * d) / 6.0;
        v0 = (11.0 * c - 7.0 * d + 2.0 * e) / 6.0;
        break;
    }

    // Smoothness indicators (Jiang and Shu, 1996)
    b2 = (13.0 / 12.0) * pow((a - 2.0 * b + c), 2) + (1.0 / 4.0) * pow((a - 4.0 * b + 3.0 * c), 2);
    b1 = (13.0 / 12.0) * pow((b - 2.0 * c + d), 2) + (1.0 / 4.0) * pow((b - d), 2);
    b0 = (13.0 / 12.0) * pow((c - 2.0 * d + e), 2) + (1.0 / 4.0) * pow((3.0 * c - 4.0 * d + e), 2);

    // WENO weights
    a0 = d0 / pow((epsilon + b0), p);
    a1 = d1 / pow((epsilon + b1), p);
    a2 = d2 / pow((epsilon + b2), p);

    sum = a0 + a1 + a2;

    // Check for degenerate case where sum is too small
    if (sum < epsilon)
    {
        // Fall back to simple average of the three candidate values
        U = (v0 + v1 + v2) / 3.0;
    }
    else
    {
        w0 = a0 / sum;
        w1 = a1 / sum;
        w2 = a2 / sum;
        U = w0 * v0 + w1 * v1 + w2 * v2;
    }

    // Final validation
    if (!std::isfinite(U))
    {
        std::cout << "Warning: Non-finite result in 3D WENO reconstruction, using central value" << std::endl;
        U = c;
    }
}

/**
 * @brief Compute left and right eigenvector matrices for 3D characteristic decomposition
 * @param CellNo Current cell index
 * @param Neighbour Neighbor cell index
 * @param Face_No Face index (0-5 for hexahedral cells)
 * @param L Left eigenvector matrix (25 elements, 5x5)
 * @param IL Right eigenvector matrix (25 elements, 5x5)
 */
void Get_LR_3D(int &CellNo, int &Neighbour, const int &Face_No, V_D &L, V_D &IL)
{
    int index = Face_No * 3; // 3D normal vector index

    double dL = 0.0, uL = 0.0, vL = 0.0, wL = 0.0, aL = 0.0;
    double dR = 0.0, uR = 0.0, vR = 0.0, wR = 0.0, aR = 0.0;
    double h = 0.0, ek = 0.0, vn = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;
    double a_RL = 0.0, u_RL = 0.0, v_RL = 0.0, w_RL = 0.0, t1 = 0.0, t2 = 0.0;

    // Left state primitive variables
    dL = Primitive_Cells[CellNo][PRIM_RHO];
    uL = Primitive_Cells[CellNo][PRIM_U];
    vL = Primitive_Cells[CellNo][PRIM_V];
    wL = Primitive_Cells[CellNo][PRIM_W];
    aL = Primitive_Cells[CellNo][PRIM_C]; // Speed of sound

    // Right state primitive variables
    dR = Primitive_Cells[Neighbour][PRIM_RHO];
    uR = Primitive_Cells[Neighbour][PRIM_U];
    vR = Primitive_Cells[Neighbour][PRIM_V];
    wR = Primitive_Cells[Neighbour][PRIM_W];
    aR = Primitive_Cells[Neighbour][PRIM_C];

    // Roe averaging for 3D
    double sqrt_dL = sqrt(fmax(dL, 1e-14));
    double sqrt_dR = sqrt(fmax(dR, 1e-14));
    double denom = sqrt_dL + sqrt_dR;

    if (denom < 1e-14)
    {
        // Handle degenerate case - use simple average
        u_RL = 0.5 * (uL + uR);
        v_RL = 0.5 * (vL + vR);
        w_RL = 0.5 * (wL + wR);
        a_RL = 0.5 * (aL + aR);
    }
    else
    {
        u_RL = (uL * sqrt_dL + uR * sqrt_dR) / denom;
        v_RL = (vL * sqrt_dL + vR * sqrt_dR) / denom;
        w_RL = (wL * sqrt_dL + wR * sqrt_dR) / denom;
        a_RL = (aL * sqrt_dL + aR * sqrt_dR) / denom;
    }

    // Face normal components
    nx = Cells[CellNo].Face_Normals[index + 0];
    ny = Cells[CellNo].Face_Normals[index + 1];
    nz = Cells[CellNo].Face_Normals[index + 2];

    // Normal velocity and kinetic energy
    vn = (u_RL * nx + v_RL * ny + w_RL * nz);
    ek = 0.5 * (u_RL * u_RL + v_RL * v_RL + w_RL * w_RL);
    h = (a_RL * a_RL / gamma_M_1) + ek; // Total enthalpy

    // Add safety check for speed of sound
    if (a_RL < 1e-14)
    {
        std::cout << "Warning: Very small speed of sound in 3D WENO, a_RL = " << a_RL << std::endl;
        a_RL = 1e-14;
    }

    t1 = 0.5 / (a_RL * a_RL);
    t2 = gamma_M_1 * t1;

    // Compute tangent vectors for 3D (Gram-Schmidt orthogonalization)
    V_D t1_vec(3), t2_vec(3);

    // First tangent vector
    if (fabs(nx) < 0.9)
    {
        t1_vec[0] = 0.0;
        t1_vec[1] = nz;
        t1_vec[2] = -ny;
    }
    else
    {
        t1_vec[0] = -nz;
        t1_vec[1] = 0.0;
        t1_vec[2] = nx;
    }

    // Normalize first tangent vector
    double t1_mag = sqrt(t1_vec[0] * t1_vec[0] + t1_vec[1] * t1_vec[1] + t1_vec[2] * t1_vec[2]);
    if (t1_mag > 1e-14)
    {
        t1_vec[0] /= t1_mag;
        t1_vec[1] /= t1_mag;
        t1_vec[2] /= t1_mag;
    }

    // Second tangent vector (cross product of normal and first tangent)
    t2_vec[0] = ny * t1_vec[2] - nz * t1_vec[1];
    t2_vec[1] = nz * t1_vec[0] - nx * t1_vec[2];
    t2_vec[2] = nx * t1_vec[1] - ny * t1_vec[0];

    // 3D Left eigenvector matrix (Conservative to Characteristic)
    // Row 0: Entropy wave
    L[0] = 1.0 - 2.0 * t2 * ek;
    L[5] = 2.0 * t2 * u_RL;
    L[10] = 2.0 * t2 * v_RL;
    L[15] = 2.0 * t2 * w_RL;
    L[20] = -2.0 * t2;

    // Row 1: First shear wave
    L[1] = u_RL * t1_vec[0] + v_RL * t1_vec[1] + w_RL * t1_vec[2];
    L[6] = t1_vec[0];
    L[11] = t1_vec[1];
    L[16] = t1_vec[2];
    L[21] = 0.0;

    // Row 2: Second shear wave
    L[2] = u_RL * t2_vec[0] + v_RL * t2_vec[1] + w_RL * t2_vec[2];
    L[7] = t2_vec[0];
    L[12] = t2_vec[1];
    L[17] = t2_vec[2];
    L[22] = 0.0;

    // Row 3: Acoustic wave (negative)
    L[3] = t2 * ek - a_RL * vn * t1;
    L[8] = -t2 * u_RL + a_RL * t1 * nx;
    L[13] = -t2 * v_RL + a_RL * t1 * ny;
    L[18] = -t2 * w_RL + a_RL * t1 * nz;
    L[23] = t2;

    // Row 4: Acoustic wave (positive)
    L[4] = t2 * ek + a_RL * vn * t1;
    L[9] = -t2 * u_RL - a_RL * t1 * nx;
    L[14] = -t2 * v_RL - a_RL * t1 * ny;
    L[19] = -t2 * w_RL - a_RL * t1 * nz;
    L[24] = t2;

    // Validate L matrix
    for (int i = 0; i < 25; i++)
    {
        if (isnan(L[i]))
        {
            cout << "NaN occurred in 3D L matrix at index " << i << "\t" << L[i] << endl;
            exit(0);
        }
    }

    // 3D Right eigenvector matrix (Characteristic to Conservative)
    // Column 0: Entropy wave
    IL[0] = 1.0;
    IL[1] = u_RL;
    IL[2] = v_RL;
    IL[3] = w_RL;
    IL[4] = ek;

    // Column 1: First shear wave
    IL[5] = 0.0;
    IL[6] = t1_vec[0];
    IL[7] = t1_vec[1];
    IL[8] = t1_vec[2];
    IL[9] = u_RL * t1_vec[0] + v_RL * t1_vec[1] + w_RL * t1_vec[2];

    // Column 2: Second shear wave
    IL[10] = 0.0;
    IL[11] = t2_vec[0];
    IL[12] = t2_vec[1];
    IL[13] = t2_vec[2];
    IL[14] = u_RL * t2_vec[0] + v_RL * t2_vec[1] + w_RL * t2_vec[2];

    // Column 3: Acoustic wave (negative)
    IL[15] = 1.0;
    IL[16] = u_RL - a_RL * nx;
    IL[17] = v_RL - a_RL * ny;
    IL[18] = w_RL - a_RL * nz;
    IL[19] = h - a_RL * vn;

    // Column 4: Acoustic wave (positive)
    IL[20] = 1.0;
    IL[21] = u_RL + a_RL * nx;
    IL[22] = v_RL + a_RL * ny;
    IL[23] = w_RL + a_RL * nz;
    IL[24] = h + a_RL * vn;

    // Validate IL matrix
    for (int i = 0; i < 25; i++)
    {
        if (isnan(IL[i]))
        {
            cout << "NaN occurred in 3D IL matrix at index " << i << "\t" << IL[i] << endl;
            exit(0);
        }
    }
}

/**
 * @brief Matrix-vector multiplication for 3D (5x5 matrix with 5-element vector)
 * @param A Matrix (25 elements stored column-wise)
 * @param x Input vector (5 elements)
 * @param B Output vector (5 elements)
 */
void MatVecMul_3D(V_D &A, V_D &x, V_D &B)
{
    int size = 5; // 3D conservative variables
    for (int i = 0; i < size; i++)
    {
        B[i] = 0.0;
        for (int j = 0; j < size; j++)
        {
            B[i] += A[i + size * j] * x[j];
        }
    }
}

/**
 * @brief Get reconstructed conservative variables using WENO for 3D
 * @param Cell_No Current cell index
 * @param Face_No Face index (0-5)
 * @param i1, i2, i3, i4, i5 Stencil cell indices
 * @param LR Left/Right flag (0 for left, 1 for right)
 * @param U Output reconstructed conservative variables
 */
void Get_Reconstructed_U_3D(int &Cell_No, const int &Face_No, int &i1, int &i2, int &i3, int &i4, int &i5, int &LR, V_D &U)
{
    int Neighbour = 0;

    // Get neighbor based on face
    if (Face_No >= 0 && Face_No < NUM_FACES_3D)
    {
        Neighbour = Cells[Cell_No].Neighbours[Face_No];
    }
    else
    {
        cout << "ERROR: Invalid face number " << Face_No << " in Get_Reconstructed_U_3D" << endl;
        exit(1);
    }

    V_D U1(NUM_CONSERVATIVE_VARS, 0.0), U2(NUM_CONSERVATIVE_VARS, 0.0), U3(NUM_CONSERVATIVE_VARS, 0.0);
    V_D U4(NUM_CONSERVATIVE_VARS, 0.0), U5(NUM_CONSERVATIVE_VARS, 0.0);
    V_D W1(NUM_CONSERVATIVE_VARS, 0.0), W2(NUM_CONSERVATIVE_VARS, 0.0), W3(NUM_CONSERVATIVE_VARS, 0.0);
    V_D W4(NUM_CONSERVATIVE_VARS, 0.0), W5(NUM_CONSERVATIVE_VARS, 0.0), W(NUM_CONSERVATIVE_VARS, 0.0);

    V_D L(25, 0.0), InvL(25, 0.0); // 5x5 matrices for 3D

    // Get conservative variables from stencil cells
    Calculate_Computational_Variables_3D(i1, U1);
    Calculate_Computational_Variables_3D(i2, U2);
    Calculate_Computational_Variables_3D(i3, U3);
    Calculate_Computational_Variables_3D(i4, U4);
    Calculate_Computational_Variables_3D(i5, U5);

    if (Is_Char)
    {
        // Characteristic variable reconstruction
        Get_LR_3D(Cell_No, Neighbour, Face_No, L, InvL);
        MatVecMul_3D(L, U1, W1);
        MatVecMul_3D(L, U2, W2);
        MatVecMul_3D(L, U3, W3);
        MatVecMul_3D(L, U4, W4);
        MatVecMul_3D(L, U5, W5);

        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            WENO_Reconstruction_3D(W1[i], W2[i], W3[i], W4[i], W5[i], LR, W[i]);
        }

        MatVecMul_3D(InvL, W, U);
    }
    else
    {
        // Conservative variable reconstruction
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            WENO_Reconstruction_3D(U1[i], U2[i], U3[i], U4[i], U5[i], LR, U[i]);
        }
    }
}

/**
 * @brief WENO reconstruction for 3D hexahedral cells along specified direction
 * @param Cell_No Current cell index
 * @param Face_No Face index (0-5)
 * @param Direction Coordinate direction (0=x, 1=y, 2=z)
 * @param U_L Left reconstructed state
 * @param U_R Right reconstructed state
 */
void WENO_Reconstruction_Direction_3D(int &Cell_No, const int &Face_No, const int &Direction, V_D &U_L, V_D &U_R)
{
    // Initialize neighbor indices
    int im3 = 0, im2 = 0, im1 = 0, ip1 = 0, ip2 = 0, ip3 = 0;
    int LR = 0;

    // Get immediate neighbors
    im1 = Cells[Cell_No].Neighbours[2 * Direction];     // Negative direction neighbor
    ip1 = Cells[Cell_No].Neighbours[2 * Direction + 1]; // Positive direction neighbor

    // Build stencil in negative direction
    if (im1 >= No_Physical_Cells || im1 < 0)
    {
        im2 = im1;
        im3 = im1;
    }
    else
    {
        im2 = Cells[im1].Neighbours[2 * Direction];
        if (im2 >= No_Physical_Cells || im2 < 0)
        {
            im3 = im2;
        }
        else
        {
            im3 = Cells[im2].Neighbours[2 * Direction];
            if (im3 < 0)
                im3 = im2;
        }
    }

    // Build stencil in positive direction
    if (ip1 >= No_Physical_Cells || ip1 < 0)
    {
        ip2 = ip1;
        ip3 = ip1;
    }
    else
    {
        ip2 = Cells[ip1].Neighbours[2 * Direction + 1];
        if (ip2 >= No_Physical_Cells || ip2 < 0)
        {
            ip3 = ip2;
        }
        else
        {
            ip3 = Cells[ip2].Neighbours[2 * Direction + 1];
            if (ip3 < 0)
                ip3 = ip2;
        }
    }

    // Perform WENO reconstruction
    switch (Face_No)
    {
    case 2 * Direction: // Left face in direction
        LR = 0;
        Get_Reconstructed_U_3D(Cell_No, Face_No, im3, im2, im1, Cell_No, ip1, LR, U_L);
        LR = 1;
        Get_Reconstructed_U_3D(Cell_No, Face_No, im2, im1, Cell_No, ip1, ip2, LR, U_R);
        break;

    case 2 * Direction + 1: // Right face in direction
        LR = 0;
        Get_Reconstructed_U_3D(Cell_No, Face_No, im2, im1, Cell_No, ip1, ip2, LR, U_L);
        LR = 1;
        Get_Reconstructed_U_3D(Cell_No, Face_No, im1, Cell_No, ip1, ip2, ip3, LR, U_R);
        break;

    default:
        cout << "ERROR: Invalid face-direction combination in WENO_Reconstruction_Direction_3D" << endl;
        exit(1);
    }
}

/**
 * @brief Main WENO reconstruction function for 3D hexahedral cells
 * @param Cell_No Current cell index
 * @param Face_No Face index (0-5)
 * @param U_L Left reconstructed state
 * @param U_R Right reconstructed state
 */
void WENO_Reconstruction_3D_Main(int &Cell_No, const int &Face_No, V_D &U_L, V_D &U_R)
{
    // Determine coordinate direction based on face
    int Direction = Face_No / 2; // 0=x-direction, 1=y-direction, 2=z-direction

    WENO_Reconstruction_Direction_3D(Cell_No, Face_No, Direction, U_L, U_R);
}

/**
 * @brief Calculate face flux using WENO reconstruction for 3D
 * @param Cell_No Current cell index
 * @param N_Cell_No Neighbor cell index
 * @param Face_No Face index
 * @param Is_Wall_Face Wall boundary flag
 */
void Calculate_Face_WENO_Flux_3D(int &Cell_No, int &N_Cell_No, const int &Face_No, bool Is_Wall_Face)
{
    // Initialize flux variables
    double Rho_L = 0.0, T_L = 0.0, P_L = 0.0, u_L = 0.0, v_L = 0.0, w_L = 0.0;
    double Vdotn_L = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, Vmag_L = 0.0, C_L = 0.0;
    double Rho_R = 0.0, T_R = 0.0, P_R = 0.0, u_R = 0.0, v_R = 0.0, w_R = 0.0;
    double Vdotn_R = 0.0, Vmag_R = 0.0, C_R = 0.0, face_area = 0.0;

    V_D Flux_L(NUM_CONSERVATIVE_VARS, 0.0), Flux_R(NUM_CONSERVATIVE_VARS, 0.0);
    V_D Average_Convective_Flux(NUM_CONSERVATIVE_VARS, 0.0), Dissipative_Flux(NUM_CONSERVATIVE_VARS, 0.0);
    V_D U_L(NUM_CONSERVATIVE_VARS, 0.0), U_R(NUM_CONSERVATIVE_VARS, 0.0);
    V_D S(6, 0.0); // Wave speeds

    int index = Face_No * 3; // 3D normal vector index

    // WENO reconstruction
    WENO_Reconstruction_3D_Main(Cell_No, Face_No, U_L, U_R);

    // Calculate primitive variables from reconstructed states
    Calculate_Primitive_Variables_3D(Cell_No, U_L);
    Rho_L = Global_Primitive_3D[PRIM_RHO];
    P_L = Global_Primitive_3D[PRIM_P];
    u_L = Global_Primitive_3D[PRIM_U];
    v_L = Global_Primitive_3D[PRIM_V];
    w_L = Global_Primitive_3D[PRIM_W];
    C_L = Global_Primitive_3D[PRIM_C];

    Calculate_Primitive_Variables_3D(N_Cell_No, U_R);
    Rho_R = Global_Primitive_3D[PRIM_RHO];
    P_R = Global_Primitive_3D[PRIM_P];
    u_R = Global_Primitive_3D[PRIM_U];
    v_R = Global_Primitive_3D[PRIM_V];
    w_R = Global_Primitive_3D[PRIM_W];
    C_R = Global_Primitive_3D[PRIM_C];

    // Face geometry
    nx = Cells[Cell_No].Face_Normals[index + 0];
    ny = Cells[Cell_No].Face_Normals[index + 1];
    nz = Cells[Cell_No].Face_Normals[index + 2];
    face_area = Cells[Cell_No].Face_Areas[Face_No];

    // Normal velocities and kinetic energies
    Vdotn_L = (u_L * nx + v_L * ny + w_L * nz);
    Vdotn_R = (u_R * nx + v_R * ny + w_R * nz);
    Vmag_L = 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);
    Vmag_R = 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);

    // Wall boundary condition
    if (Is_Wall_Face)
    {
        P_R = P_L;
    }

    // Left state flux
    Flux_L[RHO_INDEX] = Rho_L * Vdotn_L * face_area;
    Flux_L[RHU_INDEX] = Rho_L * u_L * Vdotn_L * face_area + P_L * nx * face_area;
    Flux_L[RHV_INDEX] = Rho_L * v_L * Vdotn_L * face_area + P_L * ny * face_area;
    Flux_L[RHW_INDEX] = Rho_L * w_L * Vdotn_L * face_area + P_L * nz * face_area;
    Flux_L[ENERGY_INDEX] = (gamma1 * P_L + Rho_L * Vmag_L) * Vdotn_L * face_area;

    // Right state flux
    Flux_R[RHO_INDEX] = Rho_R * Vdotn_R * face_area;
    Flux_R[RHU_INDEX] = Rho_R * u_R * Vdotn_R * face_area + P_R * nx * face_area;
    Flux_R[RHV_INDEX] = Rho_R * v_R * Vdotn_R * face_area + P_R * ny * face_area;
    Flux_R[RHW_INDEX] = Rho_R * w_R * Vdotn_R * face_area + P_R * nz * face_area;
    Flux_R[ENERGY_INDEX] = (gamma1 * P_R + Rho_R * Vmag_R) * Vdotn_R * face_area;

    // Wave speed evaluation
    if (C_L < 1e-14)
        C_L = 1e-14;
    if (C_R < 1e-14)
        C_R = 1e-14;

    S[0] = fabs(Vdotn_L - C_L) * face_area;
    S[1] = fabs(Vdotn_L + C_L) * face_area;
    S[2] = fabs(Vdotn_L) * face_area;
    S[3] = fabs(Vdotn_R - C_R) * face_area;
    S[4] = fabs(Vdotn_R + C_R) * face_area;
    S[5] = fabs(Vdotn_R) * face_area;

    // Maximum wave speed
    double Max1 = max({S[0], S[1], S[2]});
    double Max2 = max({S[3], S[4], S[5]});
    double max_eigen_value = max(Max1, Max2);

    // Final flux calculation
    for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
    {
        Average_Convective_Flux[i] = 0.5 * (Flux_L[i] + Flux_R[i]);
        Dissipative_Flux[i] = 0.5 * max_eigen_value * (U_R[i] - U_L[i]);

        // Add to cell net flux
        Cells_Net_Flux[Cell_No][i] += Average_Convective_Flux[i] - Dissipative_Flux[i];
    }
}

/**
 * @brief Evaluate net flux for all cells using 3D WENO scheme
 */
void Evaluate_Cell_Net_Flux_WENO_3D()
{
    for (int Current_Cell_No = 0; Current_Cell_No < No_Physical_Cells; Current_Cell_No++)
    {
        // Initialize net flux
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            Cells_Net_Flux[Current_Cell_No][i] = 0.0;
        }

        // Process all 6 faces of hexahedral cell
        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            int N_Cell_No = Cells[Current_Cell_No].Neighbours[face];

            if (N_Cell_No >= 0) // Valid neighbor
            {
                Calculate_Face_WENO_Flux_3D(Current_Cell_No, N_Cell_No, face,
                                            Cells_Face_Boundary_Type[Current_Cell_No][face]);
            }
        }

        // Normalize by cell volume
        double inv_volume = Cells[Current_Cell_No].Inv_Volume;
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++)
        {
            Cells_Net_Flux[Current_Cell_No][i] *= inv_volume;
        }

        // Evaluate time step for this cell
        Evaluate_Time_Step_3D(Current_Cell_No);
    }
}