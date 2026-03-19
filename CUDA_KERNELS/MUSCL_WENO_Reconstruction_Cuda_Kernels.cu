/**
 * @file MUSCL_WENO_Reconstruction_Cuda_Kernels.cu
 * @brief CUDA implementation of MUSCL and WENO reconstruction schemes
 * @author AI Assistant
 * @date 2026-01-15
 * @version 1.0
 * 
 * @details Implements:
 * - MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
 * - WENO5 (5th-order Weighted Essentially Non-Oscillatory)
 */

#include <cuda_runtime.h>
#include <math.h>

//=============================================================================
// MUSCL RECONSTRUCTION KERNEL
//=============================================================================

/**
 * @brief MUSCL reconstruction with slope limiting
 * 
 * Computes 2nd-order accurate left and right states at cell faces
 * using linear reconstruction with slope limiters.
 * 
 * Q_L = Q_i + 0.5 * phi_L * (Q_i - Q_{i-1})
 * Q_R = Q_j - 0.5 * phi_R * (Q_{j+1} - Q_j)
 */
__global__ void muscl_reconstruction_kernel(
    const double* U_cells,
    const int* cell_neighbors,
    const double* cell_distances,
    const double* limited_gradients,  // From limiter kernels
    double* Q_left,   // Output: left states at faces
    double* Q_right,  // Output: right states at faces
    int num_cells,
    double kappa      // MUSCL parameter: -1 (2nd upwind), 0 (Fromm), 0.5 (Quick), 1 (central)
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;
    
    for (int face = 0; face < 4; face++) {
        int neighbor_idx = cell_neighbors[cell_idx * 4 + face];
        if (neighbor_idx < 0 || neighbor_idx >= num_cells) {
            // Boundary face - use cell-centered values
            for (int var = 0; var < 4; var++) {
                int face_offset = (cell_idx * 4 + face) * 4 + var;
                Q_left[face_offset] = U_cells[cell_idx * 4 + var];
                Q_right[face_offset] = U_cells[cell_idx * 4 + var];
            }
            continue;
        }
        
        double distance = cell_distances[cell_idx * 4 + face];
        
        for (int var = 0; var < 4; var++) {
            // Get cell values
            double U_i = U_cells[cell_idx * 4 + var];
            double U_j = U_cells[neighbor_idx * 4 + var];
            
            // Get limited gradients (phi already applied)
            int grad_offset_L = (cell_idx * 4 + face) * 4 + var;
            int grad_offset_R = (neighbor_idx * 4 + face) * 4 + var;
            
            double phi_L = limited_gradients[grad_offset_L];
            double phi_R = limited_gradients[grad_offset_R];
            
            // MUSCL reconstruction
            // Left state (from current cell)
            double slope_L = (U_j - U_i) / distance;
            double Q_L = U_i + 0.5 * kappa * phi_L * slope_L * distance;
            
            // Right state (from neighbor cell)
            double slope_R = (U_i - U_j) / distance;
            double Q_R = U_j + 0.5 * kappa * phi_R * slope_R * distance;
            
            // Store reconstructed values
            int face_offset = (cell_idx * 4 + face) * 4 + var;
            Q_left[face_offset] = Q_L;
            Q_right[face_offset] = Q_R;
        }
    }
}

//=============================================================================
// WENO5 RECONSTRUCTION KERNEL
//=============================================================================

/**
 * @brief Device function for WENO5 smoothness indicator
 */
__device__ double weno5_smoothness_indicator_device(
    double u1, double u2, double u3
) {
    double term1 = u1 - 2.0 * u2 + u3;
    double term2 = u1 - 4.0 * u2 + 3.0 * u3;
    return (13.0 / 12.0) * term1 * term1 + 0.25 * term2 * term2;
}

/**
 * @brief WENO5 reconstruction kernel
 * 
 * Computes 5th-order accurate WENO reconstruction at cell faces.
 * Uses 5-point stencil with 3 candidate polynomials.
 * 
 * @param U_cells Conservative variables
 * @param cell_neighbors Neighbor indices (extended for 5-point stencil)
 * @param Q_left Output: left reconstructed states
 * @param Q_right Output: right reconstructed states
 * @param num_cells Total number of cells
 * @param epsilon WENO epsilon parameter (typically 1e-6)
 * @param p Power for nonlinear weights (typically 2)
 */
/**
 * @brief WENO5 reconstruction kernel (corrected)
 *
 * The left-biased reconstruction at face i+1/2 uses the 5-point stencil
 * {i-2, i-1, i, i+1, i+2} centered on cell i.
 * The right-biased reconstruction at face i+1/2 uses the shifted stencil
 * {i-1, i, i+1, i+2, i+3} centered on cell i+1.
 *
 * Stencil layout per face: 6 neighbor indices
 *   [0]=i-2, [1]=i-1, [2]=i+1, [3]=i+2, [4]=i+3
 * with cell_idx itself serving as i.
 */
__global__ void weno5_reconstruction_kernel(
    const double* U_cells,
    const int* cell_neighbors,
    double* Q_left,
    double* Q_right,
    int num_cells,
    double epsilon,
    int p
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    for (int face = 0; face < 4; face++) {
        int stencil_base = (cell_idx * 4 + face) * 5;
        int idx_m2 = cell_neighbors[stencil_base + 0];  // i-2
        int idx_m1 = cell_neighbors[stencil_base + 1];  // i-1
        int idx_0  = cell_idx;                           // i
        int idx_p1 = cell_neighbors[stencil_base + 2];  // i+1
        int idx_p2 = cell_neighbors[stencil_base + 3];  // i+2
        int idx_p3 = cell_neighbors[stencil_base + 4];  // i+3

        bool valid_left  = (idx_m2 >= 0 && idx_m2 < num_cells &&
                            idx_m1 >= 0 && idx_m1 < num_cells &&
                            idx_p1 >= 0 && idx_p1 < num_cells &&
                            idx_p2 >= 0 && idx_p2 < num_cells);
        bool valid_right = (idx_m1 >= 0 && idx_m1 < num_cells &&
                            idx_p1 >= 0 && idx_p1 < num_cells &&
                            idx_p2 >= 0 && idx_p2 < num_cells &&
                            idx_p3 >= 0 && idx_p3 < num_cells);

        if (!valid_left && !valid_right) {
            for (int var = 0; var < 4; var++) {
                int face_offset = (cell_idx * 4 + face) * 4 + var;
                Q_left[face_offset]  = U_cells[cell_idx * 4 + var];
                Q_right[face_offset] = (idx_p1 >= 0 && idx_p1 < num_cells)
                    ? U_cells[idx_p1 * 4 + var]
                    : U_cells[cell_idx * 4 + var];
            }
            continue;
        }

        for (int var = 0; var < 4; var++) {
            double u_0  = U_cells[idx_0 * 4 + var];

            // ============ LEFT STATE (u^- at i+1/2) ============
            // Stencil: {i-2, i-1, i, i+1, i+2}
            double Q_L;
            if (valid_left) {
                double u_m2 = U_cells[idx_m2 * 4 + var];
                double u_m1 = U_cells[idx_m1 * 4 + var];
                double u_p1 = U_cells[idx_p1 * 4 + var];
                double u_p2 = U_cells[idx_p2 * 4 + var];

                double v0_L = ( 2.0 * u_0  + 5.0 * u_p1 -       u_p2) / 6.0;
                double v1_L = (      -u_m1 + 5.0 * u_0  + 2.0 * u_p1) / 6.0;
                double v2_L = ( 2.0 * u_m2 - 7.0 * u_m1 + 11.0 * u_0) / 6.0;

                double beta0_L = weno5_smoothness_indicator_device(u_0,  u_p1, u_p2);
                double beta1_L = weno5_smoothness_indicator_device(u_m1, u_0,  u_p1);
                double beta2_L = weno5_smoothness_indicator_device(u_m2, u_m1, u_0);

                double alpha0_L = (3.0 / 10.0) / pow(epsilon + beta0_L, p);
                double alpha1_L = (6.0 / 10.0) / pow(epsilon + beta1_L, p);
                double alpha2_L = (1.0 / 10.0) / pow(epsilon + beta2_L, p);
                double sum_L = alpha0_L + alpha1_L + alpha2_L;

                Q_L = (alpha0_L * v0_L + alpha1_L * v1_L + alpha2_L * v2_L) / sum_L;
                if (!isfinite(Q_L)) Q_L = u_0;
            } else {
                Q_L = u_0;
            }

            // ============ RIGHT STATE (u^+ at i+1/2) ============
            // Stencil: {i-1, i, i+1, i+2, i+3}
            double Q_R;
            if (valid_right) {
                double u_m1 = U_cells[idx_m1 * 4 + var];
                double u_p1 = U_cells[idx_p1 * 4 + var];
                double u_p2 = U_cells[idx_p2 * 4 + var];
                double u_p3 = U_cells[idx_p3 * 4 + var];

                double v0_R = (11.0 * u_p1 - 7.0 * u_p2 + 2.0 * u_p3) / 6.0;
                double v1_R = ( 2.0 * u_0  + 5.0 * u_p1 -       u_p2) / 6.0;
                double v2_R = (      -u_m1 + 5.0 * u_0  + 2.0 * u_p1) / 6.0;

                double beta0_R = weno5_smoothness_indicator_device(u_p1, u_p2, u_p3);
                double beta1_R = weno5_smoothness_indicator_device(u_0,  u_p1, u_p2);
                double beta2_R = weno5_smoothness_indicator_device(u_m1, u_0,  u_p1);

                double alpha0_R = (1.0 / 10.0) / pow(epsilon + beta0_R, p);
                double alpha1_R = (6.0 / 10.0) / pow(epsilon + beta1_R, p);
                double alpha2_R = (3.0 / 10.0) / pow(epsilon + beta2_R, p);
                double sum_R = alpha0_R + alpha1_R + alpha2_R;

                Q_R = (alpha0_R * v0_R + alpha1_R * v1_R + alpha2_R * v2_R) / sum_R;
                if (!isfinite(Q_R)) Q_R = u_p1;
            } else {
                Q_R = (idx_p1 >= 0 && idx_p1 < num_cells)
                    ? U_cells[idx_p1 * 4 + var] : u_0;
            }

            int face_offset = (cell_idx * 4 + face) * 4 + var;
            Q_left[face_offset]  = Q_L;
            Q_right[face_offset] = Q_R;
        }
    }
}

//=============================================================================
// CHARACTERISTIC-BASED WENO5 KERNEL (More Accurate)
//=============================================================================

/**
 * @brief WENO5 reconstruction in characteristic space
 * 
 * Performs WENO reconstruction on characteristic variables for better
 * accuracy near discontinuities. Requires transformation matrices.
 */
/**
 * @brief WENO5 reconstruction in characteristic space (corrected)
 *
 * Left stencil: {i-2, i-1, i, i+1, i+2}
 * Right stencil: {i-1, i, i+1, i+2, i+3}
 * Stencil layout: [0]=i-2, [1]=i-1, [2]=i+1, [3]=i+2, [4]=i+3
 */
__global__ void weno5_characteristic_reconstruction_kernel(
    const double* U_cells,
    const int* cell_neighbors,
    const double* L_matrices,
    const double* R_matrices,
    double* Q_left,
    double* Q_right,
    int num_cells,
    double epsilon,
    int p
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    for (int face = 0; face < 4; face++) {
        int stencil_base = (cell_idx * 4 + face) * 5;
        int idx_m2 = cell_neighbors[stencil_base + 0];
        int idx_m1 = cell_neighbors[stencil_base + 1];
        int idx_0  = cell_idx;
        int idx_p1 = cell_neighbors[stencil_base + 2];
        int idx_p2 = cell_neighbors[stencil_base + 3];
        int idx_p3 = cell_neighbors[stencil_base + 4];

        bool valid_left  = (idx_m2 >= 0 && idx_m2 < num_cells &&
                            idx_m1 >= 0 && idx_m1 < num_cells &&
                            idx_p1 >= 0 && idx_p1 < num_cells &&
                            idx_p2 >= 0 && idx_p2 < num_cells);
        bool valid_right = (idx_m1 >= 0 && idx_m1 < num_cells &&
                            idx_p1 >= 0 && idx_p1 < num_cells &&
                            idx_p2 >= 0 && idx_p2 < num_cells &&
                            idx_p3 >= 0 && idx_p3 < num_cells);

        if (!valid_left && !valid_right) {
            int face_offset = (cell_idx * 4 + face) * 4;
            for (int var = 0; var < 4; var++) {
                Q_left[face_offset + var]  = U_cells[cell_idx * 4 + var];
                Q_right[face_offset + var] = (idx_p1 >= 0 && idx_p1 < num_cells)
                    ? U_cells[idx_p1 * 4 + var]
                    : U_cells[cell_idx * 4 + var];
            }
            continue;
        }

        int matrix_offset = (cell_idx * 4 + face) * 16;

        // Transform all stencil points to characteristic space
        double w_m2[4] = {0}, w_m1[4] = {0}, w_0[4] = {0};
        double w_p1[4] = {0}, w_p2[4] = {0}, w_p3[4] = {0};

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double L_ij = L_matrices[matrix_offset + i * 4 + j];
                w_0[i]  += L_ij * U_cells[idx_0 * 4 + j];
                if (valid_left) {
                    w_m2[i] += L_ij * U_cells[idx_m2 * 4 + j];
                }
                if (valid_left || valid_right) {
                    w_m1[i] += L_ij * U_cells[idx_m1 * 4 + j];
                    w_p1[i] += L_ij * U_cells[idx_p1 * 4 + j];
                    w_p2[i] += L_ij * U_cells[idx_p2 * 4 + j];
                }
                if (valid_right) {
                    w_p3[i] += L_ij * U_cells[idx_p3 * 4 + j];
                }
            }
        }

        double w_L[4], w_R[4];

        for (int var = 0; var < 4; var++) {
            // ============ LEFT STATE ============
            // Stencil: {w_m2, w_m1, w_0, w_p1, w_p2}
            if (valid_left) {
                double v0_L = ( 2.0 * w_0[var]  + 5.0 * w_p1[var] -       w_p2[var]) / 6.0;
                double v1_L = (      -w_m1[var] + 5.0 * w_0[var]  + 2.0 * w_p1[var]) / 6.0;
                double v2_L = ( 2.0 * w_m2[var] - 7.0 * w_m1[var] + 11.0 * w_0[var]) / 6.0;

                double b0_L = weno5_smoothness_indicator_device(w_0[var],  w_p1[var], w_p2[var]);
                double b1_L = weno5_smoothness_indicator_device(w_m1[var], w_0[var],  w_p1[var]);
                double b2_L = weno5_smoothness_indicator_device(w_m2[var], w_m1[var], w_0[var]);

                double a0_L = (3.0 / 10.0) / pow(epsilon + b0_L, p);
                double a1_L = (6.0 / 10.0) / pow(epsilon + b1_L, p);
                double a2_L = (1.0 / 10.0) / pow(epsilon + b2_L, p);
                double sum_L = a0_L + a1_L + a2_L;

                w_L[var] = (a0_L * v0_L + a1_L * v1_L + a2_L * v2_L) / sum_L;
            } else {
                w_L[var] = w_0[var];
            }

            // ============ RIGHT STATE ============
            // Stencil: {w_m1, w_0, w_p1, w_p2, w_p3}
            if (valid_right) {
                double v0_R = (11.0 * w_p1[var] - 7.0 * w_p2[var] + 2.0 * w_p3[var]) / 6.0;
                double v1_R = ( 2.0 * w_0[var]  + 5.0 * w_p1[var] -       w_p2[var]) / 6.0;
                double v2_R = (      -w_m1[var] + 5.0 * w_0[var]  + 2.0 * w_p1[var]) / 6.0;

                double b0_R = weno5_smoothness_indicator_device(w_p1[var], w_p2[var], w_p3[var]);
                double b1_R = weno5_smoothness_indicator_device(w_0[var],  w_p1[var], w_p2[var]);
                double b2_R = weno5_smoothness_indicator_device(w_m1[var], w_0[var],  w_p1[var]);

                double a0_R = (1.0 / 10.0) / pow(epsilon + b0_R, p);
                double a1_R = (6.0 / 10.0) / pow(epsilon + b1_R, p);
                double a2_R = (3.0 / 10.0) / pow(epsilon + b2_R, p);
                double sum_R = a0_R + a1_R + a2_R;

                w_R[var] = (a0_R * v0_R + a1_R * v1_R + a2_R * v2_R) / sum_R;
            } else {
                w_R[var] = w_p1[var];
            }
        }

        // Transform back to conservative variables
        double QL[4] = {0}, QR[4] = {0};
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double R_ij = R_matrices[matrix_offset + i * 4 + j];
                QL[i] += R_ij * w_L[j];
                QR[i] += R_ij * w_R[j];
            }
        }

        int face_offset = (cell_idx * 4 + face) * 4;
        for (int var = 0; var < 4; var++) {
            Q_left[face_offset + var]  = isfinite(QL[var]) ? QL[var] : U_cells[cell_idx * 4 + var];
            Q_right[face_offset + var] = isfinite(QR[var]) ? QR[var]
                : ((idx_p1 >= 0 && idx_p1 < num_cells) ? U_cells[idx_p1 * 4 + var] : U_cells[cell_idx * 4 + var]);
        }
    }
}
