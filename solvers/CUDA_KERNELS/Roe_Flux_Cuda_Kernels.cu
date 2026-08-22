/**
 * @file Roe_Flux_Cuda_Kernels.cu
 * @brief CUDA implementation of Roe's approximate Riemann solver
 * @author AI Assistant
 * @date 2026-01-15
 * @version 1.0
 * 
 * @details Implements the classical Roe flux scheme with:
 * - Roe-averaged states (density-weighted averaging)
 * - Eigenvalue/eigenvector decomposition
 * - Characteristic variable projection
 * - Entropy fix for sonic rarefactions
 * 
 * Reference: Roe, P.L. "Approximate Riemann Solvers, Parameter Vectors,
 *            and Difference Schemes" JCP 43, 357-372 (1981)
 */

#include <cuda_runtime.h>
#include <math.h>

//=============================================================================
// DEVICE UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Compute Roe-averaged quantities
 */
__device__ void compute_roe_averages_device(
    double rho_L, double u_L, double v_L, double H_L,
    double rho_R, double u_R, double v_R, double H_R,
    double gamma,
    double& roe_u, double& roe_v, double& roe_H, double& roe_a
) {
    double sqrt_rho_L = sqrt(rho_L);
    double sqrt_rho_R = sqrt(rho_R);
    double term = 1.0 / (sqrt_rho_L + sqrt_rho_R);
    
    // Density-weighted averaging
    roe_u = term * (u_L * sqrt_rho_L + u_R * sqrt_rho_R);
    roe_v = term * (v_L * sqrt_rho_L + v_R * sqrt_rho_R);
    roe_H = term * (H_L * sqrt_rho_L + H_R * sqrt_rho_R);
    
    // Sound speed from enthalpy
    double roe_Vmag = 0.5 * (roe_u * roe_u + roe_v * roe_v);
    roe_a = sqrt(fmax((gamma - 1.0) * (roe_H - roe_Vmag), 1e-14));
}

/**
 * @brief Apply Harten's entropy fix
 */
__device__ double apply_entropy_fix_device(double lambda, double a) {
    double delta = 0.1 * a;
    if (fabs(lambda) < delta) {
        return 0.5 * (lambda * lambda / delta + delta);
    }
    return fabs(lambda);
}

//=============================================================================
// ROE FLUX KERNEL - FIRST ORDER
//=============================================================================

/**
 * @brief First-order Roe flux kernel
 * 
 * Computes Roe's approximate Riemann flux for each cell face.
 * 
 * @param U_cells Conservative variables [rho, rho*u, rho*v, E] per cell
 * @param P_cells Primitive variables [rho, u, v, P, T, a] per cell
 * @param face_neighbors Neighbor cell indices for each face
 * @param face_normals Normal vectors [nx, ny] for each face (8 per cell)
 * @param face_areas Face areas
 * @param flux_output Output flux array [4 components × num_faces]
 * @param num_cells Total number of cells
 * @param gamma Ratio of specific heats
 */
__global__ void roe_flux_kernel(
    const double* U_cells,
    const double* P_cells,
    const int* face_neighbors,
    const double* face_normals,
    const double* face_areas,
    double* flux_output,
    int num_cells,
    double gamma
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;
    
    // Process each face
    for (int face = 0; face < 4; face++) {
        int neighbor_idx = face_neighbors[cell_idx * 4 + face];
        
        // Skip invalid neighbors
        if (neighbor_idx < 0 || neighbor_idx >= num_cells) continue;
        
        // Extract face geometry
        int normal_offset = cell_idx * 8 + face * 2;
        double nx = face_normals[normal_offset + 0];
        double ny = face_normals[normal_offset + 1];
        double area = face_areas[cell_idx * 4 + face];
        
        // Normalize face normal
        double normal_mag = sqrt(nx * nx + ny * ny);
        if (normal_mag < 1e-14) continue;
        nx /= normal_mag;
        ny /= normal_mag;
        
        // Extract left state (current cell)
        int L_offset = cell_idx * 6;
        double rho_L = P_cells[L_offset + 0];
        double u_L = P_cells[L_offset + 1];
        double v_L = P_cells[L_offset + 2];
        double P_L = P_cells[L_offset + 3];
        double a_L = P_cells[L_offset + 5];
        
        // Validate left state
        if (rho_L <= 0.0 || P_L <= 0.0 || a_L <= 0.0) continue;
        
        double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
        double H_L = (a_L * a_L) / (gamma - 1.0) + Vmag_L;
        
        // Extract right state (neighbor cell)
        int R_offset = neighbor_idx * 6;
        double rho_R = P_cells[R_offset + 0];
        double u_R = P_cells[R_offset + 1];
        double v_R = P_cells[R_offset + 2];
        double P_R = P_cells[R_offset + 3];
        double a_R = P_cells[R_offset + 5];
        
        // Validate right state
        if (rho_R <= 0.0 || P_R <= 0.0 || a_R <= 0.0) continue;
        
        double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);
        double H_R = (a_R * a_R) / (gamma - 1.0) + Vmag_R;
        
        // Compute Roe averages
        double roe_u, roe_v, roe_H, roe_a;
        compute_roe_averages_device(rho_L, u_L, v_L, H_L,
                                   rho_R, u_R, v_R, H_R,
                                   gamma, roe_u, roe_v, roe_H, roe_a);
        
        double sqrt_rho_L = sqrt(rho_L);
        double sqrt_rho_R = sqrt(rho_R);
        double roe_rho = sqrt(rho_L * rho_R);
        
        // Normal and tangential velocities
        double Un = nx * roe_u + ny * roe_v;
        double Ut = -ny * roe_u + nx * roe_v;
        
        // State differences
        double du = u_R - u_L;
        double dv = v_R - v_L;
        double dP = P_R - P_L;
        double drho = rho_R - rho_L;
        double dUn = nx * du + ny * dv;
        double dUt = -ny * du + nx * dv;
        
        // Wave strengths (characteristic variable coefficients)
        double a2_inv = 1.0 / (roe_a * roe_a);
        double alpha_1 = 0.5 * a2_inv * (dP - roe_rho * roe_a * dUn);
        double alpha_2 = drho - a2_inv * dP;
        double alpha_3 = roe_rho * dUt;
        double alpha_4 = 0.5 * a2_inv * (dP + roe_rho * roe_a * dUn);
        
        // Eigenvalues with entropy fix
        double lambda_1 = apply_entropy_fix_device(Un - roe_a, roe_a);
        double lambda_2 = fabs(Un);
        double lambda_3 = fabs(Un);
        double lambda_4 = apply_entropy_fix_device(Un + roe_a, roe_a);
        
        // Right eigenvectors
        // R1: lambda = u - a
        double R11 = 1.0;
        double R12 = roe_u - nx * roe_a;
        double R13 = roe_v - ny * roe_a;
        double R14 = roe_H - roe_a * Un;
        
        // R2: lambda = u (entropy)
        double R21 = 1.0;
        double R22 = roe_u;
        double R23 = roe_v;
        double R24 = 0.5 * (roe_u * roe_u + roe_v * roe_v);
        
        // R3: lambda = u (shear)
        double R31 = 0.0;
        double R32 = -ny;
        double R33 = nx;
        double R34 = Ut;
        
        // R4: lambda = u + a
        double R41 = 1.0;
        double R42 = roe_u + nx * roe_a;
        double R43 = roe_v + ny * roe_a;
        double R44 = roe_H + roe_a * Un;
        
        // Compute dissipative flux: D = 0.5 * sum(|lambda_k| * alpha_k * R_k)
        double D1 = 0.5 * (lambda_1 * alpha_1 * R11 + 
                          lambda_2 * alpha_2 * R21 + 
                          lambda_3 * alpha_3 * R31 + 
                          lambda_4 * alpha_4 * R41);
        
        double D2 = 0.5 * (lambda_1 * alpha_1 * R12 + 
                          lambda_2 * alpha_2 * R22 + 
                          lambda_3 * alpha_3 * R32 + 
                          lambda_4 * alpha_4 * R42);
        
        double D3 = 0.5 * (lambda_1 * alpha_1 * R13 + 
                          lambda_2 * alpha_2 * R23 + 
                          lambda_3 * alpha_3 * R33 + 
                          lambda_4 * alpha_4 * R43);
        
        double D4 = 0.5 * (lambda_1 * alpha_1 * R14 + 
                          lambda_2 * alpha_2 * R24 + 
                          lambda_3 * alpha_3 * R34 + 
                          lambda_4 * alpha_4 * R44);
        
        // Compute convective fluxes
        double Un_L = nx * u_L + ny * v_L;
        double Un_R = nx * u_R + ny * v_R;
        
        double F_L[4], F_R[4];
        F_L[0] = rho_L * Un_L;
        F_L[1] = rho_L * u_L * Un_L + P_L * nx;
        F_L[2] = rho_L * v_L * Un_L + P_L * ny;
        F_L[3] = (U_cells[cell_idx * 4 + 3] + P_L) * Un_L;
        
        F_R[0] = rho_R * Un_R;
        F_R[1] = rho_R * u_R * Un_R + P_R * nx;
        F_R[2] = rho_R * v_R * Un_R + P_R * ny;
        F_R[3] = (U_cells[neighbor_idx * 4 + 3] + P_R) * Un_R;
        
        // Final Roe flux: F_roe = 0.5*(F_L + F_R) - D
        int flux_offset = (cell_idx * 4 + face) * 4;
        flux_output[flux_offset + 0] = 0.5 * (F_L[0] + F_R[0] - D1) * area;
        flux_output[flux_offset + 1] = 0.5 * (F_L[1] + F_R[1] - D2) * area;
        flux_output[flux_offset + 2] = 0.5 * (F_L[2] + F_R[2] - D3) * area;
        flux_output[flux_offset + 3] = 0.5 * (F_L[3] + F_R[3] - D4) * area;
    }
}

//=============================================================================
// ROE FLUX KERNEL - SECOND ORDER (with reconstruction)
//=============================================================================

/**
 * @brief Second-order Roe flux kernel with MUSCL reconstruction
 * 
 * @param U_cells Conservative variables
 * @param P_cells Primitive variables
 * @param face_neighbors Neighbor indices
 * @param face_normals Face normal vectors
 * @param face_areas Face areas
 * @param limited_gradients Pre-computed limited gradients from limiter kernel
 * @param cell_centers Cell centroid coordinates
 * @param face_centers Face centroid coordinates
 * @param flux_output Output flux array
 * @param num_cells Total number of cells
 * @param gamma Ratio of specific heats
 */
__global__ void roe_flux_2nd_order_kernel(
    const double* U_cells,
    const double* P_cells,
    const int* face_neighbors,
    const double* face_normals,
    const double* face_areas,
    const double* limited_gradients,
    const double* cell_centers,
    const double* face_centers,
    double* flux_output,
    int num_cells,
    double gamma
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;
    
    // Get cell center
    double xc = cell_centers[cell_idx * 2 + 0];
    double yc = cell_centers[cell_idx * 2 + 1];
    
    // Process each face
    for (int face = 0; face < 4; face++) {
        int neighbor_idx = face_neighbors[cell_idx * 4 + face];
        if (neighbor_idx < 0 || neighbor_idx >= num_cells) continue;
        
        // Face geometry
        int normal_offset = cell_idx * 8 + face * 2;
        double nx = face_normals[normal_offset + 0];
        double ny = face_normals[normal_offset + 1];
        double area = face_areas[cell_idx * 4 + face];
        
        double normal_mag = sqrt(nx * nx + ny * ny);
        if (normal_mag < 1e-14) continue;
        nx /= normal_mag;
        ny /= normal_mag;
        
        // Get face center
        int face_center_offset = (cell_idx * 4 + face) * 2;
        double xf = face_centers[face_center_offset + 0];
        double yf = face_centers[face_center_offset + 1];
        
        // Distance from cell center to face center
        double dx_L = xf - xc;
        double dy_L = yf - yc;
        
        // Get neighbor cell center
        double xc_R = cell_centers[neighbor_idx * 2 + 0];
        double yc_R = cell_centers[neighbor_idx * 2 + 1];
        double dx_R = xf - xc_R;
        double dy_R = yf - yc_R;
        
        // Reconstruct left state using gradients
        int grad_offset_L = (cell_idx * 4 + face) * 4;
        double U_L[4];
        for (int var = 0; var < 4; var++) {
            double gradient = limited_gradients[grad_offset_L + var];
            U_L[var] = U_cells[cell_idx * 4 + var] + gradient * dx_L;
        }
        
        // Reconstruct right state
        int grad_offset_R = (neighbor_idx * 4 + face) * 4;
        double U_R[4];
        for (int var = 0; var < 4; var++) {
            double gradient = limited_gradients[grad_offset_R + var];
            U_R[var] = U_cells[neighbor_idx * 4 + var] - gradient * dx_R;
        }
        
        // Convert to primitive variables
        double rho_L = fmax(U_L[0], 1e-14);
        double u_L = U_L[1] / rho_L;
        double v_L = U_L[2] / rho_L;
        double E_L = U_L[3];
        double P_L = fmax((gamma - 1.0) * (E_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L)), 1e-14);
        double a_L = sqrt(gamma * P_L / rho_L);
        double H_L = (E_L + P_L) / rho_L;
        
        double rho_R = fmax(U_R[0], 1e-14);
        double u_R = U_R[1] / rho_R;
        double v_R = U_R[2] / rho_R;
        double E_R = U_R[3];
        double P_R = fmax((gamma - 1.0) * (E_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R)), 1e-14);
        double a_R = sqrt(gamma * P_R / rho_R);
        double H_R = (E_R + P_R) / rho_R;
        
        // Compute Roe averages
        double roe_u, roe_v, roe_H, roe_a;
        compute_roe_averages_device(rho_L, u_L, v_L, H_L,
                                   rho_R, u_R, v_R, H_R,
                                   gamma, roe_u, roe_v, roe_H, roe_a);
        
        double roe_rho = sqrt(rho_L * rho_R);
        double Un = nx * roe_u + ny * roe_v;
        double Ut = -ny * roe_u + nx * roe_v;
        
        // State differences
        double du = u_R - u_L;
        double dv = v_R - v_L;
        double dP = P_R - P_L;
        double drho = rho_R - rho_L;
        double dUn = nx * du + ny * dv;
        double dUt = -ny * du + nx * dv;
        
        // Wave strengths
        double a2_inv = 1.0 / (roe_a * roe_a);
        double alpha_1 = 0.5 * a2_inv * (dP - roe_rho * roe_a * dUn);
        double alpha_2 = drho - a2_inv * dP;
        double alpha_3 = roe_rho * dUt;
        double alpha_4 = 0.5 * a2_inv * (dP + roe_rho * roe_a * dUn);
        
        // Eigenvalues with entropy fix
        double lambda_1 = apply_entropy_fix_device(Un - roe_a, roe_a);
        double lambda_2 = fabs(Un);
        double lambda_3 = fabs(Un);
        double lambda_4 = apply_entropy_fix_device(Un + roe_a, roe_a);
        
        // Right eigenvectors (same as first order)
        double R11 = 1.0, R12 = roe_u - nx * roe_a, R13 = roe_v - ny * roe_a, R14 = roe_H - roe_a * Un;
        double R21 = 1.0, R22 = roe_u, R23 = roe_v, R24 = 0.5 * (roe_u * roe_u + roe_v * roe_v);
        double R31 = 0.0, R32 = -ny, R33 = nx, R34 = Ut;
        double R41 = 1.0, R42 = roe_u + nx * roe_a, R43 = roe_v + ny * roe_a, R44 = roe_H + roe_a * Un;
        
        // Dissipative flux
        double D1 = 0.5 * (lambda_1 * alpha_1 * R11 + lambda_2 * alpha_2 * R21 + lambda_3 * alpha_3 * R31 + lambda_4 * alpha_4 * R41);
        double D2 = 0.5 * (lambda_1 * alpha_1 * R12 + lambda_2 * alpha_2 * R22 + lambda_3 * alpha_3 * R32 + lambda_4 * alpha_4 * R42);
        double D3 = 0.5 * (lambda_1 * alpha_1 * R13 + lambda_2 * alpha_2 * R23 + lambda_3 * alpha_3 * R33 + lambda_4 * alpha_4 * R43);
        double D4 = 0.5 * (lambda_1 * alpha_1 * R14 + lambda_2 * alpha_2 * R24 + lambda_3 * alpha_3 * R34 + lambda_4 * alpha_4 * R44);
        
        // Convective fluxes with reconstructed states
        double Un_L = nx * u_L + ny * v_L;
        double Un_R = nx * u_R + ny * v_R;
        
        double F_L[4], F_R[4];
        F_L[0] = rho_L * Un_L;
        F_L[1] = rho_L * u_L * Un_L + P_L * nx;
        F_L[2] = rho_L * v_L * Un_L + P_L * ny;
        F_L[3] = (E_L + P_L) * Un_L;
        
        F_R[0] = rho_R * Un_R;
        F_R[1] = rho_R * u_R * Un_R + P_R * nx;
        F_R[2] = rho_R * v_R * Un_R + P_R * ny;
        F_R[3] = (E_R + P_R) * Un_R;
        
        // Final flux
        int flux_offset = (cell_idx * 4 + face) * 4;
        flux_output[flux_offset + 0] = 0.5 * (F_L[0] + F_R[0] - D1) * area;
        flux_output[flux_offset + 1] = 0.5 * (F_L[1] + F_R[1] - D2) * area;
        flux_output[flux_offset + 2] = 0.5 * (F_L[2] + F_R[2] - D3) * area;
        flux_output[flux_offset + 3] = 0.5 * (F_L[3] + F_R[3] - D4) * area;
    }
}
