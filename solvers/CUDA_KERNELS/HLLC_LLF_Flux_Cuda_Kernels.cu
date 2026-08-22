/**
 * @file HLLC_Flux_Cuda_Kernels.cu
 * @brief CUDA implementation of HLLC (Harten-Lax-van Leer-Contact) Riemann solver
 * @author AI Assistant
 * @date 2026-01-15
 * @version 1.0
 * 
 * @details Implements the HLLC flux scheme which resolves contact discontinuities
 * exactly while being more robust than Roe scheme (no carbuncle phenomenon).
 * 
 * Reference: Toro, E.F. "Riemann Solvers and Numerical Methods for Fluid Dynamics"
 *            3rd Edition, Chapter 10
 */

#include <cuda_runtime.h>
#include <math.h>

//=============================================================================
// DEVICE UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Estimate wave speeds for HLLC solver
 */
__device__ void estimate_wave_speeds_device(
    double rho_L, double u_L, double v_L, double P_L, double a_L,
    double rho_R, double u_R, double v_R, double P_R, double a_R,
    double nx, double ny,
    double& S_L, double& S_R, double& S_star
) {
    // Normal velocities
    double un_L = nx * u_L + ny * v_L;
    double un_R = nx * u_R + ny * v_R;
    
    // Roe-averaged sound speed
    double sqrt_rho_L = sqrt(rho_L);
    double sqrt_rho_R = sqrt(rho_R);
    double a_roe = (a_L * sqrt_rho_L + a_R * sqrt_rho_R) / (sqrt_rho_L + sqrt_rho_R);
    double u_roe = (u_L * sqrt_rho_L + u_R * sqrt_rho_R) / (sqrt_rho_L + sqrt_rho_R);
    double v_roe = (v_L * sqrt_rho_L + v_R * sqrt_rho_R) / (sqrt_rho_L + sqrt_rho_R);
    double un_roe = nx * u_roe + ny * v_roe;
    
    // Wave speed estimates (Davis)
    S_L = fmin(un_L - a_L, un_roe - a_roe);
    S_R = fmax(un_R + a_R, un_roe + a_roe);
    
    // Contact wave speed (middle wave)
    double numer = P_R - P_L + rho_L * un_L * (S_L - un_L) - rho_R * un_R * (S_R - un_R);
    double denom = rho_L * (S_L - un_L) - rho_R * (S_R - un_R);
    
    if (fabs(denom) > 1e-14) {
        S_star = numer / denom;
    } else {
        S_star = 0.5 * (un_L + un_R);
    }
}

/**
 * @brief Compute star region state for HLLC
 */
__device__ void compute_star_state_device(
    double rho, double u, double v, double E, double P,
    double nx, double ny, double S_K, double S_star,
    double U_star[4]
) {
    double un = nx * u + ny * v;
    
    double factor = rho * (S_K - un) / (S_K - S_star);
    
    U_star[0] = factor;
    U_star[1] = factor * (u + (S_star - un) * nx);
    U_star[2] = factor * (v + (S_star - un) * ny);
    U_star[3] = factor * (E / rho + (S_star - un) * (S_star + P / (rho * (S_K - un))));
}

//=============================================================================
// HLLC FLUX KERNEL
//=============================================================================

/**
 * @brief HLLC flux kernel
 * 
 * Computes HLLC Riemann flux with exact contact discontinuity resolution.
 * More robust than Roe scheme (no carbuncle, no entropy fix needed).
 */
__global__ void hllc_flux_kernel(
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
        
        // Left state
        int L_idx = cell_idx * 6;
        double rho_L = P_cells[L_idx + 0];
        double u_L = P_cells[L_idx + 1];
        double v_L = P_cells[L_idx + 2];
        double P_L = P_cells[L_idx + 3];
        double a_L = P_cells[L_idx + 5];
        
        if (rho_L <= 0.0 || P_L <= 0.0 || a_L <= 0.0) continue;
        
        int L_cons = cell_idx * 4;
        double U_L[4] = {U_cells[L_cons + 0], U_cells[L_cons + 1], 
                         U_cells[L_cons + 2], U_cells[L_cons + 3]};
        
        // Right state
        int R_idx = neighbor_idx * 6;
        double rho_R = P_cells[R_idx + 0];
        double u_R = P_cells[R_idx + 1];
        double v_R = P_cells[R_idx + 2];
        double P_R = P_cells[R_idx + 3];
        double a_R = P_cells[R_idx + 5];
        
        if (rho_R <= 0.0 || P_R <= 0.0 || a_R <= 0.0) continue;
        
        int R_cons = neighbor_idx * 4;
        double U_R[4] = {U_cells[R_cons + 0], U_cells[R_cons + 1],
                         U_cells[R_cons + 2], U_cells[R_cons + 3]};
        
        // Estimate wave speeds
        double S_L, S_R, S_star;
        estimate_wave_speeds_device(rho_L, u_L, v_L, P_L, a_L,
                                   rho_R, u_R, v_R, P_R, a_R,
                                   nx, ny, S_L, S_R, S_star);
        
        // Compute flux based on wave speeds
        double flux[4];
        
        if (S_L >= 0.0) {
            // Supersonic left
            double un_L = nx * u_L + ny * v_L;
            flux[0] = rho_L * un_L;
            flux[1] = rho_L * u_L * un_L + P_L * nx;
            flux[2] = rho_L * v_L * un_L + P_L * ny;
            flux[3] = (U_L[3] + P_L) * un_L;
            
        } else if (S_R <= 0.0) {
            // Supersonic right
            double un_R = nx * u_R + ny * v_R;
            flux[0] = rho_R * un_R;
            flux[1] = rho_R * u_R * un_R + P_R * nx;
            flux[2] = rho_R * v_R * un_R + P_R * ny;
            flux[3] = (U_R[3] + P_R) * un_R;
            
        } else if (S_star >= 0.0) {
            // Left star region
            double un_L = nx * u_L + ny * v_L;
            double F_L[4];
            F_L[0] = rho_L * un_L;
            F_L[1] = rho_L * u_L * un_L + P_L * nx;
            F_L[2] = rho_L * v_L * un_L + P_L * ny;
            F_L[3] = (U_L[3] + P_L) * un_L;
            
            double U_star_L[4];
            compute_star_state_device(rho_L, u_L, v_L, U_L[3], P_L,
                                     nx, ny, S_L, S_star, U_star_L);
            
            for (int i = 0; i < 4; i++) {
                flux[i] = F_L[i] + S_L * (U_star_L[i] - U_L[i]);
            }
            
        } else {
            // Right star region
            double un_R = nx * u_R + ny * v_R;
            double F_R[4];
            F_R[0] = rho_R * un_R;
            F_R[1] = rho_R * u_R * un_R + P_R * nx;
            F_R[2] = rho_R * v_R * un_R + P_R * ny;
            F_R[3] = (U_R[3] + P_R) * un_R;
            
            double U_star_R[4];
            compute_star_state_device(rho_R, u_R, v_R, U_R[3], P_R,
                                     nx, ny, S_R, S_star, U_star_R);
            
            for (int i = 0; i < 4; i++) {
                flux[i] = F_R[i] + S_R * (U_star_R[i] - U_R[i]);
            }
        }
        
        // Store flux
        int flux_offset = (cell_idx * 4 + face) * 4;
        for (int i = 0; i < 4; i++) {
            flux_output[flux_offset + i] = flux[i] * area;
        }
    }
}

//=============================================================================
// LLF (LOCAL LAX-FRIEDRICHS) FLUX KERNEL
//=============================================================================

/**
 * @brief LLF (Local Lax-Friedrichs) flux kernel
 * 
 * Simplest Riemann solver - very robust, good fallback.
 * F_LLF = 0.5*(F_L + F_R) - 0.5*lambda_max*(U_R - U_L)
 */
__global__ void llf_flux_kernel(
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
        
        // Left state
        int L_idx = cell_idx * 6;
        double rho_L = P_cells[L_idx + 0];
        double u_L = P_cells[L_idx + 1];
        double v_L = P_cells[L_idx + 2];
        double P_L = P_cells[L_idx + 3];
        double a_L = P_cells[L_idx + 5];
        
        if (rho_L <= 0.0 || P_L <= 0.0 || a_L <= 0.0) continue;
        
        double un_L = nx * u_L + ny * v_L;
        
        int L_cons = cell_idx * 4;
        double U_L[4] = {U_cells[L_cons + 0], U_cells[L_cons + 1], 
                         U_cells[L_cons + 2], U_cells[L_cons + 3]};
        
        // Right state
        int R_idx = neighbor_idx * 6;
        double rho_R = P_cells[R_idx + 0];
        double u_R = P_cells[R_idx + 1];
        double v_R = P_cells[R_idx + 2];
        double P_R = P_cells[R_idx + 3];
        double a_R = P_cells[R_idx + 5];
        
        if (rho_R <= 0.0 || P_R <= 0.0 || a_R <= 0.0) continue;
        
        double un_R = nx * u_R + ny * v_R;
        
        int R_cons = neighbor_idx * 4;
        double U_R[4] = {U_cells[R_cons + 0], U_cells[R_cons + 1],
                         U_cells[R_cons + 2], U_cells[R_cons + 3]};
        
        // Maximum wave speed
        double lambda_max = fmax(fabs(un_L) + a_L, fabs(un_R) + a_R);
        
        // Convective fluxes
        double F_L[4], F_R[4];
        F_L[0] = rho_L * un_L;
        F_L[1] = rho_L * u_L * un_L + P_L * nx;
        F_L[2] = rho_L * v_L * un_L + P_L * ny;
        F_L[3] = (U_L[3] + P_L) * un_L;
        
        F_R[0] = rho_R * un_R;
        F_R[1] = rho_R * u_R * un_R + P_R * nx;
        F_R[2] = rho_R * v_R * un_R + P_R * ny;
        F_R[3] = (U_R[3] + P_R) * un_R;
        
        // LLF flux: 0.5*(F_L + F_R) - 0.5*lambda_max*(U_R - U_L)
        int flux_offset = (cell_idx * 4 + face) * 4;
        for (int i = 0; i < 4; i++) {
            flux_output[flux_offset + i] = 0.5 * (F_L[i] + F_R[i] - lambda_max * (U_R[i] - U_L[i])) * area;
        }
    }
}
