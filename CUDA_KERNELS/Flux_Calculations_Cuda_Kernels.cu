// Flux_Calculations_Cuda_Kernels.cu  
// Priority 1: Complete Flux calculation CUDA kernels (AUSM, Van Leer, Central)
#include <cuda_runtime.h>
#include <math.h>

// Constants from the original code
#define gamma 1.4
#define R 287.5
#define gamma_R 402.5

// Device functions for flux calculations
__device__ double evaluate_mach_plus(double M) {
    if (M >= 1.0) return M;
    else if (M > -1.0 && M < 1.0) return 0.25 * (M + 1.0) * (M + 1.0);
    else return 0.0;
}

__device__ double evaluate_mach_minus(double M) {
    if (M >= 1.0) return 0.0;
    else if (M > -1.0 && M < 1.0) return -0.25 * (M - 1.0) * (M - 1.0);
    else return M;
}

__device__ double evaluate_pressure_plus(double P, double M) {
    if (M >= 1.0) return P;
    else if (M > -1.0 && M < 1.0) return (P * 0.25) * (M + 1.0) * (M + 1.0) * (2.0 - M);
    else return 0.0;
}

__device__ double evaluate_pressure_minus(double P, double M) {
    if (M >= 1.0) return 0.0;
    else if (M > -1.0 && M < 1.0) return (P * 0.25) * (M - 1.0) * (M - 1.0) * (2.0 + M);
    else return P;
}

// Complete AUSM flux kernel
__global__ void ausm_flux_kernel(
    const double* U_left,      // [num_faces * 5] - Left state conservative variables
    const double* U_right,     // [num_faces * 5] - Right state conservative variables  
    const double* face_normals, // [num_faces * 3] - Face normal vectors
    const double* face_areas,   // [num_faces] - Face areas
    double* fluxes,            // [num_faces * 5] - Output fluxes
    int num_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_faces) return;

    // Left state primitive variables
    double rho_L = U_left[idx * 5 + 0];
    double u_L = U_left[idx * 5 + 1] / rho_L;
    double v_L = U_left[idx * 5 + 2] / rho_L;
    double w_L = U_left[idx * 5 + 3] / rho_L;
    double E_L = U_left[idx * 5 + 4];
    double P_L = (gamma - 1.0) * (E_L - 0.5 * rho_L * (u_L*u_L + v_L*v_L + w_L*w_L));
    double T_L = P_L / (rho_L * R);
    double a_L = sqrt(gamma * R * T_L);

    // Right state primitive variables  
    double rho_R = U_right[idx * 5 + 0];
    double u_R = U_right[idx * 5 + 1] / rho_R;
    double v_R = U_right[idx * 5 + 2] / rho_R;
    double w_R = U_right[idx * 5 + 3] / rho_R;
    double E_R = U_right[idx * 5 + 4];
    double P_R = (gamma - 1.0) * (E_R - 0.5 * rho_R * (u_R*u_R + v_R*v_R + w_R*w_R));
    double T_R = P_R / (rho_R * R);
    double a_R = sqrt(gamma * R * T_R);

    // Face normal components
    double nx = face_normals[idx * 3 + 0];
    double ny = face_normals[idx * 3 + 1]; 
    double nz = face_normals[idx * 3 + 2];
    double area = face_areas[idx];

    // Normal velocities
    double V_L = u_L * nx + v_L * ny + w_L * nz;
    double V_R = u_R * nx + v_R * ny + w_R * nz;

    // Mach numbers
    double M_L = V_L / a_L;
    double M_R = V_R / a_R;

    // AUSM split Mach numbers and pressures
    double M_L_plus = evaluate_mach_plus(M_L);
    double M_R_minus = evaluate_mach_minus(M_R);
    double P_L_plus = evaluate_pressure_plus(P_L, M_L);
    double P_R_minus = evaluate_pressure_minus(P_R, M_R);

    // Interface Mach number and pressure
    double M_interface = M_L_plus + M_R_minus;
    double P_interface = P_L_plus + P_R_minus;

    // AUSM fluxes
    double H_L = (E_L + P_L) / rho_L;
    double H_R = (E_R + P_R) / rho_R;

    if (M_interface > 0) {
        // Use left state
        fluxes[idx * 5 + 0] = area * rho_L * a_L * M_interface;
        fluxes[idx * 5 + 1] = area * (rho_L * a_L * M_interface * u_L + P_interface * nx);
        fluxes[idx * 5 + 2] = area * (rho_L * a_L * M_interface * v_L + P_interface * ny);
        fluxes[idx * 5 + 3] = area * (rho_L * a_L * M_interface * w_L + P_interface * nz);
        fluxes[idx * 5 + 4] = area * rho_L * a_L * M_interface * H_L;
    } else {
        // Use right state
        fluxes[idx * 5 + 0] = area * rho_R * a_R * M_interface;
        fluxes[idx * 5 + 1] = area * (rho_R * a_R * M_interface * u_R + P_interface * nx);
        fluxes[idx * 5 + 2] = area * (rho_R * a_R * M_interface * v_R + P_interface * ny);
        fluxes[idx * 5 + 3] = area * (rho_R * a_R * M_interface * w_R + P_interface * nz);
        fluxes[idx * 5 + 4] = area * rho_R * a_R * M_interface * H_R;
    }
}

// Complete Van Leer flux splitting kernel
__global__ void vanleer_flux_kernel(
    const double* U_left,      // [num_faces * 5] - Left state
    const double* U_right,     // [num_faces * 5] - Right state
    const double* face_normals, // [num_faces * 3] - Face normals
    const double* face_areas,   // [num_faces] - Face areas
    double* fluxes,            // [num_faces * 5] - Output fluxes
    int num_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_faces) return;

    // Left state primitive variables
    double rho_L = U_left[idx * 5 + 0];
    double u_L = U_left[idx * 5 + 1] / rho_L;
    double v_L = U_left[idx * 5 + 2] / rho_L;
    double w_L = U_left[idx * 5 + 3] / rho_L;
    double E_L = U_left[idx * 5 + 4];
    double P_L = (gamma - 1.0) * (E_L - 0.5 * rho_L * (u_L*u_L + v_L*v_L + w_L*w_L));
    double a_L = sqrt(gamma * P_L / rho_L);

    // Right state primitive variables
    double rho_R = U_right[idx * 5 + 0];
    double u_R = U_right[idx * 5 + 1] / rho_R;
    double v_R = U_right[idx * 5 + 2] / rho_R;
    double w_R = U_right[idx * 5 + 3] / rho_R;
    double E_R = U_right[idx * 5 + 4];
    double P_R = (gamma - 1.0) * (E_R - 0.5 * rho_R * (u_R*u_R + v_R*v_R + w_R*w_R));
    double a_R = sqrt(gamma * P_R / rho_R);

    // Face normal and area
    double nx = face_normals[idx * 3 + 0];
    double ny = face_normals[idx * 3 + 1];
    double nz = face_normals[idx * 3 + 2];
    double area = face_areas[idx];

    // Normal velocities and Mach numbers
    double V_L = u_L * nx + v_L * ny + w_L * nz;
    double V_R = u_R * nx + v_R * ny + w_R * nz;
    double M_L = V_L / a_L;
    double M_R = V_R / a_R;

    // Van Leer flux splitting
    double flux_plus[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double flux_minus[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    // Left state (plus) flux
    if (M_L >= 1.0) {
        flux_plus[0] = rho_L * V_L;
        flux_plus[1] = rho_L * V_L * u_L + P_L * nx;
        flux_plus[2] = rho_L * V_L * v_L + P_L * ny;
        flux_plus[3] = rho_L * V_L * w_L + P_L * nz;
        flux_plus[4] = rho_L * V_L * (E_L + P_L) / rho_L;
    } else if (M_L > -1.0) {
        double rho_plus = rho_L * a_L * 0.25 * (M_L + 1.0) * (M_L + 1.0);
        double temp = (gamma - 1.0) * V_L + 2.0 * a_L;
        flux_plus[0] = rho_plus;
        flux_plus[1] = rho_plus * (nx * (-V_L + 2.0*a_L)/gamma + u_L);
        flux_plus[2] = rho_plus * (ny * (-V_L + 2.0*a_L)/gamma + v_L);
        flux_plus[3] = rho_plus * (nz * (-V_L + 2.0*a_L)/gamma + w_L);
        flux_plus[4] = 0.5 * rho_plus * (temp*temp/(gamma*gamma-1) + (u_L*u_L + v_L*v_L + w_L*w_L - V_L*V_L));
    }

    // Right state (minus) flux
    if (M_R <= -1.0) {
        flux_minus[0] = rho_R * V_R;
        flux_minus[1] = rho_R * V_R * u_R + P_R * nx;
        flux_minus[2] = rho_R * V_R * v_R + P_R * ny;
        flux_minus[3] = rho_R * V_R * w_R + P_R * nz;
        flux_minus[4] = rho_R * V_R * (E_R + P_R) / rho_R;
    } else if (M_R < 1.0) {
        double rho_minus = -rho_R * a_R * 0.25 * (M_R - 1.0) * (M_R - 1.0);
        double temp = (gamma - 1.0) * V_R - 2.0 * a_R;
        flux_minus[0] = rho_minus;
        flux_minus[1] = rho_minus * (nx * (-V_R - 2.0*a_R)/gamma + u_R);
        flux_minus[2] = rho_minus * (ny * (-V_R - 2.0*a_R)/gamma + v_R);
        flux_minus[3] = rho_minus * (nz * (-V_R - 2.0*a_R)/gamma + w_R);
        flux_minus[4] = 0.5 * rho_minus * (temp*temp/(gamma*gamma-1) + (u_R*u_R + v_R*v_R + w_R*w_R - V_R*V_R));
    }

    // Combine fluxes
    for (int i = 0; i < 5; i++) {
        fluxes[idx * 5 + i] = area * (flux_plus[i] + flux_minus[i]);
    }
}

// Central difference flux kernel
__global__ void central_flux_kernel(
    const double* U_left,
    const double* U_right,
    const double* face_normals,
    const double* face_areas,
    double* fluxes,
    int num_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_faces) return;

    // Average states
    double rho_avg = 0.5 * (U_left[idx * 5 + 0] + U_right[idx * 5 + 0]);
    double rhou_avg = 0.5 * (U_left[idx * 5 + 1] + U_right[idx * 5 + 1]);
    double rhov_avg = 0.5 * (U_left[idx * 5 + 2] + U_right[idx * 5 + 2]);
    double rhow_avg = 0.5 * (U_left[idx * 5 + 3] + U_right[idx * 5 + 3]);
    double E_avg = 0.5 * (U_left[idx * 5 + 4] + U_right[idx * 5 + 4]);

    double u_avg = rhou_avg / rho_avg;
    double v_avg = rhov_avg / rho_avg;
    double w_avg = rhow_avg / rho_avg;
    double P_avg = (gamma - 1.0) * (E_avg - 0.5 * rho_avg * (u_avg*u_avg + v_avg*v_avg + w_avg*w_avg));

    // Face properties
    double nx = face_normals[idx * 3 + 0];
    double ny = face_normals[idx * 3 + 1];
    double nz = face_normals[idx * 3 + 2];
    double area = face_areas[idx];
    double V_n = u_avg * nx + v_avg * ny + w_avg * nz;

    // Central fluxes
    fluxes[idx * 5 + 0] = area * rho_avg * V_n;
    fluxes[idx * 5 + 1] = area * (rho_avg * V_n * u_avg + P_avg * nx);
    fluxes[idx * 5 + 2] = area * (rho_avg * V_n * v_avg + P_avg * ny);
    fluxes[idx * 5 + 3] = area * (rho_avg * V_n * w_avg + P_avg * nz);
    fluxes[idx * 5 + 4] = area * rho_avg * V_n * (E_avg + P_avg) / rho_avg;
}
