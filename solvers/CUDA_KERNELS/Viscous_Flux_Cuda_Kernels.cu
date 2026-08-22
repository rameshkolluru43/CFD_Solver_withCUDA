// Viscous_Flux_Cuda_Kernels.cu
// Priority 2: Complete viscous flux and stress tensor CUDA kernels
#include <cuda_runtime.h>
#include <math.h>

// Constants from original code
#define gamma 1.4
#define R 287.5
#define V_S_T 110.4
#define V_C1 1.45793265452e-06
#define T_S_T 194.4  
#define T_C1 0.0025001353447
#define Pr 0.72

// Device function for Sutherland's law viscosity
__device__ double calculate_viscosity(double T) {
    return V_C1 * pow(T, 1.5) / (T + V_S_T);
}

// Device function for thermal conductivity
__device__ double calculate_thermal_conductivity(double T, double mu) {
    return T_C1 * pow(T, 1.5) / (T + T_S_T);
    // Alternative: return mu * gamma * R / ((gamma - 1.0) * Pr);
}

// Complete viscous stress tensor calculation kernel
__global__ void calculate_viscous_stress_kernel(
    const double* velocity_gradients, // [num_cells * 9] - Velocity gradient tensor
    const double* temperature_gradients, // [num_cells * 3] - Temperature gradients
    const double* primitive_vars,      // [num_cells * 5] - Primitive variables (rho,u,v,w,P)
    const double* face_normals,        // [num_faces * 3] - Face normal vectors
    const double* face_areas,          // [num_faces] - Face areas
    const int* face_cells,             // [num_faces * 2] - Left and right cell indices
    double* viscous_fluxes,           // [num_faces * 5] - Output viscous fluxes
    int num_faces
) {
    int face_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (face_idx >= num_faces) return;

    int left_cell = face_cells[face_idx * 2 + 0];
    int right_cell = face_cells[face_idx * 2 + 1];
    
    // Face normal and area
    double nx = face_normals[face_idx * 3 + 0];
    double ny = face_normals[face_idx * 3 + 1];
    double nz = face_normals[face_idx * 3 + 2];
    double area = face_areas[face_idx];

    // Average velocity gradients at face
    double dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    double dTdx, dTdy, dTdz;
    
    if (right_cell >= 0) {
        // Internal face - average gradients
        dudx = 0.5 * (velocity_gradients[left_cell * 9 + 0] + velocity_gradients[right_cell * 9 + 0]);
        dudy = 0.5 * (velocity_gradients[left_cell * 9 + 1] + velocity_gradients[right_cell * 9 + 1]);
        dudz = 0.5 * (velocity_gradients[left_cell * 9 + 2] + velocity_gradients[right_cell * 9 + 2]);
        dvdx = 0.5 * (velocity_gradients[left_cell * 9 + 3] + velocity_gradients[right_cell * 9 + 3]);
        dvdy = 0.5 * (velocity_gradients[left_cell * 9 + 4] + velocity_gradients[right_cell * 9 + 4]);
        dvdz = 0.5 * (velocity_gradients[left_cell * 9 + 5] + velocity_gradients[right_cell * 9 + 5]);
        dwdx = 0.5 * (velocity_gradients[left_cell * 9 + 6] + velocity_gradients[right_cell * 9 + 6]);
        dwdy = 0.5 * (velocity_gradients[left_cell * 9 + 7] + velocity_gradients[right_cell * 9 + 7]);
        dwdz = 0.5 * (velocity_gradients[left_cell * 9 + 8] + velocity_gradients[right_cell * 9 + 8]);
        
        dTdx = 0.5 * (temperature_gradients[left_cell * 3 + 0] + temperature_gradients[right_cell * 3 + 0]);
        dTdy = 0.5 * (temperature_gradients[left_cell * 3 + 1] + temperature_gradients[right_cell * 3 + 1]);
        dTdz = 0.5 * (temperature_gradients[left_cell * 3 + 2] + temperature_gradients[right_cell * 3 + 2]);
    } else {
        // Boundary face - use left cell gradients
        dudx = velocity_gradients[left_cell * 9 + 0];
        dudy = velocity_gradients[left_cell * 9 + 1];
        dudz = velocity_gradients[left_cell * 9 + 2];
        dvdx = velocity_gradients[left_cell * 9 + 3];
        dvdy = velocity_gradients[left_cell * 9 + 4];
        dvdz = velocity_gradients[left_cell * 9 + 5];
        dwdx = velocity_gradients[left_cell * 9 + 6];
        dwdy = velocity_gradients[left_cell * 9 + 7];
        dwdz = velocity_gradients[left_cell * 9 + 8];
        
        dTdx = temperature_gradients[left_cell * 3 + 0];
        dTdy = temperature_gradients[left_cell * 3 + 1];
        dTdz = temperature_gradients[left_cell * 3 + 2];
    }

    // Get average temperature and calculate transport properties
    double T_avg = 0.5 * (primitive_vars[left_cell * 5 + 4] / (primitive_vars[left_cell * 5 + 0] * R));
    if (right_cell >= 0) {
        T_avg += 0.5 * (primitive_vars[right_cell * 5 + 4] / (primitive_vars[right_cell * 5 + 0] * R));
    }
    
    double mu = calculate_viscosity(T_avg);
    double k = calculate_thermal_conductivity(T_avg, mu);
    double lambda = -2.0/3.0 * mu; // Bulk viscosity

    // Calculate velocity divergence
    double div_v = dudx + dvdy + dwdz;

    // Calculate stress tensor components
    double Txx = mu * (2.0 * dudx + lambda * div_v);
    double Tyy = mu * (2.0 * dvdy + lambda * div_v);
    double Tzz = mu * (2.0 * dwdz + lambda * div_v);
    double Txy = mu * (dudy + dvdx);
    double Txz = mu * (dudz + dwdx);
    double Tyz = mu * (dvdz + dwdy);

    // Calculate heat flux components
    double qx = -k * dTdx;
    double qy = -k * dTdy;
    double qz = -k * dTdz;

    // Calculate viscous flux components in face normal direction
    double tau_nn = Txx * nx * nx + Tyy * ny * ny + Tzz * nz * nz +
                    2.0 * (Txy * nx * ny + Txz * nx * nz + Tyz * ny * nz);
    double tau_nx = Txx * nx + Txy * ny + Txz * nz;
    double tau_ny = Txy * nx + Tyy * ny + Tyz * nz;
    double tau_nz = Txz * nx + Tyz * ny + Tzz * nz;
    double q_n = qx * nx + qy * ny + qz * nz;

    // Get average velocities
    double u_avg = 0.5 * (primitive_vars[left_cell * 5 + 1]);
    double v_avg = 0.5 * (primitive_vars[left_cell * 5 + 2]);
    double w_avg = 0.5 * (primitive_vars[left_cell * 5 + 3]);
    
    if (right_cell >= 0) {
        u_avg += 0.5 * primitive_vars[right_cell * 5 + 1];
        v_avg += 0.5 * primitive_vars[right_cell * 5 + 2];
        w_avg += 0.5 * primitive_vars[right_cell * 5 + 3];
    }

    // Assemble viscous fluxes
    viscous_fluxes[face_idx * 5 + 0] = 0.0; // No mass flux from viscous terms
    viscous_fluxes[face_idx * 5 + 1] = area * tau_nx;
    viscous_fluxes[face_idx * 5 + 2] = area * tau_ny;
    viscous_fluxes[face_idx * 5 + 3] = area * tau_nz;
    viscous_fluxes[face_idx * 5 + 4] = area * (tau_nx * u_avg + tau_ny * v_avg + tau_nz * w_avg + q_n);
}

// Kernel for calculating cell-wise viscous flux by summing face contributions
__global__ void accumulate_cell_viscous_flux_kernel(
    const double* face_viscous_fluxes, // [num_faces * 5] - Face viscous fluxes
    const int* cell_faces,             // [num_cells * max_faces] - Face indices for each cell
    const int* num_faces_per_cell,     // [num_cells] - Number of faces per cell
    const int* face_orientations,      // [num_cells * max_faces] - Face orientation (+1 or -1)
    double* cell_viscous_fluxes,       // [num_cells * 5] - Output cell viscous fluxes
    int num_cells,
    int max_faces
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    int num_faces = num_faces_per_cell[cell_idx];
    
    // Initialize cell viscous flux
    for (int var = 0; var < 5; var++) {
        cell_viscous_fluxes[cell_idx * 5 + var] = 0.0;
    }

    // Sum contributions from all faces
    for (int f = 0; f < num_faces; f++) {
        int face_id = cell_faces[cell_idx * max_faces + f];
        int orientation = face_orientations[cell_idx * max_faces + f];
        
        if (face_id >= 0) {
            for (int var = 0; var < 5; var++) {
                cell_viscous_fluxes[cell_idx * 5 + var] += 
                    orientation * face_viscous_fluxes[face_id * 5 + var];
            }
        }
    }
}

// Simple viscous dissipation kernel for artificial viscosity
__global__ void artificial_viscosity_kernel(
    const double* U_cells,         // [num_cells * 5] - Conservative variables
    const double* pressure,        // [num_cells] - Pressure field
    const int* cell_neighbors,     // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,      // [num_cells] - Number of neighbors
    double* artificial_viscosity,  // [num_cells * 5] - Output artificial viscosity
    double C2,                     // Artificial viscosity coefficient
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    double P_center = pressure[idx];
    
    for (int var = 0; var < 5; var++) {
        artificial_viscosity[idx * 5 + var] = 0.0;
    }

    // Calculate pressure-based artificial viscosity
    double max_pressure_jump = 0.0;
    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = cell_neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double P_neigh = pressure[neigh_id];
            double pressure_jump = fabs(P_neigh - P_center) / (P_neigh + P_center + 1e-12);
            max_pressure_jump = fmax(max_pressure_jump, pressure_jump);
        }
    }

    // Apply artificial viscosity if pressure jump is significant
    if (max_pressure_jump > 0.1) {
        double eps = C2 * max_pressure_jump;
        
        for (int n = 0; n < num_neigh; n++) {
            int neigh_id = cell_neighbors[idx * max_neighbors + n];
            if (neigh_id >= 0 && neigh_id < num_cells) {
                for (int var = 0; var < 5; var++) {
                    artificial_viscosity[idx * 5 + var] += 
                        eps * (U_cells[neigh_id * 5 + var] - U_cells[idx * 5 + var]);
                }
            }
        }
    }
}
