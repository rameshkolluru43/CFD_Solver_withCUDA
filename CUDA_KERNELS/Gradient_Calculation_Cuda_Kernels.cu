// Gradient_Calculation_Cuda_Kernels.cu
// Priority 1: Complete Gradient calculation CUDA kernels for Del operations
// Note: We locally define gamma to avoid collision with the math library gamma() function.
#include <cuda_runtime.h>
#include <math.h>
#ifndef gamma
#define gamma 1.4
#endif

// Kernel for calculating Q gradients using Green's theorem
__global__ void calculate_Q_gradients_kernel(
    const double* U_cells,        // [num_cells * 5] - Conservative variables
    const double* face_QnDS,      // [num_faces * 5] - Face Q*n*dS values
    const int* cell_faces,        // [num_cells * max_faces] - Face indices for each cell
    const int* num_faces_per_cell,// [num_cells] - Number of faces per cell
    const double* cell_volumes,   // [num_cells] - Cell volumes
    double* Q_gradients,          // [num_cells * 5 * 3] - Output gradients (x,y,z for each variable)
    int num_cells,
    int max_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    double inv_vol = 1.0 / cell_volumes[idx];
    int num_faces = num_faces_per_cell[idx];
    
    // Initialize gradients to zero
    for (int var = 0; var < 5; var++) {
        for (int dir = 0; dir < 3; dir++) {
            Q_gradients[idx * 15 + var * 3 + dir] = 0.0;
        }
    }

    // Sum face contributions using Green's theorem
    for (int f = 0; f < num_faces; f++) {
        int face_id = cell_faces[idx * max_faces + f];
        if (face_id >= 0) {
            for (int var = 0; var < 5; var++) {
                for (int dir = 0; dir < 3; dir++) {
                    Q_gradients[idx * 15 + var * 3 + dir] += 
                        face_QnDS[face_id * 15 + var * 3 + dir];
                }
            }
        }
    }

    // Apply inverse volume
    for (int var = 0; var < 5; var++) {
        for (int dir = 0; dir < 3; dir++) {
            Q_gradients[idx * 15 + var * 3 + dir] *= inv_vol;
        }
    }
}

// Kernel for calculating velocity gradients at cell center
__global__ void calculate_velocity_gradients_kernel(
    const double* U_cells,        // [num_cells * 5] - Conservative variables
    const double* face_undS,      // [num_faces * 3] - Face u*n*dS values
    const double* face_vndS,      // [num_faces * 3] - Face v*n*dS values  
    const double* face_wndS,      // [num_faces * 3] - Face w*n*dS values
    const double* face_TndS,      // [num_faces * 3] - Face T*n*dS values
    const int* cell_faces,        // [num_cells * max_faces] - Face indices
    const int* num_faces_per_cell,// [num_cells] - Number of faces per cell
    const double* cell_volumes,   // [num_cells] - Cell volumes
    double* velocity_gradients,   // [num_cells * 9] - u,v,w gradients (3x3 tensor)
    double* temperature_gradients,// [num_cells * 3] - Temperature gradients
    int num_cells,
    int max_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    double inv_vol = 1.0 / cell_volumes[idx];
    int num_faces = num_faces_per_cell[idx];

    // Initialize gradients
    for (int i = 0; i < 9; i++) velocity_gradients[idx * 9 + i] = 0.0;
    for (int i = 0; i < 3; i++) temperature_gradients[idx * 3 + i] = 0.0;

    // Sum face contributions
    for (int f = 0; f < num_faces; f++) {
        int face_id = cell_faces[idx * max_faces + f];
        if (face_id >= 0) {
            // Velocity gradients
            for (int dir = 0; dir < 3; dir++) {
                velocity_gradients[idx * 9 + 0 * 3 + dir] += face_undS[face_id * 3 + dir];
                velocity_gradients[idx * 9 + 1 * 3 + dir] += face_vndS[face_id * 3 + dir];
                velocity_gradients[idx * 9 + 2 * 3 + dir] += face_wndS[face_id * 3 + dir];
            }
            
            // Temperature gradients
            for (int dir = 0; dir < 3; dir++) {
                temperature_gradients[idx * 3 + dir] += face_TndS[face_id * 3 + dir];
            }
        }
    }

    // Apply inverse volume
    for (int i = 0; i < 9; i++) velocity_gradients[idx * 9 + i] *= inv_vol;
    for (int i = 0; i < 3; i++) temperature_gradients[idx * 3 + i] *= inv_vol;
}

// Kernel for calculating second-order dissipation (Del2_Q)
__global__ void calculate_del2_Q_kernel(
    const double* U_cells,        // [num_cells * 5] - Conservative variables
    const int* cell_neighbors,    // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,     // [num_cells] - Number of neighbors per cell
    double* Del2_Q,              // [num_cells * 5 * 3] - Second-order differences
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    
    // Calculate Del2_Q for each direction and variable
    for (int var = 0; var < 5; var++) {
        for (int dir = 0; dir < 3; dir++) {
            double sum = 0.0;
            int valid_neighbors = 0;
            
            // Central difference approximation
            for (int n = 0; n < num_neigh; n++) {
                int neigh_id = cell_neighbors[idx * max_neighbors + n];
                if (neigh_id >= 0 && neigh_id < num_cells) {
                    sum += U_cells[neigh_id * 5 + var] - U_cells[idx * 5 + var];
                    valid_neighbors++;
                }
            }
            
            if (valid_neighbors > 0) {
                Del2_Q[idx * 15 + var * 3 + dir] = sum / valid_neighbors;
            } else {
                Del2_Q[idx * 15 + var * 3 + dir] = 0.0;
            }
        }
    }
}

// Kernel for calculating fourth-order dissipation (Del4_Q)
__global__ void calculate_del4_Q_kernel(
    const double* Del2_Q,         // [num_cells * 5 * 3] - Second-order differences
    const int* cell_neighbors,    // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,     // [num_cells] - Number of neighbors per cell
    double* Del4_Q,              // [num_cells * 5 * 3] - Fourth-order differences
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    
    // Calculate Del4_Q as second difference of Del2_Q
    for (int var = 0; var < 5; var++) {
        for (int dir = 0; dir < 3; dir++) {
            double sum = 0.0;
            int valid_neighbors = 0;
            
            for (int n = 0; n < num_neigh; n++) {
                int neigh_id = cell_neighbors[idx * max_neighbors + n];
                if (neigh_id >= 0 && neigh_id < num_cells) {
                    sum += Del2_Q[neigh_id * 15 + var * 3 + dir] - Del2_Q[idx * 15 + var * 3 + dir];
                    valid_neighbors++;
                }
            }
            
            if (valid_neighbors > 0) {
                Del4_Q[idx * 15 + var * 3 + dir] = sum / valid_neighbors;
            } else {
                Del4_Q[idx * 15 + var * 3 + dir] = 0.0;
            }
        }
    }
}

// Kernel for calculating sixth-order dissipation (Del6_Q) 
__global__ void calculate_del6_Q_kernel(
    const double* Del4_Q,         // [num_cells * 5 * 3] - Fourth-order differences
    const int* cell_neighbors,    // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,     // [num_cells] - Number of neighbors per cell
    double* Del6_Q,              // [num_cells * 5 * 3] - Sixth-order differences
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    
    // Calculate Del6_Q as second difference of Del4_Q
    for (int var = 0; var < 5; var++) {
        for (int dir = 0; dir < 3; dir++) {
            double sum = 0.0;
            int valid_neighbors = 0;
            
            for (int n = 0; n < num_neigh; n++) {
                int neigh_id = cell_neighbors[idx * max_neighbors + n];
                if (neigh_id >= 0 && neigh_id < num_cells) {
                    sum += Del4_Q[neigh_id * 15 + var * 3 + dir] - Del4_Q[idx * 15 + var * 3 + dir];
                    valid_neighbors++;
                }
            }
            
            if (valid_neighbors > 0) {
                Del6_Q[idx * 15 + var * 3 + dir] = sum / valid_neighbors;
            } else {
                Del6_Q[idx * 15 + var * 3 + dir] = 0.0;
            }
        }
    }
}

// Kernel for calculating primitive variable gradients
__global__ void calculate_primitive_gradients_kernel(
    const double* U_cells,          // [num_cells * 5] - Conservative variables
    const double* Q_gradients,      // [num_cells * 5 * 3] - Conservative gradients
    double* primitive_gradients,    // [num_cells * 5 * 3] - Primitive gradients (rho,u,v,w,P)
    int num_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    // Get conservative variables
    double rho = U_cells[idx * 5 + 0];
    double rhou = U_cells[idx * 5 + 1];
    double rhov = U_cells[idx * 5 + 2];
    double rhow = U_cells[idx * 5 + 3];
    double E = U_cells[idx * 5 + 4];

    // Calculate primitive variables
    double u = rhou / rho;
    double v = rhov / rho;
    double w = rhow / rho;
    double P = (gamma - 1.0) * (E - 0.5 * rho * (u*u + v*v + w*w));

    double inv_rho = 1.0 / rho;
    double inv_rho2 = inv_rho * inv_rho;

    // Transform gradients from conservative to primitive
    for (int dir = 0; dir < 3; dir++) {
        double drho_dx = Q_gradients[idx * 15 + 0 * 3 + dir];
        double drhou_dx = Q_gradients[idx * 15 + 1 * 3 + dir];
        double drhov_dx = Q_gradients[idx * 15 + 2 * 3 + dir];
        double drhow_dx = Q_gradients[idx * 15 + 3 * 3 + dir];
        double dE_dx = Q_gradients[idx * 15 + 4 * 3 + dir];

        // Density gradient
        primitive_gradients[idx * 15 + 0 * 3 + dir] = drho_dx;

        // Velocity gradients
        primitive_gradients[idx * 15 + 1 * 3 + dir] = inv_rho * (drhou_dx - u * drho_dx);
        primitive_gradients[idx * 15 + 2 * 3 + dir] = inv_rho * (drhov_dx - v * drho_dx);
        primitive_gradients[idx * 15 + 3 * 3 + dir] = inv_rho * (drhow_dx - w * drho_dx);

        // Pressure gradient
        double du_dx = primitive_gradients[idx * 15 + 1 * 3 + dir];
        double dv_dx = primitive_gradients[idx * 15 + 2 * 3 + dir];
        double dw_dx = primitive_gradients[idx * 15 + 3 * 3 + dir];

        primitive_gradients[idx * 15 + 4 * 3 + dir] = (gamma - 1.0) * 
            (dE_dx - 0.5 * (drho_dx * (u*u + v*v + w*w) + 2.0 * rho * (u*du_dx + v*dv_dx + w*dw_dx)));
    }
}