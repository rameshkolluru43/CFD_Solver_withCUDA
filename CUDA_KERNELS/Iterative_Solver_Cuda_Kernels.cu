// Iterative_Solver_Cuda_Kernels.cu
// Priority 2: Complete iterative solver CUDA kernels based on actual solver implementations
#include <cuda_runtime.h>
#include <math.h>

// Complete Poisson equation solver using Jacobi iteration
__global__ void poisson_jacobi_kernel(
    const double* phi_old,          // [num_cells] - Current solution
    const double* source_term,      // [num_cells] - Source term (RHS)
    const double* coefficients,     // [num_cells * max_neighbors] - Neighbor coefficients
    const int* neighbors,           // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,       // [num_cells] - Number of neighbors per cell
    const double* face_areas,       // [num_faces] - Face areas
    const double* distances,        // [num_cells * max_neighbors] - Distances to neighbors
    const double* cell_volumes,     // [num_cells] - Cell volumes
    double* phi_new,               // [num_cells] - Updated solution
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    double diagonal_coeff = 0.0;
    double off_diagonal_sum = 0.0;
    double volume = cell_volumes[idx];

    // Calculate coefficient matrix elements and RHS
    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double area = face_areas[idx * max_neighbors + n];
            double distance = distances[idx * max_neighbors + n];
            
            if (distance > 1e-12) {
                double coeff = area / distance;
                diagonal_coeff += coeff;
                off_diagonal_sum += coeff * phi_old[neigh_id];
            }
        }
    }

    // Jacobi update: phi_new = (RHS - off_diagonal_terms) / diagonal_coeff
    if (diagonal_coeff > 1e-12) {
        phi_new[idx] = (source_term[idx] * volume + off_diagonal_sum) / diagonal_coeff;
    } else {
        phi_new[idx] = phi_old[idx];
    }
}

// Gauss-Seidel iteration kernel (requires careful memory access)
__global__ void poisson_gauss_seidel_kernel(
    double* phi,                    // [num_cells] - Solution array (updated in-place)
    const double* source_term,      // [num_cells] - Source term
    const double* coefficients,     // [num_cells * max_neighbors] - Coefficients
    const int* neighbors,           // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,       // [num_cells] - Number of neighbors
    const double* face_areas,       // Face areas
    const double* distances,        // Distances to neighbors
    const double* cell_volumes,     // Cell volumes
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    double diagonal_coeff = 0.0;
    double off_diagonal_sum = 0.0;
    double volume = cell_volumes[idx];

    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double area = face_areas[idx * max_neighbors + n];
            double distance = distances[idx * max_neighbors + n];
            
            if (distance > 1e-12) {
                double coeff = area / distance;
                diagonal_coeff += coeff;
                off_diagonal_sum += coeff * phi[neigh_id]; // Uses latest values
            }
        }
    }

    if (diagonal_coeff > 1e-12) {
        phi[idx] = (source_term[idx] * volume + off_diagonal_sum) / diagonal_coeff;
    }
}

// Conjugate Gradient solver - Matrix-vector multiplication
__global__ void cg_matrix_vector_mult_kernel(
    const double* vector,           // [num_cells] - Input vector
    const double* coefficients,     // [num_cells * max_neighbors] - Matrix coefficients
    const int* neighbors,           // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,       // [num_cells] - Number of neighbors
    const double* face_areas,       // Face areas
    const double* distances,        // Distances
    double* result,                 // [num_cells] - A * vector
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    double diagonal_coeff = 0.0;
    double off_diagonal_sum = 0.0;

    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double area = face_areas[idx * max_neighbors + n];
            double distance = distances[idx * max_neighbors + n];
            
            if (distance > 1e-12) {
                double coeff = area / distance;
                diagonal_coeff += coeff;
                off_diagonal_sum += coeff * vector[neigh_id];
            }
        }
    }

    result[idx] = diagonal_coeff * vector[idx] - off_diagonal_sum;
}

// SOR (Successive Over-Relaxation) kernel
__global__ void poisson_sor_kernel(
    double* phi,                    // [num_cells] - Solution array
    const double* phi_old,          // [num_cells] - Previous iteration
    const double* source_term,      // [num_cells] - Source term
    const int* neighbors,           // [num_cells * max_neighbors] - Neighbor indices
    const int* num_neighbors,       // [num_cells] - Number of neighbors
    const double* face_areas,       // Face areas
    const double* distances,        // Distances
    const double* cell_volumes,     // Cell volumes
    double omega,                   // Relaxation parameter (1.0-2.0)
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    double diagonal_coeff = 0.0;
    double off_diagonal_sum = 0.0;
    double volume = cell_volumes[idx];

    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double area = face_areas[idx * max_neighbors + n];
            double distance = distances[idx * max_neighbors + n];
            
            if (distance > 1e-12) {
                double coeff = area / distance;
                diagonal_coeff += coeff;
                off_diagonal_sum += coeff * phi[neigh_id];
            }
        }
    }

    if (diagonal_coeff > 1e-12) {
        double phi_jacobi = (source_term[idx] * volume + off_diagonal_sum) / diagonal_coeff;
        phi[idx] = (1.0 - omega) * phi_old[idx] + omega * phi_jacobi;
    }
}

// Residual calculation kernel
__global__ void calculate_residual_kernel(
    const double* phi,              // [num_cells] - Current solution
    const double* source_term,      // [num_cells] - Source term
    const int* neighbors,           // [num_cells * max_neighbors] - Neighbors
    const int* num_neighbors,       // [num_cells] - Number of neighbors
    const double* face_areas,       // Face areas
    const double* distances,        // Distances
    const double* cell_volumes,     // Cell volumes
    double* residual,              // [num_cells] - Output residual
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_neigh = num_neighbors[idx];
    double diagonal_term = 0.0;
    double off_diagonal_term = 0.0;
    double volume = cell_volumes[idx];

    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double area = face_areas[idx * max_neighbors + n];
            double distance = distances[idx * max_neighbors + n];
            
            if (distance > 1e-12) {
                double coeff = area / distance;
                diagonal_term += coeff * phi[idx];
                off_diagonal_term += coeff * phi[neigh_id];
            }
        }
    }

    // Residual = RHS - A*phi
    residual[idx] = source_term[idx] * volume - (diagonal_term - off_diagonal_term);
}

// BiCGSTAB solver kernels
__global__ void bicgstab_vector_operations_kernel(
    const double* x,     // [num_cells] - Input vector 1
    const double* y,     // [num_cells] - Input vector 2  
    const double* z,     // [num_cells] - Input vector 3 (optional)
    double* result,      // [num_cells] - Output vector
    double alpha,        // Scalar coefficient 1
    double beta,         // Scalar coefficient 2
    int operation,       // 0: result = alpha*x + beta*y
                        // 1: result = x + alpha*y
                        // 2: result = x - alpha*y
                        // 3: result = alpha*x + beta*y + z
    int num_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    switch (operation) {
        case 0: // result = alpha*x + beta*y
            result[idx] = alpha * x[idx] + beta * y[idx];
            break;
        case 1: // result = x + alpha*y
            result[idx] = x[idx] + alpha * y[idx];
            break;
        case 2: // result = x - alpha*y
            result[idx] = x[idx] - alpha * y[idx];
            break;
        case 3: // result = alpha*x + beta*y + z
            result[idx] = alpha * x[idx] + beta * y[idx] + z[idx];
            break;
    }
}

// Momentum equation solver kernel (semi-implicit)
__global__ void momentum_solver_kernel(
    const double* u_old,            // [num_cells] - Previous velocity component
    const double* pressure_grad,    // [num_cells] - Pressure gradient
    const double* viscous_term,     // [num_cells] - Viscous term
    const double* source_term,      // [num_cells] - Source/body force term
    const double* density,          // [num_cells] - Density field
    double* u_new,                 // [num_cells] - Updated velocity component
    double dt,                     // Time step
    double under_relaxation,       // Under-relaxation factor
    int num_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    double rho = density[idx];
    if (rho > 1e-12) {
        // Semi-implicit momentum update
        double u_explicit = u_old[idx] + dt * (
            -pressure_grad[idx] / rho +
            viscous_term[idx] / rho +
            source_term[idx]
        );
        
        // Apply under-relaxation
        u_new[idx] = under_relaxation * u_explicit + (1.0 - under_relaxation) * u_old[idx];
    } else {
        u_new[idx] = u_old[idx];
    }
}

// Pressure correction equation solver (SIMPLE algorithm)
__global__ void pressure_correction_kernel(
    const double* u_star,           // [num_cells] - Intermediate velocity (u-component)
    const double* v_star,           // [num_cells] - Intermediate velocity (v-component)  
    const double* w_star,           // [num_cells] - Intermediate velocity (w-component)
    const double* density,          // [num_cells] - Density field
    const int* neighbors,           // [num_cells * max_neighbors] - Neighbors
    const int* num_neighbors,       // [num_cells] - Number of neighbors
    const double* face_areas,       // Face areas
    const double* distances,        // Distances
    const double* cell_volumes,     // Cell volumes
    double* pressure_correction,   // [num_cells] - Pressure correction
    int num_cells,
    int max_neighbors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    // Calculate mass source term from velocity divergence
    double mass_source = 0.0;
    int num_neigh = num_neighbors[idx];
    
    for (int n = 0; n < num_neigh; n++) {
        int neigh_id = neighbors[idx * max_neighbors + n];
        if (neigh_id >= 0 && neigh_id < num_cells) {
            double area = face_areas[idx * max_neighbors + n];
            
            // Face velocity (average)
            double u_face = 0.5 * (u_star[idx] + u_star[neigh_id]);
            double v_face = 0.5 * (v_star[idx] + v_star[neigh_id]);
            double w_face = 0.5 * (w_star[idx] + w_star[neigh_id]);
            double rho_face = 0.5 * (density[idx] + density[neigh_id]);
            
            // Mass flux through face (simplified)
            mass_source += rho_face * (u_face + v_face + w_face) * area;
        }
    }

    // Store mass source as initial pressure correction RHS
    pressure_correction[idx] = -mass_source;
}

// Atomic add for double (pre-sm_60 architectures) using CAS loop
__device__ inline double atomicAdd_double(double* address, double val) {
#if __CUDA_ARCH__ >= 600
    return atomicAdd(address, val);
#else
    unsigned long long int* addr_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *addr_as_ull, assumed;
    do {
        assumed = old;
        double updated = val + __longlong_as_double(assumed);
        old = atomicCAS(addr_as_ull, assumed, __double_as_longlong(updated));
    } while (assumed != old);
    return __longlong_as_double(old);
#endif
}

// Atomic max for double via CAS loop
__device__ inline void atomicMax_double(double* address, double val) {
    unsigned long long int* addr_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *addr_as_ull, assumed;
    while (__longlong_as_double(old) < val) {
        assumed = old;
        old = atomicCAS(addr_as_ull, assumed, __double_as_longlong(val));
        if (assumed == old) break; // success
    }
}

// Convergence check kernel
__global__ void convergence_check_kernel(
    const double* residual,         // [num_cells] - Residual vector
    double* max_residual,          // [1] - Output maximum residual
    double* sum_residual,          // [1] - Output sum of residuals
    int num_cells
) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Load data
    double local_max = 0.0;
    double local_sum = 0.0;
    
    if (idx < num_cells) {
        double abs_res = fabs(residual[idx]);
        local_max = abs_res;
        local_sum = abs_res;
    }
    
    sdata[tid] = local_max;
    sdata[tid + blockDim.x] = local_sum;
    __syncthreads();
    
    // Reduction for max
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
            sdata[tid + blockDim.x] += sdata[tid + s + blockDim.x];
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        atomicMax_double(max_residual, sdata[0]);
        atomicAdd_double(sum_residual, sdata[blockDim.x]);
    }
}
