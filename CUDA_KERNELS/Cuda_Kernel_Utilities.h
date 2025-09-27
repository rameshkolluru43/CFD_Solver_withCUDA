// Cuda_Kernel_Utilities.h
// Comprehensive CUDA kernel declarations and utility functions
#ifndef CUDA_KERNEL_UTILITIES_H
#define CUDA_KERNEL_UTILITIES_H

#include <cuda_runtime.h>

// ===== TIME INTEGRATION KERNELS =====
__global__ void rk4_complete_kernel(double *U_cells, double *U_cells_temp, double *Net_Flux,
                                    double *Viscous_Flux, double *Cell_Volumes, double *Del2_Q,
                                    double *Del4_Q, double *Cell_Diagonal, double dt, double lambda,
                                    int num_cells, int stage);
__global__ void tvd_rk3_kernel(double *U_cells, double *U_temp1, double *U_temp2, double *Net_Flux,
                               double *Viscous_Flux, double *Cell_Volumes, double dt, int num_cells, int stage);
__global__ void euler_time_step_kernel(const double *U_old, const double *flux, double *U_new,
                                       const double *volumes, double dt, int num_cells);

// ===== FLUX CALCULATION KERNELS =====
// AUSM flux scheme
__global__ void ausm_flux_kernel(const double *U_left, const double *U_right, const double *face_normals,
                                 const double *face_areas, double *fluxes, int num_faces);
// Van Leer flux splitting
__global__ void vanleer_flux_kernel(const double *U_left, const double *U_right, const double *face_normals,
                                    const double *face_areas, double *fluxes, int num_faces);
// Central difference flux
__global__ void central_flux_kernel(const double *U_left, const double *U_right, const double *face_normals,
                                    const double *face_areas, double *fluxes, int num_faces);

// ===== GRADIENT CALCULATION KERNELS =====
__global__ void calculate_Q_gradients_kernel(const double *U_cells, const double *face_QnDS,
                                             const int *cell_faces, const int *num_faces_per_cell,
                                             const double *cell_volumes, double *Q_gradients,
                                             int num_cells, int max_faces);
__global__ void calculate_velocity_gradients_kernel(const double *U_cells, const double *face_undS,
                                                    const double *face_vndS, const double *face_wndS,
                                                    const double *face_TndS, const int *cell_faces,
                                                    const int *num_faces_per_cell, const double *cell_volumes,
                                                    double *velocity_gradients, double *temperature_gradients,
                                                    int num_cells, int max_faces);
__global__ void calculate_del2_Q_kernel(const double *U_cells, const int *cell_neighbors,
                                        const int *num_neighbors, double *Del2_Q, int num_cells, int max_neighbors);
__global__ void calculate_del4_Q_kernel(const double *Del2_Q, const int *cell_neighbors,
                                        const int *num_neighbors, double *Del4_Q, int num_cells, int max_neighbors);
__global__ void calculate_del6_Q_kernel(const double *Del4_Q, const int *cell_neighbors,
                                        const int *num_neighbors, double *Del6_Q, int num_cells, int max_neighbors);
__global__ void calculate_primitive_gradients_kernel(const double *U_cells, const double *Q_gradients,
                                                     double *primitive_gradients, int num_cells);

// ===== VISCOUS FLUX KERNELS =====
__global__ void calculate_viscous_stress_kernel(const double *velocity_gradients, const double *temperature_gradients,
                                                const double *primitive_vars, const double *face_normals,
                                                const double *face_areas, const int *face_cells,
                                                double *viscous_fluxes, int num_faces);
__global__ void accumulate_cell_viscous_flux_kernel(const double *face_viscous_fluxes, const int *cell_faces,
                                                    const int *num_faces_per_cell, const int *face_orientations,
                                                    double *cell_viscous_fluxes, int num_cells, int max_faces);
__global__ void artificial_viscosity_kernel(const double *U_cells, const double *pressure,
                                            const int *cell_neighbors, const int *num_neighbors,
                                            double *artificial_viscosity, double C2, int num_cells, int max_neighbors);

// ===== GEOMETRY CALCULATION KERNELS =====
__global__ void calculate_cell_volumes_kernel(const double *cell_coords, const double *face_centers,
                                              const double *face_areas, const double *face_normals,
                                              const int *cell_faces, const int *num_faces_per_cell,
                                              double *cell_volumes, int num_cells, int max_faces);
__global__ void calculate_quad_face_properties_kernel(const double *face_coords, double *face_areas,
                                                      double *face_normals, double *face_centers, int num_faces);
__global__ void calculate_triangle_face_properties_kernel(const double *face_coords, double *face_areas,
                                                          double *face_normals, double *face_centers, int num_faces);
__global__ void calculate_cell_centers_kernel(const double *cell_coords, const int *num_points_per_cell,
                                              double *cell_centers, int num_cells, int max_points);
__global__ void calculate_cell_distances_kernel(const double *cell_centers, const int *cell_faces,
                                                const int *face_cells, const int *num_faces_per_cell,
                                                double *cell_distances, int num_cells, int max_faces);
__global__ void calculate_hex_volumes_kernel(const double *cell_coords, double *cell_volumes, int num_cells);
__global__ void calculate_mesh_quality_kernel(const double *cell_volumes, const double *face_areas,
                                              const int *cell_faces, const int *num_faces_per_cell,
                                              double *aspect_ratios, double *skewness, int num_cells, int max_faces);

// ===== GRID PROCESSING KERNELS =====
// Include Grid CUDA kernel declarations
#include "Grid_Cuda_Kernels.h"

// ===== LINEAR ALGEBRA KERNELS =====
__global__ void vector_add_kernel(const double *a, const double *b, double *c, int n);
__global__ void vector_subtract_kernel(const double *a, const double *b, double *c, int n);
__global__ void vector_scalar_mult_kernel(const double *a, double scalar, double *result, int n);
__global__ void vector_dot_product_kernel(const double *a, const double *b, double *partial_results, int n);
__global__ void vector_cross_product_kernel(const double *a, const double *b, double *result, int num_vectors);
__global__ void vector_magnitude_kernel(const double *vectors, double *magnitudes, int num_vectors);
__global__ void vector_normalize_kernel(const double *vectors, double *unit_vectors, int num_vectors);
__global__ void matrix_vector_mult_kernel(const double *matrix, const double *vector, double *result,
                                          int num_systems, int n);
__global__ void vector_distance_kernel(const double *points_a, const double *points_b, double *distances, int num_pairs);
__global__ void vector_interpolation_kernel(const double *vec_a, const double *vec_b, const double *weights,
                                            double *result, int num_vectors);
__global__ void vector_component_ops_kernel(const double *vectors, double *components, int operation, int num_vectors);
__global__ void tensor_operations_kernel(const double *tensors, double *results, int operation, int num_tensors);
__global__ void advanced_reduction_kernel(const double *input, double *output, int operation, int n);

// ===== ITERATIVE SOLVER KERNELS =====
__global__ void poisson_jacobi_kernel(const double *phi_old, const double *source_term, const double *coefficients,
                                      const int *neighbors, const int *num_neighbors, const double *face_areas,
                                      const double *distances, const double *cell_volumes, double *phi_new,
                                      int num_cells, int max_neighbors);
__global__ void poisson_gauss_seidel_kernel(double *phi, const double *source_term, const double *coefficients,
                                            const int *neighbors, const int *num_neighbors, const double *face_areas,
                                            const double *distances, const double *cell_volumes,
                                            int num_cells, int max_neighbors);
__global__ void cg_matrix_vector_mult_kernel(const double *vector, const double *coefficients, const int *neighbors,
                                             const int *num_neighbors, const double *face_areas, const double *distances,
                                             double *result, int num_cells, int max_neighbors);
__global__ void poisson_sor_kernel(double *phi, const double *phi_old, const double *source_term,
                                   const int *neighbors, const int *num_neighbors, const double *face_areas,
                                   const double *distances, const double *cell_volumes, double omega,
                                   int num_cells, int max_neighbors);
__global__ void calculate_residual_kernel(const double *phi, const double *source_term, const int *neighbors,
                                          const int *num_neighbors, const double *face_areas, const double *distances,
                                          const double *cell_volumes, double *residual, int num_cells, int max_neighbors);
__global__ void bicgstab_vector_operations_kernel(const double *x, const double *y, const double *z, double *result,
                                                  double alpha, double beta, int operation, int num_cells);
__global__ void momentum_solver_kernel(const double *u_old, const double *pressure_grad, const double *viscous_term,
                                       const double *source_term, const double *density, double *u_new,
                                       double dt, double under_relaxation, int num_cells);
__global__ void pressure_correction_kernel(const double *u_star, const double *v_star, const double *w_star,
                                           const double *density, const int *neighbors, const int *num_neighbors,
                                           const double *face_areas, const double *distances, const double *cell_volumes,
                                           double *pressure_correction, int num_cells, int max_neighbors);
__global__ void convergence_check_kernel(const double *residual, double *max_residual, double *sum_residual, int num_cells);

// ===== ADVANCED OPTIMIZATION KERNELS =====
__global__ void shared_memory_reduction_kernel(const double *input, double *output, int n);

// ===== UTILITY FUNCTIONS =====
inline void cudaCheckError(cudaError_t err, const char *msg)
{
    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s: %s\n", msg, cudaGetErrorString(err));
        exit(1);
    }
}

inline void cudaCheckKernel(const char *kernel_name)
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("CUDA Kernel Error in %s: %s\n", kernel_name, cudaGetErrorString(err));
        exit(1);
    }
}

// Macro for convenient error checking
#define CUDA_CHECK(call) cudaCheckError(call, #call)
#define CUDA_CHECK_KERNEL(kernel) cudaCheckKernel(#kernel)

// Common kernel launch configuration helper
inline dim3 calculateGridSize(int num_elements, int block_size = 256)
{
    return dim3((num_elements + block_size - 1) / block_size, 1, 1);
}

inline dim3 calculateBlockSize(int block_size = 256)
{
    return dim3(block_size, 1, 1);
}

#endif // CUDA_KERNEL_UTILITIES_H
