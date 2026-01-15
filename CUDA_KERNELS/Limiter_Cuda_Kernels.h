/**
 * @file Limiter_Cuda_Kernels.h
 * @brief Header file for CUDA slope limiter kernels
 * @author AI Assistant
 * @date 2026-01-14
 * @version 1.0
 */

#ifndef LIMITER_CUDA_KERNELS_H
#define LIMITER_CUDA_KERNELS_H

#include <cuda_runtime.h>

//=============================================================================
// LIMITER TYPE ENUMERATION
//=============================================================================

/**
 * @brief Enumeration for limiter types
 */
enum LimiterType
{
    LIMITER_MINMOD = 0,
    LIMITER_VANLEER = 1,
    LIMITER_SUPERBEE = 2,
    LIMITER_VANALBADA = 3,
    LIMITER_VENKATAKRISHNAN = 4,
    LIMITER_BARTH_JESPERSEN = 5
};

//=============================================================================
// CUDA KERNEL DECLARATIONS
//=============================================================================

/**
 * @brief MinMod limiter kernel
 */
__global__ void minmod_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    double *limited_gradients,
    int num_cells,
    double limiter_zeta,
    int use_three_arg);

/**
 * @brief Van Leer limiter kernel
 */
__global__ void vanleer_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    double *limited_gradients,
    int num_cells,
    double limiter_zeta);

/**
 * @brief Superbee limiter kernel
 */
__global__ void superbee_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    double *limited_gradients,
    int num_cells,
    double limiter_zeta);

/**
 * @brief Van Albada limiter kernel
 */
__global__ void vanalbada_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    double *limited_gradients,
    int num_cells,
    double limiter_zeta);

/**
 * @brief Venkatakrishnan limiter kernel
 */
__global__ void venkatakrishnan_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    const double *cell_volumes,
    double *limited_gradients,
    int num_cells,
    double K_constant);

/**
 * @brief Barth-Jespersen limiter kernel
 */
__global__ void barth_jespersen_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    const double *unlimited_gradients,
    double *limited_gradients,
    int num_cells);

/**
 * @brief Generic limiter kernel with runtime selection
 */
__global__ void generic_limiter_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    double *limited_gradients,
    int num_cells,
    double limiter_zeta,
    int limiter_type);

//=============================================================================
// HOST WRAPPER FUNCTIONS
//=============================================================================

/**
 * @brief Host wrapper for MinMod limiter
 */
cudaError_t launch_minmod_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    double *d_limited_gradients,
    int num_cells,
    double limiter_zeta = 1.0,
    bool use_three_arg = false,
    int block_size = 256);

/**
 * @brief Host wrapper for Van Leer limiter
 */
cudaError_t launch_vanleer_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    double *d_limited_gradients,
    int num_cells,
    double limiter_zeta = 1.0,
    int block_size = 256);

/**
 * @brief Host wrapper for Superbee limiter
 */
cudaError_t launch_superbee_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    double *d_limited_gradients,
    int num_cells,
    double limiter_zeta = 1.0,
    int block_size = 256);

/**
 * @brief Host wrapper for Van Albada limiter
 */
cudaError_t launch_vanalbada_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    double *d_limited_gradients,
    int num_cells,
    double limiter_zeta = 1.0,
    int block_size = 256);

/**
 * @brief Host wrapper for Venkatakrishnan limiter
 */
cudaError_t launch_venkatakrishnan_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    const double *d_cell_volumes,
    double *d_limited_gradients,
    int num_cells,
    double K_constant = 1.0,
    int block_size = 256);

/**
 * @brief Host wrapper for Barth-Jespersen limiter
 */
cudaError_t launch_barth_jespersen_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    const double *d_unlimited_gradients,
    double *d_limited_gradients,
    int num_cells,
    int block_size = 256);

/**
 * @brief Host wrapper for generic limiter with runtime selection
 */
cudaError_t launch_generic_limiter(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    double *d_limited_gradients,
    int num_cells,
    LimiterType limiter_type,
    double limiter_zeta = 1.0,
    int block_size = 256);

/**
 * @brief Utility function to get limiter name string
 */
const char *get_limiter_name(LimiterType limiter_type);

#endif // LIMITER_CUDA_KERNELS_H
