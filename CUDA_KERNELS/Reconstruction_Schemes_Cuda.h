/**
 * @file Reconstruction_Schemes_Cuda.h
 * @brief Header for MUSCL and WENO reconstruction schemes
 * @date 2026-01-15
 */

#ifndef RECONSTRUCTION_SCHEMES_CUDA_H
#define RECONSTRUCTION_SCHEMES_CUDA_H

#include <cuda_runtime.h>

//=============================================================================
// KERNEL DECLARATIONS - MUSCL
//=============================================================================

__global__ void muscl_reconstruction_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *cell_distances,
    const double *limited_gradients,
    double *Q_left,
    double *Q_right,
    int num_cells,
    double kappa);

//=============================================================================
// KERNEL DECLARATIONS - WENO5
//=============================================================================

__global__ void weno5_reconstruction_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    double *Q_left,
    double *Q_right,
    int num_cells,
    double epsilon,
    int p);

__global__ void weno5_characteristic_reconstruction_kernel(
    const double *U_cells,
    const int *cell_neighbors,
    const double *L_matrices,
    const double *R_matrices,
    double *Q_left,
    double *Q_right,
    int num_cells,
    double epsilon,
    int p);

//=============================================================================
// HOST WRAPPER FUNCTIONS
//=============================================================================

cudaError_t launch_muscl_reconstruction(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_cell_distances,
    const double *d_limited_gradients,
    double *d_Q_left,
    double *d_Q_right,
    int num_cells,
    double kappa = 0.5, // QUICK scheme by default
    int block_size = 256);

cudaError_t launch_weno5_reconstruction(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    double *d_Q_left,
    double *d_Q_right,
    int num_cells,
    double epsilon = 1e-6,
    int p = 2,
    int block_size = 256);

cudaError_t launch_weno5_characteristic_reconstruction(
    const double *d_U_cells,
    const int *d_cell_neighbors,
    const double *d_L_matrices,
    const double *d_R_matrices,
    double *d_Q_left,
    double *d_Q_right,
    int num_cells,
    double epsilon = 1e-6,
    int p = 2,
    int block_size = 256);

#endif // RECONSTRUCTION_SCHEMES_CUDA_H
