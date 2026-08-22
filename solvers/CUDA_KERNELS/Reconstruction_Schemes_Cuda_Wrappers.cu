/**
 * @file Reconstruction_Schemes_Cuda_Wrappers.cu
 * @brief Host wrapper implementations for MUSCL and WENO reconstruction
 * @date 2026-01-15
 */

#include "Reconstruction_Schemes_Cuda.h"
#include <stdio.h>

//=============================================================================
// MUSCL RECONSTRUCTION WRAPPERS
//=============================================================================

cudaError_t launch_muscl_reconstruction(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    const double* d_limited_gradients,
    double* d_Q_left,
    double* d_Q_right,
    int num_cells,
    double kappa,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    muscl_reconstruction_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_limited_gradients,
        d_Q_left,
        d_Q_right,
        num_cells,
        kappa
    );
    
    return cudaGetLastError();
}

//=============================================================================
// WENO5 RECONSTRUCTION WRAPPERS
//=============================================================================

cudaError_t launch_weno5_reconstruction(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    double* d_Q_left,
    double* d_Q_right,
    int num_cells,
    double epsilon,
    int p,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    weno5_reconstruction_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_Q_left,
        d_Q_right,
        num_cells,
        epsilon,
        p
    );
    
    return cudaGetLastError();
}

cudaError_t launch_weno5_characteristic_reconstruction(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_L_matrices,
    const double* d_R_matrices,
    double* d_Q_left,
    double* d_Q_right,
    int num_cells,
    double epsilon,
    int p,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    weno5_characteristic_reconstruction_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_L_matrices,
        d_R_matrices,
        d_Q_left,
        d_Q_right,
        num_cells,
        epsilon,
        p
    );
    
    return cudaGetLastError();
}
