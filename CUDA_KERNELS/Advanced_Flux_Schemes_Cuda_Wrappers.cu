/**
 * @file Advanced_Flux_Schemes_Cuda_Wrappers.cu
 * @brief Host wrapper implementations for Roe, HLLC, and LLF flux schemes
 * @date 2026-01-15
 */

#include "Advanced_Flux_Schemes_Cuda.h"
#include <stdio.h>

//=============================================================================
// ROE FLUX WRAPPERS
//=============================================================================

cudaError_t launch_roe_flux(
    const double* d_U_cells,
    const double* d_P_cells,
    const int* d_face_neighbors,
    const double* d_face_normals,
    const double* d_face_areas,
    double* d_flux_output,
    int num_cells,
    double gamma,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    roe_flux_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_P_cells,
        d_face_neighbors,
        d_face_normals,
        d_face_areas,
        d_flux_output,
        num_cells,
        gamma
    );
    
    return cudaGetLastError();
}

cudaError_t launch_roe_flux_2nd_order(
    const double* d_U_cells,
    const double* d_P_cells,
    const int* d_face_neighbors,
    const double* d_face_normals,
    const double* d_face_areas,
    const double* d_limited_gradients,
    const double* d_cell_centers,
    const double* d_face_centers,
    double* d_flux_output,
    int num_cells,
    double gamma,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    roe_flux_2nd_order_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_P_cells,
        d_face_neighbors,
        d_face_normals,
        d_face_areas,
        d_limited_gradients,
        d_cell_centers,
        d_face_centers,
        d_flux_output,
        num_cells,
        gamma
    );
    
    return cudaGetLastError();
}

//=============================================================================
// HLLC FLUX WRAPPERS
//=============================================================================

cudaError_t launch_hllc_flux(
    const double* d_U_cells,
    const double* d_P_cells,
    const int* d_face_neighbors,
    const double* d_face_normals,
    const double* d_face_areas,
    double* d_flux_output,
    int num_cells,
    double gamma,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    hllc_flux_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_P_cells,
        d_face_neighbors,
        d_face_normals,
        d_face_areas,
        d_flux_output,
        num_cells,
        gamma
    );
    
    return cudaGetLastError();
}

//=============================================================================
// LLF FLUX WRAPPERS
//=============================================================================

cudaError_t launch_llf_flux(
    const double* d_U_cells,
    const double* d_P_cells,
    const int* d_face_neighbors,
    const double* d_face_normals,
    const double* d_face_areas,
    double* d_flux_output,
    int num_cells,
    double gamma,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    llf_flux_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_P_cells,
        d_face_neighbors,
        d_face_normals,
        d_face_areas,
        d_flux_output,
        num_cells,
        gamma
    );
    
    return cudaGetLastError();
}
