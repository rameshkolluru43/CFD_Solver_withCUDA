/**
 * @file Limiter_Cuda_Wrappers.cu
 * @brief Host wrapper implementations for limiter CUDA kernels
 * @author AI Assistant
 * @date 2026-01-14
 */

#include "Limiter_Cuda_Kernels.h"
#include <stdio.h>

//=============================================================================
// HOST WRAPPER IMPLEMENTATIONS
//=============================================================================

cudaError_t launch_minmod_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    double* d_limited_gradients,
    int num_cells,
    double limiter_zeta,
    bool use_three_arg,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    minmod_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_limited_gradients,
        num_cells,
        limiter_zeta,
        use_three_arg ? 1 : 0
    );
    
    return cudaGetLastError();
}

cudaError_t launch_vanleer_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    double* d_limited_gradients,
    int num_cells,
    double limiter_zeta,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    vanleer_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_limited_gradients,
        num_cells,
        limiter_zeta
    );
    
    return cudaGetLastError();
}

cudaError_t launch_superbee_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    double* d_limited_gradients,
    int num_cells,
    double limiter_zeta,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    superbee_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_limited_gradients,
        num_cells,
        limiter_zeta
    );
    
    return cudaGetLastError();
}

cudaError_t launch_vanalbada_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    double* d_limited_gradients,
    int num_cells,
    double limiter_zeta,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    vanalbada_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_limited_gradients,
        num_cells,
        limiter_zeta
    );
    
    return cudaGetLastError();
}

cudaError_t launch_venkatakrishnan_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    const double* d_cell_volumes,
    double* d_limited_gradients,
    int num_cells,
    double K_constant,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    venkatakrishnan_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_cell_volumes,
        d_limited_gradients,
        num_cells,
        K_constant
    );
    
    return cudaGetLastError();
}

cudaError_t launch_barth_jespersen_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    const double* d_unlimited_gradients,
    double* d_limited_gradients,
    int num_cells,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    barth_jespersen_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_unlimited_gradients,
        d_limited_gradients,
        num_cells
    );
    
    return cudaGetLastError();
}

cudaError_t launch_generic_limiter(
    const double* d_U_cells,
    const int* d_cell_neighbors,
    const double* d_cell_distances,
    double* d_limited_gradients,
    int num_cells,
    LimiterType limiter_type,
    double limiter_zeta,
    int block_size
) {
    if (num_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_cells + block_size - 1) / block_size;
    
    generic_limiter_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_cell_neighbors,
        d_cell_distances,
        d_limited_gradients,
        num_cells,
        limiter_zeta,
        static_cast<int>(limiter_type)
    );
    
    return cudaGetLastError();
}

const char* get_limiter_name(LimiterType limiter_type) {
    switch (limiter_type) {
        case LIMITER_MINMOD:
            return "MinMod";
        case LIMITER_VANLEER:
            return "Van Leer";
        case LIMITER_SUPERBEE:
            return "Superbee";
        case LIMITER_VANALBADA:
            return "Van Albada";
        case LIMITER_VENKATAKRISHNAN:
            return "Venkatakrishnan";
        case LIMITER_BARTH_JESPERSEN:
            return "Barth-Jespersen";
        default:
            return "Unknown";
    }
}
