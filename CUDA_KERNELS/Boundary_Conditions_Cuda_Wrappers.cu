/**
 * @file Boundary_Conditions_Cuda_Wrappers.cu
 * @brief Host wrapper implementations for boundary condition CUDA kernels
 * @author AI Assistant
 * @date 2026-01-14
 */

#include "Boundary_Conditions_Cuda_Kernels.h"
#include <stdio.h>

//=============================================================================
// HOST WRAPPER IMPLEMENTATIONS
//=============================================================================

cudaError_t launch_subsonic_inlet_bc(
    double* d_U_cells,
    const int* d_inlet_cell_list,
    const double* d_face_normals,
    const InletCondition_GPU& bc_data,
    int num_inlet_cells,
    double gamma,
    int block_size
) {
    if (num_inlet_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_inlet_cells + block_size - 1) / block_size;
    
    subsonic_inlet_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_inlet_cell_list,
        d_face_normals,
        bc_data,
        num_inlet_cells,
        gamma
    );
    
    return cudaGetLastError();
}

cudaError_t launch_supersonic_inlet_bc(
    double* d_U_cells,
    const int* d_inlet_cell_list,
    const InletCondition_GPU& bc_data,
    int num_inlet_cells,
    double gamma,
    int block_size
) {
    if (num_inlet_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_inlet_cells + block_size - 1) / block_size;
    
    supersonic_inlet_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_inlet_cell_list,
        bc_data,
        num_inlet_cells,
        gamma
    );
    
    return cudaGetLastError();
}

cudaError_t launch_subsonic_exit_bc(
    double* d_U_cells,
    const int* d_exit_cell_list,
    const double* d_face_normals,
    const ExitCondition_GPU& bc_data,
    int num_exit_cells,
    double gamma,
    int block_size
) {
    if (num_exit_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_exit_cells + block_size - 1) / block_size;
    
    subsonic_exit_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_exit_cell_list,
        d_face_normals,
        bc_data,
        num_exit_cells,
        gamma
    );
    
    return cudaGetLastError();
}

cudaError_t launch_supersonic_exit_bc(
    double* d_U_cells,
    const int* d_exit_cell_list,
    int num_exit_cells,
    int block_size
) {
    if (num_exit_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_exit_cells + block_size - 1) / block_size;
    
    supersonic_exit_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_exit_cell_list,
        num_exit_cells
    );
    
    return cudaGetLastError();
}

cudaError_t launch_viscous_wall_bc(
    double* d_U_cells,
    const int* d_wall_cell_list,
    int num_wall_cells,
    int block_size
) {
    if (num_wall_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_wall_cells + block_size - 1) / block_size;
    
    viscous_wall_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_wall_cell_list,
        num_wall_cells
    );
    
    return cudaGetLastError();
}

cudaError_t launch_inviscid_wall_bc(
    double* d_U_cells,
    const int* d_wall_cell_list,
    const double* d_face_normals,
    int num_wall_cells,
    int block_size
) {
    if (num_wall_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_wall_cells + block_size - 1) / block_size;
    
    inviscid_wall_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_wall_cell_list,
        d_face_normals,
        num_wall_cells
    );
    
    return cudaGetLastError();
}

cudaError_t launch_symmetry_bc(
    double* d_U_cells,
    const int* d_symmetry_cell_list,
    const double* d_face_normals,
    int num_symmetry_cells,
    int block_size
) {
    if (num_symmetry_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_symmetry_cells + block_size - 1) / block_size;
    
    symmetry_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_symmetry_cell_list,
        d_face_normals,
        num_symmetry_cells
    );
    
    return cudaGetLastError();
}

cudaError_t launch_farfield_bc(
    double* d_U_cells,
    const int* d_farfield_cell_list,
    const double* d_face_normals,
    const FarfieldCondition_GPU& bc_data,
    int num_farfield_cells,
    double gamma,
    int block_size
) {
    if (num_farfield_cells == 0) return cudaSuccess;
    
    int num_blocks = (num_farfield_cells + block_size - 1) / block_size;
    
    farfield_bc_kernel<<<num_blocks, block_size>>>(
        d_U_cells,
        d_farfield_cell_list,
        d_face_normals,
        bc_data,
        num_farfield_cells,
        gamma
    );
    
    return cudaGetLastError();
}
