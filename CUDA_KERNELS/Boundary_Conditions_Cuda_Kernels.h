/**
 * @file Boundary_Conditions_Cuda_Kernels.h
 * @brief Header file for CUDA boundary condition kernels
 * @author AI Assistant
 * @date 2026-01-14
 * @version 1.0
 */

#ifndef BOUNDARY_CONDITIONS_CUDA_KERNELS_H
#define BOUNDARY_CONDITIONS_CUDA_KERNELS_H

#include <cuda_runtime.h>

//=============================================================================
// BOUNDARY CONDITION DATA STRUCTURES
//=============================================================================

/**
 * @brief Structure for inlet boundary conditions
 */
struct InletCondition_GPU
{
    double Rho; // Density
    double u;   // x-velocity
    double v;   // y-velocity
    double P;   // Pressure
    double T;   // Temperature
};

/**
 * @brief Structure for exit boundary conditions
 */
struct ExitCondition_GPU
{
    double P; // Prescribed pressure
};

/**
 * @brief Structure for farfield boundary conditions
 */
struct FarfieldCondition_GPU
{
    double Rho_inf;
    double u_inf;
    double v_inf;
    double P_inf;
    double Mach;
};

//=============================================================================
// CUDA KERNEL DECLARATIONS
//=============================================================================

/**
 * @brief Subsonic inlet boundary condition kernel
 */
__global__ void subsonic_inlet_bc_kernel(
    double *U_cells,
    const int *inlet_cell_list,
    const double *face_normals,
    const InletCondition_GPU bc_data,
    int num_inlet_cells,
    double gamma);

/**
 * @brief Supersonic inlet boundary condition kernel
 */
__global__ void supersonic_inlet_bc_kernel(
    double *U_cells,
    const int *inlet_cell_list,
    const InletCondition_GPU bc_data,
    int num_inlet_cells,
    double gamma);

/**
 * @brief Subsonic exit boundary condition kernel
 */
__global__ void subsonic_exit_bc_kernel(
    double *U_cells,
    const int *exit_cell_list,
    const double *face_normals,
    const ExitCondition_GPU bc_data,
    int num_exit_cells,
    double gamma);

/**
 * @brief Supersonic exit boundary condition kernel
 */
__global__ void supersonic_exit_bc_kernel(
    double *U_cells,
    const int *exit_cell_list,
    int num_exit_cells);

/**
 * @brief Viscous wall boundary condition kernel (no-slip)
 */
__global__ void viscous_wall_bc_kernel(
    double *U_cells,
    const int *wall_cell_list,
    int num_wall_cells);

/**
 * @brief Inviscid wall boundary condition kernel (slip)
 */
__global__ void inviscid_wall_bc_kernel(
    double *U_cells,
    const int *wall_cell_list,
    const double *face_normals,
    int num_wall_cells);

/**
 * @brief Symmetry boundary condition kernel
 */
__global__ void symmetry_bc_kernel(
    double *U_cells,
    const int *symmetry_cell_list,
    const double *face_normals,
    int num_symmetry_cells);

/**
 * @brief Farfield boundary condition kernel (characteristic-based)
 */
__global__ void farfield_bc_kernel(
    double *U_cells,
    const int *farfield_cell_list,
    const double *face_normals,
    const FarfieldCondition_GPU bc_data,
    int num_farfield_cells,
    double gamma);

//=============================================================================
// HOST WRAPPER FUNCTIONS
//=============================================================================

/**
 * @brief Host wrapper for launching subsonic inlet BC
 */
cudaError_t launch_subsonic_inlet_bc(
    double *d_U_cells,
    const int *d_inlet_cell_list,
    const double *d_face_normals,
    const InletCondition_GPU &bc_data,
    int num_inlet_cells,
    double gamma,
    int block_size = 256);

/**
 * @brief Host wrapper for launching supersonic inlet BC
 */
cudaError_t launch_supersonic_inlet_bc(
    double *d_U_cells,
    const int *d_inlet_cell_list,
    const InletCondition_GPU &bc_data,
    int num_inlet_cells,
    double gamma,
    int block_size = 256);

/**
 * @brief Host wrapper for launching subsonic exit BC
 */
cudaError_t launch_subsonic_exit_bc(
    double *d_U_cells,
    const int *d_exit_cell_list,
    const double *d_face_normals,
    const ExitCondition_GPU &bc_data,
    int num_exit_cells,
    double gamma,
    int block_size = 256);

/**
 * @brief Host wrapper for launching supersonic exit BC
 */
cudaError_t launch_supersonic_exit_bc(
    double *d_U_cells,
    const int *d_exit_cell_list,
    int num_exit_cells,
    int block_size = 256);

/**
 * @brief Host wrapper for launching viscous wall BC
 */
cudaError_t launch_viscous_wall_bc(
    double *d_U_cells,
    const int *d_wall_cell_list,
    int num_wall_cells,
    int block_size = 256);

/**
 * @brief Host wrapper for launching inviscid wall BC
 */
cudaError_t launch_inviscid_wall_bc(
    double *d_U_cells,
    const int *d_wall_cell_list,
    const double *d_face_normals,
    int num_wall_cells,
    int block_size = 256);

/**
 * @brief Host wrapper for launching symmetry BC
 */
cudaError_t launch_symmetry_bc(
    double *d_U_cells,
    const int *d_symmetry_cell_list,
    const double *d_face_normals,
    int num_symmetry_cells,
    int block_size = 256);

/**
 * @brief Host wrapper for launching farfield BC
 */
cudaError_t launch_farfield_bc(
    double *d_U_cells,
    const int *d_farfield_cell_list,
    const double *d_face_normals,
    const FarfieldCondition_GPU &bc_data,
    int num_farfield_cells,
    double gamma,
    int block_size = 256);

#endif // BOUNDARY_CONDITIONS_CUDA_KERNELS_H
