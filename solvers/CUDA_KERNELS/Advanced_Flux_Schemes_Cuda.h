/**
 * @file Advanced_Flux_Schemes_Cuda.h
 * @brief Header for Roe, HLLC, and LLF flux schemes
 * @date 2026-01-15
 */

#ifndef ADVANCED_FLUX_SCHEMES_CUDA_H
#define ADVANCED_FLUX_SCHEMES_CUDA_H

#include <cuda_runtime.h>

//=============================================================================
// KERNEL DECLARATIONS - ROE FLUX
//=============================================================================

__global__ void roe_flux_kernel(
    const double *U_cells,
    const double *P_cells,
    const int *face_neighbors,
    const double *face_normals,
    const double *face_areas,
    double *flux_output,
    int num_cells,
    double gamma);

__global__ void roe_flux_2nd_order_kernel(
    const double *U_cells,
    const double *P_cells,
    const int *face_neighbors,
    const double *face_normals,
    const double *face_areas,
    const double *limited_gradients,
    const double *cell_centers,
    const double *face_centers,
    double *flux_output,
    int num_cells,
    double gamma);

//=============================================================================
// KERNEL DECLARATIONS - HLLC FLUX
//=============================================================================

__global__ void hllc_flux_kernel(
    const double *U_cells,
    const double *P_cells,
    const int *face_neighbors,
    const double *face_normals,
    const double *face_areas,
    double *flux_output,
    int num_cells,
    double gamma);

//=============================================================================
// KERNEL DECLARATIONS - LLF FLUX
//=============================================================================

__global__ void llf_flux_kernel(
    const double *U_cells,
    const double *P_cells,
    const int *face_neighbors,
    const double *face_normals,
    const double *face_areas,
    double *flux_output,
    int num_cells,
    double gamma);

//=============================================================================
// HOST WRAPPER FUNCTIONS
//=============================================================================

cudaError_t launch_roe_flux(
    const double *d_U_cells,
    const double *d_P_cells,
    const int *d_face_neighbors,
    const double *d_face_normals,
    const double *d_face_areas,
    double *d_flux_output,
    int num_cells,
    double gamma,
    int block_size = 256);

cudaError_t launch_roe_flux_2nd_order(
    const double *d_U_cells,
    const double *d_P_cells,
    const int *d_face_neighbors,
    const double *d_face_normals,
    const double *d_face_areas,
    const double *d_limited_gradients,
    const double *d_cell_centers,
    const double *d_face_centers,
    double *d_flux_output,
    int num_cells,
    double gamma,
    int block_size = 256);

cudaError_t launch_hllc_flux(
    const double *d_U_cells,
    const double *d_P_cells,
    const int *d_face_neighbors,
    const double *d_face_normals,
    const double *d_face_areas,
    double *d_flux_output,
    int num_cells,
    double gamma,
    int block_size = 256);

cudaError_t launch_llf_flux(
    const double *d_U_cells,
    const double *d_P_cells,
    const int *d_face_neighbors,
    const double *d_face_normals,
    const double *d_face_areas,
    double *d_flux_output,
    int num_cells,
    double gamma,
    int block_size = 256);

#endif // ADVANCED_FLUX_SCHEMES_CUDA_H
