// File: Matrix_Assembly_Cuda_Kernels.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-01-19
// Description: Header file for CUDA kernels for matrix assembly operations
// Author: AI Assistant

#ifndef MATRIX_ASSEMBLY_CUDA_KERNELS_H
#define MATRIX_ASSEMBLY_CUDA_KERNELS_H

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <utility>

//=============================================================================
// CUDA KERNEL DECLARATIONS
//=============================================================================

/**
 * CUDA kernel for assembling the dense matrix A in implicit time integration
 * @param d_cell_areas Device array of cell areas [No_Physical_Cells]
 * @param d_face_areas Device array of face areas [No_Physical_Cells * 4]
 * @param d_neighbors Device array of neighbor indices [No_Physical_Cells * 4]
 * @param d_conservative_vars Device array of conservative variables [No_Physical_Cells * 4]
 * @param d_matrix Device array for output matrix [4*No_Physical_Cells * 4*No_Physical_Cells]
 * @param dt Time step
 * @param No_Physical_Cells Number of physical cells
 */
__global__ void assemble_dense_matrix_kernel(
    const double* d_cell_areas,
    const double* d_face_areas,
    const int* d_neighbors,
    const double* d_conservative_vars,
    double* d_matrix,
    double dt,
    int No_Physical_Cells);

/**
 * CUDA kernel for counting non-zero entries per cell for sparse matrix assembly
 * @param d_neighbors Device array of neighbor indices [No_Physical_Cells * 4]
 * @param d_nnz_per_cell Device array for output nnz count per cell [No_Physical_Cells]
 * @param No_Physical_Cells Number of physical cells
 */
__global__ void count_nonzeros_per_cell_kernel(
    const int* d_neighbors,
    int* d_nnz_per_cell,
    int No_Physical_Cells);

/**
 * CUDA kernel for sparse matrix assembly in COO format
 * @param d_cell_areas Device array of cell areas [No_Physical_Cells]
 * @param d_face_areas Device array of face areas [No_Physical_Cells * 4]
 * @param d_neighbors Device array of neighbor indices [No_Physical_Cells * 4]
 * @param d_conservative_vars Device array of conservative variables [No_Physical_Cells * 4]
 * @param d_nnz_offsets Device array of cumulative nnz offsets [No_Physical_Cells]
 * @param d_row_indices Device array for output row indices [total_nnz]
 * @param d_col_indices Device array for output column indices [total_nnz]
 * @param d_values Device array for output matrix values [total_nnz]
 * @param dt Time step
 * @param No_Physical_Cells Number of physical cells
 */
__global__ void assemble_sparse_matrix_kernel(
    const double* d_cell_areas,
    const double* d_face_areas,
    const int* d_neighbors,
    const double* d_conservative_vars,
    const int* d_nnz_offsets,
    int* d_row_indices,
    int* d_col_indices,
    double* d_values,
    double dt,
    int No_Physical_Cells);

/**
 * Memory-coalesced version of sparse matrix assembly
 */
__global__ void assemble_sparse_matrix_coalesced_kernel(
    const double* d_cell_areas,
    const double* d_face_areas,
    const int* d_neighbors,
    const double* d_conservative_vars,
    const int* d_nnz_offsets,
    int* d_row_indices,
    int* d_col_indices,
    double* d_values,
    double dt,
    int No_Physical_Cells);

/**
 * CUDA kernel for assembling the RHS vector b = -R (negative residual)
 * @param d_net_flux Device array of net flux per cell [No_Physical_Cells * 4]
 * @param d_vector_b Device array for output RHS vector [No_Physical_Cells * 4]
 * @param No_Physical_Cells Number of physical cells
 */
__global__ void assemble_vector_b_kernel(
    const double* d_net_flux,
    double* d_vector_b,
    int No_Physical_Cells);

/**
 * Kernel to initialize matrix with zeros
 * @param matrix Device matrix to initialize
 * @param size Total number of elements in matrix
 */
__global__ void initialize_matrix_kernel(double* matrix, int size);

/**
 * Kernel to validate matrix assembly results
 * @param d_matrix Device matrix to validate
 * @param matrix_size Size of square matrix (number of rows/columns)
 * @param d_validation_results Device array for validation results [4]
 */
__global__ void validate_matrix_kernel(
    const double* d_matrix,
    int matrix_size,
    double* d_validation_results);

//=============================================================================
// HOST WRAPPER FUNCTIONS
//=============================================================================

/**
 * Host wrapper for dense matrix assembly
 * @param cell_areas Host array of cell areas
 * @param face_areas Host array of face areas (4 per cell)
 * @param neighbors Host array of neighbor indices (4 per cell)
 * @param conservative_vars Host array of conservative variables (4 per cell)
 * @param matrix Host output matrix (flattened)
 * @param dt Time step
 * @param No_Physical_Cells Number of physical cells
 * @return Execution time in milliseconds
 */
double assemble_dense_matrix_cuda(
    const std::vector<double>& cell_areas,
    const std::vector<double>& face_areas,
    const std::vector<int>& neighbors,
    const std::vector<double>& conservative_vars,
    std::vector<double>& matrix,
    double dt,
    int No_Physical_Cells);

/**
 * Host wrapper for sparse matrix assembly
 * @param cell_areas Host array of cell areas
 * @param face_areas Host array of face areas (4 per cell)
 * @param neighbors Host array of neighbor indices (4 per cell)
 * @param conservative_vars Host array of conservative variables (4 per cell)
 * @param row_indices Host output array for row indices
 * @param col_indices Host output array for column indices
 * @param values Host output array for matrix values
 * @param dt Time step
 * @param No_Physical_Cells Number of physical cells
 * @return Execution time in milliseconds
 */
double assemble_sparse_matrix_cuda(
    const std::vector<double>& cell_areas,
    const std::vector<double>& face_areas,
    const std::vector<int>& neighbors,
    const std::vector<double>& conservative_vars,
    std::vector<int>& row_indices,
    std::vector<int>& col_indices,
    std::vector<double>& values,
    double dt,
    int No_Physical_Cells);

/**
 * Host wrapper for vector b assembly
 * @param net_flux Host array of net flux per cell
 * @param vector_b Host output RHS vector
 * @param No_Physical_Cells Number of physical cells
 * @return Execution time in milliseconds
 */
double assemble_vector_b_cuda(
    const std::vector<double>& net_flux,
    std::vector<double>& vector_b,
    int No_Physical_Cells);

/**
 * Utility function to validate CUDA matrix assembly results
 * @param matrix Host matrix to validate
 * @param matrix_size Size of square matrix
 * @return Validation statistics (min, max, nnz, diagonal_sum)
 */
std::vector<double> validate_matrix_cuda(
    const std::vector<double>& matrix,
    int matrix_size);

//=============================================================================
// BENCHMARKING AND PERFORMANCE ANALYSIS
//=============================================================================

/**
 * Structure to hold performance metrics
 */
struct MatrixAssemblyPerformanceMetrics {
    double kernel_time_ms;
    double memory_transfer_time_ms;
    double total_time_ms;
    size_t memory_usage_bytes;
    double throughput_cells_per_sec;
    double speedup_vs_cpu;
    int grid_size;
    int block_size;
    std::string kernel_variant;
};

/**
 * Comprehensive performance benchmarking function
 * @param cell_areas Cell areas array
 * @param face_areas Face areas array
 * @param neighbors Neighbors array
 * @param conservative_vars Conservative variables array
 * @param dt Time step
 * @param No_Physical_Cells Number of physical cells
 * @param num_iterations Number of benchmark iterations
 * @return Performance metrics for different kernel variants
 */
std::vector<MatrixAssemblyPerformanceMetrics> benchmark_matrix_assembly_cuda(
    const std::vector<double>& cell_areas,
    const std::vector<double>& face_areas,
    const std::vector<int>& neighbors,
    const std::vector<double>& conservative_vars,
    double dt,
    int No_Physical_Cells,
    int num_iterations = 10);

/**
 * Function to get optimal CUDA launch parameters
 * @param No_Physical_Cells Number of physical cells
 * @param kernel_type Type of kernel (0=dense, 1=sparse, 2=coalesced)
 * @return Pair of (grid_size, block_size)
 */
std::pair<int, int> get_optimal_launch_params(int No_Physical_Cells, int kernel_type);

//=============================================================================
// ERROR CHECKING AND DEBUGGING
//=============================================================================

/**
 * CUDA error checking macro - expanded version for debugging
 */
#define CUDA_CHECK_DEBUG(call, file, line) \
    do { \
        cudaError_t cudaStatus = call; \
        if (cudaStatus != cudaSuccess) { \
            fprintf(stderr, "CUDA Error at %s:%d in %s:%d - %s\n", \
                    file, line, __FILE__, __LINE__, cudaGetErrorString(cudaStatus)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

/**
 * Function to print CUDA device properties relevant to matrix assembly
 */
void print_cuda_device_info_for_matrix_assembly();

/**
 * Function to estimate memory requirements for matrix assembly
 * @param No_Physical_Cells Number of physical cells
 * @param use_sparse True for sparse matrix, false for dense
 * @return Estimated memory usage in bytes
 */
size_t estimate_matrix_assembly_memory_usage(int No_Physical_Cells, bool use_sparse);

#endif // MATRIX_ASSEMBLY_CUDA_KERNELS_H