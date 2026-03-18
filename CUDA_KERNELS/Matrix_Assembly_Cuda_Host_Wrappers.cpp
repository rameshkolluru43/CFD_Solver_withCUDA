// File: Matrix_Assembly_Cuda_Host_Wrappers.cpp
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-01-19
// Description: Host wrapper functions for CUDA matrix assembly kernels
// Author: AI Assistant

#include "Matrix_Assembly_Cuda_Kernels.h"
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

//=============================================================================
// CUDA ERROR CHECKING AND UTILITIES
//=============================================================================

#define CUDA_CHECK(call)                                                                                       \
    do                                                                                                         \
    {                                                                                                          \
        cudaError_t cudaStatus = call;                                                                         \
        if (cudaStatus != cudaSuccess)                                                                         \
        {                                                                                                      \
            fprintf(stderr, "CUDA Error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(cudaStatus)); \
            exit(EXIT_FAILURE);                                                                                \
        }                                                                                                      \
    } while (0)

/**
 * Print CUDA device properties relevant to matrix assembly
 */
void print_cuda_device_info_for_matrix_assembly()
{
    int device;
    CUDA_CHECK(cudaGetDevice(&device));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device));

    std::cout << "\n=== CUDA Device Info for Matrix Assembly ===" << std::endl;
    std::cout << "Device: " << prop.name << std::endl;
    std::cout << "Compute Capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "Global Memory: " << prop.totalGlobalMem / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Shared Memory per Block: " << prop.sharedMemPerBlock / 1024 << " KB" << std::endl;
    std::cout << "Max Threads per Block: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "Max Grid Size: " << prop.maxGridSize[0] << " x " << prop.maxGridSize[1] << std::endl;
    std::cout << "Warp Size: " << prop.warpSize << std::endl;
    std::cout << "Memory Clock Rate: " << prop.memoryClockRate / 1000 << " MHz" << std::endl;
    std::cout << "Memory Bus Width: " << prop.memoryBusWidth << " bits" << std::endl;
    std::cout << "L2 Cache Size: " << prop.l2CacheSize / 1024 << " KB" << std::endl;
    std::cout << "===========================================" << std::endl;
}

/**
 * Estimate memory requirements for matrix assembly
 */
size_t estimate_matrix_assembly_memory_usage(int No_Physical_Cells, bool use_sparse)
{
    size_t total_memory = 0;

    // Input arrays
    total_memory += No_Physical_Cells * sizeof(double);     // cell_areas
    total_memory += No_Physical_Cells * 4 * sizeof(double); // face_areas
    total_memory += No_Physical_Cells * 4 * sizeof(int);    // neighbors
    total_memory += No_Physical_Cells * 4 * sizeof(double); // conservative_vars

    if (use_sparse)
    {
        // Sparse matrix format (COO)
        int avg_neighbors = 3;                                        // Average neighbors per cell in 2D
        int total_nnz = No_Physical_Cells * 16 * (1 + avg_neighbors); // 4x4 entries per cell and neighbors
        total_memory += total_nnz * sizeof(int);                      // row_indices
        total_memory += total_nnz * sizeof(int);                      // col_indices
        total_memory += total_nnz * sizeof(double);                   // values
        total_memory += No_Physical_Cells * sizeof(int);              // nnz_per_cell
        total_memory += No_Physical_Cells * sizeof(int);              // nnz_offsets
    }
    else
    {
        // Dense matrix format
        size_t matrix_size = 4 * No_Physical_Cells;
        total_memory += matrix_size * matrix_size * sizeof(double);
    }

    // Vector b
    total_memory += No_Physical_Cells * 4 * sizeof(double);

    return total_memory;
}

/**
 * Get optimal CUDA launch parameters
 */
std::pair<int, int> get_optimal_launch_params(int No_Physical_Cells, int kernel_type)
{
    int device;
    CUDA_CHECK(cudaGetDevice(&device));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device));

    int block_size;
    int grid_size;

    switch (kernel_type)
    {
    case 0:
        block_size = (256 < prop.maxThreadsPerBlock) ? 256 : prop.maxThreadsPerBlock;
        grid_size = (No_Physical_Cells + block_size - 1) / block_size;
        break;

    case 1:
        block_size = (128 < prop.maxThreadsPerBlock) ? 128 : prop.maxThreadsPerBlock;
        grid_size = (No_Physical_Cells + block_size - 1) / block_size;
        break;

    case 2:
        block_size = (256 < prop.maxThreadsPerBlock) ? 256 : prop.maxThreadsPerBlock;
        block_size = (block_size / prop.warpSize) * prop.warpSize;
        grid_size = (No_Physical_Cells * prop.warpSize + block_size - 1) / block_size;
        break;

    default:
        block_size = 256;
        grid_size = (No_Physical_Cells + block_size - 1) / block_size;
    }

    if (grid_size > prop.maxGridSize[0]) grid_size = prop.maxGridSize[0];

    return std::make_pair(grid_size, block_size);
}

//=============================================================================
// HOST WRAPPER FUNCTIONS
//=============================================================================

/**
 * Host wrapper for dense matrix assembly
 */
double assemble_dense_matrix_cuda(
    const std::vector<double> &cell_areas,
    const std::vector<double> &face_areas,
    const std::vector<int> &neighbors,
    const std::vector<double> &conservative_vars,
    std::vector<double> &matrix,
    double dt,
    int No_Physical_Cells)
{

    cudaEvent_t ev_start, ev_end, ev_k_start, ev_k_end;
    cudaEventCreate(&ev_start); cudaEventCreate(&ev_end);
    cudaEventCreate(&ev_k_start); cudaEventCreate(&ev_k_end);
    cudaEventRecord(ev_start);

    int matrix_size = 4 * No_Physical_Cells;
    matrix.resize(matrix_size * matrix_size, 0.0);

    double *d_cell_areas, *d_face_areas, *d_conservative_vars, *d_matrix;
    int *d_neighbors;

    CUDA_CHECK(cudaMalloc(&d_cell_areas, No_Physical_Cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_face_areas, No_Physical_Cells * 4 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_neighbors, No_Physical_Cells * 4 * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_conservative_vars, No_Physical_Cells * 4 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_matrix, matrix_size * matrix_size * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(d_cell_areas, cell_areas.data(), No_Physical_Cells * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_face_areas, face_areas.data(), No_Physical_Cells * 4 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_neighbors, neighbors.data(), No_Physical_Cells * 4 * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_conservative_vars, conservative_vars.data(), No_Physical_Cells * 4 * sizeof(double), cudaMemcpyHostToDevice));

    auto launch_params = get_optimal_launch_params(matrix_size * matrix_size, 0);
    initialize_matrix_kernel<<<launch_params.first, launch_params.second>>>(d_matrix, matrix_size * matrix_size);
    CUDA_CHECK(cudaGetLastError());

    launch_params = get_optimal_launch_params(No_Physical_Cells, 0);

    cudaEventRecord(ev_k_start);

    assemble_dense_matrix_kernel<<<launch_params.first, launch_params.second>>>(
        d_cell_areas, d_face_areas, d_neighbors, d_conservative_vars,
        d_matrix, dt, No_Physical_Cells);

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());

    cudaEventRecord(ev_k_end);

    CUDA_CHECK(cudaMemcpy(matrix.data(), d_matrix, matrix_size * matrix_size * sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(d_cell_areas));
    CUDA_CHECK(cudaFree(d_face_areas));
    CUDA_CHECK(cudaFree(d_neighbors));
    CUDA_CHECK(cudaFree(d_conservative_vars));
    CUDA_CHECK(cudaFree(d_matrix));

    cudaEventRecord(ev_end);
    cudaEventSynchronize(ev_end);

    float total_time_f = 0.0f, kernel_time_f = 0.0f;
    cudaEventElapsedTime(&total_time_f, ev_start, ev_end);
    cudaEventElapsedTime(&kernel_time_f, ev_k_start, ev_k_end);
    double total_time = total_time_f;
    double kernel_time = kernel_time_f;

    cudaEventDestroy(ev_start); cudaEventDestroy(ev_end);
    cudaEventDestroy(ev_k_start); cudaEventDestroy(ev_k_end);

    std::cout << "Dense Matrix Assembly CUDA Performance:" << std::endl;
    std::cout << "  Total Time: " << total_time << " ms" << std::endl;
    std::cout << "  Kernel Time: " << kernel_time << " ms" << std::endl;
    std::cout << "  Memory Transfer Time: " << (total_time - kernel_time) << " ms" << std::endl;
    std::cout << "  Matrix Size: " << matrix_size << " x " << matrix_size << std::endl;
    std::cout << "  Grid/Block Size: " << launch_params.first << "/" << launch_params.second << std::endl;

    return total_time;
}

/**
 * Host wrapper for sparse matrix assembly
 */
double assemble_sparse_matrix_cuda(
    const std::vector<double> &cell_areas,
    const std::vector<double> &face_areas,
    const std::vector<int> &neighbors,
    const std::vector<double> &conservative_vars,
    std::vector<int> &row_indices,
    std::vector<int> &col_indices,
    std::vector<double> &values,
    double dt,
    int No_Physical_Cells)
{

    cudaEvent_t ev_start, ev_end, ev_k_start, ev_k_end;
    cudaEventCreate(&ev_start); cudaEventCreate(&ev_end);
    cudaEventCreate(&ev_k_start); cudaEventCreate(&ev_k_end);
    cudaEventRecord(ev_start);

    double *d_cell_areas, *d_face_areas, *d_conservative_vars;
    int *d_neighbors;

    CUDA_CHECK(cudaMalloc(&d_cell_areas, No_Physical_Cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_face_areas, No_Physical_Cells * 4 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_neighbors, No_Physical_Cells * 4 * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_conservative_vars, No_Physical_Cells * 4 * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(d_cell_areas, cell_areas.data(), No_Physical_Cells * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_face_areas, face_areas.data(), No_Physical_Cells * 4 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_neighbors, neighbors.data(), No_Physical_Cells * 4 * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_conservative_vars, conservative_vars.data(), No_Physical_Cells * 4 * sizeof(double), cudaMemcpyHostToDevice));

    int *d_nnz_per_cell;
    CUDA_CHECK(cudaMalloc(&d_nnz_per_cell, No_Physical_Cells * sizeof(int)));

    auto launch_params = get_optimal_launch_params(No_Physical_Cells, 1);
    count_nonzeros_per_cell_kernel<<<launch_params.first, launch_params.second>>>(
        d_neighbors, d_nnz_per_cell, No_Physical_Cells);
    CUDA_CHECK(cudaGetLastError());

    // Host-side prefix sum (avoids Thrust dependency)
    std::vector<int> h_nnz_per_cell(No_Physical_Cells);
    CUDA_CHECK(cudaMemcpy(h_nnz_per_cell.data(), d_nnz_per_cell, No_Physical_Cells * sizeof(int), cudaMemcpyDeviceToHost));

    std::vector<int> h_nnz_offsets(No_Physical_Cells + 1);
    h_nnz_offsets[0] = 0;
    for (int i = 0; i < No_Physical_Cells; i++) {
        h_nnz_offsets[i + 1] = h_nnz_offsets[i] + h_nnz_per_cell[i];
    }
    int total_nnz = h_nnz_offsets[No_Physical_Cells];

    int *d_nnz_offsets;
    CUDA_CHECK(cudaMalloc(&d_nnz_offsets, (No_Physical_Cells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_nnz_offsets, h_nnz_offsets.data(), (No_Physical_Cells + 1) * sizeof(int), cudaMemcpyHostToDevice));

    row_indices.resize(total_nnz);
    col_indices.resize(total_nnz);
    values.resize(total_nnz);

    int *d_row_indices, *d_col_indices;
    double *d_values;

    CUDA_CHECK(cudaMalloc(&d_row_indices, total_nnz * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_col_indices, total_nnz * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_values, total_nnz * sizeof(double)));

    cudaEventRecord(ev_k_start);

    assemble_sparse_matrix_kernel<<<launch_params.first, launch_params.second>>>(
        d_cell_areas, d_face_areas, d_neighbors, d_conservative_vars,
        d_nnz_offsets,
        d_row_indices, d_col_indices, d_values,
        dt, No_Physical_Cells);

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());

    cudaEventRecord(ev_k_end);

    CUDA_CHECK(cudaMemcpy(row_indices.data(), d_row_indices, total_nnz * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(col_indices.data(), d_col_indices, total_nnz * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(values.data(), d_values, total_nnz * sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(d_cell_areas));
    CUDA_CHECK(cudaFree(d_face_areas));
    CUDA_CHECK(cudaFree(d_neighbors));
    CUDA_CHECK(cudaFree(d_conservative_vars));
    CUDA_CHECK(cudaFree(d_nnz_per_cell));
    CUDA_CHECK(cudaFree(d_nnz_offsets));
    CUDA_CHECK(cudaFree(d_row_indices));
    CUDA_CHECK(cudaFree(d_col_indices));
    CUDA_CHECK(cudaFree(d_values));

    cudaEventRecord(ev_end);
    cudaEventSynchronize(ev_end);

    float total_time_f = 0.0f, kernel_time_f = 0.0f;
    cudaEventElapsedTime(&total_time_f, ev_start, ev_end);
    cudaEventElapsedTime(&kernel_time_f, ev_k_start, ev_k_end);
    double total_time = total_time_f;
    double kernel_time = kernel_time_f;

    cudaEventDestroy(ev_start); cudaEventDestroy(ev_end);
    cudaEventDestroy(ev_k_start); cudaEventDestroy(ev_k_end);

    std::cout << "Sparse Matrix Assembly CUDA Performance:" << std::endl;
    std::cout << "  Total Time: " << total_time << " ms" << std::endl;
    std::cout << "  Kernel Time: " << kernel_time << " ms" << std::endl;
    std::cout << "  Memory Transfer Time: " << (total_time - kernel_time) << " ms" << std::endl;
    std::cout << "  Total Non-zeros: " << total_nnz << std::endl;
    std::cout << "  Sparsity: " << (1.0 - (double)total_nnz / ((double)4 * No_Physical_Cells * 4 * No_Physical_Cells)) * 100 << "%" << std::endl;
    std::cout << "  Grid/Block Size: " << launch_params.first << "/" << launch_params.second << std::endl;

    return total_time;
}

/**
 * Host wrapper for vector b assembly
 */
double assemble_vector_b_cuda(
    const std::vector<double> &net_flux,
    std::vector<double> &vector_b,
    int No_Physical_Cells)
{

    cudaEvent_t ev_start, ev_end, ev_k_start, ev_k_end;
    cudaEventCreate(&ev_start); cudaEventCreate(&ev_end);
    cudaEventCreate(&ev_k_start); cudaEventCreate(&ev_k_end);
    cudaEventRecord(ev_start);

    int vector_size = No_Physical_Cells * 4;
    vector_b.resize(vector_size);

    double *d_net_flux, *d_vector_b;
    CUDA_CHECK(cudaMalloc(&d_net_flux, vector_size * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_vector_b, vector_size * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_net_flux, net_flux.data(), vector_size * sizeof(double), cudaMemcpyHostToDevice));

    auto launch_params = get_optimal_launch_params(vector_size, 0);

    cudaEventRecord(ev_k_start);

    assemble_vector_b_kernel<<<launch_params.first, launch_params.second>>>(
        d_net_flux, d_vector_b, No_Physical_Cells);

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());

    cudaEventRecord(ev_k_end);

    CUDA_CHECK(cudaMemcpy(vector_b.data(), d_vector_b, vector_size * sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(d_net_flux));
    CUDA_CHECK(cudaFree(d_vector_b));

    cudaEventRecord(ev_end);
    cudaEventSynchronize(ev_end);

    float total_time_f = 0.0f, kernel_time_f = 0.0f;
    cudaEventElapsedTime(&total_time_f, ev_start, ev_end);
    cudaEventElapsedTime(&kernel_time_f, ev_k_start, ev_k_end);

    cudaEventDestroy(ev_start); cudaEventDestroy(ev_end);
    cudaEventDestroy(ev_k_start); cudaEventDestroy(ev_k_end);

    std::cout << "Vector b Assembly CUDA Performance:" << std::endl;
    std::cout << "  Total Time: " << total_time_f << " ms" << std::endl;
    std::cout << "  Kernel Time: " << kernel_time_f << " ms" << std::endl;
    std::cout << "  Vector Size: " << vector_size << std::endl;

    return total_time_f;
}

/**
 * Validate CUDA matrix assembly results
 */
std::vector<double> validate_matrix_cuda(
    const std::vector<double> &matrix,
    int matrix_size)
{

    // Allocate device memory
    double *d_matrix, *d_validation_results;

    CUDA_CHECK(cudaMalloc(&d_matrix, matrix_size * matrix_size * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_validation_results, 4 * sizeof(double)));

    // Initialize validation results
    std::vector<double> init_vals = {0.0, 0.0, 0.0, 0.0};
    CUDA_CHECK(cudaMemcpy(d_validation_results, init_vals.data(), 4 * sizeof(double), cudaMemcpyHostToDevice));

    // Copy matrix to device
    CUDA_CHECK(cudaMemcpy(d_matrix, matrix.data(), matrix_size * matrix_size * sizeof(double), cudaMemcpyHostToDevice));

    // Launch validation kernel
    int block_size = 256;
    int grid_size = (matrix_size * matrix_size + block_size - 1) / block_size;

    validate_matrix_kernel<<<grid_size, block_size>>>(
        d_matrix, matrix_size, d_validation_results);

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());

    // Copy results back
    std::vector<double> results(4);
    CUDA_CHECK(cudaMemcpy(results.data(), d_validation_results, 4 * sizeof(double), cudaMemcpyDeviceToHost));

    // Cleanup
    CUDA_CHECK(cudaFree(d_matrix));
    CUDA_CHECK(cudaFree(d_validation_results));

    return results;
}

//=============================================================================
// BENCHMARKING FUNCTIONS
//=============================================================================

/**
 * Comprehensive performance benchmarking function
 */
std::vector<MatrixAssemblyPerformanceMetrics> benchmark_matrix_assembly_cuda(
    const std::vector<double> &cell_areas,
    const std::vector<double> &face_areas,
    const std::vector<int> &neighbors,
    const std::vector<double> &conservative_vars,
    double dt,
    int No_Physical_Cells,
    int num_iterations)
{

    std::vector<MatrixAssemblyPerformanceMetrics> metrics;

    std::cout << "\n=== Matrix Assembly CUDA Benchmarking ===" << std::endl;
    std::cout << "Number of Physical Cells: " << No_Physical_Cells << std::endl;
    std::cout << "Number of Iterations: " << num_iterations << std::endl;

    // Print device info
    print_cuda_device_info_for_matrix_assembly();

    // Benchmark Dense Matrix Assembly
    {
        std::cout << "\nBenchmarking Dense Matrix Assembly..." << std::endl;

        std::vector<double> times;
        std::vector<double> matrix;

        for (int i = 0; i < num_iterations; i++)
        {
            double time = assemble_dense_matrix_cuda(
                cell_areas, face_areas, neighbors, conservative_vars,
                matrix, dt, No_Physical_Cells);
            times.push_back(time);
        }

        MatrixAssemblyPerformanceMetrics metric;
        metric.kernel_variant = "Dense Matrix Assembly";
        double min_t = times[0]; for (size_t k=1;k<times.size();k++) if(times[k]<min_t) min_t=times[k];
        metric.total_time_ms = min_t;
        metric.memory_usage_bytes = estimate_matrix_assembly_memory_usage(No_Physical_Cells, false);
        metric.throughput_cells_per_sec = No_Physical_Cells / (metric.total_time_ms / 1000.0);

        auto launch_params = get_optimal_launch_params(No_Physical_Cells, 0);
        metric.grid_size = launch_params.first;
        metric.block_size = launch_params.second;

        metrics.push_back(metric);
    }

    // Benchmark Sparse Matrix Assembly
    {
        std::cout << "\nBenchmarking Sparse Matrix Assembly..." << std::endl;

        std::vector<double> times;
        std::vector<int> row_indices, col_indices;
        std::vector<double> values;

        for (int i = 0; i < num_iterations; i++)
        {
            double time = assemble_sparse_matrix_cuda(
                cell_areas, face_areas, neighbors, conservative_vars,
                row_indices, col_indices, values, dt, No_Physical_Cells);
            times.push_back(time);
        }

        MatrixAssemblyPerformanceMetrics metric;
        metric.kernel_variant = "Sparse Matrix Assembly";
        double min_t = times[0]; for (size_t k=1;k<times.size();k++) if(times[k]<min_t) min_t=times[k];
        metric.total_time_ms = min_t;
        metric.memory_usage_bytes = estimate_matrix_assembly_memory_usage(No_Physical_Cells, true);
        metric.throughput_cells_per_sec = No_Physical_Cells / (metric.total_time_ms / 1000.0);

        auto launch_params = get_optimal_launch_params(No_Physical_Cells, 1);
        metric.grid_size = launch_params.first;
        metric.block_size = launch_params.second;

        metrics.push_back(metric);
    }

    // Benchmark Vector b Assembly
    {
        std::cout << "\nBenchmarking Vector b Assembly..." << std::endl;

        std::vector<double> times;
        std::vector<double> vector_b;
        std::vector<double> net_flux(No_Physical_Cells * 4, 1.0); // Dummy data

        for (int i = 0; i < num_iterations; i++)
        {
            double time = assemble_vector_b_cuda(net_flux, vector_b, No_Physical_Cells);
            times.push_back(time);
        }

        MatrixAssemblyPerformanceMetrics metric;
        metric.kernel_variant = "Vector b Assembly";
        double min_t = times[0]; for (size_t k=1;k<times.size();k++) if(times[k]<min_t) min_t=times[k];
        metric.total_time_ms = min_t;
        metric.memory_usage_bytes = No_Physical_Cells * 4 * 2 * sizeof(double); // Input + output
        metric.throughput_cells_per_sec = No_Physical_Cells / (metric.total_time_ms / 1000.0);

        auto launch_params = get_optimal_launch_params(No_Physical_Cells * 4, 0);
        metric.grid_size = launch_params.first;
        metric.block_size = launch_params.second;

        metrics.push_back(metric);
    }

    // Print summary
    std::cout << "\n=== Performance Summary ===" << std::endl;
    for (const auto &metric : metrics)
    {
        std::cout << metric.kernel_variant << ":" << std::endl;
        std::cout << "  Time: " << metric.total_time_ms << " ms" << std::endl;
        std::cout << "  Throughput: " << metric.throughput_cells_per_sec << " cells/sec" << std::endl;
        std::cout << "  Memory Usage: " << metric.memory_usage_bytes / (1024 * 1024) << " MB" << std::endl;
        std::cout << "  Grid/Block: " << metric.grid_size << "/" << metric.block_size << std::endl;
        std::cout << std::endl;
    }

    return metrics;
}