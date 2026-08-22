// File: Matrix_Assembly_CUDA_Integration_Example.cpp
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-01-19
// Description: Example integration of CUDA matrix assembly kernels with existing CPU code
// Author: AI Assistant

#include "Matrix_Assembly_Cuda_Kernels.h"
#include "definitions.h"
#include "Globals.h"
#include "Assemble_Matrix.h"
#include <chrono>
#include <iostream>
#include <algorithm>

//=============================================================================
// CUDA-ACCELERATED MATRIX ASSEMBLY FUNCTIONS
//=============================================================================

/**
 * CUDA-accelerated version of Assemble_A function
 * Replaces the CPU version with GPU acceleration
 */
vector<V_D> Assemble_A_CUDA(vector<V_D> &A, double &dt)
{
    // Input validation
    if (dt <= 0.0)
    {
        std::cout << "Error: Invalid time step dt = " << dt << std::endl;
        return A;
    }
    if (No_Physical_Cells <= 0)
    {
        std::cout << "Error: Invalid number of physical cells = " << No_Physical_Cells << std::endl;
        return A;
    }

    std::cout << "Using CUDA-accelerated matrix assembly..." << std::endl;

    // Prepare input data for CUDA kernels
    std::vector<double> cell_areas(No_Physical_Cells);
    std::vector<double> face_areas(No_Physical_Cells * 4);
    std::vector<int> neighbors(No_Physical_Cells * 4);
    std::vector<double> conservative_vars(No_Physical_Cells * 4);

    // Extract data from global Cells structure
    for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
    {
        // Cell areas
        cell_areas[cell_idx] = Cells[cell_idx].Area;

        // Face areas (4 per cell: left, bottom, right, top)
        for (int face = 0; face < 4; face++)
        {
            face_areas[cell_idx * 4 + face] = Cells[cell_idx].Face_Areas[face];
        }

        // Neighbors (4 per cell: left, bottom, right, top)
        for (int face = 0; face < 4; face++)
        {
            neighbors[cell_idx * 4 + face] = Cells[cell_idx].Neighbours[face];
        }

        // Conservative variables (rho, rho*u, rho*v, E)
        for (int var = 0; var < 4; var++)
        {
            conservative_vars[cell_idx * 4 + var] = U_Cells[cell_idx][var];
        }
    }

    // Prepare output matrix
    int matrix_size = 4 * No_Physical_Cells;
    A.resize(matrix_size, V_D(matrix_size, 0.0));
    std::vector<double> matrix_flat(matrix_size * matrix_size, 0.0);

    // Call CUDA matrix assembly
    auto start_time = std::chrono::high_resolution_clock::now();

    double cuda_time = assemble_dense_matrix_cuda(
        cell_areas, face_areas, neighbors, conservative_vars,
        matrix_flat, dt, No_Physical_Cells);

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    // Convert flat matrix back to 2D format
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            A[i][j] = matrix_flat[i * matrix_size + j];
        }
    }

    std::cout << "CUDA matrix assembly completed in " << total_time << " ms" << std::endl;

    return A;
}

/**
 * CUDA-accelerated version of Assemble_A1 function for sparse matrix assembly
 * Replaces the CPU version with GPU acceleration and populates global sparse matrix arrays
 */
void Assemble_A1_CUDA(double &dt)
{
    // Input validation
    if (dt <= 0.0)
    {
        std::cout << "Error: Invalid time step dt = " << dt << std::endl;
        return;
    }
    if (No_Physical_Cells <= 0)
    {
        std::cout << "Error: Invalid number of physical cells = " << No_Physical_Cells << std::endl;
        return;
    }

    std::cout << "Using CUDA-accelerated sparse matrix assembly..." << std::endl;

    // Clear existing global arrays
    row_indices.clear();
    col_indices.clear();
    values.clear();

    // Prepare input data for CUDA kernels
    std::vector<double> cell_areas(No_Physical_Cells);
    std::vector<double> face_areas(No_Physical_Cells * 4);
    std::vector<int> neighbors(No_Physical_Cells * 4);
    std::vector<double> conservative_vars(No_Physical_Cells * 4);

    // Extract data from global Cells structure
    for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
    {
        // Bounds checking
        if (cell_idx >= (int)Cells.size())
        {
            std::cout << "Error: Cell index " << cell_idx << " out of bounds" << std::endl;
            continue;
        }

        cell_areas[cell_idx] = Cells[cell_idx].Area;

        for (int face = 0; face < 4; face++)
        {
            face_areas[cell_idx * 4 + face] = Cells[cell_idx].Face_Areas[face];
            neighbors[cell_idx * 4 + face] = Cells[cell_idx].Neighbours[face];
        }

        for (int var = 0; var < 4; var++)
        {
            conservative_vars[cell_idx * 4 + var] = U_Cells[cell_idx][var];
        }
    }

    // Call CUDA sparse matrix assembly
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<int> cuda_row_indices, cuda_col_indices;
    std::vector<double> cuda_values;

    double cuda_time = assemble_sparse_matrix_cuda(
        cell_areas, face_areas, neighbors, conservative_vars,
        cuda_row_indices, cuda_col_indices, cuda_values,
        dt, No_Physical_Cells);

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    // Copy results to global arrays
    row_indices = std::move(cuda_row_indices);
    col_indices = std::move(cuda_col_indices);
    values = std::move(cuda_values);

    std::cout << "CUDA sparse matrix assembly completed in " << total_time << " ms" << std::endl;
    std::cout << "Generated " << values.size() << " non-zero entries" << std::endl;
}

/**
 * CUDA-accelerated version of Assemble_b function
 * Replaces the CPU version with GPU acceleration
 */
V_D Assemble_b_CUDA(V_D &b)
{
    // Input validation
    if (No_Physical_Cells <= 0)
    {
        std::cout << "Error: Invalid number of physical cells = " << No_Physical_Cells << std::endl;
        return b;
    }

    std::cout << "Using CUDA-accelerated vector b assembly..." << std::endl;

    // Prepare input data - flatten Cells_Net_Flux
    std::vector<double> net_flux_flat(No_Physical_Cells * 4);
    for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
    {
        for (int var = 0; var < 4; var++)
        {
            net_flux_flat[cell_idx * 4 + var] = Cells_Net_Flux[cell_idx][var];
        }
    }

    // Call CUDA vector assembly
    std::vector<double> vector_b_flat;
    double cuda_time = assemble_vector_b_cuda(net_flux_flat, vector_b_flat, No_Physical_Cells);

    // Convert back to the expected format
    b.resize(No_Physical_Cells * 4);
    for (int i = 0; i < No_Physical_Cells * 4; i++)
    {
        b[i] = vector_b_flat[i];
    }

    std::cout << "CUDA vector b assembly completed in " << cuda_time << " ms" << std::endl;

    return b;
}

//=============================================================================
// COMPARATIVE BENCHMARKING FUNCTIONS
//=============================================================================

/**
 * Function to compare CPU vs CUDA performance for matrix assembly
 */
void benchmark_matrix_assembly_cpu_vs_cuda(double dt, int num_iterations = 5)
{
    std::cout << "\n=== CPU vs CUDA Matrix Assembly Benchmark ===" << std::endl;
    std::cout << "Number of Physical Cells: " << No_Physical_Cells << std::endl;
    std::cout << "Time Step: " << dt << std::endl;
    std::cout << "Iterations: " << num_iterations << std::endl;

    // Prepare test data
    std::vector<double> cpu_times, cuda_times;

    // CPU Benchmark - Dense Matrix Assembly
    std::cout << "\nBenchmarking CPU Dense Matrix Assembly..." << std::endl;
    for (int i = 0; i < num_iterations; i++)
    {
        vector<V_D> A_cpu(4 * No_Physical_Cells, V_D(4 * No_Physical_Cells, 0.0));

        auto start = std::chrono::high_resolution_clock::now();
        A_cpu = Assemble_A(A_cpu, dt);
        auto end = std::chrono::high_resolution_clock::now();

        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cpu_times.push_back(time_ms);
        std::cout << "  Iteration " << (i + 1) << ": " << time_ms << " ms" << std::endl;
    }

    // CUDA Benchmark - Dense Matrix Assembly
    std::cout << "\nBenchmarking CUDA Dense Matrix Assembly..." << std::endl;
    for (int i = 0; i < num_iterations; i++)
    {
        vector<V_D> A_cuda(4 * No_Physical_Cells, V_D(4 * No_Physical_Cells, 0.0));

        auto start = std::chrono::high_resolution_clock::now();
        A_cuda = Assemble_A_CUDA(A_cuda, dt);
        auto end = std::chrono::high_resolution_clock::now();

        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cuda_times.push_back(time_ms);
        std::cout << "  Iteration " << (i + 1) << ": " << time_ms << " ms" << std::endl;
    }

    // Calculate statistics
    double cpu_avg = std::accumulate(cpu_times.begin(), cpu_times.end(), 0.0) / cpu_times.size();
    double cuda_avg = std::accumulate(cuda_times.begin(), cuda_times.end(), 0.0) / cuda_times.size();
    double cpu_min = *std::min_element(cpu_times.begin(), cpu_times.end());
    double cuda_min = *std::min_element(cuda_times.begin(), cuda_times.end());

    double speedup_avg = cpu_avg / cuda_avg;
    double speedup_best = cpu_min / cuda_min;

    // Print results
    std::cout << "\n=== Performance Results ===" << std::endl;
    std::cout << "CPU Average Time: " << cpu_avg << " ms" << std::endl;
    std::cout << "CUDA Average Time: " << cuda_avg << " ms" << std::endl;
    std::cout << "CPU Best Time: " << cpu_min << " ms" << std::endl;
    std::cout << "CUDA Best Time: " << cuda_min << " ms" << std::endl;
    std::cout << "Average Speedup: " << speedup_avg << "x" << std::endl;
    std::cout << "Best Speedup: " << speedup_best << "x" << std::endl;

    if (speedup_avg > 1.0)
    {
        std::cout << "✅ CUDA is " << speedup_avg << "x faster than CPU on average" << std::endl;
    }
    else
    {
        std::cout << "⚠️  CPU is " << (1.0 / speedup_avg) << "x faster than CUDA on average" << std::endl;
        std::cout << "    This may be due to small problem size or memory transfer overhead" << std::endl;
    }

    // Memory usage analysis
    size_t dense_memory = estimate_matrix_assembly_memory_usage(No_Physical_Cells, false);
    size_t sparse_memory = estimate_matrix_assembly_memory_usage(No_Physical_Cells, true);

    std::cout << "\n=== Memory Usage Analysis ===" << std::endl;
    std::cout << "Dense Matrix Memory: " << dense_memory / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Sparse Matrix Memory: " << sparse_memory / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Memory Savings with Sparse: " << (1.0 - (double)sparse_memory / dense_memory) * 100 << "%" << std::endl;
}

/**
 * Function to benchmark sparse matrix assembly
 */
void benchmark_sparse_matrix_assembly(double dt, int num_iterations = 5)
{
    std::cout << "\n=== Sparse Matrix Assembly Benchmark ===" << std::endl;

    std::vector<double> cpu_times, cuda_times;

    // CPU Sparse Matrix Assembly Benchmark
    std::cout << "\nBenchmarking CPU Sparse Matrix Assembly..." << std::endl;
    for (int i = 0; i < num_iterations; i++)
    {
        // Clear global arrays
        row_indices.clear();
        col_indices.clear();
        values.clear();

        auto start = std::chrono::high_resolution_clock::now();
        Assemble_A1(dt);
        auto end = std::chrono::high_resolution_clock::now();

        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cpu_times.push_back(time_ms);
        std::cout << "  Iteration " << (i + 1) << ": " << time_ms << " ms, "
                  << values.size() << " nnz" << std::endl;
    }

    // CUDA Sparse Matrix Assembly Benchmark
    std::cout << "\nBenchmarking CUDA Sparse Matrix Assembly..." << std::endl;
    for (int i = 0; i < num_iterations; i++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        Assemble_A1_CUDA(dt);
        auto end = std::chrono::high_resolution_clock::now();

        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cuda_times.push_back(time_ms);
        std::cout << "  Iteration " << (i + 1) << ": " << time_ms << " ms, "
                  << values.size() << " nnz" << std::endl;
    }

    // Results
    double cpu_avg = std::accumulate(cpu_times.begin(), cpu_times.end(), 0.0) / cpu_times.size();
    double cuda_avg = std::accumulate(cuda_times.begin(), cuda_times.end(), 0.0) / cuda_times.size();
    double speedup = cpu_avg / cuda_avg;

    std::cout << "\n=== Sparse Matrix Results ===" << std::endl;
    std::cout << "CPU Average: " << cpu_avg << " ms" << std::endl;
    std::cout << "CUDA Average: " << cuda_avg << " ms" << std::endl;
    std::cout << "Speedup: " << speedup << "x" << std::endl;
}

//=============================================================================
// INTEGRATION TEST FUNCTIONS
//=============================================================================

/**
 * Function to validate that CUDA and CPU results are consistent
 */
bool validate_cuda_cpu_consistency(double dt, double tolerance = 1e-10)
{
    std::cout << "\n=== Validating CUDA vs CPU Consistency ===" << std::endl;

    // Test Dense Matrix Assembly
    vector<V_D> A_cpu(4 * No_Physical_Cells, V_D(4 * No_Physical_Cells, 0.0));
    vector<V_D> A_cuda(4 * No_Physical_Cells, V_D(4 * No_Physical_Cells, 0.0));

    A_cpu = Assemble_A(A_cpu, dt);
    A_cuda = Assemble_A_CUDA(A_cuda, dt);

    // Compare matrices
    double max_diff = 0.0;
    int diff_count = 0;

    for (int i = 0; i < 4 * No_Physical_Cells; i++)
    {
        for (int j = 0; j < 4 * No_Physical_Cells; j++)
        {
            double diff = std::abs(A_cpu[i][j] - A_cuda[i][j]);
            if (diff > tolerance)
            {
                diff_count++;
                max_diff = std::max(max_diff, diff);
            }
        }
    }

    std::cout << "Dense Matrix Comparison:" << std::endl;
    std::cout << "  Max Difference: " << max_diff << std::endl;
    std::cout << "  Entries with diff > " << tolerance << ": " << diff_count << std::endl;

    bool dense_consistent = (max_diff < tolerance * 10); // Allow some numerical error

    // Test Vector b Assembly
    V_D b_cpu, b_cuda;
    b_cpu = Assemble_b(b_cpu);
    b_cuda = Assemble_b_CUDA(b_cuda);

    double max_b_diff = 0.0;
    for (int i = 0; i < (int)b_cpu.size() && i < (int)b_cuda.size(); i++)
    {
        double diff = std::abs(b_cpu[i] - b_cuda[i]);
        max_b_diff = std::max(max_b_diff, diff);
    }

    std::cout << "Vector b Comparison:" << std::endl;
    std::cout << "  Max Difference: " << max_b_diff << std::endl;

    bool vector_consistent = (max_b_diff < tolerance * 10);

    bool overall_consistent = dense_consistent && vector_consistent;

    if (overall_consistent)
    {
        std::cout << "✅ CUDA and CPU results are consistent!" << std::endl;
    }
    else
    {
        std::cout << "❌ CUDA and CPU results differ significantly!" << std::endl;
        if (!dense_consistent)
            std::cout << "   - Dense matrix assembly inconsistent" << std::endl;
        if (!vector_consistent)
            std::cout << "   - Vector b assembly inconsistent" << std::endl;
    }

    return overall_consistent;
}

//=============================================================================
// MAIN DEMONSTRATION FUNCTION
//=============================================================================

/**
 * Main function to demonstrate CUDA matrix assembly integration
 */
void demonstrate_cuda_matrix_assembly()
{
    std::cout << "\n=== CUDA Matrix Assembly Integration Demo ===" << std::endl;

    // Check if we have valid data
    if (No_Physical_Cells <= 0)
    {
        std::cout << "Error: No physical cells available. Initialize grid first." << std::endl;
        return;
    }

    if (Cells.empty())
    {
        std::cout << "Error: Cells array is empty. Initialize grid first." << std::endl;
        return;
    }

    double dt = 1e-4; // Example time step

    // Print system info
    print_cuda_device_info_for_matrix_assembly();

    // 1. Validate consistency
    if (validate_cuda_cpu_consistency(dt))
    {
        std::cout << "✅ Validation passed - proceeding with benchmarks" << std::endl;
    }
    else
    {
        std::cout << "⚠️  Validation failed - results may not be reliable" << std::endl;
    }

    // 2. Performance benchmarking
    benchmark_matrix_assembly_cpu_vs_cuda(dt, 3);
    benchmark_sparse_matrix_assembly(dt, 3);

    // 3. Memory usage analysis
    std::cout << "\n=== Memory Usage Analysis ===" << std::endl;
    size_t dense_mem = estimate_matrix_assembly_memory_usage(No_Physical_Cells, false);
    size_t sparse_mem = estimate_matrix_assembly_memory_usage(No_Physical_Cells, true);

    std::cout << "Problem Size: " << No_Physical_Cells << " cells" << std::endl;
    std::cout << "Dense Matrix Memory: " << dense_mem / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Sparse Matrix Memory: " << sparse_mem / (1024 * 1024) << " MB" << std::endl;

    // Get device memory info
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "GPU Memory Available: " << free_mem / (1024 * 1024) << " MB" << std::endl;
    std::cout << "GPU Memory Total: " << total_mem / (1024 * 1024) << " MB" << std::endl;

    if (dense_mem > free_mem)
    {
        std::cout << "⚠️  Dense matrix exceeds available GPU memory - use sparse format" << std::endl;
    }
    else
    {
        std::cout << "✅ Dense matrix fits in GPU memory" << std::endl;
    }

    std::cout << "\n=== Integration Demo Complete ===" << std::endl;
}