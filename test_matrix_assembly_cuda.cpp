// File: test_matrix_assembly_cuda.cpp
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-01-19
// Description: Test program for CUDA matrix assembly kernels
// Author: AI Assistant

#include "Matrix_Assembly_Cuda_Kernels.h"
#include "definitions.h"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>

//=============================================================================
// TEST DATA GENERATION
//=============================================================================

/**
 * Generate synthetic test data for matrix assembly
 */
struct TestData
{
    std::vector<double> cell_areas;
    std::vector<double> face_areas;
    std::vector<int> neighbors;
    std::vector<double> conservative_vars;
    int No_Physical_Cells;
    double dt;
};

TestData generate_test_data(int num_cells, bool structured_grid = true)
{
    TestData data;
    data.No_Physical_Cells = num_cells;
    data.dt = 1e-4;

    std::random_device rd;
    std::mt19937 gen(42); // Fixed seed for reproducibility
    std::uniform_real_distribution<> area_dist(0.001, 0.01);
    std::uniform_real_distribution<> face_dist(0.01, 0.1);
    std::uniform_real_distribution<> rho_dist(0.5, 2.0);
    std::uniform_real_distribution<> vel_dist(-100.0, 100.0);
    std::uniform_real_distribution<> energy_dist(1e5, 3e5);

    // Generate cell areas
    data.cell_areas.resize(num_cells);
    for (int i = 0; i < num_cells; i++)
    {
        data.cell_areas[i] = area_dist(gen);
    }

    // Generate face areas (4 per cell)
    data.face_areas.resize(num_cells * 4);
    for (int i = 0; i < num_cells * 4; i++)
    {
        data.face_areas[i] = face_dist(gen);
    }

    // Generate neighbors (structured 2D grid pattern)
    data.neighbors.resize(num_cells * 4);
    if (structured_grid)
    {
        int nx = static_cast<int>(std::sqrt(num_cells));
        int ny = (num_cells + nx - 1) / nx;

        for (int cell = 0; cell < num_cells; cell++)
        {
            int i = cell % nx;
            int j = cell / nx;

            // Left neighbor
            data.neighbors[cell * 4 + 0] = (i > 0) ? (j * nx + i - 1) : -1;
            // Bottom neighbor
            data.neighbors[cell * 4 + 1] = (j > 0) ? ((j - 1) * nx + i) : -1;
            // Right neighbor
            data.neighbors[cell * 4 + 2] = (i < nx - 1 && j * nx + i + 1 < num_cells) ? (j * nx + i + 1) : -1;
            // Top neighbor
            data.neighbors[cell * 4 + 3] = (j < ny - 1 && (j + 1) * nx + i < num_cells) ? ((j + 1) * nx + i) : -1;
        }
    }
    else
    {
        // Random connectivity
        std::uniform_int_distribution<> neighbor_dist(-1, num_cells - 1);
        for (int i = 0; i < num_cells * 4; i++)
        {
            int neighbor = neighbor_dist(gen);
            data.neighbors[i] = (neighbor < num_cells && neighbor >= 0) ? neighbor : -1;
        }
    }

    // Generate conservative variables (rho, rho*u, rho*v, E)
    data.conservative_vars.resize(num_cells * 4);
    for (int cell = 0; cell < num_cells; cell++)
    {
        double rho = rho_dist(gen);
        double u = vel_dist(gen);
        double v = vel_dist(gen);
        double E = energy_dist(gen);

        data.conservative_vars[cell * 4 + 0] = rho;
        data.conservative_vars[cell * 4 + 1] = rho * u;
        data.conservative_vars[cell * 4 + 2] = rho * v;
        data.conservative_vars[cell * 4 + 3] = E;
    }

    return data;
}

//=============================================================================
// TEST FUNCTIONS
//=============================================================================

/**
 * Test basic CUDA functionality
 */
bool test_cuda_availability()
{
    std::cout << "=== Testing CUDA Availability ===" << std::endl;

    int device_count;
    cudaError_t error = cudaGetDeviceCount(&device_count);

    if (error != cudaSuccess)
    {
        std::cout << "❌ CUDA Error: " << cudaGetErrorString(error) << std::endl;
        return false;
    }

    if (device_count == 0)
    {
        std::cout << "❌ No CUDA devices found" << std::endl;
        return false;
    }

    std::cout << "✅ Found " << device_count << " CUDA device(s)" << std::endl;

    // Print device info
    print_cuda_device_info_for_matrix_assembly();

    return true;
}

/**
 * Test dense matrix assembly
 */
bool test_dense_matrix_assembly()
{
    std::cout << "\n=== Testing Dense Matrix Assembly ===" << std::endl;

    // Generate small test data
    TestData data = generate_test_data(100); // Small problem for testing

    std::vector<double> matrix;

    try
    {
        double execution_time = assemble_dense_matrix_cuda(
            data.cell_areas, data.face_areas, data.neighbors, data.conservative_vars,
            matrix, data.dt, data.No_Physical_Cells);

        int matrix_size = 4 * data.No_Physical_Cells;

        std::cout << "✅ Dense matrix assembly completed" << std::endl;
        std::cout << "   Execution time: " << execution_time << " ms" << std::endl;
        std::cout << "   Matrix size: " << matrix_size << " x " << matrix_size << std::endl;
        std::cout << "   Matrix elements: " << matrix.size() << std::endl;

        // Basic validation
        if (matrix.size() != matrix_size * matrix_size)
        {
            std::cout << "❌ Matrix size mismatch" << std::endl;
            return false;
        }

        // Check for reasonable values
        double sum = 0.0, abs_sum = 0.0;
        int zero_count = 0;
        for (double val : matrix)
        {
            sum += val;
            abs_sum += std::abs(val);
            if (std::abs(val) < 1e-15)
                zero_count++;
        }

        std::cout << "   Matrix sum: " << sum << std::endl;
        std::cout << "   Matrix abs sum: " << abs_sum << std::endl;
        std::cout << "   Zero entries: " << zero_count << " ("
                  << (100.0 * zero_count / matrix.size()) << "%)" << std::endl;

        if (abs_sum < 1e-10)
        {
            std::cout << "⚠️  Matrix appears to be nearly zero - check implementation" << std::endl;
        }

        return true;
    }
    catch (const std::exception &e)
    {
        std::cout << "❌ Exception in dense matrix assembly: " << e.what() << std::endl;
        return false;
    }
}

/**
 * Test sparse matrix assembly
 */
bool test_sparse_matrix_assembly()
{
    std::cout << "\n=== Testing Sparse Matrix Assembly ===" << std::endl;

    TestData data = generate_test_data(200); // Medium problem for testing

    std::vector<int> row_indices, col_indices;
    std::vector<double> values;

    try
    {
        double execution_time = assemble_sparse_matrix_cuda(
            data.cell_areas, data.face_areas, data.neighbors, data.conservative_vars,
            row_indices, col_indices, values, data.dt, data.No_Physical_Cells);

        std::cout << "✅ Sparse matrix assembly completed" << std::endl;
        std::cout << "   Execution time: " << execution_time << " ms" << std::endl;
        std::cout << "   Non-zero entries: " << values.size() << std::endl;
        std::cout << "   Matrix size: " << (4 * data.No_Physical_Cells)
                  << " x " << (4 * data.No_Physical_Cells) << std::endl;

        if (row_indices.size() != values.size() || col_indices.size() != values.size())
        {
            std::cout << "❌ Sparse matrix array size mismatch" << std::endl;
            return false;
        }

        // Check sparsity
        int total_entries = 4 * data.No_Physical_Cells * 4 * data.No_Physical_Cells;
        double sparsity = 1.0 - (double)values.size() / total_entries;
        std::cout << "   Sparsity: " << (sparsity * 100) << "%" << std::endl;

        // Validate indices
        int max_row = *std::max_element(row_indices.begin(), row_indices.end());
        int max_col = *std::max_element(col_indices.begin(), col_indices.end());
        int min_row = *std::min_element(row_indices.begin(), row_indices.end());
        int min_col = *std::min_element(col_indices.begin(), col_indices.end());

        std::cout << "   Row indices range: [" << min_row << ", " << max_row << "]" << std::endl;
        std::cout << "   Col indices range: [" << min_col << ", " << max_col << "]" << std::endl;

        int matrix_dim = 4 * data.No_Physical_Cells;
        if (max_row >= matrix_dim || max_col >= matrix_dim || min_row < 0 || min_col < 0)
        {
            std::cout << "❌ Invalid indices detected" << std::endl;
            return false;
        }

        return true;
    }
    catch (const std::exception &e)
    {
        std::cout << "❌ Exception in sparse matrix assembly: " << e.what() << std::endl;
        return false;
    }
}

/**
 * Test vector b assembly
 */
bool test_vector_b_assembly()
{
    std::cout << "\n=== Testing Vector b Assembly ===" << std::endl;

    TestData data = generate_test_data(150);

    // Generate net flux data
    std::vector<double> net_flux(data.No_Physical_Cells * 4);
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> flux_dist(-1000.0, 1000.0);

    for (int i = 0; i < data.No_Physical_Cells * 4; i++)
    {
        net_flux[i] = flux_dist(gen);
    }

    std::vector<double> vector_b;

    try
    {
        double execution_time = assemble_vector_b_cuda(net_flux, vector_b, data.No_Physical_Cells);

        std::cout << "✅ Vector b assembly completed" << std::endl;
        std::cout << "   Execution time: " << execution_time << " ms" << std::endl;
        std::cout << "   Vector size: " << vector_b.size() << std::endl;

        if (vector_b.size() != data.No_Physical_Cells * 4)
        {
            std::cout << "❌ Vector size mismatch" << std::endl;
            return false;
        }

        // Validate that b = -net_flux
        double max_diff = 0.0;
        for (int i = 0; i < vector_b.size(); i++)
        {
            double diff = std::abs(vector_b[i] + net_flux[i]); // Should be b = -net_flux
            max_diff = std::max(max_diff, diff);
        }

        std::cout << "   Max difference from -net_flux: " << max_diff << std::endl;

        if (max_diff > 1e-10)
        {
            std::cout << "❌ Vector b computation incorrect" << std::endl;
            return false;
        }

        return true;
    }
    catch (const std::exception &e)
    {
        std::cout << "❌ Exception in vector b assembly: " << e.what() << std::endl;
        return false;
    }
}

/**
 * Performance scaling test
 */
void test_performance_scaling()
{
    std::cout << "\n=== Testing Performance Scaling ===" << std::endl;

    std::vector<int> problem_sizes = {100, 500, 1000, 2000, 5000};

    std::cout << "Problem Size | Dense Time (ms) | Sparse Time (ms) | Vector Time (ms)" << std::endl;
    std::cout << "-------------|-----------------|-------------------|------------------" << std::endl;

    for (int size : problem_sizes)
    {
        TestData data = generate_test_data(size);

        // Test dense matrix assembly
        std::vector<double> matrix;
        double dense_time = -1.0;
        try
        {
            dense_time = assemble_dense_matrix_cuda(
                data.cell_areas, data.face_areas, data.neighbors, data.conservative_vars,
                matrix, data.dt, data.No_Physical_Cells);
        }
        catch (...)
        {
            dense_time = -1.0;
        }

        // Test sparse matrix assembly
        std::vector<int> row_indices, col_indices;
        std::vector<double> values;
        double sparse_time = -1.0;
        try
        {
            sparse_time = assemble_sparse_matrix_cuda(
                data.cell_areas, data.face_areas, data.neighbors, data.conservative_vars,
                row_indices, col_indices, values, data.dt, data.No_Physical_Cells);
        }
        catch (...)
        {
            sparse_time = -1.0;
        }

        // Test vector assembly
        std::vector<double> net_flux(size * 4, 1.0);
        std::vector<double> vector_b;
        double vector_time = -1.0;
        try
        {
            vector_time = assemble_vector_b_cuda(net_flux, vector_b, size);
        }
        catch (...)
        {
            vector_time = -1.0;
        }

        std::cout << std::setw(12) << size << " | ";
        if (dense_time > 0)
            std::cout << std::setw(15) << std::fixed << std::setprecision(3) << dense_time;
        else
            std::cout << std::setw(15) << "FAILED";
        std::cout << " | ";
        if (sparse_time > 0)
            std::cout << std::setw(17) << std::fixed << std::setprecision(3) << sparse_time;
        else
            std::cout << std::setw(17) << "FAILED";
        std::cout << " | ";
        if (vector_time > 0)
            std::cout << std::setw(16) << std::fixed << std::setprecision(3) << vector_time;
        else
            std::cout << std::setw(16) << "FAILED";
        std::cout << std::endl;
    }
}

//=============================================================================
// MAIN TEST PROGRAM
//=============================================================================

int main()
{
    std::cout << "CUDA Matrix Assembly Kernels Test Program" << std::endl;
    std::cout << "==========================================" << std::endl;

    bool all_tests_passed = true;

    // Test 1: CUDA availability
    if (!test_cuda_availability())
    {
        std::cout << "❌ CUDA not available - cannot run tests" << std::endl;
        return 1;
    }

    // Test 2: Dense matrix assembly
    if (!test_dense_matrix_assembly())
    {
        all_tests_passed = false;
    }

    // Test 3: Sparse matrix assembly
    if (!test_sparse_matrix_assembly())
    {
        all_tests_passed = false;
    }

    // Test 4: Vector b assembly
    if (!test_vector_b_assembly())
    {
        all_tests_passed = false;
    }

    // Test 5: Performance scaling
    test_performance_scaling();

    std::cout << "\n=== Test Summary ===" << std::endl;
    if (all_tests_passed)
    {
        std::cout << "✅ All tests PASSED" << std::endl;
        std::cout << "   CUDA matrix assembly kernels are working correctly" << std::endl;
        std::cout << "   Ready for integration with CFD solver" << std::endl;
    }
    else
    {
        std::cout << "❌ Some tests FAILED" << std::endl;
        std::cout << "   Review implementation before integration" << std::endl;
    }

    return all_tests_passed ? 0 : 1;
}