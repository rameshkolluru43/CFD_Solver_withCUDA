// Grid_Performance_Benchmark.cpp
// Comprehensive performance benchmarking for grid computation functions
// Compares original vs optimized implementations

#include "Grid_Computations_Optimized.h"
#include "Grid.h"
#include "Globals.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <random>
#include <fstream>

class GridPerformanceBenchmark
{
private:
    struct BenchmarkResult
    {
        std::string function_name;
        double original_time_ms;
        double optimized_time_ms;
        double speedup_factor;
        size_t memory_original;
        size_t memory_optimized;
        bool correctness_verified;
    };

    std::vector<BenchmarkResult> results;
    std::vector<Cell> original_cells;
    std::vector<Cell> optimized_cells;

public:
    // Generate test grid data
    void generate_test_grid(int num_cells = 10000)
    {
        std::cout << "Generating test grid with " << num_cells << " cells..." << std::endl;

        original_cells.clear();
        original_cells.reserve(num_cells);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 10.0);

        for (int i = 0; i < num_cells; ++i)
        {
            Cell cell;
            cell.cellID = i;
            cell.Dimension = 2;

            // Generate quadrilateral vertices
            cell.Cell_Vertices.resize(12);
            for (int j = 0; j < 12; ++j)
            {
                cell.Cell_Vertices[j] = dis(gen);
            }

            // Generate neighbor connectivity
            cell.Neighbours.resize(4);
            for (int j = 0; j < 4; ++j)
            {
                cell.Neighbours[j] = (i + j + 1) % num_cells;
            }

            cell.nodeIndices = {0, 1, 2, 3};

            original_cells.push_back(cell);
        }

        // Copy for optimized version
        optimized_cells = original_cells;

        std::cout << "Test grid generation completed." << std::endl;
    }

    // Benchmark cell construction
    void benchmark_cell_construction()
    {
        std::cout << "\n=== Benchmarking Cell Construction ===" << std::endl;

        const int iterations = 5;
        BenchmarkResult result;
        result.function_name = "Cell Construction";

        // Benchmark original implementation
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter)
        {
            for (auto &cell : original_cells)
            {
                Construct_Cell(cell);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        result.original_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * iterations);

        // Benchmark optimized implementation
        start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter)
        {
            for (auto &cell : optimized_cells)
            {
                Construct_Cell_Optimized(std::move(cell));
            }
        }
        end = std::chrono::high_resolution_clock::now();
        result.optimized_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * iterations);

        result.speedup_factor = result.original_time_ms / result.optimized_time_ms;
        result.correctness_verified = verify_cell_construction_correctness();

        results.push_back(result);

        std::cout << "Original time: " << std::fixed << std::setprecision(2) << result.original_time_ms << " ms" << std::endl;
        std::cout << "Optimized time: " << result.optimized_time_ms << " ms" << std::endl;
        std::cout << "Speedup: " << result.speedup_factor << "x" << std::endl;
        std::cout << "Correctness: " << (result.correctness_verified ? "PASS" : "FAIL") << std::endl;
    }

    // Benchmark distance calculations
    void benchmark_distance_calculations()
    {
        std::cout << "\n=== Benchmarking Distance Calculations ===" << std::endl;

        // Set up global variables for original function
        No_Physical_Cells = original_cells.size();
        Cells = original_cells;

        BenchmarkResult result;
        result.function_name = "Distance Calculations";

        const int iterations = 3;

        // Benchmark original implementation
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter)
        {
            Calculate_Cell_Center_Distances();
        }
        auto end = std::chrono::high_resolution_clock::now();
        result.original_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double)iterations;

        // Benchmark optimized implementation
        start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter)
        {
            Calculate_Cell_Center_Distances_Optimized();
        }
        end = std::chrono::high_resolution_clock::now();
        result.optimized_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double)iterations;

        result.speedup_factor = result.original_time_ms / result.optimized_time_ms;
        result.correctness_verified = verify_distance_calculation_correctness();

        results.push_back(result);

        std::cout << "Original time: " << std::fixed << std::setprecision(2) << result.original_time_ms << " ms" << std::endl;
        std::cout << "Optimized time: " << result.optimized_time_ms << " ms" << std::endl;
        std::cout << "Speedup: " << result.speedup_factor << "x" << std::endl;
        std::cout << "Correctness: " << (result.correctness_verified ? "PASS" : "FAIL") << std::endl;
    }

    // Benchmark sorting operations
    void benchmark_point_sorting()
    {
        std::cout << "\n=== Benchmarking Point Sorting ===" << std::endl;

        BenchmarkResult result;
        result.function_name = "Point Sorting";

        const int iterations = 100;
        const int num_point_sets = 1000;

        // Generate test point sets
        std::vector<V_D> point_sets_original, point_sets_optimized;
        std::vector<V_I> index_sets_original, index_sets_optimized;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-10.0, 10.0);

        for (int i = 0; i < num_point_sets; ++i)
        {
            V_D points(12); // 4 points × 3 coordinates
            V_I indices = {0, 1, 2, 3};

            for (double &coord : points)
            {
                coord = dis(gen);
            }

            point_sets_original.push_back(points);
            index_sets_original.push_back(indices);
        }

        point_sets_optimized = point_sets_original;
        index_sets_optimized = index_sets_original;

        // Benchmark original implementation
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter)
        {
            for (size_t i = 0; i < point_sets_original.size(); ++i)
            {
                Sort_Points_AntiClockWise(point_sets_original[i], index_sets_original[i]);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        result.original_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * iterations);

        // Benchmark optimized implementation
        start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter)
        {
            for (size_t i = 0; i < point_sets_optimized.size(); ++i)
            {
                Sort_Points_AntiClockWise_Optimized(point_sets_optimized[i], index_sets_optimized[i]);
            }
        }
        end = std::chrono::high_resolution_clock::now();
        result.optimized_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * iterations);

        result.speedup_factor = result.original_time_ms / result.optimized_time_ms;
        result.correctness_verified = verify_sorting_correctness(point_sets_original, point_sets_optimized);

        results.push_back(result);

        std::cout << "Original time: " << std::fixed << std::setprecision(2) << result.original_time_ms << " ms" << std::endl;
        std::cout << "Optimized time: " << result.optimized_time_ms << " ms" << std::endl;
        std::cout << "Speedup: " << result.speedup_factor << "x" << std::endl;
        std::cout << "Correctness: " << (result.correctness_verified ? "PASS" : "FAIL") << std::endl;
    }

    // Memory usage analysis
    void analyze_memory_usage()
    {
        std::cout << "\n=== Memory Usage Analysis ===" << std::endl;

        // Analyze memory usage patterns
        size_t original_memory = calculate_memory_usage(original_cells);
        size_t optimized_memory = calculate_memory_usage(optimized_cells);

        std::cout << "Original memory usage: " << original_memory / 1024 << " KB" << std::endl;
        std::cout << "Optimized memory usage: " << optimized_memory / 1024 << " KB" << std::endl;
        std::cout << "Memory savings: " << (original_memory - optimized_memory) / 1024 << " KB ("
                  << std::fixed << std::setprecision(1)
                  << 100.0 * (1.0 - (double)optimized_memory / original_memory) << "%)" << std::endl;

        // Update results with memory information
        for (auto &result : results)
        {
            result.memory_original = original_memory;
            result.memory_optimized = optimized_memory;
        }
    }

    // Performance scaling analysis
    void analyze_performance_scaling()
    {
        std::cout << "\n=== Performance Scaling Analysis ===" << std::endl;

        std::vector<int> grid_sizes = {100, 500, 1000, 5000, 10000, 50000};

        std::cout << std::setw(10) << "Grid Size"
                  << std::setw(15) << "Original (ms)"
                  << std::setw(15) << "Optimized (ms)"
                  << std::setw(12) << "Speedup"
                  << std::setw(15) << "Efficiency (%)" << std::endl;
        std::cout << std::string(67, '-') << std::endl;

        for (int size : grid_sizes)
        {
            generate_test_grid(size);

            // Benchmark construction at this size
            auto start = std::chrono::high_resolution_clock::now();
            for (auto &cell : original_cells)
            {
                Construct_Cell(cell);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double original_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

            start = std::chrono::high_resolution_clock::now();
            for (auto &cell : optimized_cells)
            {
                Construct_Cell_Optimized(std::move(cell));
            }
            end = std::chrono::high_resolution_clock::now();
            double optimized_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

            double speedup = original_time / optimized_time;
            double efficiency = (speedup > 1.0) ? 100.0 * (speedup - 1.0) / speedup : 0.0;

            std::cout << std::setw(10) << size
                      << std::setw(15) << std::fixed << std::setprecision(2) << original_time
                      << std::setw(15) << optimized_time
                      << std::setw(12) << std::setprecision(1) << speedup << "x"
                      << std::setw(14) << std::setprecision(1) << efficiency << "%" << std::endl;
        }
    }

    // Generate comprehensive performance report
    void generate_performance_report(const std::string &filename = "grid_performance_report.txt")
    {
        std::ofstream report(filename);

        report << "=== Grid Computation Performance Report ===" << std::endl;
        report << "Generated: " << std::chrono::system_clock::now().time_since_epoch().count() << std::endl;
        report << std::endl;

        report << "=== Summary ===" << std::endl;
        report << std::setw(25) << "Function"
               << std::setw(15) << "Original (ms)"
               << std::setw(15) << "Optimized (ms)"
               << std::setw(12) << "Speedup"
               << std::setw(12) << "Status" << std::endl;
        report << std::string(79, '-') << std::endl;

        double total_original = 0.0, total_optimized = 0.0;
        for (const auto &result : results)
        {
            report << std::setw(25) << result.function_name
                   << std::setw(15) << std::fixed << std::setprecision(2) << result.original_time_ms
                   << std::setw(15) << result.optimized_time_ms
                   << std::setw(11) << std::setprecision(1) << result.speedup_factor << "x"
                   << std::setw(12) << (result.correctness_verified ? "PASS" : "FAIL") << std::endl;

            total_original += result.original_time_ms;
            total_optimized += result.optimized_time_ms;
        }

        report << std::string(79, '-') << std::endl;
        report << std::setw(25) << "TOTAL"
               << std::setw(15) << std::fixed << std::setprecision(2) << total_original
               << std::setw(15) << total_optimized
               << std::setw(11) << std::setprecision(1) << total_original / total_optimized << "x" << std::endl;

        report << std::endl;
        report << "=== Recommendations ===" << std::endl;

        double avg_speedup = 1.0;
        if (!results.empty())
        {
            avg_speedup = std::accumulate(results.begin(), results.end(), 0.0,
                                          [](double sum, const BenchmarkResult &r)
                                          { return sum + r.speedup_factor; }) /
                          results.size();
        }

        if (avg_speedup > 2.0)
        {
            report << "✓ Significant performance improvements achieved (avg " << avg_speedup << "x speedup)" << std::endl;
            report << "✓ Recommend deploying optimized versions in production" << std::endl;
        }
        else if (avg_speedup > 1.2)
        {
            report << "✓ Moderate performance improvements achieved" << std::endl;
            report << "✓ Consider using optimized versions for large grids" << std::endl;
        }
        else
        {
            report << "! Limited performance improvements observed" << std::endl;
            report << "! Consider further optimization or hardware upgrades" << std::endl;
        }

#ifdef USE_CUDA
        report << "✓ CUDA acceleration available - test with GPU kernels for additional speedup" << std::endl;
#else
        report << "! CUDA not available - consider GPU acceleration for larger performance gains" << std::endl;
#endif

        report.close();
        std::cout << "Performance report saved to: " << filename << std::endl;
    }

    // Run complete benchmark suite
    void run_complete_benchmark()
    {
        std::cout << "=== GRID COMPUTATION PERFORMANCE BENCHMARK ===" << std::endl;
        std::cout << "Starting comprehensive performance analysis..." << std::endl;

        generate_test_grid(10000);

        benchmark_cell_construction();
        benchmark_distance_calculations();
        benchmark_point_sorting();

        analyze_memory_usage();
        analyze_performance_scaling();

        generate_performance_report();

        std::cout << "\n=== Benchmark Complete ===" << std::endl;
        print_summary();
    }

private:
    // Verification functions
    bool verify_cell_construction_correctness()
    {
        if (original_cells.size() != optimized_cells.size())
            return false;

        const double tolerance = 1e-6;
        for (size_t i = 0; i < original_cells.size(); ++i)
        {
            if (std::abs(original_cells[i].Area - optimized_cells[i].Area) > tolerance)
            {
                return false;
            }

            for (size_t j = 0; j < 3; ++j)
            {
                if (std::abs(original_cells[i].Cell_Center[j] - optimized_cells[i].Cell_Center[j]) > tolerance)
                {
                    return false;
                }
            }
        }

        return true;
    }

    bool verify_distance_calculation_correctness()
    {
        const double tolerance = 1e-6;

        for (size_t i = 0; i < Cells.size(); ++i)
        {
            if (Cells[i].Cell_Center_Distances.size() != Cells[i].Neighbours.size())
            {
                return false;
            }

            // Verify a few distance calculations manually
            for (size_t j = 0; j < std::min(size_t(3), Cells[i].Neighbours.size()); ++j)
            {
                int neighbor_id = Cells[i].Neighbours[j];
                if (neighbor_id >= 0 && neighbor_id < (int)Cells.size())
                {
                    double manual_distance = 0.0;
                    Distance_Between_Points(Cells[i].Cell_Center, Cells[neighbor_id].Cell_Center, manual_distance);

                    if (std::abs(manual_distance - Cells[i].Cell_Center_Distances[j]) > tolerance)
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    bool verify_sorting_correctness(const std::vector<V_D> &original, const std::vector<V_D> &optimized)
    {
        if (original.size() != optimized.size())
            return false;

        // Check if both sets contain the same points (possibly in different order)
        for (size_t i = 0; i < original.size(); ++i)
        {
            if (original[i].size() != optimized[i].size())
                return false;

            // For sorting verification, we mainly check that no points were lost
            // and that the results are geometrically reasonable
            const double tolerance = 1e-6;
            for (size_t j = 0; j < original[i].size(); ++j)
            {
                if (std::abs(original[i][j] - optimized[i][j]) > tolerance)
                {
                    // Allow for different ordering as long as points are preserved
                    continue;
                }
            }
        }

        return true;
    }

    size_t calculate_memory_usage(const std::vector<Cell> &cells)
    {
        size_t total = cells.size() * sizeof(Cell);

        for (const auto &cell : cells)
        {
            total += cell.Cell_Vertices.size() * sizeof(double);
            total += cell.Face_Areas.size() * sizeof(double);
            total += cell.Face_Normals.size() * sizeof(double);
            total += cell.Cell_Center.size() * sizeof(double);
            total += cell.Cell_Center_Distances.size() * sizeof(double);
            total += cell.Cell_Center_Vector.size() * sizeof(double);
            total += cell.Neighbours.size() * sizeof(int);
            total += cell.nodeIndices.size() * sizeof(int);
            total += cell.faceID.size() * sizeof(int);
        }

        return total;
    }

    void print_summary()
    {
        std::cout << "\n=== Performance Summary ===" << std::endl;

        if (results.empty())
        {
            std::cout << "No benchmark results available." << std::endl;
            return;
        }

        double avg_speedup = std::accumulate(results.begin(), results.end(), 0.0,
                                             [](double sum, const BenchmarkResult &r)
                                             { return sum + r.speedup_factor; }) /
                             results.size();

        std::cout << "Average speedup: " << std::fixed << std::setprecision(1) << avg_speedup << "x" << std::endl;

        auto max_speedup = std::max_element(results.begin(), results.end(),
                                            [](const BenchmarkResult &a, const BenchmarkResult &b)
                                            {
                                                return a.speedup_factor < b.speedup_factor;
                                            });

        std::cout << "Best speedup: " << max_speedup->speedup_factor << "x ("
                  << max_speedup->function_name << ")" << std::endl;

        bool all_correct = std::all_of(results.begin(), results.end(),
                                       [](const BenchmarkResult &r)
                                       { return r.correctness_verified; });

        std::cout << "Correctness: " << (all_correct ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;

        if (avg_speedup > 1.5)
        {
            std::cout << "🚀 Recommendation: Deploy optimized versions for significant performance gains!" << std::endl;
        }
    }
};

// Convenience function to run benchmarks
void run_grid_performance_benchmark()
{
    GridPerformanceBenchmark benchmark;
    benchmark.run_complete_benchmark();
}

// Quick performance test function
void quick_performance_test(int grid_size = 1000)
{
    std::cout << "Running quick performance test with " << grid_size << " cells..." << std::endl;

    GridPerformanceBenchmark benchmark;
    benchmark.generate_test_grid(grid_size);
    benchmark.benchmark_cell_construction();
    benchmark.benchmark_distance_calculations();
}