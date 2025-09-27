// Grid_Computations_Optimized.h
// Header for optimized grid computation functions
// Provides performance-enhanced alternatives to standard grid operations

#ifndef GRID_COMPUTATIONS_OPTIMIZED_H
#define GRID_COMPUTATIONS_OPTIMIZED_H

#include "definitions.h"
#include "Grid.h"
#include "Globals.h"
#include <vector>
#include <array>
#include <chrono>

// ===== PERFORMANCE TRACKING =====
namespace GridPerformance
{
    void start_timer();
    void end_timer(const std::string &operation_name);
    void print_performance_summary();
}

// ===== OPTIMIZED MATHEMATICAL OPERATIONS =====

// Fast mathematical operations
inline float fast_inv_sqrt(float x);
inline double distance_squared(const V_D &P1, const V_D &P2);

// SIMD-optimized operations (when available)
#ifdef __AVX2__
inline void simd_cross_product(const double *a, const double *b, double *result);
#endif

// Fast trigonometric approximations
inline double fast_atan2(double y, double x);

// ===== OPTIMIZED GRID CONSTRUCTION FUNCTIONS =====

// Memory optimization
void optimize_memory_allocation(int num_physical_cells, int num_ghost_cells);

// Enhanced cell construction with move semantics
void Construct_Cell_Optimized(Cell &&Grid_Cell);

// Optimized face construction
void Construct_Face_Optimized(Cell &Grid_Cell);

// Parallel distance calculations
void Calculate_Cell_Center_Distances_Optimized();

// Enhanced centroid computation
void Compute_Centroid_Optimized(const V_D &vertices, V_D &centroid);

// Fast anti-clockwise sorting
void Sort_Points_AntiClockWise_Optimized(V_D &points, V_I &indices);

// ===== CUDA ACCELERATION WRAPPERS =====
#ifdef USE_CUDA

// Hybrid CPU/GPU grid construction
bool Construct_Cells_Hybrid(bool use_gpu = true);

// GPU-accelerated distance calculations
void Calculate_Distances_GPU();

#endif

// ===== PERFORMANCE ANALYSIS =====

// Grid performance analysis
void analyze_grid_performance();

// Benchmark comparison functions
void benchmark_grid_functions();

// ===== INTEGRATION INTERFACE =====

// Drop-in replacement namespace for existing grid functions
namespace OptimizedGrid
{

    // Initialize optimized grid processing
    void Initialize();

    // Process all cells with optimal performance
    void ProcessCells();

    // Finalize and report performance
    void Finalize();

    // Configuration options
    struct Config
    {
        bool use_gpu = true;                     // Enable CUDA acceleration when available
        bool use_parallel_cpu = true;            // Enable CPU parallelization
        bool use_simd = true;                    // Enable SIMD operations
        bool enable_performance_tracking = true; // Track and report performance
        int gpu_threshold = 1000;                // Minimum cells to use GPU
        int parallel_threshold = 100;            // Minimum cells for CPU parallelization
    };

    extern Config config;

    // Performance comparison utilities
    void compare_with_original();
    void memory_usage_analysis();
    void grid_quality_metrics();
}

// ===== OPTIMIZATION RECOMMENDATIONS =====

namespace GridOptimization
{

    // Analyze current grid and suggest optimizations
    struct OptimizationReport
    {
        bool recommend_gpu;
        bool recommend_memory_reallocation;
        bool recommend_data_structure_changes;
        double estimated_speedup;
        size_t memory_savings_bytes;
        std::vector<std::string> recommendations;
    };

    OptimizationReport analyze_optimization_potential();
    void apply_recommended_optimizations(const OptimizationReport &report);

    // Grid-specific optimizations
    void optimize_cell_connectivity();
    void optimize_memory_layout();
    void optimize_computation_order();
}

// ===== COMPATIBILITY LAYER =====

// Macro definitions for easy switching between original and optimized versions
#ifdef USE_OPTIMIZED_GRID
#define CONSTRUCT_CELL(cell) Construct_Cell_Optimized(std::move(cell))
#define CALCULATE_DISTANCES() Calculate_Cell_Center_Distances_Optimized()
#define SORT_POINTS(points, indices) Sort_Points_AntiClockWise_Optimized(points, indices)
#define COMPUTE_CENTROID(vertices, centroid) Compute_Centroid_Optimized(vertices, centroid)
#else
#define CONSTRUCT_CELL(cell) Construct_Cell(cell)
#define CALCULATE_DISTANCES() Calculate_Cell_Center_Distances()
#define SORT_POINTS(points, indices) Sort_Points_AntiClockWise(points, indices)
#define COMPUTE_CENTROID(vertices, centroid) Compute_Centroid(vertices, centroid)
#endif

// ===== PERFORMANCE CONSTANTS =====

namespace GridConstants
{
    // Performance thresholds
    constexpr int MIN_CELLS_FOR_GPU = 1000;
    constexpr int MIN_CELLS_FOR_PARALLEL = 100;
    constexpr int SIMD_ALIGNMENT = 32; // AVX2 alignment

    // Memory optimization constants
    constexpr double MEMORY_GROWTH_FACTOR = 1.5;
    constexpr size_t CACHE_LINE_SIZE = 64;

    // Numerical precision constants
    constexpr double FAST_MATH_EPSILON = 1e-6;
    constexpr double ANGLE_EPSILON = 1e-9;
}

// ===== ERROR HANDLING =====

enum class GridOptimizationError
{
    NONE = 0,
    INSUFFICIENT_MEMORY,
    CUDA_NOT_AVAILABLE,
    SIMD_NOT_SUPPORTED,
    INVALID_GRID_DATA,
    PERFORMANCE_DEGRADATION
};

class GridOptimizationException : public std::exception
{
private:
    GridOptimizationError error_type;
    std::string message;

public:
    GridOptimizationException(GridOptimizationError type, const std::string &msg)
        : error_type(type), message(msg) {}

    const char *what() const noexcept override
    {
        return message.c_str();
    }

    GridOptimizationError get_error_type() const
    {
        return error_type;
    }
};

#endif // GRID_COMPUTATIONS_OPTIMIZED_H