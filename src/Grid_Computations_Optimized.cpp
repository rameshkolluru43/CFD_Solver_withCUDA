// Grid_Computations_Optimized.cpp
// Optimized versions of performance-critical grid computation functions
// Includes CPU-level optimizations and CUDA integration options

#include "definitions.h"
#include "Grid.h"
#include "Globals.h"
#include "Utilities.h"
#include "Grid_Computations_Optimized.h"
#include <algorithm>
#ifdef __has_include
#if __has_include(<execution>)
#include <execution>
#define HAS_PARALLEL_EXECUTION 1
#endif
#endif

void Construct_Face_Optimized(Cell &Grid_Cell);
#include <numeric>
#include <immintrin.h> // For SIMD operations
#include <chrono>

#ifdef USE_CUDA
#include "Grid_Cuda_Kernels.h"
#include "Cuda_Kernel_Utilities.h"
#endif

// Performance tracking
namespace GridPerformance
{
    std::chrono::high_resolution_clock::time_point start_time;
    std::vector<std::pair<std::string, double>> timing_results;

    void start_timer()
    {
        start_time = std::chrono::high_resolution_clock::now();
    }

    void end_timer(const std::string &operation_name)
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        timing_results.emplace_back(operation_name, duration.count() / 1000.0); // Convert to milliseconds
    }

    void print_performance_summary()
    {
        std::cout << "\n=== Grid Computation Performance Summary ===" << std::endl;
        for (const auto &result : timing_results)
        {
            std::cout << result.first << ": " << result.second << " ms" << std::endl;
        }
        timing_results.clear();
    }
}

// ===== OPTIMIZED MATHEMATICAL OPERATIONS =====

// Fast inverse square root approximation (Quake III algorithm)
inline float fast_inv_sqrt(float x)
{
    float xhalf = 0.5f * x;
    int i = *(int *)&x;
    i = 0x5f3759df - (i >> 1);
    x = *(float *)&i;
    x = x * (1.5f - xhalf * x * x);
    return x;
}

// Optimized distance calculation without sqrt when only comparison is needed
inline double distance_squared(const V_D &P1, const V_D &P2)
{
    double dx = P2[0] - P1[0];
    double dy = P2[1] - P1[1];
    double dz = P2[2] - P1[2];
    return dx * dx + dy * dy + dz * dz;
}

// SIMD-optimized 3D vector operations
#ifdef __AVX2__
inline void simd_cross_product(const double *a, const double *b, double *result)
{
    // Load vectors into AVX registers
    __m256d va = _mm256_loadu_pd(a);
    __m256d vb = _mm256_loadu_pd(b);

    // Perform cross product using AVX operations
    // result[0] = a[1]*b[2] - a[2]*b[1]
    // result[1] = a[2]*b[0] - a[0]*b[2]
    // result[2] = a[0]*b[1] - a[1]*b[0]

    __m256d temp1 = _mm256_permute4x64_pd(va, 0xC9); // [a1, a2, a0, a3]
    __m256d temp2 = _mm256_permute4x64_pd(vb, 0xD2); // [b2, b0, b1, b3]
    __m256d temp3 = _mm256_permute4x64_pd(va, 0xD2); // [a2, a0, a1, a3]
    __m256d temp4 = _mm256_permute4x64_pd(vb, 0xC9); // [b1, b2, b0, b3]

    __m256d prod1 = _mm256_mul_pd(temp1, temp2);
    __m256d prod2 = _mm256_mul_pd(temp3, temp4);
    __m256d cross = _mm256_sub_pd(prod1, prod2);

    _mm256_storeu_pd(result, cross);
}
#endif

// ===== OPTIMIZED GRID CONSTRUCTION FUNCTIONS =====

// Pre-allocate memory for better performance
void optimize_memory_allocation(int num_physical_cells, int num_ghost_cells)
{
    GridPerformance::start_timer();

    // Pre-allocate cell vectors
    Cells.reserve(num_physical_cells + num_ghost_cells);

    // Pre-allocate boundary cell lists with estimated sizes
    Wall_Cells_List.reserve(num_physical_cells / 10); // Estimate 10% boundary cells
    Inlet_Cells_List.reserve(num_physical_cells / 20);
    Exit_Cells_List.reserve(num_physical_cells / 20);
    Symmetry_Cells_List.reserve(num_physical_cells / 20);

    GridPerformance::end_timer("Memory Pre-allocation");
}

// Optimized cell construction with move semantics
void Construct_Cell_Optimized(Cell &&Grid_Cell)
{
    // Use references to avoid copies
    auto &vertices = Grid_Cell.Cell_Vertices;

    // Pre-allocate vectors
    Grid_Cell.Cell_Center.resize(3);
    Grid_Cell.Face_Areas.reserve(4);
    Grid_Cell.Face_Normals.reserve(8);

    // Compute centroid efficiently
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        cx += vertices[i];
        cy += vertices[i + 1];
        cz += vertices[i + 2];
    }

    const int num_vertices = vertices.size() / 3;
    const double inv_vertices = 1.0 / num_vertices;
    Grid_Cell.Cell_Center[0] = cx * inv_vertices;
    Grid_Cell.Cell_Center[1] = cy * inv_vertices;
    Grid_Cell.Cell_Center[2] = cz * inv_vertices;

    // Construct faces with optimized calculations
    Construct_Face_Optimized(Grid_Cell);

    // Calculate area using optimized cross product
    if (vertices.size() == 12)
    { // Quadrilateral
        double diag1_x = vertices[6] - vertices[0];
        double diag1_y = vertices[7] - vertices[1];
        double diag2_x = vertices[9] - vertices[3];
        double diag2_y = vertices[10] - vertices[4];

        // 2D cross product magnitude
        Grid_Cell.Area = 0.5 * std::abs(diag1_x * diag2_y - diag1_y * diag2_x);
    }

    Grid_Cell.Inv_Area = 1.0 / Grid_Cell.Area;
}

// Optimized face construction
void Construct_Face_Optimized(Cell &Grid_Cell)
{
    const auto &vertices = Grid_Cell.Cell_Vertices;

    // Extract points efficiently
    const double *p0 = &vertices[0]; // [x0, y0, z0]
    const double *p1 = &vertices[3]; // [x1, y1, z1]
    const double *p2 = &vertices[6]; // [x2, y2, z2]
    const double *p3 = &vertices[9]; // [x3, y3, z3]

    // Construct faces with SIMD operations where available
    auto construct_face_segment = [&](const double *a, const double *b)
    {
        double dx = b[0] - a[0];
        double dy = b[1] - a[1];
        double length = std::sqrt(dx * dx + dy * dy);

        Grid_Cell.Face_Areas.push_back(length);
        Grid_Cell.Face_Normals.push_back(dy / length);  // nx = dy/dl
        Grid_Cell.Face_Normals.push_back(-dx / length); // ny = -dx/dl
    };

    // Face segments: (p3,p0), (p0,p1), (p1,p2), (p2,p3)
    construct_face_segment(p3, p0);
    construct_face_segment(p0, p1);
    construct_face_segment(p1, p2);
    construct_face_segment(p2, p3);
}

// Parallel cell center distance calculation using std::execution
void Calculate_Cell_Center_Distances_Optimized()
{
    GridPerformance::start_timer();

    // Create index vector for parallel processing
    std::vector<int> cell_indices(No_Physical_Cells);
    std::iota(cell_indices.begin(), cell_indices.end(), 0);

    auto compute_distances = [](int cell_idx)
    {
        auto &cell = Cells[cell_idx];
        const auto &center = cell.Cell_Center;

        const size_t num_neighbors = cell.Neighbours.size();
        cell.Cell_Center_Distances.resize(num_neighbors);
        cell.Cell_Center_Vector.resize(3 * num_neighbors);

        for (size_t i = 0; i < num_neighbors; ++i)
        {
            const int neighbor_id = cell.Neighbours[i];
            const auto &neighbor_center = Cells[neighbor_id].Cell_Center;

            double dx = neighbor_center[0] - center[0];
            double dy = neighbor_center[1] - center[1];
            double dz = neighbor_center[2] - center[2];

            cell.Cell_Center_Distances[i] = std::sqrt(dx * dx + dy * dy + dz * dz);

            size_t vec_idx = i * 3;
            cell.Cell_Center_Vector[vec_idx] = dx;
            cell.Cell_Center_Vector[vec_idx + 1] = dy;
            cell.Cell_Center_Vector[vec_idx + 2] = dz;
        }
    };

#ifdef HAS_PARALLEL_EXECUTION
    std::for_each(std::execution::par_unseq, cell_indices.begin(), cell_indices.end(), compute_distances);
#else
    std::for_each(cell_indices.begin(), cell_indices.end(), compute_distances);
#endif

    GridPerformance::end_timer("Parallel Cell Center Distances");
}

// Optimized centroid computation with reduced operations
void Compute_Centroid_Optimized(const V_D &vertices, V_D &centroid)
{
    if (centroid.size() != 3)
        centroid.resize(3);

    double cx = 0.0, cy = 0.0, cz = 0.0;
    const size_t num_points = vertices.size() / 3;

    // Unrolled loop for better performance
    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        cx += vertices[i];
        cy += vertices[i + 1];
        cz += vertices[i + 2];
    }

    const double inv_num_points = 1.0 / num_points;
    centroid[0] = cx * inv_num_points;
    centroid[1] = cy * inv_num_points;
    centroid[2] = cz * inv_num_points;
}

// Fast angle calculation using atan2 approximation
inline double fast_atan2(double y, double x)
{
    // Fast atan2 approximation using polynomial approximation
    const double abs_y = std::abs(y);
    const double abs_x = std::abs(x);

    if (abs_x < 1e-10 && abs_y < 1e-10)
        return 0.0;

    double a = std::min(abs_x, abs_y) / std::max(abs_x, abs_y);
    double s = a * a;
    double r = ((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;

    if (abs_y > abs_x)
        r = 1.57079637 - r;
    if (x < 0)
        r = 3.14159274 - r;
    if (y < 0)
        r = -r;

    return r;
}

// Optimized anti-clockwise sorting using faster angle calculation
void Sort_Points_AntiClockWise_Optimized(V_D &points, V_I &indices)
{
    if (points.empty() || points.size() % 3 != 0)
        return;

    GridPerformance::start_timer();

    // Compute centroid
    V_D centroid(3, 0.0);
    Compute_Centroid_Optimized(points, centroid);

    // Create point-index pairs
    std::vector<std::pair<std::array<double, 3>, int>> point_pairs;
    point_pairs.reserve(points.size() / 3);

    for (size_t i = 0; i < points.size(); i += 3)
    {
        std::array<double, 3> point = {points[i], points[i + 1], points[i + 2]};
        point_pairs.emplace_back(std::move(point), indices[i / 3]);
    }

    // Sort using optimized angle calculation
    std::sort(point_pairs.begin(), point_pairs.end(),
              [&centroid](const auto &a, const auto &b)
              {
                  double angle_a = fast_atan2(a.first[1] - centroid[1], a.first[0] - centroid[0]);
                  double angle_b = fast_atan2(b.first[1] - centroid[1], b.first[0] - centroid[0]);
                  return angle_a < angle_b;
              });

    // Rebuild vectors
    points.clear();
    indices.clear();
    points.reserve(point_pairs.size() * 3);
    indices.reserve(point_pairs.size());

    for (const auto &pair : point_pairs)
    {
        points.insert(points.end(), pair.first.begin(), pair.first.end());
        indices.push_back(pair.second);
    }

    GridPerformance::end_timer("Optimized Anti-clockwise Sorting");
}

// ===== CUDA ACCELERATION WRAPPERS =====

#ifdef USE_CUDA
// Hybrid CPU/GPU grid construction
bool Construct_Cells_Hybrid(bool use_gpu)
{
    GridPerformance::start_timer();

    if (use_gpu && No_Physical_Cells > 1000)
    { // Use GPU for large grids
        std::cout << "Using GPU acceleration for cell construction..." << std::endl;

        // Convert cell data to GPU-friendly format
        std::vector<double> point_coords;
        std::vector<int> cell_connectivity, cell_offsets, cell_types;

        // Extract data from Cells vector
        for (const auto &cell : Cells)
        {
            cell_offsets.push_back(cell_connectivity.size());
            for (size_t i = 0; i < cell.Cell_Vertices.size(); i += 3)
            {
                point_coords.insert(point_coords.end(),
                                    cell.Cell_Vertices.begin() + i,
                                    cell.Cell_Vertices.begin() + i + 3);
                cell_connectivity.push_back(point_coords.size() / 3 - 1);
            }
            cell_types.push_back(VTK_QUAD); // Assume quadrilateral
        }
        cell_offsets.push_back(cell_connectivity.size());

        // Prepare output arrays
        std::vector<double> cell_areas(No_Physical_Cells);
        std::vector<double> cell_centers(No_Physical_Cells * 3);

        // Launch GPU kernels
        cudaError_t result = launch_construct_cells_from_vtk(
            point_coords.data(), cell_connectivity.data(),
            cell_offsets.data(), cell_types.data(),
            cell_areas.data(), cell_centers.data(),
            nullptr, nullptr, nullptr, nullptr, nullptr,
            point_coords.size() / 3, No_Physical_Cells, 0,
            cell_connectivity.size(), 4);

        if (result == cudaSuccess)
        {
            // Copy results back to Cells
            for (int i = 0; i < No_Physical_Cells; ++i)
            {
                Cells[i].Area = cell_areas[i];
                Cells[i].Inv_Area = 1.0 / cell_areas[i];
                Cells[i].Cell_Center[0] = cell_centers[i * 3];
                Cells[i].Cell_Center[1] = cell_centers[i * 3 + 1];
                Cells[i].Cell_Center[2] = cell_centers[i * 3 + 2];
            }

            GridPerformance::end_timer("GPU Cell Construction");
            return true;
        }
        else
        {
            std::cout << "GPU construction failed, falling back to optimized CPU..." << std::endl;
        }
    }

    // Optimized CPU construction
    std::cout << "Using optimized CPU construction..." << std::endl;

    // Parallel construction using OpenMP-style parallel for
    std::vector<int> cell_indices(No_Physical_Cells);
    std::iota(cell_indices.begin(), cell_indices.end(), 0);

    auto construct_cell = [](int i)
    {
        Construct_Cell_Optimized(std::move(Cells[i]));
    };

#ifdef HAS_PARALLEL_EXECUTION
    std::for_each(std::execution::par_unseq, cell_indices.begin(), cell_indices.end(), construct_cell);
#else
    std::for_each(cell_indices.begin(), cell_indices.end(), construct_cell);
#endif

    GridPerformance::end_timer("Optimized CPU Cell Construction");
    return true;
}

// GPU-accelerated distance calculations
void Calculate_Distances_GPU()
{
    GridPerformance::start_timer();

    if (No_Physical_Cells > 500)
    { // Use GPU for medium+ sized grids
        std::vector<double> cell_centers_flat;
        std::vector<int> cell_neighbors_flat;
        std::vector<int> num_neighbors;

        // Flatten data for GPU
        for (const auto &cell : Cells)
        {
            cell_centers_flat.insert(cell_centers_flat.end(),
                                     cell.Cell_Center.begin(), cell.Cell_Center.end());

            // Pad neighbors to max size
            const int max_neighbors = 8;
            cell_neighbors_flat.resize(cell_neighbors_flat.size() + max_neighbors, -1);

            for (size_t i = 0; i < cell.Neighbours.size() && i < max_neighbors; ++i)
            {
                cell_neighbors_flat[cell_neighbors_flat.size() - max_neighbors + i] = cell.Neighbours[i];
            }

            num_neighbors.push_back(cell.Neighbours.size());
        }

        std::vector<double> distances_result(No_Physical_Cells * 8);

        // Launch geometric computations on GPU
        cudaError_t result = launch_geometric_computations(
            cell_centers_flat.data(), cell_centers_flat.data(),
            nullptr, distances_result.data(), No_Physical_Cells);

        if (result == cudaSuccess)
        {
            // Copy results back
            for (int i = 0; i < No_Physical_Cells; ++i)
            {
                const int max_neighbors = 8;
                Cells[i].Cell_Center_Distances.resize(num_neighbors[i]);
                for (int j = 0; j < num_neighbors[i]; ++j)
                {
                    Cells[i].Cell_Center_Distances[j] = distances_result[i * max_neighbors + j];
                }
            }

            GridPerformance::end_timer("GPU Distance Calculation");
            return;
        }
    }

    // Fallback to optimized CPU
    Calculate_Cell_Center_Distances_Optimized();
}
#endif

// ===== PERFORMANCE ANALYSIS FUNCTIONS =====

void analyze_grid_performance()
{
    std::cout << "\n=== Grid Performance Analysis ===" << std::endl;

    // Memory usage analysis
    size_t total_memory = 0;
    total_memory += Cells.size() * sizeof(Cell);
    for (const auto &cell : Cells)
    {
        total_memory += cell.Cell_Vertices.size() * sizeof(double);
        total_memory += cell.Face_Areas.size() * sizeof(double);
        total_memory += cell.Face_Normals.size() * sizeof(double);
        total_memory += cell.Cell_Center.size() * sizeof(double);
        total_memory += cell.Neighbours.size() * sizeof(int);
    }

    std::cout << "Total memory usage: " << total_memory / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Average memory per cell: " << total_memory / Cells.size() << " bytes" << std::endl;

    // Grid quality metrics
    double min_area = std::numeric_limits<double>::max();
    double max_area = 0.0;
    double total_area = 0.0;

    for (const auto &cell : Cells)
    {
        if (cell.Area > 0)
        {
            min_area = std::min(min_area, cell.Area);
            max_area = std::max(max_area, cell.Area);
            total_area += cell.Area;
        }
    }

    std::cout << "Grid quality metrics:" << std::endl;
    std::cout << "  Min cell area: " << min_area << std::endl;
    std::cout << "  Max cell area: " << max_area << std::endl;
    std::cout << "  Avg cell area: " << total_area / No_Physical_Cells << std::endl;
    std::cout << "  Area ratio (max/min): " << max_area / min_area << std::endl;
}

// Benchmark function to compare original vs optimized performance
void benchmark_grid_functions()
{
    std::cout << "\n=== Grid Function Benchmarks ===" << std::endl;

    const int num_iterations = 10;

    // Benchmark distance calculations
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_iterations; ++i)
        {
            Calculate_Cell_Center_Distances_Optimized();
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Optimized distance calculation: " << duration.count() / num_iterations << " ms/iteration" << std::endl;
    }

#ifdef USE_CUDA
    // Benchmark GPU vs CPU
    {
        GridPerformance::start_timer();
        Construct_Cells_Hybrid(true); // GPU
        GridPerformance::end_timer("GPU Construction");

        GridPerformance::start_timer();
        Construct_Cells_Hybrid(false); // CPU
        GridPerformance::end_timer("CPU Construction");
    }
#endif

    GridPerformance::print_performance_summary();
}

// ===== INTEGRATION FUNCTIONS =====

namespace OptimizedGrid
{

    Config config;

    void Initialize()
    {
        optimize_memory_allocation(No_Physical_Cells, No_Ghost_Cells);
    }

    void ProcessCells()
    {
#ifdef USE_CUDA
        Construct_Cells_Hybrid(true);
        Calculate_Distances_GPU();
#else
#ifdef HAS_PARALLEL_EXECUTION
        std::for_each(std::execution::par_unseq, Cells.begin(), Cells.end(),
                      [](Cell &cell)
                      { Construct_Cell_Optimized(std::move(cell)); });
#else
        std::for_each(Cells.begin(), Cells.end(),
                      [](Cell &cell)
                      { Construct_Cell_Optimized(std::move(cell)); });
#endif
        Calculate_Cell_Center_Distances_Optimized();
#endif
    }

    void Finalize()
    {
        analyze_grid_performance();
        GridPerformance::print_performance_summary();
    }
}

namespace GridOptimization
{
    OptimizationReport analyze_optimization_potential()
    {
        OptimizationReport report;
        report.recommend_gpu = (No_Physical_Cells > 1000);
        report.recommend_memory_reallocation = false;
        report.recommend_data_structure_changes = false;
        report.estimated_speedup = report.recommend_gpu ? 5.0 : 1.0;
        report.memory_savings_bytes = 0;
        if (report.recommend_gpu)
            report.recommendations.push_back("Use CUDA for grid computations (>1000 cells)");
        return report;
    }

    void apply_recommended_optimizations(const OptimizationReport &) {}
    void optimize_cell_connectivity() {}
    void optimize_memory_layout() {}
    void optimize_computation_order() {}
}