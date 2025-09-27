// Grid_Cuda_Host_Wrappers.cpp
// Host wrapper functions for Grid CUDA kernels
// Provides easy-to-use interface for GPU-accelerated grid operations

#include "Grid_Cuda_Kernels.h"
#include "Cuda_Kernel_Utilities.h"
#include <iostream>
#include <vector>

// ===== HOST WRAPPER IMPLEMENTATIONS =====

// Host wrapper for VTK grid construction
cudaError_t launch_construct_cells_from_vtk(
    const double *h_point_coords, const int *h_cell_connectivity,
    const int *h_cell_offsets, const int *h_cell_types,
    double *h_cell_areas, double *h_cell_centers,
    double *h_face_areas, double *h_face_normals, double *h_face_centers,
    int *h_cell_face_connectivity, int *h_face_cell_connectivity,
    int num_points, int num_cells, int num_faces,
    int total_connectivity, int max_faces_per_cell)
{
    // Allocate device memory
    double *d_point_coords, *d_cell_areas, *d_cell_centers;
    double *d_face_areas, *d_face_normals, *d_face_centers;
    int *d_cell_connectivity, *d_cell_offsets, *d_cell_types;
    int *d_cell_face_connectivity, *d_face_cell_connectivity;

    // Point coordinates
    CUDA_CHECK(cudaMalloc(&d_point_coords, num_points * 3 * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_point_coords, h_point_coords,
                          num_points * 3 * sizeof(double), cudaMemcpyHostToDevice));

    // Cell connectivity data
    CUDA_CHECK(cudaMalloc(&d_cell_connectivity, total_connectivity * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cell_connectivity, h_cell_connectivity,
                          total_connectivity * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_offsets, (num_cells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cell_offsets, h_cell_offsets,
                          (num_cells + 1) * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_types, num_cells * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cell_types, h_cell_types,
                          num_cells * sizeof(int), cudaMemcpyHostToDevice));

    // Output arrays
    CUDA_CHECK(cudaMalloc(&d_cell_areas, num_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_cell_centers, num_cells * 3 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_face_areas, num_faces * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_face_normals, num_faces * 3 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_face_centers, num_faces * 3 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_cell_face_connectivity, num_cells * max_faces_per_cell * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_face_cell_connectivity, num_faces * 2 * sizeof(int)));

    // Launch kernel
    dim3 gridSize = calculateGridSizeForGrid(num_cells);
    dim3 blockSize = calculateBlockSize();

    construct_cells_from_vtk_kernel<<<gridSize, blockSize>>>(
        d_point_coords, d_cell_connectivity, d_cell_offsets, d_cell_types,
        d_cell_areas, d_cell_centers, d_face_areas, d_face_normals, d_face_centers,
        d_cell_face_connectivity, d_face_cell_connectivity,
        num_cells, max_faces_per_cell);

    CUDA_CHECK_KERNEL("construct_cells_from_vtk_kernel");

    // Copy results back to host
    CUDA_CHECK(cudaMemcpy(h_cell_areas, d_cell_areas,
                          num_cells * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_cell_centers, d_cell_centers,
                          num_cells * 3 * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_face_areas, d_face_areas,
                          num_faces * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_face_normals, d_face_normals,
                          num_faces * 3 * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_face_centers, d_face_centers,
                          num_faces * 3 * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_cell_face_connectivity, d_cell_face_connectivity,
                          num_cells * max_faces_per_cell * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_face_cell_connectivity, d_face_cell_connectivity,
                          num_faces * 2 * sizeof(int), cudaMemcpyDeviceToHost));

    // Cleanup device memory
    cudaFree(d_point_coords);
    cudaFree(d_cell_connectivity);
    cudaFree(d_cell_offsets);
    cudaFree(d_cell_types);
    cudaFree(d_cell_areas);
    cudaFree(d_cell_centers);
    cudaFree(d_face_areas);
    cudaFree(d_face_normals);
    cudaFree(d_face_centers);
    cudaFree(d_cell_face_connectivity);
    cudaFree(d_face_cell_connectivity);

    return cudaSuccess;
}

// Host wrapper for neighbor identification
cudaError_t launch_identify_neighbors(
    const int *h_face_cell_connectivity, const int *h_cell_face_connectivity,
    const int *h_num_faces_per_cell, int *h_cell_neighbors, int *h_num_neighbors_per_cell,
    int num_cells, int num_faces, int max_faces_per_cell, int max_neighbors)
{
    // Allocate device memory
    int *d_face_cell_connectivity, *d_cell_face_connectivity, *d_num_faces_per_cell;
    int *d_cell_neighbors, *d_num_neighbors_per_cell;

    CUDA_CHECK(cudaMalloc(&d_face_cell_connectivity, num_faces * 2 * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_face_cell_connectivity, h_face_cell_connectivity,
                          num_faces * 2 * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_face_connectivity, num_cells * max_faces_per_cell * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cell_face_connectivity, h_cell_face_connectivity,
                          num_cells * max_faces_per_cell * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_num_faces_per_cell, num_cells * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_num_faces_per_cell, h_num_faces_per_cell,
                          num_cells * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_neighbors, num_cells * max_neighbors * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_num_neighbors_per_cell, num_cells * sizeof(int)));

    // Launch kernel
    dim3 gridSize = calculateGridSizeForGrid(num_cells);
    dim3 blockSize = calculateBlockSize();

    identify_neighbors_kernel<<<gridSize, blockSize>>>(
        d_face_cell_connectivity, d_cell_face_connectivity, d_num_faces_per_cell,
        d_cell_neighbors, d_num_neighbors_per_cell,
        num_cells, num_faces, max_faces_per_cell, max_neighbors);

    CUDA_CHECK_KERNEL("identify_neighbors_kernel");

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_cell_neighbors, d_cell_neighbors,
                          num_cells * max_neighbors * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_num_neighbors_per_cell, d_num_neighbors_per_cell,
                          num_cells * sizeof(int), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_face_cell_connectivity);
    cudaFree(d_cell_face_connectivity);
    cudaFree(d_num_faces_per_cell);
    cudaFree(d_cell_neighbors);
    cudaFree(d_num_neighbors_per_cell);

    return cudaSuccess;
}

// Host wrapper for geometric computations
cudaError_t launch_geometric_computations(
    const double *h_vec_a, const double *h_vec_b,
    double *h_cross_result, double *h_dot_result,
    int num_vectors)
{
    // Allocate device memory
    double *d_vec_a, *d_vec_b, *d_cross_result, *d_dot_result;

    CUDA_CHECK(cudaMalloc(&d_vec_a, num_vectors * 3 * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_vec_a, h_vec_a,
                          num_vectors * 3 * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_vec_b, num_vectors * 3 * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_vec_b, h_vec_b,
                          num_vectors * 3 * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cross_result, num_vectors * 3 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_dot_result, num_vectors * sizeof(double)));

    // Launch kernels
    dim3 gridSize = calculateGridSizeForGrid(num_vectors);
    dim3 blockSize = calculateBlockSize();

    cross_product_3d_kernel<<<gridSize, blockSize>>>(d_vec_a, d_vec_b, d_cross_result, num_vectors);
    CUDA_CHECK_KERNEL("cross_product_3d_kernel");

    dot_product_3d_kernel<<<gridSize, blockSize>>>(d_vec_a, d_vec_b, d_dot_result, num_vectors);
    CUDA_CHECK_KERNEL("dot_product_3d_kernel");

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_cross_result, d_cross_result,
                          num_vectors * 3 * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_dot_result, d_dot_result,
                          num_vectors * sizeof(double), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_vec_a);
    cudaFree(d_vec_b);
    cudaFree(d_cross_result);
    cudaFree(d_dot_result);

    return cudaSuccess;
}

// Host wrapper for boundary classification
cudaError_t launch_classify_boundaries(
    const int *h_face_cell_connectivity, const double *h_face_centers,
    int *h_boundary_markers, double *h_boundary_distances,
    int num_faces, double domain_bounds[6], double tolerance)
{
    // Allocate device memory
    int *d_face_cell_connectivity, *d_boundary_markers;
    double *d_face_centers, *d_boundary_distances;

    CUDA_CHECK(cudaMalloc(&d_face_cell_connectivity, num_faces * 2 * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_face_cell_connectivity, h_face_cell_connectivity,
                          num_faces * 2 * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_face_centers, num_faces * 3 * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_face_centers, h_face_centers,
                          num_faces * 3 * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_boundary_markers, num_faces * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_boundary_distances, num_faces * sizeof(double)));

    // Launch kernel
    dim3 gridSize = calculateGridSizeForGrid(num_faces);
    dim3 blockSize = calculateBlockSize();

    classify_boundary_faces_kernel<<<gridSize, blockSize>>>(
        d_face_cell_connectivity, d_face_centers,
        d_boundary_markers, d_boundary_distances, num_faces,
        domain_bounds[0], domain_bounds[1], domain_bounds[2],
        domain_bounds[3], domain_bounds[4], domain_bounds[5], tolerance);

    CUDA_CHECK_KERNEL("classify_boundary_faces_kernel");

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_boundary_markers, d_boundary_markers,
                          num_faces * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_boundary_distances, d_boundary_distances,
                          num_faces * sizeof(double), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_face_cell_connectivity);
    cudaFree(d_face_centers);
    cudaFree(d_boundary_markers);
    cudaFree(d_boundary_distances);

    return cudaSuccess;
}

// Host wrapper for grid quality assessment
cudaError_t launch_grid_quality_assessment(
    const double *h_cell_areas, const double *h_cell_centers,
    const int *h_cell_neighbors, const int *h_num_neighbors_per_cell,
    double *h_aspect_ratios, double *h_skewness_metrics, double *h_orthogonality_metrics,
    int num_cells, int max_neighbors)
{
    // Allocate device memory
    double *d_cell_areas, *d_cell_centers;
    int *d_cell_neighbors, *d_num_neighbors_per_cell;
    double *d_aspect_ratios, *d_skewness_metrics, *d_orthogonality_metrics;

    CUDA_CHECK(cudaMalloc(&d_cell_areas, num_cells * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_cell_areas, h_cell_areas,
                          num_cells * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_centers, num_cells * 3 * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_cell_centers, h_cell_centers,
                          num_cells * 3 * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_neighbors, num_cells * max_neighbors * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cell_neighbors, h_cell_neighbors,
                          num_cells * max_neighbors * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_num_neighbors_per_cell, num_cells * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_num_neighbors_per_cell, h_num_neighbors_per_cell,
                          num_cells * sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_aspect_ratios, num_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_skewness_metrics, num_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_orthogonality_metrics, num_cells * sizeof(double)));

    // Launch kernel
    dim3 gridSize = calculateGridSizeForGrid(num_cells);
    dim3 blockSize = calculateBlockSize();

    calculate_grid_quality_kernel<<<gridSize, blockSize>>>(
        d_cell_areas, d_cell_centers, d_cell_neighbors, d_num_neighbors_per_cell,
        d_aspect_ratios, d_skewness_metrics, d_orthogonality_metrics,
        num_cells, max_neighbors);

    CUDA_CHECK_KERNEL("calculate_grid_quality_kernel");

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_aspect_ratios, d_aspect_ratios,
                          num_cells * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_skewness_metrics, d_skewness_metrics,
                          num_cells * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_orthogonality_metrics, d_orthogonality_metrics,
                          num_cells * sizeof(double), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_cell_areas);
    cudaFree(d_cell_centers);
    cudaFree(d_cell_neighbors);
    cudaFree(d_num_neighbors_per_cell);
    cudaFree(d_aspect_ratios);
    cudaFree(d_skewness_metrics);
    cudaFree(d_orthogonality_metrics);

    return cudaSuccess;
}

// Host wrapper for grid transformations
cudaError_t launch_grid_transformations(
    double *h_point_coords, int num_points,
    double scale_factors[3], double translation[3], const double *rotation_matrix)
{
    // Allocate device memory
    double *d_point_coords, *d_rotation_matrix;

    CUDA_CHECK(cudaMalloc(&d_point_coords, num_points * 3 * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_point_coords, h_point_coords,
                          num_points * 3 * sizeof(double), cudaMemcpyHostToDevice));

    dim3 gridSize = calculateGridSizeForGrid(num_points);
    dim3 blockSize = calculateBlockSize();

    // Apply scaling if needed
    if (scale_factors[0] != 1.0 || scale_factors[1] != 1.0 || scale_factors[2] != 1.0)
    {
        scale_grid_kernel<<<gridSize, blockSize>>>(
            d_point_coords, scale_factors[0], scale_factors[1], scale_factors[2], num_points);
        CUDA_CHECK_KERNEL("scale_grid_kernel");
    }

    // Apply rotation if needed
    if (rotation_matrix != nullptr)
    {
        CUDA_CHECK(cudaMalloc(&d_rotation_matrix, 9 * sizeof(double)));
        CUDA_CHECK(cudaMemcpy(d_rotation_matrix, rotation_matrix,
                              9 * sizeof(double), cudaMemcpyHostToDevice));

        rotate_grid_kernel<<<gridSize, blockSize>>>(d_point_coords, d_rotation_matrix, num_points);
        CUDA_CHECK_KERNEL("rotate_grid_kernel");

        cudaFree(d_rotation_matrix);
    }

    // Apply translation if needed
    if (translation[0] != 0.0 || translation[1] != 0.0 || translation[2] != 0.0)
    {
        translate_grid_kernel<<<gridSize, blockSize>>>(
            d_point_coords, translation[0], translation[1], translation[2], num_points);
        CUDA_CHECK_KERNEL("translate_grid_kernel");
    }

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_point_coords, d_point_coords,
                          num_points * 3 * sizeof(double), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_point_coords);

    return cudaSuccess;
}

// Host wrapper for refinement marking
cudaError_t launch_mark_refinement_cells(
    const double *h_solution_gradients, const double *h_cell_areas,
    int *h_refinement_markers, double gradient_threshold,
    int num_cells, int num_vars)
{
    // Allocate device memory
    double *d_solution_gradients, *d_cell_areas;
    int *d_refinement_markers;

    CUDA_CHECK(cudaMalloc(&d_solution_gradients, num_cells * num_vars * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_solution_gradients, h_solution_gradients,
                          num_cells * num_vars * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_cell_areas, num_cells * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_cell_areas, h_cell_areas,
                          num_cells * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_refinement_markers, num_cells * sizeof(int)));

    // Launch kernel
    dim3 gridSize = calculateGridSizeForGrid(num_cells);
    dim3 blockSize = calculateBlockSize();

    mark_refinement_cells_kernel<<<gridSize, blockSize>>>(
        d_solution_gradients, d_cell_areas, d_refinement_markers,
        gradient_threshold, num_cells, num_vars);

    CUDA_CHECK_KERNEL("mark_refinement_cells_kernel");

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_refinement_markers, d_refinement_markers,
                          num_cells * sizeof(int), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_solution_gradients);
    cudaFree(d_cell_areas);
    cudaFree(d_refinement_markers);

    return cudaSuccess;
}

// ===== UTILITY FUNCTION IMPLEMENTATIONS =====

// Estimate memory requirements for grid operations
size_t estimate_grid_memory_usage(int num_points, int num_cells, int num_faces,
                                  int max_faces_per_cell, int max_neighbors)
{
    size_t total_memory = 0;

    // Point coordinates
    total_memory += num_points * 3 * sizeof(double);

    // Cell data
    total_memory += num_cells * sizeof(double);                   // Areas
    total_memory += num_cells * 3 * sizeof(double);               // Centers
    total_memory += num_cells * max_faces_per_cell * sizeof(int); // Face connectivity
    total_memory += num_cells * max_neighbors * sizeof(int);      // Neighbors

    // Face data
    total_memory += num_faces * sizeof(double);     // Areas
    total_memory += num_faces * 3 * sizeof(double); // Centers
    total_memory += num_faces * 3 * sizeof(double); // Normals
    total_memory += num_faces * 2 * sizeof(int);    // Cell connectivity

    // Quality metrics
    total_memory += num_cells * 3 * sizeof(double); // Aspect ratio, skewness, orthogonality

    return total_memory;
}

// Validate grid data consistency
bool validate_grid_data(int num_points, int num_cells, int num_faces,
                        const int *cell_connectivity, const int *cell_offsets)
{
    if (num_points <= 0 || num_cells <= 0 || num_faces <= 0)
    {
        std::cerr << "Invalid grid dimensions" << std::endl;
        return false;
    }

    if (cell_connectivity == nullptr || cell_offsets == nullptr)
    {
        std::cerr << "Null connectivity data" << std::endl;
        return false;
    }

    // Check connectivity indices are within bounds
    for (int i = 0; i < num_cells; i++)
    {
        int start = cell_offsets[i];
        int end = cell_offsets[i + 1];

        for (int j = start; j < end; j++)
        {
            if (cell_connectivity[j] < 0 || cell_connectivity[j] >= num_points)
            {
                std::cerr << "Invalid point index in connectivity: " << cell_connectivity[j] << std::endl;
                return false;
            }
        }
    }

    return true;
}