// Grid_Cuda_Kernels.h
// Header file for Grid CUDA Kernels
// Declares GPU-accelerated grid processing functions

#ifndef GRID_CUDA_KERNELS_H
#define GRID_CUDA_KERNELS_H

#include <cuda_runtime.h>
#include "../include/definitions.h"

// ===== GRID CONSTRUCTION KERNEL DECLARATIONS =====

// Parallel cell construction from VTK grid data
__global__ void construct_cells_from_vtk_kernel(
    const double *point_coords,   // [num_points * 3] - VTK point coordinates
    const int *cell_connectivity, // [total_connectivity] - VTK cell connectivity data
    const int *cell_offsets,      // [num_cells + 1] - Offsets into connectivity array
    const int *cell_types,        // [num_cells] - VTK cell types
    double *cell_areas,           // [num_cells] - Output cell areas/volumes
    double *cell_centers,         // [num_cells * 3] - Output cell centers
    double *face_areas,           // [num_faces] - Output face areas
    double *face_normals,         // [num_faces * 3] - Output face normals
    double *face_centers,         // [num_faces * 3] - Output face centers
    int *cell_face_connectivity,  // [num_cells * max_faces] - Cell to face connectivity
    int *face_cell_connectivity,  // [num_faces * 2] - Face to cell connectivity
    int num_cells,
    int max_faces_per_cell);

// Parallel neighbor identification
__global__ void identify_neighbors_kernel(
    const int *face_cell_connectivity, // [num_faces * 2] - Face to cell connectivity
    const int *cell_face_connectivity, // [num_cells * max_faces] - Cell to face connectivity
    const int *num_faces_per_cell,     // [num_cells] - Number of faces per cell
    int *cell_neighbors,               // [num_cells * max_neighbors] - Output neighbors
    int *num_neighbors_per_cell,       // [num_cells] - Output neighbor counts
    int num_cells,
    int num_faces,
    int max_faces_per_cell,
    int max_neighbors);

// ===== GEOMETRIC COMPUTATION KERNEL DECLARATIONS =====

// Parallel 3D cross product calculations
__global__ void cross_product_3d_kernel(
    const double *vec_a, // [num_vectors * 3] - First vectors
    const double *vec_b, // [num_vectors * 3] - Second vectors
    double *result,      // [num_vectors * 3] - Output cross products
    int num_vectors);

// Parallel 3D dot product calculations
__global__ void dot_product_3d_kernel(
    const double *vec_a, // [num_vectors * 3] - First vectors
    const double *vec_b, // [num_vectors * 3] - Second vectors
    double *result,      // [num_vectors] - Output dot products
    int num_vectors);

// Calculate face area vectors (area * normal)
__global__ void calculate_face_area_vectors_kernel(
    const double *face_coords,    // [num_faces * max_vertices * 3] - Face vertex coordinates
    const int *face_num_vertices, // [num_faces] - Number of vertices per face
    double *face_area_vectors,    // [num_faces * 3] - Output area vectors
    int num_faces,
    int max_vertices_per_face);

// ===== BOUNDARY CLASSIFICATION KERNEL DECLARATIONS =====

// Classify boundary faces and calculate distances
__global__ void classify_boundary_faces_kernel(
    const int *face_cell_connectivity, // [num_faces * 2] - Face to cell connectivity
    const double *face_centers,        // [num_faces * 3] - Face centers
    int *boundary_markers,             // [num_faces] - Output boundary markers
    double *boundary_distances,        // [num_faces] - Distance to nearest boundary
    int num_faces,
    double domain_xmin, double domain_xmax,
    double domain_ymin, double domain_ymax,
    double domain_zmin, double domain_zmax,
    double boundary_tolerance);

// ===== GRID QUALITY ASSESSMENT KERNEL DECLARATIONS =====

// Calculate comprehensive grid quality metrics
__global__ void calculate_grid_quality_kernel(
    const double *cell_areas,          // [num_cells] - Cell areas
    const double *cell_centers,        // [num_cells * 3] - Cell centers
    const int *cell_neighbors,         // [num_cells * max_neighbors] - Cell neighbors
    const int *num_neighbors_per_cell, // [num_cells] - Number of neighbors
    double *aspect_ratios,             // [num_cells] - Output aspect ratios
    double *skewness_metrics,          // [num_cells] - Output skewness
    double *orthogonality_metrics,     // [num_cells] - Output orthogonality
    int num_cells,
    int max_neighbors);

// ===== GRID TRANSFORMATION KERNEL DECLARATIONS =====

// Scale grid coordinates
__global__ void scale_grid_kernel(
    double *point_coords, // [num_points * 3] - Point coordinates to scale
    double scale_x,       // X-direction scaling factor
    double scale_y,       // Y-direction scaling factor
    double scale_z,       // Z-direction scaling factor
    int num_points);

// Translate grid coordinates
__global__ void translate_grid_kernel(
    double *point_coords, // [num_points * 3] - Point coordinates to translate
    double offset_x,      // X-direction offset
    double offset_y,      // Y-direction offset
    double offset_z,      // Z-direction offset
    int num_points);

// Rotate grid coordinates
__global__ void rotate_grid_kernel(
    double *point_coords,          // [num_points * 3] - Point coordinates to rotate
    const double *rotation_matrix, // [9] - 3x3 rotation matrix (row-major)
    int num_points);

// ===== GRID REFINEMENT KERNEL DECLARATIONS =====

// Mark cells for refinement based on solution gradients
__global__ void mark_refinement_cells_kernel(
    const double *solution_gradients, // [num_cells * num_vars] - Solution gradients
    const double *cell_areas,         // [num_cells] - Cell areas
    int *refinement_markers,          // [num_cells] - Output refinement markers
    double gradient_threshold,        // Threshold for refinement
    int num_cells,
    int num_vars);

// ===== HOST WRAPPER FUNCTIONS FOR GRID OPERATIONS =====

// Host wrapper for VTK grid construction
cudaError_t launch_construct_cells_from_vtk(
    const double *h_point_coords, const int *h_cell_connectivity,
    const int *h_cell_offsets, const int *h_cell_types,
    double *h_cell_areas, double *h_cell_centers,
    double *h_face_areas, double *h_face_normals, double *h_face_centers,
    int *h_cell_face_connectivity, int *h_face_cell_connectivity,
    int num_points, int num_cells, int num_faces,
    int total_connectivity, int max_faces_per_cell);

// Host wrapper for neighbor identification
cudaError_t launch_identify_neighbors(
    const int *h_face_cell_connectivity, const int *h_cell_face_connectivity,
    const int *h_num_faces_per_cell, int *h_cell_neighbors, int *h_num_neighbors_per_cell,
    int num_cells, int num_faces, int max_faces_per_cell, int max_neighbors);

// Host wrapper for geometric computations
cudaError_t launch_geometric_computations(
    const double *h_vec_a, const double *h_vec_b,
    double *h_cross_result, double *h_dot_result,
    int num_vectors);

// Host wrapper for boundary classification
cudaError_t launch_classify_boundaries(
    const int *h_face_cell_connectivity, const double *h_face_centers,
    int *h_boundary_markers, double *h_boundary_distances,
    int num_faces, double domain_bounds[6], double tolerance);

// Host wrapper for grid quality assessment
cudaError_t launch_grid_quality_assessment(
    const double *h_cell_areas, const double *h_cell_centers,
    const int *h_cell_neighbors, const int *h_num_neighbors_per_cell,
    double *h_aspect_ratios, double *h_skewness_metrics, double *h_orthogonality_metrics,
    int num_cells, int max_neighbors);

// Host wrapper for grid transformations
cudaError_t launch_grid_transformations(
    double *h_point_coords, int num_points,
    double scale_factors[3], double translation[3], const double *rotation_matrix);

// Host wrapper for refinement marking
cudaError_t launch_mark_refinement_cells(
    const double *h_solution_gradients, const double *h_cell_areas,
    int *h_refinement_markers, double gradient_threshold,
    int num_cells, int num_vars);

// ===== UTILITY FUNCTIONS =====

// Calculate optimal grid and block sizes for grid operations
inline dim3 calculateGridSizeForGrid(int num_elements, int block_size = 256)
{
    return dim3((num_elements + block_size - 1) / block_size, 1, 1);
}

// Estimate memory requirements for grid operations
size_t estimate_grid_memory_usage(int num_points, int num_cells, int num_faces,
                                  int max_faces_per_cell, int max_neighbors);

// Validate grid data consistency
bool validate_grid_data(int num_points, int num_cells, int num_faces,
                        const int *cell_connectivity, const int *cell_offsets);

#endif // GRID_CUDA_KERNELS_H