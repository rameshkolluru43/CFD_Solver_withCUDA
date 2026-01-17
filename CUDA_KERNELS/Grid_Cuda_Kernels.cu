// Grid_Cuda_Kernels.cu
// CUDA Kernels for Grid Generation and Processing Operations
// Implements GPU-accelerated versions of Grid.h functions

#include <cuda_runtime.h>
#include <math.h>
#include "Cuda_Kernel_Utilities.h"
#include "../include/definitions.h"

// ===== GRID CONSTRUCTION KERNELS =====

// Kernel for parallel cell construction from VTK grid data
__global__ void construct_cells_from_vtk_kernel(
    const double* point_coords,      // [num_points * 3] - VTK point coordinates
    const int* cell_connectivity,    // [total_connectivity] - VTK cell connectivity data
    const int* cell_offsets,         // [num_cells + 1] - Offsets into connectivity array
    const int* cell_types,          // [num_cells] - VTK cell types
    double* cell_areas,             // [num_cells] - Output cell areas/volumes
    double* cell_centers,           // [num_cells * 3] - Output cell centers
    double* face_areas,             // [num_faces] - Output face areas
    double* face_normals,           // [num_faces * 3] - Output face normals
    double* face_centers,           // [num_faces * 3] - Output face centers
    int* cell_face_connectivity,    // [num_cells * max_faces] - Cell to face connectivity
    int* face_cell_connectivity,    // [num_faces * 2] - Face to cell connectivity
    int num_cells,
    int max_faces_per_cell
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    // Get cell connectivity data
    int start_idx = cell_offsets[cell_idx];
    int end_idx = cell_offsets[cell_idx + 1];
    int num_vertices = end_idx - start_idx;
    
    // Calculate cell center as average of vertices
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int v = 0; v < num_vertices; v++) {
        int vertex_id = cell_connectivity[start_idx + v];
        cx += point_coords[vertex_id * 3 + 0];
        cy += point_coords[vertex_id * 3 + 1];
        cz += point_coords[vertex_id * 3 + 2];
    }
    cell_centers[cell_idx * 3 + 0] = cx / num_vertices;
    cell_centers[cell_idx * 3 + 1] = cy / num_vertices;
    cell_centers[cell_idx * 3 + 2] = cz / num_vertices;

    // Calculate cell area/volume based on cell type
    int cell_type = cell_types[cell_idx];
    double area_volume = 0.0;

    if (cell_type == VTK_QUAD && num_vertices == 4) {
        // Quadrilateral area calculation
        double x1 = point_coords[cell_connectivity[start_idx + 0] * 3 + 0];
        double y1 = point_coords[cell_connectivity[start_idx + 0] * 3 + 1];
        double x2 = point_coords[cell_connectivity[start_idx + 1] * 3 + 0];
        double y2 = point_coords[cell_connectivity[start_idx + 1] * 3 + 1];
        double x3 = point_coords[cell_connectivity[start_idx + 2] * 3 + 0];
        double y3 = point_coords[cell_connectivity[start_idx + 2] * 3 + 1];
        double x4 = point_coords[cell_connectivity[start_idx + 3] * 3 + 0];
        double y4 = point_coords[cell_connectivity[start_idx + 3] * 3 + 1];

        // Shoelace formula for quadrilateral
        area_volume = 0.5 * fabs((x1*y2 - x2*y1) + (x2*y3 - x3*y2) + 
                                 (x3*y4 - x4*y3) + (x4*y1 - x1*y4));
    }
    else if (cell_type == VTK_TRIANGLE && num_vertices == 3) {
        // Triangle area calculation
        double x1 = point_coords[cell_connectivity[start_idx + 0] * 3 + 0];
        double y1 = point_coords[cell_connectivity[start_idx + 0] * 3 + 1];
        double x2 = point_coords[cell_connectivity[start_idx + 1] * 3 + 0];
        double y2 = point_coords[cell_connectivity[start_idx + 1] * 3 + 1];
        double x3 = point_coords[cell_connectivity[start_idx + 2] * 3 + 0];
        double y3 = point_coords[cell_connectivity[start_idx + 2] * 3 + 1];

        // Cross product magnitude / 2
        area_volume = 0.5 * fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    }

    cell_areas[cell_idx] = area_volume;
}

// Kernel for identifying neighboring cells
__global__ void identify_neighbors_kernel(
    const int* face_cell_connectivity,  // [num_faces * 2] - Face to cell connectivity
    const int* cell_face_connectivity,  // [num_cells * max_faces] - Cell to face connectivity
    const int* num_faces_per_cell,      // [num_cells] - Number of faces per cell
    int* cell_neighbors,                // [num_cells * max_neighbors] - Output neighbors
    int* num_neighbors_per_cell,        // [num_cells] - Output neighbor counts
    int num_cells,
    int num_faces,
    int max_faces_per_cell,
    int max_neighbors
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    int neighbor_count = 0;
    int cell_num_faces = num_faces_per_cell[cell_idx];  // Renamed to avoid shadowing parameter

    // Check all faces of this cell
    for (int f = 0; f < cell_num_faces && f < max_faces_per_cell; f++) {
        int face_id = cell_face_connectivity[cell_idx * max_faces_per_cell + f];
        if (face_id >= 0 && face_id < num_faces) {
            // Get cells connected to this face
            int left_cell = face_cell_connectivity[face_id * 2 + 0];
            int right_cell = face_cell_connectivity[face_id * 2 + 1];
            
            // The neighbor is the other cell (not this one)
            int neighbor = (left_cell == cell_idx) ? right_cell : left_cell;
            
            if (neighbor >= 0 && neighbor < num_cells && neighbor_count < max_neighbors) {
                cell_neighbors[cell_idx * max_neighbors + neighbor_count] = neighbor;
                neighbor_count++;
            }
        }
    }

    num_neighbors_per_cell[cell_idx] = neighbor_count;
}

// ===== GEOMETRIC COMPUTATION KERNELS =====

// Kernel for parallel cross product calculations
__global__ void cross_product_3d_kernel(
    const double* vec_a,    // [num_vectors * 3] - First vectors
    const double* vec_b,    // [num_vectors * 3] - Second vectors
    double* result,         // [num_vectors * 3] - Output cross products
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors) return;

    double ax = vec_a[idx * 3 + 0], ay = vec_a[idx * 3 + 1], az = vec_a[idx * 3 + 2];
    double bx = vec_b[idx * 3 + 0], by = vec_b[idx * 3 + 1], bz = vec_b[idx * 3 + 2];

    result[idx * 3 + 0] = ay * bz - az * by;
    result[idx * 3 + 1] = az * bx - ax * bz;
    result[idx * 3 + 2] = ax * by - ay * bx;
}

// Kernel for parallel dot product calculations
__global__ void dot_product_3d_kernel(
    const double* vec_a,    // [num_vectors * 3] - First vectors
    const double* vec_b,    // [num_vectors * 3] - Second vectors
    double* result,         // [num_vectors] - Output dot products
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors) return;

    double ax = vec_a[idx * 3 + 0], ay = vec_a[idx * 3 + 1], az = vec_a[idx * 3 + 2];
    double bx = vec_b[idx * 3 + 0], by = vec_b[idx * 3 + 1], bz = vec_b[idx * 3 + 2];

    result[idx] = ax * bx + ay * by + az * bz;
}

// Kernel for calculating face area vectors (area * normal)
__global__ void calculate_face_area_vectors_kernel(
    const double* face_coords,      // [num_faces * max_vertices * 3] - Face vertex coordinates
    const int* face_num_vertices,   // [num_faces] - Number of vertices per face
    double* face_area_vectors,      // [num_faces * 3] - Output area vectors
    int num_faces,
    int max_vertices_per_face
) {
    int face_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (face_idx >= num_faces) return;

    int num_vertices = face_num_vertices[face_idx];
    
    if (num_vertices == 3) {
        // Triangle: cross product of two edges
        double x1 = face_coords[face_idx * max_vertices_per_face * 3 + 0];
        double y1 = face_coords[face_idx * max_vertices_per_face * 3 + 1];
        double z1 = face_coords[face_idx * max_vertices_per_face * 3 + 2];
        double x2 = face_coords[face_idx * max_vertices_per_face * 3 + 3];
        double y2 = face_coords[face_idx * max_vertices_per_face * 3 + 4];
        double z2 = face_coords[face_idx * max_vertices_per_face * 3 + 5];
        double x3 = face_coords[face_idx * max_vertices_per_face * 3 + 6];
        double y3 = face_coords[face_idx * max_vertices_per_face * 3 + 7];
        double z3 = face_coords[face_idx * max_vertices_per_face * 3 + 8];

        // Edge vectors
        double e1x = x2 - x1, e1y = y2 - y1, e1z = z2 - z1;
        double e2x = x3 - x1, e2y = y3 - y1, e2z = z3 - z1;

        // Cross product (gives 2 * area vector)
        double ax = e1y * e2z - e1z * e2y;
        double ay = e1z * e2x - e1x * e2z;
        double az = e1x * e2y - e1y * e2x;

        face_area_vectors[face_idx * 3 + 0] = 0.5 * ax;
        face_area_vectors[face_idx * 3 + 1] = 0.5 * ay;
        face_area_vectors[face_idx * 3 + 2] = 0.5 * az;
    }
    else if (num_vertices == 4) {
        // Quadrilateral: cross product of diagonals
        double x1 = face_coords[face_idx * max_vertices_per_face * 3 + 0];
        double y1 = face_coords[face_idx * max_vertices_per_face * 3 + 1];
        double z1 = face_coords[face_idx * max_vertices_per_face * 3 + 2];
        double x2 = face_coords[face_idx * max_vertices_per_face * 3 + 3];
        double y2 = face_coords[face_idx * max_vertices_per_face * 3 + 4];
        double z2 = face_coords[face_idx * max_vertices_per_face * 3 + 5];
        double x3 = face_coords[face_idx * max_vertices_per_face * 3 + 6];
        double y3 = face_coords[face_idx * max_vertices_per_face * 3 + 7];
        double z3 = face_coords[face_idx * max_vertices_per_face * 3 + 8];
        double x4 = face_coords[face_idx * max_vertices_per_face * 3 + 9];
        double y4 = face_coords[face_idx * max_vertices_per_face * 3 + 10];
        double z4 = face_coords[face_idx * max_vertices_per_face * 3 + 11];

        // Diagonal vectors
        double d1x = x3 - x1, d1y = y3 - y1, d1z = z3 - z1;
        double d2x = x4 - x2, d2y = y4 - y2, d2z = z4 - z2;

        // Cross product
        double ax = d1y * d2z - d1z * d2y;
        double ay = d1z * d2x - d1x * d2z;
        double az = d1x * d2y - d1y * d2x;

        face_area_vectors[face_idx * 3 + 0] = 0.5 * ax;
        face_area_vectors[face_idx * 3 + 1] = 0.5 * ay;
        face_area_vectors[face_idx * 3 + 2] = 0.5 * az;
    }
}

// ===== BOUNDARY CLASSIFICATION KERNELS =====

// Kernel for classifying boundary faces
__global__ void classify_boundary_faces_kernel(
    const int* face_cell_connectivity,  // [num_faces * 2] - Face to cell connectivity
    const double* face_centers,         // [num_faces * 3] - Face centers
    int* boundary_markers,              // [num_faces] - Output boundary markers
    double* boundary_distances,         // [num_faces] - Distance to nearest boundary
    int num_faces,
    double domain_xmin, double domain_xmax,
    double domain_ymin, double domain_ymax,
    double domain_zmin, double domain_zmax,
    double boundary_tolerance
) {
    int face_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (face_idx >= num_faces) return;

    // Check if face is on boundary (has only one neighboring cell)
    int left_cell = face_cell_connectivity[face_idx * 2 + 0];
    int right_cell = face_cell_connectivity[face_idx * 2 + 1];
    
    bool is_boundary = (left_cell < 0 || right_cell < 0);
    
    if (is_boundary) {
        double fx = face_centers[face_idx * 3 + 0];
        double fy = face_centers[face_idx * 3 + 1];
        double fz = face_centers[face_idx * 3 + 2];

        // Classify boundary type based on position
        if (fabs(fx - domain_xmin) < boundary_tolerance) {
            boundary_markers[face_idx] = BOUNDARY_XMIN;  // Left boundary
        }
        else if (fabs(fx - domain_xmax) < boundary_tolerance) {
            boundary_markers[face_idx] = BOUNDARY_XMAX;  // Right boundary
        }
        else if (fabs(fy - domain_ymin) < boundary_tolerance) {
            boundary_markers[face_idx] = BOUNDARY_YMIN;  // Bottom boundary
        }
        else if (fabs(fy - domain_ymax) < boundary_tolerance) {
            boundary_markers[face_idx] = BOUNDARY_YMAX;  // Top boundary
        }
        else if (fabs(fz - domain_zmin) < boundary_tolerance) {
            boundary_markers[face_idx] = BOUNDARY_ZMIN;  // Front boundary
        }
        else if (fabs(fz - domain_zmax) < boundary_tolerance) {
            boundary_markers[face_idx] = BOUNDARY_ZMAX;  // Back boundary
        }
        else {
            boundary_markers[face_idx] = BOUNDARY_WALL;   // Wall/obstacle boundary
        }

        // Calculate distance to nearest domain boundary
        double dist_x = fmin(fabs(fx - domain_xmin), fabs(fx - domain_xmax));
        double dist_y = fmin(fabs(fy - domain_ymin), fabs(fy - domain_ymax));
        double dist_z = fmin(fabs(fz - domain_zmin), fabs(fz - domain_zmax));
        boundary_distances[face_idx] = fmin(fmin(dist_x, dist_y), dist_z);
    }
    else {
        boundary_markers[face_idx] = BOUNDARY_NONE;  // Interior face
        boundary_distances[face_idx] = INFINITY;     // Far from boundary
    }
}

// ===== GRID QUALITY ASSESSMENT KERNELS =====

// Kernel for calculating grid quality metrics
__global__ void calculate_grid_quality_kernel(
    const double* cell_areas,          // [num_cells] - Cell areas
    const double* cell_centers,        // [num_cells * 3] - Cell centers
    const int* cell_neighbors,         // [num_cells * max_neighbors] - Cell neighbors
    const int* num_neighbors_per_cell, // [num_cells] - Number of neighbors
    double* aspect_ratios,             // [num_cells] - Output aspect ratios
    double* skewness_metrics,          // [num_cells] - Output skewness
    double* orthogonality_metrics,     // [num_cells] - Output orthogonality
    int num_cells,
    int max_neighbors
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    double area = cell_areas[cell_idx];
    double cx = cell_centers[cell_idx * 3 + 0];
    double cy = cell_centers[cell_idx * 3 + 1];
    double cz = cell_centers[cell_idx * 3 + 2];

    int num_neighbors = num_neighbors_per_cell[cell_idx];
    
    // Calculate characteristic length scales
    double min_distance = INFINITY;
    double max_distance = 0.0;
    double sum_distances = 0.0;
    
    for (int n = 0; n < num_neighbors && n < max_neighbors; n++) {
        int neighbor_id = cell_neighbors[cell_idx * max_neighbors + n];
        if (neighbor_id >= 0 && neighbor_id < num_cells) {
            double nx = cell_centers[neighbor_id * 3 + 0];
            double ny = cell_centers[neighbor_id * 3 + 1];
            double nz = cell_centers[neighbor_id * 3 + 2];
            
            double distance = sqrt((nx - cx)*(nx - cx) + (ny - cy)*(ny - cy) + (nz - cz)*(nz - cz));
            min_distance = fmin(min_distance, distance);
            max_distance = fmax(max_distance, distance);
            sum_distances += distance;
        }
    }

    // Aspect ratio: ratio of max to min neighbor distance
    if (min_distance > 0 && min_distance != INFINITY) {
        aspect_ratios[cell_idx] = max_distance / min_distance;
    } else {
        aspect_ratios[cell_idx] = 1.0;
    }

    // Skewness: deviation from regular spacing
    if (num_neighbors > 0) {
        double avg_distance = sum_distances / num_neighbors;
        double variance = 0.0;
        
        for (int n = 0; n < num_neighbors && n < max_neighbors; n++) {
            int neighbor_id = cell_neighbors[cell_idx * max_neighbors + n];
            if (neighbor_id >= 0 && neighbor_id < num_cells) {
                double nx = cell_centers[neighbor_id * 3 + 0];
                double ny = cell_centers[neighbor_id * 3 + 1];
                double nz = cell_centers[neighbor_id * 3 + 2];
                
                double distance = sqrt((nx - cx)*(nx - cx) + (ny - cy)*(ny - cy) + (nz - cz)*(nz - cz));
                variance += (distance - avg_distance) * (distance - avg_distance);
            }
        }
        
        if (avg_distance > 0) {
            skewness_metrics[cell_idx] = sqrt(variance / num_neighbors) / avg_distance;
        } else {
            skewness_metrics[cell_idx] = 0.0;
        }
    } else {
        skewness_metrics[cell_idx] = 0.0;
    }

    // Orthogonality: measure of grid alignment (simplified)
    orthogonality_metrics[cell_idx] = 1.0 / (1.0 + skewness_metrics[cell_idx]);
}

// ===== GRID TRANSFORMATION KERNELS =====

// Kernel for scaling grid coordinates
__global__ void scale_grid_kernel(
    double* point_coords,    // [num_points * 3] - Point coordinates to scale
    double scale_x,          // X-direction scaling factor
    double scale_y,          // Y-direction scaling factor
    double scale_z,          // Z-direction scaling factor
    int num_points
) {
    int point_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_idx >= num_points) return;

    point_coords[point_idx * 3 + 0] *= scale_x;
    point_coords[point_idx * 3 + 1] *= scale_y;
    point_coords[point_idx * 3 + 2] *= scale_z;
}

// Kernel for translating grid coordinates
__global__ void translate_grid_kernel(
    double* point_coords,    // [num_points * 3] - Point coordinates to translate
    double offset_x,         // X-direction offset
    double offset_y,         // Y-direction offset
    double offset_z,         // Z-direction offset
    int num_points
) {
    int point_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_idx >= num_points) return;

    point_coords[point_idx * 3 + 0] += offset_x;
    point_coords[point_idx * 3 + 1] += offset_y;
    point_coords[point_idx * 3 + 2] += offset_z;
}

// Kernel for rotating grid coordinates
__global__ void rotate_grid_kernel(
    double* point_coords,       // [num_points * 3] - Point coordinates to rotate
    const double* rotation_matrix, // [9] - 3x3 rotation matrix (row-major)
    int num_points
) {
    int point_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_idx >= num_points) return;

    double x = point_coords[point_idx * 3 + 0];
    double y = point_coords[point_idx * 3 + 1];
    double z = point_coords[point_idx * 3 + 2];

    // Matrix multiplication: new_coords = rotation_matrix * old_coords
    point_coords[point_idx * 3 + 0] = rotation_matrix[0] * x + rotation_matrix[1] * y + rotation_matrix[2] * z;
    point_coords[point_idx * 3 + 1] = rotation_matrix[3] * x + rotation_matrix[4] * y + rotation_matrix[5] * z;
    point_coords[point_idx * 3 + 2] = rotation_matrix[6] * x + rotation_matrix[7] * y + rotation_matrix[8] * z;
}

// ===== GRID REFINEMENT KERNELS =====

// Kernel for marking cells for refinement based on solution gradients
__global__ void mark_refinement_cells_kernel(
    const double* solution_gradients,  // [num_cells * num_vars] - Solution gradients
    const double* cell_areas,          // [num_cells] - Cell areas
    int* refinement_markers,           // [num_cells] - Output refinement markers
    double gradient_threshold,         // Threshold for refinement
    int num_cells,
    int num_vars
) {
    int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_idx >= num_cells) return;

    double max_gradient = 0.0;
    double area = cell_areas[cell_idx];
    double characteristic_length = sqrt(area);

    // Find maximum normalized gradient magnitude
    for (int var = 0; var < num_vars; var++) {
        double grad_magnitude = fabs(solution_gradients[cell_idx * num_vars + var]);
        double normalized_gradient = grad_magnitude * characteristic_length;
        max_gradient = fmax(max_gradient, normalized_gradient);
    }

    // Mark for refinement if gradient exceeds threshold
    refinement_markers[cell_idx] = (max_gradient > gradient_threshold) ? 1 : 0;
}
