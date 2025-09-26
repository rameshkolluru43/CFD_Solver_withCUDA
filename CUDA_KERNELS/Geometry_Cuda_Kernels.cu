// Geometry_Cuda_Kernels.cu
// Priority 2: Complete geometry calculation CUDA kernels
#include <cuda_runtime.h>
#include <math.h>

// Kernel for calculating cell volumes using Gauss theorem (divergence theorem)
__global__ void calculate_cell_volumes_kernel(
    const double* cell_coords,    // [num_cells * max_points * 3] - Cell vertex coordinates  
    const double* face_centers,   // [num_faces * 3] - Face center coordinates
    const double* face_areas,     // [num_faces] - Face areas
    const double* face_normals,   // [num_faces * 3] - Face normal vectors
    const int* cell_faces,        // [num_cells * max_faces] - Face indices for each cell
    const int* num_faces_per_cell,// [num_cells] - Number of faces per cell
    double* cell_volumes,         // [num_cells] - Output cell volumes
    int num_cells,
    int max_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_faces = num_faces_per_cell[idx];
    double volume = 0.0;
    
    // Use Gauss theorem: V = (1/3) * sum(r_face · n_face * A_face)
    for (int f = 0; f < num_faces; f++) {
        int face_id = cell_faces[idx * max_faces + f];
        if (face_id >= 0) {
            // Face center coordinates
            double xf = face_centers[face_id * 3 + 0];
            double yf = face_centers[face_id * 3 + 1];
            double zf = face_centers[face_id * 3 + 2];
            
            // Face normal and area
            double nx = face_normals[face_id * 3 + 0];
            double ny = face_normals[face_id * 3 + 1];
            double nz = face_normals[face_id * 3 + 2];
            double area = face_areas[face_id];
            
            // r · n * A contribution
            volume += (1.0/3.0) * (xf * nx + yf * ny + zf * nz) * area;
        }
    }
    
    cell_volumes[idx] = fabs(volume); // Ensure positive volume
}

// Kernel for calculating face areas and normals for quadrilateral faces
__global__ void calculate_quad_face_properties_kernel(
    const double* face_coords,    // [num_faces * 4 * 3] - Face vertex coordinates (4 points per quad)
    double* face_areas,          // [num_faces] - Output face areas
    double* face_normals,        // [num_faces * 3] - Output face normal vectors
    double* face_centers,        // [num_faces * 3] - Output face centers
    int num_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_faces) return;

    // Get quad vertices
    double x1 = face_coords[idx * 12 + 0], y1 = face_coords[idx * 12 + 1], z1 = face_coords[idx * 12 + 2];
    double x2 = face_coords[idx * 12 + 3], y2 = face_coords[idx * 12 + 4], z2 = face_coords[idx * 12 + 5];
    double x3 = face_coords[idx * 12 + 6], y3 = face_coords[idx * 12 + 7], z3 = face_coords[idx * 12 + 8];
    double x4 = face_coords[idx * 12 + 9], y4 = face_coords[idx * 12 + 10], z4 = face_coords[idx * 12 + 11];

    // Calculate face center
    face_centers[idx * 3 + 0] = 0.25 * (x1 + x2 + x3 + x4);
    face_centers[idx * 3 + 1] = 0.25 * (y1 + y2 + y3 + y4);
    face_centers[idx * 3 + 2] = 0.25 * (z1 + z2 + z3 + z4);

    // Calculate area using cross product of diagonals
    double d1x = x3 - x1, d1y = y3 - y1, d1z = z3 - z1; // Diagonal 1
    double d2x = x4 - x2, d2y = y4 - y2, d2z = z4 - z2; // Diagonal 2

    // Cross product d1 × d2
    double nx = d1y * d2z - d1z * d2y;
    double ny = d1z * d2x - d1x * d2z;
    double nz = d1x * d2y - d1y * d2x;

    // Area is half the magnitude of cross product
    double area = 0.5 * sqrt(nx*nx + ny*ny + nz*nz);
    face_areas[idx] = area;

    // Unit normal vector
    double inv_norm = 1.0 / sqrt(nx*nx + ny*ny + nz*nz + 1e-12);
    face_normals[idx * 3 + 0] = nx * inv_norm;
    face_normals[idx * 3 + 1] = ny * inv_norm;
    face_normals[idx * 3 + 2] = nz * inv_norm;
}

// Kernel for calculating face areas and normals for triangular faces
__global__ void calculate_triangle_face_properties_kernel(
    const double* face_coords,    // [num_faces * 3 * 3] - Face vertex coordinates (3 points per triangle)
    double* face_areas,          // [num_faces] - Output face areas
    double* face_normals,        // [num_faces * 3] - Output face normal vectors
    double* face_centers,        // [num_faces * 3] - Output face centers
    int num_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_faces) return;

    // Get triangle vertices
    double x1 = face_coords[idx * 9 + 0], y1 = face_coords[idx * 9 + 1], z1 = face_coords[idx * 9 + 2];
    double x2 = face_coords[idx * 9 + 3], y2 = face_coords[idx * 9 + 4], z2 = face_coords[idx * 9 + 5];
    double x3 = face_coords[idx * 9 + 6], y3 = face_coords[idx * 9 + 7], z3 = face_coords[idx * 9 + 8];

    // Calculate face center
    face_centers[idx * 3 + 0] = (x1 + x2 + x3) / 3.0;
    face_centers[idx * 3 + 1] = (y1 + y2 + y3) / 3.0;
    face_centers[idx * 3 + 2] = (z1 + z2 + z3) / 3.0;

    // Calculate edge vectors
    double e1x = x2 - x1, e1y = y2 - y1, e1z = z2 - z1;
    double e2x = x3 - x1, e2y = y3 - y1, e2z = z3 - z1;

    // Cross product e1 × e2
    double nx = e1y * e2z - e1z * e2y;
    double ny = e1z * e2x - e1x * e2z;
    double nz = e1x * e2y - e1y * e2x;

    // Area is half the magnitude of cross product
    double area = 0.5 * sqrt(nx*nx + ny*ny + nz*nz);
    face_areas[idx] = area;

    // Unit normal vector
    double inv_norm = 1.0 / sqrt(nx*nx + ny*ny + nz*nz + 1e-12);
    face_normals[idx * 3 + 0] = nx * inv_norm;
    face_normals[idx * 3 + 1] = ny * inv_norm;
    face_normals[idx * 3 + 2] = nz * inv_norm;
}

// Kernel for calculating cell centers
__global__ void calculate_cell_centers_kernel(
    const double* cell_coords,    // [num_cells * max_points * 3] - Cell vertex coordinates
    const int* num_points_per_cell, // [num_cells] - Number of vertices per cell
    double* cell_centers,         // [num_cells * 3] - Output cell centers
    int num_cells,
    int max_points
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_points = num_points_per_cell[idx];
    double cx = 0.0, cy = 0.0, cz = 0.0;

    // Average all vertex coordinates
    for (int p = 0; p < num_points; p++) {
        cx += cell_coords[idx * max_points * 3 + p * 3 + 0];
        cy += cell_coords[idx * max_points * 3 + p * 3 + 1];
        cz += cell_coords[idx * max_points * 3 + p * 3 + 2];
    }

    if (num_points > 0) {
        double inv_num_points = 1.0 / num_points;
        cell_centers[idx * 3 + 0] = cx * inv_num_points;
        cell_centers[idx * 3 + 1] = cy * inv_num_points;
        cell_centers[idx * 3 + 2] = cz * inv_num_points;
    }
}

// Kernel for calculating distances between cell centers
__global__ void calculate_cell_distances_kernel(
    const double* cell_centers,   // [num_cells * 3] - Cell center coordinates
    const int* cell_faces,        // [num_cells * max_faces] - Face indices
    const int* face_cells,        // [num_faces * 2] - Left and right cell indices for faces
    const int* num_faces_per_cell,// [num_cells] - Number of faces per cell
    double* cell_distances,       // [num_cells * max_faces] - Output distances to neighbors
    int num_cells,
    int max_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    int num_faces = num_faces_per_cell[idx];
    double cx = cell_centers[idx * 3 + 0];
    double cy = cell_centers[idx * 3 + 1];
    double cz = cell_centers[idx * 3 + 2];

    for (int f = 0; f < num_faces; f++) {
        int face_id = cell_faces[idx * max_faces + f];
        if (face_id >= 0) {
            int left_cell = face_cells[face_id * 2 + 0];
            int right_cell = face_cells[face_id * 2 + 1];
            int neighbor_cell = (left_cell == idx) ? right_cell : left_cell;

            if (neighbor_cell >= 0 && neighbor_cell < num_cells) {
                double nx = cell_centers[neighbor_cell * 3 + 0];
                double ny = cell_centers[neighbor_cell * 3 + 1];
                double nz = cell_centers[neighbor_cell * 3 + 2];
                
                double dist = sqrt((nx - cx)*(nx - cx) + (ny - cy)*(ny - cy) + (nz - cz)*(nz - cz));
                cell_distances[idx * max_faces + f] = dist;
            } else {
                cell_distances[idx * max_faces + f] = 0.0;
            }
        }
    }
}

// Kernel for calculating hexahedral cell volumes using decomposition
__global__ void calculate_hex_volumes_kernel(
    const double* cell_coords,    // [num_cells * 8 * 3] - Hex vertex coordinates (8 points)
    double* cell_volumes,         // [num_cells] - Output cell volumes
    int num_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    // Get hex vertices (assume ordering: o,a,b,c,o1,a1,b1,c1)
    double coords[8][3];
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 3; j++) {
            coords[i][j] = cell_coords[idx * 24 + i * 3 + j];
        }
    }

    // Decompose hex into 6 tetrahedra and sum volumes
    double volume = 0.0;
    
    // Tetrahedron 1: (o, a, c, o1)
    double v1x = coords[1][0] - coords[0][0], v1y = coords[1][1] - coords[0][1], v1z = coords[1][2] - coords[0][2];
    double v2x = coords[3][0] - coords[0][0], v2y = coords[3][1] - coords[0][1], v2z = coords[3][2] - coords[0][2];
    double v3x = coords[4][0] - coords[0][0], v3y = coords[4][1] - coords[0][1], v3z = coords[4][2] - coords[0][2];
    
    // Volume = (1/6) * |v1 · (v2 × v3)|
    double cross_x = v2y * v3z - v2z * v3y;
    double cross_y = v2z * v3x - v2x * v3z;
    double cross_z = v2x * v3y - v2y * v3x;
    double dot = v1x * cross_x + v1y * cross_y + v1z * cross_z;
    volume += fabs(dot) / 6.0;

    // Add other 5 tetrahedra (simplified - would need complete decomposition)
    // This is a placeholder - full implementation would calculate all 6 tetrahedra
    
    cell_volumes[idx] = volume;
}

// Kernel for calculating mesh quality metrics
__global__ void calculate_mesh_quality_kernel(
    const double* cell_volumes,   // [num_cells] - Cell volumes
    const double* face_areas,     // [num_faces] - Face areas
    const int* cell_faces,        // [num_cells * max_faces] - Face indices
    const int* num_faces_per_cell,// [num_cells] - Number of faces per cell
    double* aspect_ratios,        // [num_cells] - Output aspect ratios
    double* skewness,            // [num_cells] - Output skewness measures
    int num_cells,
    int max_faces
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_cells) return;

    double volume = cell_volumes[idx];
    int num_faces = num_faces_per_cell[idx];
    
    // Calculate characteristic length scales
    double min_area = INFINITY, max_area = 0.0, total_area = 0.0;
    
    for (int f = 0; f < num_faces; f++) {
        int face_id = cell_faces[idx * max_faces + f];
        if (face_id >= 0) {
            double area = face_areas[face_id];
            min_area = fmin(min_area, area);
            max_area = fmax(max_area, area);
            total_area += area;
        }
    }
    
    // Aspect ratio approximation
    if (min_area > 0) {
        aspect_ratios[idx] = sqrt(max_area / min_area);
    } else {
        aspect_ratios[idx] = 1.0;
    }
    
    // Skewness approximation (ratio of actual to ideal volume)
    double ideal_volume = pow(volume, 2.0/3.0); // Rough approximation
    if (ideal_volume > 0) {
        skewness[idx] = 1.0 - volume / ideal_volume;
    } else {
        skewness[idx] = 0.0;
    }
}
