// Linear_Algebra_Cuda_Kernels.cu
// Priority 2: Complete linear algebra CUDA kernels based on Vector operations
#include <cuda_runtime.h>
#include <math.h>

// Vector addition kernel
__global__ void vector_add_kernel(
    const double* a,
    const double* b,
    double* c,
    int n
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    c[idx] = a[idx] + b[idx];
}

// Vector subtraction kernel  
__global__ void vector_subtract_kernel(
    const double* a,
    const double* b,
    double* c,
    int n
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    c[idx] = a[idx] - b[idx];
}

// Vector scalar multiplication kernel
__global__ void vector_scalar_mult_kernel(
    const double* a,
    double scalar,
    double* result,
    int n
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    result[idx] = scalar * a[idx];
}

// Vector dot product kernel with reduction
__global__ void vector_dot_product_kernel(
    const double* a,
    const double* b,
    double* partial_results,
    int n
) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Load data into shared memory
    sdata[tid] = (idx < n) ? a[idx] * b[idx] : 0.0;
    __syncthreads();
    
    // Parallel reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    // Write result for this block
    if (tid == 0) {
        partial_results[blockIdx.x] = sdata[0];
    }
}

// Vector cross product kernel (for 3D vectors)
__global__ void vector_cross_product_kernel(
    const double* a,    // [n/3] - 3D vectors stored as [x1,y1,z1,x2,y2,z2,...]
    const double* b,    // [n/3] - 3D vectors
    double* result,     // [n/3] - Result 3D vectors
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors) return;
    
    // Get vector components
    double ax = a[idx * 3 + 0], ay = a[idx * 3 + 1], az = a[idx * 3 + 2];
    double bx = b[idx * 3 + 0], by = b[idx * 3 + 1], bz = b[idx * 3 + 2];
    
    // Calculate cross product: a × b
    result[idx * 3 + 0] = ay * bz - az * by; // x component
    result[idx * 3 + 1] = az * bx - ax * bz; // y component  
    result[idx * 3 + 2] = ax * by - ay * bx; // z component
}

// Vector magnitude kernel
__global__ void vector_magnitude_kernel(
    const double* vectors,  // [num_vectors * 3] - 3D vectors
    double* magnitudes,     // [num_vectors] - Output magnitudes
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors) return;
    
    double x = vectors[idx * 3 + 0];
    double y = vectors[idx * 3 + 1];
    double z = vectors[idx * 3 + 2];
    
    magnitudes[idx] = sqrt(x*x + y*y + z*z);
}

// Vector normalization kernel
__global__ void vector_normalize_kernel(
    const double* vectors,      // [num_vectors * 3] - Input 3D vectors
    double* unit_vectors,       // [num_vectors * 3] - Output unit vectors
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors) return;
    
    double x = vectors[idx * 3 + 0];
    double y = vectors[idx * 3 + 1];
    double z = vectors[idx * 3 + 2];
    
    double magnitude = sqrt(x*x + y*y + z*z);
    
    if (magnitude > 1e-12) {
        double inv_mag = 1.0 / magnitude;
        unit_vectors[idx * 3 + 0] = x * inv_mag;
        unit_vectors[idx * 3 + 1] = y * inv_mag;
        unit_vectors[idx * 3 + 2] = z * inv_mag;
    } else {
        unit_vectors[idx * 3 + 0] = 0.0;
        unit_vectors[idx * 3 + 1] = 0.0;
        unit_vectors[idx * 3 + 2] = 0.0;
    }
}

// Matrix-vector multiplication kernel (for small matrices)
__global__ void matrix_vector_mult_kernel(
    const double* matrix,   // [num_systems * n * n] - Square matrices
    const double* vector,   // [num_systems * n] - Input vectors
    double* result,         // [num_systems * n] - Output vectors
    int num_systems,
    int n
) {
    int sys_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (sys_idx >= num_systems) return;
    
    // For each row of the matrix
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += matrix[sys_idx * n * n + i * n + j] * vector[sys_idx * n + j];
        }
        result[sys_idx * n + i] = sum;
    }
}

// Vector distance calculation kernel
__global__ void vector_distance_kernel(
    const double* points_a,  // [num_pairs * 3] - First set of 3D points
    const double* points_b,  // [num_pairs * 3] - Second set of 3D points
    double* distances,       // [num_pairs] - Output distances
    int num_pairs
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_pairs) return;
    
    double dx = points_b[idx * 3 + 0] - points_a[idx * 3 + 0];
    double dy = points_b[idx * 3 + 1] - points_a[idx * 3 + 1];
    double dz = points_b[idx * 3 + 2] - points_a[idx * 3 + 2];
    
    distances[idx] = sqrt(dx*dx + dy*dy + dz*dz);
}

// Vector interpolation kernel (linear interpolation)
__global__ void vector_interpolation_kernel(
    const double* vec_a,     // [num_vectors * 3] - First vectors
    const double* vec_b,     // [num_vectors * 3] - Second vectors  
    const double* weights,   // [num_vectors] - Interpolation weights (0 to 1)
    double* result,          // [num_vectors * 3] - Interpolated vectors
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors) return;
    
    double w = weights[idx];
    double inv_w = 1.0 - w;
    
    result[idx * 3 + 0] = inv_w * vec_a[idx * 3 + 0] + w * vec_b[idx * 3 + 0];
    result[idx * 3 + 1] = inv_w * vec_a[idx * 3 + 1] + w * vec_b[idx * 3 + 1];
    result[idx * 3 + 2] = inv_w * vec_a[idx * 3 + 2] + w * vec_b[idx * 3 + 2];
}

// Vector component-wise operations kernel
__global__ void vector_component_ops_kernel(
    const double* vectors,   // [num_vectors * 3] - Input vectors
    double* components,      // [num_vectors * 3] - Output components
    int operation,           // 0=absolute value, 1=square, 2=sqrt
    int num_vectors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_vectors * 3) return;
    
    double value = vectors[idx];
    
    switch (operation) {
        case 0: // Absolute value
            components[idx] = fabs(value);
            break;
        case 1: // Square
            components[idx] = value * value;
            break;
        case 2: // Square root (of absolute value)
            components[idx] = sqrt(fabs(value));
            break;
        default:
            components[idx] = value;
    }
}

// Tensor operations kernel (for 3x3 tensors)
__global__ void tensor_operations_kernel(
    const double* tensors,   // [num_tensors * 9] - 3x3 tensors
    double* results,         // [num_tensors * 9] - Output tensors
    int operation,           // 0=transpose, 1=trace, 2=determinant
    int num_tensors
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_tensors) return;
    
    // Load tensor components
    double T[9];
    for (int i = 0; i < 9; i++) {
        T[i] = tensors[idx * 9 + i];
    }
    
    if (operation == 0) { // Transpose
        results[idx * 9 + 0] = T[0]; results[idx * 9 + 1] = T[3]; results[idx * 9 + 2] = T[6];
        results[idx * 9 + 3] = T[1]; results[idx * 9 + 4] = T[4]; results[idx * 9 + 5] = T[7];
        results[idx * 9 + 6] = T[2]; results[idx * 9 + 7] = T[5]; results[idx * 9 + 8] = T[8];
    } else if (operation == 1) { // Trace (store in first component)
        results[idx * 9 + 0] = T[0] + T[4] + T[8];
    } else if (operation == 2) { // Determinant (store in first component)
        results[idx * 9 + 0] = T[0]*(T[4]*T[8] - T[5]*T[7]) - T[1]*(T[3]*T[8] - T[5]*T[6]) + T[2]*(T[3]*T[7] - T[4]*T[6]);
    }
}

// Advanced reduction kernel with multiple operations
__global__ void advanced_reduction_kernel(
    const double* input,
    double* output,
    int operation,  // 0=sum, 1=max, 2=min, 3=sum_of_squares
    int n
) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Load data with operation
    double value = 0.0;
    if (idx < n) {
        switch (operation) {
            case 0: value = input[idx]; break;                    // Sum
            case 1: value = input[idx]; break;                    // Max
            case 2: value = input[idx]; break;                    // Min
            case 3: value = input[idx] * input[idx]; break;       // Sum of squares
        }
    }
    sdata[tid] = value;
    __syncthreads();
    
    // Reduction with appropriate operation
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            switch (operation) {
                case 0: case 3: sdata[tid] += sdata[tid + s]; break;        // Sum operations
                case 1: sdata[tid] = fmax(sdata[tid], sdata[tid + s]); break; // Max
                case 2: sdata[tid] = fmin(sdata[tid], sdata[tid + s]); break; // Min
            }
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        output[blockIdx.x] = sdata[0];
    }
}
