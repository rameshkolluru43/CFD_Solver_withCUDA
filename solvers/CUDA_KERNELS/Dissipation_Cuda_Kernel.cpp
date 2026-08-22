#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>

__device__ double Sign(double value) {
    return (value >= 0) ? 1.0 : -1.0;
}

__global__ void Condition_For_MOVERS_CUDA(double* d_U, double* d_F, double* L_Max, double* L_Min, double* Alpha, int n) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= n) return;

    double epsilon = 1e-10;
    double local_Alpha = 0.0;
    
    if (fabs(d_F[idx]) < epsilon && fabs(d_U[idx]) > epsilon) {
        local_Alpha = 0.0;
    } else if (fabs(d_F[idx]) < epsilon && fabs(d_U[idx]) < epsilon) {
        local_Alpha = L_Min[idx];
    } else if (fabs(d_F[idx]) > epsilon && fabs(d_U[idx]) > epsilon) {
        local_Alpha = fabs(d_F[idx] / d_U[idx]);
        if (local_Alpha >= L_Max[idx]) {
            local_Alpha = Sign(local_Alpha) * L_Max[idx];
        } else if (local_Alpha <= L_Min[idx]) {
            local_Alpha = Sign(local_Alpha) * L_Min[idx];
        }
    } else {
        local_Alpha = L_Min[idx];
    }

    Alpha[idx] = local_Alpha;
}

int main() {
    const int n = 1024;
    double *h_d_U, *h_d_F, *h_L_Max, *h_L_Min, *h_Alpha;
    double *d_d_U, *d_d_F, *d_L_Max, *d_L_Min, *d_Alpha;

    // Allocate memory on host
    h_d_U = (double*)malloc(n * sizeof(double));
    h_d_F = (double*)malloc(n * sizeof(double));
    h_L_Max = (double*)malloc(n * sizeof(double));
    h_L_Min = (double*)malloc(n * sizeof(double));
    h_Alpha = (double*)malloc(n * sizeof(double));

    // Initialize data
    for (int i = 0; i < n; i++) {
        h_d_U[i] = 1.0; // Replace with real values
        h_d_F[i] = 0.5; // Replace with real values
        h_L_Max[i] = 2.0; // Replace with real values
        h_L_Min[i] = 0.1; // Replace with real values
    }

    // Allocate memory on device
    cudaMalloc((void**)&d_d_U, n * sizeof(double));
    cudaMalloc((void**)&d_d_F, n * sizeof(double));
    cudaMalloc((void**)&d_L_Max, n * sizeof(double));
    cudaMalloc((void**)&d_L_Min, n * sizeof(double));
    cudaMalloc((void**)&d_Alpha, n * sizeof(double));

    // Copy data from host to device
    cudaMemcpy(d_d_U, h_d_U, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_d_F, h_d_F, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L_Max, h_L_Max, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L_Min, h_L_Min, n * sizeof(double), cudaMemcpyHostToDevice);

    // Launch kernel
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;
    Condition_For_MOVERS_CUDA<<<numBlocks, blockSize>>>(d_d_U, d_d_F, d_L_Max, d_L_Min, d_Alpha, n);

    // Copy results from device to host
    cudaMemcpy(h_Alpha, d_Alpha, n * sizeof(double), cudaMemcpyDeviceToHost);

    // Free memory
    cudaFree(d_d_U);
    cudaFree(d_d_F);
    cudaFree(d_L_Max);
    cudaFree(d_L_Min);
    cudaFree(d_Alpha);
    free(h_d_U);
    free(h_d_F);
    free(h_L_Max);
    free(h_L_Min);
    free(h_Alpha);

    return 0;
}
