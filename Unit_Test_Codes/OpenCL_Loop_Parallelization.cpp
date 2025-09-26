#include <CL/cl.h>
#include <iostream>
#include <vector>

const char* kernelSource = R"(
__kernel void vector_add(__global const float* A, __global const float* B, __global float* C, const int n) {
    int i = get_global_id(0);
    if (i < n) {
        C[i] = A[i] + B[i];
    }
}
)";

int main() {
    const int n = 1024;
    std::vector<float> A(n, 1.0f), B(n, 2.0f), C(n);

    // Step 1: Set up OpenCL environment
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    cl_kernel kernel;

    cl_int err = clGetPlatformIDs(1, &platform, nullptr);
    err |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
    context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
    queue = clCreateCommandQueue(context, device, 0, &err);

    // Step 2: Compile the kernel
    program = clCreateProgramWithSource(context, 1, &kernelSource, nullptr, &err);
    err = clBuildProgram(program, 1, &device, nullptr, nullptr, nullptr);
    kernel = clCreateKernel(program, "vector_add", &err);

    // Step 3: Allocate device memory and transfer data
    cl_mem bufferA = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, n * sizeof(float), A.data(), &err);
    cl_mem bufferB = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, n * sizeof(float), B.data(), &err);
    cl_mem bufferC = clCreateBuffer(context, CL_MEM_WRITE_ONLY, n * sizeof(float), nullptr, &err);

    // Step 4: Set kernel arguments and launch
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufferA);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufferB);
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &bufferC);
    err |= clSetKernelArg(kernel, 3, sizeof(int), &n);

    size_t globalWorkSize = n;
    err = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &globalWorkSize, nullptr, 0, nullptr, nullptr);

    // Step 5: Retrieve results
    err = clEnqueueReadBuffer(queue, bufferC, CL_TRUE, 0, n * sizeof(float), C.data(), 0, nullptr, nullptr);

    // Verify results
    for (int i = 0; i < n; ++i) {
        if (C[i] != 3.0f) {
            std::cerr << "Error at index " << i << ": " << C[i] << std::endl;
            return -1;
        }
    }
    std::cout << "Computation successful!" << std::endl;

    // Step 6: Clean up
    clReleaseMemObject(bufferA);
    clReleaseMemObject(bufferB);
    clReleaseMemObject(bufferC);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);

    return 0;
}
