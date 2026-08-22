#include "../Basic_Files/definitions.h"
#include <cuda_runtime.h>

// CUDA kernel for Explicit_Method
/**
 * CUDA kernel for updating cell values using the explicit method.
 *
 * @param d_Cells_Inv_Area Pointer to device array containing inverse cell areas.
 * @param d_Cells_Net_Flux Pointer to device array containing net flux values for each cell.
 * @param d_Cells_Viscous_Flux Pointer to device array containing viscous flux values for each cell.
 * @param d_Cells_DelU Pointer to device array where the updated cell values will be stored.
 * @param Min_dt Minimum time step size for the simulation.
 * @param Is_Viscous_Wall Boolean flag indicating whether viscous wall conditions are applied.
 * @param No_Physical_Cells Total number of physical cells in the simulation.
 */
#define NUM_COMPONENTS 4 // Configurable parameter for the number of components per cell

__global__ void Explicit_Method_Kernel(
    double *d_Cells_Inv_Area, double *d_Cells_Net_Flux, double *d_Cells_Viscous_Flux,
    double *d_Cells_DelU, double Min_dt, bool Is_Viscous_Wall, size_t No_Physical_Cells)
{
    size_t Cell_Index = blockIdx.x * blockDim.x + threadIdx.x;

    if (Cell_Index < No_Physical_Cells)
    {
        double inv_Area = d_Cells_Inv_Area[Cell_Index];

        if (Is_Viscous_Wall)
        {
            // Apply preconditioned viscous and inviscid fluxes when viscous fluxes are enabled
            for (int i = 0; i < NUM_COMPONENTS; i++)
            {
                int index = Cell_Index * NUM_COMPONENTS + i;
                d_Cells_DelU[index] = -Min_dt * inv_Area *
                                      (d_Cells_Net_Flux[index] - d_Cells_Viscous_Flux[index]);
            }
        }
        else
        {
            // Apply preconditioned inviscid fluxes only
            for (int i = 0; i < NUM_COMPONENTS; i++)
            {
                int index = Cell_Index * NUM_COMPONENTS + i;
                d_Cells_DelU[index] = -Min_dt * inv_Area * d_Cells_Net_Flux[index];
            }
        }
    }
}

void PrecomputeFluxes()
{
    // Check if the Weighted Essentially Non-Oscillatory (WENO) method is enabled.
    // WENO is a high-order numerical method used for solving hyperbolic partial differential equations.
    if (Is_WENO)
    {
        Evaluate_Cell_Net_Flux_WENO();
    }
    else
    {
        // If WENO is not enabled, check if the second-order method is used.
        // The second-order method provides a balance between accuracy and computational cost.
        if (Is_Second_Order)
            Evaluate_Cell_Net_Flux_2O();
        else
            // If neither WENO nor second-order methods are enabled, fall back to the first-order method.
            // The first-order method is less accurate but computationally cheaper.
            Evaluate_Cell_Net_Flux_1O();
    }
}

void AllocateAndCopyToDevice(double *&d_Cells_Inv_Area, double *&d_Cells_Net_Flux, double *&d_Cells_Viscous_Flux, double *&d_Cells_DelU, const vector<double> &h_inv_Area)
{
    if (No_Physical_Cells == 0)
    {
        fprintf(stderr, "Error: No_Physical_Cells is zero. Memory allocation aborted.\n");
        exit(EXIT_FAILURE);
    }

    if (cudaMalloc(&d_Cells_Inv_Area, No_Physical_Cells * sizeof(double)) != cudaSuccess ||
        cudaMalloc(&d_Cells_Net_Flux, No_Physical_Cells * 4 * sizeof(double)) != cudaSuccess ||
        cudaMalloc(&d_Cells_Viscous_Flux, No_Physical_Cells * 4 * sizeof(double)) != cudaSuccess)
    {
        throw std::runtime_error("cudaMalloc failed: Unable to allocate device memory.");
    }

    if (No_Physical_Cells > 0)
    {
        cudaMemcpy(d_Cells_Inv_Area, h_inv_Area.data(), No_Physical_Cells * sizeof(double), cudaMemcpyHostToDevice);
    }
    else
    {
        throw std::invalid_argument("Error: No_Physical_Cells is zero or invalid.");
    }

    cudaMemcpy(d_Cells_Inv_Area, h_inv_Area.data(), No_Physical_Cells * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Cells_Net_Flux, Cells_Net_Flux, No_Physical_Cells * 4 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Cells_Viscous_Flux, Cells_Viscous_Flux, No_Physical_Cells * 4 * sizeof(double), cudaMemcpyHostToDevice);
}

void LaunchKernel(double *d_Cells_Inv_Area, double *d_Cells_Net_Flux, double *d_Cells_Viscous_Flux, double *d_Cells_DelU)
{
    size_t threadsPerBlock = 256;
    size_t maxBlocksPerGrid = (size_t)cudaDeviceProp().maxGridSize[0];
    size_t blocksPerGrid = (No_Physical_Cells + threadsPerBlock - 1) / threadsPerBlock;
    if (blocksPerGrid > maxBlocksPerGrid)
    {
        fprintf(stderr, "Error: blocksPerGrid exceeds the maximum supported grid size.\n");
        exit(EXIT_FAILURE);
    }
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    ComputeTimeSteps();
    Min_dt = get_Min_dt();
    if (Min_dt == 0.0)
    {
        fprintf(stderr, "Error: Min_dt is zero. Cannot proceed with the simulation.\n");
        exit(EXIT_FAILURE);
    }
    if (Min_dt < 0.0)
    {
        fprintf(stderr, "Error: Min_dt is negative. Cannot proceed with the simulation.\n");
        exit(EXIT_FAILURE);
    }

    Explicit_Method_Kernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_Cells_Inv_Area, d_Cells_Net_Flux, d_Cells_Viscous_Flux, d_Cells_DelU,
        Min_dt, Is_Viscous_Wall, No_Physical_Cells);

    cudaDeviceSynchronize();

    if (d_Cells_Inv_Area)
        cudaFree(d_Cells_Inv_Area);
    if (d_Cells_Net_Flux)
        cudaFree(d_Cells_Net_Flux);
    if (d_Cells_Viscous_Flux)
        cudaFree(d_Cells_Viscous_Flux);
    if (d_Cells_DelU)
        cudaFree(d_Cells_DelU);
    cudaFree(d_Cells_Viscous_Flux);
    cudaFree(d_Cells_DelU);
}

void Explicit_Method()
{
    PrecomputeFluxes();
    vector<double> h_inv_Area(No_Physical_Cells);
    std::transform(Cells.begin(), Cells.begin() + No_Physical_Cells, h_inv_Area.begin(),
                   [](const auto &cell)
                   { return cell.Inv_Area; });

    double *d_Cells_Inv_Area, *d_Cells_Net_Flux, *d_Cells_Viscous_Flux, *d_Cells_DelU;

    AllocateAndCopyToDevice(d_Cells_Inv_Area, d_Cells_Net_Flux, d_Cells_Viscous_Flux, d_Cells_DelU, h_inv_Area);

    LaunchKernel(d_Cells_Inv_Area, d_Cells_Net_Flux, d_Cells_Viscous_Flux, d_Cells_DelU);

    FreeDeviceMemory(d_Cells_Inv_Area, d_Cells_Net_Flux, d_Cells_Viscous_Flux, d_Cells_DelU);
}