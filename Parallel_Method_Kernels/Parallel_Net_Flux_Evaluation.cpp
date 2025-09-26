#include <cuda_runtime.h>
#include <iostream>


// Net_Flux.cu
// CUDA implementation of net flux evaluation for all cells

#include "Flux.h"
#include "definitions.h"
#include "Globals.h"
#include "Timestep.h"

__device__ double d_Cells_Net_Flux[MAX_CELLS][4];
__device__ Cell d_Cells[MAX_CELLS];
__device__ int d_Cells_Face_Boundary_Type[MAX_CELLS][4];
__device__ double d_Average_Convective_Flux[4];
__device__ double d_Dissipative_Flux[4];

// Host function to launch CUDA kernels
void Evaluate_Cell_Net_Flux_1O_CUDA_Launch(int No_Physical_Cells, int Dissipation_Type, float (*Cells_Net_Flux)[4])
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (No_Physical_Cells + threadsPerBlock - 1) / threadsPerBlock;

    Evaluate_Cell_Net_Flux_1O_CUDA<<<blocksPerGrid, threadsPerBlock>>>(No_Physical_Cells, Dissipation_Type, Cells_Net_Flux);
    cudaDeviceSynchronize();
}

void Evaluate_Cell_Net_Flux_2O_CUDA_Launch(int No_Physical_Cells, int Dissipation_Type, float (*Cells_Net_Flux)[4])
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (No_Physical_Cells + threadsPerBlock - 1) / threadsPerBlock;

    Evaluate_Cell_Net_Flux_2O_CUDA<<<blocksPerGrid, threadsPerBlock>>>(No_Physical_Cells, Dissipation_Type, Cells_Net_Flux);
    cudaDeviceSynchronize();
}

__device__ void LLF_Device(const int &i, int &j, const int &f)
{
    // Example LLF flux computation (simplified)
    // Fetch state vectors from Cells[i] and Cells[j]
    double rhoL = d_Cells[i].U[0];
    double rhoR = d_Cells[j].U[0];

    double max_lambda = max(abs(d_Cells[i].a), abs(d_Cells[j].a));
    for (int k = 0; k < 4; ++k)
    {
        d_Dissipative_Flux[k] = 0.5 * max_lambda * (d_Cells[j].U[k] - d_Cells[i].U[k]);
    }
}

__device__ void Calculate_Flux_For_All_Faces_Device(int Current_Cell_No, int Dissipation_Type)
{
    int Neighbours[4];
    for (int f = 0; f < 4; ++f)
        Neighbours[f] = d_Cells[Current_Cell_No].Neighbours[f];

    for (int Face_No = 0; Face_No < 4; ++Face_No)
    {
        int Neigh = Neighbours[Face_No];
        int is_boundary = d_Cells_Face_Boundary_Type[Current_Cell_No][Face_No];

        Calculate_Face_Average_Flux_Device(Current_Cell_No, Neigh, Face_No, is_boundary);

        switch (Dissipation_Type)
        {
        case 1:
            LLF_Device(Current_Cell_No, Neigh, Face_No);
            break;
        // Add other dissipation types here
        default:
            LLF_Device(Current_Cell_No, Neigh, Face_No);
        }

        for (int i = 0; i < 4; ++i)
        {
            d_Cells_Net_Flux[Current_Cell_No][i] +=
                d_Average_Convective_Flux[i] - d_Dissipative_Flux[i];
        }
    }
}



__global__ void Net_Flux_Kernel(int No_Cells, int Dissipation_Type)
{
    int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_id >= No_Cells)
        return;

    for (int i = 0; i < 4; ++i)
        d_Cells_Net_Flux[cell_id][i] = 0.0;

    Calculate_Flux_For_All_Faces_Device(cell_id, Dissipation_Type);
    Evaluate_Time_Step_Device(cell_id);
}

void Launch_Flux_Kernel(Cell *h_Cells, int (*h_BC_Type)[4], int No_Cells, int Diss_Type)
{
    Cell *d_Cells_host;
    int (*d_BC_Type)[4];

    cudaMalloc((void **)&d_Cells_host, sizeof(Cell) * No_Cells);
    cudaMemcpyToSymbol(d_Cells, h_Cells, sizeof(Cell) * No_Cells);

    cudaMalloc((void **)&d_BC_Type, sizeof(int) * No_Cells * 4);
    cudaMemcpyToSymbol(d_Cells_Face_Boundary_Type, h_BC_Type, sizeof(int) * No_Cells * 4);

    int threads = 256;
    int blocks = (No_Cells + threads - 1) / threads;

    Net_Flux_Kernel<<<blocks, threads>>>(No_Cells, Diss_Type);
    cudaDeviceSynchronize();

    cudaFree(d_Cells_host);
    cudaFree(d_BC_Type);
}