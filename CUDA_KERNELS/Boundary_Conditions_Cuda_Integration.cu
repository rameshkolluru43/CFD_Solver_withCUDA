/**
 * @file Boundary_Conditions_Cuda_Integration.cu
 * @brief Viscous no-slip wall BC on ghost cells (matches Viscous_Wall_Boundary_Condition in src/Wall_Boundary_Conditions.cpp).
 */

#include <cuda_runtime.h>
#include <vector>
#include <cstdio>

#include "../include/definitions.h"
#include "../include/Globals.h"
#include "../include/Boundary_Conditions.h"
#include "../include/Primitive_Computational.h"
#include "../include/Utilities.h"
#include "Boundary_Conditions_Cuda_Integration.h"

namespace
{
constexpr int kNv = 4;

__global__ void viscous_wall_ghost_kernel(
    double *U,
    const int *wall_list,
    int n_triplets)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_triplets)
        return;

    const int cell = wall_list[idx * 3 + 0];
    const int ghost = wall_list[idx * 3 + 2];

    const int cb = cell * kNv;
    const int gb = ghost * kNv;

    U[gb + 0] = U[cb + 0];
    U[gb + 1] = -U[cb + 1];
    U[gb + 2] = -U[cb + 2];
    U[gb + 3] = U[cb + 3];
}

struct WallBcBuffers
{
    int triplet_cap = 0;
    int cell_cap = 0;
    double *d_U = nullptr;
    int *d_wall_list = nullptr;

    void free_all()
    {
        if (d_U)
            cudaFree(d_U);
        if (d_wall_list)
            cudaFree(d_wall_list);
        *this = WallBcBuffers{};
    }
};

WallBcBuffers g_wall_bc;

bool ensure_wall_buffers(int n_cells, int n_triplets)
{
    if (n_cells <= g_wall_bc.cell_cap && n_triplets <= g_wall_bc.triplet_cap)
        return true;

    g_wall_bc.free_all();
    g_wall_bc.cell_cap = n_cells;
    g_wall_bc.triplet_cap = n_triplets;

    if (cudaMalloc(&g_wall_bc.d_U, static_cast<size_t>(n_cells) * kNv * sizeof(double)) != cudaSuccess)
        return false;
    if (cudaMalloc(&g_wall_bc.d_wall_list, static_cast<size_t>(n_triplets) * 3 * sizeof(int)) != cudaSuccess)
    {
        cudaFree(g_wall_bc.d_U);
        g_wall_bc.d_U = nullptr;
        return false;
    }
    return true;
}

} // namespace

bool Apply_Viscous_Wall_Boundary_Condition_CUDA()
{
    if (!Is_Viscous_Wall || Wall_Cells_List.empty())
        return false;

    const int n_triplets = static_cast<int>(Wall_Cells_List.size() / 3);
    if (n_triplets <= 0)
        return false;

    const int n_cells = static_cast<int>(U_Cells.size());
    if (n_cells <= 0)
        return false;

    if (!ensure_wall_buffers(n_cells, n_triplets))
    {
        fprintf(stderr, "Viscous wall BC CUDA: buffer allocation failed\n");
        return false;
    }

    std::vector<double> h_U(static_cast<size_t>(n_cells) * kNv, 0.0);
    for (int c = 0; c < n_cells; ++c)
    {
        if (c < static_cast<int>(U_Cells.size()) && U_Cells[c].size() >= static_cast<size_t>(kNv))
        {
            for (int v = 0; v < kNv; ++v)
                h_U[c * kNv + v] = U_Cells[c][v];
        }
    }

  cudaMemcpy(g_wall_bc.d_U, h_U.data(), h_U.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(g_wall_bc.d_wall_list, Wall_Cells_List.data(),
             static_cast<size_t>(n_triplets) * 3 * sizeof(int), cudaMemcpyHostToDevice);

    const int threads = 256;
    const int blocks = (n_triplets + threads - 1) / threads;
    viscous_wall_ghost_kernel<<<blocks, threads>>>(g_wall_bc.d_U, g_wall_bc.d_wall_list, n_triplets);

    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Viscous wall BC CUDA kernel failed: %s\n", cudaGetErrorString(err));
        return false;
    }

    cudaMemcpy(h_U.data(), g_wall_bc.d_U, h_U.size() * sizeof(double), cudaMemcpyDeviceToHost);

    for (int c = 0; c < n_cells; ++c)
    {
        if (c < static_cast<int>(U_Cells.size()) && U_Cells[c].size() >= static_cast<size_t>(kNv))
        {
            for (int v = 0; v < kNv; ++v)
                U_Cells[c][v] = h_U[c * kNv + v];
        }
    }

    for (int t = 0; t < n_triplets; ++t)
    {
        const int ghost = Wall_Cells_List[t * 3 + 2];
        if (ghost < 0 || ghost >= n_cells)
            continue;
        Calculate_Primitive_Variables(ghost, U_Cells[ghost]);
        if (ghost < static_cast<int>(Primitive_Cells.size()))
        {
            Vector_Reset(Primitive_Cells[ghost]);
            for (unsigned int i = 0; i < Global_Primitive.size(); ++i)
                Primitive_Cells[ghost][i] = Global_Primitive[i];
        }
    }

    return true;
}
