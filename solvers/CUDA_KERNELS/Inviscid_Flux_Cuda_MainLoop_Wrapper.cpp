#include <cuda_runtime.h>

#include <iostream>
#include <vector>

#include "Flux.h"
#include "Grid.h"

extern "C" cudaError_t launch_inviscid_flux_1o_mainloop_cuda(const int *d_leaf_cells,
                                                             int num_leaf_cells,
                                                             const double *d_u_cells,
                                                             const double *d_prim_cells,
                                                             const int *d_face_offsets,
                                                             const int *d_neighbours,
                                                             const double *d_face_normals,
                                                             const double *d_face_areas,
                                                             const unsigned char *d_boundary_faces,
                                                             const double *d_cell_areas,
                                                             int dissipation_type,
                                                             int is_movers_1,
                                                             int enable_entropy_fix,
                                                             double cfl,
                                                             double *d_net_flux,
                                                             double *d_cell_dt);

namespace
{
constexpr int kNumEq = 4;
constexpr int kPrimStride = 6;

bool cuda_check(cudaError_t status, const char *where)
{
    if (status == cudaSuccess)
        return true;

    static bool warned = false;
    if (!warned)
    {
        std::cerr << "CUDA inviscid flux disabled at " << where << ": "
                  << cudaGetErrorString(status) << ". Falling back to CPU flux." << std::endl;
        warned = true;
    }
    return false;
}

template <typename T>
bool copy_to_device(T **device_ptr, const std::vector<T> &host, const char *where)
{
    *device_ptr = nullptr;
    if (host.empty())
        return true;

    if (!cuda_check(cudaMalloc(reinterpret_cast<void **>(device_ptr), host.size() * sizeof(T)), where))
        return false;
    return cuda_check(cudaMemcpy(*device_ptr, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice), where);
}

void free_all(int *d_leaf, int *d_offsets, int *d_neigh,
              double *d_u, double *d_prim, double *d_normals, double *d_face_areas,
              double *d_cell_area, unsigned char *d_boundary, double *d_net, double *d_dt)
{
    cudaFree(d_leaf);
    cudaFree(d_offsets);
    cudaFree(d_neigh);
    cudaFree(d_u);
    cudaFree(d_prim);
    cudaFree(d_normals);
    cudaFree(d_face_areas);
    cudaFree(d_cell_area);
    cudaFree(d_boundary);
    cudaFree(d_net);
    cudaFree(d_dt);
}
} // namespace

bool Evaluate_Cell_Net_Flux_1O_CUDA()
{
    if (NUM_FLUX_COMPONENTS != kNumEq)
        return false;
    if (!(Dissipation_Type == 1 || Dissipation_Type == 2 || Dissipation_Type == 4 || Dissipation_Type == 5))
        return false;
    if (Cells.empty() || U_Cells.empty() || Primitive_Cells.empty())
        return false;

    V_I leaf_cells;
    Build_Leaf_Cell_List(leaf_cells);
    if (leaf_cells.empty())
        return true;

    const int n_cells = static_cast<int>(Cells.size());
    std::vector<double> h_u(n_cells * kNumEq, 0.0);
    std::vector<double> h_prim(n_cells * kPrimStride, 0.0);
    std::vector<double> h_cell_area(n_cells, 0.0);
    std::vector<int> h_offsets(n_cells + 1, 0);
    std::vector<int> h_neigh;
    std::vector<double> h_normals;
    std::vector<double> h_face_areas;
    std::vector<unsigned char> h_boundary;

    for (int c = 0; c < n_cells; ++c)
    {
        if (c < static_cast<int>(U_Cells.size()))
        {
            for (int eq = 0; eq < kNumEq && eq < static_cast<int>(U_Cells[c].size()); ++eq)
                h_u[c * kNumEq + eq] = U_Cells[c][eq];
        }

        if (c < static_cast<int>(Primitive_Cells.size()))
        {
            if (Primitive_Cells[c].size() > 0)
                h_prim[c * kPrimStride + 0] = Primitive_Cells[c][0];
            if (Primitive_Cells[c].size() > 1)
                h_prim[c * kPrimStride + 1] = Primitive_Cells[c][1];
            if (Primitive_Cells[c].size() > 2)
                h_prim[c * kPrimStride + 2] = Primitive_Cells[c][2];
            if (Primitive_Cells[c].size() > 4)
                h_prim[c * kPrimStride + 4] = Primitive_Cells[c][4];
            if (Primitive_Cells[c].size() > 5)
                h_prim[c * kPrimStride + 5] = Primitive_Cells[c][5];
        }
        h_cell_area[c] = Cells[c].Area;

        h_offsets[c] = static_cast<int>(h_neigh.size());
        const int n_faces = (Cells[c].numFaces > 0) ? Cells[c].numFaces : static_cast<int>(Cells[c].Face_Areas.size());
        for (int f = 0; f < n_faces; ++f)
        {
            const int neigh = (f < static_cast<int>(Cells[c].Neighbours.size())) ? Cells[c].Neighbours[f] : -1;
            h_neigh.push_back(neigh);
            h_face_areas.push_back((f < static_cast<int>(Cells[c].Face_Areas.size())) ? Cells[c].Face_Areas[f] : 0.0);
            h_normals.push_back((2 * f + 0 < static_cast<int>(Cells[c].Face_Normals.size())) ? Cells[c].Face_Normals[2 * f + 0] : 0.0);
            h_normals.push_back((2 * f + 1 < static_cast<int>(Cells[c].Face_Normals.size())) ? Cells[c].Face_Normals[2 * f + 1] : 0.0);
            const bool is_boundary = c < static_cast<int>(Cells_Face_Boundary_Type.size()) &&
                                     f < static_cast<int>(Cells_Face_Boundary_Type[c].size()) &&
                                     Cells_Face_Boundary_Type[c][f];
            h_boundary.push_back(static_cast<unsigned char>(is_boundary ? 1 : 0));
        }
    }
    h_offsets[n_cells] = static_cast<int>(h_neigh.size());

    std::vector<double> h_net(n_cells * kNumEq, 0.0);
    std::vector<double> h_dt(n_cells, 0.0);

    int *d_leaf = nullptr;
    int *d_offsets = nullptr;
    int *d_neigh = nullptr;
    double *d_u = nullptr;
    double *d_prim = nullptr;
    double *d_normals = nullptr;
    double *d_face_areas = nullptr;
    double *d_cell_area = nullptr;
    unsigned char *d_boundary = nullptr;
    double *d_net = nullptr;
    double *d_dt = nullptr;

    bool ok = copy_to_device(&d_leaf, leaf_cells, "leaf cells") &&
              copy_to_device(&d_offsets, h_offsets, "face offsets") &&
              copy_to_device(&d_neigh, h_neigh, "neighbours") &&
              copy_to_device(&d_u, h_u, "U cells") &&
              copy_to_device(&d_prim, h_prim, "primitive cells") &&
              copy_to_device(&d_normals, h_normals, "face normals") &&
              copy_to_device(&d_face_areas, h_face_areas, "face areas") &&
              copy_to_device(&d_cell_area, h_cell_area, "cell areas") &&
              copy_to_device(&d_boundary, h_boundary, "boundary faces") &&
              cuda_check(cudaMalloc(reinterpret_cast<void **>(&d_net), h_net.size() * sizeof(double)), "net flux allocation") &&
              cuda_check(cudaMalloc(reinterpret_cast<void **>(&d_dt), h_dt.size() * sizeof(double)), "dt allocation");

    if (ok)
    {
        cudaMemset(d_net, 0, h_net.size() * sizeof(double));
        cudaMemset(d_dt, 0, h_dt.size() * sizeof(double));
        ok = cuda_check(launch_inviscid_flux_1o_mainloop_cuda(d_leaf, static_cast<int>(leaf_cells.size()),
                                                              d_u, d_prim, d_offsets, d_neigh, d_normals,
                                                              d_face_areas, d_boundary, d_cell_area,
                                                              Dissipation_Type, Is_MOVERS_1 ? 1 : 0,
                                                              Enable_Entropy_Fix ? 1 : 0, CFL, d_net, d_dt),
                        "kernel launch") &&
             cuda_check(cudaDeviceSynchronize(), "kernel synchronize") &&
             cuda_check(cudaMemcpy(h_net.data(), d_net, h_net.size() * sizeof(double), cudaMemcpyDeviceToHost), "net flux copy") &&
             cuda_check(cudaMemcpy(h_dt.data(), d_dt, h_dt.size() * sizeof(double), cudaMemcpyDeviceToHost), "dt copy");
    }

    free_all(d_leaf, d_offsets, d_neigh, d_u, d_prim, d_normals, d_face_areas, d_cell_area, d_boundary, d_net, d_dt);

    if (!ok)
        return false;

    for (int cell : leaf_cells)
    {
        if (cell < 0 || cell >= n_cells || cell >= static_cast<int>(Cells_Net_Flux.size()))
            continue;
        for (int eq = 0; eq < kNumEq && eq < static_cast<int>(Cells_Net_Flux[cell].size()); ++eq)
            Cells_Net_Flux[cell][eq] = h_net[cell * kNumEq + eq];
        Cells[cell].del_t = h_dt[cell];
    }
    return true;
}
