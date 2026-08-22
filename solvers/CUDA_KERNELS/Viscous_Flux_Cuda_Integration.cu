/**
 * @file Viscous_Flux_Cuda_Integration.cu
 * @brief 2D viscous flux on GPU: Green–Gauss + face-normal correction (α-damped),
 *        μ-scaled heat conductivity. Used by host Evaluate_Viscous_Fluxes when
 *        the resident NS stepper is not active.
 */

#include <cuda_runtime.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>

#include "../include/definitions.h"
#include "../include/Globals.h"
#include "Viscous_Flux_Cuda_Integration.h"

namespace
{
constexpr int kMaxFaces = 8;
constexpr int kNv = 4;
constexpr int kPrimStride = 9; // rho,u,v,T,P,c,pad,pad,mu

__device__ double prim_at(const double *prim, int cell, int comp)
{
    return prim[cell * kPrimStride + comp];
}

__device__ double d_sutherland_mu(double T)
{
    const double term1 = 110.4 / 288.15;
    if (!(T > 1e-12) || !isfinite(T))
        return 1.0;
    const double term2 = T + term1;
    return pow(T, 1.5) * ((1.0 + term1) / term2);
}

__device__ void green_gauss_grad_cell(
    int cell,
    int grad_comp,
    const int *num_faces,
    const int *neighbours,
    const double *face_normals,
    const double *face_areas,
    double inv_area,
    const double *prim,
    int n_cells,
    double &gx,
    double &gy)
{
    gx = 0.0;
    gy = 0.0;
    if (inv_area <= 0.0)
        return;

    const double phi_c = prim_at(prim, cell, grad_comp);
    const int nf = num_faces[cell];

    for (int f = 0; f < nf && f < kMaxFaces; ++f)
    {
        const int neigh = neighbours[cell * kMaxFaces + f];
        double phi_n = phi_c;
        if (neigh >= 0 && neigh < n_cells)
            phi_n = prim_at(prim, neigh, grad_comp);

        const int nb = (cell * kMaxFaces + f) * 2;
        const double nx = face_normals[nb + 0];
        const double ny = face_normals[nb + 1];
        const double dl = face_areas[cell * kMaxFaces + f];
        const double phi_face = 0.5 * (phi_c + phi_n);
        gx += phi_face * nx * dl;
        gy += phi_face * ny * dl;
    }
    gx *= inv_area;
    gy *= inv_area;
}

__device__ void face_grad_corrected(
    int cell,
    int neigh,
    int grad_comp,
    const int *num_faces,
    const int *neighbours,
    const double *face_normals,
    const double *face_areas,
    const double *inv_area,
    const double *prim,
    const double *centers,
    int n_cells,
    double alpha,
    double &gx,
    double &gy)
{
    double gxc = 0.0, gyc = 0.0, gxn = 0.0, gyn = 0.0;
    green_gauss_grad_cell(cell, grad_comp, num_faces, neighbours, face_normals, face_areas,
                          inv_area[cell], prim, n_cells, gxc, gyc);
    if (neigh >= 0 && neigh < n_cells)
        green_gauss_grad_cell(neigh, grad_comp, num_faces, neighbours, face_normals, face_areas,
                              inv_area[neigh], prim, n_cells, gxn, gyn);
    else
    {
        gxn = gxc;
        gyn = gyc;
    }
    gx = 0.5 * (gxc + gxn);
    gy = 0.5 * (gyc + gyn);

    if (neigh < 0 || neigh >= n_cells || centers == nullptr)
        return;

    const double dx = centers[neigh * 2 + 0] - centers[cell * 2 + 0];
    const double dy = centers[neigh * 2 + 1] - centers[cell * 2 + 1];
    const double ds2 = dx * dx + dy * dy;
    if (!(ds2 > 1e-30))
        return;
    const double ds = sqrt(ds2);
    const double ex = dx / ds;
    const double ey = dy / ds;
    const double dphi_exact = (prim_at(prim, neigh, grad_comp) - prim_at(prim, cell, grad_comp)) / ds;
    const double dphi_approx = gx * ex + gy * ey;
    const double corr = alpha * (dphi_exact - dphi_approx);
    gx += corr * ex;
    gy += corr * ey;
}

__device__ void viscous_flux_on_face(
    int cell,
    int neigh,
    int face,
    const int *num_faces,
    const int *neighbours,
    const double *face_normals,
    const double *face_areas,
    const double *inv_area,
    const double *prim,
    const double *centers,
    int n_cells,
    double inv_re,
    double k1_base,
    double alpha,
    double &f0,
    double &f1,
    double &f2,
    double &f3)
{
    f0 = f1 = f2 = f3 = 0.0;
    if (neigh < 0 || neigh >= n_cells)
        return;

    const double nx = face_normals[(cell * kMaxFaces + face) * 2 + 0];
    const double ny = face_normals[(cell * kMaxFaces + face) * 2 + 1];
    const double dl = face_areas[cell * kMaxFaces + face];
    if (dl <= 0.0)
        return;

    const double u_L = prim_at(prim, cell, 1);
    const double v_L = prim_at(prim, cell, 2);
    const double u_R = prim_at(prim, neigh, 1);
    const double v_R = prim_at(prim, neigh, 2);
    double mu = 0.5 * (prim_at(prim, cell, 8) + prim_at(prim, neigh, 8));
    if (!(mu > 0.0) || !isfinite(mu))
    {
        const double T_L = prim_at(prim, cell, 3);
        const double T_R = prim_at(prim, neigh, 3);
        mu = 0.5 * (d_sutherland_mu(T_L) + d_sutherland_mu(T_R));
    }

    const double v1 = 0.5 * (u_L + u_R);
    const double v2 = 0.5 * (v_L + v_R);

    double u11, u12, u21, u22, tgx, tgy;
    face_grad_corrected(cell, neigh, 1, num_faces, neighbours, face_normals, face_areas, inv_area,
                        prim, centers, n_cells, alpha, u11, u12);
    face_grad_corrected(cell, neigh, 2, num_faces, neighbours, face_normals, face_areas, inv_area,
                        prim, centers, n_cells, alpha, u21, u22);
    face_grad_corrected(cell, neigh, 3, num_faces, neighbours, face_normals, face_areas, inv_area,
                        prim, centers, n_cells, alpha, tgx, tgy);

    const double T11 = (2.0 / 3.0) * mu * inv_re * (2.0 * u11 - u22);
    const double T12 = mu * inv_re * (u12 + u21);
    const double T21 = T12;
    const double T22 = (2.0 / 3.0) * mu * inv_re * (2.0 * u22 - u11);
    const double k1 = k1_base * mu;
    const double Qx = k1 * tgx;
    const double Qy = k1 * tgy;

    f0 = 0.0;
    f1 = (T11 * nx + T21 * ny) * dl;
    f2 = (T12 * nx + T22 * ny) * dl;
    f3 = ((T11 * v1 + T12 * v2 + Qx) * nx + (T21 * v1 + T22 * v2 + Qy) * ny) * dl;
}

__global__ void evaluate_viscous_flux_kernel(
    const int *num_faces,
    const int *neighbours,
    const double *face_normals,
    const double *face_areas,
    const double *inv_area,
    const double *prim,
    const double *centers,
    double *viscous_flux,
    int n_cells,
    int n_physical,
    double inv_re,
    double k1,
    double alpha)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_physical)
        return;

    for (int v = 0; v < kNv; ++v)
        viscous_flux[cell * kNv + v] = 0.0;

    const int nf = num_faces[cell];
    for (int face = 0; face < nf && face < kMaxFaces; ++face)
    {
        const int neigh = neighbours[cell * kMaxFaces + face];
        double f0, f1, f2, f3;
        viscous_flux_on_face(cell, neigh, face, num_faces, neighbours, face_normals, face_areas,
                             inv_area, prim, centers, n_cells, inv_re, k1, alpha, f0, f1, f2, f3);
        viscous_flux[cell * kNv + 0] += f0;
        viscous_flux[cell * kNv + 1] += f1;
        viscous_flux[cell * kNv + 2] += f2;
        viscous_flux[cell * kNv + 3] += f3;
    }
}

struct ViscousCudaBuffers
{
    int cap = 0;
    int *d_num_faces = nullptr;
    int *d_neighbours = nullptr;
    double *d_face_normals = nullptr;
    double *d_face_areas = nullptr;
    double *d_inv_area = nullptr;
    double *d_prim = nullptr;
    double *d_viscous = nullptr;
    double *d_centers = nullptr;

    void free_all()
    {
        if (d_num_faces)
            cudaFree(d_num_faces);
        if (d_neighbours)
            cudaFree(d_neighbours);
        if (d_face_normals)
            cudaFree(d_face_normals);
        if (d_face_areas)
            cudaFree(d_face_areas);
        if (d_inv_area)
            cudaFree(d_inv_area);
        if (d_prim)
            cudaFree(d_prim);
        if (d_viscous)
            cudaFree(d_viscous);
        if (d_centers)
            cudaFree(d_centers);
        *this = ViscousCudaBuffers{};
    }
};

ViscousCudaBuffers g_visc;

bool ensure_visc_buffers(int n_cells)
{
    if (n_cells <= g_visc.cap)
        return true;
    g_visc.free_all();
    g_visc.cap = n_cells;

    auto chk = [](cudaError_t e, const char *msg) {
        if (e != cudaSuccess)
        {
            fprintf(stderr, "Viscous flux CUDA: %s (%s)\n", msg, cudaGetErrorString(e));
            return false;
        }
        return true;
    };

    if (!chk(cudaMalloc(&g_visc.d_num_faces, n_cells * sizeof(int)), "num_faces"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_neighbours, n_cells * kMaxFaces * sizeof(int)), "neighbours"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_face_normals, n_cells * kMaxFaces * 2 * sizeof(double)), "normals"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_face_areas, n_cells * kMaxFaces * sizeof(double)), "areas"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_inv_area, n_cells * sizeof(double)), "inv_area"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_prim, n_cells * kPrimStride * sizeof(double)), "prim"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_viscous, n_cells * kNv * sizeof(double)), "viscous"))
        return false;
    if (!chk(cudaMalloc(&g_visc.d_centers, n_cells * 2 * sizeof(double)), "centers"))
        return false;
    return true;
}

} // namespace

bool Evaluate_Viscous_Fluxes_CUDA()
{
    if (!Is_Viscous_Wall)
        return false;

    const int n_cells = static_cast<int>(Cells.size());
    if (n_cells <= 0 || No_Physical_Cells <= 0)
        return false;

    if (!ensure_visc_buffers(n_cells))
        return false;

    std::vector<int> h_num_faces(n_cells, 0);
    std::vector<int> h_neighbours(n_cells * kMaxFaces, -1);
    std::vector<double> h_normals(n_cells * kMaxFaces * 2, 0.0);
    std::vector<double> h_areas(n_cells * kMaxFaces, 0.0);
    std::vector<double> h_inv_area(n_cells, 0.0);
    std::vector<double> h_prim(n_cells * kPrimStride, 0.0);
    std::vector<double> h_centers(n_cells * 2, 0.0);

    for (int c = 0; c < n_cells; ++c)
    {
        const int nf = (Cells[c].numFaces > 0)
                           ? Cells[c].numFaces
                           : static_cast<int>(Cells[c].Face_Areas.size());
        h_num_faces[c] = std::min(nf, kMaxFaces);
        h_inv_area[c] = Cells[c].Inv_Area;
        if (Cells[c].Cell_Center.size() >= 2)
        {
            h_centers[c * 2 + 0] = Cells[c].Cell_Center[0];
            h_centers[c * 2 + 1] = Cells[c].Cell_Center[1];
        }

        for (int f = 0; f < h_num_faces[c]; ++f)
        {
            if (f < static_cast<int>(Cells[c].Neighbours.size()))
                h_neighbours[c * kMaxFaces + f] = Cells[c].Neighbours[f];
            const int nb = (c * kMaxFaces + f) * 2;
            if (f * 2 + 1 < static_cast<int>(Cells[c].Face_Normals.size()))
            {
                h_normals[nb + 0] = Cells[c].Face_Normals[f * 2 + 0];
                h_normals[nb + 1] = Cells[c].Face_Normals[f * 2 + 1];
            }
            if (f < static_cast<int>(Cells[c].Face_Areas.size()))
                h_areas[c * kMaxFaces + f] = Cells[c].Face_Areas[f];
        }

        if (c < static_cast<int>(Primitive_Cells.size()))
        {
            const auto &p = Primitive_Cells[c];
            if (p.size() > 0)
                h_prim[c * kPrimStride + 0] = p[0];
            if (p.size() > 1)
                h_prim[c * kPrimStride + 1] = p[1];
            if (p.size() > 2)
                h_prim[c * kPrimStride + 2] = p[2];
            if (p.size() > 3)
                h_prim[c * kPrimStride + 3] = p[3];
            if (p.size() > 4)
                h_prim[c * kPrimStride + 4] = p[4];
            if (p.size() > 5)
                h_prim[c * kPrimStride + 5] = p[5];
            if (p.size() > 8)
                h_prim[c * kPrimStride + 8] = p[8];
            else if (p.size() > 3)
            {
                const double T = p[3];
                const double term1 = T_S_MU_CONST / T_REF_CONST;
                const double term2 = T + term1;
                h_prim[c * kPrimStride + 8] = (term2 > 1e-14)
                                                  ? pow(T, 1.5) * ((1.0 + term1) / term2)
                                                  : 0.0;
            }
        }
    }

    cudaMemcpy(g_visc.d_num_faces, h_num_faces.data(), n_cells * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_visc.d_neighbours, h_neighbours.data(), n_cells * kMaxFaces * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_visc.d_face_normals, h_normals.data(), n_cells * kMaxFaces * 2 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_visc.d_face_areas, h_areas.data(), n_cells * kMaxFaces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_visc.d_inv_area, h_inv_area.data(), n_cells * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_visc.d_prim, h_prim.data(), n_cells * kPrimStride * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_visc.d_centers, h_centers.data(), n_cells * 2 * sizeof(double), cudaMemcpyHostToDevice);

    const int threads = 256;
    const int blocks = (No_Physical_Cells + threads - 1) / threads;
    constexpr double kAlpha = 0.75;

    evaluate_viscous_flux_kernel<<<blocks, threads>>>(
        g_visc.d_num_faces, g_visc.d_neighbours, g_visc.d_face_normals, g_visc.d_face_areas,
        g_visc.d_inv_area, g_visc.d_prim, g_visc.d_centers, g_visc.d_viscous,
        n_cells, No_Physical_Cells, Inv_Re, K1, kAlpha);

    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Viscous flux CUDA kernel failed: %s\n", cudaGetErrorString(err));
        return false;
    }

    std::vector<double> h_visc(n_cells * kNv, 0.0);
    cudaMemcpy(h_visc.data(), g_visc.d_viscous, n_cells * kNv * sizeof(double), cudaMemcpyDeviceToHost);

    const int nv = std::min(kNv, NUM_FLUX_COMPONENTS);
    for (int c = 0; c < No_Physical_Cells; ++c)
    {
        if (c >= static_cast<int>(Cells_Viscous_Flux.size()))
            continue;
        for (int v = 0; v < nv; ++v)
            Cells_Viscous_Flux[c][v] = h_visc[c * kNv + v];
    }

    return true;
}
