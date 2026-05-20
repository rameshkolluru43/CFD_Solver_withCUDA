/**
 * @file MOVERS_RICCA_Flux_Cuda_Kernels.cu
 * @brief Device net-flux kernels for MOVERS, RICCA, and MOVERS_NWSC (1st order).
 */

#include <cuda_runtime.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>

#include "../include/definitions.h"
#include "../include/Globals.h"
#include "../include/Grid.h"
#include "MOVERS_RICCA_Flux_Cuda.h"

namespace
{
constexpr int kMaxFaces = 8;
constexpr int kNv = 4;

__device__ double d_sign(double v)
{
    return (v >= 0.0) ? 1.0 : -1.0;
}

__device__ void d_max3(double a, double b, double c, double &out)
{
    out = fmax(a, fmax(b, c));
}

__device__ void d_min3(double a, double b, double c, double &out)
{
    out = fmin(a, fmin(b, c));
}

__device__ void condition_for_movers(double d_U, double d_F, double L_Max, double L_Min, double &Alpha)
{
    const double epsilon = 1e-10;
    if (fabs(d_F) < epsilon && fabs(d_U) > epsilon)
        Alpha = 0.0;
    else if (fabs(d_F) < epsilon && fabs(d_U) < epsilon)
        Alpha = L_Min;
    else if (fabs(d_F) > epsilon && fabs(d_U) > epsilon)
    {
        Alpha = fabs(d_F / d_U);
        if (Alpha >= L_Max)
            Alpha = d_sign(Alpha) * L_Max;
        else if (Alpha <= L_Min)
            Alpha = d_sign(Alpha) * L_Min;
    }
    else
        Alpha = L_Min;
}

__device__ void entropy_fix(double &Alpha, double L_Max)
{
    const double k = 1.0;
    const double delta = k * L_Max;
    if (fabs(Alpha) < delta)
        Alpha = (Alpha * Alpha + delta * delta) / (2.0 * delta);
    else
        Alpha = L_Max;
}

__device__ void condition_for_ricca(double d_U, double d_F, double Vn_L, double Vn_R, double dP,
                                    double Rho_I, double P_I, double &Alpha, double gamma)
{
    const double epsilon = 1e-8;
    Alpha = 0.0;
    if (dP < epsilon)
        dP = 0.0;
    const double Vmax = fmax(fabs(Vn_L), fabs(Vn_R));
    if (fabs(d_F) < epsilon && fabs(d_U) < epsilon)
        Alpha = 0.5 * (fabs(Vn_L) + fabs(Vn_R));
    else
        Alpha = Vmax + d_sign(dP) * sqrt(gamma * P_I / Rho_I);
}

__device__ void condition_for_movers_nwsc(double d_U, double d_F, double &Alpha)
{
    const double epsilon = 1e-8;
    if (fabs(d_F) <= epsilon && fabs(d_U) >= epsilon)
        Alpha = 0.0;
    else
        Alpha = d_sign(d_U) * fabs(d_F);
}

__device__ void face_average_convective(
    double rho_L, double u_L, double v_L, double P_L,
    double rho_R, double u_R, double v_R, double P_R,
    double nx, double ny, double dl, bool is_wall, double gamma1,
    double out_avg[kNv])
{
    if (is_wall)
        P_R = P_L;

    const double Vdotn_L = u_L * nx + v_L * ny;
    const double Vdotn_R = u_R * nx + v_R * ny;
    const double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    const double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    const double Vdotn_L_dl = Vdotn_L * dl;
    const double Vdotn_R_dl = Vdotn_R * dl;

    const double Flux_L[kNv] = {
        rho_L * Vdotn_L_dl,
        rho_L * u_L * Vdotn_L_dl + P_L * nx * dl,
        rho_L * v_L * Vdotn_L_dl + P_L * ny * dl,
        (gamma1 * P_L + rho_L * Vmag_L) * Vdotn_L_dl};

    const double Flux_R[kNv] = {
        rho_R * Vdotn_R_dl,
        rho_R * u_R * Vdotn_R_dl + P_R * nx * dl,
        rho_R * v_R * Vdotn_R_dl + P_R * ny * dl,
        (gamma1 * P_R + rho_R * Vmag_R) * Vdotn_R_dl};

    for (int i = 0; i < kNv; ++i)
        out_avg[i] = 0.5 * (Flux_L[i] + Flux_R[i]);
}

__device__ void dissipative_movers(
    int cell, int neigh,
    const double *d_prim, const double *d_U,
    double nx, double ny, double dl,
    int dissipation_type, bool is_movers_1, bool entropy_fix_on, double gamma,
    double out_diss[kNv])
{
    const int Lp = cell * 6;
    const int Rp = neigh * 6;
    const double rho_L = d_prim[Lp + 0];
    const double u_L = d_prim[Lp + 1];
    const double v_L = d_prim[Lp + 2];
    const double P_L = d_prim[Lp + 4];
    const double C_L = d_prim[Lp + 5];

    const double rho_R = d_prim[Rp + 0];
    const double u_R = d_prim[Rp + 1];
    const double v_R = d_prim[Rp + 2];
    const double P_R = d_prim[Rp + 4];
    const double C_R = d_prim[Rp + 5];

    const int Lc = cell * kNv;
    const int Rc = neigh * kNv;
    const double d_Uv[kNv] = {
        d_U[Rc + 0] - d_U[Lc + 0],
        d_U[Rc + 1] - d_U[Lc + 1],
        d_U[Rc + 2] - d_U[Lc + 2],
        d_U[Rc + 3] - d_U[Lc + 3]};

    const double Vdotn_L = u_L * nx + v_L * ny;
    const double Vdotn_R = u_R * nx + v_R * ny;
    const double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    const double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    const double d_F[kNv] = {
        (rho_R * Vdotn_R - rho_L * Vdotn_L) * dl,
        (rho_R * u_R * Vdotn_R + P_R * nx - rho_L * u_L * Vdotn_L + P_L * nx) * dl,
        (rho_R * v_R * Vdotn_R + P_R * ny - rho_L * v_L * Vdotn_L + P_L * ny) * dl,
        ((((P_R / (gamma - 1.0)) + rho_R * Vmag_R) + P_R) * Vdotn_R -
         (((P_L / (gamma - 1.0)) + rho_L * Vmag_L) + P_L) * Vdotn_L) * dl};

    double S[6];
    S[0] = fabs(Vdotn_L - C_L) * dl;
    S[1] = fabs(Vdotn_L + C_L) * dl;
    S[2] = fabs(Vdotn_L) * dl;
    S[3] = fabs(Vdotn_R - C_R) * dl;
    S[4] = fabs(Vdotn_R + C_R) * dl;
    S[5] = fabs(Vdotn_R) * dl;

    double Max1, Max2, Min1, Min2, Lambda_Max, Lambda_Min;
    d_max3(S[0], S[1], S[2], Max1);
    d_max3(S[3], S[4], S[5], Max2);
    Lambda_Max = fmax(Max1, Max2);
    d_min3(S[0], S[1], S[2], Min1);
    d_min3(S[3], S[4], S[5], Min2);
    Lambda_Min = fmax(Min1, Min2);

    double alpha[kNv];
    condition_for_movers(d_Uv[3], d_F[3], Lambda_Max, Lambda_Min, alpha[3]);
    alpha[0] = alpha[3];
    alpha[1] = alpha[3];
    alpha[2] = alpha[3];

    if (!is_movers_1)
    {
        condition_for_movers(d_Uv[0], d_F[0], Lambda_Max, Lambda_Min, alpha[0]);
        condition_for_movers(d_Uv[1], d_F[1], Lambda_Max, Lambda_Min, alpha[1]);
        condition_for_movers(d_Uv[2], d_F[2], Lambda_Max, Lambda_Min, alpha[2]);
    }

    if (entropy_fix_on)
    {
        entropy_fix(alpha[0], Lambda_Max);
        entropy_fix(alpha[1], Lambda_Max);
        entropy_fix(alpha[2], Lambda_Max);
        entropy_fix(alpha[3], Lambda_Max);
    }

    for (int i = 0; i < kNv; ++i)
        out_diss[i] = 0.5 * alpha[i] * d_Uv[i];
}

__device__ void dissipative_ricca(
    int cell, int neigh,
    const double *d_prim,
    double nx, double ny, double dl, double gamma,
    double out_diss[kNv])
{
    const int Lp = cell * 6;
    const int Rp = neigh * 6;
    const double rho_L = d_prim[Lp + 0];
    const double u_L = d_prim[Lp + 1];
    const double v_L = d_prim[Lp + 2];
    const double P_L = d_prim[Lp + 4];

    const double rho_R = d_prim[Rp + 0];
    const double u_R = d_prim[Rp + 1];
    const double v_R = d_prim[Rp + 2];
    const double P_R = d_prim[Rp + 4];

    const double Vdotn_L = u_L * nx + v_L * ny;
    const double Vdotn_R = u_R * nx + v_R * ny;
    const double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    const double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    double d_Uv[kNv] = {
        rho_R - rho_L,
        rho_R * u_R - rho_L * u_L,
        rho_R * v_R - rho_L * v_L,
        ((P_R / (gamma - 1.0)) + 0.5 * rho_R * (u_R * u_R + v_R * v_R)) -
            ((P_L / (gamma - 1.0)) + 0.5 * rho_L * (u_L * u_L + v_L * v_L))};

    const double d_F[kNv] = {
        (rho_R * Vdotn_R - rho_L * Vdotn_L) * dl,
        (rho_R * u_R * Vdotn_R + P_R * nx - rho_L * u_L * Vdotn_L + P_L * nx) * dl,
        (rho_R * v_R * Vdotn_R + P_R * ny - rho_L * v_L * Vdotn_L + P_L * ny) * dl,
        ((((P_R / (gamma - 1.0)) + rho_R * Vmag_R) + P_R) * Vdotn_R -
         (((P_L / (gamma - 1.0)) + rho_L * Vmag_L) + P_L) * Vdotn_L) * dl};

    double dP = fabs(P_L - P_R);
    const double P_I = 0.5 * (P_L + P_R);
    const double Rho_I = 0.5 * (rho_L + rho_R);

    for (int i = 0; i < kNv; ++i)
    {
        double alpha = 0.0;
        condition_for_ricca(d_Uv[i], d_F[i], Vdotn_L, Vdotn_R, dP, Rho_I, P_I, alpha, gamma);
        out_diss[i] = 0.5 * alpha * dl * d_Uv[i];
    }
}

__device__ void dissipative_movers_nwsc(
    int cell, int neigh,
    const double *d_prim, const double *d_U,
    double nx, double ny, double dl,
    double out_diss[kNv])
{
    const int Lp = cell * 6;
    const int Rp = neigh * 6;
    const double rho_L = d_prim[Lp + 0];
    const double u_L = d_prim[Lp + 1];
    const double v_L = d_prim[Lp + 2];
    const double P_L = d_prim[Lp + 4];

    const double rho_R = d_prim[Rp + 0];
    const double u_R = d_prim[Rp + 1];
    const double v_R = d_prim[Rp + 2];
    const double P_R = d_prim[Rp + 4];

    const int Lc = cell * kNv;
    const int Rc = neigh * kNv;
    const double d_Uv[kNv] = {
        d_U[Rc + 0] - d_U[Lc + 0],
        d_U[Rc + 1] - d_U[Lc + 1],
        d_U[Rc + 2] - d_U[Lc + 2],
        d_U[Rc + 3] - d_U[Lc + 3]};

    const double Vdotn_L = u_L * nx + v_L * ny;
    const double Vdotn_R = u_R * nx + v_R * ny;
    const double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    const double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    const double d_F[kNv] = {
        (rho_R * Vdotn_R - rho_L * Vdotn_L) * dl,
        (rho_R * u_R * Vdotn_R + P_R * nx - rho_L * u_L * Vdotn_L + P_L * nx) * dl,
        (rho_R * v_R * Vdotn_R + P_R * ny - rho_L * v_L * Vdotn_L + P_L * ny) * dl,
        ((((P_R / (gamma - 1.0)) + rho_R * Vmag_R) + P_R) * Vdotn_R -
         (((P_L / (gamma - 1.0)) + rho_L * Vmag_L) + P_L) * Vdotn_L) * dl};

    double alpha[kNv];
    for (int i = 0; i < kNv; ++i)
        condition_for_movers_nwsc(d_Uv[i], d_F[i], alpha[i]);

    const double d_P = P_R - P_L;
    const double P_I = 0.5 * (P_R + P_L);
    const double Alpha_P = 0.5 * (fabs(Vdotn_R * dl) + fabs(Vdotn_L * dl));
    const double beta = 0.21;
    const double Sensor = (fabs(P_I) > 1e-14) ? beta * fabs(d_P / (2.0 * P_I)) : 0.0;

    for (int i = 0; i < kNv; ++i)
        out_diss[i] = 0.5 * (Sensor * alpha[i] + Alpha_P * d_Uv[i]);
}

__global__ void net_flux_movers_ricca_kernel(
    const int *leaf_cells,
    int n_leaf,
    const int *d_num_faces,
    const int *d_neighbours,
    const double *d_face_normals,
    const double *d_face_areas,
    const int *d_wall_face,
    const double *d_prim,
    const double *d_U,
    double *d_net_flux,
    int n_cells,
    int dissipation_type,
    bool is_movers_1,
    bool entropy_fix_on,
    double gamma,
    double gamma1)
{
    const int li = blockIdx.x * blockDim.x + threadIdx.x;
    if (li >= n_leaf)
        return;

    const int cell = leaf_cells[li];
    if (cell < 0 || cell >= n_cells)
        return;

    for (int v = 0; v < kNv; ++v)
        d_net_flux[cell * kNv + v] = 0.0;

    const int nf = d_num_faces[cell];
    for (int face = 0; face < nf && face < kMaxFaces; ++face)
    {
        const int neigh = d_neighbours[cell * kMaxFaces + face];
        if (neigh < 0 || neigh >= n_cells)
            continue;

        const int norm_base = (cell * kMaxFaces + face) * 2;
        const double nx = d_face_normals[norm_base + 0];
        const double ny = d_face_normals[norm_base + 1];
        const double dl = d_face_areas[cell * kMaxFaces + face];
        if (dl <= 0.0)
            continue;

        const bool is_wall = d_wall_face[cell * kMaxFaces + face] != 0;

        const int Lp = cell * 6;
        const int Rp = neigh * 6;
        double rho_L = d_prim[Lp + 0], u_L = d_prim[Lp + 1], v_L = d_prim[Lp + 2], P_L = d_prim[Lp + 4];
        double rho_R = d_prim[Rp + 0], u_R = d_prim[Rp + 1], v_R = d_prim[Rp + 2], P_R = d_prim[Rp + 4];

        double avg[kNv], diss[kNv];
        face_average_convective(rho_L, u_L, v_L, P_L, rho_R, u_R, v_R, P_R, nx, ny, dl, is_wall, gamma1, avg);

        if (dissipation_type == 2)
            dissipative_movers(cell, neigh, d_prim, d_U, nx, ny, dl, dissipation_type, is_movers_1, entropy_fix_on, gamma, diss);
        else if (dissipation_type == 4)
            dissipative_ricca(cell, neigh, d_prim, nx, ny, dl, gamma, diss);
        else if (dissipation_type == 5)
            dissipative_movers_nwsc(cell, neigh, d_prim, d_U, nx, ny, dl, diss);
        else
            for (int v = 0; v < kNv; ++v)
                diss[v] = 0.0;

        for (int v = 0; v < kNv; ++v)
            d_net_flux[cell * kNv + v] += avg[v] - diss[v];
    }
}

struct CudaFluxBuffers
{
    int cell_cap = 0;
    int leaf_cap = 0;
    int *d_num_faces = nullptr;
    int *d_neighbours = nullptr;
    double *d_face_normals = nullptr;
    double *d_face_areas = nullptr;
    int *d_wall_face = nullptr;
    double *d_prim = nullptr;
    double *d_U = nullptr;
    double *d_net_flux = nullptr;
    int *d_leaf = nullptr;

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
        if (d_wall_face)
            cudaFree(d_wall_face);
        if (d_prim)
            cudaFree(d_prim);
        if (d_U)
            cudaFree(d_U);
        if (d_net_flux)
            cudaFree(d_net_flux);
        if (d_leaf)
            cudaFree(d_leaf);
        *this = CudaFluxBuffers{};
    }
};

CudaFluxBuffers g_buf;

bool ensure_capacity(int n_cells, int n_leaf)
{
    if (n_cells <= g_buf.cell_cap && n_leaf <= g_buf.leaf_cap)
        return true;

    g_buf.free_all();
    g_buf.cell_cap = n_cells;
    g_buf.leaf_cap = n_leaf;

    auto chk = [](cudaError_t e, const char *msg) {
        if (e != cudaSuccess)
        {
            fprintf(stderr, "MOVERS/RICCA CUDA: %s (%s)\n", msg, cudaGetErrorString(e));
            return false;
        }
        return true;
    };

    if (!chk(cudaMalloc(&g_buf.d_num_faces, n_cells * sizeof(int)), "malloc num_faces"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_neighbours, n_cells * kMaxFaces * sizeof(int)), "malloc neighbours"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_face_normals, n_cells * kMaxFaces * 2 * sizeof(double)), "malloc normals"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_face_areas, n_cells * kMaxFaces * sizeof(double)), "malloc areas"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_wall_face, n_cells * kMaxFaces * sizeof(int)), "malloc wall"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_prim, n_cells * 6 * sizeof(double)), "malloc prim"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_U, n_cells * kNv * sizeof(double)), "malloc U"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_net_flux, n_cells * kNv * sizeof(double)), "malloc net_flux"))
        return false;
    if (!chk(cudaMalloc(&g_buf.d_leaf, n_leaf * sizeof(int)), "malloc leaf"))
        return false;
    return true;
}

} // namespace

bool Evaluate_Cell_Net_Flux_CUDA_Movers_Ricca(bool second_order)
{
    if (second_order)
        return false;

    if (Dissipation_Type != 2 && Dissipation_Type != 4 && Dissipation_Type != 5)
        return false;

    const int n_cells = static_cast<int>(Cells.size());
    if (n_cells <= 0)
        return false;

    V_I leafCells;
    Build_Leaf_Cell_List(leafCells);
    const int n_leaf = static_cast<int>(leafCells.size());
    if (n_leaf <= 0)
        return false;

    if (!ensure_capacity(n_cells, n_leaf))
        return false;

    std::vector<int> h_num_faces(n_cells, 0);
    std::vector<int> h_neighbours(n_cells * kMaxFaces, -1);
    std::vector<double> h_normals(n_cells * kMaxFaces * 2, 0.0);
    std::vector<double> h_areas(n_cells * kMaxFaces, 0.0);
    std::vector<int> h_wall(n_cells * kMaxFaces, 0);
    std::vector<double> h_prim(n_cells * 6, 0.0);
    std::vector<double> h_U(n_cells * kNv, 0.0);

    for (int c = 0; c < n_cells; ++c)
    {
        const int nf = (Cells[c].numFaces > 0)
                           ? Cells[c].numFaces
                           : static_cast<int>(Cells[c].Face_Areas.size());
        h_num_faces[c] = std::min(nf, kMaxFaces);

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

            bool wall = false;
            if (c < static_cast<int>(Cells_Face_Boundary_Type.size()) &&
                f < static_cast<int>(Cells_Face_Boundary_Type[c].size()))
                wall = Cells_Face_Boundary_Type[c][f];
            h_wall[c * kMaxFaces + f] = wall ? 1 : 0;
        }

        if (c < static_cast<int>(Primitive_Cells.size()) && Primitive_Cells[c].size() >= 6)
        {
            h_prim[c * 6 + 0] = Primitive_Cells[c][0];
            h_prim[c * 6 + 1] = Primitive_Cells[c][1];
            h_prim[c * 6 + 2] = Primitive_Cells[c][2];
            h_prim[c * 6 + 4] = Primitive_Cells[c][4];
            h_prim[c * 6 + 5] = Primitive_Cells[c][5];
        }
        if (c < static_cast<int>(U_Cells.size()) && U_Cells[c].size() >= static_cast<size_t>(kNv))
        {
            for (int v = 0; v < kNv; ++v)
                h_U[c * kNv + v] = U_Cells[c][v];
        }
    }

    cudaMemcpy(g_buf.d_num_faces, h_num_faces.data(), n_cells * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_neighbours, h_neighbours.data(), n_cells * kMaxFaces * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_face_normals, h_normals.data(), n_cells * kMaxFaces * 2 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_face_areas, h_areas.data(), n_cells * kMaxFaces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_wall_face, h_wall.data(), n_cells * kMaxFaces * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_prim, h_prim.data(), n_cells * 6 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_U, h_U.data(), n_cells * kNv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_buf.d_leaf, leafCells.data(), n_leaf * sizeof(int), cudaMemcpyHostToDevice);

    const int threads = 256;
    const int blocks = (n_leaf + threads - 1) / threads;

    net_flux_movers_ricca_kernel<<<blocks, threads>>>(
        g_buf.d_leaf, n_leaf,
        g_buf.d_num_faces, g_buf.d_neighbours,
        g_buf.d_face_normals, g_buf.d_face_areas, g_buf.d_wall_face,
        g_buf.d_prim, g_buf.d_U, g_buf.d_net_flux,
        n_cells, Dissipation_Type,
        Is_MOVERS_1, Enable_Entropy_Fix,
        GAMMA_CONST, GAMMA1_CONST);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "MOVERS/RICCA kernel launch failed: %s\n", cudaGetErrorString(err));
        return false;
    }
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "MOVERS/RICCA kernel sync failed: %s\n", cudaGetErrorString(err));
        return false;
    }

    std::vector<double> h_net(n_cells * kNv, 0.0);
    cudaMemcpy(h_net.data(), g_buf.d_net_flux, n_cells * kNv * sizeof(double), cudaMemcpyDeviceToHost);

    const int nv = std::min(kNv, NUM_FLUX_COMPONENTS);
    for (int idx = 0; idx < n_leaf; ++idx)
    {
        const int c = leafCells[idx];
        if (c < 0 || c >= static_cast<int>(Cells_Net_Flux.size()))
            continue;
        for (int v = 0; v < nv; ++v)
            Cells_Net_Flux[c][v] = h_net[c * kNv + v];
        Evaluate_Time_Step(c);
    }

    return true;
}
