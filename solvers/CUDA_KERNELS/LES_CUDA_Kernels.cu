/**
 * @file LES_CUDA_Kernels.cu
 * @brief CUDA kernels for Large Eddy Simulation (LES) SGS models
 * @author CFD Solver Team
 * @date 2026
 *
 * GPU-accelerated implementations of:
 * - Smagorinsky model
 * - Dynamic Smagorinsky model
 * - WALE model
 * - Vreman model
 * - Sigma model
 */

#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <vector>

#include "../include/definitions.h"
#include "../include/Globals.h"
#include "../include/LES_Models.h"

namespace
{

//=============================================================================
// CUDA CONSTANTS AND HELPERS
//=============================================================================

constexpr int kMaxFaces = 8;
constexpr int kMaxNeighbors = 8;

__device__ bool cuda_ok_device(cudaError_t err)
{
    return err == cudaSuccess;
}

//=============================================================================
// DEVICE CONSTANTS
//=============================================================================

__constant__ double d_Cs_Smagorinsky;
__constant__ double d_Cw_WALE;
__constant__ double d_Cv_Vreman;
__constant__ double d_Csigma;
__constant__ double d_A_plus;
__constant__ int d_les_model;

//=============================================================================
// DEVICE HELPER FUNCTIONS
//=============================================================================

__device__ void d_calculate_velocity_gradient(
    int cell,
    const double* d_primitive,
    const double* d_centroids,
    const int* d_neighbors,
    int num_neighbors,
    int n_cells,
    double G[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            G[i][j] = 0.0;

    if (num_neighbors < 2)
        return;

    double u0 = d_primitive[cell * 9 + 1];
    double v0 = d_primitive[cell * 9 + 2];
    double x0 = d_centroids[cell * 2 + 0];
    double y0 = d_centroids[cell * 2 + 1];

    double sum_wxx = 0.0, sum_wyy = 0.0, sum_wxy = 0.0;
    double sum_wudx = 0.0, sum_wudy = 0.0;
    double sum_wvdx = 0.0, sum_wvdy = 0.0;

    for (int n = 0; n < num_neighbors; ++n)
    {
        int neigh = d_neighbors[cell * kMaxNeighbors + n];
        if (neigh < 0 || neigh >= n_cells)
            continue;

        double xn = d_centroids[neigh * 2 + 0];
        double yn = d_centroids[neigh * 2 + 1];
        double dx = xn - x0;
        double dy = yn - y0;
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq < 1e-28)
            continue;

        double w = 1.0 / dist_sq;
        double un = d_primitive[neigh * 9 + 1];
        double vn = d_primitive[neigh * 9 + 2];
        double du = un - u0;
        double dv = vn - v0;

        sum_wxx += w * dx * dx;
        sum_wyy += w * dy * dy;
        sum_wxy += w * dx * dy;
        sum_wudx += w * du * dx;
        sum_wudy += w * du * dy;
        sum_wvdx += w * dv * dx;
        sum_wvdy += w * dv * dy;
    }

    double det = sum_wxx * sum_wyy - sum_wxy * sum_wxy;
    if (fabs(det) > 1e-14)
    {
        double inv_det = 1.0 / det;
        G[0][0] = (sum_wyy * sum_wudx - sum_wxy * sum_wudy) * inv_det;
        G[0][1] = (sum_wxx * sum_wudy - sum_wxy * sum_wudx) * inv_det;
        G[1][0] = (sum_wyy * sum_wvdx - sum_wxy * sum_wvdy) * inv_det;
        G[1][1] = (sum_wxx * sum_wvdy - sum_wxy * sum_wvdx) * inv_det;
    }
}

__device__ void d_calculate_strain_rate(const double G[3][3], double S[3][3])
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            S[i][j] = 0.5 * (G[i][j] + G[j][i]);
        }
    }
}

__device__ double d_strain_rate_magnitude(const double S[3][3])
{
    double S_sq = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            S_sq += S[i][j] * S[i][j];
    return sqrt(2.0 * S_sq);
}

__device__ double d_van_driest_damping(double y_plus)
{
    return 1.0 - exp(-y_plus / d_A_plus);
}

//=============================================================================
// SMAGORINSKY MODEL KERNEL
//=============================================================================

__global__ void smagorinsky_kernel(
    const int* leaf_cells,
    int n_leaf,
    const double* d_primitive,
    const double* d_centroids,
    const int* d_neighbors,
    const int* d_num_neighbors,
    const double* d_filter_width,
    const double* d_y_plus,
    double* d_nu_sgs,
    double* d_mu_sgs,
    int n_cells,
    bool use_damping)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_leaf)
        return;

    int cell = leaf_cells[idx];
    if (cell < 0 || cell >= n_cells)
        return;

    int num_neighbors = d_num_neighbors[cell];
    
    double G[3][3], S[3][3];
    d_calculate_velocity_gradient(cell, d_primitive, d_centroids, d_neighbors, 
                                   num_neighbors, n_cells, G);
    d_calculate_strain_rate(G, S);
    
    double S_mag = d_strain_rate_magnitude(S);
    double delta = d_filter_width[cell];
    double Cs = d_Cs_Smagorinsky;

    if (use_damping)
    {
        double y_plus = d_y_plus[cell];
        if (y_plus < 100.0)
        {
            Cs *= d_van_driest_damping(y_plus);
        }
    }

    double nu_sgs = Cs * Cs * delta * delta * S_mag;
    d_nu_sgs[cell] = nu_sgs;

    double rho = d_primitive[cell * 9 + 0];
    if (rho < 1e-14)
        rho = 1e-14;
    d_mu_sgs[cell] = rho * nu_sgs;
}

//=============================================================================
// WALE MODEL KERNEL
//=============================================================================

__global__ void wale_kernel(
    const int* leaf_cells,
    int n_leaf,
    const double* d_primitive,
    const double* d_centroids,
    const int* d_neighbors,
    const int* d_num_neighbors,
    const double* d_filter_width,
    double* d_nu_sgs,
    double* d_mu_sgs,
    int n_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_leaf)
        return;

    int cell = leaf_cells[idx];
    if (cell < 0 || cell >= n_cells)
        return;

    int num_neighbors = d_num_neighbors[cell];

    double G[3][3], S[3][3];
    d_calculate_velocity_gradient(cell, d_primitive, d_centroids, d_neighbors,
                                   num_neighbors, n_cells, G);
    d_calculate_strain_rate(G, S);

    double g2[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                g2[i][j] += G[i][k] * G[k][j];

    double trace_g2 = g2[0][0] + g2[1][1] + g2[2][2];

    double Sd[3][3];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Sd[i][j] = 0.5 * (g2[i][j] + g2[j][i]);
            if (i == j)
                Sd[i][j] -= trace_g2 / 3.0;
        }
    }

    double Sd_sq = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            Sd_sq += Sd[i][j] * Sd[i][j];
    double Sd_mag = sqrt(Sd_sq);

    double S_mag = d_strain_rate_magnitude(S);

    double delta = d_filter_width[cell];
    double Cw = d_Cw_WALE;

    double Sd_52 = pow(Sd_mag, 2.5);
    double denom = pow(S_mag, 5.0) + Sd_52 * Sd_52;
    
    double nu_sgs = 0.0;
    if (denom > 1e-20)
    {
        nu_sgs = Cw * Cw * delta * delta * Sd_52 * Sd_52 / sqrt(denom);
    }

    d_nu_sgs[cell] = nu_sgs;

    double rho = d_primitive[cell * 9 + 0];
    if (rho < 1e-14)
        rho = 1e-14;
    d_mu_sgs[cell] = rho * nu_sgs;
}

//=============================================================================
// VREMAN MODEL KERNEL
//=============================================================================

__global__ void vreman_kernel(
    const int* leaf_cells,
    int n_leaf,
    const double* d_primitive,
    const double* d_centroids,
    const int* d_neighbors,
    const int* d_num_neighbors,
    const double* d_filter_width,
    double* d_nu_sgs,
    double* d_mu_sgs,
    int n_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_leaf)
        return;

    int cell = leaf_cells[idx];
    if (cell < 0 || cell >= n_cells)
        return;

    int num_neighbors = d_num_neighbors[cell];

    double G[3][3];
    d_calculate_velocity_gradient(cell, d_primitive, d_centroids, d_neighbors,
                                   num_neighbors, n_cells, G);

    double delta = d_filter_width[cell];

    double beta[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                beta[i][j] += delta * delta * G[k][i] * G[k][j];

    double B_beta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
                    beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
                    beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];
    B_beta = fmax(0.0, B_beta);

    double alpha_sq = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            alpha_sq += G[i][j] * G[i][j];

    double nu_sgs = 0.0;
    if (alpha_sq > 1e-14)
    {
        nu_sgs = d_Cv_Vreman * sqrt(B_beta / alpha_sq);
    }

    d_nu_sgs[cell] = nu_sgs;

    double rho = d_primitive[cell * 9 + 0];
    if (rho < 1e-14)
        rho = 1e-14;
    d_mu_sgs[cell] = rho * nu_sgs;
}

//=============================================================================
// SIGMA MODEL KERNEL
//=============================================================================

__global__ void sigma_kernel(
    const int* leaf_cells,
    int n_leaf,
    const double* d_primitive,
    const double* d_centroids,
    const int* d_neighbors,
    const int* d_num_neighbors,
    const double* d_filter_width,
    double* d_nu_sgs,
    double* d_mu_sgs,
    int n_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_leaf)
        return;

    int cell = leaf_cells[idx];
    if (cell < 0 || cell >= n_cells)
        return;

    int num_neighbors = d_num_neighbors[cell];

    double G[3][3];
    d_calculate_velocity_gradient(cell, d_primitive, d_centroids, d_neighbors,
                                   num_neighbors, n_cells, G);

    double GtG[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                GtG[i][j] += G[k][i] * G[k][j];

    double I1 = GtG[0][0] + GtG[1][1] + GtG[2][2];
    double I2 = GtG[0][0] * GtG[1][1] + GtG[0][0] * GtG[2][2] + GtG[1][1] * GtG[2][2] -
                GtG[0][1] * GtG[0][1] - GtG[0][2] * GtG[0][2] - GtG[1][2] * GtG[1][2];

    double p = I1 / 3.0;
    double r = I2 / 3.0 - p * p;

    double sigma_sq[3] = {p, p, p};

    if (fabs(r) > 1e-14)
    {
        double I3 = GtG[0][0] * (GtG[1][1] * GtG[2][2] - GtG[1][2] * GtG[1][2]) -
                    GtG[0][1] * (GtG[0][1] * GtG[2][2] - GtG[1][2] * GtG[0][2]) +
                    GtG[0][2] * (GtG[0][1] * GtG[1][2] - GtG[1][1] * GtG[0][2]);
        double q = (p * p * p - I1 * p * p / 2.0 + I2 * p / 2.0 - I3 / 2.0);

        double sqrtr = sqrt(-r);
        double cos_arg = q / (sqrtr * sqrtr * sqrtr);
        cos_arg = fmax(-1.0, fmin(1.0, cos_arg));
        double phi = acos(cos_arg) / 3.0;

        sigma_sq[0] = p + 2.0 * sqrtr * cos(phi);
        sigma_sq[1] = p + 2.0 * sqrtr * cos(phi + 2.0 * M_PI / 3.0);
        sigma_sq[2] = p + 2.0 * sqrtr * cos(phi + 4.0 * M_PI / 3.0);
    }

    for (int i = 0; i < 3; ++i)
        sigma_sq[i] = fmax(0.0, sigma_sq[i]);

    if (sigma_sq[0] < sigma_sq[1])
    {
        double tmp = sigma_sq[0];
        sigma_sq[0] = sigma_sq[1];
        sigma_sq[1] = tmp;
    }
    if (sigma_sq[1] < sigma_sq[2])
    {
        double tmp = sigma_sq[1];
        sigma_sq[1] = sigma_sq[2];
        sigma_sq[2] = tmp;
    }
    if (sigma_sq[0] < sigma_sq[1])
    {
        double tmp = sigma_sq[0];
        sigma_sq[0] = sigma_sq[1];
        sigma_sq[1] = tmp;
    }

    double sigma1 = sqrt(sigma_sq[0]);
    double sigma2 = sqrt(sigma_sq[1]);
    double sigma3 = sqrt(sigma_sq[2]);

    double delta = d_filter_width[cell];
    double nu_sgs = 0.0;

    if (sigma1 > 1e-14)
    {
        double D_sigma = sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3) / (sigma1 * sigma1);
        nu_sgs = d_Csigma * d_Csigma * delta * delta * D_sigma;
        nu_sgs = fmax(0.0, nu_sgs);
    }

    d_nu_sgs[cell] = nu_sgs;

    double rho = d_primitive[cell * 9 + 0];
    if (rho < 1e-14)
        rho = 1e-14;
    d_mu_sgs[cell] = rho * nu_sgs;
}

//=============================================================================
// EFFECTIVE VISCOSITY UPDATE KERNEL
//=============================================================================

__global__ void update_effective_viscosity_kernel(
    const int* leaf_cells,
    int n_leaf,
    const double* d_primitive,
    const double* d_mu_sgs,
    double* d_mu_effective,
    int n_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_leaf)
        return;

    int cell = leaf_cells[idx];
    if (cell < 0 || cell >= n_cells)
        return;

    double mu_lam = d_primitive[cell * 9 + 8];
    if (mu_lam < 1e-14)
        mu_lam = 1e-14;

    d_mu_effective[cell] = mu_lam + d_mu_sgs[cell];
}

//=============================================================================
// GPU BUFFER STRUCTURE
//=============================================================================

struct LES_GPU_Buffers
{
    int* d_leaf_cells;
    double* d_primitive;
    double* d_centroids;
    int* d_neighbors;
    int* d_num_neighbors;
    double* d_filter_width;
    double* d_y_plus;
    double* d_nu_sgs;
    double* d_mu_sgs;
    double* d_mu_effective;
    int n_cells;
    int n_leaf;
    bool initialized;

    LES_GPU_Buffers() : d_leaf_cells(nullptr), d_primitive(nullptr),
                        d_centroids(nullptr), d_neighbors(nullptr),
                        d_num_neighbors(nullptr), d_filter_width(nullptr),
                        d_y_plus(nullptr), d_nu_sgs(nullptr),
                        d_mu_sgs(nullptr), d_mu_effective(nullptr),
                        n_cells(0), n_leaf(0), initialized(false) {}
};

LES_GPU_Buffers g_les_buf;

bool cuda_ok(cudaError_t err, const char* msg)
{
    if (err != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error (%s): %s\n", msg, cudaGetErrorString(err));
        return false;
    }
    return true;
}

} // anonymous namespace

//=============================================================================
// PUBLIC GPU INTERFACE FUNCTIONS
//=============================================================================

bool LES_GPU_Available()
{
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    return (err == cudaSuccess && device_count > 0 && Is_LES && Use_GPU_LES);
}

void Initialize_LES_GPU()
{
    if (!LES_GPU_Available())
        return;

    int n_cells = Total_No_Cells;
    g_les_buf.n_cells = n_cells;

    cudaMalloc(&g_les_buf.d_primitive, n_cells * 9 * sizeof(double));
    cudaMalloc(&g_les_buf.d_centroids, n_cells * 2 * sizeof(double));
    cudaMalloc(&g_les_buf.d_neighbors, n_cells * kMaxNeighbors * sizeof(int));
    cudaMalloc(&g_les_buf.d_num_neighbors, n_cells * sizeof(int));
    cudaMalloc(&g_les_buf.d_filter_width, n_cells * sizeof(double));
    cudaMalloc(&g_les_buf.d_y_plus, n_cells * sizeof(double));
    cudaMalloc(&g_les_buf.d_nu_sgs, n_cells * sizeof(double));
    cudaMalloc(&g_les_buf.d_mu_sgs, n_cells * sizeof(double));
    cudaMalloc(&g_les_buf.d_mu_effective, n_cells * sizeof(double));

    std::vector<double> h_centroids(n_cells * 2);
    std::vector<int> h_neighbors(n_cells * kMaxNeighbors, -1);
    std::vector<int> h_num_neighbors(n_cells);
    std::vector<double> h_filter_width(n_cells);

    for (int i = 0; i < n_cells; ++i)
    {
        h_centroids[i * 2 + 0] = Cells[i].Cell_Center[0];
        h_centroids[i * 2 + 1] = Cells[i].Cell_Center[1];

        int num_n = static_cast<int>(Cells[i].Neighbours.size());
        h_num_neighbors[i] = std::min(num_n, kMaxNeighbors);

        for (int n = 0; n < h_num_neighbors[i]; ++n)
        {
            h_neighbors[i * kMaxNeighbors + n] = Cells[i].Neighbours[n];
        }

        h_filter_width[i] = les_vars[i].delta;
    }

    cudaMemcpy(g_les_buf.d_centroids, h_centroids.data(), n_cells * 2 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_les_buf.d_neighbors, h_neighbors.data(), n_cells * kMaxNeighbors * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_les_buf.d_num_neighbors, h_num_neighbors.data(), n_cells * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_les_buf.d_filter_width, h_filter_width.data(), n_cells * sizeof(double), cudaMemcpyHostToDevice);

    double h_Cs = Cs_Smagorinsky;
    double h_Cw = Cw_WALE;
    double h_Cv = Cv_Vreman;
    double h_Csig = Csigma;
    double h_Aplus = LES_Constants::A_plus;
    int h_model = static_cast<int>(current_les_model);

    cudaMemcpyToSymbol(d_Cs_Smagorinsky, &h_Cs, sizeof(double));
    cudaMemcpyToSymbol(d_Cw_WALE, &h_Cw, sizeof(double));
    cudaMemcpyToSymbol(d_Cv_Vreman, &h_Cv, sizeof(double));
    cudaMemcpyToSymbol(d_Csigma, &h_Csig, sizeof(double));
    cudaMemcpyToSymbol(d_A_plus, &h_Aplus, sizeof(double));
    cudaMemcpyToSymbol(d_les_model, &h_model, sizeof(int));

    g_les_buf.initialized = true;
    std::cout << "LES GPU buffers initialized for " << n_cells << " cells" << std::endl;
}

void Finalize_LES_GPU()
{
    if (!g_les_buf.initialized)
        return;

    cudaFree(g_les_buf.d_leaf_cells);
    cudaFree(g_les_buf.d_primitive);
    cudaFree(g_les_buf.d_centroids);
    cudaFree(g_les_buf.d_neighbors);
    cudaFree(g_les_buf.d_num_neighbors);
    cudaFree(g_les_buf.d_filter_width);
    cudaFree(g_les_buf.d_y_plus);
    cudaFree(g_les_buf.d_nu_sgs);
    cudaFree(g_les_buf.d_mu_sgs);
    cudaFree(g_les_buf.d_mu_effective);

    g_les_buf = LES_GPU_Buffers();
}

void Copy_LES_Data_To_Device()
{
    if (!g_les_buf.initialized)
        return;

    int n_cells = g_les_buf.n_cells;

    std::vector<double> h_primitive(n_cells * 9);
    std::vector<double> h_y_plus(n_cells);

    for (int i = 0; i < n_cells; ++i)
    {
        for (int v = 0; v < 9; ++v)
        {
            h_primitive[i * 9 + v] = Primitive_Cells[i][v];
        }
        h_y_plus[i] = les_vars[i].y_plus;
    }

    cudaMemcpy(g_les_buf.d_primitive, h_primitive.data(), n_cells * 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(g_les_buf.d_y_plus, h_y_plus.data(), n_cells * sizeof(double), cudaMemcpyHostToDevice);
}

void Copy_LES_Data_From_Device()
{
    if (!g_les_buf.initialized)
        return;

    int n_cells = g_les_buf.n_cells;

    std::vector<double> h_nu_sgs(n_cells);
    std::vector<double> h_mu_sgs(n_cells);
    std::vector<double> h_mu_effective(n_cells);

    cudaMemcpy(h_nu_sgs.data(), g_les_buf.d_nu_sgs, n_cells * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mu_sgs.data(), g_les_buf.d_mu_sgs, n_cells * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mu_effective.data(), g_les_buf.d_mu_effective, n_cells * sizeof(double), cudaMemcpyDeviceToHost);

    for (int i = 0; i < n_cells; ++i)
    {
        les_vars[i].nu_sgs = h_nu_sgs[i];
        les_vars[i].mu_sgs = h_mu_sgs[i];
        les_vars[i].mu_effective = h_mu_effective[i];
    }
}

bool Calculate_SGS_Viscosity_GPU()
{
    if (!g_les_buf.initialized || !LES_GPU_Available())
        return false;

    Copy_LES_Data_To_Device();

    V_I leafCells;
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
            leafCells.push_back(i);
    }
    int n_leaf = static_cast<int>(leafCells.size());

    if (g_les_buf.d_leaf_cells == nullptr || g_les_buf.n_leaf != n_leaf)
    {
        cudaFree(g_les_buf.d_leaf_cells);
        cudaMalloc(&g_les_buf.d_leaf_cells, n_leaf * sizeof(int));
        g_les_buf.n_leaf = n_leaf;
    }
    cudaMemcpy(g_les_buf.d_leaf_cells, leafCells.data(), n_leaf * sizeof(int), cudaMemcpyHostToDevice);

    const int threads = 256;
    const int blocks = (n_leaf + threads - 1) / threads;

    switch (current_les_model)
    {
    case LES_SGS_Model::SMAGORINSKY:
        smagorinsky_kernel<<<blocks, threads>>>(
            g_les_buf.d_leaf_cells, n_leaf,
            g_les_buf.d_primitive, g_les_buf.d_centroids,
            g_les_buf.d_neighbors, g_les_buf.d_num_neighbors,
            g_les_buf.d_filter_width, g_les_buf.d_y_plus,
            g_les_buf.d_nu_sgs, g_les_buf.d_mu_sgs,
            g_les_buf.n_cells, Use_Wall_Damping);
        break;

    case LES_SGS_Model::WALE:
        wale_kernel<<<blocks, threads>>>(
            g_les_buf.d_leaf_cells, n_leaf,
            g_les_buf.d_primitive, g_les_buf.d_centroids,
            g_les_buf.d_neighbors, g_les_buf.d_num_neighbors,
            g_les_buf.d_filter_width,
            g_les_buf.d_nu_sgs, g_les_buf.d_mu_sgs,
            g_les_buf.n_cells);
        break;

    case LES_SGS_Model::VREMAN:
        vreman_kernel<<<blocks, threads>>>(
            g_les_buf.d_leaf_cells, n_leaf,
            g_les_buf.d_primitive, g_les_buf.d_centroids,
            g_les_buf.d_neighbors, g_les_buf.d_num_neighbors,
            g_les_buf.d_filter_width,
            g_les_buf.d_nu_sgs, g_les_buf.d_mu_sgs,
            g_les_buf.n_cells);
        break;

    case LES_SGS_Model::SIGMA:
        sigma_kernel<<<blocks, threads>>>(
            g_les_buf.d_leaf_cells, n_leaf,
            g_les_buf.d_primitive, g_les_buf.d_centroids,
            g_les_buf.d_neighbors, g_les_buf.d_num_neighbors,
            g_les_buf.d_filter_width,
            g_les_buf.d_nu_sgs, g_les_buf.d_mu_sgs,
            g_les_buf.n_cells);
        break;

    default:
        return false;
    }

    if (!cuda_ok(cudaGetLastError(), "LES SGS kernel launch"))
        return false;

    update_effective_viscosity_kernel<<<blocks, threads>>>(
        g_les_buf.d_leaf_cells, n_leaf,
        g_les_buf.d_primitive, g_les_buf.d_mu_sgs,
        g_les_buf.d_mu_effective, g_les_buf.n_cells);

    if (!cuda_ok(cudaDeviceSynchronize(), "LES kernel sync"))
        return false;

    Copy_LES_Data_From_Device();

    return true;
}

bool Calculate_Smagorinsky_GPU()
{
    if (current_les_model != LES_SGS_Model::SMAGORINSKY)
        return false;
    return Calculate_SGS_Viscosity_GPU();
}

bool Calculate_WALE_GPU()
{
    if (current_les_model != LES_SGS_Model::WALE)
        return false;
    return Calculate_SGS_Viscosity_GPU();
}

bool Calculate_Vreman_GPU()
{
    if (current_les_model != LES_SGS_Model::VREMAN)
        return false;
    return Calculate_SGS_Viscosity_GPU();
}

bool Calculate_Dynamic_Smagorinsky_GPU()
{
    return false;
}
