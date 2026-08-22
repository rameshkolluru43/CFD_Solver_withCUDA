#include <cuda_runtime.h>
#include <math.h>

namespace
{
constexpr int kNumEq = 4;
constexpr int kPrimStride = 6;
constexpr double kGamma = 1.4;
constexpr double kGamma1 = kGamma / (kGamma - 1.0);

__device__ double device_sign(const double value)
{
    return (value > 0.0) - (value < 0.0);
}

__device__ double max3(const double a, const double b, const double c)
{
    return fmax(a, fmax(b, c));
}

__device__ double min3(const double a, const double b, const double c)
{
    return fmin(a, fmin(b, c));
}

__device__ void entropy_fix(double &alpha, const double lambda_max)
{
    const double delta = lambda_max;
    if (fabs(alpha) < delta && delta > 0.0)
        alpha = (alpha * alpha + delta * delta) / (2.0 * delta);
    else
        alpha = lambda_max;
}

__device__ double movers_alpha(const double d_u, const double d_f, const double lambda_max, const double lambda_min)
{
    const double eps = 1e-10;
    if (fabs(d_f) < eps && fabs(d_u) > eps)
        return 0.0;
    if (fabs(d_f) < eps && fabs(d_u) < eps)
        return lambda_min;
    if (fabs(d_f) > eps && fabs(d_u) > eps)
    {
        double alpha = fabs(d_f / d_u);
        if (alpha >= lambda_max)
            alpha = device_sign(alpha) * lambda_max;
        else if (alpha <= lambda_min)
            alpha = device_sign(alpha) * lambda_min;
        return alpha;
    }
    return lambda_min;
}

__device__ double movers_nwsc_alpha(const double d_u, const double d_f)
{
    const double eps = 1e-8;
    if (fabs(d_f) <= eps && fabs(d_u) >= eps)
        return 0.0;
    return device_sign(d_u) * fabs(d_f);
}

__device__ double ricca_alpha(const double d_u, const double d_f, const double vn_l, const double vn_r,
                              double d_p, const double rho_i, const double p_i)
{
    const double eps = 1e-8;
    if (d_p < eps)
        d_p = 0.0;

    if (fabs(d_f) < eps && fabs(d_u) < eps)
        return 0.5 * (fabs(vn_l) + fabs(vn_r));

    const double acoustic = (rho_i > eps && p_i > eps) ? sqrt(kGamma * p_i / rho_i) : 0.0;
    return fmax(fabs(vn_l), fabs(vn_r)) + device_sign(d_p) * acoustic;
}

__global__ void inviscid_flux_1o_kernel(const int *leaf_cells,
                                        int num_leaf_cells,
                                        const double *u_cells,
                                        const double *prim_cells,
                                        const int *face_offsets,
                                        const int *neighbours,
                                        const double *face_normals,
                                        const double *face_areas,
                                        const unsigned char *boundary_faces,
                                        const double *cell_areas,
                                        int dissipation_type,
                                        int is_movers_1,
                                        int enable_entropy_fix,
                                        double cfl,
                                        double *net_flux,
                                        double *cell_dt)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_leaf_cells)
        return;

    const int cell = leaf_cells[tid];
    const int f0 = face_offsets[cell];
    const int f1 = face_offsets[cell + 1];

    double local_net[kNumEq] = {0.0, 0.0, 0.0, 0.0};
    const double rho_l = prim_cells[cell * kPrimStride + 0];
    const double u_l = prim_cells[cell * kPrimStride + 1];
    const double v_l = prim_cells[cell * kPrimStride + 2];
    const double p_l_base = prim_cells[cell * kPrimStride + 4];
    const double a_l = prim_cells[cell * kPrimStride + 5];
    double dt_denom = 0.0;

    for (int face_pos = f0; face_pos < f1; ++face_pos)
    {
        const int neigh = neighbours[face_pos];
        if (neigh < 0)
            continue;

        const double nx = face_normals[2 * face_pos + 0];
        const double ny = face_normals[2 * face_pos + 1];
        const double dl = face_areas[face_pos];

        const double rho_r = prim_cells[neigh * kPrimStride + 0];
        const double u_r = prim_cells[neigh * kPrimStride + 1];
        const double v_r = prim_cells[neigh * kPrimStride + 2];
        double p_l = p_l_base;
        double p_r = prim_cells[neigh * kPrimStride + 4];
        const double a_r = prim_cells[neigh * kPrimStride + 5];

        const double vn_l = u_l * nx + v_l * ny;
        const double vn_r = u_r * nx + v_r * ny;
        const double vmag_l = 0.5 * (u_l * u_l + v_l * v_l);
        const double vmag_r = 0.5 * (u_r * u_r + v_r * v_r);
        dt_denom += (fabs(vn_l) + a_l) * dl;

        if (boundary_faces[face_pos])
            p_r = p_l;

        double flux_l[kNumEq];
        double flux_r[kNumEq];
        const double vn_l_dl = vn_l * dl;
        const double vn_r_dl = vn_r * dl;
        flux_l[0] = rho_l * vn_l_dl;
        flux_l[1] = rho_l * u_l * vn_l_dl + p_l * nx * dl;
        flux_l[2] = rho_l * v_l * vn_l_dl + p_l * ny * dl;
        flux_l[3] = (kGamma1 * p_l + rho_l * vmag_l) * vn_l_dl;

        flux_r[0] = rho_r * vn_r_dl;
        flux_r[1] = rho_r * u_r * vn_r_dl + p_r * nx * dl;
        flux_r[2] = rho_r * v_r * vn_r_dl + p_r * ny * dl;
        flux_r[3] = (kGamma1 * p_r + rho_r * vmag_r) * vn_r_dl;

        double d_u[kNumEq];
        double d_f[kNumEq];
        for (int eq = 0; eq < kNumEq; ++eq)
        {
            d_u[eq] = u_cells[neigh * kNumEq + eq] - u_cells[cell * kNumEq + eq];
            d_f[eq] = flux_r[eq] - flux_l[eq];
        }

        const double s0 = fabs(vn_l - a_l) * dl;
        const double s1 = fabs(vn_l + a_l) * dl;
        const double s2 = fabs(vn_l) * dl;
        const double s3 = fabs(vn_r - a_r) * dl;
        const double s4 = fabs(vn_r + a_r) * dl;
        const double s5 = fabs(vn_r) * dl;
        const double lambda_max = fmax(max3(s0, s1, s2), max3(s3, s4, s5));
        const double lambda_min = fmax(min3(s0, s1, s2), min3(s3, s4, s5));

        double diss[kNumEq] = {0.0, 0.0, 0.0, 0.0};
        if (dissipation_type == 1)
        {
            for (int eq = 0; eq < kNumEq; ++eq)
                diss[eq] = 0.5 * lambda_max * d_u[eq];
        }
        else if (dissipation_type == 2)
        {
            double alpha[kNumEq];
            alpha[3] = movers_alpha(d_u[3], d_f[3], lambda_max, lambda_min);
            alpha[0] = alpha[1] = alpha[2] = alpha[3];
            if (!is_movers_1)
            {
                alpha[0] = movers_alpha(d_u[0], d_f[0], lambda_max, lambda_min);
                alpha[1] = movers_alpha(d_u[1], d_f[1], lambda_max, lambda_min);
                alpha[2] = movers_alpha(d_u[2], d_f[2], lambda_max, lambda_min);
            }
            if (enable_entropy_fix)
                for (int eq = 0; eq < kNumEq; ++eq)
                    entropy_fix(alpha[eq], lambda_max);
            for (int eq = 0; eq < kNumEq; ++eq)
                diss[eq] = 0.5 * alpha[eq] * d_u[eq];
        }
        else if (dissipation_type == 4)
        {
            const double d_p = fabs(p_l - p_r);
            const double p_i = 0.5 * (p_l + p_r);
            const double rho_i = 0.5 * (rho_l + rho_r);
            for (int eq = 0; eq < kNumEq; ++eq)
            {
                const double alpha = ricca_alpha(d_u[eq], d_f[eq], vn_l, vn_r, d_p, rho_i, p_i);
                diss[eq] = 0.5 * alpha * dl * d_u[eq];
            }
        }
        else if (dissipation_type == 5)
        {
            const double p_i = 0.5 * (p_l + p_r);
            const double d_p_signed = p_r - p_l;
            const double alpha_p = 0.5 * (fabs(vn_r * dl) + fabs(vn_l * dl));
            const double sensor = (fabs(p_i) > 1e-12) ? 0.21 * fabs(d_p_signed / (2.0 * p_i)) : 0.0;
            for (int eq = 0; eq < kNumEq; ++eq)
            {
                const double alpha = movers_nwsc_alpha(d_u[eq], d_f[eq]);
                diss[eq] = 0.5 * (sensor * alpha + alpha_p * d_u[eq]);
            }
        }

        for (int eq = 0; eq < kNumEq; ++eq)
            local_net[eq] += 0.5 * (flux_l[eq] + flux_r[eq]) - diss[eq];
    }

    for (int eq = 0; eq < kNumEq; ++eq)
        net_flux[cell * kNumEq + eq] = local_net[eq];

    const double area = cell_areas[cell];
    cell_dt[cell] = (dt_denom > 0.0 && isfinite(dt_denom) && area > 0.0) ? cfl * area / dt_denom : 0.0;
}
} // namespace

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
                                                             double *d_cell_dt)
{
    const int block = 128;
    const int grid = (num_leaf_cells + block - 1) / block;
    inviscid_flux_1o_kernel<<<grid, block>>>(d_leaf_cells, num_leaf_cells, d_u_cells, d_prim_cells,
                                             d_face_offsets, d_neighbours, d_face_normals, d_face_areas,
                                             d_boundary_faces, d_cell_areas, dissipation_type, is_movers_1,
                                             enable_entropy_fix, cfl, d_net_flux, d_cell_dt);
    return cudaGetLastError();
}
