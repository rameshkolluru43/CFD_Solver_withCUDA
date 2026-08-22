/**
 * @file MOVERS_RICCA_Flux_Cuda_Kernels.cu
 * @brief Device net-flux kernels for MOVERS, RICCA, and MOVERS_NWSC (1st order).
 */

#include <cuda_runtime.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>

#include "../include/definitions.h"
#include "../include/Globals.h"
#include "../include/Grid.h"
#include "../include/Timestep.h"
#include "../include/MPI_Utils.h"
#include "MOVERS_RICCA_Flux_Cuda.h"
#include <cmath>

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
    /* Slip wall: ghost velocity from reflection (host BC); do not zero V·n. */
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

/* -------------------- WENO5 (host WENO2D-compatible) -------------------- */

__device__ double d_weno5_scalar(double a, double b, double c, double d, double e, int shift)
{
    if (!isfinite(a) || !isfinite(b) || !isfinite(c) || !isfinite(d) || !isfinite(e))
        return c;
    const double epsilon = 1e-6;
    double d0, d1, d2, v0, v1, v2;
    if (shift == 0)
    {
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        v2 = (2.0 * a - 7.0 * b + 11.0 * c) / 6.0;
        v1 = (-b + 5.0 * c + 2.0 * d) / 6.0;
        v0 = (2.0 * c + 5.0 * d - e) / 6.0;
    }
    else
    {
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        v2 = (-a + 5.0 * b + 2.0 * c) / 6.0;
        v1 = (2.0 * b + 5.0 * c - d) / 6.0;
        v0 = (11.0 * c - 7.0 * d + 2.0 * e) / 6.0;
    }
    const double t2 = a - 2.0 * b + c;
    const double t2b = a - 4.0 * b + 3.0 * c;
    const double b2 = (13.0 / 12.0) * t2 * t2 + 0.25 * t2b * t2b;
    const double t1 = b - 2.0 * c + d;
    const double t1b = b - d;
    const double b1 = (13.0 / 12.0) * t1 * t1 + 0.25 * t1b * t1b;
    const double t0 = c - 2.0 * d + e;
    const double t0b = 3.0 * c - 4.0 * d + e;
    const double b0 = (13.0 / 12.0) * t0 * t0 + 0.25 * t0b * t0b;
    const double s0 = epsilon + b0;
    const double s1 = epsilon + b1;
    const double s2 = epsilon + b2;
    const double a0 = d0 / (s0 * s0);
    const double a1 = d1 / (s1 * s1);
    const double a2 = d2 / (s2 * s2);
    const double sum = a0 + a1 + a2;
    double U;
    if (sum < epsilon)
        U = (v0 + v1 + v2) / 3.0;
    else
        U = (a0 * v0 + a1 * v1 + a2 * v2) / sum;
    return isfinite(U) ? U : c;
}

__device__ void d_slip_reflect_U(const double U_in[kNv], double nx, double ny, double gamma, double U_out[kNv])
{
    const double rho = U_in[0];
    if (!(rho > 1e-14) || !isfinite(rho))
    {
        for (int v = 0; v < kNv; ++v)
            U_out[v] = U_in[v];
        return;
    }
    const double u = U_in[1] / rho;
    const double vv = U_in[2] / rho;
    const double un = u * nx + vv * ny;
    const double ug = u - 2.0 * un * nx;
    const double vg = vv - 2.0 * un * ny;
    double p = (gamma - 1.0) * (U_in[3] - 0.5 * rho * (u * u + vv * vv));
    if (!(p > 0.0) || !isfinite(p))
        p = 1e-12;
    U_out[0] = rho;
    U_out[1] = rho * ug;
    U_out[2] = rho * vg;
    U_out[3] = p / (gamma - 1.0) + rho * 0.5 * (ug * ug + vg * vg);
}

__device__ int d_interior_face_dir(int wallFace)
{
    switch (wallFace)
    {
    case 0: return 2;
    case 1: return 3;
    case 2: return 0;
    case 3: return 1;
    default: return -1;
    }
}

__device__ void d_load_cell_U(int cell, const double *d_U, int n_cells, double U[kNv])
{
    if (cell < 0 || cell >= n_cells)
    {
        for (int v = 0; v < kNv; ++v)
            U[v] = 0.0;
        return;
    }
    for (int v = 0; v < kNv; ++v)
        U[v] = d_U[cell * kNv + v];
}

__device__ void d_load_offset_state(int cell0, int negFace, int posFace, int offset,
                                   const int *d_neighbours, const int *d_wall_face,
                                   const double *d_face_normals, const double *d_U,
                                   int n_phys, int n_cells, double gamma, double U[kNv])
{
    if (offset == 0)
    {
        d_load_cell_U(cell0, d_U, n_cells, U);
        return;
    }
    const int dir = (offset > 0) ? posFace : negFace;
    const int nsteps = (offset > 0) ? offset : -offset;
    int cur = cell0;
    for (int s = 0; s < nsteps; ++s)
    {
        int nxt = -1;
        if (cur >= 0 && cur < n_phys)
            nxt = d_neighbours[cur * kMaxFaces + dir];
        if (nxt < 0)
        {
            d_load_cell_U(cur >= 0 ? cur : cell0, d_U, n_cells, U);
            return;
        }
        if (nxt >= n_phys)
        {
            const int ghost_depth = nsteps - s;
            const int bdry_cell = cur;
            const bool is_wall = (bdry_cell >= 0 && bdry_cell < n_cells &&
                                  d_wall_face[bdry_cell * kMaxFaces + dir] != 0);
            if (!is_wall)
            {
                if (ghost_depth == 1)
                    d_load_cell_U(nxt, d_U, n_cells, U);
                else
                    d_load_cell_U(bdry_cell, d_U, n_cells, U);
                return;
            }
            const int idir = d_interior_face_dir(dir);
            int refl_src = bdry_cell;
            for (int d = 0; d < ghost_depth - 1; ++d)
            {
                const int n = (refl_src >= 0 && refl_src < n_phys)
                                  ? d_neighbours[refl_src * kMaxFaces + idir]
                                  : -1;
                if (n < 0 || n >= n_phys)
                    break;
                refl_src = n;
            }
            double U_src[kNv];
            d_load_cell_U(refl_src, d_U, n_cells, U_src);
            double nx = d_face_normals[(bdry_cell * kMaxFaces + dir) * 2 + 0];
            double ny = d_face_normals[(bdry_cell * kMaxFaces + dir) * 2 + 1];
            const double nmag = sqrt(nx * nx + ny * ny);
            if (nmag > 1e-14)
                d_slip_reflect_U(U_src, nx / nmag, ny / nmag, gamma, U);
            else
                d_load_cell_U(nxt, d_U, n_cells, U);
            return;
        }
        cur = nxt;
    }
    d_load_cell_U(cur, d_U, n_cells, U);
}

__device__ void d_mat4_vec(const double A[16], const double x[kNv], double y[kNv])
{
    for (int i = 0; i < kNv; ++i)
    {
        double s = 0.0;
        for (int j = 0; j < kNv; ++j)
            s += A[i + kNv * j] * x[j];
        y[i] = s;
    }
}

/** Roe-averaged L / L^{-1} for face-normal characteristic WENO (host Get_LR). */
__device__ bool d_get_LR(double dL, double uL, double vL, double aL,
                         double dR, double uR, double vR, double aR,
                         double nx, double ny, double gamma,
                         double L[16], double IL[16])
{
    const double sqrt_dL = sqrt(fmax(dL, 1e-14));
    const double sqrt_dR = sqrt(fmax(dR, 1e-14));
    const double denom = sqrt_dL + sqrt_dR;
    double u_RL, v_RL, a_RL;
    if (denom < 1e-14)
    {
        u_RL = 0.5 * (uL + uR);
        v_RL = 0.5 * (vL + vR);
        a_RL = 0.5 * (aL + aR);
    }
    else
    {
        u_RL = (uL * sqrt_dL + uR * sqrt_dR) / denom;
        v_RL = (vL * sqrt_dL + vR * sqrt_dR) / denom;
        a_RL = (aL * sqrt_dL + aR * sqrt_dR) / denom;
    }
    if (!(a_RL > 1e-14) || !isfinite(a_RL))
        a_RL = 1e-14;

    const double gamma_m1 = gamma - 1.0;
    const double vn = u_RL * nx + v_RL * ny;
    const double ek = 0.5 * (u_RL * u_RL + v_RL * v_RL);
    const double h = (a_RL * a_RL / gamma_m1) + ek;
    const double t1 = 0.5 / (a_RL * a_RL);
    const double t2 = gamma_m1 * t1;

    L[0] = 1.0 - 2.0 * t2 * ek;
    L[4] = 2.0 * t2 * u_RL;
    L[8] = 2.0 * t2 * v_RL;
    L[12] = -2.0 * t2;
    L[1] = v_RL * nx - u_RL * ny;
    L[5] = ny;
    L[9] = -nx;
    L[13] = 0.0;
    L[2] = t2 * ek - a_RL * vn * t1;
    L[6] = -t2 * u_RL + a_RL * t1 * nx;
    L[10] = -t2 * v_RL + a_RL * t1 * ny;
    L[14] = t2;
    L[3] = t2 * ek + a_RL * vn * t1;
    L[7] = -t2 * u_RL - a_RL * t1 * nx;
    L[11] = -t2 * v_RL - a_RL * t1 * ny;
    L[15] = t2;

    IL[0] = 1.0;
    IL[4] = 0.0;
    IL[8] = 1.0;
    IL[12] = 1.0;
    IL[1] = u_RL;
    IL[5] = ny;
    IL[9] = u_RL + a_RL * nx;
    IL[13] = u_RL - a_RL * nx;
    IL[2] = v_RL;
    IL[6] = -nx;
    IL[10] = v_RL + a_RL * ny;
    IL[14] = v_RL - a_RL * ny;
    IL[3] = ek;
    IL[7] = u_RL * ny - v_RL * nx;
    IL[11] = h + a_RL * vn;
    IL[15] = h - a_RL * vn;

    for (int i = 0; i < 16; ++i)
    {
        if (!isfinite(L[i]) || !isfinite(IL[i]))
            return false;
    }
    return true;
}

__device__ void d_load_prim6(int cell, const double *d_prim, int n_cells,
                             double &rho, double &u, double &v, double &a)
{
    if (cell < 0 || cell >= n_cells || d_prim == nullptr)
    {
        rho = 1.0;
        u = 0.0;
        v = 0.0;
        a = 1.0;
        return;
    }
    rho = d_prim[cell * 6 + 0];
    u = d_prim[cell * 6 + 1];
    v = d_prim[cell * 6 + 2];
    a = d_prim[cell * 6 + 5];
    if (!(rho > 1e-14) || !isfinite(rho))
        rho = 1e-14;
    if (!isfinite(u))
        u = 0.0;
    if (!isfinite(v))
        v = 0.0;
    if (!(a > 1e-14) || !isfinite(a))
        a = 1e-14;
}

__device__ void d_weno_from_offset(int cell, int face, int neigh, int negFace, int posFace,
                                  int o1, int o2, int o3, int o4, int o5, int LR,
                                  const int *d_neighbours, const int *d_wall_face,
                                  const double *d_face_normals, const double *d_U,
                                  const double *d_prim,
                                  int n_phys, int n_cells, double gamma, bool is_char,
                                  double Uout[kNv])
{
    double U1[kNv], U2[kNv], U3[kNv], U4[kNv], U5[kNv];
    d_load_offset_state(cell, negFace, posFace, o1, d_neighbours, d_wall_face, d_face_normals, d_U, n_phys, n_cells, gamma, U1);
    d_load_offset_state(cell, negFace, posFace, o2, d_neighbours, d_wall_face, d_face_normals, d_U, n_phys, n_cells, gamma, U2);
    d_load_offset_state(cell, negFace, posFace, o3, d_neighbours, d_wall_face, d_face_normals, d_U, n_phys, n_cells, gamma, U3);
    d_load_offset_state(cell, negFace, posFace, o4, d_neighbours, d_wall_face, d_face_normals, d_U, n_phys, n_cells, gamma, U4);
    d_load_offset_state(cell, negFace, posFace, o5, d_neighbours, d_wall_face, d_face_normals, d_U, n_phys, n_cells, gamma, U5);

    if (!is_char)
    {
        for (int v = 0; v < kNv; ++v)
            Uout[v] = d_weno5_scalar(U1[v], U2[v], U3[v], U4[v], U5[v], LR);
        return;
    }

    double dL, uL, vL, aL, dR, uR, vR, aR;
    d_load_prim6(cell, d_prim, n_cells, dL, uL, vL, aL);
    d_load_prim6(neigh, d_prim, n_cells, dR, uR, vR, aR);
    double nx = d_face_normals[(cell * kMaxFaces + face) * 2 + 0];
    double ny = d_face_normals[(cell * kMaxFaces + face) * 2 + 1];
    double Lmat[16], ILmat[16];
    if (!d_get_LR(dL, uL, vL, aL, dR, uR, vR, aR, nx, ny, gamma, Lmat, ILmat))
    {
        for (int v = 0; v < kNv; ++v)
            Uout[v] = d_weno5_scalar(U1[v], U2[v], U3[v], U4[v], U5[v], LR);
        return;
    }

    double W1[kNv], W2[kNv], W3[kNv], W4[kNv], W5[kNv], W[kNv];
    d_mat4_vec(Lmat, U1, W1);
    d_mat4_vec(Lmat, U2, W2);
    d_mat4_vec(Lmat, U3, W3);
    d_mat4_vec(Lmat, U4, W4);
    d_mat4_vec(Lmat, U5, W5);
    for (int v = 0; v < kNv; ++v)
        W[v] = d_weno5_scalar(W1[v], W2[v], W3[v], W4[v], W5[v], LR);
    d_mat4_vec(ILmat, W, Uout);
    for (int v = 0; v < kNv; ++v)
    {
        if (!isfinite(Uout[v]))
            Uout[v] = U3[v];
    }
}

__device__ double d_pressure_from_U(const double U[kNv], double gamma)
{
    const double rho = U[0];
    if (!(rho > 1e-14) || !isfinite(rho))
        return -1.0;
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    return (gamma - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v));
}

__device__ double d_blend_pressure(const double U0[kNv], const double Uw[kNv],
                                  double t, double gamma)
{
    double U[kNv];
    for (int i = 0; i < kNv; ++i)
        U[i] = (1.0 - t) * U0[i] + t * Uw[i];
    return d_pressure_from_U(U, gamma);
}

/**
 * Zhang–Shu positivity limiter (convex blend):
 *   U* = (1-θ) U_1O + θ U_WENO,  θ ∈ [0,1]
 * Largest θ such that ρ(U*)≥ε and p(U*)≥ε.
 */
__device__ double d_positivity_theta(const double U0[kNv], const double Uw[kNv],
                                    double gamma, double eps)
{
    double theta = 1.0;

    if (!(Uw[0] > eps) || !isfinite(Uw[0]))
    {
        if (U0[0] > eps && isfinite(U0[0]) && (U0[0] - Uw[0]) > 1e-30)
            theta = fmin(theta, (U0[0] - eps) / (U0[0] - Uw[0]));
        else
            theta = 0.0;
    }
    theta = fmax(0.0, fmin(1.0, theta));

    if (!(d_blend_pressure(U0, Uw, theta, gamma) > eps))
    {
        double lo = 0.0, hi = theta;
        for (int it = 0; it < 14; ++it)
        {
            const double mid = 0.5 * (lo + hi);
            if (d_blend_pressure(U0, Uw, mid, gamma) > eps)
                lo = mid;
            else
                hi = mid;
        }
        theta = lo;
    }
    if (!isfinite(theta) || theta < 0.0)
        theta = 0.0;
    if (theta > 1.0)
        theta = 1.0;
    return theta;
}

__device__ void d_positivity_limit_face(const double U_cell[kNv], double U_face[kNv],
                                       double gamma, double eps)
{
    const double th = d_positivity_theta(U_cell, U_face, gamma, eps);
    if (th >= 1.0 - 1e-15)
        return;
    for (int i = 0; i < kNv; ++i)
        U_face[i] = (1.0 - th) * U_cell[i] + th * U_face[i];
}

__device__ bool d_weno_wall_face(int cell, int face,
                                const int *d_neighbours, const int *d_wall_face,
                                const double *d_face_normals, const double *d_U,
                                int n_phys, int n_cells, double gamma,
                                double U_L[kNv], double U_R[kNv])
{
    const int idir = d_interior_face_dir(face);
    if (idir < 0 || cell < 0 || cell >= n_phys)
        return false;
    double nx = d_face_normals[(cell * kMaxFaces + face) * 2 + 0];
    double ny = d_face_normals[(cell * kMaxFaces + face) * 2 + 1];
    const double nmag = sqrt(nx * nx + ny * ny);
    if (!(nmag > 1e-14))
        return false;
    nx /= nmag;
    ny /= nmag;

    double C0[kNv], C1[kNv], C2[kNv];
    d_load_cell_U(cell, d_U, n_cells, C0);
    const int c1 = d_neighbours[cell * kMaxFaces + idir];
    if (c1 < 0 || c1 >= n_phys)
        return false;
    d_load_cell_U(c1, d_U, n_cells, C1);
    int c2 = d_neighbours[c1 * kMaxFaces + idir];
    if (c2 < 0 || c2 >= n_phys)
    {
        c2 = c1;
        for (int v = 0; v < kNv; ++v)
            C2[v] = C1[v];
    }
    else
        d_load_cell_U(c2, d_U, n_cells, C2);

    double G1[kNv], G2[kNv], G3[kNv];
    d_slip_reflect_U(C0, nx, ny, gamma, G1);
    d_slip_reflect_U(C1, nx, ny, gamma, G2);
    d_slip_reflect_U(C2, nx, ny, gamma, G3);

    double U_a[kNv], U_b[kNv];
    for (int v = 0; v < kNv; ++v)
    {
        U_a[v] = d_weno5_scalar(G3[v], G2[v], G1[v], C0[v], C1[v], 0);
        U_b[v] = d_weno5_scalar(G2[v], G1[v], C0[v], C1[v], C2[v], 1);
    }
    if (face == 0 || face == 1)
    {
        for (int v = 0; v < kNv; ++v)
        {
            U_L[v] = U_b[v];
            U_R[v] = U_a[v];
        }
    }
    else
    {
        for (int v = 0; v < kNv; ++v)
        {
            U_L[v] = U_a[v];
            U_R[v] = U_b[v];
        }
    }
    const double pL = d_pressure_from_U(U_L, gamma);
    if (isfinite(pL) && pL > 0.0 && U_R[0] > 1e-14)
    {
        const double u = U_R[1] / U_R[0];
        const double vv = U_R[2] / U_R[0];
        U_R[3] = pL / (gamma - 1.0) + U_R[0] * 0.5 * (u * u + vv * vv);
    }
    return isfinite(U_L[0]) && isfinite(U_R[0]) && U_L[0] > 0.0 && U_R[0] > 0.0;
}

__device__ void d_weno_interior_face(int cell, int face,
                                    const int *d_neighbours, const int *d_wall_face,
                                    const double *d_face_normals, const double *d_U,
                                    const double *d_prim,
                                    int n_phys, int n_cells, double gamma, bool is_char,
                                    double U_L[kNv], double U_R[kNv])
{
    const int neigh = d_neighbours[cell * kMaxFaces + face];
    switch (face)
    {
    case 0:
        d_weno_from_offset(cell, face, neigh, 0, 2, -3, -2, -1, 0, 1, 0, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_L);
        d_weno_from_offset(cell, face, neigh, 0, 2, -2, -1, 0, 1, 2, 1, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_R);
        break;
    case 1:
        d_weno_from_offset(cell, face, neigh, 1, 3, -3, -2, -1, 0, 1, 0, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_L);
        d_weno_from_offset(cell, face, neigh, 1, 3, -2, -1, 0, 1, 2, 1, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_R);
        break;
    case 2:
        d_weno_from_offset(cell, face, neigh, 0, 2, -2, -1, 0, 1, 2, 0, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_L);
        d_weno_from_offset(cell, face, neigh, 0, 2, -1, 0, 1, 2, 3, 1, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_R);
        break;
    case 3:
        d_weno_from_offset(cell, face, neigh, 1, 3, -2, -1, 0, 1, 2, 0, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_L);
        d_weno_from_offset(cell, face, neigh, 1, 3, -1, 0, 1, 2, 3, 1, d_neighbours, d_wall_face, d_face_normals, d_U, d_prim, n_phys, n_cells, gamma, is_char, U_R);
        break;
    default:
        d_load_cell_U(cell, d_U, n_cells, U_L);
        d_load_cell_U(cell, d_U, n_cells, U_R);
        break;
    }
    if (face == 0 || face == 1)
    {
        for (int v = 0; v < kNv; ++v)
        {
            const double t = U_L[v];
            U_L[v] = U_R[v];
            U_R[v] = t;
        }
    }
}

__global__ void net_flux_weno_ricca_kernel(
    const int *leaf_cells,
    int n_leaf,
    const int *d_num_faces,
    const int *d_neighbours,
    const double *d_face_normals,
    const double *d_face_areas,
    const int *d_wall_face,
    const double *d_U,
    const double *d_prim,
    double *d_net_flux,
    int n_cells,
    int n_phys,
    int dissipation_type,
    bool is_char,
    double gamma,
    double gamma1)
{
    const int li = blockIdx.x * blockDim.x + threadIdx.x;
    if (li >= n_leaf)
        return;
    const int cell = leaf_cells[li];
    if (cell < 0 || cell >= n_phys)
        return;

    for (int v = 0; v < kNv; ++v)
        d_net_flux[cell * kNv + v] = 0.0;

    const int nf = d_num_faces[cell];
    if (nf != 4)
        return;

    for (int face = 0; face < 4; ++face)
    {
        const int neigh = d_neighbours[cell * kMaxFaces + face];
        if (neigh < 0 || neigh >= n_cells)
            continue;

        const bool is_wall = d_wall_face[cell * kMaxFaces + face] != 0;
        double U_L[kNv], U_R[kNv];
        double U_cell[kNv], U_neigh[kNv];
        d_load_cell_U(cell, d_U, n_cells, U_cell);
        d_load_cell_U(neigh, d_U, n_cells, U_neigh);

        /* Host WENO2D: walls / invalid neighbours use 1O cell averages (no WENO). */
        if (is_wall || neigh >= n_phys)
        {
            for (int v = 0; v < kNv; ++v)
            {
                U_L[v] = U_cell[v];
                U_R[v] = U_neigh[v];
            }
        }
        else
        {
            d_weno_interior_face(cell, face, d_neighbours, d_wall_face, d_face_normals, d_U,
                                 d_prim, n_phys, n_cells, gamma, is_char, U_L, U_R);
        }

        /* Zhang–Shu: convex-blend WENO faces toward 1O cell averages if ρ/p ≤ ε. */
        {
            constexpr double pos_eps = 1e-8;
            d_positivity_limit_face(U_cell, U_L, gamma, pos_eps);
            d_positivity_limit_face(U_neigh, U_R, gamma, pos_eps);
            const double pL = d_pressure_from_U(U_L, gamma);
            const double pR = d_pressure_from_U(U_R, gamma);
            if (!(U_L[0] > pos_eps) || !(U_R[0] > pos_eps) || !(pL > pos_eps) || !(pR > pos_eps) ||
                !isfinite(pL) || !isfinite(pR))
            {
                for (int v = 0; v < kNv; ++v)
                {
                    U_L[v] = U_cell[v];
                    U_R[v] = U_neigh[v];
                }
            }
        }

        double rho_L = U_L[0], u_L = U_L[1] / fmax(U_L[0], 1e-14), v_L = U_L[2] / fmax(U_L[0], 1e-14);
        double rho_R = U_R[0], u_R = U_R[1] / fmax(U_R[0], 1e-14), v_R = U_R[2] / fmax(U_R[0], 1e-14);
        double P_L = d_pressure_from_U(U_L, gamma);
        double P_R = d_pressure_from_U(U_R, gamma);
        if (!(P_L > 1e-12))
            P_L = 1e-12;
        if (!(P_R > 1e-12))
            P_R = 1e-12;

        const double nx = d_face_normals[(cell * kMaxFaces + face) * 2 + 0];
        const double ny = d_face_normals[(cell * kMaxFaces + face) * 2 + 1];
        const double dl = d_face_areas[cell * kMaxFaces + face];
        if (!(dl > 0.0))
            continue;

        if (is_wall)
            P_R = P_L;

        double avg[kNv], diss[kNv];
        face_average_convective(rho_L, u_L, v_L, P_L, rho_R, u_R, v_R, P_R, nx, ny, dl, is_wall, gamma1, avg);

        if (dissipation_type == 4)
        {
            /* RICCA on reconstructed face conservatives (host WENO path). */
            const double Vdotn_L = u_L * nx + v_L * ny;
            const double Vdotn_R = u_R * nx + v_R * ny;
            const double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
            const double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);
            double d_Uv[kNv] = {U_R[0] - U_L[0], U_R[1] - U_L[1], U_R[2] - U_L[2], U_R[3] - U_L[3]};
            /* Match host WENO2D / 1O RICCA: dF = F_R - F_L with +P on both sides. */
            const double d_F[kNv] = {
                (rho_R * Vdotn_R - rho_L * Vdotn_L) * dl,
                (rho_R * u_R * Vdotn_R + P_R * nx - rho_L * u_L * Vdotn_L + P_L * nx) * dl,
                (rho_R * v_R * Vdotn_R + P_R * ny - rho_L * v_L * Vdotn_L + P_L * ny) * dl,
                ((((P_R / (gamma - 1.0)) + rho_R * Vmag_R) + P_R) * Vdotn_R -
                 (((P_L / (gamma - 1.0)) + rho_L * Vmag_L) + P_L) * Vdotn_L) *
                    dl};
            double dP = fabs(P_L - P_R);
            const double P_I = 0.5 * (P_L + P_R);
            const double Rho_I = 0.5 * (rho_L + rho_R);
            for (int i = 0; i < kNv; ++i)
            {
                double alpha = 0.0;
                condition_for_ricca(d_Uv[i], d_F[i], Vdotn_L, Vdotn_R, dP, Rho_I, P_I, alpha, gamma);
                diss[i] = 0.5 * alpha * dl * d_Uv[i];
                if (!isfinite(diss[i]))
                    diss[i] = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < kNv; ++i)
                diss[i] = 0.0;
        }

        for (int v = 0; v < kNv; ++v)
            d_net_flux[cell * kNv + v] += avg[v] - diss[v];
    }
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

__device__ void conserv_to_prim(const double *U, double *prim, double gamma)
{
    double rho = U[0];
    if (rho < 1e-12)
        rho = 1e-12;
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    const double ke = 0.5 * (u * u + v * v);
    double p = (gamma - 1.0) * (U[3] - rho * ke);
    if (p < 1e-12)
        p = 1e-12;
    prim[0] = rho;
    prim[1] = u;
    prim[2] = v;
    prim[3] = 0.0;
    prim[4] = p;
    prim[5] = sqrt(gamma * p / rho);
}

/* Repair conservatives in-place so prim and U stay thermodynamically consistent.
   Blind p-clamp in conserv_to_prim alone leaves U[3] negative-internal-energy and
   the next WENO step amplifies until |V|→∞ and dt collapses. */
__device__ void repair_conserved(double *U, double gamma)
{
    const double rho_floor = 1e-8;
    const double p_floor = 1e-8;
    double rho = U[0];
    const bool rho_bad = !isfinite(rho) || rho < rho_floor;
    if (rho_bad)
    {
        /* Tiny/NaN density with leftover energy → a→∞; reset to a calm floor state. */
        U[0] = rho_floor;
        U[1] = 0.0;
        U[2] = 0.0;
        U[3] = p_floor / (gamma - 1.0);
        return;
    }
    double u = U[1] / rho;
    double v = U[2] / rho;
    if (!isfinite(u))
        u = 0.0;
    if (!isfinite(v))
        v = 0.0;
    /* Cap extreme velocities that already indicate a failed update. */
    const double v_cap = 20.0;
    const double speed = sqrt(u * u + v * v);
    if (speed > v_cap)
    {
        const double s = v_cap / speed;
        u *= s;
        v *= s;
    }
    U[1] = rho * u;
    U[2] = rho * v;
    double ke = 0.5 * (u * u + v * v);
    double p = (gamma - 1.0) * (U[3] - rho * ke);
    if (!isfinite(p) || p < p_floor)
    {
        p = p_floor;
        /* Near-vacuum: drop unsupported KE so a stays finite. */
        if (ke > 2.0)
        {
            u = 0.0;
            v = 0.0;
            ke = 0.0;
            U[1] = 0.0;
            U[2] = 0.0;
        }
    }
    U[3] = p / (gamma - 1.0) + rho * ke;
}

__global__ void inlet_ghost_kernel(const int *list, int n, double *U, double *prim,
                                   double rho, double u, double v, double p, double gamma)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    const int g = list[i * 3 + 2];
    const double ke = 0.5 * (u * u + v * v);
    U[g * kNv + 0] = rho;
    U[g * kNv + 1] = rho * u;
    U[g * kNv + 2] = rho * v;
    U[g * kNv + 3] = p / (gamma - 1.0) + rho * ke;
    conserv_to_prim(&U[g * kNv], &prim[g * 6], gamma);
}

__global__ void exit_ghost_kernel(const int *list, int n, double *U, double *prim, double gamma)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    const int c = list[i * 3 + 0];
    const int g = list[i * 3 + 2];
    for (int v = 0; v < kNv; ++v)
        U[g * kNv + v] = U[c * kNv + v];
    conserv_to_prim(&U[g * kNv], &prim[g * 6], gamma);
}

__global__ void wall_ghost_kernel(const int *list, int n, double *U, double *prim,
                                  const double *face_normals, double gamma)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    const int c = list[i * 3 + 0];
    const int f = list[i * 3 + 1];
    const int g = list[i * 3 + 2];
    const double n1 = face_normals[(c * kMaxFaces + f) * 2 + 0];
    const double n2 = face_normals[(c * kMaxFaces + f) * 2 + 1];
    const double rho = prim[c * 6 + 0];
    const double u = prim[c * 6 + 1];
    const double v = prim[c * 6 + 2];
    const double p = prim[c * 6 + 4];
    const double un = u * n1 + v * n2;
    const double ug = u - 2.0 * un * n1;
    const double vg = v - 2.0 * un * n2;
    const double ke = 0.5 * (ug * ug + vg * vg);
    U[g * kNv + 0] = rho;
    U[g * kNv + 1] = rho * ug;
    U[g * kNv + 2] = rho * vg;
    U[g * kNv + 3] = p / (gamma - 1.0) + rho * ke;
    conserv_to_prim(&U[g * kNv], &prim[g * 6], gamma);
}

__global__ void cell_dt_kernel(const int *leaf, int n_leaf, const int *num_faces,
                               const double *face_normals, const double *face_areas,
                               const double *prim, const double *inv_area_unused,
                               const double *cell_area, double *dt_out, double cfl)
{
    (void)inv_area_unused;
    const int li = blockIdx.x * blockDim.x + threadIdx.x;
    if (li >= n_leaf)
        return;
    const int c = leaf[li];
    const double u = prim[c * 6 + 1];
    const double v = prim[c * 6 + 2];
    const double a = prim[c * 6 + 5];
    const int nf = num_faces[c];
    double denom = 0.0;
    for (int f = 0; f < nf && f < kMaxFaces; ++f)
    {
        const double nx = face_normals[(c * kMaxFaces + f) * 2 + 0];
        const double ny = face_normals[(c * kMaxFaces + f) * 2 + 1];
        const double dl = face_areas[c * kMaxFaces + f];
        denom += (fabs(u * nx + v * ny) + a) * dl;
    }
    dt_out[li] = (denom > 1e-30) ? (cfl * cell_area[c] / denom) : 0.0;
}

__global__ void update_and_error_kernel(const int *leaf, int n_leaf, double *U, double *prim,
                                        const double *net_flux, const double *inv_area,
                                        double min_dt, double *err_sq, double gamma)
{
    const int li = blockIdx.x * blockDim.x + threadIdx.x;
    if (li >= n_leaf)
        return;
    const int c = leaf[li];
    const double s = -min_dt * inv_area[c];
    double local_e[kNv] = {0.0, 0.0, 0.0, 0.0};
    double Uold[kNv], Utry[kNv];
    for (int v = 0; v < kNv; ++v)
        Uold[v] = U[c * kNv + v];
    for (int v = 0; v < kNv; ++v)
    {
        const double dU = s * net_flux[c * kNv + v];
        Utry[v] = Uold[v] + dU;
    }
    /* Zhang–Shu cell update: U ← (1-θ)Uold + θ Utry with ρ,p ≥ ε (θ≤1). */
    double th = 1.0;
    {
        constexpr double pos_eps = 1e-8;
        for (int v = 0; v < kNv; ++v)
        {
            if (!isfinite(Utry[v]))
            {
                th = 0.0;
                break;
            }
        }
        if (th > 0.0)
            th = d_positivity_theta(Uold, Utry, gamma, pos_eps);
        for (int v = 0; v < kNv; ++v)
            U[c * kNv + v] = (1.0 - th) * Uold[v] + th * Utry[v];
    }
    /* Residual from *applied* update (θ·dU), not the rejected trial. */
    for (int v = 0; v < kNv; ++v)
    {
        const double dU_app = th * (Utry[v] - Uold[v]);
        const double denom = Uold[v];
        double t = (!isfinite(dU_app) || !isfinite(denom)) ? 0.0
                                                           : (fabs(denom) < 1e-12 ? fabs(dU_app) : fabs(dU_app) / fabs(denom));
        local_e[v] = t * t;
    }
    repair_conserved(&U[c * kNv], gamma);
    /* Last-resort safety net if repair still leaves a runaway / vacuum state. */
    {
        const double rho = U[c * kNv + 0];
        bool bad = !isfinite(rho) || !(rho > 1e-10);
        if (!bad)
        {
            const double u = U[c * kNv + 1] / rho;
            const double v = U[c * kNv + 2] / rho;
            const double ke = 0.5 * (u * u + v * v);
            const double p = (gamma - 1.0) * (U[c * kNv + 3] - rho * ke);
            const double a2 = (rho > 0.0) ? (gamma * p / rho) : 0.0;
            const double speed = sqrt(u * u + v * v);
            bad = !isfinite(u) || !isfinite(v) || !isfinite(p) || !isfinite(a2) ||
                  !(p > 1e-10) || speed > 20.0 || a2 > 20.0 * 20.0;
        }
        if (bad)
        {
            for (int v = 0; v < kNv; ++v)
                U[c * kNv + v] = Uold[v];
            repair_conserved(&U[c * kNv], gamma);
        }
    }
    conserv_to_prim(&U[c * kNv], &prim[c * 6], gamma);
    for (int v = 0; v < kNv; ++v)
        atomicAdd(&err_sq[v], local_e[v]);
}

struct CudaFluxBuffers
{
    int cell_cap = 0;
    int leaf_cap = 0;
    int n_cells = 0;
    int n_leaf = 0;
    int n_phys = 0;
    bool geom_ready = false;
    bool resident = false;

    int *d_num_faces = nullptr;
    int *d_neighbours = nullptr;
    double *d_face_normals = nullptr;
    double *d_face_areas = nullptr;
    int *d_wall_face = nullptr;
    double *d_prim = nullptr;
    double *d_U = nullptr;
    double *d_net_flux = nullptr;
    int *d_leaf = nullptr;
    double *d_inv_area = nullptr;
    double *d_cell_area = nullptr;
    double *d_dt = nullptr;
    double *d_err_sq = nullptr;
    int *d_inlet = nullptr;
    int *d_exit = nullptr;
    int *d_wall = nullptr;
    int n_inlet = 0, n_exit = 0, n_wall = 0;
    std::vector<int> h_leaf;

    void free_all()
    {
        cudaFree(d_num_faces);
        cudaFree(d_neighbours);
        cudaFree(d_face_normals);
        cudaFree(d_face_areas);
        cudaFree(d_wall_face);
        cudaFree(d_prim);
        cudaFree(d_U);
        cudaFree(d_net_flux);
        cudaFree(d_leaf);
        cudaFree(d_inv_area);
        cudaFree(d_cell_area);
        cudaFree(d_dt);
        cudaFree(d_err_sq);
        cudaFree(d_inlet);
        cudaFree(d_exit);
        cudaFree(d_wall);
        *this = CudaFluxBuffers{};
    }
};

CudaFluxBuffers g_buf;

bool cuda_ok(cudaError_t e, const char *msg)
{
    if (e != cudaSuccess)
    {
        fprintf(stderr, "MOVERS/RICCA CUDA: %s (%s)\n", msg, cudaGetErrorString(e));
        return false;
    }
    return true;
}

bool ensure_capacity(int n_cells, int n_leaf)
{
    if (n_cells <= g_buf.cell_cap && n_leaf <= g_buf.leaf_cap && g_buf.d_num_faces)
        return true;

    g_buf.free_all();
    g_buf.cell_cap = n_cells;
    g_buf.leaf_cap = n_leaf;

    return cuda_ok(cudaMalloc(&g_buf.d_num_faces, n_cells * sizeof(int)), "malloc num_faces") &&
           cuda_ok(cudaMalloc(&g_buf.d_neighbours, n_cells * kMaxFaces * sizeof(int)), "malloc neighbours") &&
           cuda_ok(cudaMalloc(&g_buf.d_face_normals, n_cells * kMaxFaces * 2 * sizeof(double)), "malloc normals") &&
           cuda_ok(cudaMalloc(&g_buf.d_face_areas, n_cells * kMaxFaces * sizeof(double)), "malloc areas") &&
           cuda_ok(cudaMalloc(&g_buf.d_wall_face, n_cells * kMaxFaces * sizeof(int)), "malloc wall") &&
           cuda_ok(cudaMalloc(&g_buf.d_prim, n_cells * 6 * sizeof(double)), "malloc prim") &&
           cuda_ok(cudaMalloc(&g_buf.d_U, n_cells * kNv * sizeof(double)), "malloc U") &&
           cuda_ok(cudaMalloc(&g_buf.d_net_flux, n_cells * kNv * sizeof(double)), "malloc net_flux") &&
           cuda_ok(cudaMalloc(&g_buf.d_leaf, n_leaf * sizeof(int)), "malloc leaf") &&
           cuda_ok(cudaMalloc(&g_buf.d_inv_area, n_cells * sizeof(double)), "malloc inv_area") &&
           cuda_ok(cudaMalloc(&g_buf.d_cell_area, n_cells * sizeof(double)), "malloc cell_area") &&
           cuda_ok(cudaMalloc(&g_buf.d_dt, n_leaf * sizeof(double)), "malloc dt") &&
           cuda_ok(cudaMalloc(&g_buf.d_err_sq, 4 * sizeof(double)), "malloc err");
}

void pack_geometry_host(int n_cells,
                        std::vector<int> &h_num_faces,
                        std::vector<int> &h_neighbours,
                        std::vector<double> &h_normals,
                        std::vector<double> &h_areas,
                        std::vector<int> &h_wall,
                        std::vector<double> &h_inv_area,
                        std::vector<double> &h_cell_area)
{
    h_num_faces.assign(n_cells, 0);
    h_neighbours.assign(n_cells * kMaxFaces, -1);
    h_normals.assign(n_cells * kMaxFaces * 2, 0.0);
    h_areas.assign(n_cells * kMaxFaces, 0.0);
    h_wall.assign(n_cells * kMaxFaces, 0);
    h_inv_area.assign(n_cells, 0.0);
    h_cell_area.assign(n_cells, 0.0);

    for (int c = 0; c < n_cells; ++c)
    {
        const int nf = (Cells[c].numFaces > 0)
                           ? Cells[c].numFaces
                           : static_cast<int>(Cells[c].Face_Areas.size());
        h_num_faces[c] = std::min(nf, kMaxFaces);
        h_cell_area[c] = Cells[c].Area;
        h_inv_area[c] = Cells[c].Inv_Area;

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
    }
}

void pack_state_host(int n_cells, std::vector<double> &h_prim, std::vector<double> &h_U)
{
    h_prim.assign(n_cells * 6, 0.0);
    h_U.assign(n_cells * kNv, 0.0);
    for (int c = 0; c < n_cells; ++c)
    {
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
}

bool upload_geometry_once(int n_cells, int n_leaf, const V_I &leafCells)
{
    if (g_buf.geom_ready && g_buf.n_cells == n_cells && g_buf.n_leaf == n_leaf)
        return true;

    std::vector<int> h_num_faces, h_neighbours, h_wall;
    std::vector<double> h_normals, h_areas, h_inv_area, h_cell_area;
    pack_geometry_host(n_cells, h_num_faces, h_neighbours, h_normals, h_areas, h_wall, h_inv_area, h_cell_area);

    if (!cuda_ok(cudaMemcpy(g_buf.d_num_faces, h_num_faces.data(), n_cells * sizeof(int), cudaMemcpyHostToDevice), "H2D num_faces") ||
        !cuda_ok(cudaMemcpy(g_buf.d_neighbours, h_neighbours.data(), n_cells * kMaxFaces * sizeof(int), cudaMemcpyHostToDevice), "H2D neigh") ||
        !cuda_ok(cudaMemcpy(g_buf.d_face_normals, h_normals.data(), n_cells * kMaxFaces * 2 * sizeof(double), cudaMemcpyHostToDevice), "H2D normals") ||
        !cuda_ok(cudaMemcpy(g_buf.d_face_areas, h_areas.data(), n_cells * kMaxFaces * sizeof(double), cudaMemcpyHostToDevice), "H2D areas") ||
        !cuda_ok(cudaMemcpy(g_buf.d_wall_face, h_wall.data(), n_cells * kMaxFaces * sizeof(int), cudaMemcpyHostToDevice), "H2D wall") ||
        !cuda_ok(cudaMemcpy(g_buf.d_inv_area, h_inv_area.data(), n_cells * sizeof(double), cudaMemcpyHostToDevice), "H2D inv_area") ||
        !cuda_ok(cudaMemcpy(g_buf.d_cell_area, h_cell_area.data(), n_cells * sizeof(double), cudaMemcpyHostToDevice), "H2D cell_area") ||
        !cuda_ok(cudaMemcpy(g_buf.d_leaf, leafCells.data(), n_leaf * sizeof(int), cudaMemcpyHostToDevice), "H2D leaf"))
        return false;

    g_buf.n_cells = n_cells;
    g_buf.n_leaf = n_leaf;
    g_buf.geom_ready = true;
    return true;
}

bool launch_flux_kernel(int n_cells, int n_leaf)
{
    const int threads = 256;
    const int blocks = (n_leaf + threads - 1) / threads;
    if (Is_WENO && Dissipation_Type == 4)
    {
        net_flux_weno_ricca_kernel<<<blocks, threads>>>(
            g_buf.d_leaf, n_leaf,
            g_buf.d_num_faces, g_buf.d_neighbours,
            g_buf.d_face_normals, g_buf.d_face_areas, g_buf.d_wall_face,
            g_buf.d_U, g_buf.d_prim, g_buf.d_net_flux,
            n_cells, g_buf.n_phys, Dissipation_Type,
            Is_Char, GAMMA_CONST, GAMMA1_CONST);
    }
    else
    {
        net_flux_movers_ricca_kernel<<<blocks, threads>>>(
            g_buf.d_leaf, n_leaf,
            g_buf.d_num_faces, g_buf.d_neighbours,
            g_buf.d_face_normals, g_buf.d_face_areas, g_buf.d_wall_face,
            g_buf.d_prim, g_buf.d_U, g_buf.d_net_flux,
            n_cells, Dissipation_Type,
            Is_MOVERS_1, Enable_Entropy_Fix,
            GAMMA_CONST, GAMMA1_CONST);
    }
    return cuda_ok(cudaGetLastError(), "flux launch") &&
           cuda_ok(cudaDeviceSynchronize(), "flux sync");
}

} // namespace

bool Evaluate_Cell_Net_Flux_CUDA_Movers_Ricca(bool second_order)
{
    if (second_order)
        return false;
    if (Dissipation_Type != 2 && Dissipation_Type != 4 && Dissipation_Type != 5)
        return false;
    if (g_buf.resident)
        return false; /* resident path owns device state */

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
    if (!upload_geometry_once(n_cells, n_leaf, leafCells))
        return false;

    std::vector<double> h_prim, h_U;
    pack_state_host(n_cells, h_prim, h_U);
    if (!cuda_ok(cudaMemcpy(g_buf.d_prim, h_prim.data(), n_cells * 6 * sizeof(double), cudaMemcpyHostToDevice), "H2D prim") ||
        !cuda_ok(cudaMemcpy(g_buf.d_U, h_U.data(), n_cells * kNv * sizeof(double), cudaMemcpyHostToDevice), "H2D U"))
        return false;

    if (!launch_flux_kernel(n_cells, n_leaf))
        return false;

    std::vector<double> h_net(n_cells * kNv, 0.0);
    if (!cuda_ok(cudaMemcpy(h_net.data(), g_buf.d_net_flux, n_cells * kNv * sizeof(double), cudaMemcpyDeviceToHost), "D2H net"))
        return false;

    const int nv = std::min(kNv, NUM_FLUX_COMPONENTS);
    for (int idx = 0; idx < n_leaf; ++idx)
    {
        int c = leafCells[idx];
        if (c < 0 || c >= static_cast<int>(Cells_Net_Flux.size()))
            continue;
        for (int v = 0; v < nv; ++v)
            Cells_Net_Flux[c][v] = h_net[c * kNv + v];
        Evaluate_Time_Step(c);
    }
    return true;
}

bool Resident_GPU_Explicit_Available()
{
    /* 1O MOVERS/RICCA/MOVERS_NWSC, or WENO5+RICCA (Dissipation_Type==4). */
    if (Is_Second_Order || Is_Viscous || Is_Viscous_Wall || Time_Accurate || Is_Implicit_Method)
        return false;
    if (NUM_FLUX_COMPONENTS != 4 || Cells.empty())
        return false;
    if (Is_WENO)
        return (Dissipation_Type == 4); /* GPU WENO5 + RICCA */
    return (Dissipation_Type == 2 || Dissipation_Type == 4 || Dissipation_Type == 5);
}

bool Resident_GPU_Explicit_Init()
{
    if (!Resident_GPU_Explicit_Available())
        return false;

    const int n_cells = static_cast<int>(Cells.size());
    V_I leafCells;
    Build_Leaf_Cell_List(leafCells);
    const int n_leaf = static_cast<int>(leafCells.size());
    if (n_leaf <= 0)
        return false;

    g_buf.free_all();
    if (!ensure_capacity(n_cells, n_leaf))
        return false;
    if (!upload_geometry_once(n_cells, n_leaf, leafCells))
        return false;

    std::vector<double> h_prim, h_U;
    pack_state_host(n_cells, h_prim, h_U);
    if (!cuda_ok(cudaMemcpy(g_buf.d_prim, h_prim.data(), n_cells * 6 * sizeof(double), cudaMemcpyHostToDevice), "init prim") ||
        !cuda_ok(cudaMemcpy(g_buf.d_U, h_U.data(), n_cells * kNv * sizeof(double), cudaMemcpyHostToDevice), "init U"))
        return false;

    auto upload_bc = [](const V_I &src, int *&dst, int &n) -> bool {
        n = static_cast<int>(src.size() / 3);
        if (n <= 0)
        {
            dst = nullptr;
            return true;
        }
        if (!cuda_ok(cudaMalloc(&dst, src.size() * sizeof(int)), "malloc bc"))
            return false;
        return cuda_ok(cudaMemcpy(dst, src.data(), src.size() * sizeof(int), cudaMemcpyHostToDevice), "H2D bc");
    };
    if (!upload_bc(Inlet_Cells_List, g_buf.d_inlet, g_buf.n_inlet) ||
        !upload_bc(Exit_Cells_List, g_buf.d_exit, g_buf.n_exit) ||
        !upload_bc(Wall_Cells_List, g_buf.d_wall, g_buf.n_wall))
        return false;

    g_buf.n_phys = No_Physical_Cells;
    g_buf.h_leaf.assign(leafCells.begin(), leafCells.end());
    g_buf.resident = true;
    if (CFD_MPI_Is_Root())
        std::printf("Resident GPU explicit solver ON (%s%s cells=%d leaf=%d phys=%d inlet=%d exit=%d wall=%d)\n",
                    Is_WENO ? "WENO+RICCA" : "1O",
                    (Is_WENO && Is_Char) ? "+Char" : "",
                    n_cells, n_leaf, g_buf.n_phys, g_buf.n_inlet, g_buf.n_exit, g_buf.n_wall);
    return true;
}

static void diagnose_dt_collapse(const std::vector<double> &h_dt, double min_dt, double ref_dt, const char *tag)
{
    const int n_leaf = static_cast<int>(h_dt.size());
    if (n_leaf <= 0 || g_buf.h_leaf.size() != static_cast<size_t>(n_leaf))
        return;

    int imin = 0;
    double max_dt = h_dt[0];
    int n_small = 0;
    const double thresh = (ref_dt > 0.0) ? 0.1 * ref_dt : 1e-6;
    for (int i = 0; i < n_leaf; ++i)
    {
        if (h_dt[i] < h_dt[imin])
            imin = i;
        if (h_dt[i] > max_dt)
            max_dt = h_dt[i];
        if (h_dt[i] < thresh)
            ++n_small;
    }
    const int c = g_buf.h_leaf[imin];

    std::vector<double> U4(4), prim6(6);
    cudaMemcpy(U4.data(), g_buf.d_U + c * kNv, 4 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(prim6.data(), g_buf.d_prim + c * 6, 6 * sizeof(double), cudaMemcpyDeviceToHost);

    const double rho = U4[0];
    const double u = (fabs(rho) > 1e-30) ? U4[1] / rho : 0.0;
    const double v = (fabs(rho) > 1e-30) ? U4[2] / rho : 0.0;
    const double ke = 0.5 * (u * u + v * v);
    const double praw = (GAMMA_CONST - 1.0) * (U4[3] - rho * ke);
    const double speed = sqrt(u * u + v * v);

    /* Global pathology counts (sample all leaf cells' raw pressure from U). */
    std::vector<double> h_U(g_buf.n_cells * kNv);
    std::vector<double> h_prim(g_buf.n_cells * 6);
    cudaMemcpy(h_U.data(), g_buf.d_U, g_buf.n_cells * kNv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_prim.data(), g_buf.d_prim, g_buf.n_cells * 6 * sizeof(double), cudaMemcpyDeviceToHost);

    int n_neg_p = 0, n_tiny_rho = 0, n_fast = 0, n_clamped_p = 0;
    double max_speed = 0.0, max_a = 0.0, min_praw = 1e300;
    int c_fast = -1, c_neg = -1;
    for (int li = 0; li < n_leaf; ++li)
    {
        const int ci = g_buf.h_leaf[li];
        const double r = h_U[ci * kNv + 0];
        const double uu = (fabs(r) > 1e-30) ? h_U[ci * kNv + 1] / r : 0.0;
        const double vv = (fabs(r) > 1e-30) ? h_U[ci * kNv + 2] / r : 0.0;
        const double pr = (GAMMA_CONST - 1.0) * (h_U[ci * kNv + 3] - 0.5 * r * (uu * uu + vv * vv));
        const double spd = sqrt(uu * uu + vv * vv);
        const double a = h_prim[ci * 6 + 5];
        if (!(pr > 0.0) || !std::isfinite(pr))
        {
            ++n_neg_p;
            if (c_neg < 0)
                c_neg = ci;
        }
        if (!(r > 1e-6) || !std::isfinite(r))
            ++n_tiny_rho;
        if (pr < 1e-12 && r > 1e-6)
            ++n_clamped_p; /* would hit conserv_to_prim floor */
        if (spd + a > max_speed)
        {
            max_speed = spd + a;
            c_fast = ci;
        }
        if (a > max_a)
            max_a = a;
        if (pr < min_praw)
            min_praw = pr;
        if (spd + a > 50.0)
            ++n_fast;
    }

    std::cerr << std::scientific;
    std::cerr << "[GPU-DT-DEBUG] " << tag
              << " min_dt=" << min_dt << " ref_dt=" << ref_dt
              << " max_dt=" << max_dt << " n_dt<0.1*ref=" << n_small << "\n";
    std::cerr << "[GPU-DT-DEBUG] worst-dt cell=" << c
              << " U=[" << U4[0] << "," << U4[1] << "," << U4[2] << "," << U4[3] << "]"
              << " praw=" << praw << " |V|=" << speed
              << " prim(rho,u,v,p,a)=[" << prim6[0] << "," << prim6[1] << "," << prim6[2]
              << "," << prim6[4] << "," << prim6[5] << "]\n";
    std::cerr << "[GPU-DT-DEBUG] counts: neg_p=" << n_neg_p
              << " tiny_rho=" << n_tiny_rho
              << " would_clamp_p=" << n_clamped_p
              << " fast(|V|+a>50)=" << n_fast
              << " min_praw=" << min_praw
              << " max(|V|+a)=" << max_speed << " @cell " << c_fast
              << " max_a=" << max_a
              << " first_neg_p_cell=" << c_neg << "\n";
}

bool Resident_GPU_Explicit_Step(double &min_dt_out, double err_out[4])
{
    if (!g_buf.resident)
        return false;

    const int n_cells = g_buf.n_cells;
    const int n_leaf = g_buf.n_leaf;
    const int threads = 256;

    if (!Is_Inlet_SubSonic && g_buf.n_inlet > 0)
    {
        const int blocks = (g_buf.n_inlet + threads - 1) / threads;
        inlet_ghost_kernel<<<blocks, threads>>>(g_buf.d_inlet, g_buf.n_inlet, g_buf.d_U, g_buf.d_prim,
                                                inletCond.Rho, inletCond.u, inletCond.v, inletCond.P, GAMMA_CONST);
        if (!cuda_ok(cudaGetLastError(), "inlet bc"))
            return false;
    }
    if (!Is_Exit_SubSonic && g_buf.n_exit > 0)
    {
        const int blocks = (g_buf.n_exit + threads - 1) / threads;
        exit_ghost_kernel<<<blocks, threads>>>(g_buf.d_exit, g_buf.n_exit, g_buf.d_U, g_buf.d_prim, GAMMA_CONST);
        if (!cuda_ok(cudaGetLastError(), "exit bc"))
            return false;
    }
    if (g_buf.n_wall > 0)
    {
        const int blocks = (g_buf.n_wall + threads - 1) / threads;
        wall_ghost_kernel<<<blocks, threads>>>(g_buf.d_wall, g_buf.n_wall, g_buf.d_U, g_buf.d_prim,
                                               g_buf.d_face_normals, GAMMA_CONST);
        if (!cuda_ok(cudaGetLastError(), "wall bc"))
            return false;
    }

    if (!launch_flux_kernel(n_cells, n_leaf))
        return false;

    {
        const int blocks = (n_leaf + threads - 1) / threads;
        cell_dt_kernel<<<blocks, threads>>>(g_buf.d_leaf, n_leaf, g_buf.d_num_faces,
                                            g_buf.d_face_normals, g_buf.d_face_areas,
                                            g_buf.d_prim, g_buf.d_inv_area, g_buf.d_cell_area,
                                            g_buf.d_dt, CFL);
        if (!cuda_ok(cudaGetLastError(), "dt launch"))
            return false;
    }

    std::vector<double> h_dt(n_leaf);
    if (!cuda_ok(cudaMemcpy(h_dt.data(), g_buf.d_dt, n_leaf * sizeof(double), cudaMemcpyDeviceToHost), "D2H dt"))
        return false;
    double min_dt = h_dt[0];
    for (int i = 1; i < n_leaf; ++i)
        if (h_dt[i] < min_dt)
            min_dt = h_dt[i];
    if (!(min_dt > 0.0) || !std::isfinite(min_dt))
    {
        diagnose_dt_collapse(h_dt, min_dt, -1.0, "nonfinite/nonpositive min_dt");
        return false;
    }
    /* GPU WENO can produce non-physical states that drive dt → 0 without NaN.
       Force host fallback before the solution is destroyed. */
    static double s_ref_dt = -1.0;
    static bool s_warned = false;
    if (s_ref_dt < 0.0 && min_dt > 1e-8)
        s_ref_dt = min_dt;
    if (s_ref_dt > 0.0 && min_dt < 0.5 * s_ref_dt && !s_warned)
    {
        diagnose_dt_collapse(h_dt, min_dt, s_ref_dt, "early-warning dt<0.5*ref");
        s_warned = true;
    }
    if (min_dt < 1e-12 || (s_ref_dt > 0.0 && min_dt < 1e-4 * s_ref_dt))
    {
        diagnose_dt_collapse(h_dt, min_dt, s_ref_dt, "COLLAPSE");
        std::cerr << "Resident_GPU_Explicit_Step: min_dt collapsed (" << min_dt
                  << "); requesting host fallback.\n";
        s_ref_dt = -1.0;
        s_warned = false;
        return false;
    }
    min_dt_out = min_dt;
    Min_dt = min_dt;

    double zeros[4] = {0, 0, 0, 0};
    if (!cuda_ok(cudaMemcpy(g_buf.d_err_sq, zeros, 4 * sizeof(double), cudaMemcpyHostToDevice), "clear err"))
        return false;
    {
        const int blocks = (n_leaf + threads - 1) / threads;
        update_and_error_kernel<<<blocks, threads>>>(g_buf.d_leaf, n_leaf, g_buf.d_U, g_buf.d_prim,
                                                     g_buf.d_net_flux, g_buf.d_inv_area, min_dt,
                                                     g_buf.d_err_sq, GAMMA_CONST);
        if (!cuda_ok(cudaGetLastError(), "update launch") ||
            !cuda_ok(cudaDeviceSynchronize(), "update sync"))
            return false;
    }
    double err_sq[4];
    if (!cuda_ok(cudaMemcpy(err_sq, g_buf.d_err_sq, 4 * sizeof(double), cudaMemcpyDeviceToHost), "D2H err"))
        return false;
    for (int i = 0; i < 4; ++i)
        err_out[i] = std::sqrt(err_sq[i]);
    Error[0] = err_out[0];
    Error[1] = err_out[1];
    Error[2] = err_out[2];
    Error[3] = err_out[3];
    return true;
}

bool Resident_GPU_Explicit_Download_Host()
{
    if (!g_buf.resident)
        return false;
    const int n_cells = g_buf.n_cells;
    std::vector<double> h_U(n_cells * kNv), h_prim(n_cells * 6);
    if (!cuda_ok(cudaMemcpy(h_U.data(), g_buf.d_U, n_cells * kNv * sizeof(double), cudaMemcpyDeviceToHost), "D2H U") ||
        !cuda_ok(cudaMemcpy(h_prim.data(), g_buf.d_prim, n_cells * 6 * sizeof(double), cudaMemcpyDeviceToHost), "D2H prim"))
        return false;

    for (int c = 0; c < n_cells; ++c)
    {
        if (c < static_cast<int>(U_Cells.size()))
            for (int v = 0; v < kNv && v < static_cast<int>(U_Cells[c].size()); ++v)
                U_Cells[c][v] = h_U[c * kNv + v];
        if (c < static_cast<int>(Primitive_Cells.size()) && Primitive_Cells[c].size() >= 8)
        {
            const double rho = h_prim[c * 6 + 0];
            const double u = h_prim[c * 6 + 1];
            const double v = h_prim[c * 6 + 2];
            const double p = h_prim[c * 6 + 4];
            const double a = h_prim[c * 6 + 5];
            const double speed = std::sqrt(u * u + v * v);
            const double mach = (a > 1e-14) ? (speed / a) : 0.0;
            double T = 0.0;
            if (rho > 1e-14)
            {
                if (Non_Dimensional_Form)
                    T = p / (rho * R_ref);
                else
                    T = p / (rho * R_GC);
            }
            Primitive_Cells[c][0] = rho;
            Primitive_Cells[c][1] = u;
            Primitive_Cells[c][2] = v;
            Primitive_Cells[c][3] = T;
            Primitive_Cells[c][4] = p;
            Primitive_Cells[c][5] = a;
            Primitive_Cells[c][6] = h_U[c * kNv + 3]; /* rho*E */
            Primitive_Cells[c][7] = mach;
            if (Primitive_Cells[c].size() > 10)
                Primitive_Cells[c][10] = p * std::pow(1.0 + 0.5 * GAMMA_M_1_CONST * mach * mach,
                                                      GAMMA_CONST / GAMMA_M_1_CONST);
        }
    }
    return true;
}

void Resident_GPU_Explicit_Shutdown()
{
    g_buf.free_all();
}
