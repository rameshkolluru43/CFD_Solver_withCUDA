#include "Euler3D_Cuda.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

constexpr int kNv = 5;
constexpr int kFaces = 6;
constexpr int kBlockSize = 256;
constexpr double kGamma = 1.4;
constexpr double kTiny = 1.0e-14;

#define CUDA_CHECK(call)                                                            \
    do {                                                                            \
        const cudaError_t err_ = (call);                                             \
        if (err_ != cudaSuccess) {                                                   \
            std::fprintf(stderr, "Euler3D CUDA error at %s:%d: %s\n", __FILE__,     \
                         __LINE__, cudaGetErrorString(err_));                         \
            std::abort();                                                           \
        }                                                                           \
    } while (0)

void require_size(std::size_t actual, std::size_t expected, const char *name)
{
    if (actual != expected) {
        std::fprintf(stderr, "Euler3D CUDA: %s has size %zu, expected %zu\n",
                     name, actual, expected);
        std::abort();
    }
}

__device__ __forceinline__ bool admissible(const double U[kNv])
{
    if (!isfinite(U[0]) || !(U[0] > kTiny))
        return false;
    const double ke = 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0];
    const double p = (kGamma - 1.0) * (U[4] - ke);
    return isfinite(p) && p > kTiny;
}

__device__ __forceinline__ void fill_prim(const double U[kNv], double q[6])
{
    const double rho = fmax(U[0], kTiny);
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    const double w = U[3] / rho;
    const double p = fmax((kGamma - 1.0) *
                              (U[4] - 0.5 * rho * (u * u + v * v + w * w)),
                          kTiny);
    q[0] = rho;
    q[1] = u;
    q[2] = v;
    q[3] = w;
    q[4] = p;
    q[5] = sqrt(kGamma * p / rho);
}

__device__ __forceinline__ void from_prim(double rho, double u, double v, double w,
                                          double p, double U[kNv])
{
    rho = fmax(rho, kTiny);
    p = fmax(p, kTiny);
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = rho * w;
    U[4] = p / (kGamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

__device__ __forceinline__ void physical_flux(const double U[kNv], const double q[6],
                                               double nx, double ny, double nz,
                                               double F[kNv])
{
    const double un = q[1] * nx + q[2] * ny + q[3] * nz;
    F[0] = U[0] * un;
    F[1] = U[1] * un + q[4] * nx;
    F[2] = U[2] * un + q[4] * ny;
    F[3] = U[3] * un + q[4] * nz;
    F[4] = (U[4] + q[4]) * un;
}

__device__ __forceinline__ void boundary_ghost(const double UL[kNv], int bc,
                                                double nx, double ny, double nz,
                                                const double *fs, const double *ps,
                                                double UR[kNv])
{
    if (bc == 1 || bc == 2) {
        double q[6];
        fill_prim(UL, q);
        if (bc == 1) {
            const double un = q[1] * nx + q[2] * ny + q[3] * nz;
            from_prim(q[0], q[1] - 2.0 * un * nx, q[2] - 2.0 * un * ny,
                      q[3] - 2.0 * un * nz, q[4], UR);
        } else {
            from_prim(q[0], -q[1], -q[2], -q[3], q[4], UR);
        }
    } else if (bc == 4) {
        from_prim(fs[0], fs[1], fs[2], fs[3], fs[4], UR);
    } else if (bc == 5) {
        from_prim(ps[0], ps[1], ps[2], ps[3], ps[4], UR);
    } else {
#pragma unroll
        for (int v = 0; v < kNv; ++v)
            UR[v] = UL[v];
    }
}

__device__ __forceinline__ double weno5(double v0, double v1, double v2, double v3,
                                        double v4, bool use_z)
{
    const double p0 = (2.0 * v0 - 7.0 * v1 + 11.0 * v2) / 6.0;
    const double p1 = (-v1 + 5.0 * v2 + 2.0 * v3) / 6.0;
    const double p2 = (2.0 * v2 + 5.0 * v3 - v4) / 6.0;
    const double a = v0 - 2.0 * v1 + v2;
    const double b = v0 - 4.0 * v1 + 3.0 * v2;
    const double c = v1 - 2.0 * v2 + v3;
    const double d = v1 - v3;
    const double e = v2 - 2.0 * v3 + v4;
    const double f = 3.0 * v2 - 4.0 * v3 + v4;
    const double b0 = (13.0 / 12.0) * a * a + 0.25 * b * b;
    const double b1 = (13.0 / 12.0) * c * c + 0.25 * d * d;
    const double b2 = (13.0 / 12.0) * e * e + 0.25 * f * f;
    double a0, a1, a2;
    if (use_z) {
        constexpr double eps = 1.0e-40;
        const double tau = fabs(b0 - b2);
        const double r0 = tau / (eps + b0);
        const double r1 = tau / (eps + b1);
        const double r2 = tau / (eps + b2);
        a0 = 0.1 * (1.0 + r0 * r0);
        a1 = 0.6 * (1.0 + r1 * r1);
        a2 = 0.3 * (1.0 + r2 * r2);
    } else {
        constexpr double eps = 1.0e-6;
        a0 = 0.1 / ((eps + b0) * (eps + b0));
        a1 = 0.6 / ((eps + b1) * (eps + b1));
        a2 = 0.3 / ((eps + b2) * (eps + b2));
    }
    const double sum = a0 + a1 + a2;
    if (!(sum > 0.0) || !isfinite(sum))
        return v2;
    const double out = (a0 * p0 + a1 * p1 + a2 * p2) / sum;
    return isfinite(out) ? out : v2;
}

__device__ __forceinline__ int line_cell(int cell, int face, int offset,
                                         int nx, int ny, int nz)
{
    const int plane = nx * ny;
    int k = cell / plane;
    const int rem = cell - k * plane;
    int j = rem / nx;
    int i = rem - j * nx;
    const int step = (face & 1) ? offset : -offset;
    if (face < 2)
        i += step;
    else if (face < 4)
        j += step;
    else
        k += step;
    if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz)
        return -1;
    return i + j * nx + k * plane;
}

__device__ __forceinline__ void build_char_matrices(const double qL[6],
                                                     const double qR[6],
                                                     double nx, double ny, double nz,
                                                     double L[25], double R[25])
{
    const double rL = sqrt(fmax(qL[0], kTiny));
    const double rR = sqrt(fmax(qR[0], kTiny));
    const double den = fmax(rL + rR, kTiny);
    const double u = (rL * qL[1] + rR * qR[1]) / den;
    const double v = (rL * qL[2] + rR * qR[2]) / den;
    const double w = (rL * qL[3] + rR * qR[3]) / den;
    const double sound = fmax((rL * qL[5] + rR * qR[5]) / den, kTiny);
    const double un = u * nx + v * ny + w * nz;
    const double ek = 0.5 * (u * u + v * v + w * w);
    const double H = sound * sound / (kGamma - 1.0) + ek;
    const double t1 = 0.5 / (sound * sound);
    const double t2 = (kGamma - 1.0) * t1;

    double tx, ty, tz;
    if (fabs(nx) > fabs(nz)) {
        const double inv = rsqrt(nx * nx + ny * ny + 1.0e-30);
        tx = -ny * inv;
        ty = nx * inv;
        tz = 0.0;
    } else {
        const double inv = rsqrt(ny * ny + nz * nz + 1.0e-30);
        tx = 0.0;
        ty = -nz * inv;
        tz = ny * inv;
    }
    const double sx = ny * tz - nz * ty;
    const double sy = nz * tx - nx * tz;
    const double sz = nx * ty - ny * tx;
    const double ut = u * tx + v * ty + w * tz;
    const double us = u * sx + v * sy + w * sz;
    const double cols[5][5] = {
        {1.0, u - sound * nx, v - sound * ny, w - sound * nz, H - sound * un},
        {1.0, u, v, w, ek},
        {0.0, tx, ty, tz, ut},
        {0.0, sx, sy, sz, us},
        {1.0, u + sound * nx, v + sound * ny, w + sound * nz, H + sound * un}};
#pragma unroll
    for (int col = 0; col < 5; ++col)
#pragma unroll
        for (int row = 0; row < 5; ++row)
            R[row * 5 + col] = cols[col][row];

    L[0] = t2 * ek + t1 * sound * un;
    L[1] = -t2 * u - t1 * sound * nx;
    L[2] = -t2 * v - t1 * sound * ny;
    L[3] = -t2 * w - t1 * sound * nz;
    L[4] = t2;
    L[5] = 1.0 - 2.0 * t2 * ek;
    L[6] = 2.0 * t2 * u;
    L[7] = 2.0 * t2 * v;
    L[8] = 2.0 * t2 * w;
    L[9] = -2.0 * t2;
    L[10] = -ut;
    L[11] = tx;
    L[12] = ty;
    L[13] = tz;
    L[14] = 0.0;
    L[15] = -us;
    L[16] = sx;
    L[17] = sy;
    L[18] = sz;
    L[19] = 0.0;
    L[20] = t2 * ek - t1 * sound * un;
    L[21] = -t2 * u + t1 * sound * nx;
    L[22] = -t2 * v + t1 * sound * ny;
    L[23] = -t2 * w + t1 * sound * nz;
    L[24] = t2;
}

__device__ __forceinline__ void reconstruct_weno(
    const double *__restrict__ U, int cell, int face, int nx, int ny, int nz,
    double fnx, double fny, double fnz, bool use_char, bool use_z, int weno_hybrid,
    double UL[kNv], double UR[kNv])
{
    int ids[6];
#pragma unroll
    for (int s = 0; s < 6; ++s) {
        ids[s] = line_cell(cell, face, s - 2, nx, ny, nz);
        if (ids[s] < 0)
            return;
    }
    double sten[6][kNv];
#pragma unroll
    for (int s = 0; s < 6; ++s)
#pragma unroll
        for (int v = 0; v < kNv; ++v)
            sten[s][v] = U[ids[s] * kNv + v];

    if (!use_char) {
#pragma unroll
        for (int v = 0; v < kNv; ++v) {
            UL[v] = weno5(sten[0][v], sten[1][v], sten[2][v],
                          sten[3][v], sten[4][v], use_z);
            UR[v] = weno5(sten[5][v], sten[4][v], sten[3][v],
                          sten[2][v], sten[1][v], use_z);
        }
    } else {
        double qL[6], qR[6], L[25], R[25], chars[6][kNv];
        fill_prim(sten[2], qL);
        fill_prim(sten[3], qR);
        build_char_matrices(qL, qR, fnx, fny, fnz, L, R);
#pragma unroll
        for (int s = 0; s < 6; ++s)
#pragma unroll
            for (int row = 0; row < kNv; ++row) {
                double sum = 0.0;
#pragma unroll
                for (int col = 0; col < kNv; ++col)
                    sum += L[row * kNv + col] * sten[s][col];
                chars[s][row] = sum;
            }
        double wl[kNv], wr[kNv];
#pragma unroll
        for (int v = 0; v < kNv; ++v) {
            wl[v] = weno5(chars[0][v], chars[1][v], chars[2][v],
                          chars[3][v], chars[4][v], use_z);
            wr[v] = weno5(chars[5][v], chars[4][v], chars[3][v],
                          chars[2][v], chars[1][v], use_z);
        }
#pragma unroll
        for (int row = 0; row < kNv; ++row) {
            double sl = 0.0, sr = 0.0;
#pragma unroll
            for (int col = 0; col < kNv; ++col) {
                sl += R[row * kNv + col] * wl[col];
                sr += R[row * kNv + col] * wr[col];
            }
            UL[row] = sl;
            UR[row] = sr;
        }
    }
    if (!admissible(UL) || !admissible(UR)) {
#pragma unroll
        for (int v = 0; v < kNv; ++v) {
            UL[v] = sten[2][v];
            UR[v] = sten[3][v];
        }
        return;
    }

    if (weno_hybrid > 0) {
        double qL[6], qR[6];
        fill_prim(UL, qL);
        fill_prim(UR, qR);
        const double unL = qL[1] * fnx + qL[2] * fny + qL[3] * fnz;
        const double unR = qR[1] * fnx + qR[2] * fny + qR[3] * fnz;
        const double a_ref = 0.5 * (qL[5] + qR[5]) + 1.0;
        const double phi_P =
            fabs(qR[4] - qL[4]) / (0.5 * (fabs(qL[4]) + fabs(qR[4])) + 1.0);
        const double phi_rho =
            fabs(qR[0] - qL[0]) / (0.5 * (fabs(qL[0]) + fabs(qR[0])) + 1.0);
        const double sL = qL[4] / pow(fmax(qL[0], kTiny), kGamma);
        const double sR = qR[4] / pow(fmax(qR[0], kTiny), kGamma);
        const double phi_S =
            fabs(sR - sL) / (0.5 * (fabs(sL) + fabs(sR)) + 1.0);
        const double phi_u = fabs(unR - unL) / a_ref;
        double phi = phi_P;
        if (weno_hybrid == 2) {
            phi = fmax(phi_P, phi_u);
            if (phi_P < 0.05)
                phi = fmax(phi, fmax(phi_rho, phi_S));
        }
#pragma unroll
        for (int v = 0; v < kNv; ++v) {
            const double C = (v >= 1 && v <= 3) ? 120.0 : 50.0;
            const double theta = 1.0 / (1.0 + C * phi);
            UL[v] = theta * UL[v] + (1.0 - theta) * sten[2][v];
            UR[v] = theta * UR[v] + (1.0 - theta) * sten[3][v];
            const double lo = fmin(fmin(sten[1][v], sten[2][v]),
                                   fmin(sten[3][v], sten[4][v]));
            const double hi = fmax(fmax(sten[1][v], sten[2][v]),
                                   fmax(sten[3][v], sten[4][v]));
            UL[v] = fmin(hi, fmax(lo, UL[v]));
            UR[v] = fmin(hi, fmax(lo, UR[v]));
        }
        if (!admissible(UL) || !admissible(UR)) {
#pragma unroll
            for (int v = 0; v < kNv; ++v) {
                UL[v] = sten[2][v];
                UR[v] = sten[3][v];
            }
        }
    }
}

__device__ __forceinline__ void flux_llf(const double UL[kNv], const double UR[kNv],
                                         double nx, double ny, double nz, double area,
                                         double flux[kNv])
{
    double qL[6], qR[6], FL[kNv], FR[kNv];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    physical_flux(UL, qL, nx, ny, nz, FL);
    physical_flux(UR, qR, nx, ny, nz, FR);
    const double unL = qL[1] * nx + qL[2] * ny + qL[3] * nz;
    const double unR = qR[1] * nx + qR[2] * ny + qR[3] * nz;
    const double speed = fmax(fabs(unL) + qL[5], fabs(unR) + qR[5]);
#pragma unroll
    for (int v = 0; v < kNv; ++v)
        flux[v] = 0.5 * (FL[v] + FR[v] - speed * (UR[v] - UL[v])) * area;
}

__device__ __forceinline__ double ricca_phi_rh(const double qL[6], const double qR[6],
                                               double unL, double unR)
{
    const double massL = qL[0] * unL;
    const double massR = qR[0] * unR;
    const double FnL = massL * unL + qL[4];
    const double FnR = massR * unR + qR[4];
    const double dr = qR[0] - qL[0];
    const double s = (fabs(dr) > kTiny) ? (massR - massL) / dr : 0.5 * (unL + unR);
    return fabs(FnR - FnL - s * (massR - massL)) /
           (0.5 * (fabs(FnL) + fabs(FnR)) + 1.0);
}

__device__ __forceinline__ void flux_ricca(const double UL[kNv], const double UR[kNv],
                                           double nx, double ny, double nz, double area,
                                           int sensor_type, double rh_threshold,
                                           double flux[kNv])
{
    double qL[6], qR[6], FL[kNv], FR[kNv];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    physical_flux(UL, qL, nx, ny, nz, FL);
    physical_flux(UR, qR, nx, ny, nz, FR);
    const double unL = qL[1] * nx + qL[2] * ny + qL[3] * nz;
    const double unR = qR[1] * nx + qR[2] * ny + qR[3] * nz;
    const double a_contact = fmax(fabs(unL), fabs(unR));
    const double a_acoustic = fmax(fabs(unL) + qL[5], fabs(unR) + qR[5]);
    const double phi_P =
        fabs(qR[4] - qL[4]) / (0.5 * (fabs(qL[4]) + fabs(qR[4])) + 1.0);
    bool use_acoustic = false;
    if (sensor_type == 1) {
        const double phi_RH = ricca_phi_rh(qL, qR, unL, unR);
        use_acoustic = (phi_P > 1.0e-8) || (phi_RH > rh_threshold);
    } else {
        use_acoustic = (phi_P > 1.0e-8);
    }
    const double alpha = use_acoustic ? a_acoustic : a_contact;
#pragma unroll
    for (int v = 0; v < kNv; ++v)
        flux[v] = 0.5 * (FL[v] + FR[v]) * area - 0.5 * alpha * area * (UR[v] - UL[v]);
}

__device__ __forceinline__ double entropy_abs(double x, double delta)
{
    const double ax = fabs(x);
    return ax < delta ? 0.5 * (ax * ax / delta + delta) : ax;
}

__device__ __forceinline__ void flux_roe(const double UL[kNv], const double UR[kNv],
                                         double nx, double ny, double nz, double area,
                                         double flux[kNv])
{
    double qL[6], qR[6], FL[kNv], FR[kNv];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    physical_flux(UL, qL, nx, ny, nz, FL);
    physical_flux(UR, qR, nx, ny, nz, FR);
    double tx, ty, tz;
    if (fabs(nx) > fabs(nz)) {
        const double inv = rsqrt(nx * nx + ny * ny + 1.0e-30);
        tx = -ny * inv; ty = nx * inv; tz = 0.0;
    } else {
        const double inv = rsqrt(ny * ny + nz * nz + 1.0e-30);
        tx = 0.0; ty = -nz * inv; tz = ny * inv;
    }
    const double sx = ny * tz - nz * ty;
    const double sy = nz * tx - nx * tz;
    const double sz = nx * ty - ny * tx;
    const double rl = sqrt(qL[0]), rr = sqrt(qR[0]), den = rl + rr;
    const double u = (rl * qL[1] + rr * qR[1]) / den;
    const double v = (rl * qL[2] + rr * qR[2]) / den;
    const double w = (rl * qL[3] + rr * qR[3]) / den;
    const double HL = (UL[4] + qL[4]) / qL[0];
    const double HR = (UR[4] + qR[4]) / qR[0];
    const double H = (rl * HL + rr * HR) / den;
    const double vel2 = u * u + v * v + w * w;
    const double sound = sqrt(fmax((kGamma - 1.0) * (H - 0.5 * vel2), kTiny));
    const double rho = rl * rr;
    const double un = u * nx + v * ny + w * nz;
    const double ut = u * tx + v * ty + w * tz;
    const double us = u * sx + v * sy + w * sz;
    const double dr = qR[0] - qL[0];
    const double dun = (qR[1] - qL[1]) * nx + (qR[2] - qL[2]) * ny +
                       (qR[3] - qL[3]) * nz;
    const double dut = (qR[1] - qL[1]) * tx + (qR[2] - qL[2]) * ty +
                       (qR[3] - qL[3]) * tz;
    const double dus = (qR[1] - qL[1]) * sx + (qR[2] - qL[2]) * sy +
                       (qR[3] - qL[3]) * sz;
    const double dp = qR[4] - qL[4], a2 = sound * sound;
    const double alpha[kNv] = {(dp - rho * sound * dun) / (2.0 * a2),
                               dr - dp / a2, rho * dut, rho * dus,
                               (dp + rho * sound * dun) / (2.0 * a2)};
    const double delta = fmax(0.1 * sound, kTiny);
    const double lambda[kNv] = {entropy_abs(un - sound, delta),
                                 entropy_abs(un, delta), entropy_abs(un, delta),
                                 entropy_abs(un, delta), entropy_abs(un + sound, delta)};
    const double wave[kNv][kNv] = {
        {1.0, u - sound * nx, v - sound * ny, w - sound * nz, H - sound * un},
        {1.0, u, v, w, 0.5 * vel2},
        {0.0, tx, ty, tz, ut},
        {0.0, sx, sy, sz, us},
        {1.0, u + sound * nx, v + sound * ny, w + sound * nz, H + sound * un}};
#pragma unroll
    for (int n = 0; n < kNv; ++n) {
        double diss = 0.0;
#pragma unroll
        for (int m = 0; m < kNv; ++m)
            diss += lambda[m] * alpha[m] * wave[m][n];
        flux[n] = (0.5 * (FL[n] + FR[n]) - 0.5 * diss) * area;
    }
}

__device__ __forceinline__ void vanleer_split(const double U[kNv], const double q[6],
                                               double nx, double ny, double nz,
                                               bool plus, double F[kNv])
{
    const double un = q[1] * nx + q[2] * ny + q[3] * nz;
    const double M = un / q[5];
    if ((plus && M >= 1.0) || (!plus && M <= -1.0)) {
        physical_flux(U, q, nx, ny, nz, F);
        return;
    }
    if ((plus && M <= -1.0) || (!plus && M >= 1.0)) {
#pragma unroll
        for (int v = 0; v < kNv; ++v)
            F[v] = 0.0;
        return;
    }
    const double sign = plus ? 1.0 : -1.0;
    const double ms = M + sign;
    const double mass = sign * 0.25 * q[0] * q[5] * ms * ms;
    const double ps = 0.25 * q[4] * ms * ms * (2.0 - sign * M);
    const double vt2 = q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - un * un;
    const double hb = (kGamma - 1.0) * un + sign * 2.0 * q[5];
    const double hs = hb * hb / (2.0 * (kGamma * kGamma - 1.0)) + 0.5 * vt2;
    F[0] = mass;
    F[1] = mass * q[1] + ps * nx;
    F[2] = mass * q[2] + ps * ny;
    F[3] = mass * q[3] + ps * nz;
    F[4] = mass * hs;
}

__device__ __forceinline__ void flux_vanleer(const double UL[kNv], const double UR[kNv],
                                              double nx, double ny, double nz, double area,
                                              double flux[kNv])
{
    double qL[6], qR[6], plus[kNv], minus[kNv];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    vanleer_split(UL, qL, nx, ny, nz, true, plus);
    vanleer_split(UR, qR, nx, ny, nz, false, minus);
#pragma unroll
    for (int v = 0; v < kNv; ++v)
        flux[v] = (plus[v] + minus[v]) * area;
}

__global__ void net_flux_kernel(const double *__restrict__ U, double *__restrict__ residual,
                                const int *__restrict__ neigh,
                                const double *__restrict__ normals,
                                const double *__restrict__ areas,
                                const int *__restrict__ bc, const double *__restrict__ fs,
                                const double *__restrict__ ps,
                                int n_cells, int nx, int ny, int nz,
                                int scheme, int order, bool weno_char, bool weno_z,
                                int ricca_sensor, double ricca_rh_thresh, int weno_hybrid)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;
    double center[kNv], sum[kNv] = {0.0, 0.0, 0.0, 0.0, 0.0};
#pragma unroll
    for (int v = 0; v < kNv; ++v)
        center[v] = U[cell * kNv + v];

#pragma unroll
    for (int face = 0; face < kFaces; ++face) {
        const int fi = cell * kFaces + face;
        const double fnx = normals[fi * 3];
        const double fny = normals[fi * 3 + 1];
        const double fnz = normals[fi * 3 + 2];
        const int nb = neigh[fi];
        double UL[kNv], UR[kNv], flux[kNv];
#pragma unroll
        for (int v = 0; v < kNv; ++v)
            UL[v] = center[v];
        if (nb >= 0) {
#pragma unroll
            for (int v = 0; v < kNv; ++v)
                UR[v] = U[nb * kNv + v];
            if (order == 2)
                reconstruct_weno(U, cell, face, nx, ny, nz, fnx, fny, fnz,
                                 weno_char, weno_z, weno_hybrid, UL, UR);
        } else {
            boundary_ghost(UL, bc[fi], fnx, fny, fnz, fs, ps, UR);
        }
        if (scheme == 1)
            flux_ricca(UL, UR, fnx, fny, fnz, areas[fi], ricca_sensor, ricca_rh_thresh,
                       flux);
        else if (scheme == 2)
            flux_roe(UL, UR, fnx, fny, fnz, areas[fi], flux);
        else if (scheme == 3)
            flux_vanleer(UL, UR, fnx, fny, fnz, areas[fi], flux);
        else
            flux_llf(UL, UR, fnx, fny, fnz, areas[fi], flux);
#pragma unroll
        for (int v = 0; v < kNv; ++v)
            sum[v] += flux[v];
    }
#pragma unroll
    for (int v = 0; v < kNv; ++v)
        residual[cell * kNv + v] = sum[v];
}

__global__ void update_kernel(double *__restrict__ U, const double *__restrict__ residual,
                              const double *__restrict__ volume, double dt, int n_cells)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;
    double old[kNv], next[kNv];
    const double scale = dt / volume[cell];
#pragma unroll
    for (int v = 0; v < kNv; ++v) {
        old[v] = U[cell * kNv + v];
        next[v] = old[v] - scale * residual[cell * kNv + v];
    }
    const bool keep = admissible(next);
#pragma unroll
    for (int v = 0; v < kNv; ++v)
        U[cell * kNv + v] = keep ? next[v] : old[v];
}

__global__ void local_dt_kernel(const double *__restrict__ U,
                                const double *__restrict__ volume,
                                const double *__restrict__ normals,
                                const double *__restrict__ areas,
                                double *__restrict__ local_dt, double cfl, int n_cells)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;
    double q[6];
    fill_prim(U + cell * kNv, q);
    double denom = 0.0;
#pragma unroll
    for (int face = 0; face < kFaces; ++face) {
        const int fi = cell * kFaces + face;
        const double un = q[1] * normals[fi * 3] + q[2] * normals[fi * 3 + 1] +
                          q[3] * normals[fi * 3 + 2];
        denom += (fabs(un) + q[5]) * areas[fi];
    }
    local_dt[cell] = denom > 1.0e-30 ? cfl * volume[cell] / denom : 1.0e300;
}

constexpr int kNGrad = 4;
constexpr int kNDim = 3;
constexpr int kBCNoslip = 2;

__device__ __forceinline__ double primitive_value(const double *U, int cell, int variable)
{
    double q[6];
    fill_prim(U + cell * kNv, q);
    if (variable < 3)
        return q[variable + 1];
    return q[4] / q[0];
}

__device__ __forceinline__ double neighbour_distance(const double *xc, const double *yc,
                                                     const double *zc, int cell, int neighbour,
                                                     int face, double fallback)
{
    if (neighbour < 0)
        return fallback;
    if (face < 2)
        return fabs(xc[neighbour] - xc[cell]);
    if (face < 4)
        return fabs(yc[neighbour] - yc[cell]);
    return fabs(zc[neighbour] - zc[cell]);
}

__device__ __forceinline__ double directional_gradient(const double *U, const int *neigh,
                                                       const int *bc, const double *xc,
                                                       const double *yc, const double *zc,
                                                       int cell, int variable, int neg_face,
                                                       int pos_face, double spacing)
{
    const int negative = neigh[cell * kFaces + neg_face];
    const int positive = neigh[cell * kFaces + pos_face];

    if (negative >= 0 && positive >= 0) {
        const double dn = neighbour_distance(xc, yc, zc, cell, negative, neg_face, spacing);
        const double dp = neighbour_distance(xc, yc, zc, cell, positive, pos_face, spacing);
        return (primitive_value(U, positive, variable) - primitive_value(U, negative, variable)) /
               fmax(dn + dp, 1.0e-14);
    }

    const double center = primitive_value(U, cell, variable);
    if (negative < 0) {
        if (bc[cell * kFaces + neg_face] == kBCNoslip) {
            const double boundary = variable < 3 ? 0.0 : center;
            const double h = (positive >= 0)
                                 ? neighbour_distance(xc, yc, zc, cell, positive, pos_face, spacing)
                                 : spacing;
            return (center - boundary) / fmax(0.5 * h, 1.0e-14);
        }
        if (positive >= 0)
            return (primitive_value(U, positive, variable) - center) /
                   fmax(neighbour_distance(xc, yc, zc, cell, positive, pos_face, spacing), 1.0e-14);
        return 0.0;
    }

    if (bc[cell * kFaces + pos_face] == kBCNoslip) {
        const double boundary = variable < 3 ? 0.0 : center;
        const double h = neighbour_distance(xc, yc, zc, cell, negative, neg_face, spacing);
        return (boundary - center) / fmax(0.5 * h, 1.0e-14);
    }
    return (center - primitive_value(U, negative, variable)) /
           fmax(neighbour_distance(xc, yc, zc, cell, negative, neg_face, spacing), 1.0e-14);
}

__global__ void viscous_grad_kernel(const double *__restrict__ U, const int *__restrict__ neigh,
                                    const int *__restrict__ bc, const double *__restrict__ xc,
                                    const double *__restrict__ yc, const double *__restrict__ zc,
                                    double *__restrict__ grad, double dx, double dy, double dz,
                                    int n_cells)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;
    for (int variable = 0; variable < kNGrad; ++variable) {
        grad[cell * 12 + variable * kNDim] =
            directional_gradient(U, neigh, bc, xc, yc, zc, cell, variable, 0, 1, dx);
        grad[cell * 12 + variable * kNDim + 1] =
            directional_gradient(U, neigh, bc, xc, yc, zc, cell, variable, 2, 3, dy);
        grad[cell * 12 + variable * kNDim + 2] =
            directional_gradient(U, neigh, bc, xc, yc, zc, cell, variable, 4, 5, dz);
    }
}

__device__ __forceinline__ double d_strain_mag(const double g[12])
{
    const double sxx = g[0], syy = g[4], szz = g[8];
    const double sxy = 0.5 * (g[1] + g[3]);
    const double sxz = 0.5 * (g[2] + g[6]);
    const double syz = 0.5 * (g[5] + g[7]);
    return sqrt(fmax(2.0 * (sxx * sxx + syy * syy + szz * szz) +
                         4.0 * (sxy * sxy + sxz * sxz + syz * syz),
                     0.0));
}

__device__ __forceinline__ double d_wale_sd(const double g[12])
{
    const double G[3][3] = {{g[0], g[1], g[2]}, {g[3], g[4], g[5]}, {g[6], g[7], g[8]}};
    double G2[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                G2[i][j] += G[i][k] * G[k][j];
    const double tr = G2[0][0] + G2[1][1] + G2[2][2];
    double Sd2 = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double sd = 0.5 * (G2[i][j] + G2[j][i]);
            if (i == j)
                sd -= tr / 3.0;
            Sd2 += sd * sd;
        }
    return sqrt(fmax(Sd2, 0.0));
}

__device__ __forceinline__ double d_vreman(const double g[12], double delta)
{
    const double a[3][3] = {{g[0], g[3], g[6]}, {g[1], g[4], g[7]}, {g[2], g[5], g[8]}};
    double alpha2 = 0.0;
    double beta[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int m = 0; m < 3; ++m) {
                s += a[m][i] * a[m][j];
                alpha2 += a[i][j] * a[i][j];
            }
            beta[i][j] = delta * delta * s;
        }
    const double B = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
                     beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
                     beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];
    if (alpha2 < 1.0e-30 || B < 0.0)
        return 0.0;
    return sqrt(B / alpha2);
}

__global__ void sgs_viscosity_kernel(const double *__restrict__ U, const double *__restrict__ grad,
                                     const double *__restrict__ volume,
                                     const double *__restrict__ wall_dist,
                                     double *__restrict__ mu_sgs, int turb, double mu_lam,
                                     double Cs, double Cw, double Cv, double kappa, double Aplus,
                                     double Uinf, int n_cells)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;
    if (turb <= 1) {
        mu_sgs[cell] = 0.0;
        return;
    }
    double q[6];
    fill_prim(U + cell * kNv, q);
    const double rho = fmax(q[0], 1.0e-14);
    const double nu = mu_lam / rho;
    const double g[12] = {grad[cell * 12 + 0],  grad[cell * 12 + 1],  grad[cell * 12 + 2],
                          grad[cell * 12 + 3],  grad[cell * 12 + 4],  grad[cell * 12 + 5],
                          grad[cell * 12 + 6],  grad[cell * 12 + 7],  grad[cell * 12 + 8],
                          grad[cell * 12 + 9],  grad[cell * 12 + 10], grad[cell * 12 + 11]};
    const double S = d_strain_mag(g);
    const double delta = pow(fmax(volume[cell], 1.0e-30), 1.0 / 3.0);
    const double yw = fmax(wall_dist[cell], 0.0);
    double mut = 0.0;
    if (turb == 2) {
        const double ut = 0.04 * fmax(fabs(Uinf), 1.0e-6);
        const double yp = yw * ut / fmax(nu, 1.0e-14);
        const double mix = kappa * yw * (1.0 - exp(-yp / fmax(Aplus, 1.0)));
        mut = rho * mix * mix * S;
    } else if (turb == 3) {
        const double ut = 0.04 * fmax(fabs(Uinf), 1.0e-6);
        const double yp = yw * ut / fmax(nu, 1.0e-14);
        const double fvd = 1.0 - exp(-yp / fmax(Aplus, 1.0));
        const double ls = Cs * delta * fvd;
        mut = rho * ls * ls * S;
    } else if (turb == 4) {
        const double Sd = d_wale_sd(g);
        const double den = pow(S, 2.5) + pow(Sd, 1.25);
        const double cw = Cw * delta;
        mut = den > 1.0e-30 ? rho * cw * cw * pow(Sd, 1.5) / den : 0.0;
    } else if (turb == 5) {
        mut = rho * Cv * d_vreman(g, delta);
    }
    mu_sgs[cell] = mut;
}

__global__ void viscous_flux_kernel(const double *__restrict__ U, const double *__restrict__ grad,
                                    const int *__restrict__ neigh, const double *__restrict__ normals,
                                    const double *__restrict__ areas,
                                    const double *__restrict__ mu_sgs,
                                    double *__restrict__ residual, double mu_lam, double cp,
                                    double Pr, double Pr_t, int n_cells)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;

    for (int face = 0; face < kFaces; ++face) {
        const int fi = cell * kFaces + face;
        const int neighbour = neigh[fi];
        const bool wall = neighbour < 0;

        double fg[12];
#pragma unroll
        for (int component = 0; component < 12; ++component) {
            fg[component] = grad[cell * 12 + component];
            if (!wall)
                fg[component] = 0.5 * (fg[component] + grad[neighbour * 12 + component]);
        }

        const double ux = fg[0], uy = fg[1], uz = fg[2];
        const double vx = fg[3], vy = fg[4], vz = fg[5];
        const double wx = fg[6], wy = fg[7], wz = fg[8];
            const double mu_t_f =
                wall ? mu_sgs[cell] : 0.5 * (mu_sgs[cell] + mu_sgs[neighbour]);
            const double mu = mu_lam + mu_t_f;
            const double kappa =
                (Pr > 0.0 ? mu_lam * cp / Pr : 0.0) + (Pr_t > 0.0 ? mu_t_f * cp / Pr_t : 0.0);
            const double two_thirds_mu = (2.0 / 3.0) * mu;
        const double tau_xx = two_thirds_mu * (2.0 * ux - vy - wz);
        const double tau_yy = two_thirds_mu * (2.0 * vy - ux - wz);
        const double tau_zz = two_thirds_mu * (2.0 * wz - ux - vy);
        const double tau_xy = mu * (uy + vx);
        const double tau_xz = mu * (uz + wx);
        const double tau_yz = mu * (vz + wy);
        const double qx = -kappa * fg[9];
        const double qy = -kappa * fg[10];
        const double qz = -kappa * fg[11];

        double q[6];
        fill_prim(U + cell * kNv, q);
        double u = q[1], v = q[2], w = q[3];
        if (wall) {
            u = v = w = 0.0;
        } else {
            double qn[6];
            fill_prim(U + neighbour * kNv, qn);
            u = 0.5 * (u + qn[1]);
            v = 0.5 * (v + qn[2]);
            w = 0.5 * (w + qn[3]);
        }

        const double nx = normals[fi * 3];
        const double ny = normals[fi * 3 + 1];
        const double nz = normals[fi * 3 + 2];
        const double area = areas[fi];
        const double traction_x = tau_xx * nx + tau_xy * ny + tau_xz * nz;
        const double traction_y = tau_xy * nx + tau_yy * ny + tau_yz * nz;
        const double traction_z = tau_xz * nx + tau_yz * ny + tau_zz * nz;
        const double energy_flux =
            u * traction_x + v * traction_y + w * traction_z - (qx * nx + qy * ny + qz * nz);

        residual[cell * kNv + 1] -= traction_x * area;
        residual[cell * kNv + 2] -= traction_y * area;
        residual[cell * kNv + 3] -= traction_z * area;
        residual[cell * kNv + 4] -= energy_flux * area;
    }
}

__global__ void viscous_dt_kernel(const double *__restrict__ U, const double *__restrict__ volume,
                                  const double *__restrict__ normals,
                                  const double *__restrict__ areas, const double *__restrict__ xc,
                                  const double *__restrict__ yc, const double *__restrict__ zc,
                                  double *__restrict__ local_dt, double cfl, double mu, double dx,
                                  double dy, double dz, int nx, int ny, int nz, int n_cells)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= n_cells)
        return;
    double q[6];
    fill_prim(U + cell * kNv, q);
    double denom = 0.0;
#pragma unroll
    for (int face = 0; face < kFaces; ++face) {
        const int fi = cell * kFaces + face;
        const double un = q[1] * normals[fi * 3] + q[2] * normals[fi * 3 + 1] +
                          q[3] * normals[fi * 3 + 2];
        denom += (fabs(un) + q[5]) * areas[fi];
    }
    double dt = denom > 1.0e-30 ? cfl * volume[cell] / denom : 1.0e300;
    if (mu > 0.0) {
        const int i = cell % nx;
        const int k = cell / (nx * ny);
        const int j = (cell / nx) % ny;
        double hx = dx, hy = dy, hz = dz;
        if (i + 1 < nx)
            hx = fabs(xc[cell + 1] - xc[cell]);
        else if (i > 0)
            hx = fabs(xc[cell] - xc[cell - 1]);
        if (j + 1 < ny)
            hy = fabs(yc[cell + nx] - yc[cell]);
        else if (j > 0)
            hy = fabs(yc[cell] - yc[cell - nx]);
        if (k + 1 < nz)
            hz = fabs(zc[cell + nx * ny] - zc[cell]);
        else if (k > 0)
            hz = fabs(zc[cell] - zc[cell - nx * ny]);
        const double h2 = fmin(hx * hx, fmin(hy * hy, hz * hz));
        const double nu = mu / fmax(q[0], 1.0e-14);
        dt = fmin(dt, 0.25 * cfl * h2 / fmax(nu, 1.0e-30));
    }
    local_dt[cell] = dt;
}

} // namespace

namespace euler3d_cuda {

void release(DeviceBundle &dev)
{
#define FREE_MEMBER(member)             \
    do {                                \
        if (dev.member != nullptr) {    \
            CUDA_CHECK(cudaFree(dev.member)); \
            dev.member = nullptr;       \
        }                               \
    } while (0)
    FREE_MEMBER(d_U);
    FREE_MEMBER(d_R);
    FREE_MEMBER(d_volume);
    FREE_MEMBER(d_nxyz);
    FREE_MEMBER(d_area);
    FREE_MEMBER(d_neigh);
    FREE_MEMBER(d_bc);
    FREE_MEMBER(d_freestream);
    FREE_MEMBER(d_postshock);
    FREE_MEMBER(d_xc);
    FREE_MEMBER(d_yc);
    FREE_MEMBER(d_zc);
    FREE_MEMBER(d_grad);
    FREE_MEMBER(d_wall_dist);
    FREE_MEMBER(d_mu_sgs);
#undef FREE_MEMBER
    dev.n_cells = dev.nx = dev.ny = dev.nz = 0;
    dev.dx = dev.dy = dev.dz = 0.0;
}

void allocate(DeviceBundle &dev, int n_cells)
{
    if (n_cells <= 0) {
        std::fprintf(stderr, "Euler3D CUDA: n_cells must be positive\n");
        std::abort();
    }
    release(dev);
    dev.n_cells = n_cells;
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_U),
                          sizeof(double) * n_cells * kNv));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_R),
                          sizeof(double) * n_cells * kNv));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_volume),
                          sizeof(double) * n_cells));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_nxyz),
                          sizeof(double) * n_cells * kFaces * 3));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_area),
                          sizeof(double) * n_cells * kFaces));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_neigh),
                          sizeof(int) * n_cells * kFaces));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_bc),
                          sizeof(int) * n_cells * kFaces));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_freestream),
                          sizeof(double) * kNv));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_postshock),
                          sizeof(double) * kNv));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_xc), sizeof(double) * n_cells));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_yc), sizeof(double) * n_cells));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_zc), sizeof(double) * n_cells));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_grad),
                          sizeof(double) * n_cells * 12));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_wall_dist),
                          sizeof(double) * n_cells));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&dev.d_mu_sgs),
                          sizeof(double) * n_cells));
}

void set_grid_size(DeviceBundle &dev, int nx, int ny, int nz)
{
    if (nx <= 0 || ny <= 0 || nz <= 0 ||
        static_cast<long long>(nx) * ny * nz != dev.n_cells) {
        std::fprintf(stderr, "Euler3D CUDA: grid %dx%dx%d does not match %d cells\n",
                     nx, ny, nz, dev.n_cells);
        std::abort();
    }
    dev.nx = nx;
    dev.ny = ny;
    dev.nz = nz;
}

void upload_mesh(DeviceBundle &dev, const std::vector<double> &volume,
                 const std::vector<int> &neigh, const std::vector<double> &nxyz,
                 const std::vector<double> &area, const std::vector<int> &bc_type,
                 const std::vector<double> &xc, const std::vector<double> &yc,
                 const std::vector<double> &zc, const std::vector<double> &wall_dist,
                 double dx, double dy, double dz)
{
    const std::size_t cells = static_cast<std::size_t>(dev.n_cells);
    require_size(volume.size(), cells, "volume");
    require_size(neigh.size(), cells * kFaces, "neigh");
    require_size(nxyz.size(), cells * kFaces * 3, "nxyz");
    require_size(area.size(), cells * kFaces, "area");
    require_size(bc_type.size(), cells * kFaces, "bc_type");
    require_size(xc.size(), cells, "xc");
    require_size(yc.size(), cells, "yc");
    require_size(zc.size(), cells, "zc");
    require_size(wall_dist.size(), cells, "wall_dist");
    CUDA_CHECK(cudaMemcpy(dev.d_volume, volume.data(), sizeof(double) * cells,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_neigh, neigh.data(), sizeof(int) * cells * kFaces,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_nxyz, nxyz.data(), sizeof(double) * cells * kFaces * 3,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_area, area.data(), sizeof(double) * cells * kFaces,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_bc, bc_type.data(), sizeof(int) * cells * kFaces,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_xc, xc.data(), sizeof(double) * cells, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_yc, yc.data(), sizeof(double) * cells, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_zc, zc.data(), sizeof(double) * cells, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev.d_wall_dist, wall_dist.data(), sizeof(double) * cells,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(dev.d_mu_sgs, 0, sizeof(double) * cells));
    dev.dx = dx;
    dev.dy = dy;
    dev.dz = dz;
}

void set_freestream(DeviceBundle &dev, double rho, double u, double v, double w, double p)
{
    const double fs[kNv] = {rho, u, v, w, p};
    CUDA_CHECK(cudaMemcpy(dev.d_freestream, fs, sizeof(fs), cudaMemcpyHostToDevice));
}

void set_postshock(DeviceBundle &dev, double rho, double u, double v, double w, double p)
{
    const double ps[kNv] = {rho, u, v, w, p};
    CUDA_CHECK(cudaMemcpy(dev.d_postshock, ps, sizeof(ps), cudaMemcpyHostToDevice));
}

void upload_U(DeviceBundle &dev, const std::vector<double> &U)
{
    const std::size_t count = static_cast<std::size_t>(dev.n_cells) * kNv;
    require_size(U.size(), count, "U");
    CUDA_CHECK(cudaMemcpy(dev.d_U, U.data(), sizeof(double) * count,
                          cudaMemcpyHostToDevice));
}

void download_U(DeviceBundle &dev, std::vector<double> &U)
{
    const std::size_t count = static_cast<std::size_t>(dev.n_cells) * kNv;
    U.resize(count);
    CUDA_CHECK(cudaMemcpy(U.data(), dev.d_U, sizeof(double) * count,
                          cudaMemcpyDeviceToHost));
}

void upload_R(DeviceBundle &dev, const std::vector<double> &R)
{
    const std::size_t count = static_cast<std::size_t>(dev.n_cells) * kNv;
    require_size(R.size(), count, "R");
    CUDA_CHECK(cudaMemcpy(dev.d_R, R.data(), sizeof(double) * count,
                          cudaMemcpyHostToDevice));
}

void launch_net_flux(DeviceBundle &dev, int scheme, int order,
                     bool weno_char, bool weno_z,
                     int ricca_sensor, double ricca_rh_thresh, int weno_hybrid)
{
    if (dev.nx <= 0 || dev.ny <= 0 || dev.nz <= 0) {
        std::fprintf(stderr, "Euler3D CUDA: call set_grid_size before launch_net_flux\n");
        std::abort();
    }
    const int blocks = (dev.n_cells + kBlockSize - 1) / kBlockSize;
    net_flux_kernel<<<blocks, kBlockSize>>>(
        dev.d_U, dev.d_R, dev.d_neigh, dev.d_nxyz, dev.d_area, dev.d_bc,
        dev.d_freestream, dev.d_postshock, dev.n_cells, dev.nx, dev.ny, dev.nz,
        scheme, order, weno_char, weno_z, ricca_sensor, ricca_rh_thresh, weno_hybrid);
    CUDA_CHECK(cudaGetLastError());
}

void launch_update(DeviceBundle &dev, double dt)
{
    const int blocks = (dev.n_cells + kBlockSize - 1) / kBlockSize;
    update_kernel<<<blocks, kBlockSize>>>(dev.d_U, dev.d_R, dev.d_volume, dt, dev.n_cells);
    CUDA_CHECK(cudaGetLastError());
}

void launch_viscous_flux(DeviceBundle &dev, double mu, double Pr, int turb_model,
                         double Cs, double Cw, double Cv, double Pr_t, double kappa,
                         double A_plus, double Uinf)
{
    if (mu <= 0.0)
        return;
    const double cp = kGamma / (kGamma - 1.0);
    const int blocks = (dev.n_cells + kBlockSize - 1) / kBlockSize;
    viscous_grad_kernel<<<blocks, kBlockSize>>>(dev.d_U, dev.d_neigh, dev.d_bc, dev.d_xc,
                                                dev.d_yc, dev.d_zc, dev.d_grad, dev.dx,
                                                dev.dy, dev.dz, dev.n_cells);
    CUDA_CHECK(cudaGetLastError());
    sgs_viscosity_kernel<<<blocks, kBlockSize>>>(
        dev.d_U, dev.d_grad, dev.d_volume, dev.d_wall_dist, dev.d_mu_sgs, turb_model, mu,
        Cs, Cw, Cv, kappa, A_plus, Uinf, dev.n_cells);
    CUDA_CHECK(cudaGetLastError());
    viscous_flux_kernel<<<blocks, kBlockSize>>>(dev.d_U, dev.d_grad, dev.d_neigh, dev.d_nxyz,
                                                dev.d_area, dev.d_mu_sgs, dev.d_R, mu, cp, Pr,
                                                Pr_t, dev.n_cells);
    CUDA_CHECK(cudaGetLastError());
}

double launch_min_dt(DeviceBundle &dev, double cfl, bool viscous, double mu)
{
    const std::size_t count = static_cast<std::size_t>(dev.n_cells);
    double *d_dt = nullptr;
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_dt), sizeof(double) * count));
    const int blocks = (dev.n_cells + kBlockSize - 1) / kBlockSize;
    if (viscous && mu > 0.0 && dev.nx > 0) {
        viscous_dt_kernel<<<blocks, kBlockSize>>>(
            dev.d_U, dev.d_volume, dev.d_nxyz, dev.d_area, dev.d_xc, dev.d_yc, dev.d_zc,
            d_dt, cfl, mu, dev.dx, dev.dy, dev.dz, dev.nx, dev.ny, dev.nz, dev.n_cells);
    } else {
        local_dt_kernel<<<blocks, kBlockSize>>>(dev.d_U, dev.d_volume, dev.d_nxyz,
                                                dev.d_area, d_dt, cfl, dev.n_cells);
    }
    CUDA_CHECK(cudaGetLastError());
    std::vector<double> host_dt(count);
    CUDA_CHECK(cudaMemcpy(host_dt.data(), d_dt, sizeof(double) * count,
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_dt));
    return *std::min_element(host_dt.begin(), host_dt.end());
}

} // namespace euler3d_cuda
