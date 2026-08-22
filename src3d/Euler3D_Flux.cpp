/**
 * @file Euler3D_Flux.cpp
 * @brief Host inviscid fluxes and structured-grid reconstruction for Euler3D.
 */
#include "Euler3D.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace euler3d
{
namespace
{

constexpr double kTiny = 1.0e-14;

inline void fill_prim(const double U[5], double prim[6])
{
    const double rho = std::max(U[0], kTiny);
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    const double w = U[3] / rho;
    const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
    const double p = std::max((kGamma - 1.0) * (U[4] - kinetic), kTiny);
    prim[0] = rho;
    prim[1] = u;
    prim[2] = v;
    prim[3] = w;
    prim[4] = p;
    prim[5] = std::sqrt(kGamma * p / rho);
}

inline void U_from_prim(double rho, double u, double v, double w, double p, double U[5])
{
    rho = std::max(rho, kTiny);
    p = std::max(p, kTiny);
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = rho * w;
    U[4] = p / (kGamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

inline void physical_flux(const double U[5], const double p[6],
                          double nx, double ny, double nz, double F[5])
{
    const double un = p[1] * nx + p[2] * ny + p[3] * nz;
    F[0] = U[0] * un;
    F[1] = U[1] * un + p[4] * nx;
    F[2] = U[2] * un + p[4] * ny;
    F[3] = U[3] * un + p[4] * nz;
    F[4] = (U[4] + p[4]) * un;
}

inline void slip_reflect(const double U_in[5], double nx, double ny, double nz,
                         double U_out[5])
{
    double q[6];
    fill_prim(U_in, q);
    const double un = q[1] * nx + q[2] * ny + q[3] * nz;
    U_from_prim(q[0], q[1] - 2.0 * un * nx, q[2] - 2.0 * un * ny,
                q[3] - 2.0 * un * nz, q[4], U_out);
}

inline void no_slip_reflect(const double U_in[5], double U_out[5])
{
    double q[6];
    fill_prim(U_in, q);
    U_from_prim(q[0], -q[1], -q[2], -q[3], q[4], U_out);
}

inline void extrapolate_ghost(const double U_in[5], double U_out[5])
{
    std::copy(U_in, U_in + 5, U_out);
}

template <class T>
auto cfg_rho(const T &cfg, int) -> decltype(cfg.freestream_rho)
{
    return cfg.freestream_rho;
}
template <class T>
double cfg_rho(const T &, long) { return 1.0; }

template <class T>
auto cfg_u(const T &cfg, int) -> decltype(cfg.freestream_u)
{
    return cfg.freestream_u;
}
template <class T>
double cfg_u(const T &, long) { return 0.0; }

template <class T>
auto cfg_v(const T &cfg, int) -> decltype(cfg.freestream_v)
{
    return cfg.freestream_v;
}
template <class T>
double cfg_v(const T &, long) { return 0.0; }

template <class T>
auto cfg_w(const T &cfg, int) -> decltype(cfg.freestream_w)
{
    return cfg.freestream_w;
}
template <class T>
double cfg_w(const T &, long) { return 0.0; }

template <class T>
auto cfg_p(const T &cfg, int) -> decltype(cfg.freestream_p)
{
    return cfg.freestream_p;
}
template <class T>
double cfg_p(const T &, long) { return 1.0; }

template <class T>
auto cfg_order(const T &cfg, int) -> decltype(cfg.order)
{
    return cfg.order;
}
template <class T>
int cfg_order(const T &, long) { return 0; }

template <class T>
auto boundary_type(const T &m, int face, int) -> decltype(m.bc_type.size(), int())
{
    return face < static_cast<int>(m.bc_type.size()) ? m.bc_type[face] : 3;
}
template <class T>
int boundary_type(const T &m, int face, long)
{
    return face < static_cast<int>(m.wall_face.size()) && m.wall_face[face] ? 1 : 3;
}

inline void freestream_ghost(const Config &cfg, double U_out[5])
{
    U_from_prim(cfg_rho(cfg, 0), cfg_u(cfg, 0), cfg_v(cfg, 0), cfg_w(cfg, 0),
                cfg_p(cfg, 0), U_out);
}

inline void postshock_ghost(const Config &cfg, double U_out[5])
{
    U_from_prim(cfg.post_rho, cfg.post_u, cfg.post_v, cfg.post_w, cfg.post_p, U_out);
}

inline void boundary_ghost(const double U_in[5], const Mesh &m, const Config &cfg,
                           int c, int f, double nx, double ny, double nz, double U_out[5])
{
    const int type = boundary_type(m, c * kFaces + f, 0);
    if (type == 1)
        slip_reflect(U_in, nx, ny, nz, U_out);
    else if (type == 2)
        no_slip_reflect(U_in, U_out);
    else if (type == 4)
        freestream_ghost(cfg, U_out);
    else if (type == 5)
        postshock_ghost(cfg, U_out);
    else
        extrapolate_ghost(U_in, U_out);
}

inline int cell_id(int i, int j, int k, const Mesh &m)
{
    if (i < 0 || i >= m.nx || j < 0 || j >= m.ny || k < 0 || k >= m.nz)
        return -1;
    return i + j * m.nx + k * m.nx * m.ny;
}

inline void cell_ijk(int c, const Mesh &m, int &i, int &j, int &k)
{
    const int plane = m.nx * m.ny;
    k = c / plane;
    const int r = c - k * plane;
    j = r / m.nx;
    i = r - j * m.nx;
}

inline int line_cell(int c, int f, int offset, const Mesh &m)
{
    int i, j, k;
    cell_ijk(c, m, i, j, k);
    const int sign = (f & 1) ? 1 : -1;
    if (f < 2)
        i += sign * offset;
    else if (f < 4)
        j += sign * offset;
    else
        k += sign * offset;
    return cell_id(i, j, k, m);
}

inline double minmod(double a, double b)
{
    if (a * b <= 0.0)
        return 0.0;
    return std::copysign(std::min(std::fabs(a), std::fabs(b)), a);
}

inline void primitive_at(const State &st, int c, double q[6])
{
    if (static_cast<int>(st.prim.size()) >= (c + 1) * 6)
        std::copy(&st.prim[c * 6], &st.prim[c * 6 + 6], q);
    else
        fill_prim(&st.U[c * 5], q);
}

inline bool admissible(const double U[5])
{
    if (!(U[0] > kTiny))
        return false;
    const double kinetic = 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0];
    return (kGamma - 1.0) * (U[4] - kinetic) > kTiny;
}

inline double weno5_js(double v0, double v1, double v2, double v3, double v4)
{
    /* Jiang–Shu WENO5 for u_{i+1/2}^- from stencil {i-2..i+2}. */
    const double d0 = 0.1, d1 = 0.6, d2 = 0.3;
    const double p0 = (2.0 * v0 - 7.0 * v1 + 11.0 * v2) / 6.0;
    const double p1 = (-v1 + 5.0 * v2 + 2.0 * v3) / 6.0;
    const double p2 = (2.0 * v2 + 5.0 * v3 - v4) / 6.0;
    const double b0 = (13.0 / 12.0) * std::pow(v0 - 2.0 * v1 + v2, 2) +
                      0.25 * std::pow(v0 - 4.0 * v1 + 3.0 * v2, 2);
    const double b1 = (13.0 / 12.0) * std::pow(v1 - 2.0 * v2 + v3, 2) +
                      0.25 * std::pow(v1 - v3, 2);
    const double b2 = (13.0 / 12.0) * std::pow(v2 - 2.0 * v3 + v4, 2) +
                      0.25 * std::pow(3.0 * v2 - 4.0 * v3 + v4, 2);
    const double eps = 1.0e-6;
    const double a0 = d0 / ((eps + b0) * (eps + b0));
    const double a1 = d1 / ((eps + b1) * (eps + b1));
    const double a2 = d2 / ((eps + b2) * (eps + b2));
    const double sum = a0 + a1 + a2;
    if (!(sum > eps))
        return v2;
    return (a0 * p0 + a1 * p1 + a2 * p2) / sum;
}

inline double weno5_z(double v0, double v1, double v2, double v3, double v4)
{
    /* Borges WENO-Z for u_{i+1/2}^- from stencil {i-2..i+2}. */
    const double d0 = 0.1, d1 = 0.6, d2 = 0.3;
    const double p0 = (2.0 * v0 - 7.0 * v1 + 11.0 * v2) / 6.0;
    const double p1 = (-v1 + 5.0 * v2 + 2.0 * v3) / 6.0;
    const double p2 = (2.0 * v2 + 5.0 * v3 - v4) / 6.0;
    const double b0 = (13.0 / 12.0) * std::pow(v0 - 2.0 * v1 + v2, 2) +
                      0.25 * std::pow(v0 - 4.0 * v1 + 3.0 * v2, 2);
    const double b1 = (13.0 / 12.0) * std::pow(v1 - 2.0 * v2 + v3, 2) +
                      0.25 * std::pow(v1 - v3, 2);
    const double b2 = (13.0 / 12.0) * std::pow(v2 - 2.0 * v3 + v4, 2) +
                      0.25 * std::pow(3.0 * v2 - 4.0 * v3 + v4, 2);
    const double eps = 1.0e-40;
    const double tau5 = std::fabs(b0 - b2);
    const double a0 = d0 * (1.0 + std::pow(tau5 / (eps + b0), 2));
    const double a1 = d1 * (1.0 + std::pow(tau5 / (eps + b1), 2));
    const double a2 = d2 * (1.0 + std::pow(tau5 / (eps + b2), 2));
    const double sum = a0 + a1 + a2;
    if (!(sum > 0.0) || !std::isfinite(sum))
        return v2;
    const double out = (a0 * p0 + a1 * p1 + a2 * p2) / sum;
    return std::isfinite(out) ? out : v2;
}

inline double weno5_scalar(double v0, double v1, double v2, double v3, double v4, bool use_z)
{
    return use_z ? weno5_z(v0, v1, v2, v3, v4) : weno5_js(v0, v1, v2, v3, v4);
}

/** Roe-averaged L,R at face; char = L*U, U = R*char (row-major 5×5). */
inline void build_char_matrices(const double qL[6], const double qR[6], double nx, double ny,
                                double nz, double L[25], double R[25])
{
    const double rL = std::sqrt(std::max(qL[0], kTiny));
    const double rR = std::sqrt(std::max(qR[0], kTiny));
    const double den = std::max(rL + rR, kTiny);
    const double u = (rL * qL[1] + rR * qR[1]) / den;
    const double v = (rL * qL[2] + rR * qR[2]) / den;
    const double w = (rL * qL[3] + rR * qR[3]) / den;
    const double a = std::max((rL * qL[5] + rR * qR[5]) / den, kTiny);
    const double un = u * nx + v * ny + w * nz;
    const double vel2 = u * u + v * v + w * w;
    const double ek = 0.5 * vel2;
    const double H = a * a / (kGamma - 1.0) + ek;
    const double gm1 = kGamma - 1.0;
    const double a2 = a * a;
    const double t1 = 0.5 / a2;
    const double t2 = gm1 * t1;

    double tx, ty, tz;
    if (std::fabs(nx) > std::fabs(nz))
    {
        const double inv = 1.0 / std::sqrt(nx * nx + ny * ny + 1e-30);
        tx = -ny * inv;
        ty = nx * inv;
        tz = 0.0;
    }
    else
    {
        const double inv = 1.0 / std::sqrt(ny * ny + nz * nz + 1e-30);
        tx = 0.0;
        ty = -nz * inv;
        tz = ny * inv;
    }
    const double sx = ny * tz - nz * ty;
    const double sy = nz * tx - nx * tz;
    const double sz = nx * ty - ny * tx;
    const double ut = u * tx + v * ty + w * tz;
    const double us = u * sx + v * sy + w * sz;

    /* Columns of R: acoustic-, entropy, shear-t, shear-s, acoustic+ */
    const double cols[5][5] = {
        {1.0, u - a * nx, v - a * ny, w - a * nz, H - a * un},
        {1.0, u, v, w, ek},
        {0.0, tx, ty, tz, ut},
        {0.0, sx, sy, sz, us},
        {1.0, u + a * nx, v + a * ny, w + a * nz, H + a * un}};
    for (int col = 0; col < 5; ++col)
        for (int row = 0; row < 5; ++row)
            R[row * 5 + col] = cols[col][row];

    /* Rows of L */
    L[0 * 5 + 0] = t2 * ek + t1 * a * un;
    L[0 * 5 + 1] = -t2 * u - t1 * a * nx;
    L[0 * 5 + 2] = -t2 * v - t1 * a * ny;
    L[0 * 5 + 3] = -t2 * w - t1 * a * nz;
    L[0 * 5 + 4] = t2;

    L[1 * 5 + 0] = 1.0 - 2.0 * t2 * ek;
    L[1 * 5 + 1] = 2.0 * t2 * u;
    L[1 * 5 + 2] = 2.0 * t2 * v;
    L[1 * 5 + 3] = 2.0 * t2 * w;
    L[1 * 5 + 4] = -2.0 * t2;

    L[2 * 5 + 0] = -ut;
    L[2 * 5 + 1] = tx;
    L[2 * 5 + 2] = ty;
    L[2 * 5 + 3] = tz;
    L[2 * 5 + 4] = 0.0;

    L[3 * 5 + 0] = -us;
    L[3 * 5 + 1] = sx;
    L[3 * 5 + 2] = sy;
    L[3 * 5 + 3] = sz;
    L[3 * 5 + 4] = 0.0;

    L[4 * 5 + 0] = t2 * ek - t1 * a * un;
    L[4 * 5 + 1] = -t2 * u + t1 * a * nx;
    L[4 * 5 + 2] = -t2 * v + t1 * a * ny;
    L[4 * 5 + 3] = -t2 * w + t1 * a * nz;
    L[4 * 5 + 4] = t2;
}

inline void mat5_vec(const double M[25], const double x[5], double y[5])
{
    for (int i = 0; i < 5; ++i)
    {
        y[i] = 0.0;
        for (int j = 0; j < 5; ++j)
            y[i] += M[i * 5 + j] * x[j];
    }
}

inline void weno5_face_states(const State &st, const Mesh &m, const Config &cfg, int c, int f,
                              int nb, double UL[5], double UR[5])
{
    int ids[6];
    for (int s = -2; s <= 3; ++s)
    {
        ids[s + 2] = line_cell(c, f, s, m);
        if (ids[s + 2] < 0)
            return;
    }

    double Usten[6][5];
    for (int s = 0; s < 6; ++s)
        for (int v = 0; v < 5; ++v)
            Usten[s][v] = st.U[ids[s] * 5 + v];

    const bool use_z = cfg.weno_z;
    auto recon_LR = [&](const double sten[6][5], double outL[5], double outR[5]) {
        for (int v = 0; v < 5; ++v)
        {
            outL[v] = weno5_scalar(sten[0][v], sten[1][v], sten[2][v], sten[3][v], sten[4][v],
                                   use_z);
            /* Right state: mirror stencil into the same left-biased operator. */
            outR[v] = weno5_scalar(sten[5][v], sten[4][v], sten[3][v], sten[2][v], sten[1][v],
                                   use_z);
        }
    };

    if (!cfg.weno_char)
    {
        recon_LR(Usten, UL, UR);
    }
    else
    {
        double qL[6], qR[6], Lmat[25], Rmat[25];
        fill_prim(Usten[2], qL);
        fill_prim(Usten[3], qR);
        const int ni = (c * kFaces + f) * 3;
        build_char_matrices(qL, qR, m.nxyz[ni], m.nxyz[ni + 1], m.nxyz[ni + 2], Lmat, Rmat);

        double Wsten[6][5];
        for (int s = 0; s < 6; ++s)
            mat5_vec(Lmat, Usten[s], Wsten[s]);

        double WL[5], WR[5];
        recon_LR(Wsten, WL, WR);
        mat5_vec(Rmat, WL, UL);
        mat5_vec(Rmat, WR, UR);
    }

    if (!admissible(UL) || !admissible(UR) || !std::isfinite(UL[0]) || !std::isfinite(UR[0]))
    {
        std::copy(&st.U[c * 5], &st.U[c * 5 + 5], UL);
        std::copy(&st.U[nb * 5], &st.U[nb * 5 + 5], UR);
        return;
    }

    /* Paper WENO hybrid: blend toward FO near discontinuities.
       Indicators are face-based (reconstructed L/R). Combo uses shock
       (p, u_n) always; density/entropy only when φ_P is small (contacts),
       so post-shock entropy noise does not freeze the whole field to FO. */
    if (cfg.weno_hybrid > 0)
    {
        const double *UfoL = Usten[2];
        const double *UfoR = Usten[3];
        double qL[6], qR[6];
        fill_prim(UL, qL);
        fill_prim(UR, qR);
        const int ni = (c * kFaces + f) * 3;
        const double nx = m.nxyz[ni], ny = m.nxyz[ni + 1], nz = m.nxyz[ni + 2];
        const double unL = qL[1] * nx + qL[2] * ny + qL[3] * nz;
        const double unR = qR[1] * nx + qR[2] * ny + qR[3] * nz;
        const double a_ref = 0.5 * (qL[5] + qR[5]) + 1.0;
        const double phi_P =
            std::fabs(qR[4] - qL[4]) / (0.5 * (std::fabs(qL[4]) + std::fabs(qR[4])) + 1.0);
        const double phi_rho =
            std::fabs(qR[0] - qL[0]) / (0.5 * (std::fabs(qL[0]) + std::fabs(qR[0])) + 1.0);
        const double sL = qL[4] / std::pow(std::max(qL[0], kTiny), kGamma);
        const double sR = qR[4] / std::pow(std::max(qR[0], kTiny), kGamma);
        const double phi_S =
            std::fabs(sR - sL) / (0.5 * (std::fabs(sL) + std::fabs(sR)) + 1.0);
        const double phi_u = std::fabs(unR - unL) / a_ref;
        double phi = phi_P;
        if (cfg.weno_hybrid == 2)
        {
            phi = std::max(phi_P, phi_u);
            if (phi_P < 0.05)
                phi = std::max({phi, phi_rho, phi_S});
        }

        for (int v = 0; v < 5; ++v)
        {
            const double C = (v >= 1 && v <= 3) ? 120.0 : 50.0;
            const double theta = 1.0 / (1.0 + C * phi);
            UL[v] = theta * UL[v] + (1.0 - theta) * UfoL[v];
            UR[v] = theta * UR[v] + (1.0 - theta) * UfoR[v];
            const double lo = std::min({Usten[1][v], Usten[2][v], Usten[3][v], Usten[4][v]});
            const double hi = std::max({Usten[1][v], Usten[2][v], Usten[3][v], Usten[4][v]});
            UL[v] = std::min(hi, std::max(lo, UL[v]));
            UR[v] = std::min(hi, std::max(lo, UR[v]));
        }
        if (!admissible(UL) || !admissible(UR))
        {
            std::copy(UfoL, UfoL + 5, UL);
            std::copy(UfoR, UfoR + 5, UR);
        }
    }
}

void face_flux_llf(const double UL[5], const double UR[5],
                   double nx, double ny, double nz, double area, double flux[5])
{
    double qL[6], qR[6], FL[5], FR[5];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    physical_flux(UL, qL, nx, ny, nz, FL);
    physical_flux(UR, qR, nx, ny, nz, FR);
    const double unL = qL[1] * nx + qL[2] * ny + qL[3] * nz;
    const double unR = qR[1] * nx + qR[2] * ny + qR[3] * nz;
    const double smax = std::max(std::fabs(unL) + qL[5], std::fabs(unR) + qR[5]);
    for (int v = 0; v < 5; ++v)
        flux[v] = 0.5 * (FL[v] + FR[v] - smax * (UR[v] - UL[v])) * area;
}

/** Normalized Rankine–Hugoniot residual on the normal momentum flux. */
inline double ricca_phi_rh(const double qL[6], const double qR[6], double unL, double unR)
{
    const double massL = qL[0] * unL;
    const double massR = qR[0] * unR;
    const double FnL = massL * unL + qL[4];
    const double FnR = massR * unR + qR[4];
    const double dr = qR[0] - qL[0];
    const double s = (std::fabs(dr) > kTiny) ? (massR - massL) / dr : 0.5 * (unL + unR);
    const double resid = std::fabs(FnR - FnL - s * (massR - massL));
    return resid / (0.5 * (std::fabs(FnL) + std::fabs(FnR)) + 1.0);
}

/**
 * Paper RICCA flux:
 *   F = 1/2 (FL+FR) - 1/2 α (UR-UL) dℓ
 * α_contact = max(|VnL|,|VnR|)
 * α_acoustic = max(|VnL|+aL, |VnR|+aR)   [single-phase; multiphase c(1-c)|Vrn| = 0]
 * Legacy: pressure jump → acoustic on all single-phase equations (ρ is gas mass).
 *          no pressure jump → contact speed (keeps contacts sharp).
 * RICCA-RH: binary σ from φ_P or φ_RH.
 */
void face_flux_ricca(const double UL[5], const double UR[5],
                     double nx, double ny, double nz, double area, double flux[5],
                     int sensor_type, double rh_threshold)
{
    double qL[6], qR[6], FL[5], FR[5];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    physical_flux(UL, qL, nx, ny, nz, FL);
    physical_flux(UR, qR, nx, ny, nz, FR);
    const double unL = qL[1] * nx + qL[2] * ny + qL[3] * nz;
    const double unR = qR[1] * nx + qR[2] * ny + qR[3] * nz;
    const double a_contact = std::max(std::fabs(unL), std::fabs(unR));
    const double a_acoustic =
        std::max(std::fabs(unL) + qL[5], std::fabs(unR) + qR[5]);
    const double phi_P =
        std::fabs(qR[4] - qL[4]) / (0.5 * (std::fabs(qL[4]) + std::fabs(qR[4])) + 1.0);

    bool use_acoustic = false;
    if (sensor_type == 1)
    {
        const double phi_RH = ricca_phi_rh(qL, qR, unL, unR);
        use_acoustic = (phi_P > 1.0e-8) || (phi_RH > rh_threshold);
    }
    else
    {
        use_acoustic = (phi_P > 1.0e-8);
    }

    const double alpha = use_acoustic ? a_acoustic : a_contact;
    for (int v = 0; v < 5; ++v)
        flux[v] = 0.5 * (FL[v] + FR[v]) * area - 0.5 * alpha * area * (UR[v] - UL[v]);
}

inline double entropy_abs(double lambda, double delta)
{
    const double a = std::fabs(lambda);
    return a < delta ? 0.5 * (a * a / delta + delta) : a;
}

void face_flux_roe(const double UL[5], const double UR[5],
                   double nx, double ny, double nz, double area, double flux[5])
{
    double qL[6], qR[6], FL[5], FR[5];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    physical_flux(UL, qL, nx, ny, nz, FL);
    physical_flux(UR, qR, nx, ny, nz, FR);

    double tx, ty, tz;
    if (std::fabs(nx) > std::fabs(nz))
    {
        const double inv = 1.0 / std::sqrt(nx * nx + ny * ny);
        tx = -ny * inv;
        ty = nx * inv;
        tz = 0.0;
    }
    else
    {
        const double inv = 1.0 / std::sqrt(ny * ny + nz * nz);
        tx = 0.0;
        ty = -nz * inv;
        tz = ny * inv;
    }
    const double sx = ny * tz - nz * ty;
    const double sy = nz * tx - nx * tz;
    const double sz = nx * ty - ny * tx;

    const double rL = std::sqrt(qL[0]), rR = std::sqrt(qR[0]);
    const double den = rL + rR;
    const double u = (rL * qL[1] + rR * qR[1]) / den;
    const double v = (rL * qL[2] + rR * qR[2]) / den;
    const double w = (rL * qL[3] + rR * qR[3]) / den;
    const double HL = (UL[4] + qL[4]) / qL[0];
    const double HR = (UR[4] + qR[4]) / qR[0];
    const double H = (rL * HL + rR * HR) / den;
    const double vel2 = u * u + v * v + w * w;
    const double a = std::sqrt(std::max((kGamma - 1.0) * (H - 0.5 * vel2), kTiny));
    const double rho = rL * rR;
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
    const double dp = qR[4] - qL[4];
    const double a2 = a * a;
    const double alpha[5] = {
        (dp - rho * a * dun) / (2.0 * a2),
        dr - dp / a2,
        rho * dut,
        rho * dus,
        (dp + rho * a * dun) / (2.0 * a2)};
    const double delta = std::max(0.1 * a, kTiny);
    const double lambda[5] = {
        entropy_abs(un - a, delta), entropy_abs(un, delta), entropy_abs(un, delta),
        entropy_abs(un, delta), entropy_abs(un + a, delta)};
    const double wave[5][5] = {
        {1.0, u - a * nx, v - a * ny, w - a * nz, H - a * un},
        {1.0, u, v, w, 0.5 * vel2},
        {0.0, tx, ty, tz, ut},
        {0.0, sx, sy, sz, us},
        {1.0, u + a * nx, v + a * ny, w + a * nz, H + a * un}};
    for (int n = 0; n < 5; ++n)
    {
        double diss = 0.0;
        for (int k = 0; k < 5; ++k)
            diss += lambda[k] * alpha[k] * wave[k][n];
        flux[n] = (0.5 * (FL[n] + FR[n]) - 0.5 * diss) * area;
    }
}

void vanleer_split(const double U[5], const double q[6],
                   double nx, double ny, double nz, bool plus, double F[5])
{
    const double un = q[1] * nx + q[2] * ny + q[3] * nz;
    const double M = un / q[5];
    const bool all = plus ? (M >= 1.0) : (M <= -1.0);
    const bool none = plus ? (M <= -1.0) : (M >= 1.0);
    if (all)
    {
        physical_flux(U, q, nx, ny, nz, F);
        return;
    }
    if (none)
    {
        std::fill(F, F + 5, 0.0);
        return;
    }
    const double s = plus ? 1.0 : -1.0;
    const double mass = s * 0.25 * q[0] * q[5] * (M + s) * (M + s);
    const double ps = 0.25 * q[4] * (M + s) * (M + s) * (2.0 - s * M);
    const double vt2 = q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - un * un;
    const double hs = std::pow((kGamma - 1.0) * un + s * 2.0 * q[5], 2) /
                          (2.0 * (kGamma * kGamma - 1.0)) +
                      0.5 * vt2;
    F[0] = mass;
    F[1] = mass * q[1] + ps * nx;
    F[2] = mass * q[2] + ps * ny;
    F[3] = mass * q[3] + ps * nz;
    F[4] = mass * hs;
}

void face_flux_vanleer(const double UL[5], const double UR[5],
                       double nx, double ny, double nz, double area, double flux[5])
{
    double qL[6], qR[6], Fp[5], Fm[5];
    fill_prim(UL, qL);
    fill_prim(UR, qR);
    vanleer_split(UL, qL, nx, ny, nz, true, Fp);
    vanleer_split(UR, qR, nx, ny, nz, false, Fm);
    for (int v = 0; v < 5; ++v)
        flux[v] = (Fp[v] + Fm[v]) * area;
}

} // namespace

void fill_prim_cell(const double U[5], double prim[6])
{
    fill_prim(U, prim);
}

void U_from_prim_vals(double rho, double u, double v, double w, double p, double U[5])
{
    U_from_prim(rho, u, v, w, p, U);
}

void get_LR_states(const State &st, const Mesh &m, const Config &cfg,
                   int c, int f, double UL[5], double UR[5])
{
    std::copy(&st.U[c * 5], &st.U[c * 5 + 5], UL);
    const int nb = m.neigh[c * kFaces + f];
    const int ni = (c * kFaces + f) * 3;
    const double nx = m.nxyz[ni], ny = m.nxyz[ni + 1], nz = m.nxyz[ni + 2];
    if (nb < 0)
    {
        boundary_ghost(UL, m, cfg, c, f, nx, ny, nz, UR);
        return;
    }
    std::copy(&st.U[nb * 5], &st.U[nb * 5 + 5], UR);

    const int order = cfg_order(cfg, 0);
    if (order == 1)
    {
        const int cm1 = line_cell(c, f, -1, m);
        const int cp1 = line_cell(c, f, 1, m);
        const int cp2 = line_cell(c, f, 2, m);
        if (cm1 < 0 || cp1 < 0 || cp2 < 0)
            return;
        double qm1[6], q0[6], q1[6], q2[6], qL[5], qR[5];
        primitive_at(st, cm1, qm1);
        primitive_at(st, c, q0);
        primitive_at(st, cp1, q1);
        primitive_at(st, cp2, q2);
        for (int v = 0; v < 5; ++v)
        {
            qL[v] = q0[v] + 0.5 * minmod(q0[v] - qm1[v], q1[v] - q0[v]);
            qR[v] = q1[v] - 0.5 * minmod(q1[v] - q0[v], q2[v] - q1[v]);
        }
        if (qL[0] <= kTiny || qL[4] <= kTiny || qR[0] <= kTiny || qR[4] <= kTiny)
            return;
        U_from_prim(qL[0], qL[1], qL[2], qL[3], qL[4], UL);
        U_from_prim(qR[0], qR[1], qR[2], qR[3], qR[4], UR);
    }
    else if (order == 2)
    {
        weno5_face_states(st, m, cfg, c, f, nb, UL, UR);
    }
}

void face_flux(const double UL[5], const double UR[5],
               double nx, double ny, double nz, double area, int scheme, double flux[5])
{
    Config cfg;
    cfg.scheme = scheme;
    face_flux(UL, UR, nx, ny, nz, area, cfg, flux);
}

void face_flux(const double UL[5], const double UR[5],
               double nx, double ny, double nz, double area, const Config &cfg, double flux[5])
{
    if (cfg.scheme == 1)
        face_flux_ricca(UL, UR, nx, ny, nz, area, flux, cfg.ricca_sensor,
                        cfg.ricca_rh_threshold);
    else if (cfg.scheme == 2)
        face_flux_roe(UL, UR, nx, ny, nz, area, flux);
    else if (cfg.scheme == 3)
        face_flux_vanleer(UL, UR, nx, ny, nz, area, flux);
    else
        face_flux_llf(UL, UR, nx, ny, nz, area, flux);
}

void net_flux_host(const State &st, const Mesh &m, const Config &cfg, std::vector<double> &R)
{
    R.assign(m.n_cells * kNv, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int c = 0; c < m.n_cells; ++c)
        for (int f = 0; f < kFaces; ++f)
        {
            const int fi = c * kFaces + f;
            const int ni = fi * 3;
            double UL[5], UR[5], flux[5];
            get_LR_states(st, m, cfg, c, f, UL, UR);
            face_flux(UL, UR, m.nxyz[ni], m.nxyz[ni + 1], m.nxyz[ni + 2], m.area[fi], cfg,
                      flux);
            for (int v = 0; v < kNv; ++v)
                R[c * kNv + v] += flux[v];
        }
}

} // namespace euler3d
