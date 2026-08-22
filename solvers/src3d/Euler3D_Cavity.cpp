/**
 * @file Euler3D_Cavity.cpp
 * @brief Four-block Cartesian open cavity (90° lips).
 *
 * Blocks packed into one ijk array:
 *   bay     [xs,xe] x [-D, 0]
 *   over    [xs,xe] x [ 0, Ly]
 *   fore    [0, xs] x [ 0, Ly]
 *   aft     [xe,Lx] x [ 0, Ly]
 * Cells under the plate (x<xs or x>xe, y<0) are solid (inactive).
 */
#include "Euler3D.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace euler3d
{
namespace
{

inline int nid(int i, int j, int k, int npx, int npy)
{
    return i + j * npx + k * npx * npy;
}

inline int cid(int i, int j, int k, int nx, int ny)
{
    return i + j * nx + k * nx * ny;
}

double cluster_right(double u, double p)
{
    u = std::min(1.0, std::max(0.0, u));
    p = std::max(p, 1.0);
    return 1.0 - std::pow(1.0 - u, p);
}

double cluster_left(double u, double p)
{
    u = std::min(1.0, std::max(0.0, u));
    p = std::max(p, 1.0);
    return std::pow(u, p);
}

double cluster_both(double u, double p)
{
    u = std::min(1.0, std::max(0.0, u));
    p = std::max(p, 1.0);
    const double a = std::pow(u, p);
    const double b = std::pow(1.0 - u, p);
    return a / std::max(a + b, 1.0e-30);
}

void make_x_nodes(int npx, double Lx, double xs, double xe, double p, std::vector<double> &xn)
{
    xn.assign(static_cast<std::size_t>(npx), 0.0);
    const int nx = npx - 1;
    int n1 = std::max(4, static_cast<int>(std::lround(nx * (xs / Lx))));
    int n3 = std::max(4, static_cast<int>(std::lround(nx * ((Lx - xe) / Lx))));
    int n2 = nx - n1 - n3;
    if (n2 < 8)
    {
        n2 = 8;
        n1 = std::max(4, (nx - n2) / 2);
        n3 = nx - n1 - n2;
    }
    xn[0] = 0.0;
    for (int k = 1; k <= n1; ++k)
        xn[k] = xs * cluster_right(static_cast<double>(k) / n1, p);
    xn[n1] = xs;
    for (int k = 1; k <= n2; ++k)
        xn[n1 + k] = xs + (xe - xs) * cluster_both(static_cast<double>(k) / n2, p);
    xn[n1 + n2] = xe;
    for (int k = 1; k <= n3; ++k)
        xn[n1 + n2 + k] = xe + (Lx - xe) * cluster_left(static_cast<double>(k) / n3, p);
    xn[nx] = Lx;
    for (int k = 1; k < npx; ++k)
        if (!(xn[k] > xn[k - 1]))
            xn[k] = xn[k - 1] + 1.0e-12 * Lx;
}

int parse_bc_local(const std::string &s)
{
    if (s == "noslip" || s == "no_slip" || s == "wall_viscous")
        return BC_NOSLIP;
    if (s == "extrapolate" || s == "outlet" || s == "outflow")
        return BC_EXTRAP;
    if (s == "freestream" || s == "inlet" || s == "farfield")
        return BC_FREESTREAM;
    return BC_SLIP;
}

} // namespace

Mesh make_cavity_mesh(const Config &cfg)
{
    if (cfg.nx < 8 || cfg.ny < 12 || cfg.nz < 1)
        throw std::runtime_error("Cavity mesh needs nx>=8, ny>=12, nz>=1");

    Mesh m;
    m.nx = cfg.nx;
    m.ny = cfg.ny;
    m.nz = cfg.nz;
    m.Lx = cfg.Lx;
    m.Ly = cfg.Ly;
    m.Lz = cfg.Lz;
    m.n_cells = m.nx * m.ny * m.nz;
    m.dx = cfg.Lx / m.nx;
    m.dy = (cfg.Ly + cfg.cavity_depth) / m.ny;
    m.dz = cfg.Lz / m.nz;

    const double xs = std::min(cfg.cavity_x_start, cfg.cavity_x_end);
    const double xe = std::max(cfg.cavity_x_start, cfg.cavity_x_end);
    const double D = std::max(cfg.cavity_depth, 1.0e-8);
    if (!(xs > 0.0) || !(xe < cfg.Lx) || !(xe > xs))
        throw std::runtime_error("Cavity_X_Start/End must satisfy 0 < start < end < Lx");

    const int n_bay = std::max(12, m.ny / 4);
    const int n_up = m.ny - n_bay;
    if (n_up < 8)
        throw std::runtime_error("Cavity ny too small for bay+upper blocks");

    const int npx = m.nx + 1, npy = m.ny + 1, npz = m.nz + 1;
    std::vector<double> xnode;
    make_x_nodes(npx, cfg.Lx, xs, xe, std::max(cfg.cavity_x_cluster, 1.0), xnode);

    std::vector<double> ynode(static_cast<std::size_t>(npy), 0.0);
    for (int j = 0; j <= n_bay; ++j)
    {
        const double u = static_cast<double>(j) / n_bay;
        ynode[j] = -D + D * 0.5 * (1.0 - std::cos(3.141592653589793 * u));
    }
    ynode[n_bay] = 0.0;
    const double p_y = std::max(cfg.y_stretch, 1.0);
    for (int j = 1; j <= n_up; ++j)
    {
        const double u = static_cast<double>(j) / n_up;
        ynode[n_bay + j] = cfg.Ly * std::pow(u, 1.0 / p_y);
    }
    ynode[npy - 1] = cfg.Ly;

    m.has_nodes = true;
    m.nodes.assign(static_cast<std::size_t>(npx) * npy * npz * 3, 0.0);
    for (int k = 0; k < npz; ++k)
        for (int j = 0; j < npy; ++j)
            for (int i = 0; i < npx; ++i)
            {
                const int id = nid(i, j, k, npx, npy);
                m.nodes[id * 3 + 0] = xnode[i];
                m.nodes[id * 3 + 1] = ynode[j];
                m.nodes[id * 3 + 2] = cfg.Lz * (static_cast<double>(k) / m.nz);
            }

    m.xc.resize(m.n_cells);
    m.yc.resize(m.n_cells);
    m.zc.resize(m.n_cells);
    m.volume.assign(m.n_cells, 0.0);
    m.neigh.assign(m.n_cells * kFaces, -1);
    m.nxyz.assign(m.n_cells * kFaces * 3, 0.0);
    m.area.assign(m.n_cells * kFaces, 0.0);
    m.wall_face.assign(m.n_cells * kFaces, 0);
    m.bc_type.assign(m.n_cells * kFaces, BC_INTERIOR);
    m.active.assign(m.n_cells, 1);

    const int bc_face[kFaces] = {
        parse_bc_local(cfg.bc_xmin), parse_bc_local(cfg.bc_xmax), parse_bc_local(cfg.bc_ymin),
        parse_bc_local(cfg.bc_ymax), parse_bc_local(cfg.bc_zmin), parse_bc_local(cfg.bc_zmax)};
    const double nrm[kFaces][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                                   {0, 1, 0},  {0, 0, -1}, {0, 0, 1}};

    int n_fluid = 0;
    for (int k = 0; k < m.nz; ++k)
        for (int j = 0; j < m.ny; ++j)
            for (int i = 0; i < m.nx; ++i)
            {
                const int c = cid(i, j, k, m.nx, m.ny);
                const double x0 = xnode[i], x1 = xnode[i + 1];
                const double y0 = ynode[j], y1 = ynode[j + 1];
                const double z0 = cfg.Lz * k / m.nz;
                const double z1 = cfg.Lz * (k + 1) / m.nz;
                m.xc[c] = 0.5 * (x0 + x1);
                m.yc[c] = 0.5 * (y0 + y1);
                m.zc[c] = 0.5 * (z0 + z1);
                m.volume[c] = std::max((x1 - x0) * (y1 - y0) * (z1 - z0), 1.0e-30);

                const bool in_bay_x = (x0 >= xs - 1.0e-12 && x1 <= xe + 1.0e-12);
                const bool below = (y1 <= 1.0e-12);
                m.active[c] = (below && !in_bay_x) ? 0 : 1;
                if (m.active[c])
                    ++n_fluid;

                const double face_area[kFaces] = {(y1 - y0) * (z1 - z0), (y1 - y0) * (z1 - z0),
                                                  (x1 - x0) * (z1 - z0), (x1 - x0) * (z1 - z0),
                                                  (x1 - x0) * (y1 - y0), (x1 - x0) * (y1 - y0)};
                const int ni[kFaces][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k},
                                           {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
                for (int f = 0; f < kFaces; ++f)
                {
                    m.area[c * kFaces + f] = face_area[f];
                    m.nxyz[(c * kFaces + f) * 3 + 0] = nrm[f][0];
                    m.nxyz[(c * kFaces + f) * 3 + 1] = nrm[f][1];
                    m.nxyz[(c * kFaces + f) * 3 + 2] = nrm[f][2];
                    const int ii = ni[f][0], jj = ni[f][1], kk = ni[f][2];
                    const bool inb = (ii >= 0 && ii < m.nx && jj >= 0 && jj < m.ny && kk >= 0 &&
                                      kk < m.nz);
                    int nb = inb ? cid(ii, jj, kk, m.nx, m.ny) : -1;
                    if (nb >= 0 && m.active[nb] == 0)
                        nb = -1;
                    if (!m.active[c])
                    {
                        m.neigh[c * kFaces + f] = -1;
                        m.wall_face[c * kFaces + f] = 1;
                        m.bc_type[c * kFaces + f] = BC_NOSLIP;
                        continue;
                    }
                    if (nb >= 0)
                    {
                        m.neigh[c * kFaces + f] = nb;
                        continue;
                    }
                    m.neigh[c * kFaces + f] = -1;
                    m.wall_face[c * kFaces + f] = 1;
                    int bc = BC_NOSLIP;
                    if (!inb)
                    {
                        if (f == 0)
                            bc = bc_face[0];
                        else if (f == 1)
                            bc = bc_face[1];
                        else if (f == 3)
                            bc = bc_face[3];
                        else if (f == 4)
                            bc = bc_face[4];
                        else if (f == 5)
                            bc = bc_face[5];
                        else
                            bc = BC_NOSLIP; /* ymin of bay or plate */
                    }
                    m.bc_type[c * kFaces + f] = bc;
                }
            }

    /* Two-pass: neighbour active flags are complete after the first fill. */
    for (int c = 0; c < m.n_cells; ++c)
    {
        if (!m.active[c])
            continue;
        for (int f = 0; f < kFaces; ++f)
        {
            int nb = m.neigh[c * kFaces + f];
            if (nb >= 0 && !m.active[nb])
            {
                m.neigh[c * kFaces + f] = -1;
                m.wall_face[c * kFaces + f] = 1;
                m.bc_type[c * kFaces + f] = BC_NOSLIP;
            }
        }
    }

    std::cout << "Cavity 90° blocks: D=" << D << " L=" << (xe - xs) << " L/D=" << (xe - xs) / D
              << "  bay_j=" << n_bay << "  fluid=" << n_fluid << "/" << m.n_cells << "\n";
    return m;
}

} // namespace euler3d
