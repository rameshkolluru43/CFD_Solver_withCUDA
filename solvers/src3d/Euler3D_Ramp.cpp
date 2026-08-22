/**
 * @file Euler3D_Ramp.cpp
 * @brief Body-fitted 15° compression-ramp mesh extruded in z (matches 2D ramp layout).
 *
 * Domain: x∈[0,Lx], wall y_w(x), top y=Ly, span z∈[0,Lz].
 * Flat floor to x_ramp, then 15° ramp to x_plateau, then flat plateau.
 */
#include "Euler3D.hpp"

#include <algorithm>
#include <cmath>
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

inline void cross(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

inline double norm3(const double v[3])
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline void sub3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

/** Quad area and unit outward normal from four corners (ccw when viewed from outside). */
void quad_area_normal(const double p0[3], const double p1[3], const double p2[3],
                      const double p3[3], double n_unit[3], double &area)
{
    double d1[3], d2[3], c1[3], c2[3];
    sub3(p2, p0, d1);
    sub3(p3, p1, d2);
    cross(d1, d2, c1);
    /* alternate triangulation average */
    sub3(p1, p0, d1);
    sub3(p2, p0, d2);
    cross(d1, d2, c1);
    sub3(p2, p0, d1);
    sub3(p3, p0, d2);
    cross(d1, d2, c2);
    double n[3] = {0.5 * (c1[0] + c2[0]), 0.5 * (c1[1] + c2[1]), 0.5 * (c1[2] + c2[2])};
    area = norm3(n);
    if (area < 1e-30)
    {
        n_unit[0] = n_unit[1] = n_unit[2] = 0.0;
        area = 0.0;
        return;
    }
    n_unit[0] = n[0] / area;
    n_unit[1] = n[1] / area;
    n_unit[2] = n[2] / area;
}

double hex_volume(const double p[8][3])
{
    /* Decompose into 5 or 6 tets from vertex 0 — use average of two decompositions. */
    auto tet = [](const double a[3], const double b[3], const double c[3], const double d[3]) {
        double ab[3], ac[3], ad[3], cr[3];
        sub3(b, a, ab);
        sub3(c, a, ac);
        sub3(d, a, ad);
        cross(ab, ac, cr);
        return std::fabs(cr[0] * ad[0] + cr[1] * ad[1] + cr[2] * ad[2]) / 6.0;
    };
    const double v1 = tet(p[0], p[1], p[3], p[4]) + tet(p[1], p[2], p[3], p[6]) +
                      tet(p[1], p[3], p[4], p[6]) + tet(p[1], p[4], p[5], p[6]) +
                      tet(p[3], p[4], p[6], p[7]);
    return v1;
}

double y_wall(double x, double x_ramp, double x_plateau, double angle_deg)
{
    const double ang = angle_deg * 3.14159265358979323846 / 180.0;
    const double slope = std::tan(ang);
    if (x <= x_ramp)
        return 0.0;
    if (x <= x_plateau)
        return slope * (x - x_ramp);
    return slope * (x_plateau - x_ramp);
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

Mesh make_ramp15_mesh(const Config &cfg)
{
    if (cfg.nx < 4 || cfg.ny < 4 || cfg.nz < 1)
        throw std::runtime_error("Ramp mesh needs nx,ny >= 4 and nz >= 1");

    Mesh m;
    m.nx = cfg.nx;
    m.ny = cfg.ny;
    m.nz = cfg.nz;
    m.Lx = cfg.Lx;
    m.Ly = cfg.Ly;
    m.Lz = cfg.Lz;
    m.n_cells = m.nx * m.ny * m.nz;
    m.dx = cfg.Lx / m.nx;
    m.dy = cfg.Ly / m.ny;
    m.dz = cfg.Lz / m.nz;

    const int npx = m.nx + 1, npy = m.ny + 1, npz = m.nz + 1;
    const int n_nodes = npx * npy * npz;
    m.nodes.assign(static_cast<std::size_t>(n_nodes) * 3, 0.0);
    m.has_nodes = true;

    const double x_ramp = cfg.ramp_x_start;
    const double x_plateau = cfg.ramp_x_end;
    const double ang = cfg.ramp_angle_deg;

    for (int k = 0; k < npz; ++k)
        for (int j = 0; j < npy; ++j)
            for (int i = 0; i < npx; ++i)
            {
                const double x = cfg.Lx * (static_cast<double>(i) / m.nx);
                const double z = cfg.Lz * (static_cast<double>(k) / m.nz);
                const double yw = y_wall(x, x_ramp, x_plateau, ang);
                const double eta = static_cast<double>(j) / m.ny;
                const double y = yw + eta * (cfg.Ly - yw);
                const int id = nid(i, j, k, npx, npy);
                m.nodes[id * 3 + 0] = x;
                m.nodes[id * 3 + 1] = y;
                m.nodes[id * 3 + 2] = z;
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

    const int bc_face[kFaces] = {
        parse_bc_local(cfg.bc_xmin), parse_bc_local(cfg.bc_xmax), parse_bc_local(cfg.bc_ymin),
        parse_bc_local(cfg.bc_ymax), parse_bc_local(cfg.bc_zmin), parse_bc_local(cfg.bc_zmax)};

    auto node = [&](int i, int j, int k, double out[3]) {
        const int id = nid(i, j, k, npx, npy);
        out[0] = m.nodes[id * 3 + 0];
        out[1] = m.nodes[id * 3 + 1];
        out[2] = m.nodes[id * 3 + 2];
    };

    for (int k = 0; k < m.nz; ++k)
        for (int j = 0; j < m.ny; ++j)
            for (int i = 0; i < m.nx; ++i)
            {
                const int c = cid(i, j, k, m.nx, m.ny);
                double p[8][3];
                node(i, j, k, p[0]);
                node(i + 1, j, k, p[1]);
                node(i + 1, j + 1, k, p[2]);
                node(i, j + 1, k, p[3]);
                node(i, j, k + 1, p[4]);
                node(i + 1, j, k + 1, p[5]);
                node(i + 1, j + 1, k + 1, p[6]);
                node(i, j + 1, k + 1, p[7]);

                m.volume[c] = hex_volume(p);
                if (!(m.volume[c] > 1e-30))
                    throw std::runtime_error("Non-positive hex volume in ramp mesh");

                m.xc[c] = m.yc[c] = m.zc[c] = 0.0;
                for (int v = 0; v < 8; ++v)
                {
                    m.xc[c] += p[v][0];
                    m.yc[c] += p[v][1];
                    m.zc[c] += p[v][2];
                }
                m.xc[c] *= 0.125;
                m.yc[c] *= 0.125;
                m.zc[c] *= 0.125;

                /* Face corners ordered for outward normal (right-hand rule). */
                const int face_nodes[kFaces][4] = {
                    {0, 3, 7, 4}, /* -x */
                    {1, 2, 6, 5}, /* +x */
                    {0, 1, 5, 4}, /* -y */
                    {3, 2, 6, 7}, /* +y */
                    {0, 1, 2, 3}, /* -z */
                    {4, 5, 6, 7}  /* +z */
                };
                /* Fix -x / -y / -z winding if needed by flipping against cell center. */
                for (int f = 0; f < kFaces; ++f)
                {
                    double n[3], area;
                    quad_area_normal(p[face_nodes[f][0]], p[face_nodes[f][1]],
                                     p[face_nodes[f][2]], p[face_nodes[f][3]], n, area);
                    double mid[3] = {0, 0, 0};
                    for (int q = 0; q < 4; ++q)
                        for (int d = 0; d < 3; ++d)
                            mid[d] += p[face_nodes[f][q]][d];
                    mid[0] *= 0.25;
                    mid[1] *= 0.25;
                    mid[2] *= 0.25;
                    const double outward =
                        n[0] * (mid[0] - m.xc[c]) + n[1] * (mid[1] - m.yc[c]) +
                        n[2] * (mid[2] - m.zc[c]);
                    if (outward < 0.0)
                    {
                        n[0] = -n[0];
                        n[1] = -n[1];
                        n[2] = -n[2];
                    }
                    m.area[c * kFaces + f] = area;
                    m.nxyz[(c * kFaces + f) * 3 + 0] = n[0];
                    m.nxyz[(c * kFaces + f) * 3 + 1] = n[1];
                    m.nxyz[(c * kFaces + f) * 3 + 2] = n[2];
                }

                const int ni[kFaces][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k},
                                           {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
                for (int f = 0; f < kFaces; ++f)
                {
                    const int ii = ni[f][0], jj = ni[f][1], kk = ni[f][2];
                    if (ii >= 0 && ii < m.nx && jj >= 0 && jj < m.ny && kk >= 0 && kk < m.nz)
                        m.neigh[c * kFaces + f] = cid(ii, jj, kk, m.nx, m.ny);
                    else
                    {
                        m.wall_face[c * kFaces + f] = 1;
                        m.bc_type[c * kFaces + f] = bc_face[f];
                    }
                }
            }

    return m;
}

} // namespace euler3d
