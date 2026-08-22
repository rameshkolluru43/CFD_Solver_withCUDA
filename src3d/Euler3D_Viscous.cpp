/**
 * @file Euler3D_Viscous.cpp
 * @brief Host laminar Navier--Stokes viscous fluxes for a Cartesian mesh.
 */
#include "Euler3D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace euler3d
{
namespace
{

constexpr int kNGrad = 4; // u, v, w, T
constexpr int kNDim = 3;

using Gradient = std::array<double, kNGrad * kNDim>;

inline double primitive_value(const State &st, int cell, int variable)
{
    if (variable < 3)
        return st.prim[static_cast<std::size_t>(cell) * 6 + variable + 1];

    const double rho =
        std::max(st.prim[static_cast<std::size_t>(cell) * 6], 1.0e-14);
    const double pressure =
        st.prim[static_cast<std::size_t>(cell) * 6 + 4];
    return pressure / rho; // R_gas = 1
}

inline double wall_value(int variable, double cell_value)
{
    // Stationary no-slip wall for velocity; adiabatic wall for temperature.
    return variable < 3 ? 0.0 : cell_value;
}

inline double neighbour_distance(const Mesh &m, int cell, int neighbour, int face,
                                 double fallback)
{
    if (neighbour < 0)
        return fallback;
    if (face < 2)
        return std::fabs(m.xc[neighbour] - m.xc[cell]);
    if (face < 4)
        return std::fabs(m.yc[neighbour] - m.yc[cell]);
    return std::fabs(m.zc[neighbour] - m.zc[cell]);
}

inline double directional_gradient(const State &st, const Mesh &m, int cell,
                                   int variable, int negative_face,
                                   int positive_face, double spacing)
{
    const int negative = m.neigh[cell * kFaces + negative_face];
    const int positive = m.neigh[cell * kFaces + positive_face];

    auto bc_of = [&](int face) {
        const int idx = cell * kFaces + face;
        return idx < static_cast<int>(m.bc_type.size()) ? m.bc_type[idx] : BC_EXTRAP;
    };

    if (negative >= 0 && positive >= 0)
    {
        const double dn = neighbour_distance(m, cell, negative, negative_face, spacing);
        const double dp = neighbour_distance(m, cell, positive, positive_face, spacing);
        return (primitive_value(st, positive, variable) -
                primitive_value(st, negative, variable)) /
               std::max(dn + dp, 1.0e-14);
    }

    const double center = primitive_value(st, cell, variable);
    if (negative < 0)
    {
        if (bc_of(negative_face) == BC_NOSLIP)
        {
            const double boundary = wall_value(variable, center);
            const double h = (positive >= 0)
                                 ? neighbour_distance(m, cell, positive, positive_face, spacing)
                                 : spacing;
            return (center - boundary) / std::max(0.5 * h, 1.0e-14);
        }
        if (positive >= 0)
            return (primitive_value(st, positive, variable) - center) /
                   std::max(neighbour_distance(m, cell, positive, positive_face, spacing),
                            1.0e-14);
        return 0.0;
    }

    if (bc_of(positive_face) == BC_NOSLIP)
    {
        const double boundary = wall_value(variable, center);
        const double h = neighbour_distance(m, cell, negative, negative_face, spacing);
        return (boundary - center) / std::max(0.5 * h, 1.0e-14);
    }
    return (center - primitive_value(st, negative, variable)) /
           std::max(neighbour_distance(m, cell, negative, negative_face, spacing), 1.0e-14);
}

std::vector<Gradient> compute_cell_gradients(const State &st, const Mesh &m)
{
    std::vector<Gradient> gradients(static_cast<std::size_t>(m.n_cells));

    for (int cell = 0; cell < m.n_cells; ++cell)
    {
        for (int variable = 0; variable < kNGrad; ++variable)
        {
            gradients[cell][variable * kNDim] =
                directional_gradient(st, m, cell, variable, 0, 1, m.dx);
            gradients[cell][variable * kNDim + 1] =
                directional_gradient(st, m, cell, variable, 2, 3, m.dy);
            gradients[cell][variable * kNDim + 2] =
                directional_gradient(st, m, cell, variable, 4, 5, m.dz);
        }
    }

    return gradients;
}

} // namespace

void add_viscous_flux_host(const State &st, const Mesh &m, const Config &cfg,
                           std::vector<double> &R)
{
    if (!cfg.viscous)
        return;

    const std::size_t residual_size =
        static_cast<std::size_t>(m.n_cells) * kNv;
    if (R.size() < residual_size)
        R.resize(residual_size, 0.0);

    const std::vector<Gradient> gradients = compute_cell_gradients(st, m);
    const double mu_lam = cfg.mu;
    const double cp = kGamma / (kGamma - 1.0);
    std::vector<double> mu_sgs(static_cast<std::size_t>(m.n_cells), 0.0);
    if (cfg.turb_model >= 2)
    {
        for (int cell = 0; cell < m.n_cells; ++cell)
        {
            double g[12];
            for (int i = 0; i < 12; ++i)
                g[i] = gradients[cell][static_cast<std::size_t>(i)];
            const double rho = std::max(st.prim[static_cast<std::size_t>(cell) * 6], 1.0e-14);
            const double delta = std::cbrt(std::max(m.volume[cell], 1.0e-30));
            const double yw =
                cell < static_cast<int>(m.wall_dist.size()) ? m.wall_dist[cell] : m.yc[cell];
            mu_sgs[cell] = mu_sgs_from_grad(cfg.turb_model, rho, g, delta, yw, mu_lam, cfg);
        }
    }

    for (int cell = 0; cell < m.n_cells; ++cell)
    {
        for (int face = 0; face < kFaces; ++face)
        {
            const int face_index = cell * kFaces + face;
            const int neighbour = m.neigh[face_index];
            const bool wall = neighbour < 0;

            Gradient face_gradient = gradients[cell];
            if (!wall)
            {
                for (int component = 0; component < kNGrad * kNDim; ++component)
                    face_gradient[component] =
                        0.5 * (gradients[cell][component] +
                               gradients[neighbour][component]);
            }

            const double ux = face_gradient[0], uy = face_gradient[1];
            const double uz = face_gradient[2];
            const double vx = face_gradient[3], vy = face_gradient[4];
            const double vz = face_gradient[5];
            const double wx = face_gradient[6], wy = face_gradient[7];
            const double wz = face_gradient[8];

            const double mu_t_f =
                wall ? mu_sgs[cell]
                     : 0.5 * (mu_sgs[cell] + mu_sgs[static_cast<std::size_t>(neighbour)]);
            const double mu = mu_lam + mu_t_f;
            const double thermal_conductivity =
                (cfg.Pr > 0.0 ? mu_lam * cp / cfg.Pr : 0.0) +
                (cfg.Pr_t > 0.0 ? mu_t_f * cp / cfg.Pr_t : 0.0);
            const double two_thirds_mu = (2.0 / 3.0) * mu;
            const double tau_xx = two_thirds_mu * (2.0 * ux - vy - wz);
            const double tau_yy = two_thirds_mu * (2.0 * vy - ux - wz);
            const double tau_zz = two_thirds_mu * (2.0 * wz - ux - vy);
            const double tau_xy = mu * (uy + vx);
            const double tau_xz = mu * (uz + wx);
            const double tau_yz = mu * (vz + wy);

            const double qx = -thermal_conductivity * face_gradient[9];
            const double qy = -thermal_conductivity * face_gradient[10];
            const double qz = -thermal_conductivity * face_gradient[11];

            double u = primitive_value(st, cell, 0);
            double v = primitive_value(st, cell, 1);
            double w = primitive_value(st, cell, 2);
            if (wall)
            {
                u = 0.0;
                v = 0.0;
                w = 0.0;
            }
            else
            {
                u = 0.5 * (u + primitive_value(st, neighbour, 0));
                v = 0.5 * (v + primitive_value(st, neighbour, 1));
                w = 0.5 * (w + primitive_value(st, neighbour, 2));
            }

            const double nx = m.nxyz[face_index * 3];
            const double ny = m.nxyz[face_index * 3 + 1];
            const double nz = m.nxyz[face_index * 3 + 2];
            const double area = m.area[face_index];

            const double traction_x =
                tau_xx * nx + tau_xy * ny + tau_xz * nz;
            const double traction_y =
                tau_xy * nx + tau_yy * ny + tau_yz * nz;
            const double traction_z =
                tau_xz * nx + tau_yz * ny + tau_zz * nz;
            const double energy_flux =
                u * traction_x + v * traction_y + w * traction_z -
                (qx * nx + qy * ny + qz * nz);

            // R stores outward inviscid flux.  Navier--Stokes contributes
            // minus the outward viscous flux.
            R[static_cast<std::size_t>(cell) * kNv + 1] -= traction_x * area;
            R[static_cast<std::size_t>(cell) * kNv + 2] -= traction_y * area;
            R[static_cast<std::size_t>(cell) * kNv + 3] -= traction_z * area;
            R[static_cast<std::size_t>(cell) * kNv + 4] -= energy_flux * area;
        }
    }
}

} // namespace euler3d
