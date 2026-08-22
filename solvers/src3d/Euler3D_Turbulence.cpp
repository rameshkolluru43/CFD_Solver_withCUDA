/**
 * @file Euler3D_Turbulence.cpp
 * @brief RANS mixing-length, LES SGS (Smagorinsky/WALE/Vreman), and DNS (mu_sgs=0).
 */
#include "Euler3D.hpp"

#include <algorithm>
#include <cmath>

namespace euler3d
{
namespace
{

inline double strain_rate_mag(const double g[12])
{
    const double sxx = g[0], syy = g[4], szz = g[8];
    const double sxy = 0.5 * (g[1] + g[3]);
    const double sxz = 0.5 * (g[2] + g[6]);
    const double syz = 0.5 * (g[5] + g[7]);
    const double twoSijSij =
        2.0 * (sxx * sxx + syy * syy + szz * szz) + 4.0 * (sxy * sxy + sxz * sxz + syz * syz);
    return std::sqrt(std::max(twoSijSij, 0.0));
}

inline double wale_sd_mag(const double g[12])
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
        for (int j = 0; j < 3; ++j)
        {
            double sd = 0.5 * (G2[i][j] + G2[j][i]);
            if (i == j)
                sd -= tr / 3.0;
            Sd2 += sd * sd;
        }
    return std::sqrt(std::max(Sd2, 0.0));
}

inline double vreman_op(const double g[12], double delta)
{
    const double a[3][3] = {{g[0], g[3], g[6]}, {g[1], g[4], g[7]}, {g[2], g[5], g[8]}};
    double alpha2 = 0.0;
    double beta[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            double s = 0.0;
            for (int m = 0; m < 3; ++m)
            {
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
    return std::sqrt(B / alpha2);
}

} // namespace

double mu_sgs_from_grad(int turb_model, double rho, const double g[12], double delta,
                        double y_wall, double mu_lam, const Config &cfg)
{
    if (turb_model <= 1)
        return 0.0;
    rho = std::max(rho, 1.0e-14);
    const double nu_lam = mu_lam / rho;
    const double S = strain_rate_mag(g);
    delta = std::max(delta, 1.0e-14);
    y_wall = std::max(y_wall, 0.0);

    if (turb_model == 2)
    {
        const double u_tau = 0.04 * std::max(std::fabs(cfg.freestream_u), 1.0e-6);
        const double yplus = y_wall * u_tau / std::max(nu_lam, 1.0e-14);
        const double damp = 1.0 - std::exp(-yplus / std::max(cfg.A_plus, 1.0));
        const double mix = cfg.kappa_von * y_wall * damp;
        return rho * mix * mix * S;
    }
    if (turb_model == 3)
    {
        const double u_tau = 0.04 * std::max(std::fabs(cfg.freestream_u), 1.0e-6);
        const double yplus = y_wall * u_tau / std::max(nu_lam, 1.0e-14);
        const double fvd = 1.0 - std::exp(-yplus / std::max(cfg.A_plus, 1.0));
        const double ls = cfg.Cs * delta * fvd;
        return rho * ls * ls * S;
    }
    if (turb_model == 4)
    {
        const double Sd = wale_sd_mag(g);
        const double num = std::pow(Sd, 1.5);
        const double den = std::pow(S, 2.5) + std::pow(Sd, 1.25);
        if (den < 1.0e-30)
            return 0.0;
        const double cw = cfg.Cw * delta;
        return rho * cw * cw * num / std::max(den, 1.0e-30);
    }
    if (turb_model == 5)
        return rho * cfg.Cv * vreman_op(g, delta);
    return 0.0;
}

} // namespace euler3d
