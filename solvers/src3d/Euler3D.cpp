/**
 * @file Euler3D.cpp
 * @brief Mesh, IC, time step, update, VTK, and host driver for full 3D solver.
 */
#include "Euler3D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace euler3d
{
namespace
{

inline int cid(int i, int j, int k, int nx, int ny)
{
    return i + j * nx + k * nx * ny;
}

int parse_bc(const std::string &s)
{
    if (s == "noslip" || s == "no_slip" || s == "wall_viscous")
        return BC_NOSLIP;
    if (s == "extrapolate" || s == "outlet" || s == "outflow")
        return BC_EXTRAP;
    if (s == "freestream" || s == "inlet" || s == "farfield")
        return BC_FREESTREAM;
    if (s == "postshock" || s == "shock" || s == "rh")
        return BC_POSTSHOCK;
    return BC_SLIP;
}

int parse_scheme(const std::string &s)
{
    if (s == "RICCA" || s == "ricca")
        return 1;
    if (s == "ROE" || s == "roe")
        return 2;
    if (s == "VANLEER" || s == "VanLeer" || s == "vanleer")
        return 3;
    return 0;
}

int parse_order(const std::string &s)
{
    if (s == "2O" || s == "MUSCL" || s == "2" || s == "Second")
        return 1;
    if (s == "WENO" || s == "WENO5" || s == "3" || s == "High")
        return 2;
    return 0;
}

int parse_ricca_sensor(const std::string &s)
{
    if (s == "rh" || s == "RH" || s == "ricca_rh" || s == "RICCA_RH" || s == "1")
        return 1;
    return 0; // legacy pressure
}

int parse_weno_hybrid(const std::string &s)
{
    if (s == "off" || s == "none" || s == "0" || s == "false")
        return 0;
    if (s == "pressure" || s == "old" || s == "1")
        return 1;
    return 2; // combo (paper default)
}

} // namespace

const char *scheme_name(int scheme)
{
    switch (scheme)
    {
    case 1:
        return "RICCA";
    case 2:
        return "ROE";
    case 3:
        return "VANLEER";
    default:
        return "LLF";
    }
}

const char *ricca_sensor_name(int s)
{
    return (s == 1) ? "RICCA-RH" : "legacy-pressure";
}

const char *weno_hybrid_name(int h)
{
    switch (h)
    {
    case 0:
        return "off";
    case 1:
        return "pressure";
    default:
        return "combo";
    }
}

const char *turb_name(int model)
{
    switch (model)
    {
    case 1:
        return "DNS";
    case 2:
        return "RANS-mixing";
    case 3:
        return "LES-Smagorinsky";
    case 4:
        return "LES-WALE";
    case 5:
        return "LES-Vreman";
    default:
        return "laminar";
    }
}

int parse_turb(const std::string &s)
{
    if (s == "dns" || s == "DNS")
        return 1;
    if (s == "rans" || s == "RANS" || s == "mixing" || s == "mixing_length")
        return 2;
    if (s == "les_smagorinsky" || s == "smagorinsky" || s == "LES_SMAG")
        return 3;
    if (s == "les_wale" || s == "wale" || s == "WALE")
        return 4;
    if (s == "les_vreman" || s == "vreman" || s == "VREMAN")
        return 5;
    if (s == "les" || s == "LES")
        return 4;
    return 0;
}

const char *order_name(int order)
{
    switch (order)
    {
    case 1:
        return "MUSCL";
    case 2:
        return "WENO5";
    default:
        return "1O";
    }
}

Config load_config(const std::string &json_path)
{
    Config cfg;
    std::ifstream in(json_path);
    if (!in)
        throw std::runtime_error("Cannot open config: " + json_path);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string s = ss.str();

    auto grab_int = [&](const char *key, int def) {
        const std::string k = std::string("\"") + key + "\"";
        auto pos = s.find(k);
        if (pos == std::string::npos)
            return def;
        pos = s.find(':', pos);
        if (pos == std::string::npos)
            return def;
        return std::atoi(s.c_str() + pos + 1);
    };
    auto grab_dbl = [&](const char *key, double def) {
        const std::string k = std::string("\"") + key + "\"";
        auto pos = s.find(k);
        if (pos == std::string::npos)
            return def;
        pos = s.find(':', pos);
        if (pos == std::string::npos)
            return def;
        return std::atof(s.c_str() + pos + 1);
    };
    auto grab_bool = [&](const char *key, bool def) {
        const std::string k = std::string("\"") + key + "\"";
        auto pos = s.find(k);
        if (pos == std::string::npos)
            return def;
        pos = s.find(':', pos);
        if (pos == std::string::npos)
            return def;
        const auto sub = s.substr(pos, 24);
        if (sub.find("true") != std::string::npos)
            return true;
        if (sub.find("false") != std::string::npos)
            return false;
        return def;
    };
    auto grab_str = [&](const char *key, const std::string &def) {
        const std::string k = std::string("\"") + key + "\"";
        auto pos = s.find(k);
        if (pos == std::string::npos)
            return def;
        pos = s.find(':', pos);
        if (pos == std::string::npos)
            return def;
        pos = s.find('"', pos + 1);
        if (pos == std::string::npos)
            return def;
        auto end = s.find('"', pos + 1);
        if (end == std::string::npos)
            return def;
        return s.substr(pos + 1, end - pos - 1);
    };

    cfg.nx = grab_int("nx", cfg.nx);
    cfg.ny = grab_int("ny", cfg.ny);
    cfg.nz = grab_int("nz", cfg.nz);
    cfg.Lx = grab_dbl("Lx", cfg.Lx);
    cfg.Ly = grab_dbl("Ly", cfg.Ly);
    cfg.Lz = grab_dbl("Lz", cfg.Lz);
    cfg.cfl = grab_dbl("CFL", cfg.cfl);
    cfg.max_iter = grab_int("Total_Iterations", cfg.max_iter);
    cfg.print_every = grab_int("Print_Every", cfg.print_every);
    cfg.use_cuda = grab_bool("Use_CUDA", cfg.use_cuda);
    cfg.viscous = grab_bool("Viscous", cfg.viscous);
    cfg.mu = grab_dbl("Mu", cfg.mu);
    cfg.Pr = grab_dbl("Pr", cfg.Pr);
    cfg.turb_model = parse_turb(grab_str("Turbulence", "laminar"));
    cfg.Cs = grab_dbl("Cs", cfg.Cs);
    cfg.Cw = grab_dbl("Cw", cfg.Cw);
    cfg.Cv = grab_dbl("Cv", cfg.Cv);
    cfg.Pr_t = grab_dbl("Pr_t", cfg.Pr_t);
    cfg.les_perturb = grab_dbl("LES_Perturb", cfg.les_perturb);
    if (cfg.turb_model >= 1)
        cfg.viscous = true;
    cfg.freestream_rho = grab_dbl("Freestream_Rho", cfg.freestream_rho);
    cfg.freestream_u = grab_dbl("Freestream_U", cfg.freestream_u);
    cfg.freestream_v = grab_dbl("Freestream_V", cfg.freestream_v);
    cfg.freestream_w = grab_dbl("Freestream_W", cfg.freestream_w);
    cfg.freestream_p = grab_dbl("Freestream_P", cfg.freestream_p);
    cfg.ramp_angle_deg = grab_dbl("Ramp_Angle", cfg.ramp_angle_deg);
    cfg.ramp_x_start = grab_dbl("Ramp_X_Start", cfg.ramp_x_start);
    cfg.ramp_x_end = grab_dbl("Ramp_X_End", cfg.ramp_x_end);
    cfg.shock_theta_deg = grab_dbl("Shock_Theta", cfg.shock_theta_deg);
    cfg.shock_x_top = grab_dbl("Shock_X_Top", cfg.shock_x_top);
    cfg.y_stretch = grab_dbl("Y_Stretch", cfg.y_stretch);
    cfg.cavity_x_start = grab_dbl("Cavity_X_Start", cfg.cavity_x_start);
    cfg.cavity_x_end = grab_dbl("Cavity_X_End", cfg.cavity_x_end);
    cfg.cavity_depth = grab_dbl("Cavity_Depth", cfg.cavity_depth);
    cfg.cavity_lip_width = grab_dbl("Cavity_Lip_Width", cfg.cavity_lip_width);
    cfg.cavity_x_cluster = grab_dbl("Cavity_X_Cluster", cfg.cavity_x_cluster);
    cfg.case_name = grab_str("Case", cfg.case_name);
    cfg.vtk_out = grab_str("VTK_Out", cfg.vtk_out);
    cfg.restart = grab_str("Restart", cfg.restart);
    cfg.bc_xmin = grab_str("BC_Xmin", cfg.bc_xmin);
    cfg.bc_xmax = grab_str("BC_Xmax", cfg.bc_xmax);
    cfg.bc_ymin = grab_str("BC_Ymin", cfg.bc_ymin);
    cfg.bc_ymax = grab_str("BC_Ymax", cfg.bc_ymax);
    cfg.bc_zmin = grab_str("BC_Zmin", cfg.bc_zmin);
    cfg.bc_zmax = grab_str("BC_Zmax", cfg.bc_zmax);
    cfg.scheme = parse_scheme(grab_str("Scheme", "LLF"));
    cfg.order = parse_order(grab_str("Order", "1O"));
    cfg.weno_char = grab_bool("Is_Char", cfg.weno_char);
    {
        const std::string wt = grab_str("WENO_Type", cfg.weno_z ? "Z" : "JS");
        cfg.weno_z = !(wt == "JS" || wt == "JiangShu" || wt == "0");
    }
    cfg.ricca_sensor = parse_ricca_sensor(grab_str("RICCA_Sensor", "legacy"));
    cfg.ricca_rh_threshold = grab_dbl("RICCA_RH_Threshold", cfg.ricca_rh_threshold);
    cfg.weno_hybrid = parse_weno_hybrid(grab_str("WENO_Hybrid", "combo"));
    fill_oblique_postshock(cfg);
    return cfg;
}

Mesh make_cartesian_mesh(const Config &cfg)
{
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

    m.xc.resize(m.n_cells);
    m.yc.resize(m.n_cells);
    m.zc.resize(m.n_cells);
    m.volume.assign(m.n_cells, m.dx * m.dy * m.dz);
    m.neigh.assign(m.n_cells * kFaces, -1);
    m.nxyz.assign(m.n_cells * kFaces * 3, 0.0);
    m.area.assign(m.n_cells * kFaces, 0.0);
    m.wall_face.assign(m.n_cells * kFaces, 0);
    m.bc_type.assign(m.n_cells * kFaces, BC_INTERIOR);

    const int bc_face[kFaces] = {parse_bc(cfg.bc_xmin), parse_bc(cfg.bc_xmax),
                                 parse_bc(cfg.bc_ymin), parse_bc(cfg.bc_ymax),
                                 parse_bc(cfg.bc_zmin), parse_bc(cfg.bc_zmax)};

    const double face_area[kFaces] = {m.dy * m.dz, m.dy * m.dz, m.dx * m.dz,
                                      m.dx * m.dz, m.dx * m.dy, m.dx * m.dy};
    const double nrm[kFaces][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                                   {0, 1, 0},  {0, 0, -1}, {0, 0, 1}};

    for (int k = 0; k < m.nz; ++k)
        for (int j = 0; j < m.ny; ++j)
            for (int i = 0; i < m.nx; ++i)
            {
                const int c = cid(i, j, k, m.nx, m.ny);
                m.xc[c] = (i + 0.5) * m.dx;
                m.yc[c] = (j + 0.5) * m.dy;
                m.zc[c] = (k + 0.5) * m.dz;

                const int ni[kFaces][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k},
                                           {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};

                for (int f = 0; f < kFaces; ++f)
                {
                    m.area[c * kFaces + f] = face_area[f];
                    m.nxyz[(c * kFaces + f) * 3 + 0] = nrm[f][0];
                    m.nxyz[(c * kFaces + f) * 3 + 1] = nrm[f][1];
                    m.nxyz[(c * kFaces + f) * 3 + 2] = nrm[f][2];

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

Mesh make_mesh_for_case(const Config &cfg)
{
    Mesh m;
    if (cfg.case_name == "ramp_15" || cfg.case_name == "ramp15" ||
        cfg.case_name == "Ramp_15_Degree")
        m = make_ramp15_mesh(cfg);
    else if (cfg.case_name == "swbli_plate" || cfg.case_name == "SWBLI" ||
             cfg.case_name == "flat_plate_swbli")
        m = make_swbli_plate_mesh(cfg);
    else if (cfg.case_name == "cavity" || cfg.case_name == "open_cavity" ||
             cfg.case_name == "supersonic_cavity")
        m = make_cavity_mesh(cfg);
    else
        m = make_cartesian_mesh(cfg);
    m.wall_dist = m.yc;
    if (cfg.case_name == "cavity" || cfg.case_name == "open_cavity" ||
        cfg.case_name == "supersonic_cavity")
    {
        const double xs = std::min(cfg.cavity_x_start, cfg.cavity_x_end);
        const double xe = std::max(cfg.cavity_x_start, cfg.cavity_x_end);
        const double D = std::max(cfg.cavity_depth, 1.0e-8);
        for (int c = 0; c < m.n_cells; ++c)
        {
            if (!cell_active(m, c))
            {
                m.wall_dist[c] = 0.0;
                continue;
            }
            const double x = m.xc[c];
            const double y = m.yc[c];
            if (y >= 0.0 && (x < xs || x > xe))
                m.wall_dist[c] = y; /* plate */
            else
            {
                double d = y + D; /* floor */
                d = std::min(d, std::max(x - xs, 0.0));
                d = std::min(d, std::max(xe - x, 0.0));
                if (y >= 0.0)
                    d = std::min(d, y); /* mouth */
                m.wall_dist[c] = std::max(d, 0.0);
            }
        }
    }
    return m;
}

void init_case(State &st, const Mesh &m, const Config &cfg)
{
    if (cfg.case_name == "sod_x")
        init_sod_x(st, m);
    else if (cfg.case_name == "cavity" || cfg.case_name == "open_cavity" ||
             cfg.case_name == "supersonic_cavity")
    {
        init_freestream(st, m, cfg);
        for (int c = 0; c < m.n_cells; ++c)
        {
            if (m.yc[c] >= 0.0)
                continue;
            U_from_prim_vals(cfg.freestream_rho, 0.0, 0.0, 0.0, cfg.freestream_p,
                             &st.U[c * kNv]);
        }
        prim_from_U(st, m);
    }
    else
        init_freestream(st, m, cfg);
}

void init_sod_x(State &st, const Mesh &m)
{
    st.U.assign(m.n_cells * kNv, 0.0);
    st.prim.assign(m.n_cells * 6, 0.0);
    for (int c = 0; c < m.n_cells; ++c)
    {
        const bool left = m.xc[c] < 0.5 * m.Lx;
        const double rho = left ? 1.0 : 0.125;
        const double p = left ? 1.0 : 0.1;
        U_from_prim_vals(rho, 0.0, 0.0, 0.0, p, &st.U[c * kNv]);
    }
    prim_from_U(st, m);
}

void init_freestream(State &st, const Mesh &m, const Config &cfg)
{
    st.U.assign(m.n_cells * kNv, 0.0);
    st.prim.assign(m.n_cells * 6, 0.0);
    for (int c = 0; c < m.n_cells; ++c)
        U_from_prim_vals(cfg.freestream_rho, cfg.freestream_u, cfg.freestream_v,
                         cfg.freestream_w, cfg.freestream_p, &st.U[c * kNv]);
    prim_from_U(st, m);
}

void prim_from_U(State &st, const Mesh &m)
{
    for (int c = 0; c < m.n_cells; ++c)
        fill_prim_cell(&st.U[c * kNv], &st.prim[c * 6]);
}

void apply_bc(State &st, const Mesh &m)
{
    (void)st;
    (void)m;
}

double compute_dt(const State &st, const Mesh &m, const Config &cfg)
{
    double dt_min = 1e300;
#ifdef _OPENMP
#pragma omp parallel for reduction(min : dt_min) schedule(static)
#endif
    for (int c = 0; c < m.n_cells; ++c)
    {
        if (!cell_active(m, c))
            continue;
        const double *p = &st.prim[c * 6];
        double denom = 0.0;
        for (int f = 0; f < kFaces; ++f)
        {
            const double nx = m.nxyz[(c * kFaces + f) * 3 + 0];
            const double ny = m.nxyz[(c * kFaces + f) * 3 + 1];
            const double nz = m.nxyz[(c * kFaces + f) * 3 + 2];
            const double a = m.area[c * kFaces + f];
            const double un = p[1] * nx + p[2] * ny + p[3] * nz;
            denom += (std::fabs(un) + p[5]) * a;
        }
        double dt = (denom > 1e-30) ? cfg.cfl * m.volume[c] / denom : 1e300;
        if (cfg.viscous && cfg.mu > 0.0)
        {
            double hx = m.dx, hy = m.dy, hz = m.dz;
            const int i = c % m.nx;
            const int k = c / (m.nx * m.ny);
            const int j = (c / m.nx) % m.ny;
            if (i + 1 < m.nx)
                hx = std::fabs(m.xc[c + 1] - m.xc[c]);
            else if (i > 0)
                hx = std::fabs(m.xc[c] - m.xc[c - 1]);
            if (j + 1 < m.ny)
                hy = std::fabs(m.yc[c + m.nx] - m.yc[c]);
            else if (j > 0)
                hy = std::fabs(m.yc[c] - m.yc[c - m.nx]);
            if (k + 1 < m.nz)
                hz = std::fabs(m.zc[c + m.nx * m.ny] - m.zc[c]);
            else if (k > 0)
                hz = std::fabs(m.zc[c] - m.zc[c - m.nx * m.ny]);
            const double h2 = std::min({hx * hx, hy * hy, hz * hz});
            const double nu = cfg.mu / std::max(p[0], 1e-14);
            dt = std::min(dt, 0.25 * cfg.cfl * h2 / std::max(nu, 1e-30));
        }
        dt_min = std::min(dt_min, dt);
    }
    return dt_min;
}

void update_explicit(State &st, const Mesh &m, const std::vector<double> &R, double dt,
                     double err[5])
{
    double err_acc[5] = {0, 0, 0, 0, 0};
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        double local_err[5] = {0, 0, 0, 0, 0};
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int c = 0; c < m.n_cells; ++c)
        {
            if (!cell_active(m, c))
                continue;
            double Uold[kNv];
            for (int v = 0; v < kNv; ++v)
                Uold[v] = st.U[c * kNv + v];
            for (int v = 0; v < kNv; ++v)
            {
                const double dU = -dt / m.volume[c] * R[c * kNv + v];
                const double denom = std::fabs(Uold[v]) < 1e-12 ? 1.0 : std::fabs(Uold[v]);
                local_err[v] += (dU / denom) * (dU / denom);
                st.U[c * kNv + v] = Uold[v] + dU;
            }
            double prim[6];
            fill_prim_cell(&st.U[c * kNv], prim);
            if (!(prim[0] > 1e-12) || !(prim[4] > 1e-12))
            {
                for (int v = 0; v < kNv; ++v)
                    st.U[c * kNv + v] = Uold[v];
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (int v = 0; v < kNv; ++v)
                err_acc[v] += local_err[v];
        }
    }
    for (int v = 0; v < kNv; ++v)
        err[v] = std::sqrt(err_acc[v]);
    prim_from_U(st, m);
}

bool write_vtk(const std::string &path, const Mesh &m, const State &st)
{
    std::ofstream out(path);
    if (!out)
        return false;

    if (m.has_nodes && static_cast<int>(m.nodes.size()) == (m.nx + 1) * (m.ny + 1) * (m.nz + 1) * 3)
    {
        const int npx = m.nx + 1, npy = m.ny + 1, npz = m.nz + 1;
        const int n_nodes = npx * npy * npz;
        out << "# vtk DataFile Version 3.0\nEuler3D\nASCII\nDATASET STRUCTURED_GRID\n";
        out << "DIMENSIONS " << npx << " " << npy << " " << npz << "\n";
        out << "POINTS " << n_nodes << " double\n";
        for (int n = 0; n < n_nodes; ++n)
            out << m.nodes[n * 3 + 0] << " " << m.nodes[n * 3 + 1] << " " << m.nodes[n * 3 + 2]
                << "\n";
        out << "CELL_DATA " << m.n_cells << "\n";
    }
    else
    {
        const int npx = m.nx + 1, npy = m.ny + 1, npz = m.nz + 1;
        out << "# vtk DataFile Version 3.0\nEuler3D\nASCII\nDATASET STRUCTURED_POINTS\n";
        out << "DIMENSIONS " << npx << " " << npy << " " << npz << "\n";
        out << "ORIGIN 0 0 0\n";
        out << "SPACING " << m.dx << " " << m.dy << " " << m.dz << "\n";
        out << "CELL_DATA " << m.n_cells << "\n";
    }

    out << "SCALARS Density double 1\nLOOKUP_TABLE default\n";
    for (int c = 0; c < m.n_cells; ++c)
        out << st.prim[c * 6 + 0] << "\n";
    out << "SCALARS Pressure double 1\nLOOKUP_TABLE default\n";
    for (int c = 0; c < m.n_cells; ++c)
        out << st.prim[c * 6 + 4] << "\n";
    out << "SCALARS Mach double 1\nLOOKUP_TABLE default\n";
    for (int c = 0; c < m.n_cells; ++c)
    {
        const double u = st.prim[c * 6 + 1];
        const double v = st.prim[c * 6 + 2];
        const double w = st.prim[c * 6 + 3];
        const double a = std::max(st.prim[c * 6 + 5], 1e-14);
        out << std::sqrt(u * u + v * v + w * w) / a << "\n";
    }
    out << "VECTORS Velocity double\n";
    for (int c = 0; c < m.n_cells; ++c)
        out << st.prim[c * 6 + 1] << " " << st.prim[c * 6 + 2] << " " << st.prim[c * 6 + 3]
            << "\n";
    return true;
}

bool read_vtk_restart(const std::string &path, const Mesh &m, State &st)
{
    std::ifstream in(path);
    if (!in)
        return false;

    std::string line;
    int n_cells = -1;
    while (std::getline(in, line))
    {
        std::istringstream hs(line);
        std::string tag;
        hs >> tag;
        if (tag == "CELL_DATA")
        {
            hs >> n_cells;
            break;
        }
    }
    if (n_cells <= 0)
        return false;
    const int nxy = m.nx * m.ny;
    if (nxy <= 0 || n_cells % nxy != 0 || m.n_cells % nxy != 0)
        return false;
    const int nz_src = n_cells / nxy;
    const int nz_dst = m.n_cells / nxy;
    if (nz_dst != m.nz)
        return false;

    st.U.assign(static_cast<std::size_t>(m.n_cells) * kNv, 0.0);
    st.prim.assign(static_cast<std::size_t>(m.n_cells) * 6, 0.0);

    auto skip_to = [&](const std::string &key) {
        while (std::getline(in, line))
        {
            if (line.compare(0, key.size(), key) == 0)
                return true;
        }
        return false;
    };
    auto read_n = [&](std::vector<double> &dst, int n) {
        dst.assign(static_cast<std::size_t>(n), 0.0);
        int got = 0;
        while (got < n && in >> dst[static_cast<std::size_t>(got)])
            ++got;
        return got == n;
    };

    if (!skip_to("SCALARS Density"))
        return false;
    std::getline(in, line);
    std::vector<double> rho, pressure, u, v, w;
    if (!read_n(rho, n_cells))
        return false;
    if (!skip_to("SCALARS Pressure"))
        return false;
    std::getline(in, line);
    if (!read_n(pressure, n_cells))
        return false;
    if (!skip_to("VECTORS Velocity"))
        return false;
    u.assign(n_cells, 0.0);
    v.assign(n_cells, 0.0);
    w.assign(n_cells, 0.0);
    for (int c = 0; c < n_cells; ++c)
    {
        if (!(in >> u[c] >> v[c] >> w[c]))
            return false;
    }
    for (int c = 0; c < m.n_cells; ++c)
    {
        const int k = c / nxy;
        const int ij = c - k * nxy;
        const int k_src = k % nz_src;
        const int src = ij + k_src * nxy;
        U_from_prim_vals(rho[src], u[src], v[src], w[src], pressure[src], &st.U[c * kNv]);
    }
    prim_from_U(st, m);
    if (nz_src != nz_dst)
        std::cout << "Restart: tiled z  " << nz_src << " -> " << nz_dst << "\n";
    return true;
}

void apply_les_trip(State &st, const Mesh &m, const Config &cfg)
{
    if (cfg.turb_model < 3 || cfg.les_perturb <= 0.0)
        return;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    const double amp = cfg.les_perturb * std::max(std::fabs(cfg.freestream_u), 1.0e-6);
    for (int c = 0; c < m.n_cells; ++c)
    {
        if (!cell_active(m, c))
            continue;
        double q[6];
        fill_prim_cell(&st.U[c * kNv], q);
        const double wall = (c < static_cast<int>(m.wall_dist.size())) ? m.wall_dist[c] : m.yc[c];
        const double damp = std::min(1.0, wall / 0.02);
        const double du = amp * damp * dist(rng);
        const double dv = amp * damp * dist(rng);
        const double dw = amp * damp * dist(rng);
        U_from_prim_vals(q[0], q[1] + du, q[2] + dv, q[3] + dw, q[4], &st.U[c * kNv]);
    }
    prim_from_U(st, m);
    std::cout << "LES trip: " << cfg.les_perturb << " |Uinf| (damped at wall)\n";
}

int run_host(const Config &cfg)
{
    std::cout << "=== 3D HOST " << scheme_name(cfg.scheme) << "/" << order_name(cfg.order)
              << (cfg.order == 2 ? (cfg.weno_char ? "+Char" : "+Comp") : "")
              << (cfg.order == 2 ? (cfg.weno_z ? "+WENOZ" : "+WENOJS") : "")
              << (cfg.order == 2
                      ? (std::string(" hybrid=") + weno_hybrid_name(cfg.weno_hybrid))
                      : "")
              << (cfg.scheme == 1
                      ? (std::string(" ricca=") + ricca_sensor_name(cfg.ricca_sensor))
                      : "")
              << (cfg.viscous ? " +NS" : "") << "  turb=" << turb_name(cfg.turb_model)
              << "  " << cfg.nx << "x" << cfg.ny << "x" << cfg.nz
              << "  case=" << cfg.case_name << " ===\n";
    Mesh m = make_mesh_for_case(cfg);
    State st;
    init_case(st, m, cfg);
    if (!cfg.restart.empty())
    {
        if (!read_vtk_restart(cfg.restart, m, st))
            throw std::runtime_error("Failed to restart from " + cfg.restart);
        std::cout << "Restart from " << cfg.restart << "\n";
    }
    apply_les_trip(st, m, cfg);

    std::vector<double> R;
    double err[5];
    for (int it = 1; it <= cfg.max_iter; ++it)
    {
        const double dt = compute_dt(st, m, cfg);
        net_flux_host(st, m, cfg, R);
        add_viscous_flux_host(st, m, cfg, R);
        update_explicit(st, m, R, dt, err);
        if (it % cfg.print_every == 0 || it == 1 || it == cfg.max_iter)
        {
            double pmin = 1e300, pmax = -1e300, mmin = 1e300, mmax = -1e300;
            for (int c = 0; c < m.n_cells; ++c)
            {
                pmin = std::min(pmin, st.prim[c * 6 + 4]);
                pmax = std::max(pmax, st.prim[c * 6 + 4]);
                const double spd = std::sqrt(st.prim[c * 6 + 1] * st.prim[c * 6 + 1] +
                                            st.prim[c * 6 + 2] * st.prim[c * 6 + 2] +
                                            st.prim[c * 6 + 3] * st.prim[c * 6 + 3]);
                const double mach = spd / std::max(st.prim[c * 6 + 5], 1e-14);
                mmin = std::min(mmin, mach);
                mmax = std::max(mmax, mach);
            }
            std::printf("%8d  dt=%.3e  err(rho)=%.3e  p[%.3e,%.3e]  M[%.3f,%.3f]\n", it, dt,
                        err[0], pmin, pmax, mmin, mmax);
            std::fflush(stdout);
        }
    }
    if (!write_vtk(cfg.vtk_out, m, st))
        std::cerr << "Warning: failed to write VTK " << cfg.vtk_out << "\n";
    else
        std::cout << "Wrote " << cfg.vtk_out << "\n";
    return 0;
}

} // namespace euler3d
