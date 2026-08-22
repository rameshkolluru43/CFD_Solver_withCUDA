/**
 * @file Main_3D_CUDA.cu
 * @brief Resident-GPU driver for the 3D Euler/NS solver (Use_CUDA).
 */
#include "Euler3D.hpp"
#include "Euler3D_Cuda.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace
{

void ensure_parent_dirs(const std::string &path)
{
    const std::string::size_type pos = path.find_last_of('/');
    if (pos == std::string::npos)
        return;

    const std::string dir = path.substr(0, pos);
    std::string current;
    for (std::size_t i = 0; i < dir.size(); ++i)
    {
        current.push_back(dir[i]);
        if (dir[i] == '/' || i + 1 == dir.size())
        {
            if (!current.empty() && current != "/")
                mkdir(current.c_str(), 0755);
        }
    }
}

void print_progress(int iteration, double dt, const euler3d::State &state,
                    const euler3d::Mesh &mesh)
{
    double pmin = 1.0e300;
    double pmax = -1.0e300;
    for (int c = 0; c < mesh.n_cells; ++c)
    {
        if (!euler3d::cell_active(mesh, c))
            continue;
        const double p = state.prim[c * 6 + 4];
        pmin = std::min(pmin, p);
        pmax = std::max(pmax, p);
    }
    std::printf("%8d  dt=%.3e  p[min,max]=[%.3e,%.3e]\n",
                iteration, dt, pmin, pmax);
    std::fflush(stdout);
}

} // namespace

namespace euler3d
{

int run_cuda(const Config &cfg)
{
    std::cout << "=== 3D CUDA " << scheme_name(cfg.scheme) << "/" << order_name(cfg.order)
              << (cfg.order == 2 ? (cfg.weno_char ? "+Char" : "+Comp") : "")
              << (cfg.order == 2 ? (cfg.weno_z ? "+WENOZ" : "+WENOJS") : "")
              << (cfg.order == 2
                      ? (std::string(" hybrid=") + weno_hybrid_name(cfg.weno_hybrid))
                      : "")
              << (cfg.scheme == 1
                      ? (std::string(" ricca=") + ricca_sensor_name(cfg.ricca_sensor))
                      : "")
              << (cfg.viscous ? " +NS" : "") << "  turb=" << turb_name(cfg.turb_model)
              << "  " << cfg.nx << "x" << cfg.ny << "x"
              << cfg.nz << " (resident GPU) ===\n";
    std::fflush(stdout);

    Mesh mesh = make_mesh_for_case(cfg);
    State state;
    init_case(state, mesh, cfg);
    if (!cfg.restart.empty())
    {
        if (!read_vtk_restart(cfg.restart, mesh, state))
            throw std::runtime_error("Failed to restart from " + cfg.restart);
        std::cout << "Restart from " << cfg.restart << "\n";
        std::fflush(stdout);
    }
    apply_les_trip(state, mesh, cfg);

    euler3d_cuda::DeviceBundle device;
    euler3d_cuda::allocate(device, mesh.n_cells);
    euler3d_cuda::set_grid_size(device, mesh.nx, mesh.ny, mesh.nz);
    euler3d_cuda::upload_mesh(device, mesh.volume, mesh.neigh, mesh.nxyz, mesh.area,
                              mesh.bc_type, mesh.xc, mesh.yc, mesh.zc, mesh.wall_dist,
                              mesh.active, mesh.dx, mesh.dy, mesh.dz);
    euler3d_cuda::set_freestream(device, cfg.freestream_rho, cfg.freestream_u,
                                 cfg.freestream_v, cfg.freestream_w, cfg.freestream_p);
    euler3d_cuda::set_postshock(device, cfg.post_rho, cfg.post_u, cfg.post_v,
                                cfg.post_w, cfg.post_p);
    euler3d_cuda::upload_U(device, state.U);

    for (int iteration = 1; iteration <= cfg.max_iter; ++iteration)
    {
        const double dt = euler3d_cuda::launch_min_dt(device, cfg.cfl, cfg.viscous, cfg.mu);
        euler3d_cuda::launch_net_flux(device, cfg.scheme, cfg.order,
                                      cfg.weno_char, cfg.weno_z,
                                      cfg.ricca_sensor, cfg.ricca_rh_threshold,
                                      cfg.weno_hybrid);
        if (cfg.viscous)
            euler3d_cuda::launch_viscous_flux(device, cfg.mu, cfg.Pr, cfg.turb_model,
                                              cfg.Cs, cfg.Cw, cfg.Cv, cfg.Pr_t,
                                              cfg.kappa_von, cfg.A_plus,
                                              cfg.freestream_u);

        euler3d_cuda::launch_update(device, dt);
        const bool report = iteration == 1 || iteration == cfg.max_iter ||
                            (cfg.print_every > 0 && iteration % cfg.print_every == 0);
        if (report)
        {
            euler3d_cuda::download_U(device, state.U);
            prim_from_U(state, mesh);
            print_progress(iteration, dt, state, mesh);
        }
    }

    euler3d_cuda::download_U(device, state.U);
    prim_from_U(state, mesh);
    euler3d_cuda::release(device);

    if (!write_vtk(cfg.vtk_out, mesh, state))
        std::cerr << "Warning: failed to write VTK " << cfg.vtk_out << "\n";
    else
        std::cout << "Wrote " << cfg.vtk_out << "\n";
    return 0;
}

} // namespace euler3d

int main(int argc, char **argv)
{
    const std::string config_path =
        argc > 1 ? argv[1] : "input/json_Files/Run_3D_Sod_LLF_smoke_cuda.json";
    try
    {
        const euler3d::Config config = euler3d::load_config(config_path);
        ensure_parent_dirs(config.vtk_out);
        return euler3d::run_cuda(config);
    }
    catch (const std::exception &error)
    {
        std::cerr << "Fatal: " << error.what() << "\n";
        return 1;
    }
}
