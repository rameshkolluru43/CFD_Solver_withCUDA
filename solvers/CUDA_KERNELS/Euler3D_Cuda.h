#ifndef CFD_EULER3D_CUDA_H
#define CFD_EULER3D_CUDA_H

#include <vector>

namespace euler3d_cuda {

struct DeviceBundle {
    int n_cells = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    double *d_U = nullptr, *d_R = nullptr, *d_volume = nullptr;
    double *d_nxyz = nullptr, *d_area = nullptr;
    int *d_neigh = nullptr, *d_bc = nullptr;
    double *d_freestream = nullptr;
    double *d_postshock = nullptr;
    double *d_xc = nullptr, *d_yc = nullptr, *d_zc = nullptr;
    double *d_grad = nullptr; /* n_cells * 12: du/dx..dT/dz */
    double *d_wall_dist = nullptr;
    double *d_mu_sgs = nullptr;
    int *d_active = nullptr;
};

void allocate(DeviceBundle &dev, int n_cells);
void release(DeviceBundle &dev);
void set_grid_size(DeviceBundle &dev, int nx, int ny, int nz);
void upload_mesh(DeviceBundle &dev, const std::vector<double> &volume,
                 const std::vector<int> &neigh, const std::vector<double> &nxyz,
                 const std::vector<double> &area, const std::vector<int> &bc_type,
                 const std::vector<double> &xc, const std::vector<double> &yc,
                 const std::vector<double> &zc, const std::vector<double> &wall_dist,
                 const std::vector<int> &active,
                 double dx, double dy, double dz);
void set_freestream(DeviceBundle &dev, double rho, double u, double v, double w, double p);
void set_postshock(DeviceBundle &dev, double rho, double u, double v, double w, double p);
void upload_U(DeviceBundle &dev, const std::vector<double> &U);
void download_U(DeviceBundle &dev, std::vector<double> &U);
void upload_R(DeviceBundle &dev, const std::vector<double> &R);
void launch_net_flux(DeviceBundle &dev, int scheme, int order,
                     bool weno_char, bool weno_z,
                     int ricca_sensor, double ricca_rh_thresh, int weno_hybrid);
void launch_update(DeviceBundle &dev, double dt);
void launch_viscous_flux(DeviceBundle &dev, double mu, double Pr, int turb_model,
                         double Cs, double Cw, double Cv, double Pr_t, double kappa,
                         double A_plus, double Uinf);
double launch_min_dt(DeviceBundle &dev, double cfl, bool viscous = false, double mu = 0.0);

} // namespace euler3d_cuda

#endif
