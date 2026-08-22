/**
 * @file Euler3D.hpp
 * @brief Full 3D compressible Euler/NS solver API (host + CUDA).
 *
 * Conservatives: U = [ρ, ρu, ρv, ρw, ρE]  (kNv = 5)
 * Hex faces: 0:-x 1:+x 2:-y 3:+y 4:-z 5:+z
 *
 * Schemes: 0=LLF, 1=RICCA, 2=ROE, 3=VANLEER
 * Order:   0=1O,  1=MUSCL/2O, 2=WENO5
 * BC type: 0=interior, 1=slip, 2=noslip, 3=extrapolate, 4=freestream
 */
#ifndef CFD_EULER3D_HPP
#define CFD_EULER3D_HPP

#include <string>
#include <vector>

namespace euler3d
{

constexpr int kNv = 5;
constexpr int kFaces = 6;
constexpr double kGamma = 1.4;

enum BCType : int
{
    BC_INTERIOR = 0,
    BC_SLIP = 1,
    BC_NOSLIP = 2,
    BC_EXTRAP = 3,
    BC_FREESTREAM = 4,
    BC_POSTSHOCK = 5 /* Rankine–Hugoniot state (shock generator) */
};

struct Config
{
    int nx = 64;
    int ny = 16;
    int nz = 16;
    double Lx = 1.0;
    double Ly = 0.25;
    double Lz = 0.25;
    double cfl = 0.3;
    int max_iter = 500;
    int print_every = 50;
    int scheme = 0; // LLF/RICCA/ROE/VANLEER
    int order = 0;  // 1O / MUSCL / WENO5
    bool weno_char = true; // characteristic projection (most stable)
    bool weno_z = true;    // Borges WENO-Z weights (vs Jiang–Shu)
    /* RICCA flux sensor (paper): 0=legacy pressure α_contact/α_acoustic, 1=RICCA-RH */
    int ricca_sensor = 0;
    double ricca_rh_threshold = 0.15;
    /* WENO reconstruction hybrid (paper): 0=off, 1=pressure-only, 2=combo (default) */
    int weno_hybrid = 2;
    bool use_cuda = false;
    bool viscous = false;
    double mu = 1.0e-4;
    double Pr = 0.72;
    /* 0=laminar, 1=DNS (molecular NS), 2=RANS mixing-length, 3=LES Smagorinsky,
       4=LES WALE, 5=LES Vreman */
    int turb_model = 0;
    double Cs = 0.10;
    double Cw = 0.50;
    double Cv = 0.07;
    double Pr_t = 0.90;
    double kappa_von = 0.41;
    double A_plus = 26.0;
    double les_perturb = 0.0; /* 3D trip amplitude as a fraction of |U∞| */
    double freestream_rho = 1.0;
    double freestream_u = 0.0;
    double freestream_v = 0.0;
    double freestream_w = 0.0;
    double freestream_p = 1.0;
    /* 15° ramp geometry (used when Case == ramp_15) */
    double ramp_angle_deg = 15.0;
    double ramp_x_start = 0.5;
    double ramp_x_end = 1.0;
    /* Flat-plate SWBLI (Case == swbli_plate) */
    double shock_theta_deg = 8.0; /* flow deflection of the incident shock */
    double shock_x_top = 0.6;     /* x where post-shock BC starts on ymax */
    double y_stretch = 1.12;      /* wall-normal geometric growth (>1 clusters at y=0) */
    double post_rho = 1.0;
    double post_u = 0.0;
    double post_v = 0.0;
    double post_w = 0.0;
    double post_p = 1.0;
    // Per-domain-face BC names: slip | noslip | extrapolate | freestream | postshock
    std::string bc_xmin = "slip";
    std::string bc_xmax = "slip";
    std::string bc_ymin = "slip";
    std::string bc_ymax = "slip";
    std::string bc_zmin = "slip";
    std::string bc_zmax = "slip";
    std::string case_name = "sod_x";
    std::string vtk_out = "2D_Euler_Solutions/Box_3D/sod_3d.vtk";
    std::string restart; /* optional ASCII VTK with cell Density/Pressure/Velocity */
};

struct Mesh
{
    int nx = 0, ny = 0, nz = 0;
    int n_cells = 0;
    double dx = 0, dy = 0, dz = 0;
    double Lx = 0, Ly = 0, Lz = 0;
    std::vector<double> xc, yc, zc;
    std::vector<double> wall_dist;
    std::vector<double> volume;
    std::vector<int> neigh;
    std::vector<double> nxyz;
    std::vector<double> area;
    std::vector<int> wall_face;
    std::vector<int> bc_type; // [c*6+f]
    bool has_nodes = false;
    std::vector<double> nodes; // [(nx+1)*(ny+1)*(nz+1)*3] if has_nodes
};

struct State
{
    std::vector<double> U;
    std::vector<double> prim; // rho,u,v,w,p,a
};

Config load_config(const std::string &json_path);
Mesh make_cartesian_mesh(const Config &cfg);
Mesh make_ramp15_mesh(const Config &cfg);
Mesh make_swbli_plate_mesh(const Config &cfg);
void fill_oblique_postshock(Config &cfg);
Mesh make_mesh_for_case(const Config &cfg);
void init_case(State &st, const Mesh &m, const Config &cfg);
void init_sod_x(State &st, const Mesh &m);
void init_freestream(State &st, const Mesh &m, const Config &cfg);
void prim_from_U(State &st, const Mesh &m);
void apply_bc(State &st, const Mesh &m);
double compute_dt(const State &st, const Mesh &m, const Config &cfg);

void fill_prim_cell(const double U[5], double prim[6]);
void U_from_prim_vals(double rho, double u, double v, double w, double p, double U[5]);
void get_LR_states(const State &st, const Mesh &m, const Config &cfg, int c, int f,
                   double UL[5], double UR[5]);
void face_flux(const double UL[5], const double UR[5], double nx, double ny, double nz,
               double area, int scheme, double flux[5]);
void face_flux(const double UL[5], const double UR[5], double nx, double ny, double nz,
               double area, const Config &cfg, double flux[5]);
void net_flux_host(const State &st, const Mesh &m, const Config &cfg, std::vector<double> &R);
void add_viscous_flux_host(const State &st, const Mesh &m, const Config &cfg,
                           std::vector<double> &R);
double mu_sgs_from_grad(int turb_model, double rho, const double g[12], double delta,
                        double y_wall, double mu_lam, const Config &cfg);

void update_explicit(State &st, const Mesh &m, const std::vector<double> &R, double dt,
                     double err[5]);
bool write_vtk(const std::string &path, const Mesh &m, const State &st);
bool read_vtk_restart(const std::string &path, const Mesh &m, State &st);
void apply_les_trip(State &st, const Mesh &m, const Config &cfg);
int run_host(const Config &cfg);
int run_cuda(const Config &cfg);

const char *scheme_name(int scheme);
const char *order_name(int order);
const char *ricca_sensor_name(int sensor);
const char *weno_hybrid_name(int hybrid);
const char *turb_name(int model);

} // namespace euler3d

#endif
