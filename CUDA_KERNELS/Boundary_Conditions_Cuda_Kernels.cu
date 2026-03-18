/**
 * @file Boundary_Conditions_Cuda_Kernels.cu
 * @brief CUDA kernels for applying boundary conditions in CFD solver
 *
 * Implements GPU-accelerated boundary condition kernels for the 2D Compressible
 * Navier-Stokes CFD solver. All kernels modify U_cells in-place using a
 * ghost-cell-free API.
 *
 * U_cells layout: flat array [num_cells * 4] of conservative variables
 *   per cell: [rho, rho*u, rho*v, E]
 *
 * face_normals layout: flat array [num_cells * 8] with 4 faces * 2 components
 *   per cell: [nx0, ny0, nx1, ny1, nx2, ny2, nx3, ny3]
 *   Boundary face is assumed to be face 0.
 */

#include "Boundary_Conditions_Cuda_Kernels.h"
#include <math.h>

#define GAMMA_VAL 1.4
#define GAMMA_M1  0.4
#define R_GC      287.0

//=============================================================================
// DEVICE UTILITY FUNCTIONS
//=============================================================================

__device__ void calculate_primitive_device(
    const double* U,
    double* primitive)
{
    double rho = U[0];
    if (rho < 1e-12) rho = 1e-12;

    double u = U[1] / rho;
    double v = U[2] / rho;
    double E = U[3];

    double kinetic_energy = 0.5 * rho * (u * u + v * v);
    double p = GAMMA_M1 * (E - kinetic_energy);
    if (p < 1e-12) p = 1e-12;

    double T = p / (rho * R_GC);

    primitive[0] = rho;
    primitive[1] = u;
    primitive[2] = v;
    primitive[3] = p;
    primitive[4] = T;
}

__device__ void calculate_conservative_device(
    const double* primitive,
    double* U)
{
    double rho = primitive[0];
    double u   = primitive[1];
    double v   = primitive[2];
    double p   = primitive[3];

    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = p / GAMMA_M1 + 0.5 * rho * (u * u + v * v);
}

//=============================================================================
// INLET BOUNDARY CONDITIONS
//=============================================================================

__global__ void subsonic_inlet_bc_kernel(
    double *U_cells,
    const int *inlet_cell_list,
    const double *face_normals,
    const InletCondition_GPU bc_data,
    int num_inlet_cells,
    double gamma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_inlet_cells) return;

    int cell = inlet_cell_list[idx];
    int base = cell * 4;

    double prim[5];
    calculate_primitive_device(&U_cells[base], prim);

    double u_ext = prim[1];
    double v_ext = prim[2];

    double rho_bc = bc_data.Rho;
    double p_bc   = bc_data.P;
    double v_mag_sq = u_ext * u_ext + v_ext * v_ext;

    U_cells[base + 0] = rho_bc;
    U_cells[base + 1] = rho_bc * u_ext;
    U_cells[base + 2] = rho_bc * v_ext;
    U_cells[base + 3] = p_bc / GAMMA_M1 + 0.5 * rho_bc * v_mag_sq;
}

__global__ void supersonic_inlet_bc_kernel(
    double *U_cells,
    const int *inlet_cell_list,
    const InletCondition_GPU bc_data,
    int num_inlet_cells,
    double gamma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_inlet_cells) return;

    int cell = inlet_cell_list[idx];
    int base = cell * 4;

    double rho = bc_data.Rho;
    double u   = bc_data.u;
    double v   = bc_data.v;
    double p   = bc_data.P;
    double v_mag_sq = u * u + v * v;

    U_cells[base + 0] = rho;
    U_cells[base + 1] = rho * u;
    U_cells[base + 2] = rho * v;
    U_cells[base + 3] = p / GAMMA_M1 + 0.5 * rho * v_mag_sq;
}

//=============================================================================
// EXIT BOUNDARY CONDITIONS
//=============================================================================

__global__ void subsonic_exit_bc_kernel(
    double *U_cells,
    const int *exit_cell_list,
    const double *face_normals,
    const ExitCondition_GPU bc_data,
    int num_exit_cells,
    double gamma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_exit_cells) return;

    int cell = exit_cell_list[idx];
    int base = cell * 4;

    double prim[5];
    calculate_primitive_device(&U_cells[base], prim);

    double rho = prim[0];
    double u   = prim[1];
    double v   = prim[2];
    double p_bc = bc_data.P;
    double v_mag_sq = u * u + v * v;

    U_cells[base + 0] = rho;
    U_cells[base + 1] = rho * u;
    U_cells[base + 2] = rho * v;
    U_cells[base + 3] = p_bc / GAMMA_M1 + 0.5 * rho * v_mag_sq;
}

__global__ void supersonic_exit_bc_kernel(
    double *U_cells,
    const int *exit_cell_list,
    int num_exit_cells)
{
    /* Zero-order extrapolation: interior values are already in U_cells. No modification needed. */
    (void)U_cells;
    (void)exit_cell_list;
    (void)num_exit_cells;
}

//=============================================================================
// WALL BOUNDARY CONDITIONS
//=============================================================================

__global__ void viscous_wall_bc_kernel(
    double *U_cells,
    const int *wall_cell_list,
    int num_wall_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_wall_cells) return;

    int cell = wall_cell_list[idx];
    int base = cell * 4;

    U_cells[base + 1] = -U_cells[base + 1];
    U_cells[base + 2] = -U_cells[base + 2];
}

__global__ void inviscid_wall_bc_kernel(
    double *U_cells,
    const int *wall_cell_list,
    const double *face_normals,
    int num_wall_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_wall_cells) return;

    int cell = wall_cell_list[idx];
    int base = cell * 4;

    double rho = U_cells[base + 0];
    if (rho < 1e-12) rho = 1e-12;

    double u = U_cells[base + 1] / rho;
    double v = U_cells[base + 2] / rho;

    double nx = face_normals[cell * 8 + 0];
    double ny = face_normals[cell * 8 + 1];

    double v_dot_n = u * nx + v * ny;

    double u_ref = u - 2.0 * v_dot_n * nx;
    double v_ref = v - 2.0 * v_dot_n * ny;

    U_cells[base + 1] = rho * u_ref;
    U_cells[base + 2] = rho * v_ref;
}

//=============================================================================
// SYMMETRY BOUNDARY CONDITION
//=============================================================================

__global__ void symmetry_bc_kernel(
    double *U_cells,
    const int *symmetry_cell_list,
    const double *face_normals,
    int num_symmetry_cells)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_symmetry_cells) return;

    int cell = symmetry_cell_list[idx];
    int base = cell * 4;

    double rho = U_cells[base + 0];
    if (rho < 1e-12) rho = 1e-12;

    double u = U_cells[base + 1] / rho;
    double v = U_cells[base + 2] / rho;

    double nx = face_normals[cell * 8 + 0];
    double ny = face_normals[cell * 8 + 1];

    double v_dot_n = u * nx + v * ny;

    double u_ref = u - 2.0 * v_dot_n * nx;
    double v_ref = v - 2.0 * v_dot_n * ny;

    U_cells[base + 1] = rho * u_ref;
    U_cells[base + 2] = rho * v_ref;
}

//=============================================================================
// FARFIELD BOUNDARY CONDITION
//=============================================================================

__global__ void farfield_bc_kernel(
    double *U_cells,
    const int *farfield_cell_list,
    const double *face_normals,
    const FarfieldCondition_GPU bc_data,
    int num_farfield_cells,
    double gamma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_farfield_cells) return;

    int cell = farfield_cell_list[idx];
    int base = cell * 4;

    double prim[5];
    calculate_primitive_device(&U_cells[base], prim);

    double rho_i = prim[0];
    double u_i   = prim[1];
    double v_i   = prim[2];
    double p_i   = prim[3];

    double rho_inf = bc_data.Rho_inf;
    double u_inf   = bc_data.u_inf;
    double v_inf   = bc_data.v_inf;
    double p_inf   = bc_data.P_inf;

    double nx = face_normals[cell * 8 + 0];
    double ny = face_normals[cell * 8 + 1];

    double a_i   = sqrt(GAMMA_VAL * p_i / rho_i);
    double a_inf = sqrt(GAMMA_VAL * p_inf / rho_inf);

    double v_n_i   = u_i   * nx + v_i   * ny;
    double v_n_inf = u_inf * nx + v_inf * ny;

    double R_plus  = v_n_i   + 2.0 * a_i   / GAMMA_M1;
    double R_minus = v_n_inf - 2.0 * a_inf / GAMMA_M1;

    double v_n_b = 0.5 * (R_plus + R_minus);
    double a_b   = 0.25 * GAMMA_M1 * (R_plus - R_minus);
    if (a_b < 1e-12) a_b = 1e-12;

    double rho_b, u_b, v_b, p_b;

    if (v_n_b >= 0.0) {
        /* Outflow: entropy from interior */
        double s_i = p_i / pow(rho_i, GAMMA_VAL);
        rho_b = pow(a_b * a_b / (GAMMA_VAL * s_i), 1.0 / GAMMA_M1);
        u_b = u_i + (v_n_b - v_n_i) * nx;
        v_b = v_i + (v_n_b - v_n_i) * ny;
        p_b = rho_b * a_b * a_b / GAMMA_VAL;
    } else {
        /* Inflow: entropy from freestream */
        double s_inf = p_inf / pow(rho_inf, GAMMA_VAL);
        rho_b = pow(a_b * a_b / (GAMMA_VAL * s_inf), 1.0 / GAMMA_M1);
        u_b = u_inf + (v_n_b - v_n_inf) * nx;
        v_b = v_inf + (v_n_b - v_n_inf) * ny;
        p_b = rho_b * a_b * a_b / GAMMA_VAL;
    }

    double v_mag_sq = u_b * u_b + v_b * v_b;

    U_cells[base + 0] = rho_b;
    U_cells[base + 1] = rho_b * u_b;
    U_cells[base + 2] = rho_b * v_b;
    U_cells[base + 3] = p_b / GAMMA_M1 + 0.5 * rho_b * v_mag_sq;
}
