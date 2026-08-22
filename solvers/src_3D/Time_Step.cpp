/**
 * @file Time_Step.cpp
 * @brief 3D Time stepping calculations for explicit and implicit schemes
 *
 * This module implements 3D time step calculations based on CFL stability conditions
 * for both inviscid (Euler) and viscous (Navier-Stokes) equations. The time step
 * is computed considering all 6 faces of hexahedral cells and 3D eigenvalue systems.
 *
 * Key Features:
 * - 3D CFL stability analysis for hexahedral cells
 * - Inviscid time stepping (Euler equations)
 * - Viscous time stepping (Navier-Stokes equations)
 * - Multiple viscous time step formulations
 * - Adaptive time stepping capability
 * - Stability enforcement for 3D flows
 *
 * Mathematical Framework:
 * - 3D CFL condition: dt = CFL * min(dx,dy,dz) / (|u|+|v|+|w|+c)
 * - Volume-based time step: dt = CFL * V / Σ(|λ_i| * A_i)
 * - Viscous stability: Additional constraint from diffusion terms
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"

/**
 * @brief Calculate viscous time step using method 1 (neighbor-averaged approach)
 *
 * This method computes the time step for viscous flows by considering both
 * convective and diffusive stability constraints. It uses neighbor-averaged
 * properties for viscous terms calculation.
 *
 * @param Cell_No Current cell index
 *
 * Mathematical formulation:
 * - Convective constraint: Λ_c = Σ(|u_n| + c) * A_face
 * - Viscous constraint: Λ_v = Σ(max(4μ/3ρ, γμ/ρPr) * A²/V)
 * - Time step: dt = CFL * V / (Λ_c + C * Λ_v)
 */
void Viscous_Time_Step_1_3D(int &Cell_No)
{
    // 3D face normal components and areas
    double nx[NUM_FACES_3D], ny[NUM_FACES_3D], nz[NUM_FACES_3D], area[NUM_FACES_3D];

    // Current cell properties
    double u0, v0, w0, C0, rho0, mu0;

    // Neighbor properties
    double u_neighbor[NUM_FACES_3D], v_neighbor[NUM_FACES_3D], w_neighbor[NUM_FACES_3D];
    double C_neighbor[NUM_FACES_3D], rho_neighbor[NUM_FACES_3D], mu_neighbor[NUM_FACES_3D];

    // Stability eigenvalues
    double Lambda_C = 0.0, Lambda_V = 0.0;
    double C_viscous = 4.0; // Viscous stability constant

    // Get neighbor cells for all 6 faces
    int neighbors[NUM_FACES_3D];
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        neighbors[face] = Cells[Cell_No].Neighbours[face];
    }

    // Extract face properties (normals and areas)
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int index = face * 3; // 3D normal components
        nx[face] = Cells[Cell_No].Face_Normals[index + 0];
        ny[face] = Cells[Cell_No].Face_Normals[index + 1];
        nz[face] = Cells[Cell_No].Face_Normals[index + 2];
        area[face] = Cells[Cell_No].Face_Areas[face];
    }

    // Current cell primitive variables
    u0 = Primitive_Cells_3D[Cell_No][1];   // u-velocity
    v0 = Primitive_Cells_3D[Cell_No][2];   // v-velocity
    w0 = Primitive_Cells_3D[Cell_No][3];   // w-velocity (3D extension)
    C0 = Primitive_Cells_3D[Cell_No][5];   // Speed of sound (adjusted index for 3D)
    rho0 = Primitive_Cells_3D[Cell_No][0]; // Density
    mu0 = Primitive_Cells_3D[Cell_No][8];  // Dynamic viscosity (adjusted index)

    // Neighbor primitive variables
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor = neighbors[face];
        u_neighbor[face] = Primitive_Cells_3D[neighbor][1];
        v_neighbor[face] = Primitive_Cells_3D[neighbor][2];
        w_neighbor[face] = Primitive_Cells_3D[neighbor][3]; // 3D w-velocity
        C_neighbor[face] = Primitive_Cells_3D[neighbor][5];
        rho_neighbor[face] = Primitive_Cells_3D[neighbor][0];
        mu_neighbor[face] = Primitive_Cells_3D[neighbor][8];
    }

    // Average viscosity from all neighbors (3D stencil)
    double mu_avg = mu0;
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        mu_avg += mu_neighbor[face];
    }
    mu_avg /= (NUM_FACES_3D + 1); // 7-point stencil average

    // Calculate convective eigenvalue (sum over all 6 faces)
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Face-normal velocity (3D)
        double u_n_avg = 0.5 * ((u0 * nx[face] + v0 * ny[face] + w0 * nz[face]) +
                                (u_neighbor[face] * nx[face] + v_neighbor[face] * ny[face] + w_neighbor[face] * nz[face]));

        // Average acoustic speed
        double c_avg = 0.5 * (C0 + C_neighbor[face]);

        // Convective eigenvalue contribution
        Lambda_C += (fabs(u_n_avg) + c_avg) * area[face];
    }

    // Calculate viscous eigenvalue (sum over all 6 faces)
    double inv_volume = Cells[Cell_No].Inv_Volume; // 1/Volume for 3D

    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Average density for this face
        double rho_avg = 0.5 * (rho0 + rho_neighbor[face]);

        // Viscous diffusion coefficients
        double visc_coeff = max(4.0 / (3.0 * rho_avg), gamma / rho_avg);

        // Viscous eigenvalue contribution (area²/volume scaling)
        Lambda_V += visc_coeff * mu_avg * Inv_Pr * area[face] * area[face] * inv_volume;
    }

    // Final time step calculation
    Cells[Cell_No].del_t = CFL * Cells[Cell_No].Volume / (Lambda_C + C_viscous * Lambda_V);

    // Debug output for problematic cells
    if (isnan(Cells[Cell_No].del_t) || Cells[Cell_No].del_t <= 0.0)
    {
        cout << "3D Viscous Time Step Error in Cell " << Cell_No << endl;
        cout << "Volume: " << Cells[Cell_No].Volume << " Lambda_C: " << Lambda_C
             << " Lambda_V: " << Lambda_V << " dt: " << Cells[Cell_No].del_t << endl;
        exit(1);
    }
}

/**
 * @brief Calculate viscous time step using method 2 (cell-centered approach)
 *
 * This method uses only cell-centered properties to compute viscous time step,
 * avoiding neighbor property averaging. More stable for highly stretched grids.
 *
 * @param Cell_No Current cell index
 */
void Viscous_Time_Step_2_3D(int &Cell_No)
{
    // Face properties
    double area[NUM_FACES_3D];
    double avg_area_x, avg_area_y, avg_area_z; // Average areas in each direction

    // Current cell properties
    double u0, v0, w0, C0, rho0, P0, mu0;
    double inv_volume, C_viscous = 2.0;

    // Extract face areas
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        area[face] = Cells[Cell_No].Face_Areas[face];
    }

    // Average face areas in each coordinate direction
    avg_area_x = 0.5 * (area[Face_0] + area[Face_1]); // Left/Right faces (x-direction)
    avg_area_y = 0.5 * (area[Face_2] + area[Face_3]); // Bottom/Top faces (y-direction)
    avg_area_z = 0.5 * (area[Face_4] + area[Face_5]); // Back/Front faces (z-direction)

    // Current cell primitive variables
    u0 = Primitive_Cells_3D[Cell_No][1];
    v0 = Primitive_Cells_3D[Cell_No][2];
    w0 = Primitive_Cells_3D[Cell_No][3]; // 3D w-velocity
    C0 = Primitive_Cells_3D[Cell_No][5];
    rho0 = Primitive_Cells_3D[Cell_No][0];
    P0 = Primitive_Cells_3D[Cell_No][4];
    mu0 = Primitive_Cells_3D[Cell_No][8];

    inv_volume = Cells[Cell_No].Inv_Volume;

    // Convective eigenvalue (3D extension)
    double Lambda_C = ((fabs(u0) + C0) * avg_area_x +
                       (fabs(v0) + C0) * avg_area_y +
                       (fabs(w0) + C0) * avg_area_z); // 3D extension

    // Viscous eigenvalue (area-squared sum)
    double area_squared_sum = 0.0;
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        area_squared_sum += area[face] * area[face];
    }

    // Viscous diffusion coefficient
    double visc_coeff = max(4.0 / (3.0 * rho0), gamma / rho0);
    double Lambda_V = visc_coeff * (mu0 / Pr) * area_squared_sum * inv_volume;

    // Final time step
    Cells[Cell_No].del_t = CFL * Cells[Cell_No].Volume / (Lambda_C + C_viscous * Lambda_V);

    // Stability check
    if (isnan(Cells[Cell_No].del_t) || Cells[Cell_No].del_t <= 0.0)
    {
        cout << "3D Viscous Time Step 2 Error in Cell " << Cell_No << endl;
        cout << "u0=" << u0 << " v0=" << v0 << " w0=" << w0 << " C0=" << C0 << endl;
        cout << "Lambda_C=" << Lambda_C << " Lambda_V=" << Lambda_V << endl;
        exit(1);
    }
}

/**
 * @brief Calculate viscous time step using method 3 (Reynolds number based)
 *
 * This method uses local Reynolds number to determine the appropriate time step
 * for viscous flows. Particularly useful for high Reynolds number flows.
 *
 * @param Cell_No Current cell index
 */
void Viscous_Time_Step_3_3D(int &Cell_No)
{
    // 3D face properties
    double nx[NUM_FACES_3D], ny[NUM_FACES_3D], nz[NUM_FACES_3D], area[NUM_FACES_3D];

    // Current cell properties
    double u0, v0, w0, C0, rho0, mu0;
    double inv_volume, tau = 1.0;

    // 3D averaged properties
    double avg_nx_x, avg_ny_x, avg_nz_x, avg_area_x; // x-direction averages
    double avg_nx_y, avg_ny_y, avg_nz_y, avg_area_y; // y-direction averages
    double avg_nx_z, avg_ny_z, avg_nz_z, avg_area_z; // z-direction averages

    // Extract face properties
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int index = face * 3;
        nx[face] = Cells[Cell_No].Face_Normals[index + 0];
        ny[face] = Cells[Cell_No].Face_Normals[index + 1];
        nz[face] = Cells[Cell_No].Face_Normals[index + 2];
        area[face] = Cells[Cell_No].Face_Areas[face];
    }

    // Average normal and area components (3D)
    avg_nx_x = 0.5 * (nx[Face_1] - nx[Face_0]); // Right - Left
    avg_ny_x = 0.5 * (ny[Face_1] - ny[Face_0]);
    avg_nz_x = 0.5 * (nz[Face_1] - nz[Face_0]);
    avg_area_x = 0.5 * (area[Face_0] + area[Face_1]);

    avg_nx_y = 0.5 * (nx[Face_3] - nx[Face_2]); // Top - Bottom
    avg_ny_y = 0.5 * (ny[Face_3] - ny[Face_2]);
    avg_nz_y = 0.5 * (nz[Face_3] - nz[Face_2]);
    avg_area_y = 0.5 * (area[Face_2] + area[Face_3]);

    avg_nx_z = 0.5 * (nx[Face_5] - nx[Face_4]); // Front - Back
    avg_ny_z = 0.5 * (ny[Face_5] - ny[Face_4]);
    avg_nz_z = 0.5 * (nz[Face_5] - nz[Face_4]);
    avg_area_z = 0.5 * (area[Face_4] + area[Face_5]);

    // Current cell properties
    u0 = Primitive_Cells_3D[Cell_No][1];
    v0 = Primitive_Cells_3D[Cell_No][2];
    w0 = Primitive_Cells_3D[Cell_No][3]; // 3D w-velocity
    C0 = Primitive_Cells_3D[Cell_No][5];
    rho0 = Primitive_Cells_3D[Cell_No][0];
    mu0 = Primitive_Cells_3D[Cell_No][8];

    inv_volume = Cells[Cell_No].Inv_Volume;

    // 3D eigenvalues in each coordinate direction
    double Lambda_x = (fabs(u0 * avg_nx_x + v0 * avg_ny_x + w0 * avg_nz_x) + C0) * avg_area_x;
    double Lambda_y = (fabs(u0 * avg_nx_y + v0 * avg_ny_y + w0 * avg_nz_y) + C0) * avg_area_y;
    double Lambda_z = (fabs(u0 * avg_nx_z + v0 * avg_ny_z + w0 * avg_nz_z) + C0) * avg_area_z; // 3D extension

    // 3D Reynolds numbers in each direction
    double Re_x = (rho0 * fabs(u0) * avg_area_x) / mu0;
    double Re_y = (rho0 * fabs(v0) * avg_area_y) / mu0;
    double Re_z = (rho0 * fabs(w0) * avg_area_z) / mu0; // 3D extension

    // Minimum Reynolds number (most restrictive)
    double Re_min = min({Re_x, Re_y, Re_z});

    // Stability terms
    double term1 = 1.0 + 2.0 / Re_min;
    double term2 = (Lambda_x + Lambda_y + Lambda_z) * inv_volume; // 3D sum

    // Final time step
    Cells[Cell_No].del_t = tau / (term1 * term2);

    // Stability check
    if (isnan(Cells[Cell_No].del_t) || Cells[Cell_No].del_t <= 0.0)
    {
        cout << "3D Viscous Time Step 3 Error in Cell " << Cell_No << endl;
        cout << "Re_min=" << Re_min << " term1=" << term1 << " term2=" << term2 << endl;
        exit(1);
    }
}

/**
 * @brief Calculate inviscid time step for 3D Euler equations
 *
 * This function computes the time step for inviscid flows based on the 3D CFL
 * stability condition. It considers all 6 faces of hexahedral cells and the
 * 3D eigenvalue system of the Euler equations.
 *
 * @param Cell_No Current cell index
 *
 * Mathematical formulation:
 * - 3D eigenvalues: λ = u*nx + v*ny + w*nz ± c
 * - Stability: dt = CFL * V / Σ(max|λ_i| * A_i)
 */
void Inviscid_Time_Step_3D(int &Cell_No)
{
    // 3D face properties
    double nx[NUM_FACES_3D], ny[NUM_FACES_3D], nz[NUM_FACES_3D], area[NUM_FACES_3D];

    // Current cell properties
    double u0, v0, w0, C0;
    double inv_volume;

    // 3D averaged properties for directional eigenvalues
    double avg_nx_x, avg_ny_x, avg_nz_x, avg_area_x; // x-direction
    double avg_nx_y, avg_ny_y, avg_nz_y, avg_area_y; // y-direction
    double avg_nx_z, avg_ny_z, avg_nz_z, avg_area_z; // z-direction

    // Get neighbors for all 6 faces
    int neighbors[NUM_FACES_3D];
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        neighbors[face] = Cells[Cell_No].Neighbours[face];
    }

    // Extract face properties
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int index = face * 3; // 3D normal components
        nx[face] = Cells[Cell_No].Face_Normals[index + 0];
        ny[face] = Cells[Cell_No].Face_Normals[index + 1];
        nz[face] = Cells[Cell_No].Face_Normals[index + 2];
        area[face] = Cells[Cell_No].Face_Areas[face];
    }

    // Average normal components and areas for each coordinate direction
    avg_nx_x = 0.5 * (nx[Face_1] - nx[Face_0]); // Right - Left faces
    avg_ny_x = 0.5 * (ny[Face_1] - ny[Face_0]);
    avg_nz_x = 0.5 * (nz[Face_1] - nz[Face_0]);
    avg_area_x = 0.5 * (area[Face_0] + area[Face_1]);

    avg_nx_y = 0.5 * (nx[Face_3] - nx[Face_2]); // Top - Bottom faces
    avg_ny_y = 0.5 * (ny[Face_3] - ny[Face_2]);
    avg_nz_y = 0.5 * (nz[Face_3] - nz[Face_2]);
    avg_area_y = 0.5 * (area[Face_2] + area[Face_3]);

    avg_nx_z = 0.5 * (nx[Face_5] - nx[Face_4]); // Front - Back faces
    avg_ny_z = 0.5 * (ny[Face_5] - ny[Face_4]);
    avg_nz_z = 0.5 * (nz[Face_5] - nz[Face_4]);
    avg_area_z = 0.5 * (area[Face_4] + area[Face_5]);

    // Current cell primitive variables
    u0 = Primitive_Cells_3D[Cell_No][1]; // u-velocity
    v0 = Primitive_Cells_3D[Cell_No][2]; // v-velocity
    w0 = Primitive_Cells_3D[Cell_No][3]; // w-velocity (3D extension)
    C0 = Primitive_Cells_3D[Cell_No][5]; // Speed of sound

    inv_volume = Cells[Cell_No].Inv_Volume;

    // 3D eigenvalues in each coordinate direction (Blazek 6.14 extended to 3D)
    double Lambda_x = (fabs(u0 * avg_nx_x + v0 * avg_ny_x + w0 * avg_nz_x) + C0) * avg_area_x;
    double Lambda_y = (fabs(u0 * avg_nx_y + v0 * avg_ny_y + w0 * avg_nz_y) + C0) * avg_area_y;
    double Lambda_z = (fabs(u0 * avg_nx_z + v0 * avg_ny_z + w0 * avg_nz_z) + C0) * avg_area_z; // 3D extension

    // Total eigenvalue (sum of all directional contributions)
    double Lambda_total = Lambda_x + Lambda_y + Lambda_z;

    // 3D CFL-based time step
    Cells[Cell_No].del_t = CFL * Cells[Cell_No].Volume / Lambda_total;

    // Stability and error checking
    if (isnan(Cells[Cell_No].del_t) || Cells[Cell_No].del_t <= 0.0)
    {
        cout << "3D Inviscid Time Step Error in Cell " << Cell_No << endl;
        cout << "Velocities: u=" << u0 << " v=" << v0 << " w=" << w0 << " c=" << C0 << endl;
        cout << "Eigenvalues: Λx=" << Lambda_x << " Λy=" << Lambda_y << " Λz=" << Lambda_z << endl;
        cout << "Volume=" << Cells[Cell_No].Volume << " dt=" << Cells[Cell_No].del_t << endl;

        // Print primitive variables for debugging
        cout << "Primitive variables: ";
        for (int i = 0; i < NUM_PRIMITIVE_3D; i++)
        {
            cout << Primitive_Cells_3D[Cell_No][i] << " ";
        }
        cout << endl;
        exit(1);
    }
}

/**
 * @brief Main time step evaluation function
 *
 * This function selects the appropriate time stepping method based on the
 * flow configuration (viscous vs inviscid) and the selected viscous method.
 *
 * @param Cell_No Current cell index
 */
void Evaluate_Time_Step_3D(int &Cell_No)
{
    if (Is_Viscous_Wall)
    {
        // Viscous flow - select method based on configuration
        switch (Viscous_Time_Case)
        {
        case 1:
            Viscous_Time_Step_1_3D(Cell_No); // Neighbor-averaged method
            break;
        case 2:
            Viscous_Time_Step_2_3D(Cell_No); // Cell-centered method
            break;
        case 3:
            Viscous_Time_Step_3_3D(Cell_No); // Reynolds number method
            break;
        default:
            cout << "Invalid Viscous_Time_Case: " << Viscous_Time_Case << endl;
            exit(1);
        }
    }
    else
    {
        // Inviscid flow
        Inviscid_Time_Step_3D(Cell_No);
    }
}

/**
 * @brief Find minimum time step across all cells for global time stepping
 *
 * This function computes the global minimum time step for explicit time
 * integration schemes. Essential for stability in global time stepping.
 *
 * @return Minimum time step across all physical cells
 */
double Find_Global_Min_Time_Step_3D()
{
    double dt_min = 1e10; // Initialize with large value

    // Loop through all physical cells
    for (int cell = 0; cell < No_Physical_Cells; cell++)
    {
        // Compute time step for this cell
        Evaluate_Time_Step_3D(cell);

        // Update minimum
        if (Cells[cell].del_t < dt_min)
        {
            dt_min = Cells[cell].del_t;
        }
    }

    // Safety check
    if (dt_min <= 0.0 || isnan(dt_min))
    {
        cout << "Error: Invalid minimum time step = " << dt_min << endl;
        exit(1);
    }

    return dt_min;
}

/**
 * @brief Adaptive time step calculation with stability monitoring
 *
 * This function implements adaptive time stepping that adjusts the CFL number
 * based on solution convergence and stability considerations.
 *
 * @param iteration Current iteration number
 * @param residual Current residual level
 * @return Adaptive time step
 */
double Adaptive_Time_Step_3D(int iteration, double residual)
{
    static double dt_previous = 0.0;
    double dt_current;

    // Base time step calculation
    dt_current = Find_Global_Min_Time_Step_3D();

    // Simple adaptive strategy
    if (iteration > 100)
    {
        // Check for residual reduction
        static double residual_previous = 1e10;

        if (residual < residual_previous * 0.95)
        {
            // Good convergence - allow CFL increase
            CFL = min(CFL * 1.02, CFL_max);
        }
        else if (residual > residual_previous * 1.05)
        {
            // Poor convergence - reduce CFL
            CFL = max(CFL * 0.95, CFL_min);
        }

        residual_previous = residual;
    }

    // Recompute with updated CFL
    dt_current = Find_Global_Min_Time_Step_3D();

    // Limit time step change rate for stability
    if (dt_previous > 0.0)
    {
        double dt_ratio = dt_current / dt_previous;
        if (dt_ratio > 1.2)
        {
            dt_current = dt_previous * 1.2; // Limit increase
        }
        else if (dt_ratio < 0.8)
        {
            dt_current = dt_previous * 0.8; // Limit decrease
        }
    }

    dt_previous = dt_current;
    return dt_current;
}