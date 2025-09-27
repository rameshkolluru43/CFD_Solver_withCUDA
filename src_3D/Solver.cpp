/**
 * @file Solver.cpp
 * @brief 3D implementation of solvers for Euler equations and Navier-Stokes equations.
 *
 * This file contains the implementation of solvers specifically designed for 3D hexahedral
 * grids. The solvers include methods for evaluating time steps, applying boundary conditions,
 * and solving the equations using explicit, implicit, or Runge-Kutta methods for 3D flows.
 *
 * ## Key Functions:
 * - **Evaluate_Time_Step_3D**: Determines the time step for a 3D hexahedral cell
 * - **Inviscid_Solver_3D**: Core logic for solving 3D Euler equations
 * - **Viscous_Solver_3D**: Core logic for solving 3D Navier-Stokes equations
 * - **Time step evaluation methods**: Various approaches for 3D viscous and inviscid flows
 *
 * ## 3D Euler Equations:
 * The 3D Euler equations describe the motion of an inviscid fluid in three dimensions:
 * ∂U/∂t + ∂F/∂x + ∂G/∂y + ∂H/∂z = 0
 * where U = [ρ, ρu, ρv, ρw, ρE]ᵀ
 *
 * ## 3D Navier-Stokes Equations:
 * For viscous flows, additional viscous flux terms are included in all three directions,
 * requiring evaluation of viscous gradients and stress tensors in 3D.
 *
 * ## 3D Time Step Calculation:
 * Time step is based on the 3D CFL condition considering all three coordinate directions:
 * Δt = CFL * min(Δx/(|u|+c), Δy/(|v|+c), Δz/(|w|+c))
 */

#include "definitions.h"
#include "Globals.h"
#include "Solver.h"
#include "Boundary_Conditions.h"
#include "Viscous_Functions.h"
#include "IO_Write.h"
#include "Utilities.h"
#include "Error_Update.h"
#include "Flux.h"

/**
 * @brief Evaluate time step for 3D hexahedral cell
 * @param Cell_No Current cell index
 *
 * Computes the time step based on CFL condition for 3D flow.
 * Considers maximum wave speeds in all three coordinate directions.
 */
void Evaluate_Time_Step_3D(const int &Cell_No)
{
    double u = 0.0, v = 0.0, w = 0.0, c = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    double dt_x = 0.0, dt_y = 0.0, dt_z = 0.0;
    double local_dt = 0.0;

    // Get primitive variables
    u = Primitive_Cells[Cell_No][PRIM_U];
    v = Primitive_Cells[Cell_No][PRIM_V];
    w = Primitive_Cells[Cell_No][PRIM_W];
    c = Primitive_Cells[Cell_No][PRIM_C]; // Speed of sound

    // Estimate characteristic lengths in each direction
    // Use face areas and cell volume to estimate grid spacing
    double cell_volume = 1.0 / Cells[Cell_No].Inv_Volume;

    // X-direction characteristic length
    double face_area_x = 0.5 * (Cells[Cell_No].Face_Areas[Face_0] +
                                Cells[Cell_No].Face_Areas[Face_1]);
    dx = cell_volume / face_area_x;

    // Y-direction characteristic length
    double face_area_y = 0.5 * (Cells[Cell_No].Face_Areas[Face_2] +
                                Cells[Cell_No].Face_Areas[Face_3]);
    dy = cell_volume / face_area_y;

    // Z-direction characteristic length
    double face_area_z = 0.5 * (Cells[Cell_No].Face_Areas[Face_4] +
                                Cells[Cell_No].Face_Areas[Face_5]);
    dz = cell_volume / face_area_z;

    // CFL-based time step in each direction
    dt_x = dx / (fabs(u) + c);
    dt_y = dy / (fabs(v) + c);
    dt_z = dz / (fabs(w) + c);

    // Take minimum for stability
    local_dt = min({dt_x, dt_y, dt_z});

    // Apply CFL number
    local_dt *= CFL;

    // Store in global arrays
    Cell_dt[Cell_No] = local_dt;

    // Update global min/max
    if (Cell_No == 0)
    {
        Min_dt = local_dt;
        Max_dt = local_dt;
    }
    else
    {
        Min_dt = min(Min_dt, local_dt);
        Max_dt = max(Max_dt, local_dt);
    }
}

/**
 * @brief Evaluate viscous time step for 3D cell (method 1)
 * @param Cell_No Current cell index
 */
void Viscous_Time_Step_1_3D(const int &Cell_No)
{
    double u = 0.0, v = 0.0, w = 0.0, c = 0.0;
    double rho = 0.0, mu = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    double dt_conv = 0.0, dt_visc = 0.0, local_dt = 0.0;

    // Get primitive variables
    rho = Primitive_Cells[Cell_No][PRIM_RHO];
    u = Primitive_Cells[Cell_No][PRIM_U];
    v = Primitive_Cells[Cell_No][PRIM_V];
    w = Primitive_Cells[Cell_No][PRIM_W];
    c = Primitive_Cells[Cell_No][PRIM_C];

    // Dynamic viscosity (can be computed from temperature/Sutherland's law)
    double T = Primitive_Cells[Cell_No][PRIM_T];
    mu = Mu_Ref * pow(T / T_Ref, 1.5) * (T_Ref + S_Temp) / (T + S_Temp);

    // Characteristic lengths
    double cell_volume = 1.0 / Cells[Cell_No].Inv_Volume;
    dx = cell_volume / (0.5 * (Cells[Cell_No].Face_Areas[Face_0] + Cells[Cell_No].Face_Areas[Face_1]));
    dy = cell_volume / (0.5 * (Cells[Cell_No].Face_Areas[Face_2] + Cells[Cell_No].Face_Areas[Face_3]));
    dz = cell_volume / (0.5 * (Cells[Cell_No].Face_Areas[Face_4] + Cells[Cell_No].Face_Areas[Face_5]));

    // Convective time step
    double dt_x = dx / (fabs(u) + c);
    double dt_y = dy / (fabs(v) + c);
    double dt_z = dz / (fabs(w) + c);
    dt_conv = min({dt_x, dt_y, dt_z});

    // Viscous time step (diffusion stability)
    double nu = mu / rho; // Kinematic viscosity
    double dt_visc_x = 0.5 * dx * dx / nu;
    double dt_visc_y = 0.5 * dy * dy / nu;
    double dt_visc_z = 0.5 * dz * dz / nu;
    dt_visc = min({dt_visc_x, dt_visc_y, dt_visc_z});

    // Take minimum of convective and viscous constraints
    local_dt = min(dt_conv, dt_visc);
    local_dt *= CFL;

    Cell_dt[Cell_No] = local_dt;

    // Update global min/max
    if (Cell_No == 0)
    {
        Min_dt = local_dt;
        Max_dt = local_dt;
    }
    else
    {
        Min_dt = min(Min_dt, local_dt);
        Max_dt = max(Max_dt, local_dt);
    }
}

/**
 * @brief Enhanced viscous time step evaluation for 3D (method 2)
 * @param Cell_No Current cell index
 */
void Viscous_Time_Step_2_3D(const int &Cell_No)
{
    double u = 0.0, v = 0.0, w = 0.0, c = 0.0;
    double rho = 0.0, mu = 0.0, lambda = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    double local_dt = 0.0;

    // Primitive variables
    rho = Primitive_Cells[Cell_No][PRIM_RHO];
    u = Primitive_Cells[Cell_No][PRIM_U];
    v = Primitive_Cells[Cell_No][PRIM_V];
    w = Primitive_Cells[Cell_No][PRIM_W];
    c = Primitive_Cells[Cell_No][PRIM_C];

    // Transport properties
    double T = Primitive_Cells[Cell_No][PRIM_T];
    mu = Mu_Ref * pow(T / T_Ref, 1.5) * (T_Ref + S_Temp) / (T + S_Temp);
    lambda = -2.0 / 3.0 * mu; // Bulk viscosity (Stokes hypothesis)

    // Grid spacing estimates
    double cell_volume = 1.0 / Cells[Cell_No].Inv_Volume;
    double characteristic_length = pow(cell_volume, 1.0 / 3.0); // Cubic root for 3D

    // More sophisticated grid spacing calculation
    dx = characteristic_length;
    dy = characteristic_length;
    dz = characteristic_length;

    // Enhanced stability analysis for 3D viscous flows
    double velocity_magnitude = sqrt(u * u + v * v + w * w);
    double dt_conv = characteristic_length / (velocity_magnitude + c);

    // Viscous diffusion time step
    double kinematic_viscosity = mu / rho;
    double dt_visc = 0.25 * characteristic_length * characteristic_length / kinematic_viscosity;

    // Heat conduction time step (if applicable)
    double thermal_diffusivity = mu * GAMMA / (rho * Pr); // Prandtl number based
    double dt_thermal = 0.25 * characteristic_length * characteristic_length / thermal_diffusivity;

    // Combined constraint
    local_dt = min({dt_conv, dt_visc, dt_thermal});
    local_dt *= CFL;

    Cell_dt[Cell_No] = local_dt;

    // Update global values
    if (Cell_No == 0)
    {
        Min_dt = local_dt;
        Max_dt = local_dt;
    }
    else
    {
        Min_dt = min(Min_dt, local_dt);
        Max_dt = max(Max_dt, local_dt);
    }
}

/**
 * @brief Advanced viscous time step for 3D with anisotropic grid consideration
 * @param Cell_No Current cell index
 */
void Viscous_Time_Step_3_3D(const int &Cell_No)
{
    // Get cell metrics
    double cell_volume = 1.0 / Cells[Cell_No].Inv_Volume;

    // Anisotropic grid spacing calculation
    V_D dx_vec(3), dy_vec(3), dz_vec(3);

    // Approximate grid vectors from face centers
    for (int dir = 0; dir < 3; dir++)
    {
        // X-direction
        dx_vec[dir] = Cells[Cell_No].Face_Centers[3 * Face_1 + dir] -
                      Cells[Cell_No].Face_Centers[3 * Face_0 + dir];
        // Y-direction
        dy_vec[dir] = Cells[Cell_No].Face_Centers[3 * Face_3 + dir] -
                      Cells[Cell_No].Face_Centers[3 * Face_2 + dir];
        // Z-direction
        dz_vec[dir] = Cells[Cell_No].Face_Centers[3 * Face_5 + dir] -
                      Cells[Cell_No].Face_Centers[3 * Face_4 + dir];
    }

    // Grid spacings
    double dx = MAGNITUDE_3D(dx_vec);
    double dy = MAGNITUDE_3D(dy_vec);
    double dz = MAGNITUDE_3D(dz_vec);

    // Primitive variables
    double rho = Primitive_Cells[Cell_No][PRIM_RHO];
    double u = Primitive_Cells[Cell_No][PRIM_U];
    double v = Primitive_Cells[Cell_No][PRIM_V];
    double w = Primitive_Cells[Cell_No][PRIM_W];
    double c = Primitive_Cells[Cell_No][PRIM_C];
    double T = Primitive_Cells[Cell_No][PRIM_T];

    // Transport properties
    double mu = Mu_Ref * pow(T / T_Ref, 1.5) * (T_Ref + S_Temp) / (T + S_Temp);
    double nu = mu / rho;

    // Directional time step constraints
    V_D dt_conv(3), dt_visc(3);

    dt_conv[0] = dx / (fabs(u) + c);
    dt_conv[1] = dy / (fabs(v) + c);
    dt_conv[2] = dz / (fabs(w) + c);

    dt_visc[0] = 0.5 * dx * dx / nu;
    dt_visc[1] = 0.5 * dy * dy / nu;
    dt_visc[2] = 0.5 * dz * dz / nu;

    // Find most restrictive constraint
    double local_dt = min({dt_conv[0], dt_conv[1], dt_conv[2],
                           dt_visc[0], dt_visc[1], dt_visc[2]});

    local_dt *= CFL;
    Cell_dt[Cell_No] = local_dt;

    // Update global values
    if (Cell_No == 0)
    {
        Min_dt = local_dt;
        Max_dt = local_dt;
    }
    else
    {
        Min_dt = min(Min_dt, local_dt);
        Max_dt = max(Max_dt, local_dt);
    }
}

/**
 * @brief 3D Inviscid solver for Euler equations
 * @param Error_Filename Error output file
 * @param Sol_Filename Solution output file
 * @return Success flag
 */
bool Inviscid_Solver_3D(string &Error_Filename, string &Sol_Filename)
{
    try
    {
        double start_time, end_time;
        int Solution_Data_Type = 1;
        iterations = 0;
        Total_Time = 0.0;

        cout << "Using 3D Inviscid Solver" << endl;
        cout << scientific << setprecision(6);

        // Print header for 3D output
        cout << setw(10) << "Iter"
             << setw(15) << "dt"
             << setw(15) << "Rho_Error"
             << setw(15) << "Rho_u_Error"
             << setw(15) << "Rho_v_Error"
             << setw(15) << "Rho_w_Error"
             << setw(15) << "Rho_Et_Error"
             << setw(20) << "Wall_Clock_Time"
             << setw(15) << "Total_Time" << endl;

        do
        {
            int timer = clock();

            // Apply 3D boundary conditions
            Apply_Boundary_Conditions_3D();

            // Time integration method selection
            if (Time_Accurate)
            {
                Runge_Kutta_Method_3D();
            }
            else
            {
                if (Is_Implicit_Method)
                {
                    Implicit_Method_3D();
                }
                else
                {
                    Explicit_Method_3D();
                }
            }

            Total_Time += Min_dt;
            iterations++;

            // Error estimation for 3D
            Estimate_Error_3D();

            // Update conservative variables
            Update_3D();

            // Time-dependent termination check
            if ((Total_Time >= Terminating_Time) && (Is_Time_Dependent))
            {
                cout << setw(10) << iterations
                     << setw(15) << Min_dt
                     << setw(15) << Error[0]
                     << setw(15) << Error[1]
                     << setw(15) << Error[2]
                     << setw(15) << Error[3]
                     << setw(15) << Error[4]
                     << setw(20) << (timer / CLOCKS_PER_SEC)
                     << setw(15) << Total_Time << endl;

                Write_Error_File_3D(Error_Filename);
                Write_Solution_3D(Sol_Filename, Solution_Data_Type);
                Append_Solution_3D(Solution_File, Final_Solution_File);
                break;
            }

            // Periodic output
            if (iterations % 1000 == 0)
            {
                timer = clock();
                Write_Error_File_3D(Error_Filename);
                Write_Solution_3D(Sol_Filename, Solution_Data_Type);
                Read_Write_Grid_3D(Grid_Vtk_File, Final_Solution_File);
                Append_Solution_3D(Solution_File, Final_Solution_File);

                cout << setw(10) << iterations
                     << setw(15) << Min_dt
                     << setw(15) << Error[0]
                     << setw(15) << Error[1]
                     << setw(15) << Error[2]
                     << setw(15) << Error[3]
                     << setw(15) << Error[4]
                     << setw(20) << (timer / CLOCKS_PER_SEC)
                     << setw(15) << Total_Time << endl;
            }

        } while (iterations < Total_Iterations);

        cout << "3D Inviscid solver completed successfully after " << iterations << " iterations" << endl;
        return true;
    }
    catch (const std::exception &e)
    {
        cerr << "Exception in 3D Inviscid_Solver: " << e.what() << endl;
        return false;
    }
    catch (...)
    {
        cerr << "Unknown exception occurred in 3D Inviscid_Solver" << endl;
        return false;
    }
}

/**
 * @brief 3D Viscous solver for Navier-Stokes equations
 * @param Error_Filename Error output file
 * @param Sol_Filename Solution output file
 * @return Success flag
 */
bool Viscous_Solver_3D(string &Error_Filename, string &Sol_Filename)
{
    try
    {
        int Solution_Data_Type = 1;
        iterations = 0;
        Total_Time = 0.0;

        cout << "Using 3D Viscous Solver" << endl;
        cout << setw(15) << "Min_dt"
             << setw(15) << "Iterations"
             << setw(15) << "Rho_Error"
             << setw(15) << "Rho_u_Error"
             << setw(15) << "Rho_v_Error"
             << setw(15) << "Rho_w_Error"
             << setw(15) << "Rho_Et_Error"
             << setw(20) << "Wall_Clock_Time"
             << setw(15) << "Total_Time" << endl;

        do
        {
            int timer = clock();

            // Apply 3D boundary conditions
            Apply_Boundary_Conditions_3D();

            // Evaluate 3D viscous fluxes
            Evaluate_Viscous_Fluxes_3D();

            // Evaluate 3D convective fluxes
            if (Is_Second_Order)
                Evaluate_Cell_Net_Flux_2O_3D();
            else
                Evaluate_Cell_Net_Flux_1O_3D();

            // Time integration
            if (Time_Accurate)
                Runge_Kutta_Method_3D();
            else if (Is_Implicit_Method)
                Implicit_Method_3D();
            else
                Explicit_Method_3D();

            Total_Time += Min_dt;
            iterations++;

            // Error estimation and update
            Estimate_Error_3D();
            Update_3D();

            // Time-dependent termination
            if ((Total_Time >= Terminating_Time) && (Is_Time_Dependent))
            {
                cout << setw(15) << Min_dt
                     << setw(15) << iterations
                     << setw(15) << Error[0]
                     << setw(15) << Error[1]
                     << setw(15) << Error[2]
                     << setw(15) << Error[3]
                     << setw(15) << Error[4]
                     << setw(20) << (timer / CLOCKS_PER_SEC)
                     << setw(15) << Total_Time << endl;

                Write_Error_File_3D(Error_Filename);
                Write_Solution_3D(Sol_Filename, Solution_Data_Type);
                Read_Write_Grid_3D(Grid_Vtk_File, Final_Solution_File);
                Append_Solution_3D(Solution_File, Final_Solution_File);

                if (Is_Viscous_Wall)
                {
                    Evaluate_Wall_Skin_Friction_3D();
                    Write_CF_File_3D(CF_File);
                }
                break;
            }

            // Periodic output and analysis
            if (iterations % 500 == 0)
            {
                Write_Error_File_3D(Error_Filename);
                Write_Solution_3D(Sol_Filename, Solution_Data_Type);
                Read_Write_Grid_3D(Grid_Vtk_File, Final_Solution_File);
                Append_Solution_3D(Solution_File, Final_Solution_File);

                if (Is_Viscous_Wall)
                {
                    Evaluate_Wall_Skin_Friction_3D();
                    Write_CF_File_3D(CF_File);
                }

                timer = clock();
                cout << setw(15) << Min_dt
                     << setw(15) << iterations
                     << setw(15) << Error[0]
                     << setw(15) << Error[1]
                     << setw(15) << Error[2]
                     << setw(15) << Error[3]
                     << setw(15) << Error[4]
                     << setw(20) << (timer / CLOCKS_PER_SEC)
                     << setw(15) << Total_Time << endl;
            }

        } while (iterations < Total_Iterations);

        cout << "3D Viscous solver completed successfully after " << iterations << " iterations" << endl;

        // Final skin friction evaluation
        if (Is_Viscous_Wall)
        {
            cout << "Evaluating final 3D wall skin friction coefficient" << endl;
            try
            {
                Evaluate_Wall_Skin_Friction_3D();
                Write_CF_File_3D(CF_File);
            }
            catch (const std::exception &e)
            {
                cerr << "Error during 3D wall skin friction evaluation: " << e.what() << endl;
                return false;
            }
        }

        return true;
    }
    catch (const std::exception &e)
    {
        cerr << "Exception in 3D Viscous_Solver: " << e.what() << endl;
        return false;
    }
    catch (...)
    {
        cerr << "Unknown exception occurred in 3D Viscous_Solver" << endl;
        return false;
    }
}

/**
 * @brief Evaluate time step for all cells using specified method
 */
void Evaluate_All_Time_Steps_3D()
{
    // Initialize global min/max
    Min_dt = 1e10;
    Max_dt = -1e10;

    if (Is_Viscous_Flow)
    {
        // Choose viscous time step method
        switch (Viscous_Time_Step_Method)
        {
        case 1:
            for (int i = 0; i < No_Physical_Cells; i++)
                Viscous_Time_Step_1_3D(i);
            break;
        case 2:
            for (int i = 0; i < No_Physical_Cells; i++)
                Viscous_Time_Step_2_3D(i);
            break;
        case 3:
            for (int i = 0; i < No_Physical_Cells; i++)
                Viscous_Time_Step_3_3D(i);
            break;
        default:
            for (int i = 0; i < No_Physical_Cells; i++)
                Viscous_Time_Step_1_3D(i);
            break;
        }
    }
    else
    {
        // Inviscid time step
        for (int i = 0; i < No_Physical_Cells; i++)
            Evaluate_Time_Step_3D(i);
    }
}