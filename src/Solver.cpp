/**
 * @file Solver.cpp
 * @brief Euler (`Inviscid_Solver`) and Navier–Stokes (`Viscous_Solver`) time loops.
 *
 * On `CFD_solver_gpu`, prefers `Resident_GPU_Explicit_*` when init succeeds; otherwise
 * host Explicit/RK with optional CUDA 1O inviscid kernels. Validated half-cylinder
 * WENO+RICCA NS continues use the host WENO path (see docs/HALF_CYLINDER_VALIDATION.md).
 *
 * The solvers include methods for evaluating time steps, applying boundary conditions, and solving the equations
 * using explicit, implicit, or Runge-Kutta methods. The file also includes logic for handling viscous and inviscid
 * flows, as well as error estimation and solution updates.
 *
 * ## Key Functions:
 * - **Evaluate_Time_Step**: Determines the time step for a given cell based on whether the flow is viscous or inviscid.
 * - **Inviscid_Solver**: Core logic for solving Euler equations for inviscid flows.
 * - **Viscous_Solver**: Core logic for solving Navier-Stokes equations for viscous flows.
 * - **Viscous_Time_Step_X**: Functions for evaluating time steps under different viscous flow conditions.
 * - **Inviscid_Time_Step**: Function for evaluating time steps for inviscid flows.
 *
 * ## Euler Equations:
 * The Euler equations describe the motion of an inviscid fluid and are solved using the `Inviscid_Solver` function.
 * The solver applies boundary conditions, evaluates convective fluxes, and updates the solution iteratively.
 *
 * ### Key Steps in the Euler Solver:
 * 1. **Boundary Conditions**: Apply boundary conditions to ensure physical consistency.
 * 2. **Flux Evaluation**: Compute convective fluxes for all cells.
 * 3. **Time Integration**: Use explicit or implicit methods to integrate the equations in time.
 * 4. **Error Estimation**: Estimate the error in the solution to monitor convergence.
 * 5. **Solution Update**: Update the conservative variables based on the computed fluxes.
 *
 * ### Time Step Calculation:
 * The time step is calculated based on the CFL condition, ensuring numerical stability. The maximum eigenvalues
 * in the x and y directions are used to compute the time step for each cell.
 *
 * ## Viscous Flows:
 * For viscous flows, the Navier-Stokes equations are solved using the `Viscous_Solver` function. Additional
 * viscous fluxes are computed, and the solver supports second-order accuracy.
 *
 * ### Key Features:
 * - Support for multiple viscous time step evaluation methods (`Viscous_Time_Step_1`, `Viscous_Time_Step_2`, `Viscous_Time_Step_3`).
 * - Evaluation of wall skin friction coefficients for viscous walls.
 * - Handling of both time-accurate and steady-state simulations.
 *
 * ## Output:
 * The solver outputs the following:
 * - Iteration details, including time step, errors, and wall clock time.
 * - Solution files containing the computed flow variables.
 * - Skin friction coefficients for viscous walls.
 *
 * ## Notes:
 * - The solver is designed to handle both 2D and 3D grids.
 * - The implementation assumes structured grid data with precomputed face normals and areas.
 * - The CFL number and other parameters must be set appropriately for the specific test case.
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
#include "Grid.h"
#include "AMR.hpp"
#include "MPI_Utils.h"
#include <algorithm>
#include <cmath>
#include <chrono>

#ifdef USE_CUDA
#include "../CUDA_KERNELS/MOVERS_RICCA_Flux_Cuda.h"
#endif

namespace
{
void Write_Final_Solution_Output(const string &Error_Filename, const string &Sol_Filename, const int Solution_Data_Type)
{
	if (!CFD_MPI_Is_Root())
		return;

	Write_Error_File(Error_Filename);
	Write_Solution(Sol_Filename, Solution_Data_Type);
	Read_Write_Grid(Grid_Vtk_File, Final_Solution_File);
	Append_Solution(Solution_File, Final_Solution_File);
}
}

// Core logic for solving Euler equations based on Test Case

bool Inviscid_Solver(string &Error_Filename, string &Sol_Filename)
{
	try
	{
		double start_time, end_time;
		int Solution_Data_Type = 1;
		iterations = 0;
		Total_Time = 0.0;
		// cout << "Using Inviscid Solver" << endl;
		// Set fixed floating point format and desired precision
		if (CFD_MPI_Is_Root())
			cout << scientific << setprecision(6);

		// Print header with specified widths for each column
		if (CFD_MPI_Is_Root())
		{
			cout << setw(10) << "Iter"
				 << setw(15) << "dt"
				 << setw(15) << "Rho_Error"
				 << setw(15) << "Rho_u_Error"
				 << setw(15) << "Rho_v_Error"
				 << setw(15) << "Rho_Et_Error"
				 << setw(20) << "Wall_Clock_Time"
				 << setw(15) << "Total_Time"
				 << setw(15) << "OMP_Time" << endl;
		}

#ifdef USE_CUDA
		bool use_resident_gpu = Resident_GPU_Explicit_Init();
#else
		bool use_resident_gpu = false;
#endif
		const auto wall0 = std::chrono::steady_clock::now();

		do
		{
#ifdef USE_CUDA
			if (use_resident_gpu)
			{
				double err4[4] = {0, 0, 0, 0};
				if (!Resident_GPU_Explicit_Step(Min_dt, err4))
				{
					cerr << "Resident_GPU_Explicit_Step failed; falling back to host path." << endl;
					Resident_GPU_Explicit_Shutdown();
					use_resident_gpu = false;
					Apply_Boundary_Conditions();
					Explicit_Method();
					Estimate_Error();
					Update();
				}
			}
			else
#endif
			{
				Apply_Boundary_Conditions();
				if (Time_Accurate)
					Runge_Kutta_Method();
				else if (Is_Implicit_Method)
					Implicit_Method();
				else
					Explicit_Method();
				Estimate_Error();
				Update();
			}

			Total_Time += Min_dt;
			iterations++;

			if (Enable_AMR && iterations >= AMR_Start_Iteration && iterations > 0 && iterations % AMR_Period == 0)
			{
#ifdef USE_CUDA
				if (use_resident_gpu)
					Resident_GPU_Explicit_Download_Host();
#endif
				AMR_Adaptive_Step();
			}

			const auto wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();

			if ((Total_Time >= Terminating_Time) and (Is_Time_Dependent))
			{
				if (CFD_MPI_Is_Root())
				{
#ifdef USE_CUDA
					if (use_resident_gpu)
						Resident_GPU_Explicit_Download_Host();
#endif
					cout << setw(10) << iterations
						 << setw(15) << Min_dt
						 << setw(15) << Error[0]
						 << setw(15) << Error[1]
						 << setw(15) << Error[2]
						 << setw(15) << Error[3]
						 << setw(20) << wall_s
						 << setw(15) << Total_Time
						 << endl;
					Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
				}
			}
			/* Print every 1000; write VTK/solution only every 10000 (I/O was dominating runtime). */
			if (iterations % 1000 == 0)
			{
				if (CFD_MPI_Is_Root())
				{
					cout << setw(10) << iterations
						 << setw(15) << Min_dt
						 << setw(15) << Error[0]
						 << setw(15) << Error[1]
						 << setw(15) << Error[2]
						 << setw(15) << Error[3]
						 << setw(20) << wall_s
						 << setw(15) << Total_Time
						 << endl;
				}
				if (iterations % 10000 == 0)
				{
#ifdef USE_CUDA
					if (use_resident_gpu)
						Resident_GPU_Explicit_Download_Host();
#endif
					if (CFD_MPI_Is_Root())
						Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
				}
			}

			const double res_max = std::max(std::max(Error[0], Error[1]), std::max(Error[2], Error[3]));
			constexpr double kResTol = 1.0e-14;
			if (iterations > 100 && res_max < kResTol)
			{
				if (CFD_MPI_Is_Root())
				{
					cout << "Residuals below " << kResTol << " at iteration " << iterations
						 << " (max L2 relative residual = " << res_max << ")" << endl;
#ifdef USE_CUDA
					if (use_resident_gpu)
						Resident_GPU_Explicit_Download_Host();
#endif
					Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
				}
				break;
			}
		} while (iterations < Total_Iterations);

#ifdef USE_CUDA
		if (use_resident_gpu)
			Resident_GPU_Explicit_Download_Host();
#endif
		if (CFD_MPI_Is_Root())
		{
			Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
			cout << "Inviscid solver completed successfully after " << iterations << " iterations" << endl;
			cout << "Final Min_dt=" << Min_dt
				 << " errors(rho,rhou,rhov,rhoEt)="
				 << Error[0] << " " << Error[1] << " " << Error[2] << " " << Error[3] << endl;
		}
#ifdef USE_CUDA
		if (use_resident_gpu)
			Resident_GPU_Explicit_Shutdown();
#endif
		return true;
	}
	catch (const std::exception &e)
	{
		cerr << "Exception in Inviscid_Solver: " << e.what() << endl;
		return false;
	}
	catch (...)
	{
		cerr << "Unknown exception occurred in Inviscid_Solver" << endl;
		return false;
	}
}

// Main story line of the NS equations based on Boundary condition type
bool Viscous_Solver(string &Error_Filename, string &Sol_Filename)
{
	try
	{
		int Solution_Data_Type = 1;
		iterations = 0;
		if (CFD_MPI_Is_Root())
		{
			cout << setw(10) << "Iter"
				 << setw(15) << "dt"
				 << setw(15) << "Rho_Error"
				 << setw(15) << "Rho_u_Error"
				 << setw(15) << "Rho_v_Error"
				 << setw(15) << "Rho_Et_Error"
				 << setw(20) << "Wall_Clock_Time"
				 << setw(15) << "Total_Time" << endl;
		}

#ifdef USE_CUDA
		bool use_resident_gpu = Resident_GPU_Explicit_Init();
#else
		bool use_resident_gpu = false;
#endif
		const auto wall0 = std::chrono::steady_clock::now();

		do
		{
#ifdef USE_CUDA
			if (use_resident_gpu)
			{
				double err4[4] = {0, 0, 0, 0};
				if (!Resident_GPU_Explicit_Step(Min_dt, err4))
				{
					cerr << "Resident_GPU_Explicit_Step (NS) failed; falling back to host path." << endl;
					Resident_GPU_Explicit_Shutdown();
					use_resident_gpu = false;
					Apply_Boundary_Conditions();
					Evaluate_Viscous_Fluxes();
					Explicit_Method();
					Estimate_Error();
					Update();
				}
			}
			else
#endif
			{
				/* Explicit_Method (and RK) recompute inviscid flux (WENO/1O/2O)
				   and combine with Cells_Viscous_Flux — do not pre-call 1O here. */
				Apply_Boundary_Conditions();
				Evaluate_Viscous_Fluxes();

				if (Time_Accurate)
					Runge_Kutta_Method();
				else if (Is_Implicit_Method)
					Implicit_Method();
				else
					Explicit_Method();
				Estimate_Error();
				Update();
			}

			Total_Time += Min_dt;
			iterations++;

			const auto wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();

			if ((Total_Time >= Terminating_Time) and (Is_Time_Dependent))
			{
#ifdef USE_CUDA
				if (use_resident_gpu)
					Resident_GPU_Explicit_Download_Host();
#endif
				if (CFD_MPI_Is_Root())
				{
					cout << setw(10) << iterations
						 << setw(15) << Min_dt
						 << setw(15) << Error[0]
						 << setw(15) << Error[1]
						 << setw(15) << Error[2]
						 << setw(15) << Error[3]
						 << setw(20) << wall_s
						 << setw(15) << Total_Time << endl;
					Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
				}
				break;
			}

			if (iterations % 500 == 0)
			{
#ifdef USE_CUDA
				if (use_resident_gpu)
					Resident_GPU_Explicit_Download_Host();
#endif
				if (CFD_MPI_Is_Root())
				{
					Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
					if (Is_Viscous_Wall)
					{
						Evaluate_Wall_Skin_Friction();
						Evaluate_Wall_Heat_Flux();
						Write_CF_File(CF_File);
						if (!QW_File.empty())
							Write_QW_File(QW_File);
					}
					cout << setw(10) << iterations
						 << setw(15) << Min_dt
						 << setw(15) << Error[0]
						 << setw(15) << Error[1]
						 << setw(15) << Error[2]
						 << setw(15) << Error[3]
						 << setw(20) << wall_s
						 << setw(15) << Total_Time << endl;
				}
			}
		} while (iterations < Total_Iterations);

#ifdef USE_CUDA
		if (use_resident_gpu)
			Resident_GPU_Explicit_Download_Host();
#endif
		if (CFD_MPI_Is_Root())
		{
			cout << "Iterations Completed \t" << iterations << endl;
			cout << "Evaluating Skin Friction Coefficient" << endl;
			try
			{
				Evaluate_Wall_Skin_Friction();
				Evaluate_Wall_Heat_Flux();
			}
			catch (const std::exception &e)
			{
				cerr << "Error during wall skin friction / heat flux evaluation: " << e.what() << endl;
				return false;
			}
			Write_Final_Solution_Output(Error_Filename, Sol_Filename, Solution_Data_Type);
			cout << "Writing Skin Friction Coefficient\t" << CF_File << endl;
			Write_CF_File(CF_File);
			if (!QW_File.empty())
			{
				cout << "Writing wall heat flux\t" << QW_File << endl;
				Write_QW_File(QW_File);
			}
			cout << "Viscous solver completed successfully after " << iterations << " iterations"
				 << (use_resident_gpu ? " (resident GPU NS)" : " (host/fallback)") << endl;
		}
#ifdef USE_CUDA
		if (use_resident_gpu)
			Resident_GPU_Explicit_Shutdown();
#endif
		return true;
	}
	catch (const std::exception &e)
	{
		cerr << "Exception in Viscous_Solver: " << e.what() << endl;
		return false;
	}
	catch (...)
	{
		cerr << "Unknown exception occurred in Viscous_Solver" << endl;
		return false;
	}
}
