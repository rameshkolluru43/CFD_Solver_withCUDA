/**
 * @file Solver.cpp
 * @brief This file contains the implementation of solvers for Euler equations and Navier-Stokes equations.
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

// Core logic for solving Euler equations based on Test Case

void Inviscid_Solver(string &Error_Filename, string &Sol_Filename)
{
	double start_time, end_time;
	int Solution_Data_Type = 1;
	iterations = 0;
	Total_Time = 0.0;
	// cout << "Using Inviscid Solver" << endl;
	// Set fixed floating point format and desired precision
	cout << scientific << setprecision(6);

	// Print header with specified widths for each column
	cout << setw(10) << "Iter"
		 << setw(15) << "dt"
		 << setw(15) << "Rho_Error"
		 << setw(15) << "Rho_u_Error"
		 << setw(15) << "Rho_v_Error"
		 << setw(15) << "Rho_Et_Error"
		 << setw(20) << "Wall_Clock_Time"
		 << setw(15) << "Total_Time"
		 << setw(15) << "OMP_Time" << endl;

	do
	{
		//		start_time = omp_get_wtime();
		int timer = clock();
		// Applies the boundary conditions for the cells

		Apply_Boundary_Conditions();
		// cout << "Applied Boundary Conditions" << endl;

		//  cout<<"In Solver solving with Implicit method\t"<<Is_Implicit_Method<<endl;
		if (Time_Accurate)
		{
			// cout << "Using Runge Kutta Method" << endl;
			Runge_Kutta_Method();
		}
		else
		{
			if (Is_Implicit_Method)
			{
				Implicit_Method();
			}
			else
			{
				Explicit_Method();
			}
		}

		Total_Time += Min_dt;

		iterations++;
		Estimate_Error();
		// cout << "Error Estimated" << endl;
		Update();
		// cout << "Updated Conservative Variables" << endl;
		if ((Total_Time >= Terminating_Time) and (Is_Time_Dependent))
		{
			// cout<<"Maximum and Minimum Time Step in iteration\t"<<Max_dt<<"\t"<<Min_dt<<endl;
			cout << setw(10) << iterations
				 << setw(15) << Min_dt
				 << setw(15) << Error[0]
				 << setw(15) << Error[1]
				 << setw(15) << Error[2]
				 << setw(15) << Error[3]
				 << setw(20) << (timer / CLOCKS_PER_SEC)
				 << setw(15) << Total_Time
				 << endl;
			Write_Error_File(Error_Filename);
			Write_Solution(Sol_Filename, Solution_Data_Type);
			Append_Solution(Solution_File, Final_Solution_File);
		}
		if (iterations % 1000 == 0)
		{
			// cout<<"Maximum and Minimum Time Step in iteration\t"<<Max_dt<<"\t"<<Min_dt<<endl;
			timer = clock();
			Write_Error_File(Error_Filename);
			Write_Solution(Sol_Filename, Solution_Data_Type);
			Read_Write_Grid(Grid_Vtk_File, Final_Solution_File);
			Append_Solution(Solution_File, Final_Solution_File);
			// cout << "Updated Solution File Sucessfully" << endl;
			cout << setw(10) << iterations
				 << setw(15) << Min_dt
				 << setw(15) << Error[0]
				 << setw(15) << Error[1]
				 << setw(15) << Error[2]
				 << setw(15) << Error[3]
				 << setw(20) << (timer / CLOCKS_PER_SEC)
				 << setw(15) << Total_Time
				 << endl;
		}
	} while (iterations < Total_Iterations);
}

// Main story line of the NS equations based on Boundary condition type
void Viscous_Solver(string &Error_Filename, string &Sol_Filename)
{
	int Solution_Data_Type = 1;
	iterations = 0;
	cout << "Min_dt\tIterations\tRho_Error\tRho_u_Error\tRho_v_Error\tRho_Et_Error\tWall_Clock_Time\tTotal_Time" << endl;
	do
	{
		int timer = clock();
		Apply_Boundary_Conditions();
		//		cout<<"Applied Boundary Conditions"<<endl;
		Evaluate_Viscous_Fluxes();
		//	cout<<"Viscous Fluxes evaluated"<<endl;
		// Evaluates the Convective fluxes for all the cells
		if (Is_Second_Order)
			Evaluate_Cell_Net_Flux_2O();
		else
			Evaluate_Cell_Net_Flux_1O();

		// Solves the first order Euler method using Explicit Scheme
		if (Time_Accurate)
			Runge_Kutta_Method();
		else if (Is_Implicit_Method)
			Implicit_Method();
		else
			Explicit_Method();
		Total_Time += Min_dt;
		Estimate_Error();
		//		cout<<"Error Estimated"<<endl;
		Update();

		if ((Total_Time >= Terminating_Time) and (Is_Time_Dependent))
		{
			//	cout<<"Minimum Time Step\tIterations\tRho_Error\tRho_u_Error\tRho_v_Error\tRho_Et_Error\t Wall Clock Time\t Total Time"<<endl;
			cout << Min_dt << "\t" << iterations << "\t Error\t" << Error[0] << "\t" << Error[1] << "\t" << Error[2] << "\t" << Error[3] << "\t" << timer / CLOCKS_PER_SEC << "\t" << Total_Time << endl;
			Write_Error_File(Error_Filename);
			Write_Solution(Sol_Filename, Solution_Data_Type);
			Read_Write_Grid(Grid_Vtk_File, Final_Solution_File);
			Append_Solution(Solution_File, Final_Solution_File);
			timer = clock();
			exit(0);
		}
		else if (iterations % 500 == 0)
		{
			//			cout<<"Maximum and Minimum Time Step in iteration in \t"<<iterations<<"\t"<<Max_dt<<"\t"<<Min_dt<<endl;
			Write_Error_File(Error_Filename);
			Write_Solution(Sol_Filename, Solution_Data_Type);
			Read_Write_Grid(Grid_Vtk_File, Final_Solution_File);
			Append_Solution(Solution_File, Final_Solution_File);
			// 			cout<<"Updated Solution File Sucessfully"<<endl;
			if (Is_Viscous_Wall)
			{
				//				cout<<"Evaluating Skin Friction Coefficient"<<endl;
				Evaluate_Wall_Skin_Friction();
				//				cout<<"Writing Skin Friction Coefficient\t"<<CF_File<<endl;
				Write_CF_File(CF_File);
			}
			timer = clock();
			cout << Min_dt << "\t" << iterations << "\t" << Error[0] << "\t" << Error[1] << "\t" << Error[2] << "\t" << Error[3] << "\t" << timer / CLOCKS_PER_SEC << "\t" << Total_Time << endl;
			//			cout<<"---------------------------------------------------------------------------"<<endl;
		}
		iterations++;
	} while (iterations < Total_Iterations);
	cout << "Iterations Completed \t" << iterations << endl;
	cout << "Evaluating Skin Friction Coefficient" << endl;
	try
	{
		Evaluate_Wall_Skin_Friction();
	}
	catch (const std::exception &e)
	{
		cerr << "Error during wall skin friction evaluation: " << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Writing Skin Friction Coefficient\t" << CF_File << endl;
	Write_CF_File(CF_File);
}
