#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Weno.h"
#include "Timestep.h"
#include "Solver.h"
#include "Viscous_Functions.h"
#include "Boundary_Conditions.h"
#include "Error_Update.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// Include CUDA matrix assembly functions if available
#ifdef USE_CUDA_MATRIX_ASSEMBLY
#include "Matrix_Assembly_Cuda_Kernels.h"

// Forward declarations for CUDA matrix assembly functions
extern vector<V_D> Assemble_A_CUDA(vector<V_D> &A, double &dt);
extern void Assemble_A1_CUDA(double &dt);
#endif

/* This function evaluates the net flux sum of viscous and inviscid fluxes for a given cell,
 the net flux also called as Residue for a given cell is evaluated here . */

void Explicit_Method()
{
	double inv_Area = 0.0;
	// cout << "Entered Explicit Method" << endl;
	//  Precompute fluxes based on the method selected (WENO, 2nd Order, or 1st Order)
	if (Is_WENO)
	{
		Evaluate_Cell_Net_Flux_WENO();
	}
	else
	{
		if (Is_Second_Order)
			Evaluate_Cell_Net_Flux_2O();
		else
			Evaluate_Cell_Net_Flux_1O();
	}
	// cout << "Evaluated Cell Net Flux" << endl;
	Min_dt = get_Min_dt();
	// cout << "Minimum Time Step\t" << Min_dt << endl;
	if (Min_dt == 0.0)
	{
		fprintf(stderr, "Error: Min_dt is zero. Cannot proceed with the simulation.\n");
		exit(EXIT_FAILURE);
	}
	if (Min_dt < 0.0)
	{
		fprintf(stderr, "Error: Min_dt is negative. Cannot proceed with the simulation.\n");
		exit(EXIT_FAILURE);
	}

	// Loop over all physical cells and apply preconditioning to the fluxes
	// #pragma omp parallel for private(Cell_Index, inv_Area)  // Enable OpenMP for CPU parallelization
	for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		// Precondition the fluxes (both inviscid and viscous) for each cell
		//        applyPreconditioner(Cell_Index, Is_Viscous_Wall);  // Preconditioning as discussed previously

		inv_Area = Cells[Cell_Index].Inv_Area;

		if (Is_Viscous_Wall)
		{
			// Apply preconditioned viscous and inviscid fluxes when viscous fluxes are enabled
			Cells_DelU[Cell_Index][0] = -Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][0] - Cells_Viscous_Flux[Cell_Index][0]);
			Cells_DelU[Cell_Index][1] = -Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][1] - Cells_Viscous_Flux[Cell_Index][1]);
			Cells_DelU[Cell_Index][2] = -Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][2] - Cells_Viscous_Flux[Cell_Index][2]);
			Cells_DelU[Cell_Index][3] = -Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][3] - Cells_Viscous_Flux[Cell_Index][3]);
		}
		else
		{
			// Apply preconditioned inviscid fluxes only
			Cells_DelU[Cell_Index][0] = -Min_dt * inv_Area * Cells_Net_Flux[Cell_Index][0];
			Cells_DelU[Cell_Index][1] = -Min_dt * inv_Area * Cells_Net_Flux[Cell_Index][1];
			Cells_DelU[Cell_Index][2] = -Min_dt * inv_Area * Cells_Net_Flux[Cell_Index][2];
			Cells_DelU[Cell_Index][3] = -Min_dt * inv_Area * Cells_Net_Flux[Cell_Index][3];
		}
	}
}

void Implicit_Method()
{
	std::cerr << "Implicit_Method: not yet fully implemented. "
	          << "Use explicit or Runge-Kutta time stepping." << std::endl;
	return;

#if 0
	const int max_newton_iterations = 50;
	const double newton_tolerance = 1e-8;
	const double linear_solver_tolerance = 1e-6;
	const int max_linear_iterations = 1000;
	const double under_relaxation = 0.8;

	bool use_cuda = false;
#ifdef USE_CUDA_MATRIX_ASSEMBLY
	if (No_Physical_Cells > 1000)
	{
		use_cuda = true;
		std::cout << "Using CUDA-accelerated implicit solver for " << No_Physical_Cells << " cells" << std::endl;
	}
#endif

	// Newton-Raphson iteration loop
	for (int newton_iter = 0; newton_iter < max_newton_iterations; newton_iter++)
	{
		std::cout << "Newton iteration " << newton_iter + 1 << "/" << max_newton_iterations << std::endl;

		// Step 1: Compute current residual R = -Net_Flux for all cells
		Evaluate_Cell_Net_Flux();

		// Store residual vector (negative of net flux)
		vector<V_D> Residual(No_Physical_Cells, V_D(4, 0.0));
		double max_residual = 0.0;

		for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
		{
			for (int var = 0; var < 4; var++)
			{
				Residual[cell_idx][var] = -Cells_Net_Flux[cell_idx][var];
				max_residual = std::max(max_residual, std::abs(Residual[cell_idx][var]));
			}
		}

		std::cout << "  Maximum residual: " << max_residual << std::endl;

		// Check convergence
		if (max_residual < newton_tolerance)
		{
			std::cout << "Implicit method converged in " << newton_iter << " Newton iterations" << std::endl;
			break;
		}

		// Step 2: Assemble Jacobian matrix A = I/dt + ∂F/∂U
		std::cout << "  Assembling Jacobian matrix..." << std::endl;

		vector<V_D> A; // Jacobian matrix

		if (use_cuda && No_Physical_Cells > 5000)
		{
#ifdef USE_CUDA_MATRIX_ASSEMBLY
			// Use sparse CUDA matrix assembly for large problems
			Assemble_A1_CUDA(dt);

			// Convert sparse format to dense for iterative solver
			// Note: This is a simplified approach. For production, use sparse iterative solvers
			int matrix_size = 4 * No_Physical_Cells;
			A.resize(matrix_size, V_D(matrix_size, 0.0));

			// Populate dense matrix from sparse arrays (if available)
			// This would need to be implemented based on global sparse matrix storage
			std::cout << "  Using CUDA sparse matrix assembly" << std::endl;
#endif
		}
		else if (use_cuda)
		{
#ifdef USE_CUDA_MATRIX_ASSEMBLY
			// Use dense CUDA matrix assembly for smaller problems
			A = Assemble_A_CUDA(A, dt);
			std::cout << "  Using CUDA dense matrix assembly" << std::endl;
#endif
		}
		else
		{
			// Use CPU matrix assembly
			A = Assemble_A(A, dt);
			std::cout << "  Using CPU matrix assembly" << std::endl;
		}

		// Step 3: Solve linear system A * ΔU = -R using iterative method
		std::cout << "  Solving linear system..." << std::endl;

		// Initialize solution vector ΔU
		vector<V_D> DeltaU(No_Physical_Cells, V_D(4, 0.0));

		// Use Jacobi iteration to solve the linear system
		vector<V_D> DeltaU_old(No_Physical_Cells, V_D(4, 0.0));

		bool linear_converged = false;
		for (int linear_iter = 0; linear_iter < max_linear_iterations; linear_iter++)
		{
			double linear_residual = 0.0;

			// Jacobi iteration: x_new = D^(-1) * (b - (L+U)*x_old)
			for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
			{
				for (int var = 0; var < 4; var++)
				{
					int global_row = cell_idx * 4 + var;
					double diagonal = A[global_row][global_row];

					if (std::abs(diagonal) > 1e-14)
					{
						double off_diagonal_sum = 0.0;

						// Sum off-diagonal terms
						for (int col = 0; col < 4 * No_Physical_Cells; col++)
						{
							if (col != global_row)
							{
								off_diagonal_sum += A[global_row][col] * DeltaU_old[col / 4][col % 4];
							}
						}

						// Jacobi update
						DeltaU[cell_idx][var] = (Residual[cell_idx][var] - off_diagonal_sum) / diagonal;

						// Calculate residual for convergence check
						double current_residual = std::abs(DeltaU[cell_idx][var] - DeltaU_old[cell_idx][var]);
						linear_residual = std::max(linear_residual, current_residual);
					}
					else
					{
						DeltaU[cell_idx][var] = 0.0;
					}
				}
			}

			// Check linear solver convergence
			if (linear_residual < linear_solver_tolerance)
			{
				std::cout << "    Linear solver converged in " << linear_iter + 1 << " iterations" << std::endl;
				linear_converged = true;
				break;
			}

			// Update old solution for next iteration
			DeltaU_old = DeltaU;

			// Print progress every 100 iterations
			if ((linear_iter + 1) % 100 == 0)
			{
				std::cout << "    Linear iteration " << linear_iter + 1 << ", residual: " << linear_residual << std::endl;
			}
		}

		if (!linear_converged)
		{
			std::cout << "    Warning: Linear solver did not converge" << std::endl;
		}

		// Step 4: Update solution with under-relaxation: U^(n+1) = U^n + α * ΔU
		std::cout << "  Updating solution..." << std::endl;

		double max_delta = 0.0;
		for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
		{
			for (int var = 0; var < 4; var++)
			{
				double delta = under_relaxation * DeltaU[cell_idx][var];
				U_Cells[cell_idx][var] += delta;
				max_delta = std::max(max_delta, std::abs(delta));

				// Ensure physical bounds
				if (var == 0) // Density
				{
					U_Cells[cell_idx][var] = std::max(U_Cells[cell_idx][var], 1e-10);
				}
				else if (var == 3) // Total energy
				{
					U_Cells[cell_idx][var] = std::max(U_Cells[cell_idx][var], 1e-10);
				}
			}
		}

		std::cout << "  Maximum solution change: " << max_delta << std::endl;

		// Step 5: Update primitive variables
		for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
		{
			Calculate_Primitive_Variables(cell_idx, U_Cells[cell_idx]);
		}

		// Step 6: Apply boundary conditions
		Boundary_Conditions();

		// Check for Newton convergence based on solution change
		if (max_delta < newton_tolerance)
		{
			std::cout << "Newton method converged based on solution change" << std::endl;
			break;
		}

		// Safety check for divergence
		if (max_residual > 1e10 || max_delta > 1e10)
		{
			std::cout << "Error: Newton iteration is diverging. Stopping." << std::endl;
			break;
		}
	}

	// Update time step for next iteration
	dt = get_Min_dt();

	std::cout << "Implicit time step completed. dt = " << dt << std::endl;
#endif
}

double get_Min_dt()
{
	// Use std::min_element to find the minimum time step in the Cells vector
	auto min_it = std::min_element(Cells.begin(), Cells.end(), [](const Cell &a, const Cell &b)
								   { return a.del_t < b.del_t; });
	Min_dt = min_it->del_t;
	//	cout << "Minimum Time Step\t" << Min_dt << endl;
	return Min_dt;
}

double get_Max_dt()
{
	// Use std::min_element to find the minimum time step in the Cells vector
	auto max_it = std::max_element(Cells.begin(), Cells.end(), [](const Cell &a, const Cell &b)
								   { return a.del_t < b.del_t; });
	Max_dt = max_it->del_t;

	return Max_dt;
}

void Runge_Kutta_Method()
{
	int Step_Case = 0, Cell_Index;
	double inv_Area = 0.0;
	// Step 1 - This is the First step of Runge-Kutta Method
	Step_Case = 1;
	if (Is_Viscous_Wall)
		Evaluate_Viscous_Fluxes();
	if (Is_WENO)
	{
		Evaluate_Cell_Net_Flux_WENO();
		//					cout<<"Complted WENO NET FLux"<<endl;
	}
	else
	{
		if (Is_Second_Order)
			Evaluate_Cell_Net_Flux_2O();
		else
			Evaluate_Cell_Net_Flux_1O();
	}

	Min_dt = get_Min_dt();

	for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{

		if (Local_Time_Stepping)
			Min_dt = Cells[Cell_Index].del_t;
		inv_Area = Cells[Cell_Index].Inv_Area;

		U_Cells_RK_1[Cell_Index][0] = U_Cells[Cell_Index][0] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][0] - Cells_Viscous_Flux[Cell_Index][0]);
		U_Cells_RK_1[Cell_Index][1] = U_Cells[Cell_Index][1] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][1] - Cells_Viscous_Flux[Cell_Index][1]);
		U_Cells_RK_1[Cell_Index][2] = U_Cells[Cell_Index][2] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][2] - Cells_Viscous_Flux[Cell_Index][2]);
		U_Cells_RK_1[Cell_Index][3] = U_Cells[Cell_Index][3] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][3] - Cells_Viscous_Flux[Cell_Index][3]);
	}

	for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		Update(Cell_Index, Step_Case); // This function updates the Primitive Variables using Updated Values of RK stage 1
	}

	// Step 2
	Step_Case = 2;

	Apply_Boundary_Conditions();

	if (Is_Viscous_Wall)
		Evaluate_Viscous_Fluxes();

	if (Is_WENO)
	{
		Evaluate_Cell_Net_Flux_WENO();
	}
	else
	{
		if (Is_Second_Order)
			Evaluate_Cell_Net_Flux_2O();
		else
			Evaluate_Cell_Net_Flux_1O();
	}

	Min_dt = get_Min_dt();

	for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		//			if(Local_Time_Stepping)
		//				Min_dt = Cells_DelT[Cell_Index];

		// inv_Area = Cells_Inv_Area[Cell_Index];
		inv_Area = Cells[Cell_Index].Inv_Area;

		U_Cells_RK_2[Cell_Index][0] = 0.25 * (U_Cells_RK_1[Cell_Index][0] + 3.0 * U_Cells[Cell_Index][0] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][0] - Cells_Viscous_Flux[Cell_Index][0]));
		U_Cells_RK_2[Cell_Index][1] = 0.25 * (U_Cells_RK_1[Cell_Index][1] + 3.0 * U_Cells[Cell_Index][1] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][1] - Cells_Viscous_Flux[Cell_Index][1]));
		U_Cells_RK_2[Cell_Index][2] = 0.25 * (U_Cells_RK_1[Cell_Index][2] + 3.0 * U_Cells[Cell_Index][2] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][2] - Cells_Viscous_Flux[Cell_Index][2]));
		U_Cells_RK_2[Cell_Index][3] = 0.25 * (U_Cells_RK_1[Cell_Index][3] + 3.0 * U_Cells[Cell_Index][3] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][3] - Cells_Viscous_Flux[Cell_Index][3]));
	}
	// This function updates the Primitive Variables using Updated Values of RK stage 2
	for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		Update(Cell_Index, Step_Case);
	}

	// Final Step
	Apply_Boundary_Conditions();

	if (Is_Viscous_Wall)
		Evaluate_Viscous_Fluxes();

	if (Is_WENO)
	{
		Evaluate_Cell_Net_Flux_WENO();
	}
	else
	{
		if (Is_Second_Order)
			Evaluate_Cell_Net_Flux_2O();
		else
			Evaluate_Cell_Net_Flux_1O();
	}

	Min_dt = get_Min_dt();

	for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		//			if(Local_Time_Stepping)
		//				Min_dt = Cells_DelT[Cell_Index];

		// inv_Area = Cells_Inv_Area[Cell_Index];
		inv_Area = Cells[Cell_Index].Inv_Area;

		Cells_DelU[Cell_Index][0] = (2.0 / 3.0) * (U_Cells_RK_2[Cell_Index][0] - U_Cells[Cell_Index][0] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][0] - Cells_Viscous_Flux[Cell_Index][0]));
		Cells_DelU[Cell_Index][1] = (2.0 / 3.0) * (U_Cells_RK_2[Cell_Index][1] - U_Cells[Cell_Index][1] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][1] - Cells_Viscous_Flux[Cell_Index][1]));
		Cells_DelU[Cell_Index][2] = (2.0 / 3.0) * (U_Cells_RK_2[Cell_Index][2] - U_Cells[Cell_Index][2] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][2] - Cells_Viscous_Flux[Cell_Index][2]));
		Cells_DelU[Cell_Index][3] = (2.0 / 3.0) * (U_Cells_RK_2[Cell_Index][3] - U_Cells[Cell_Index][3] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][3] - Cells_Viscous_Flux[Cell_Index][3]));
	}
}

void Lax_Fedrichs()
{

	int Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0, Neighbour_4 = 0; // Indicates the numbers to neighbours of the cell
	double Dt = 0.0, inv_Area = 0.0;
	;
	for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{

		Neighbour_1 = Cells[Cell_Index].Neighbours[0]; //(i-1,j,k)
		Neighbour_2 = Cells[Cell_Index].Neighbours[1]; //(i,j-1,k)
		Neighbour_3 = Cells[Cell_Index].Neighbours[2]; //(i+1,j,k)
		Neighbour_4 = Cells[Cell_Index].Neighbours[3]; //(i,j+1,k)
		Dt = get_Min_dt();

		inv_Area = Cells[Cell_Index].Inv_Area;

		U_Cells[Cell_Index][0] = 0.25 * (U_Cells[Neighbour_1][0] + U_Cells[Neighbour_2][0] + U_Cells[Neighbour_3][0] + U_Cells[Neighbour_4][0]) - Dt * inv_Area * Cells_Net_Flux[Cell_Index][0];
		U_Cells[Cell_Index][1] = 0.25 * (U_Cells[Neighbour_1][1] + U_Cells[Neighbour_2][1] + U_Cells[Neighbour_3][1] + U_Cells[Neighbour_4][1]) - Dt * inv_Area * Cells_Net_Flux[Cell_Index][1];
		U_Cells[Cell_Index][2] = 0.25 * (U_Cells[Neighbour_1][2] + U_Cells[Neighbour_2][2] + U_Cells[Neighbour_3][2] + U_Cells[Neighbour_4][2]) - Dt * inv_Area * Cells_Net_Flux[Cell_Index][2];
		U_Cells[Cell_Index][3] = 0.25 * (U_Cells[Neighbour_1][3] + U_Cells[Neighbour_2][3] + U_Cells[Neighbour_3][3] + U_Cells[Neighbour_4][3]) - Dt * inv_Area * Cells_Net_Flux[Cell_Index][3];
	}
}