#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Weno.h"
#include "Timestep.h"
#include "Solver.h"
#include "Viscous_Functions.h"
#include "Boundary_Conditions.h"
#include "Error_Update.h"

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