#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"
#include "Primitive_Computational.h"
#include "IO_Write.h"
#include "Grid.h"
#include "Directory_Files.h"
#include "Initialize.h"
#include "Directory_Files.h"
#include "MPI_Utils.h"
#include "Viscous_Functions.h"
#include "Boundary_Conditions.h"

namespace
{
	/* Fill every cell (physical + ghost) from a uniform (P,ρ,u,v) state. */
	void Fill_All_Cells_Uniform(double P, double Rho, double u, double v)
	{
		V_D V(2, 0.0);
		V[0] = u;
		V[1] = v;
		Calculate_Computational_Variables(P, V, Rho, 2);
		for (int index = 0; index < Total_No_Cells; ++index)
		{
			for (unsigned int j = 0; j < Global_U.size(); ++j)
				U_Cells[index][j] = Global_U[j];
			Calculate_Primitive_Variables(index, U_Cells[index]);
			for (unsigned int j = 0; j < Global_Primitive.size(); ++j)
				Primitive_Cells[index][j] = Global_Primitive[j];
		}
	}

	/* Keep BC helper globals consistent with inlet (used by some legacy BC paths). */
	void Sync_BC_Globals_From_Inlet()
	{
		Pressure_Static_Inlet = inletCond.P;
		Rho_Static_Inlet = inletCond.Rho;
		Inlet_Mach_No = inletCond.M;
		Temperature_Static_Inlet = inletCond.T;
	}

	/*
	 * After allocation / restart load:
	 *  - seed any non-positive ghost states from inlet freestream (never leave zeros)
	 *  - apply BCs so ghosts match inlet/exit/wall/symmetry before the first residual
	 */
	void Finalize_Initial_Field()
	{
		Sync_BC_Globals_From_Inlet();

		V_D V(2, 0.0);
		V[0] = inletCond.u;
		V[1] = inletCond.v;
		Calculate_Computational_Variables(inletCond.P, V, inletCond.Rho, 2);
		for (int index = No_Physical_Cells; index < Total_No_Cells; ++index)
		{
			const bool bad_rho = !(U_Cells[index][0] > 1.0e-14) || !std::isfinite(U_Cells[index][0]);
			if (!bad_rho)
				continue;
			for (unsigned int j = 0; j < Global_U.size(); ++j)
				U_Cells[index][j] = Global_U[j];
			Calculate_Primitive_Variables(index, U_Cells[index]);
			for (unsigned int j = 0; j < Global_Primitive.size(); ++j)
				Primitive_Cells[index][j] = Global_Primitive[j];
		}

		Apply_Boundary_Conditions();
	}
} // namespace

void Initialize_TestCase()
{
	if (CFD_MPI_Is_Root())
		std::cout << "Initializing the Test Case: " << Test_Case << " - " << Test_Case_Name << std::endl;

	/* Set Re/Pr/R_ref/K1 before primitive recovery so Sutherland μ is finite. */
	if (Is_Viscous || Is_Viscous_Wall)
	{
		if (!(Re > 0.0))
			Re = 1.0e5;
		if (!(Pr > 0.0))
			Pr = 0.72;
		Inv_Re = 1.0 / Re;
		Inv_Pr = 1.0 / Pr;
		L_ref = 1.0;
		M_ref = (inletCond.M > 0.0) ? inletCond.M : 6.0;
		Rho_ref = (inletCond.Rho > 0.0) ? inletCond.Rho : 1.4;
		Reference_Values();
		K1 = 1.0 / ((gamma - 1.0) * M_ref * M_ref * Re * Pr);
		Viscous_Time_Case = 2;
		if (CFD_MPI_Is_Root())
			cout << "Viscous refs: Re=" << Re << " Pr=" << Pr
				 << " M_ref=" << M_ref << " K1=" << K1
				 << " Tw=" << wallCond.T << endl;
	}

	if (Initialize_Type == 1)
	{
		/* Restart: load physical cells from Solution_File; do not overwrite with initCond. */
		Initialize(Solution_File);
		if (CFD_MPI_Is_Root())
			cout << "Restart mode (Initialize_Type=1): keeping Solution_File field" << endl;
	}
	else
	{
		/* Cold start: uniform field from Flow_Conditions.Initial_Conditions (all cells). */
		Initialize(Test_Case);
		Fill_All_Cells_Uniform(initCond.P, initCond.Rho, initCond.u, initCond.v);
		if (CFD_MPI_Is_Root())
			cout << "Cold start (Initialize_Type=0): uniform IC "
				 << "P=" << initCond.P << " Rho=" << initCond.Rho
				 << " u=" << initCond.u << " v=" << initCond.v
				 << " M=" << initCond.M << endl;
	}

	Finalize_Initial_Field();

	if (CFD_MPI_Is_Root())
	{
		cout << "Inlet Conditions: " << endl;
		cout << "Pressure_Static_Inlet: " << inletCond.P << endl;
		cout << "Rho_Static_Inlet: " << inletCond.Rho << endl;
		cout << "Inlet_Mach_No: " << inletCond.M << endl;
		cout << "V_1: " << inletCond.u << endl;
		cout << "V_2: " << inletCond.v << endl;
		cout << "Temperature_Static_Inlet: " << inletCond.T << endl;
	}

	Write_Solution(Initial_Solution_File, 1);
	CFD_MPI_Barrier();
	if (CFD_MPI_Is_Root())
		cout << "Initialized Solution, Identified Boundaries... Ready to solve." << std::endl;
}
