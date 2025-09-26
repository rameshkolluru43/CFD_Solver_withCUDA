#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Shock_Reflection_2D()
{

	Is_Time_Dependent = false;
	Terminating_Time = 10000.0;
	Viscous_Time_Case = 1;

	// Kind of Boundary Conditions to be implemented in the code
	Is_Viscous_Wall = false;
	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	has_Symmetry_BC = false;

	Directory_Name();
	File_Name();

	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_241_181.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_241_181.vtk";
		Error_File += "_241_181.txt";
		Initial_Solution_File += "_241_181.txt";
		Solution_File += "_241_181.txt";
		Final_Solution_File += "_241_181.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_71_61.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_71_61.vtk";
		Error_File += "_71_61.txt";
		Initial_Solution_File += "_71_61.txt";
		Solution_File += "_71_61.txt";
		Final_Solution_File += "_71_61.vtk";
		break;
	case 5:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_141_121.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_141_121.vtk";
		Error_File += "_Re100000_141_121.txt";
		Initial_Solution_File += "_Re100000_141_121.txt";
		Solution_File += "_Re100000_141_121.txt";
		Final_Solution_File += "_Re100000_141_121.vtk";
		break;
	case 6:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_281_241.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_281_241.vtk";
		Error_File += "_281_241.txt";
		Initial_Solution_File += "_281_241.txt";
		Solution_File += "_281_241.txt";
		Final_Solution_File += "_281_241.vtk";
		break;
	case 7:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_601_521.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_SG_Viscous_601_521.vtk";
		Error_File += "_601_521.txt";
		Initial_Solution_File += "_601_521.txt";
		Solution_File += "_601_521.txt";
		Final_Solution_File += "_601_521.vtk";
		break;
	case 2:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_121_41.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_121_41.vtk";
		Error_File += "_121_41.txt";
		Initial_Solution_File += "_121_41.txt";
		Solution_File += "_121_41.txt";
		Final_Solution_File += "_121_41.vtk";
		break;
	case 9:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_961_321.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_961_321.vtk";
		Error_File += "_961_321.txt";
		Initial_Solution_File += "_961_321.txt";
		Solution_File += "_961_321.txt";
		Final_Solution_File += "_961_321.vtk";
		break;
	case 10:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_481_161.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_481_161.vtk";
		Error_File += "_481_161.txt";
		Initial_Solution_File += "_481_161.txt";
		Solution_File += "_481_161.txt";
		Final_Solution_File += "_481_161.vtk";
		break;
	case 11:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_61_21.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_61_21.vtk";
		Error_File += "_61_21.txt";
		Initial_Solution_File += "_61_21.txt";
		Solution_File += "_61_21.txt";
		Final_Solution_File += "_61_21.vtk";
		break;
	case 8:
		Grid_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_241_81.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Reflection_Files/Shock_Reflection_2D_241_81.vtk";
		Error_File += "_241_81.txt";
		Initial_Solution_File += "_241_81.txt";
		Solution_File += "_241_81.txt";
		Final_Solution_File += "_241_81.vtk";
		break;

	case 3:
		Grid_File = "../Grid_Files/SWBLI/SWBLI_141_121.txt";
		Grid_Vtk_File = "../Grid_Files/SWBLI/SWBLI_141_121.vtk";
		Error_File += "_141_121.txt";
		Initial_Solution_File += "_141_121.txt";
		Solution_File += "_141_121.txt";
		Final_Solution_File += "_141_121.vtk";
		break;
	}

	// string Cell_Info_File = "../Grid_Files/Shock_Reflection_Files/Grid_Cell_Details_241_81.txt";
	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	// Write_Cell_Info(Cell_Info_File);

	V_D V(2, 0.0);

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";
	cin >> Initialize_Type;

	if (Is_Viscous_Wall)
	{
		Re = 100000;
		L_ref = 0.8; // Xsh
		M_ref = 2.15;
		Rho_ref = 1.0;
		cp_ref = 1.0 / ((gamma - 1.0) * M_ref * M_ref);
		R_ref = 1.0 / (gamma * M_ref * M_ref);
		P_ref = 100000;

		P_inf = 985.01;
		u_inf = 1.0;
		T_inf = 111.56;

		mu_ref = Rho_ref * u_inf * L_ref / Re;
		Inv_Re = 1.0 / Re;
		Pr = 0.72;
		Inv_Pr = 1.0 / Pr;
		Inlet_Mach_No = 2.15;
		V_1 = Inlet_Mach_No * sqrt(gamma * P_ref / Rho_ref);
		V_2 = 0.0;
	}
	else
	{ // Test Conditions for Invisicd Case
		P_ref = (1.0 / 1.4);
		Rho_ref = 1.0;
		V_1 = 2.9;
		V_2 = 0.0;
	}

	// cout<<Re<<"\t"<<u_ref<<"\t"<<mu_ref<<endl;
	if (Initialize_Type == 1)
	{
		Initialize(Solution_File);
	}
	else
	{
		Initialize(Test_Case);
		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (unsigned int index = 0; index < No_Physical_Cells; index++)
		{

			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;

			V[0] = V_1;
			V[1] = V_2;

			Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
			//		Print(Global_U);
			for (unsigned int j = 0; j < Global_U.size(); j++)
				U_Cells[index][j] = Global_U[j];
			Calculate_Primitive_Variables(index, U_Cells[index]);
			for (unsigned int j = 0; j < Global_Primitive.size(); j++)
				Primitive_Cells[index][j] = Global_Primitive[j];
		}
	}
	Identify_Wall_Boundary_Faces(Grid_Type);

	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);
}
