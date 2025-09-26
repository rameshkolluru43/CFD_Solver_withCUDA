#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

void Flat_Plate_Boundary_Layer()
{

	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	Is_Viscous_Wall = true;
	has_Symmetry_BC = true;
	Viscous_Time_Case = 2;
	Is_Time_Dependent = false;
	Terminating_Time = 10000;

	Grid_File = "../Grid_Files/Flat_Plate_Grid_Files/Flat_Plate_BL_21_21.txt";
	Grid_Vtk_File = "../Grid_Files/Flat_Plate_Grid_Files/Flat_Plate_BL_21_21.vtk";
	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	if (!Form_Cells(Grid_File))
	{
		cerr << "Error: Failed to form cells from grid file: " << Grid_File << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Grid_Type used \t" << Grid_Type << endl;
	//	Write_Cell_Info("../Grid_Files/Channel_Files/Channel_Details_101_101.txt");
	V_D V(2, 0.0);

	Directory_Name();
	File_Name();

	Error_File += "_101_101.txt";
	Initial_Solution_File += "_101_101.txt";
	Solution_File += "_101_101.txt";
	Final_Solution_File += "_101_101.vtk";

	// Reynolds Number Re and Prandtl Number Pr

	Re = 10000;
	Pr = 0.72;
	Inv_Pr = 1.0 / Pr;
	Inv_Re = 1.0 / Re;

	L_ref = 1.0;
	M_ref = 0.15;
	Rho_ref = 1.0;

	Reference_Values();

	P_ref = 1.0 / (gamma * M_ref * M_ref);
	cs = sqrt(gamma * P_ref / Rho_ref);
	u_inf = M_ref * cs;
	Inlet_Mach_No = 1.0;

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";

	cin >> Initialize_Type;

	if (Initialize_Type == 1)
	{
		Initialize(Solution_File);
		//		cout<<"Evaluating Skin Friction Coefficient"<<endl;
		Evaluate_Wall_Skin_Friction();
		//		cout<<"Writing Skin Friction Coefficient\t"<<CF_File<<endl;
		Write_CF_File(CF_File);
		//		cout<<"Writing Done"<<endl;
		exit(0);
	}
	else
	{
		Initialize(Test_Case);

		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < Total_No_Cells; index++)
		{
			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet) / u_inf;
			V_2 = 0.0;
			V[0] = V_1;
			V[1] = V_2;
			Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
			//		Print(Global_U);
			for (unsigned int j = 0; j < Global_U.size(); j++)
			{
				U_Cells[index][j] = Global_U[j];
			}
			Calculate_Primitive_Variables(index, U_Cells[index]);
			for (unsigned int j = 0; j < Global_Primitive.size(); j++)
			{
				Primitive_Cells[index][j] = Global_Primitive[j];
			}
		}
	}

	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);
}
