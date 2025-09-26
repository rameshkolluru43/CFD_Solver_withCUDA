#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Shock_Tube_2D()
{

	// Kind of Boundary Conditions to be implemented in the code
	Is_Viscous_Wall = true;
	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	has_Symmetry_BC = false;
	Is_Time_Dependent = true;
	Terminating_Time = 0.2;

	Viscous_Time_Case = 2;

	Directory_Name();
	File_Name();

	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_501_151.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_501_151.vtk";
		Error_File += "_501_151.txt";
		Initial_Solution_File += "_501_151.txt";
		Solution_File += "_501_151.txt";
		Final_Solution_File += "_501_151.vtk";
		break;
	case 2:
		Grid_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_1001_301.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_1001_301.vtk";
		Error_File += "_1001_301.txt";
		Initial_Solution_File += "_1001_301.txt";
		Solution_File += "_1001_301.txt";
		Final_Solution_File += "_1001_301.vtk";
		break;
	case 3:
		Grid_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_201_201.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_201_201.vtk";
		Error_File += "_201_201.txt";
		Initial_Solution_File += "_201_201.txt";
		Solution_File += "_201_201.txt";
		Final_Solution_File += "_201_201.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_101_101.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_101_101.vtk";
		Error_File += "_101_101.txt";
		Initial_Solution_File += "_101_101.txt";
		Solution_File += "_101_101.txt";
		Final_Solution_File += "_101_101.vtk";
		break;
	case 5:
		Grid_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_501_501.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Tube_Grid_Files/Shock_Tube_2D_SG_Alpha03_501_501.vtk";
		Error_File += "_501_501.txt";
		Initial_Solution_File += "_501_501.txt";
		Solution_File += "_501_501.txt";
		Final_Solution_File += "_501_501.vtk";
		break;
	case 6:
		Grid_File = "../Grid_Files/Shock_Tube_Grid_Files/SWBLI_L1_2D_513_513.txt";
		Grid_Vtk_File = "../Grid_Files/Shock_Tube_Grid_Files/SWBLI_L1_2D_513_513.vtk";
		Error_File += "_L1_2D_513_513.txt";
		Initial_Solution_File += "_L1_2D_513_513.txt";
		Solution_File += "_L1_2D_513_513.txt";
		Final_Solution_File += "_L1_2D_513_513.vtk";
		break;
	}

	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	//	Write_Cell_Info("../Grid_Files/Forward_Step_Grid_Files/Forward_Step_241_81.txt");
	V_D V(2, 0.0);
	int index = 0;

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";

	cin >> Initialize_Type;
	// Reynolds Number Re and Prandtl Number Pr
	Re = 25000;
	Pr = 0.72;
	L_ref = 1.0; // Xsh
	M_ref = 1.0;
	Inv_Re = 1.0 / Re;
	Inv_Pr = 1.0 / Pr;
	Inlet_Mach_No = 0.0;

	Reference_Values();

	double Diaph = 0.3;

	if (Initialize_Type == 1)
		Initialize(Final_Solution_File);
	else
	{
		Initialize(Test_Case);
		for (int j = 0; j < (ny_c - 1); j++)
		{
			for (int i = 0; i < (nx_c - 1); i++)
			{
				if (i <= Diaph * (nx_c - 1)) // Pre Shock Conditions Left side of Shock
				{
					index = i + j * (nx_c - 1);
					Pressure_Static_Inlet = 1.0;
					Rho_Static_Inlet = 1.0;
					V_1 = 0.75;
					V_2 = 0.0;
				}
				else // Post Shock Conditions Right Side of SHock
				{
					index = i + j * (nx_c - 1);
					Pressure_Static_Inlet = 0.1;
					Rho_Static_Inlet = 0.125;
					V_1 = 0.0;
					V_2 = 0.0;
				}

				V[0] = V_1;
				V[1] = V_2;

				Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
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
	}

	Identify_Wall_Boundary_Faces(Grid_Type);

	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);

	Read_Write_Grid(Grid_Vtk_File, Final_Solution_File);
	Append_Solution(Initial_Solution_File, Final_Solution_File);
	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
}
