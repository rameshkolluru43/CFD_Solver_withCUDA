#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Shock_Wedge_Reflection()
{
	V_D V(2, 0.0);
	int index = 0;

	Is_Viscous_Wall = false;
	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	has_Symmetry_BC = false;

	Is_Time_Dependent = true;
	Terminating_Time = 0.30;

	Directory_Name();
	File_Name();

	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_121_91.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_121_91.vtk";

		Error_File += "_121_91.txt";
		Initial_Solution_File += "_121_91.txt";
		Solution_File += "_121_91.txt";
		Final_Solution_File += "_121_91.vtk";

		break;
	case 2:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_241_181.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_241_181.vtk";

		Error_File += "_241_181.txt";
		Initial_Solution_File += "_241_181.txt";
		Solution_File += "_241_181.txt";
		Final_Solution_File += "_241_181.vtk";

		break;
	case 3:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_481_361.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_481_361.vtk";

		Error_File += "_481_361.txt";
		Initial_Solution_File += "_481_361.txt";
		Solution_File += "_481_361.txt";
		Final_Solution_File += "_481_361.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_961_721.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Shock_Ramp_961_721.vtk";

		Error_File += "_961_721.txt";
		Initial_Solution_File += "_961_721.txt";
		Solution_File += "_961_721.txt";
		Final_Solution_File += "_961_721.vtk";
		break;
	}
	/* Reads Input grid file and does the preprocessing required for grid  * Calculates the Cell Normals, Cell Areas * Checks for Grid  */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	//	Write_Cell_Info("../Grid_Files/Ramp_Grid_Files/Ramp_Grid_Details_241_81.txt");

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";

	cin >> Initialize_Type;

	if (Initialize_Type == 1)
		Initialize(Initial_Solution_File);
	else
	{
		Initialize(Test_Case);
		cout << "Initalizing Data with Inlet conditions \t" << endl;

		cout << nx_c << "\t" << ny_c << endl;
		for (int j = 0; j < (ny_c - 1); j++)
		{
			for (int i = 0; i < (nx_c - 1); i++)
			{
				index = i + j * (nx_c - 1);

				if (i < 0)
				{
					Pressure_Static_Inlet = 35.125;
					Rho_Static_Inlet = 7.208510638;
					V_1 = 4.43182;
					V_2 = 0.0;
				}
				else
				{
					Pressure_Static_Inlet = 1.0;
					Rho_Static_Inlet = 1.4;
					V_1 = 0.0;
					V_2 = 0.0;
				}

				V[0] = V_1;
				V[1] = V_2;

				Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
				for (unsigned int k = 0; k < Global_U.size(); k++)
					U_Cells[index][k] = Global_U[k];
				Calculate_Primitive_Variables(index, U_Cells[index]);
				for (unsigned int k = 0; k < Global_Primitive.size(); k++)
					Primitive_Cells[index][k] = Global_Primitive[k];
			}
		}
	} // End of Initialization else condition
	Identify_Wall_Boundary_Faces(Grid_Type);

	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);
}
