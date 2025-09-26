#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Shock_Diffraction()
{

	Is_Viscous_Wall = false;
	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	has_Symmetry_BC = false;
	Viscous_Time_Case = 1;
	Is_Time_Dependent = true;
	Terminating_Time = 0.1561;

	Directory_Name();
	File_Name();

	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_101_101.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_101_101.vtk";
		Error_File += "_101_101.txt";
		Initial_Solution_File += "_101_101.txt";
		Solution_File += "_101_101.txt";
		Final_Solution_File += "_101_101.vtk";
		break;
	case 2:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_401_401.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_401_401.vtk";
		Error_File += "_400_400.txt";
		Initial_Solution_File += "_400_400.txt";
		Solution_File += "_400_400.txt";
		Final_Solution_File += "_400_400.vtk";
		break;
	case 3:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_201_201.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_201_201.vtk";
		Error_File += "_200_200.txt";
		Initial_Solution_File += "_200_200.txt";
		Solution_File += "_200_200.txt";
		Final_Solution_File += "_200_200.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_1281_1281.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Shock_Diffraction_1281_1281.vtk";
		Error_File += "_1281_1281.txt";
		Initial_Solution_File += "_1281_1281.txt";
		Solution_File += "_1281_1281.txt";
		Final_Solution_File += "_1281_1281.vtk";
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

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";

	cin >> Initialize_Type;

	if (Initialize_Type == 1)
		Initialize(Initial_Solution_File);
	else
	{
		Initialize(Test_Case);

		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < Total_No_Cells; index++)
		{
			//				cout<<index<<endl;
			if (index <= ((nx_1 - 1) * (ny_1 - 1)))
			{
				Pressure_Static_Inlet = 30.05945;
				Rho_Static_Inlet = 7.041132896;
				V_1 = 4.0779;
				V_2 = 0.0;
				V[0] = V_1;
				V[1] = V_2;
			}
			else
			{
				Pressure_Static_Inlet = 1.0;
				Rho_Static_Inlet = 1.4;
				V_1 = 0.0;
				V_2 = 0.0;
				V[0] = V_1;
				V[1] = V_2;
			}

			Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);

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
