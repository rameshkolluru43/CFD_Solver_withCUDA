#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Forward_Facing_Step()
{

	Directory_Name();
	File_Name();
	double C = 0.0;

	Is_Viscous_Wall = false;
	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	has_Symmetry_BC = false;

	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_241_81.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_241_81.vtk";
		Error_File += "_M3_241_81.txt";
		Initial_Solution_File += "_M3_241_81.txt";
		Solution_File += "_M3_241_81.txt";
		Final_Solution_File += "_M3_241_81.vtk";
		break;
	case 2:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_121_41.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_121_41.vtk";
		Error_File += "_M3_121_41.txt";
		Initial_Solution_File += "_M3_121_41.txt";
		Solution_File += "_M3_121_41.txt";
		Final_Solution_File += "_M3_121_41.vtk";
		break;
	case 3:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_481_161.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_481_161.vtk";
		Error_File += "_M3_481_161.txt";
		Initial_Solution_File += "_M3_481_161.txt";
		Solution_File += "_M3_481_161.txt";
		Final_Solution_File += "_M3_481_161.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_961_321.txt";
		Grid_Vtk_File = "../Grid_Files/Forward_Step_Grid_Files/Forward_Step_Perfect_Gas_961_321.vtk";
		Error_File += "_M3_961_321.txt";
		Initial_Solution_File += "_M3_961_321.txt";
		Solution_File += "_M3_961_321.txt";
		Final_Solution_File += "_M3_961_321.vtk";
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
	Is_Time_Dependent = true;
	Terminating_Time = 4.0;
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

			Pressure_Static_Inlet = 1.0;
			Inlet_Mach_No = 3.0;
			Rho_Static_Inlet = 1.4;
			C = sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			//                 cout<<cp<<"\t"<<cv<<"\t"<<R<<endl;
			//                 cout<<Pressure_Static_Inlet<<"\t"<<Rho_Static_Inlet<<"\t"<<Temperature_Static_Inlet<<"\t"<<C<<endl;
			V_1 = Inlet_Mach_No * C;
			V_2 = 0.0;

			V[0] = V_1;
			V[1] = V_2;

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
