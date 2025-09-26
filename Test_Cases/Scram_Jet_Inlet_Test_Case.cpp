#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Scram_Jet_Inlet()
{
	Directory_Name();
	File_Name();
	double C = 0.0;

	Grid_File = "../Grid_Files/Scramjet/Scram_Jet_Inlet_Flow_300_100.txt";
	Error_File += "_M9p27_300_100.txt";
	Initial_Solution_File += "_M9p27_300_100.txt";
	Solution_File += "_M9p27_300_100.txt";

	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	Write_Cell_Info("../Grid_Files/Scramjet/Scram_Jet_Inlet_Flow_Details_300_100.txt");
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

			Pressure_Static_Inlet = 2188;
			Inlet_Mach_No = 9.27;
			Temperature_Static_Inlet = 222.5;
			Rho_Static_Inlet = Pressure_Static_Inlet / (R_GC * Temperature_Static_Inlet);
			C = sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			//                 cout<<cp<<"\t"<<cv<<"\t"<<R<<endl;

			V_1 = Inlet_Mach_No * C;
			V_2 = 0.0;

			V_Magnitude_Inlet = 0.5 * (V_1 * V_1 + V_2 * V_2);

			cout << Pressure_Static_Inlet << "\t" << Rho_Static_Inlet << "\t" << Temperature_Static_Inlet << "\t" << C << "\t" << V_Magnitude_Inlet << endl;
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
