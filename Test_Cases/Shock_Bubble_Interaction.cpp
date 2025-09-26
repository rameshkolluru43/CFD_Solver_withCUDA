#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Shock_Bubble_Interaction()
{

	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */

	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	Is_Viscous_Wall = false;
	has_Symmetry_BC = false;
	Viscous_Time_Case = 1;

	Grid_File = "../Grid_Files/Channel_Grid_Files/Channel_FLow_Inviscid_101_101.txt";
	Grid_Vtk_File = "../Grid_Files/Channel_Grid_Files/Channel_FLow_Inviscid_101_101.vtk";
	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	//	Write_Cell_Info("../Grid_Files/Channel_Files/Channel_Details_41_21.txt");
	V_D V(2, 0.0), P1(3, 0.0), P2(3, 0.0);

	// Bubble Coordinates
	double xb = 1.0, yb = 0.5, rb = 0.1, xc = 0.0, yc = 0.0, dist = 0.0;

	P2[0] = xb;
	P2[1] = yb;

	Is_Time_Dependent = true;
	Terminating_Time = 0.7;

	Directory_Name();
	File_Name();

	Error_File += "_141_141.txt";
	Initial_Solution_File += "_141_141.txt";
	Solution_File += "_141_141.txt";
	Final_Solution_File += "_141_141.vtk";

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";

	cin >> Initialize_Type;

	//							Print(Cells_Cell_Center);

	if (Initialize_Type == 1)
		Initialize(Solution_File);
	else
	{
		Initialize(Test_Case);

		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < No_Physical_Cells; index++)
		{
			Pressure_Static_Inlet = 1.5698 * (1.0 / 1.4);
			Rho_Static_Inlet = 1.25 * 1.3763;
			P1 = Cells_Cell_Center[index];
			// Print(P1);
			Distance_Between_Points(P1, P2, dist);
			if (dist <= rb)
			{
				Rho_Static_Inlet = 2.0;
				cout << "Cell No\t" << index << "\tdistance\t" << dist << endl;
			}

			//                 cout<<Pressure_Static_Inlet<<"\t"<<Rho_Static_Inlet<<"\t"<<Temperature_Static_Inlet<<"\t"<<C<<endl;

			V_1 = 0.0; // Inlet_Mach_No*sqrt(gamma*Pressure_Static_Inlet/Rho_Static_Inlet);
			V_2 = 0.0;
			V[0] = V_1;
			V[1] = V_2;

			Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
			//					 Print(Global_U);
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

	Identify_Wall_Boundary_Faces(Grid_Type);
	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);
	//			exit(0);
}
