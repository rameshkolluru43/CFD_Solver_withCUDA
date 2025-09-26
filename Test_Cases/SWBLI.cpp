#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void SWBLI()
{

	Is_Time_Dependent = true;
	Terminating_Time = 10000.0;

	// Conditions to determine the type and kind of flow to be solved
	Is_Viscous_Wall = true;
	Is_Inlet_SubSonic = false;
	Is_Exit_SubSonic = false;
	has_Symmetry_BC = true;
	Directory_Name();
	File_Name();

	Viscous_Time_Case = 2;
	switch (Grid_Size)
	{

	case 1:
		Grid_File = "../Grid_Files/SWBLI/SWBLI_101_51.txt";
		Grid_Vtk_File = "../Grid_Files/SWBLI/SWBLI_101_51.vtk";
		Error_File += "_Re1e5_SWBLI_101_51.txt";
		Initial_Solution_File += "_Re1e5_SWBLI_101_51.txt";
		Solution_File += "_Re1e5_SWBLI_101_51.txt";
		CF_File += "_Re1e5_SWBLI_101_51.txt";
		Final_Solution_File += "_Re1e5_SWBLI_101_51.vtk";
		break;
	case 2:
		Grid_File = "../Grid_Files/SWBLI/SWBLI_501_251.txt";
		Grid_Vtk_File = "../Grid_Files/SWBLI/SWBLI_501_251.vtk";
		Error_File += "_Zeta10_Re1e5_SWBLI_501_251.txt";
		Initial_Solution_File += "_Zeta10_Re1e5_SWBLI_501_251.txt";
		Solution_File += "_Zeta10_Re1e5_SWBLI_501_251.txt";
		CF_File += "_Zeta10_Re1e5_SWBLI_501_251.txt";
		Final_Solution_File += "_Zeta10_Re1e5_SWBLI_501_251.vtk";
		break;
	case 3:
		Grid_File = "../Grid_Files/SWBLI/SWBLI_1001_501.txt";
		Grid_Vtk_File = "../Grid_Files/SWBLI/SWBLI_1001_501.vtk";
		Error_File += "_Zeta08_Re1e5_SWBLI_1001_501.txt";
		Initial_Solution_File += "_Zeta08_Re1e5_SWBLI_1001_501.txt";
		Solution_File += "_Zeta08_Re1e5_SWBLI_1001_501.txt";
		CF_File += "_Zeta08_Re1e5_SWBLI_1001_501.txt";
		Final_Solution_File += "_Zeta08_Re1e5_SWBLI_1001_501.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/SWBLI/SWBLI_401_201.txt";
		Grid_Vtk_File = "../Grid_Files/SWBLI/SWBLI_401_201.vtk";
		Error_File += "_Zeta10_Re1e5_SWBLI_401_201.txt";
		Initial_Solution_File += "_Zeta10_Re1e5_SWBLI_401_201.txt";
		Solution_File += "_Zeta10_Re1e5_SWBLI_401_201.txt";
		CF_File += "_Zeta10_Re1e5_SWBLI_401_201.txt";
		Final_Solution_File += "_Zeta10_Re1e5_SWBLI_401_201.vtk";
		break;
	case 5:
		Grid_File = "../Grid_Files/SWBLI/SWBLI_801_401.txt";
		Grid_Vtk_File = "../Grid_Files/SWBLI/SWBLI_801_401.vtk";
		Error_File += "_Zeta10_Re1e5_SWBLI_801_401.txt";
		Initial_Solution_File += "_Zeta10_Re1e5_SWBLI_801_401.txt";
		Solution_File += "_Zeta10_Re1e5_SWBLI_801_401.txt";
		CF_File += "_Zeta10_Re1e5_SWBLI_801_401.txt";
		Final_Solution_File += "_Zeta10_Re1e5_SWBLI_801_401.vtk";
		break;
	}

	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	// Write_Cell_Info(Cell_Info_File);

	V_D V(2, 0.0);

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";
	cin >> Initialize_Type;

	Re = 100000;
	Inv_Re = 1.0 / Re;
	Pr = 0.72;
	Inv_Pr = 1.0 / Pr;
	M_ref = 2.15;
	K1 = 1.0 / ((gamma - 1.0) * M_ref * M_ref * Re * Pr);
	L_ref = 0.8;
	Inlet_Mach_No = 2.15;
	P_ref = 1.0 / (gamma * M_ref * M_ref);
	Rho_ref = 1.0;
	Reference_Values();
	//	q_inf = M_ref*sqrt(gamma*P_ref/Rho_ref);

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
		for (int index = 0; index < No_Physical_Cells; index++)
		{

			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;

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
	// Reynolds Number Re and Prandtl Number Pr

	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);
}
