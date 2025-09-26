#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

void Flow_Over_Bump()
{

	Directory_Name();
	File_Name();

	Viscous_Time_Case = 2;
	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_151_151.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_151_151.vtk";
		Error_File += "_LZ1010_SG2_151_151.txt";
		Initial_Solution_File += "_LZ1010_SG2_151_151.txt";
		Solution_File += "_LZ1010_SG2_151_151.txt";
		Final_Solution_File += "_LZ1010_SG2_151_151.vtk";
		CF_File += "_Re8000_LZ1010_SG2_151_151.txt";

		break;
	case 2:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_301_101.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_301_101.vtk";
		Error_File += "_Viscous_LZ10_SG2_301_101.txt";
		Initial_Solution_File += "_Viscous_LZ10_SG2_301_101.txt";
		Solution_File += "_Viscous_LZ10_SG2_301_101.txt";
		Final_Solution_File += "_Viscous_LZ10_SG2_301_101.vtk";
		CF_File += "_Re8000_LZ10_SG2_301_101.txt";
		break;

		/*	  case 2:
						Grid_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_161_321.txt";
						Grid_Vtk_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_161_321.vtk";
						Error_File +="_Viscous_LZ10_161_321.txt";
						Initial_Solution_File +="_Viscous_LZ10_161_321.txt";
						Solution_File +="_Viscous_LZ10_161_321.txt";
						Final_Solution_File +="_Viscous_LZ10_161_321.vtk";
						CF_File +="_Re8000_LZ10_161_321.txt"; 	*/
		break;
	case 3:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_901_301.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_901_301.vtk";
		Error_File += "_Viscous_LZ10_SG2_901_301.txt";
		Initial_Solution_File += "_Viscous_LZ10_SG2_901_301.txt";
		Solution_File += "_Viscous_LZ10_SG2_901_301.txt";
		Final_Solution_File += "_Viscous_LZ10_SG2_901_301.vtk";
		CF_File += "_Re8000_LZ10_SG2_901_301.txt";
		break;
		/*	  case 3:
						Grid_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_901_301.txt";
						Grid_Vtk_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_901_301.vtk";
						Error_File +="_Viscous_LZ10_901_301.txt";
						Initial_Solution_File +="_Viscous_LZ10_901_301.txt";
						Solution_File +="_Viscous_LZ10_901_301.txt";
						Final_Solution_File +="_Viscous_LZ10_901_301.vtk";
						CF_File +="_Re8000_LZ10_901_301.txt";
			  break;*/
	case 4:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_601_201.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_SG2_601_201.vtk";
		Error_File += "_Viscous_LZ10_SG2_601_201.txt";
		Initial_Solution_File += "_Viscous_LZ10_SG2_601_201.txt";
		Solution_File += "_Viscous_LZ10_SG2_601_201.txt";
		Final_Solution_File += "_Viscous_LZ10_SG2_601_201.vtk";
		CF_File += "_Re8000_LZ10_SG2_601_201.txt";
		break;
		/*	  case 4:
						Grid_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_641_161.txt";
						Grid_Vtk_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_641_161.vtk";
						Error_File +="_Viscous_LZ10_641_161.txt";
						Initial_Solution_File +="_Viscous_LZ1010_641_161.txt";
						Solution_File +="_Viscous_LZ1010_641_161.txt";
						Final_Solution_File +="_Viscous_LZ1010_641_161.vtk";
						CF_File +="_Re8000_LZ1010_641_161.txt";
			  break;*/
		/*	  case 4:
						Grid_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_161_641.txt";
						Grid_Vtk_File="../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_161_641.vtk";
						Error_File +="_Viscous_161_641.txt";
						Initial_Solution_File +="_Viscous_161_641.txt";
						Solution_File +="_Viscous_161_641.txt";
						Final_Solution_File +="_Viscous_161_641.vtk";
						CF_File +="_Re8000_161_641.txt"; 	*/
		break;
	case 5:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_601_201.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_601_201.vtk";
		Error_File += "_601_201.txt";
		Initial_Solution_File += "_601_201.txt";
		Solution_File += "_Inviscid_601_201.txt";
		Final_Solution_File += "_Inviscid_601_201.vtk";
		break;
	case 6:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_301_101.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_301_101.vtk";
		Error_File += "_301_101.txt";
		Initial_Solution_File += "_301_101.txt";
		Solution_File += "_Inviscid_301_101.txt";
		Final_Solution_File += "_Inviscid_301_101.vtk";
		CF_File += "_Re8000_301_101.txt";

		break;
	case 7:
		Grid_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_901_301.txt";
		Grid_Vtk_File = "../Grid_Files/Flow_Over_Cylinder_Files/Flow_Over_thickbump_901_301.vtk";
		Error_File += "_901_301.txt";
		Initial_Solution_File += "_901_301.txt";
		Solution_File += "_Inviscid_901_301.txt";
		Final_Solution_File += "_Inviscid_901_301.vtk";
		break;
	}

	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;

	//	Write_Cell_Info("../Grid_Files/Channel_Files/Channel_Details_101_101.txt");
	V_D V(2, 0.0);
	if (Is_Viscous_Wall)
	{
		Re = 8000;
		Inv_Re = 1.0 / Re;
		Pr = 0.72;
		Inv_Pr = 1.0 / Pr;
		M_ref = 1.4;
		K1 = 1.0 / ((gamma - 1.0) * M_ref * M_ref * Re * Pr);
		L_ref = 3.0;
		Inlet_Mach_No = 1.4;
		P_ref = 1.0 / (gamma * M_ref * M_ref);
		Rho_ref = 1.0;
		Reference_Values();
	}
	else
	{
		P_ref = 1.0;
		Rho_ref = 1.4;
		Inlet_Mach_No = 1.65;
	}

	cout << " Initial Data for Solution\t" << endl;
	cout << "Initialize from a file or from Zero, Enter 1 to read data from file:\t";

	cin >> Initialize_Type;

	if (Initialize_Type == 1)
	{
		Initialize(Solution_File);
		cout << "Evaluating Skin Friction Coefficient" << endl;
		Evaluate_Wall_Skin_Friction();
		cout << "Writing Skin Friction Coefficient\t" << CF_File << endl;
		Write_CF_File(CF_File);
		//		exit(0);
	}
	else
	{
		Initialize(Test_Case);
		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < Total_No_Cells; index++)
		{
			//				cout<<index<<endl;
			Pressure_Static_Inlet = 2.12 * P_ref;
			Rho_Static_Inlet = (2.12 / 1.25469387) * Rho_ref;
			//				Pressure_Static_Inlet = 1.01761659*P_ref;
			//				Rho_Static_Inlet=1.68946061*Rho_ref;

			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;

			V_Magnitude_Inlet = 0.5 * (V_1 * V_1 + V_2 * V_2);

			V[1] = V_2;
			V[0] = V_1;
			//				cout<<V_1<<"\t"<<V_2<<endl;

			Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
			//		Print(Global_U);
			for (unsigned int j = 0; j < Global_U.size(); j++)
				U_Cells[index][j] = Global_U[j];
			Calculate_Primitive_Variables(index, U_Cells[index]);
			for (unsigned int j = 0; j < Global_Primitive.size(); j++)
				Primitive_Cells[index][j] = Global_Primitive[j];
		}
	}

	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);
}
