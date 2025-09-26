#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Ramp_15_Degree()
{

	Directory_Name();
	File_Name();
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
	{
		std::cout << "Current working directory: " << cwd << std::endl;
	}
	else
	{
		perror("getcwd() error");
	}
	if (Is_Viscous_Wall)
	{
		Re = 386551.2;
		Inv_Re = 1.0 / Re;
		Pr = 0.72;
		Inv_Pr = 1.0 / Pr;
		M_ref = 11.63;
		K1 = 1.0 / ((gamma - 1.0) * M_ref * M_ref * Re * Pr);
		L_ref = 0.7;
		Inlet_Mach_No = 11.63;
		P_ref = 1.0 / (gamma * M_ref * M_ref);
		Rho_ref = 1.0;
		Reference_Values();
		Viscous_Time_Case = 2;
	}
	else
	{
		P_ref = 1.0 / 1.4;
		Rho_ref = 1.0;
		Inlet_Mach_No = 2.0;
	}
	std::string pwd = std::filesystem::current_path();
	// pwd = pwd + "/QPDE_FVM_CFD_Solver"	;
	std::string grid_dir = pwd + "/Grid_Files";
	std::cout << grid_dir << std::endl;

	double C = 0.0;
	switch (Grid_Size)
	{
	case 1:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_241_81.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_241_81.vtk";
		Error_File += "_WithTP_M2_241_81.txt";
		Initial_Solution_File += "_WithTP_M2_241_81.txt";
		Solution_File += "_WithTP_M2_241_81.txt";
		Final_Solution_File += "_WithTP_M2_241_81.vtk";
		break;
	case 2:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_121_41.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_121_41.vtk";
		Error_File += "_M2_121_41.txt";
		Initial_Solution_File += "_M2_121_41.txt";
		Solution_File += "_M2_121_41.txt";
		Final_Solution_File += "_M2_121_41.vtk";
		break;
	case 8:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_61_21.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_61_21.vtk";
		Error_File += "_M2_61_21.txt";
		Initial_Solution_File += "_M2_61_21.txt";
		Solution_File += "_M2_61_21.txt";
		Final_Solution_File += "_M2_61_21.vtk";
	case 9:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_109_37.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_109_37.vtk";
		Error_File += "_M2_109_37.txt";
		Initial_Solution_File += "_M2_109_37.txt";
		Solution_File += "_M2_109_37.txt";
		Final_Solution_File += "_M2_109_37.vtk";
		break;
	case 10:
		Grid_File = grid_dir + "/Ramp_Grid_Files/Ramp_15o_52_18.txt";
		Grid_Vtk_File = grid_dir + "/Ramp_Grid_Files/Ramp_15o_52_18.vtk";
		Error_File += "_M2_52_18.txt";
		Initial_Solution_File += "_M2_52_18.txt";
		Solution_File += "_M2_52_18.txt";
		Final_Solution_File += "_M2_52_18.vtk";
		break;
	case 3:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_481_161.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_481_161.vtk";
		Error_File += "_M2_481_161.txt";
		Initial_Solution_File += "_M2_481_161.txt";
		Solution_File += "_M2_481_161.txt";
		Final_Solution_File += "_M2_481_161.vtk";
		break;
	case 4:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_361_121.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Ramp_15o_361_121.vtk";
		Error_File += "_M2_361_121.txt";
		Initial_Solution_File += "_M2_361_121.txt";
		Solution_File += "_M2_361_121.txt";
		Final_Solution_File += "_M2_361_121.vtk";
		break;
	case 5:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Viscous_Ramp_15o_801_101.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Viscous_Ramp_15o_801_101.vtk";
		Error_File += "_M11_63_801_101.txt";
		Initial_Solution_File += "_M11_63_801_101.txt";
		Solution_File += "_M11_63_801_101.txt";
		Final_Solution_File += "_M11_63_801_101.vtk";
		CF_File += "_150_801_101.txt";
		break;
	case 6:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Viscous_Ramp_15o_101_301.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Viscous_Ramp_15o_101_301.vtk";
		Error_File += "_M11_63_101_301.txt";
		Initial_Solution_File += "_M11_63_101_301.txt";
		Solution_File += "_M11_63_101_301.txt";
		Final_Solution_File += "_M11_63_101_301.vtk";
		CF_File += "_150_101_301.txt";
		break;
	case 7:
		Grid_File = "../Grid_Files/Ramp_Grid_Files/Viscous_Ramp_15o_401_51.txt";
		Grid_Vtk_File = "../Grid_Files/Ramp_Grid_Files/Viscous_Ramp_15o_401_51.vtk";
		Error_File += "_M11_63_401_51.txt";
		Initial_Solution_File += "_M11_63_401_51.txt";
		Solution_File += "_M11_63_401_51.txt";
		Final_Solution_File += "_M11_63_401_51.vtk";
		CF_File += "_150_401_51.txt";
		break;
	}

	/* Reads Input grid file and does the preprocessing required for grid
	 * Calculates the Cell Normals, Cell Areas
	 * Checks for Grid
	 */
	Form_Cells(Grid_File);
	cout << "Grid_Type used \t" << Grid_Type << endl;
	Write_Cell_Info("../Grid_Files/Ramp_Grid_Files/Ramp_Grid_Details_481_161.txt");
	V_D V(2, 0.0);

	cout << "Initialize Type\t" << Initialize_Type << endl;
	if (Initialize_Type == 1)
	{

		Initialize(Solution_File);
		//		cout<<"Evaluating Skin Friction Coefficient"<<endl;
		//		Evaluate_Wall_Skin_Friction();
		//		cout<<"Writing Skin Friction Coefficient\t"<<CF_File<<endl;
		//		Write_CF_File(CF_File);
		//		cout<<"Writing Done"<<endl;
		//		exit(0);
	}
	else
	{
		Initialize(Test_Case);

		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < Total_No_Cells; index++)
		{

			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			C = sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			//                 cout<<cp<<"\t"<<cv<<"\t"<<R<<endl;
			//                 cout<<Pressure_Static_Inlet<<"\t"<<Rho_Static_Inlet<<"\t"<<Temperature_Static_Inlet<<"\t"<<C<<endl;
			V_1 = Inlet_Mach_No * C;
			V_2 = 0.0;

			V_Magnitude_Inlet = 0.5 * (V_1 * V_1 + V_2 * V_2);

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
