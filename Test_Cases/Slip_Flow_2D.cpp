#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Slip_Flow_2D()
{
	Directory_Name();
	File_Name();
	std::string pwd = std::filesystem::current_path();
	// pwd = pwd + "/QPDE_FVM_CFD_Solver"	;
	std::string grid_dir = pwd + "/Grid_Files";
	std::cout << grid_dir << std::endl;

	switch (Grid_Size)
	{
	case 1:
		Grid_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_129_129.txt";
		Grid_Vtk_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_129_129.vtk";
		Error_File += "_NEF_129_129.txt";
		Initial_Solution_File += "_NEF_129_129.txt";
		Solution_File += "_NEF_129_129.txt";
		Final_Solution_File += "_NEF_129_129.vtk";
		break;
	case 2:
		Grid_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_65_65.txt";
		Grid_Vtk_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_65_65.vtk";
		Error_File += "_NEF_65_65.txt";
		Initial_Solution_File += "_NEF_65_65.txt";
		Solution_File += "_NEF_65_65.txt";
		Final_Solution_File += "_NEF_65_65.vtk";
		break;
	case 3:
		Grid_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_33_33.txt";
		Grid_Vtk_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_33_33.vtk";
		Error_File += "_NEF_33_33.txt";
		Initial_Solution_File += "_NEF_33_33.txt";
		Solution_File += "_NEF_33_33.txt";
		Final_Solution_File += "_NEF_33_33.vtk";
		break;
	case 4:
		Grid_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_17_17.txt";
		Grid_Vtk_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_17_17.vtk";
		Error_File += "_NEF_17_17.txt";
		Initial_Solution_File += "_NEF_17_17.txt";
		Solution_File += "_NEF_17_17.txt";
		Final_Solution_File += "_NEF_17_17.vtk";
		break;
	case 5:
		Grid_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_9_9.txt";
		Grid_Vtk_File = grid_dir + "/Channel_Grid_Files/Slip_Flow_9_9.vtk";
		Error_File += "_NEF_9_9.txt";
		Initial_Solution_File += "_NEF_9_9.txt";
		Solution_File += "_NEF_9_9.txt";
		Final_Solution_File += "_NEF_9_9.vtk";
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

	Is_Time_Dependent = false;
	Terminating_Time = 100000.0;
	cin >> Initialize_Type;

	if (Initialize_Type == 1)
		Initialize(Initial_Solution_File);
	else
	{
		Initialize(Test_Case);
		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < No_Physical_Cells; index++)
		{
			Pressure_Static_Inlet = 1.0;
			Rho_Static_Inlet = 1.4;
			if (index < 0.5 * No_Physical_Cells)
				Inlet_Mach_No = 2.0;
			else
				Inlet_Mach_No = 3.0;

			V_1 = Inlet_Mach_No;
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

	cout << "Writing Initial Solution to file " << endl;
	Write_Solution(Initial_Solution_File, 1);

	Read_Write_Grid(Grid_Vtk_File, Final_Solution_File);
	Append_Solution(Initial_Solution_File, Final_Solution_File);
	cout << "Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve" << endl;
}
