#include "definitions.h"
#include "Globals.h"
#include "Directory_Files.h"

void Directory_Name()
{
	//	std::string pwd = std::filesystem::current_path();
	std::filesystem::path pwd = std::filesystem::current_path();
	// Move one folder back
	std::filesystem::path parent = pwd.parent_path();

	// Convert the parent path to a string
	std::string pwd1 = parent.string();

	// Print the parent directory
	//    std::cout << "Current Directory: " << pwd << std::endl;
	//    std::cout << "Parent Directory: " << pwd1 << std::endl;
	// pwd = pwd + "/QPDE_FVM_CFD_Solver";

	// cout<<"In Parent Directory\t"<<pwd1<<endl;
	Solution_File = pwd1 + "/2D_Euler_Solutions/";
	Final_Solution_File = pwd1 + "/2D_Euler_Solutions/";
	Initial_Solution_File = pwd1 + "/2D_Euler_Solutions/";
	Error_File = pwd1 + "/2D_Euler_Solutions/";
	CF_File = pwd1 + "/2D_Euler_Solutions/";

	/* Test case Numbers
	 * 	case 1: Half_Cylinder_Flow();
		case 2: Mach_Reflection();
		case 3: Ramp_15_Degree();
		case 4: Forward_Facing_Step();
		case 5: Air_Foil();
		case 6: Shock_Tube_2D();
		case 7: Slip_Flow_2D();
	*/

	switch (Test_Case)
	{
	case 1:
		Solution_File += "Half_Cylinder_Flow/";
		Initial_Solution_File += "Half_Cylinder_Flow/";
		Error_File += "Half_Cylinder_Flow/";
		Final_Solution_File += "Half_Cylinder_Flow/";
		CF_File += "Half_Cylinder_Flow/";
		break;

	case 2:
		Solution_File += "Shock_Reflection_2D/";
		Final_Solution_File += "Shock_Reflection_2D/";
		Initial_Solution_File += "Shock_Reflection_2D/";
		Error_File += "Shock_Reflection_2D/";
		break;

	case 3:
		Solution_File += "Ramp_15_Degree/";
		Final_Solution_File += "Ramp_15_Degree/";
		Initial_Solution_File += "Ramp_15_Degree/";
		Error_File += "Ramp_15_Degree/";
		CF_File += "Ramp_15_Degree/";
		break;
	case 4:
		Solution_File += "Forward_Step/";
		Final_Solution_File += "Forward_Step/";
		Initial_Solution_File += "Forward_Step/";
		Error_File += "Forward_Step/";
		break;
	case 5:
		Solution_File += "Air_Foil/";
		Final_Solution_File += "Air_Foil/";
		Initial_Solution_File += "Air_Foil/";
		Error_File += "Air_Foil/";
		break;

	case 6:
		Solution_File += "Shock_Tube_2D/";
		Final_Solution_File += "Shock_Tube_2D/";
		Initial_Solution_File += "Shock_Tube_2D/";
		Error_File += "Shock_Tube_2D/";
		break;
	case 7:
		Solution_File += "Slip_Flow_2D/";
		Final_Solution_File += "Slip_Flow_2D/";
		Initial_Solution_File += "Slip_Flow_2D/";
		Error_File += "Slip_Flow_2D/";
		break;
	case 8:
		Solution_File += "Channel_Flow_2D/";
		Final_Solution_File += "Channel_Flow_2D/";
		Initial_Solution_File += "Channel_Flow_2D/";
		Error_File += "Channel_Flow_2D/";
		break;
	case 9:
		Solution_File += "Scramjet_Inlet/";
		Final_Solution_File += "Scramjet_Inlet/";
		Initial_Solution_File += "Scramjet_Inlet/";
		Error_File += "Scramjet_Inlet/";
		break;
	case 10:
		Solution_File += "Flow_Over_Bump/";
		Final_Solution_File += "Flow_Over_Bump/";
		Initial_Solution_File += "Flow_Over_Bump/";
		Error_File += "Flow_Over_Bump/";
		break;
	case 11:
		Solution_File += "Slip_Flow_2D/";
		Final_Solution_File += "Slip_Flow_2D/";
		Initial_Solution_File += "Slip_Flow_2D/";
		Error_File += "Slip_Flow_2D/";
		break;
	case 12:
		Solution_File += "Flat_Plate_BL/";
		Final_Solution_File += "Flat_Plate_BL/";
		Initial_Solution_File += "Flat_Plate_BL/";
		Error_File += "Flat_Plate_BL/";
		break;
	case 13:
		Solution_File += "Shock_Wedge_Reflection/";
		Final_Solution_File += "Shock_Wedge_Reflection/";
		Initial_Solution_File += "Shock_Wedge_Reflection/";
		Error_File += "Shock_Wedge_Reflection/";
		break;
	case 14:
		Solution_File += "Shock_Diffraction/";
		Final_Solution_File += "Shock_Diffraction/";
		Initial_Solution_File += "Shock_Diffraction/";
		Error_File += "Shock_Diffraction/";
		break;
	case 15:
		Solution_File += "SWBLI/";
		Final_Solution_File += "SWBLI/";
		Initial_Solution_File += "SWBLI/";
		Error_File += "SWBLI/";
		CF_File += "SWBLI/";
		break;
	case 16:
		Solution_File += "Shock_Bubble_Interaction/";
		Final_Solution_File += "Shock_Bubble_Interaction/";
		Initial_Solution_File += "Shock_Bubble_Interaction/";
		Error_File += "Shock_Bubble_Interaction/";
		break;
	case 17:
		Solution_File += "Flow_Over_Bump/";
		Final_Solution_File += "Flow_Over_Bump/";
		Initial_Solution_File += "Flow_Over_Bump/";
		Error_File += "Flow_Over_Bump/";
		CF_File += "Flow_Over_Bump/";
		break;
	case 18:
		Solution_File += "OED/";
		Final_Solution_File += "OED/";
		Initial_Solution_File += "OED/";
		Error_File += "OED/";
		CF_File += "OED/";
		break;
	case 19:
		Solution_File += "Elling/";
		Final_Solution_File += "Elling/";
		Initial_Solution_File += "Elling/";
		Error_File += "Elling/";
		CF_File += "Elling/";
		break;
	}
	/* Dissipation_Type
	  1- LLF Dissipation
	  2 - PVU
	  3 - Roe
	  4 - Naveens Dissipatoin
	  5 - MOVERS
	  6 - RUSANOV
  */

	switch (Dissipation_Type)
	{
	case 1:
		Solution_File += "LLF/";
		Final_Solution_File += "LLF/";
		Initial_Solution_File += "LLF/";
		Error_File += "LLF/";
		CF_File += "LLF/";
		break;
	case 2:
		Solution_File += "MOVERS/";
		Final_Solution_File += "MOVERS/";
		Initial_Solution_File += "MOVERS/";
		Error_File += "MOVERS/";
		CF_File += "MOVERS/";
		break;
	case 3:
		Solution_File += "ROE/";
		Final_Solution_File += "ROE/";
		Initial_Solution_File += "ROE/";
		Error_File += "ROE/";
		CF_File += "ROE/";
		break;
	case 4:
		Solution_File += "RICCA/";
		Final_Solution_File += "RICCA/";
		Initial_Solution_File += "RICCA/";
		Error_File += "RICCA/";
		CF_File += "RICCA/";
		break;
	case 5:
		Solution_File += "MOVERS_NWSC/";
		Final_Solution_File += "MOVERS_NWSC/";
		Initial_Solution_File += "MOVERS_NWSC/";
		Error_File += "MOVERS_NWSC/";
		CF_File += "MOVERS_NWSC/";
		break;
	}
}

void File_Name()
{
	/*  	1- LLF Dissipation
		2 - PVU
		3 - Roe
		4 - Naveens Dissipatoin
		5 - MOVERS
		12 - RUSANOV*/
	switch (Dissipation_Type)
	{
	case 1:
		Solution_File += "Solution_LLF";
		Final_Solution_File += "Final_Solution_LLF";
		Initial_Solution_File += "Initial_Solution_LLF";
		Error_File += "Error_LLF";
		CF_File += "CF_LLF";
		break;
	case 3:
		Solution_File += "Solution_ROE";
		Final_Solution_File += "Final_Solution_ROE";
		Initial_Solution_File += "Initial_Solution_ROE";
		Error_File += "Error_ROE";
		CF_File += "CF_ROE";

		if (Enable_Entropy_Fix)
		{
			Solution_File += "_Entropy_Fix";
			Final_Solution_File += "_Entropy_Fix";
			Initial_Solution_File += "_Entropy_Fix";
			Error_File += "_Entropy_Fix";
			CF_File += "_Entropy_Fix";
		}
		break;
	case 2:
		Solution_File += "Solution_MOVERS";
		Final_Solution_File += "Final_Solution_MOVERS";
		Initial_Solution_File += "Initial_Solution_MOVERS";
		Error_File += "Error_MOVERS";
		CF_File += "CF_MOVERS";
		if (Is_MOVERS_1)
		{
			Solution_File += "_1Wave";
			Final_Solution_File += "_1Wave";
			Initial_Solution_File += "_1Wave";
			Error_File += "_1Wave";
			CF_File += "_1Wave";
		}
		else
		{
			Solution_File += "_NWave";
			Final_Solution_File += "_NWave";
			Initial_Solution_File += "_NWave";
			Error_File += "_NWave";
			CF_File += "_NWave";
		}
		if (Enable_Entropy_Fix)
		{
			Solution_File += "_Entropy_Fix";
			Final_Solution_File += "_Entropy_Fix";
			Initial_Solution_File += "_Entropy_Fix";
			Error_File += "_Entropy_Fix";
			CF_File += "_Entropy_Fix";
		}
		break;
	case 4:
		Solution_File += "Solution_RICCA";
		Final_Solution_File += "Final_Solution_RICCA";
		Initial_Solution_File += "Initial_Solution_RICCA";
		Error_File += "Error_RICCA";
		CF_File += "CF_RICCA";
		break;
	case 5:
		Solution_File += "Solution_MOVERS_NWSC";
		Final_Solution_File += "Final_Solution_MOVERS_NWSC";
		Initial_Solution_File += "Initial_Solution_MOVERS_NWSC";
		Error_File += "Error_MOVERS_NWSC";
		CF_File += "CF_MOVERS_NWSC";
		break;
	}

	if (Is_Viscous_Wall)
	{
		Solution_File += "_Viscous";
		Final_Solution_File += "_Viscous";
		Initial_Solution_File += "_Viscous";
		Error_File += "_Viscous";
		CF_File += "_Viscous";
	}
	else
	{
		Solution_File += "_Inviscid";
		Final_Solution_File += "_Inviscid";
		Initial_Solution_File += "_Inviscid";
		Error_File += "_Inviscid";
		CF_File += "_Inviscid";
	}
	if (Time_Accurate)
	{
		Solution_File += "_RK3";
		Final_Solution_File += "_RK3";
		Initial_Solution_File += "_RK3";
		Error_File += "_RK3";
		Error_File += "_RK3";
	}
	if (Local_Time_Stepping)
	{
		Solution_File += "_Local_Time_Stepping";
		Final_Solution_File += "_Local_Time_Stepping";
		Initial_Solution_File += "_Local_Time_Stepping";
		Error_File += "_Local_Time_Stepping";
		CF_File += "_Local_Time_Stepping";
	}
	else
	{
		Solution_File += "_Global_Time_Stepping";
		Final_Solution_File += "_Global_Time_Stepping";
		Initial_Solution_File += "_Global_Time_Stepping";
		Error_File += "_Global_Time_Stepping";
		CF_File += "_Global_Time_Stepping";
	}
	if (Is_Second_Order)
	{
		Solution_File += "_2Order";
		Final_Solution_File += "_2Order";
		Initial_Solution_File += "_2Order";
		Error_File += "_2Order";
		CF_File += "_2Order";
		switch (Limiter_Case)
		{
		case 0:
			Solution_File += "_MinMod_Limiter";
			Final_Solution_File += "_MinMod_Limiter";
			Initial_Solution_File += "_MinMod_Limiter";
			Error_File += "_MinMod_Limiter";
			CF_File += "_MinMod_Limiter";
			break;
		case 1:
			Solution_File += "_Superbee_Limiter";
			Final_Solution_File += "_Superbee_Limiter";
			Initial_Solution_File += "_Superbee_Limiter";
			Error_File += "_Superbee_Limiter";
			break;
		case 2:
			Solution_File += "_MCD_Limiter";
			Final_Solution_File += "_MCD_Limiter";
			Initial_Solution_File += "_MCD_Limiter";
			Error_File += "_MCD_Limiter";
			break;
		case 3:
			Solution_File += "_vanLeerSmooth_Limiter";
			Final_Solution_File += "_vanLeerSmooth_Limiter";
			Initial_Solution_File += "_vanLeerSmooth_Limiter";
			Error_File += "_vanLeerSmooth_Limiter";
			break;
		case 4:
			Solution_File += "_Log_Limiter";
			Final_Solution_File += "_Log_Limiter";
			Initial_Solution_File += "_Log_Limiter";
			Error_File += "_Log_Limiter";
			break;
		case 5:
			Solution_File += "_VenkatKrishnan_Limiter";
			Final_Solution_File += "_VenkatKrishnan_Limiter";
			Initial_Solution_File += "_VenkatKrishnan_Limiter";
			Error_File += "_VenkatKrishnan_Limiter";
			break;
		case 6:
			Solution_File += "_Scaled_Reconstruction";
			Final_Solution_File += "_Scaled_Reconstruction";
			Initial_Solution_File += "_Scaled_Reconstruction";
			Error_File += "_Scaled_Reconstruction";
			CF_File += "_Scaled_Reconstruction";
			break;
		}
	}
	else if (Is_WENO)
	{
		Solution_File += "_WENO";
		Final_Solution_File += "_WENO";
		Initial_Solution_File += "_WENO";
		Error_File += "_WENO";
		CF_File += "_WENO";
	}
	else
	{
		Solution_File += "_1Order";
		Final_Solution_File += "_1Order";
		Initial_Solution_File += "_1Order";
		Error_File += "_1Order";
		CF_File += "_1Order";
	}
}

void Write_VTK_File(const string &Op_file1, const string &Op_file2)
{
	ofstream outputfile1(Op_file1.c_str(), ios::out);
	ofstream outputfile2(Op_file2.c_str(), ios::out);
	int Cell_Index = 0, Cells_In_Plane_Core, Cells_In_Plane_Polar;
	Cells_In_Plane_Core = (nx_c - 1) * (ny_c - 1);
	Cells_In_Plane_Polar = (nx_p - 1) * (ny_p - 1);
	if ((outputfile1.is_open()) && (outputfile2.is_open()))
	{
		for (int k = 0; k < (nz_c - 1); k++)
		{
			for (int j = 0; j < (ny_c - 1); j++)
			{
				for (int i = 0; i < (nx_c - 1); i++)
				{
					Error_File += "Error_Inviscid_LLF";
					Cell_Index = i + j * (nx_c - 1) + k * (Cells_In_Plane_Core + Cells_In_Plane_Polar);
					// outputfile1<<Cells_Cell_Center[Cell_Index][0]<<"\t"<<Cells_Cell_Center[Cell_Index][1]<<"\t"<<Cells_Cell_Center[Cell_Index][2]<<endl;
					outputfile1 << Cells[Cell_Index].Cell_Center[0] << "\t" << Cells[Cell_Index].Cell_Center[1] << "\t" << Cells[Cell_Index].Cell_Center[2] << endl;
				}
			}
			for (int j = 0; j < (ny_p - 1); j++)
			{
				for (int i = 0; i < (nx_p - 1); i++)
				{
					Cell_Index = i + j * (nx_p - 1) + k * Cells_In_Plane_Polar + (k + 1) * Cells_In_Plane_Core;
					// outputfile2 << Cells_Cell_Center[Cell_Index][0] << "\t" << Cells_Cell_Center[Cell_Index][1] << "\t" << Cells_Cell_Center[Cell_Index][2] << endl;
					outputfile2 << Cells[Cell_Index].Cell_Center[0] << "\t" << Cells[Cell_Index].Cell_Center[1] << "\t" << Cells[Cell_Index].Cell_Center[2] << endl;
				}
			}
		}
	}
	else
	{
		cout << "Could not Open file for writing VTK file" << endl;
		exit(0);
	}
}

void Write_VTK_File(const string &File_Name)
{
	ofstream Outfile(File_Name.c_str(), ios::out);
	if (Outfile.is_open())
	{

		for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
			// Outfile << Cells_Cell_Center[Cell_Index][0] << "\t" << Cells_Cell_Center[Cell_Index][1] << "\t" << Cells_Cell_Center[Cell_Index][2] << "\t" << Primitive_Cells[Cell_Index][1] << "\t" << Primitive_Cells[Cell_Index][2] << "\t" << Primitive_Cells[Cell_Index][3] << "\t" << Primitive_Cells[Cell_Index][10] << endl;
			Outfile << Cells[Cell_Index].Cell_Center[0] << "\t" << Cells[Cell_Index].Cell_Center[1] << "\t" << Cells[Cell_Index].Cell_Center[2] << "\t" << Primitive_Cells[Cell_Index][1] << "\t" << Primitive_Cells[Cell_Index][2] << "\t" << Primitive_Cells[Cell_Index][3] << "\t" << Primitive_Cells[Cell_Index][10] << endl;
	}
	else
	{
		cout << "Could not write Data into VTK file\n";
		cout << File_Name << endl;
		exit(0);
	}
}

void Write_CF_File(const string &File_Name)
{
	ofstream Outfile(File_Name.c_str(), ios::out);
	V_D Cell_Center(3, 0);
	//	cout<<File_Name<<endl;
	int Cell_Index, j;
	if (Outfile.is_open())
	{
		j = 0;
		for (unsigned int i = 0; i < CF.size(); i++)
		{
			Cell_Index = Wall_Cells_List[j];
			j += 3;
			//			cout<<Cell_Index<<endl;
			//			Print(Cell_Center);
			// Cell_Center = Cells_Cell_Center[Cell_Index];
			Cell_Center = Cells[Cell_Index].Cell_Center;
			Outfile << i << "\t" << Cell_Center[0] << '\t' << CF[i] << endl;
		}
		//		cout<<"Writing IN CF FIle done"<<endl;
	}
	else
	{
		cout << "Could not write Data into CF file\n";
		cout << File_Name << endl;
		exit(0);
	}
	Outfile.close();
}

void Write_Solution(const string &Op_file, const int &type)
{
	//     string myfile = Solution_File;
	ofstream outputfile(Op_file.c_str(), ios::out);

	//           cout<<Op_file<<endl;

	if (outputfile.is_open())
	{
		// 	  cout<<"Opened Output File\n";
		if (Grid_Type == 0)
			outputfile << nx_c << "\t" << ny_c << "\t" << nz_c << endl;
		else if (Grid_Type == 1)
		{
			outputfile << nx_1 << "\t" << ny_1 << "\t" << nx_2 << "\t" << ny_2 << endl;
		}
		else
		{
			outputfile << nx_c << "\t" << ny_c << "\t" << nz_c << endl;
			outputfile << nx_p << "\t" << ny_p << "\t" << nz_p << endl;
		}
		outputfile << No_Physical_Cells << endl;
		outputfile << iterations << endl;
		for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
		{
			switch (type)
			{
			case 1: // Time step, density,Pressure,Temperature,u,v
				outputfile << Cells[Cell_Index].cellID << "\t\t" << Cells[Cell_Index].del_t << "\t" << Primitive_Cells[Cell_Index][0] << "\t" << Primitive_Cells[Cell_Index][4] << "\t" << Primitive_Cells[Cell_Index][3] << "\t\t\t\t" << Primitive_Cells[Cell_Index][1] << "\t\t\t\t" << Primitive_Cells[Cell_Index][2] << "\t\t\t\t" << Primitive_Cells[Cell_Index][7] << "\t\t\t\t" << Primitive_Cells[Cell_Index][10] << endl;
				break;
			case 2:
				// 					outputfile<<Cells_U_Grad[Cell_Index][0]<<"\t"<<Cells_U_Grad[Cell_Index][1]<<"\t"<<Cells_U_Grad[Cell_Index][2]<<endl;
				break;
			case 3:
				// 					outputfile<<Cells_V_Grad[Cell_Index][0]<<"\t"<<Cells_V_Grad[Cell_Index][1]<<"\t"<<Cells_V_Grad[Cell_Index][2]<<endl;
				break;
			case 4:
				// 					outputfile<<Cells_W_Grad[Cell_Index][0]<<"\t"<<Cells_W_Grad[Cell_Index][1]<<"\t"<<Cells_W_Grad[Cell_Index][2]<<endl;
				break;
			case 5:
				// 					outputfile<<Cells_T_Grad[Cell_Index][0]<<"\t"<<Cells_T_Grad[Cell_Index][1]<<"\t"<<Cells_T_Grad[Cell_Index][2]<<endl;
				break;
			case 6:
				// 					outputfile<<Cells_Viscous_Flux[Cell_Index][0]<<"\t"<<Cells_Viscous_Flux[Cell_Index][1]<<"\t"<<Cells_Viscous_Flux[Cell_Index][2]<<"\t"<<Cells_Viscous_Flux[Cell_Index][3]<<"\t"<<Cells_Viscous_Flux[Cell_Index][4]<<endl;
				break;
			case 7:
				// 					outputfile<<Convective_Flux_Cells[Cell_Index][0]<<"\t"<<Convective_Flux_Cells[Cell_Index][1]<<"\t"<<Convective_Flux_Cells[Cell_Index][2]<<"\t"<<Convective_Flux_Cells[Cell_Index][3]<<"\t"<<Convective_Flux_Cells[Cell_Index][4]<<endl;
				break;
			case 8:
				outputfile << Cells_Net_Flux[Cell_Index][0] << "\t" << Cells_Net_Flux[Cell_Index][1] << "\t" << Cells_Net_Flux[Cell_Index][2] << "\t" << Cells_Net_Flux[Cell_Index][3] << "\t" << Cells_Net_Flux[Cell_Index][4] << endl;
				break;
			case 9:
				outputfile << U_Cells[Cell_Index][0] << endl;
				break;
			}
		}
	}
	else
	{
		cout << "Could not write to solution file opfile\n";
		cout << Op_file << endl;
		exit(0);
	}
	outputfile.close();
}

void Write_Error_File(const string &File_Name)
{
	//	cout<<File_Name<<endl;
	ofstream outputfile(File_Name.c_str(), ios::out | ios::app);
	if (outputfile.is_open())
	{
		for (int i = 0; i < 4; i++)
			outputfile << Error[i] << "\t";
		outputfile << endl;
	}
	else
	{
		cout << "Could not write to  Error file\n";
		cout << File_Name << endl;
		exit(0);
	}
}

// Function to write a matrix to a file
void Write_A_MatrixToFile(vector<V_D> &A, const string &File_Name)
{
	// Open file in write mode
	ofstream Outfile(File_Name.c_str(), ios::out);

	// Check if the file is open
	if (!Outfile.is_open())
	{
		cerr << "Error opening file: " << File_Name << endl;
		return;
	}

	// Loop through the matrix and write each row to the file
	for (unsigned int i = 0; i < A.size(); i++)
	{
		for (unsigned int j = 0; j < A[i].size(); j++)
		{
			Outfile << A[i][j] << "\t";
		}
		Outfile << endl;
	}
	// Close the file
	Outfile.close();

	cout << "Matrix written to " << File_Name << " successfully!" << endl;
}

// Function to write a matrix to a file
void Write_b_VectorToFile(V_D &b, const string &File_Name)
{
	// Open file in write mode
	ofstream Outfile(File_Name.c_str(), ios::out);

	// Check if the file is open
	if (!Outfile.is_open())
	{
		cerr << "Error opening file: " << File_Name << endl;
		return;
	}

	for (unsigned int i = 0; i < b.size(); i++)
	{
		Outfile << b[i] << endl;
	}
	// Close the file
	Outfile.close();

	cout << "b vector written to " << File_Name << " successfully!" << endl;
}

string Get_ErrorFileName()
{
	return Error_File;
}

string Get_Final_Solution_FileName()
{
	return Final_Solution_File;
}

string Get_Initial_Solution_FileName()
{
	return Initial_Solution_File;
}

string Get_Grid_VtkFile()
{
	return Grid_Vtk_File;
}

string Get_SolutionFile()
{
	return Solution_File;
}

void CreateTestCaseDirectories(int &Test_Case)
{

	// Define paths
	filesystem::path execDirectory = "./";
	filesystem::path basePath = "../";
	filesystem::path solutionBasePath = basePath / "2D_Euler_Solutions";

	// Map test cases to folder names
	string folderName;
	switch (Test_Case)
	{
	case 1:
		folderName = "Half_Cylinder_Flow";
		break;
	case 2:
		folderName = "Shock_Reflection_2D";
		break;
	case 3:
		folderName = "Ramp_15_Degree";
		break;
	case 4:
		folderName = "Forward_Facing_Step";
		break;
	case 5:
		folderName = "Air_Foil";
		break;
	case 6:
		folderName = "Shock_Tube_2D";
		break;
	case 7:
		folderName = "Slip_Flow_2D";
		break;
	default:
		cerr << "Error: Invalid test case number." << endl;
		return;
	}

	// Full path for the test case
	filesystem::path casePath = solutionBasePath / folderName;

	// Check if the test case folder already exists
	if (filesystem::exists(casePath))
	{
		cout << "Directory for test case '" << folderName
			 << "' already exists at: " << casePath << endl;
		return;
	}

	// Create directories for the test case
	try
	{
		filesystem::create_directories(casePath / "Initial_Solution");
		filesystem::create_directories(casePath / "Final_Solution");
		filesystem::create_directories(casePath / "Error_File");
		filesystem::create_directories(casePath / "CF_File");

		cout << "Directories created successfully for test case: " << folderName << endl;
	}
	catch (const filesystem::filesystem_error &e)
	{
		cerr << "Error creating directories: " << e.what() << endl;
	}
}

void Write_Cell_Info(const string &opfile)
{
	ofstream File_Cell_Info(opfile.c_str(), ios::out);
	cout << "Writing Cell Information\t" << endl;
	cout << opfile << endl;
	int k, Cell_No = 0;
	if (File_Cell_Info.is_open())
	{
		File_Cell_Info << "Cell_No\tVolume\tArea\tNormal\tNeighbour" << endl;
		for (Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
		{
			File_Cell_Info << Cell_No << "\t";
			File_Cell_Info << Cells[Cell_No].Inv_Area << "\t";
			File_Cell_Info << Cells[Cell_No].Area << "\t";
			File_Cell_Info << Cells[Cell_No].Face_Normals[0] << "\t" << Cells[Cell_No].Face_Normals[1] << "\n";
			File_Cell_Info << Cells[Cell_No].Face_Normals[2] << "\t" << Cells[Cell_No].Face_Normals[3] << "\n";
			File_Cell_Info << Cells[Cell_No].Face_Normals[4] << "\t" << Cells[Cell_No].Face_Normals[5] << "\n";
			File_Cell_Info << Cells[Cell_No].Face_Normals[6] << "\t" << Cells[Cell_No].Face_Normals[7] << "\n";

			File_Cell_Info << Cells[Cell_No].Face_Areas[0] << "\t" << Cells[Cell_No].Face_Areas[1] << "\n";
			File_Cell_Info << Cells[Cell_No].Face_Areas[2] << "\t" << Cells[Cell_No].Face_Areas[3] << "\n";

			for (k = 0; k < 4; k++)
				File_Cell_Info << Cells[Cell_No].Neighbours[k] << "\t";
			File_Cell_Info << "\n";
		}
	}
	else
	{
		cout << "Could not write to Cell information to file\n";
	}
	cout << "Cell Information written to file\n";
	cout << opfile << endl;
	cout << "--------------------------------------------------------\n";
}

void Write_Cell_Info(const int &Cell_No)
{
	int nFaces = (Cells[Cell_No].numFaces > 0) ? Cells[Cell_No].numFaces : static_cast<int>(Cells[Cell_No].Face_Areas.size());
	cout << Cell_No << endl;
	for (int Face_No = 0; Face_No < nFaces; Face_No++)
	{
		double area = (Face_No < static_cast<int>(Cells[Cell_No].Face_Areas.size())) ? Cells[Cell_No].Face_Areas[Face_No] : 0.0;
		double n1 = (Face_No * 2 + 0 < static_cast<int>(Cells[Cell_No].Face_Normals.size())) ? Cells[Cell_No].Face_Normals[Face_No * 2 + 0] : 0.0;
		double n2 = (Face_No * 2 + 1 < static_cast<int>(Cells[Cell_No].Face_Normals.size())) ? Cells[Cell_No].Face_Normals[Face_No * 2 + 1] : 0.0;
		cout << Face_No << "\t area=" << area << "\t n=" << n1 << "\t" << n2 << endl;
	}
}
