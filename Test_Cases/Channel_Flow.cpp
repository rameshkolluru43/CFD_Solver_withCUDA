#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Channel_Flow()
{

	Grid_File = "../Grid_Files/Channel_Grid_Files/Channle_2D_SG_Alpha03_141_141.txt";
	Grid_Vtk_File = "../Grid_Files/Channel_Grid_Files/Channel_2D_SG_Alpha03_141_141.vtk";

	Form_Cells(Grid_File);

	V_D V(2, 0.0);

	Directory_Name();
	File_Name();

	InitialCondition initcond;
	readInitialConditions(InitCondFileName, initcond);

	if (Initialize_Type == 1)
		Initialize(Solution_File);
	else
	{
		Initialize(Test_Case);

		cout << "Initalizing Data with inlet conditions \t" << endl;
		for (int index = 0; index < Total_No_Cells; index++)
		{
			Pressure_Static_Inlet = initCond.P;
			Rho_Static_Inlet = initCond.Rho;
			V_1 = initCond.u;
			V_2 = initCond.v;
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
}
