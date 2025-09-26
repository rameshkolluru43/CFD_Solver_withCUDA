#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"
#include "Primitive_Computational.h"
#include "IO_Write.h"
#include "Grid.h"
#include "Directory_Files.h"
#include "Initialize.h"
#include "Directory_Files.h"

void Initialize_TestCase()
{

    cout << "Checking and creating the directories for solution files" << endl;
    createOutputDirectories();
    // Grid_File = gridFiles[0];
    cout << "Grid file to be read" << Grid_File << endl;
    Form_Cells(Grid_File);

    std::cout << "Grid_Type used: " << Grid_Type << std::endl;
    std::cout << Initialize_Type << std::endl;

    // Initialize solution
    if (Initialize_Type == 1)
    {
        Initialize(Solution_File);
    }
    else
    {
        Initialize(Test_Case);
    }
    V_D V(2, 0.0);
    double a = 0.0;
    // Set initial conditions
    // Print Initial Conditions read from json file
    cout << "Inlet Conditions: " << endl;
    cout << "Pressure_Static_Inlet: " << initCond.P << endl;
    cout << "Rho_Static_Inlet: " << initCond.Rho << endl;
    cout << "Inlet_Mach_No: " << initCond.M << endl;
    cout << "V_1: " << initCond.u << endl;
    cout << "V_2: " << initCond.v << endl;
    cout << "Temperature_Static_Inlet: " << initCond.T << endl;

    for (int index = 0; index < Total_No_Cells; ++index)
    {
        Pressure_Static_Inlet = initCond.P;
        Inlet_Mach_No = initCond.M;
        Rho_Static_Inlet = initCond.Rho;
        Temperature_Static_Inlet = initCond.T;
        a = sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);

        V[0] = initCond.u;
        V[1] = initCond.v;

        Calculate_Computational_Variables(Pressure_Static_Inlet, V, Rho_Static_Inlet, 2);
        for (unsigned int j = 0; j < Global_U.size(); ++j)
        {
            U_Cells[index][j] = Global_U[j];
        }

        Calculate_Primitive_Variables(index, U_Cells[index]);
        for (unsigned int j = 0; j < Global_Primitive.size(); ++j)
        {
            Primitive_Cells[index][j] = Global_Primitive[j];
        }
    }

    // Identify_Wall_Boundary_Faces(Grid_Type);
    Write_Solution(Initial_Solution_File, 1);
    cout << "Initialized Solution, Identified Boundaries... Ready to solve." << std::endl;
}
