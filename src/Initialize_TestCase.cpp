#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"
#include "Primitive_Computational.h"
#include "IO_Write.h"
#include "Grid.h"
#include "Directory_Files.h"
#include "Initialize.h"
#include "Directory_Files.h"
#include "MPI_Utils.h"
#include "Viscous_Functions.h"

void Initialize_TestCase()
{
    if (CFD_MPI_Is_Root())
        std::cout << "Initializing the Test Case: " << Test_Case << " - " << Test_Case_Name << std::endl;

    /* Set Re/Pr/R_ref/K1 before primitive recovery so Sutherland μ is finite. */
    if (Is_Viscous || Is_Viscous_Wall)
    {
        if (!(Re > 0.0))
            Re = 1.0e5;
        if (!(Pr > 0.0))
            Pr = 0.72;
        Inv_Re = 1.0 / Re;
        Inv_Pr = 1.0 / Pr;
        L_ref = 1.0;
        M_ref = (inletCond.M > 0.0) ? inletCond.M : 6.0;
        Rho_ref = (inletCond.Rho > 0.0) ? inletCond.Rho : 1.4;
        Reference_Values();
        K1 = 1.0 / ((gamma - 1.0) * M_ref * M_ref * Re * Pr);
        Viscous_Time_Case = 2;
        if (CFD_MPI_Is_Root())
            cout << "Viscous refs: Re=" << Re << " Pr=" << Pr
                 << " M_ref=" << M_ref << " K1=" << K1
                 << " Tw=" << wallCond.T << endl;
    }

    // Initialize solution
    if (Initialize_Type == 1)
    {
        Initialize(Solution_File);
        /* Restart: keep field from Solution_File — do not overwrite with initCond. */
    }
    else
    {
        Initialize(Test_Case);
        V_D V(2, 0.0);
        for (int index = 0; index < Total_No_Cells; ++index)
        {
            Pressure_Static_Inlet = initCond.P;
            Inlet_Mach_No = initCond.M;
            Rho_Static_Inlet = initCond.Rho;
            Temperature_Static_Inlet = initCond.T;

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
    }

    if (CFD_MPI_Is_Root())
    {
        cout << "Inlet Conditions: " << endl;
        cout << "Pressure_Static_Inlet: " << inletCond.P << endl;
        cout << "Rho_Static_Inlet: " << inletCond.Rho << endl;
        cout << "Inlet_Mach_No: " << inletCond.M << endl;
        cout << "V_1: " << inletCond.u << endl;
        cout << "V_2: " << inletCond.v << endl;
        cout << "Temperature_Static_Inlet: " << inletCond.T << endl;
    }

    Write_Solution(Initial_Solution_File, 1);
    CFD_MPI_Barrier();
    if (CFD_MPI_Is_Root())
        cout << "Initialized Solution, Identified Boundaries... Ready to solve." << std::endl;
}
