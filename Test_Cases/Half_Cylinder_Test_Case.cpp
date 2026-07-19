#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"
#include "Viscous_Functions.h"

extern "C" void Half_Cylinder_Flow()
{

    /*// Get current and parent directory paths
    Directory_Name();
    File_Name();

    std::filesystem::path pwd = std::filesystem::current_path();
    std::filesystem::path parent = pwd.parent_path();
    std::string basePath = parent.string();


    // Construct grid file paths based on Grid_Size
    switch (Grid_Size) {
        case 1:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_321_41.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_321_41.vtk";
            Error_File += "_M6_321_41.txt";
            Initial_Solution_File += "_M6_321_41.txt";
            Solution_File += "_M6_321_41.txt";
            Final_Solution_File += "_M6_321_41.vtk";
            break;
        case 2:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_241_81.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_241_81.vtk";
            Error_File += "_M6_241_81.txt";
            Initial_Solution_File += "_M6_241_81.txt";
            Solution_File += "_M6_241_81.txt";
            Final_Solution_File += "_M6_241_81.vtk";
            break;
        case 3:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_161_21.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_161_21.vtk";
            Error_File += "_M6_161_21.txt";
            Initial_Solution_File += "_M6_161_21.txt";
            Solution_File += "_M6_161_21.txt";
            Final_Solution_File += "_M6_161_21.vtk";
            break;
        case 4:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_121_41.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_121_41.vtk";
            Error_File += "_M6_121_41.txt";
            Initial_Solution_File += "_M6_121_41.txt";
            Solution_File += "_M6_121_41.txt";
            Final_Solution_File += "_M6_121_41.vtk";
            break;
        case 5:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_81_11.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_81_11.vtk";
            Error_File += "_M6_81_11.txt";
            Initial_Solution_File += "_M6_81_11.txt";
            Solution_File += "_M6_81_11.txt";
            Final_Solution_File += "_M6_81_11.vtk";
            break;
        case 6:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_61_21.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_61_21.vtk";
            Error_File += "_M6_61_21.txt";
            Initial_Solution_File += "_M6_61_21.txt";
            Solution_File += "_M6_61_21.txt";
            Final_Solution_File += "_M6_61_21.vtk";
            break;
        case 7:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_481_161.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_481_161.vtk";
            Error_File += "_M6_481_161.txt";
            Initial_Solution_File += "_M6_481_161.txt";
            Solution_File += "_M6_481_161.txt";
            Final_Solution_File += "_M6_481_161.vtk";
            break;
        case 8:
            Grid_File = basePath + "/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_61_21.txt";
            Grid_Vtk_File = basePath + "/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_61_21.vtk";
            Error_File += "_M6_61_21.txt";
            Initial_Solution_File += "_M6_61_21.txt";
            Solution_File += "_M6_61_21.txt";
            Final_Solution_File += "_M6_61_21.vtk";
            break;
        case 9:
                Grid_File = basePath + "/Grid_Files/Flow_Over_Bump/Single_Bump_61_21.txt";
                Grid_Vtk_File = basePath + "/Grid_Files/Flow_Over_Bump/Single_Bump_61_21.vtk";
                Error_File += "_Subsonic_61_21.txt";
                Initial_Solution_File += "_Subsonic_61_21.txt";
                Solution_File += "_Subsonic_61_21.txt";
                Final_Solution_File += "_Subsonic_61_21.vtk";
                break;
        default:
            std::cerr << "Invalid Grid_Size provided!" << std::endl;
            exit(EXIT_FAILURE);
    }

    // Log constructed file paths for debugging
    std::cout << "Grid File: " << Grid_File << std::endl;
    std::cout << "Grid VTK File: " << Grid_Vtk_File << std::endl;

    // Verify files exist
    if (!std::filesystem::exists(Grid_File)) {
        std::cerr << "Error: Grid file does not exist at path: " << Grid_File << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!std::filesystem::exists(Grid_Vtk_File)) {
        std::cerr << "Error: VTK file does not exist at path: " << Grid_Vtk_File << std::endl;
        exit(EXIT_FAILURE);
    }

    // Grid preprocessing*/

    cout << "Checking and creating the directories for solution files" << endl;
    createOutputDirectories();
    // Grid_File = gridFiles[0];
    cout << "Grid file to be read" << Grid_File << endl;
    Form_Cells(Grid_File);

    std::cout << "Grid_Type used: " << Grid_Type << std::endl;
    std::cout << Initialize_Type << std::endl;

    readInitialConditions(InitCondFileName, initCond);

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

    Identify_Wall_Boundary_Faces(Grid_Type);

    /* Viscous heat-transfer stack: Re/Pr from JSON (or defaults) + nondim heat-flux scale K1. */
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
        cout << "Half-cylinder viscous refs: Re=" << Re << " Pr=" << Pr
             << " M_ref=" << M_ref << " K1=" << K1
             << " Tw=" << wallCond.T << endl;
    }

    Write_Solution(Initial_Solution_File, 1);
    cout << "Initialized Solution, Identified Boundaries... Ready to solve." << std::endl;
}
