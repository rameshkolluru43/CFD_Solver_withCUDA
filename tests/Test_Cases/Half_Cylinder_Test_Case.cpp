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
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_321_41.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_321_41.vtk";
            Error_File += "_M6_321_41.txt";
            Initial_Solution_File += "_M6_321_41.txt";
            Solution_File += "_M6_321_41.txt";
            Final_Solution_File += "_M6_321_41.vtk";
            break;
        case 2:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_241_81.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_241_81.vtk";
            Error_File += "_M6_241_81.txt";
            Initial_Solution_File += "_M6_241_81.txt";
            Solution_File += "_M6_241_81.txt";
            Final_Solution_File += "_M6_241_81.vtk";
            break;
        case 3:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_161_21.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_161_21.vtk";
            Error_File += "_M6_161_21.txt";
            Initial_Solution_File += "_M6_161_21.txt";
            Solution_File += "_M6_161_21.txt";
            Final_Solution_File += "_M6_161_21.vtk";
            break;
        case 4:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_121_41.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_121_41.vtk";
            Error_File += "_M6_121_41.txt";
            Initial_Solution_File += "_M6_121_41.txt";
            Solution_File += "_M6_121_41.txt";
            Final_Solution_File += "_M6_121_41.vtk";
            break;
        case 5:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_81_11.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_81_11.vtk";
            Error_File += "_M6_81_11.txt";
            Initial_Solution_File += "_M6_81_11.txt";
            Solution_File += "_M6_81_11.txt";
            Final_Solution_File += "_M6_81_11.vtk";
            break;
        case 6:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_61_21.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_61_21.vtk";
            Error_File += "_M6_61_21.txt";
            Initial_Solution_File += "_M6_61_21.txt";
            Solution_File += "_M6_61_21.txt";
            Final_Solution_File += "_M6_61_21.vtk";
            break;
        case 7:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_481_161.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Half_Cylinder_Flow_Non_Elliptic_Grid_481_161.vtk";
            Error_File += "_M6_481_161.txt";
            Initial_Solution_File += "_M6_481_161.txt";
            Solution_File += "_M6_481_161.txt";
            Final_Solution_File += "_M6_481_161.vtk";
            break;
        case 8:
            Grid_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_61_21.txt";
            Grid_Vtk_File = basePath + "/input/Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_61_21.vtk";
            Error_File += "_M6_61_21.txt";
            Initial_Solution_File += "_M6_61_21.txt";
            Solution_File += "_M6_61_21.txt";
            Final_Solution_File += "_M6_61_21.vtk";
            break;
        case 9:
                Grid_File = basePath + "/input/Grid_Files/Flow_Over_Bump/Single_Bump_61_21.txt";
                Grid_Vtk_File = basePath + "/input/Grid_Files/Flow_Over_Bump/Single_Bump_61_21.vtk";
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

    /* Shared cold-start / restart + ghost BC finalize (same as Main). */
    Initialize_TestCase();
}
