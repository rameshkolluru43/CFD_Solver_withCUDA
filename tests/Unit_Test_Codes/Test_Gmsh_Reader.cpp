#include "../Basic_Files/definitions.h"
InletCondition inletCond;
ExitCondition exitCond;
InitialCondition initCond;
vector<string> gridFiles;
string gridDir, Test_Case_Name = "Default_Test_Case", GridVTKFile;
int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Error: No JSON input file provided." << endl;
        exit(EXIT_FAILURE);
    }
    string GridFileName = argv[1];
    // Read the grid file

    Read_GmshVTK_Grid(GridFileName);

    //  Read_Grid(GridFileName);
    for (int i = 0; i < Cells.size(); i++)
    {
        Print(Cells[i]);
    }

    /*
    // Example: a simple square in the XY plane with an extra point in the middle.
    V_D points = {
        1.0, 0.0, 0.0,  // Point A
        0.0, 1.0, 0.0,  // Point B
        -1.0, 0.0, 0.0, // Point C
        0.0, -1.0, 0.0, // Point D
    };

    std::cout << "Original Points (x, y, z):\n";
    for (size_t i = 0; i < points.size(); i += 3)
    {
        std::cout << "(" << points[i] << ", " << points[i + 1] << ", " << points[i + 2] << ")\n";
    }

    // Sort points anti-clockwise around their centroid (only using x, y)
    Sort_Points_AntiClockWise(points);

    std::cout << "\nSorted Points (anti-clockwise around centroid):\n";
    for (size_t i = 0; i < points.size(); i += 3)
    {
        std::cout << "(" << points[i] << ", " << points[i + 1] << ", " << points[i + 2] << ")\n";
    }
    Sort_Points_AntiClockWise(points);*/
    return 0;
}