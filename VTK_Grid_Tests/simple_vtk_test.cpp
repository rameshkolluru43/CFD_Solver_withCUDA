#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

/**
 * @brief Simple standalone VTK file reader test
 *
 * This program tests basic VTK file reading without external dependencies
 * to verify the Ramp_15o_52_18.vtk file can be opened and parsed.
 */

void readVTKFileBasic(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "❌ ERROR: Cannot open file: " << filename << std::endl;
        return;
    }

    std::cout << "✅ SUCCESS: File opened successfully!" << std::endl;

    std::string line;
    int lineCount = 0;
    int pointCount = 0;
    int cellCount = 0;
    bool foundPoints = false, foundCells = false;

    while (std::getline(file, line) && lineCount < 1000)
    { // Read first 1000 lines for inspection
        lineCount++;
        std::stringstream ss(line);
        std::string keyword;
        ss >> keyword;

        if (keyword == "POINTS")
        {
            ss >> pointCount;
            foundPoints = true;
            std::cout << "Found POINTS section: " << pointCount << " points" << std::endl;
        }
        else if (keyword == "CELLS")
        {
            int numEntries;
            ss >> cellCount >> numEntries;
            foundCells = true;
            std::cout << "Found CELLS section: " << cellCount << " cells, " << numEntries << " entries" << std::endl;
        }
        else if (lineCount <= 10)
        {
            std::cout << "Line " << lineCount << ": " << line << std::endl;
        }
    }

    file.close();

    std::cout << "\nFile Analysis:" << std::endl;
    std::cout << "  - Total lines read: " << lineCount << std::endl;
    std::cout << "  - Points section found: " << (foundPoints ? "Yes" : "No") << std::endl;
    std::cout << "  - Cells section found: " << (foundCells ? "Yes" : "No") << std::endl;

    if (foundPoints && foundCells)
    {
        std::cout << "  - Points: " << pointCount << std::endl;
        std::cout << "  - Cells: " << cellCount << std::endl;
        std::cout << "\n🎉 VTK file appears to be valid and readable!" << std::endl;
    }
    else
    {
        std::cout << "\n⚠️ WARNING: VTK file may be incomplete or invalid" << std::endl;
    }
}

int main()
{
    std::cout << "=== Standalone VTK File Reader Test ===" << std::endl;
    std::cout << "Testing file: Ramp_15o_52_18.vtk" << std::endl;
    std::cout << "=======================================" << std::endl;

    std::string vtk_file_path = "/Users/rameshkolluru/My_Research/CFD_Solver_withCUDA/Grid_Files/Ramp_Grid_Files/Ramp_15o_52_18.vtk";

    // Test file existence first
    std::ifstream test_file(vtk_file_path);
    if (!test_file.good())
    {
        std::cout << "❌ ERROR: VTK file not found at: " << vtk_file_path << std::endl;
        std::cout << "Please check the file path and ensure the file exists." << std::endl;
        return 1;
    }
    test_file.close();

    std::cout << "File exists at: " << vtk_file_path << std::endl;
    std::cout << "\nReading VTK file..." << std::endl;

    readVTKFileBasic(vtk_file_path);

    return 0;
}