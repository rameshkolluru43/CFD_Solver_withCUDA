#include "definitions.h"
#include "Globals.h"
#include "Grid.h"
#include "Directory_Files.h"
#include <iostream>
#include <string>

/**
 * @brief Simple test program to verify VTK grid file reading functionality
 *
 * This program specifically tests the Load_Mesh function with the Ramp_15o_52_18.vtk file
 * to ensure the mesh loading and error handling is working correctly.
 */
int main()
{
    std::cout << "=== VTK Grid File Reader Test ===" << std::endl;
    std::cout << "Testing with: Ramp_15o_52_18.vtk" << std::endl;
    std::cout << "=================================" << std::endl;

    // Set the path to the VTK file
    std::string vtk_file_path = "/Users/rameshkolluru/My_Research/CFD_Solver_withCUDA/Grid_Files/Ramp_Grid_Files/Ramp_15o_52_18.vtk";

    std::cout << "Attempting to load VTK file: " << vtk_file_path << std::endl;

    // Test the Load_Mesh function
    bool load_success = Load_Mesh(vtk_file_path);

    if (load_success)
    {
        std::cout << "\n✅ SUCCESS: VTK file loaded successfully!" << std::endl;
        std::cout << "Grid Statistics:" << std::endl;
        std::cout << "  - Number of physical cells: " << No_Physical_Cells << std::endl;
        std::cout << "  - Number of boundary cells: " << Boundary_Cells.size() << std::endl;
        std::cout << "  - Total cells: " << Cells.size() << std::endl;
        std::cout << "  - 2D Flow: " << (Is_2D_Flow ? "Yes" : "No") << std::endl;

        // Test Form_Cells function
        std::cout << "\nTesting Form_Cells function..." << std::endl;
        bool form_success = Form_Cells(vtk_file_path);

        if (form_success)
        {
            std::cout << "✅ SUCCESS: Form_Cells completed successfully!" << std::endl;
            std::cout << "Final cell count: " << Cells.size() << std::endl;
        }
        else
        {
            std::cout << "❌ ERROR: Form_Cells failed!" << std::endl;
            return 1;
        }

        // Display some cell information if we have cells
        if (!Cells.empty())
        {
            std::cout << "\nSample cell information (first cell):" << std::endl;
            const Cell &first_cell = Cells[0];
            std::cout << "  - Cell ID: " << first_cell.id << std::endl;
            std::cout << "  - Cell Type: " << first_cell.Cell_Type << std::endl;
            std::cout << "  - Number of vertices: " << first_cell.Cell_Vertices.size() / 3 << std::endl;
            std::cout << "  - Cell volume/area: " << first_cell.Vol << std::endl;
        }

        std::cout << "\n🎉 All tests passed successfully!" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "\n❌ ERROR: Failed to load VTK file!" << std::endl;
        std::cout << "Please check:" << std::endl;
        std::cout << "  1. File exists at the specified path" << std::endl;
        std::cout << "  2. File format is valid VTK" << std::endl;
        std::cout << "  3. Read permissions are available" << std::endl;
        return 1;
    }
}