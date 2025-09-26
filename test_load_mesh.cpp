#include <iostream>
#include <string>

/**
 * @brief Test program to verify the CFD solver's Load_Mesh function
 *
 * This test directly calls the Load_Mesh function with the Ramp VTK file
 * to ensure the full mesh loading pipeline works correctly.
 */

// We'll test the Load_Mesh function directly
bool Load_Mesh(const std::string &configOrMeshPath);

int main()
{
    std::cout << "=== CFD Solver Load_Mesh Function Test ===" << std::endl;
    std::cout << "Testing with: Ramp_15o_52_18.vtk" << std::endl;
    std::cout << "==========================================" << std::endl;

    std::string vtk_file_path = "/Users/rameshkolluru/My_Research/CFD_Solver_withCUDA/Grid_Files/Ramp_Grid_Files/Ramp_15o_52_18.vtk";

    std::cout << "Calling Load_Mesh function with: " << vtk_file_path << std::endl;

    bool success = Load_Mesh(vtk_file_path);

    if (success)
    {
        std::cout << "\n✅ SUCCESS: Load_Mesh function completed successfully!" << std::endl;
        std::cout << "The boolean return value indicates successful mesh loading." << std::endl;
        std::cout << "\n🎉 CFD Solver Load_Mesh test PASSED!" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "\n❌ ERROR: Load_Mesh function failed!" << std::endl;
        std::cout << "The boolean return value indicates mesh loading failure." << std::endl;
        return 1;
    }
}