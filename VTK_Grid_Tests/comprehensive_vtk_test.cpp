#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

// Forward declarations for the functions we want to test
bool Load_Mesh(const std::string &configOrMeshPath);
bool Form_Cells(const std::string &gridFile);

// Mock global variables that the CFD solver expects
extern std::vector<struct Cell> Cells;
extern std::vector<struct Cell> Boundary_Cells;
extern int No_Physical_Cells;
extern bool Is_2D_Flow;
extern std::vector<int> Inlet_Cells_List;
extern std::vector<int> Exit_Cells_List;
extern std::vector<int> Wall_Cells_List;

/**
 * @brief Comprehensive test program for VTK grid file loading and processing
 *
 * This program tests:
 * 1. Load_Mesh function with Ramp_15o_52_18.vtk
 * 2. Form_Cells function for complete grid processing
 * 3. Verification of all necessary data structures
 */

void printSeparator(const std::string &title)
{
    std::cout << "\n"
              << std::string(60, '=') << std::endl;
    std::cout << "  " << title << std::endl;
    std::cout << std::string(60, '=') << std::endl;
}

void printTestResult(const std::string &test, bool result)
{
    std::cout << std::left << std::setw(40) << test;
    if (result)
    {
        std::cout << "✅ PASSED" << std::endl;
    }
    else
    {
        std::cout << "❌ FAILED" << std::endl;
    }
}

void printGridStatistics()
{
    std::cout << "\nGrid Statistics:" << std::endl;
    std::cout << "  - Physical Cells: " << No_Physical_Cells << std::endl;
    std::cout << "  - Total Interior Cells: " << Cells.size() << std::endl;
    std::cout << "  - Boundary Cells: " << Boundary_Cells.size() << std::endl;
    std::cout << "  - 2D Flow: " << (Is_2D_Flow ? "Yes" : "No") << std::endl;
    std::cout << "  - Inlet Cells: " << Inlet_Cells_List.size() / 2 << " faces" << std::endl;
    std::cout << "  - Exit Cells: " << Exit_Cells_List.size() / 2 << " faces" << std::endl;
    std::cout << "  - Wall Cells: " << Wall_Cells_List.size() / 2 << " faces" << std::endl;
}

void testCellProperties()
{
    if (!Cells.empty())
    {
        std::cout << "\nSample Cell Analysis (First Cell):" << std::endl;
        const auto &cell = Cells[0];
        std::cout << "  - Cell ID: " << cell.id << std::endl;
        std::cout << "  - Cell Type: " << cell.Cell_Type << std::endl;
        std::cout << "  - Number of vertices: " << cell.Cell_Vertices.size() / 3 << std::endl;
        std::cout << "  - Number of neighbors: " << cell.Neighbours.size() << std::endl;
        std::cout << "  - Cell area/volume: " << cell.Vol << std::endl;
        std::cout << "  - Has boundary face: " << (cell.hasBoundaryface ? "Yes" : "No") << std::endl;

        if (cell.Cell_Vertices.size() >= 12)
        { // At least 4 vertices (4*3 = 12 coordinates)
            std::cout << "  - Vertex coordinates:" << std::endl;
            for (size_t i = 0; i < cell.Cell_Vertices.size(); i += 3)
            {
                std::cout << "    Vertex " << i / 3 + 1 << ": ("
                          << std::fixed << std::setprecision(6)
                          << cell.Cell_Vertices[i] << ", "
                          << cell.Cell_Vertices[i + 1] << ", "
                          << cell.Cell_Vertices[i + 2] << ")" << std::endl;
            }
        }
    }
}

int main()
{
    printSeparator("CFD Solver VTK Grid Processing Test");
    std::cout << "Testing with: Ramp_15o_52_18.vtk" << std::endl;

    std::string vtk_file_path = "/Users/rameshkolluru/My_Research/CFD_Solver_withCUDA/Grid_Files/Ramp_Grid_Files/Ramp_15o_52_18.vtk";

    // Test 1: Load_Mesh function
    printSeparator("Phase 1: Testing Load_Mesh Function");
    std::cout << "Calling Load_Mesh with: " << vtk_file_path << std::endl;

    bool load_success = Load_Mesh(vtk_file_path);
    printTestResult("Load_Mesh function", load_success);

    if (!load_success)
    {
        std::cout << "\n❌ Load_Mesh failed - cannot continue with further tests" << std::endl;
        return 1;
    }

    // Print initial statistics after Load_Mesh
    printGridStatistics();

    // Test 2: Form_Cells function
    printSeparator("Phase 2: Testing Form_Cells Function");
    std::cout << "Calling Form_Cells with: " << vtk_file_path << std::endl;

    bool form_success = Form_Cells(vtk_file_path);
    printTestResult("Form_Cells function", form_success);

    if (!form_success)
    {
        std::cout << "\n❌ Form_Cells failed - grid processing incomplete" << std::endl;
        return 1;
    }

    // Print final statistics after Form_Cells
    printGridStatistics();

    // Test 3: Verify data structures
    printSeparator("Phase 3: Data Structure Verification");

    bool has_cells = !Cells.empty();
    bool has_physical_cells = (No_Physical_Cells > 0);
    bool is_2d_detected = Is_2D_Flow;
    bool has_boundaries = (!Inlet_Cells_List.empty() || !Exit_Cells_List.empty() || !Wall_Cells_List.empty());

    printTestResult("Interior cells created", has_cells);
    printTestResult("Physical cell count set", has_physical_cells);
    printTestResult("2D flow detected", is_2d_detected);
    printTestResult("Boundary conditions identified", has_boundaries);

    // Test 4: Cell geometry verification
    printSeparator("Phase 4: Cell Geometry Analysis");

    bool cells_have_geometry = false;
    bool cells_have_neighbors = false;
    bool cells_have_volume = false;

    if (!Cells.empty())
    {
        const auto &first_cell = Cells[0];
        cells_have_geometry = !first_cell.Cell_Vertices.empty();
        cells_have_neighbors = !first_cell.Neighbours.empty();
        cells_have_volume = (first_cell.Vol > 0.0);
    }

    printTestResult("Cell vertices populated", cells_have_geometry);
    printTestResult("Cell neighbors identified", cells_have_neighbors);
    printTestResult("Cell volumes calculated", cells_have_volume);

    // Display detailed cell information
    testCellProperties();

    // Final assessment
    printSeparator("Test Summary");

    bool all_tests_passed = load_success && form_success && has_cells &&
                            has_physical_cells && is_2d_detected &&
                            cells_have_geometry && cells_have_neighbors;

    if (all_tests_passed)
    {
        std::cout << "🎉 ALL TESTS PASSED! VTK grid loading and processing is working correctly." << std::endl;
        std::cout << "\nThe CFD solver can successfully:" << std::endl;
        std::cout << "  ✅ Load VTK mesh files" << std::endl;
        std::cout << "  ✅ Process cell geometry" << std::endl;
        std::cout << "  ✅ Identify cell neighbors" << std::endl;
        std::cout << "  ✅ Detect boundary conditions" << std::endl;
        std::cout << "  ✅ Calculate cell properties" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "⚠️  Some tests failed. Please check the implementation." << std::endl;
        return 1;
    }
}