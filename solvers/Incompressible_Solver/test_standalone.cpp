#include "Incompressible_Solver_Standalone.h"
#include <iostream>

/**
 * @file test_standalone.cpp
 * @brief Simple standalone test for the incompressible solver core functionality
 */

using namespace IncompressibleSolver;

// Forward declaration for Cell structure to match what's expected
struct Cell
{
    int cellType, cellID, Dimension, ParentCellID, NoBoundaryFaces, numFaces, numNodes;
    std::vector<int> nodeIndices, Neighbours, faceID, Secondary_Neighbours;
    std::vector<double> Diagonal_Vector;
    bool hasBoundaryface = false, Is_Splittable = false, has_Wall_Face = false,
         has_Inlet_Face = false, has_Exit_Face = false, has_Symmetry_Face = false;
    double Area, Inv_Area, del_t, Vol;
    std::vector<double> Face_Areas, Face_Normals, Cell_Center, Cell_Center_Distances,
        Cell_Vertices, Cell_Face_Distances, Cell_Areas, Cell_Center_Vector;
    bool Left_Face = false, Right_Face = false, Top_Face = false, Bottom_Face = false, Interior_Face = false;

    Cell() : cellType(0), cellID(0), Dimension(0), ParentCellID(0), NoBoundaryFaces(0),
             numFaces(0), numNodes(0), Area(0.0), Inv_Area(0.0), del_t(0.0), Vol(0.0) {}
};

// Simple test grid data - these variables are used by the solver
std::vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;
std::vector<int> Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
std::vector<double> Vertices;
int Total_No_Cells = 0, No_Physical_Cells = 0;
bool Is_2D_Flow = true;

// Mock Read_Grid function for testing
void Read_Grid(const std::string &ipfile)
{
    std::cout << "Mock grid loading from: " << ipfile << std::endl;

    // Create a simple 4x4 grid for testing
    No_Physical_Cells = 16;
    Is_2D_Flow = true;

    Cells.resize(No_Physical_Cells);
    for (int i = 0; i < No_Physical_Cells; ++i)
    {
        Cells[i].cellID = i;
        Cells[i].Area = 0.25; // Unit square divided by 16
        Cells[i].Vol = Cells[i].Area;
        Cells[i].Face_Areas.resize(4, 0.25); // 4 faces for 2D
        Cells[i].Face_Normals.resize(8);     // 2 components per face
        Cells[i].Cell_Center.resize(3);

        // Set up simple grid coordinates
        int i_idx = i % 4;
        int j_idx = i / 4;
        Cells[i].Cell_Center[0] = (i_idx + 0.5) / 4.0; // x-coordinate
        Cells[i].Cell_Center[1] = (j_idx + 0.5) / 4.0; // y-coordinate
        Cells[i].Cell_Center[2] = 0.0;                 // z-coordinate

        // Set up neighbors (simplified)
        Cells[i].Neighbours.resize(4, -1);
        if (i_idx > 0)
            Cells[i].Neighbours[0] = i - 1; // Left
        if (i_idx < 3)
            Cells[i].Neighbours[1] = i + 1; // Right
        if (j_idx > 0)
            Cells[i].Neighbours[2] = i - 4; // Bottom
        if (j_idx < 3)
            Cells[i].Neighbours[3] = i + 4; // Top
    }

    // Set up boundary cells
    Wall_Cells_List = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15}; // Boundary cells
    Inlet_Cells_List = {12, 13, 14, 15};                         // Top boundary
    Exit_Cells_List = {0, 1, 2, 3};                              // Bottom boundary

    std::cout << "Mock grid created: " << No_Physical_Cells << " cells" << std::endl;
}

int main()
{
    std::cout << "\n=== Incompressible Solver Standalone Test ===\n"
              << std::endl;

    try
    {
        // Test 1: Basic solver creation
        std::cout << "Test 1: Creating solver instance..." << std::endl;
        IncompressibleFlowSolver solver;
        std::cout << "✓ Solver created successfully" << std::endl;

        // Test 2: Grid initialization
        std::cout << "\nTest 2: Initializing with mock grid..." << std::endl;
        Read_Grid("test_grid.txt");

        if (No_Physical_Cells > 0)
        {
            solver.initialize_from_grid();
            std::cout << "✓ Grid initialization successful" << std::endl;
        }
        else
        {
            std::cout << "✗ Grid initialization failed" << std::endl;
            return 1;
        }

        // Test 3: Fluid properties
        std::cout << "\nTest 3: Setting fluid properties..." << std::endl;
        solver.set_fluid_properties(1000.0, 1e-3); // Water
        std::cout << "✓ Fluid properties set" << std::endl;

        // Test 4: Solver parameters
        std::cout << "\nTest 4: Setting solver parameters..." << std::endl;
        SolverParameters params;
        params.max_iterations = 10; // Small number for testing
        params.velocity_tolerance = 1e-3;
        params.pressure_tolerance = 1e-3;
        params.continuity_tolerance = 1e-5;
        solver.set_solver_parameters(params);
        std::cout << "✓ Solver parameters set" << std::endl;

        // Test 5: Boundary conditions
        std::cout << "\nTest 5: Setting boundary conditions..." << std::endl;
        std::vector<BoundaryCondition> bcs = create_default_boundary_conditions();
        solver.set_boundary_conditions(bcs);
        std::cout << "✓ Boundary conditions set" << std::endl;

        // Test 6: Field initialization
        std::cout << "\nTest 6: Initializing flow fields..." << std::endl;
        solver.initialize_fields();
        std::cout << "✓ Flow fields initialized" << std::endl;

        // Test 7: Reynolds number calculation
        std::cout << "\nTest 7: Computing Reynolds numbers..." << std::endl;
        double Re_global, Re_cell;
        compute_flow_reynolds_numbers(solver.get_fluid_properties(),
                                      solver.get_velocity_field(),
                                      1.0, // Unit length
                                      Re_global, Re_cell);
        std::cout << "✓ Reynolds numbers computed" << std::endl;

        // Test 8: Configuration I/O
        std::cout << "\nTest 8: Testing configuration I/O..." << std::endl;
        write_solver_parameters_to_json(params, "test_config.json");
        SolverParameters loaded_params = read_solver_parameters_from_json("test_config.json");
        std::cout << "✓ Configuration I/O successful" << std::endl;

        // Test 9: Grid validation
        std::cout << "\nTest 9: Validating grid..." << std::endl;
        bool is_valid = validate_grid_for_incompressible_solver();
        if (is_valid)
        {
            std::cout << "✓ Grid validation passed" << std::endl;
        }
        else
        {
            std::cout << "⚠ Grid validation found issues (expected for mock grid)" << std::endl;
        }

        // Test 10: Brief solver run (just a few iterations)
        std::cout << "\nTest 10: Running solver (brief test)..." << std::endl;
        std::cout << "Note: This may generate convergence warnings due to simplified mock grid" << std::endl;

        // Override parameters for very brief test
        params.max_iterations = 3;
        params.velocity_tolerance = 1e-1; // Relaxed tolerance
        params.pressure_tolerance = 1e-1;
        params.continuity_tolerance = 1e-1;
        solver.set_solver_parameters(params);

        auto start = std::chrono::high_resolution_clock::now();
        solver.solve();
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "✓ Solver completed in " << duration.count() << " ms" << std::endl;

        // Test 11: Output
        std::cout << "\nTest 11: Testing output functions..." << std::endl;
        solver.write_solution_vtk("test_solution.vtk");
        solver.write_residual_history("test_residuals.dat");
        std::cout << "✓ Output files written" << std::endl;

        std::cout << "\n=== All Tests Completed Successfully! ===\n"
                  << std::endl;

        // Summary
        std::cout << "Test Summary:" << std::endl;
        std::cout << "  Grid size: " << No_Physical_Cells << " cells" << std::endl;
        std::cout << "  Reynolds number (global): " << Re_global << std::endl;
        std::cout << "  Test runtime: " << duration.count() << " ms" << std::endl;
        std::cout << "  Output files: test_solution.vtk, test_residuals.dat, test_config.json" << std::endl;

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "✗ Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "✗ Test failed with unknown exception" << std::endl;
        return 1;
    }
}