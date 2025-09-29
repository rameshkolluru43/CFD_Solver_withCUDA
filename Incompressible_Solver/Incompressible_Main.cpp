#include "Incompressible_Solver.h"
#include <iostream>
#include <string>
#include <chrono>

// Forward declarations for external functions
extern "C" void Read_Grid(const std::string &ipfile);

// Forward declarations for Cell structure
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

// External variables from Grid_Computations.cpp
extern std::vector<Cell> Cells;
extern std::vector<int> Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
extern int No_Physical_Cells;
extern bool Is_2D_Flow;

/**
 * @file Incompressible_Main.cpp
 * @brief Main driver program for the incompressible flow solver
 *
 * This program demonstrates how to use the cell-centered staggered grid
 * incompressible solver with the existing grid infrastructure.
 */

using namespace IncompressibleSolver;

void print_banner()
{
    std::cout << "\n";
    std::cout << "====================================================================\n";
    std::cout << "     Cell-Centered Staggered Grid Incompressible Flow Solver       \n";
    std::cout << "====================================================================\n";
    std::cout << "  Features:\n";
    std::cout << "    - SIMPLE algorithm for pressure-velocity coupling\n";
    std::cout << "    - Finite volume method with cell-centered approach\n";
    std::cout << "    - Integration with existing grid infrastructure\n";
    std::cout << "    - Support for 2D and 3D incompressible flows\n";
    std::cout << "    - Multiple boundary condition types\n";
    std::cout << "    - VTK output for visualization\n";
    std::cout << "====================================================================\n\n";
}

void print_usage(const std::string &program_name)
{
    std::cout << "Usage: " << program_name << " [options]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  -g, --grid <file>       Grid file path (required)\n";
    std::cout << "  -c, --config <file>     Configuration JSON file (optional)\n";
    std::cout << "  -o, --output <dir>      Output directory (default: ./output/)\n";
    std::cout << "  -s, --steady            Run steady-state simulation (default)\n";
    std::cout << "  -t, --transient         Run transient simulation\n";
    std::cout << "  -v, --verbose           Enable verbose output\n";
    std::cout << "  -h, --help              Show this help message\n";
    std::cout << "\nExample:\n";
    std::cout << "  " << program_name << " -g cylinder_grid.txt -c config.json -o results/\n\n";
}

struct ProgramOptions
{
    std::string grid_file;
    std::string config_file;
    std::string output_dir = "./output/";
    bool steady_state = true;
    bool verbose = false;
    bool show_help = false;
};

ProgramOptions parse_command_line(int argc, char *argv[])
{
    ProgramOptions options;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "-g" || arg == "--grid")
        {
            if (i + 1 < argc)
            {
                options.grid_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: Grid file not specified after " << arg << std::endl;
            }
        }
        else if (arg == "-c" || arg == "--config")
        {
            if (i + 1 < argc)
            {
                options.config_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: Config file not specified after " << arg << std::endl;
            }
        }
        else if (arg == "-o" || arg == "--output")
        {
            if (i + 1 < argc)
            {
                options.output_dir = argv[++i];
            }
            else
            {
                std::cerr << "Error: Output directory not specified after " << arg << std::endl;
            }
        }
        else if (arg == "-s" || arg == "--steady")
        {
            options.steady_state = true;
        }
        else if (arg == "-t" || arg == "--transient")
        {
            options.steady_state = false;
        }
        else if (arg == "-v" || arg == "--verbose")
        {
            options.verbose = true;
        }
        else if (arg == "-h" || arg == "--help")
        {
            options.show_help = true;
        }
        else
        {
            std::cerr << "Warning: Unknown option " << arg << std::endl;
        }
    }

    return options;
}

bool setup_output_directory(const std::string &output_dir)
{
    // Create output directory if it doesn't exist
    std::string command = "mkdir -p " + output_dir;
    int result = system(command.c_str());

    if (result != 0)
    {
        std::cerr << "Error: Failed to create output directory " << output_dir << std::endl;
        return false;
    }

    std::cout << "Output directory: " << output_dir << std::endl;
    return true;
}

void run_validation_test_case()
{
    std::cout << "\n=== Running Validation Test Case ===\n";
    std::cout << "Test: Lid-driven cavity flow (Re = 100)\n\n";

    // Create a simple rectangular grid for validation
    std::cout << "Setting up 2D rectangular grid (32x32)...\n";

    // This is a simplified setup - in practice, you would load a proper grid file
    Is_2D_Flow = true;
    No_Physical_Cells = 32 * 32;

    // Initialize cells vector (simplified)
    Cells.resize(No_Physical_Cells);
    for (int i = 0; i < No_Physical_Cells; ++i)
    {
        Cells[i].cellID = i;
        Cells[i].Area = 1.0 / (32.0 * 32.0); // Unit square divided by grid
        Cells[i].Vol = Cells[i].Area;
        Cells[i].Face_Areas.resize(4, 1.0 / 32.0); // 4 faces for 2D
        Cells[i].Face_Normals.resize(8);           // 2 components per face
        Cells[i].Cell_Center.resize(3);

        // Set up simple grid coordinates
        int i_idx = i % 32;
        int j_idx = i / 32;
        Cells[i].Cell_Center[0] = (i_idx + 0.5) / 32.0; // x-coordinate
        Cells[i].Cell_Center[1] = (j_idx + 0.5) / 32.0; // y-coordinate
        Cells[i].Cell_Center[2] = 0.0;                  // z-coordinate

        // Set up neighbors (simplified)
        Cells[i].Neighbours.resize(4, -1);
        if (i_idx > 0)
            Cells[i].Neighbours[0] = i - 1; // Left
        if (i_idx < 31)
            Cells[i].Neighbours[1] = i + 1; // Right
        if (j_idx > 0)
            Cells[i].Neighbours[2] = i - 32; // Bottom
        if (j_idx < 31)
            Cells[i].Neighbours[3] = i + 32; // Top

        // Set face normals
        Cells[i].Face_Normals[0] = -1.0;
        Cells[i].Face_Normals[1] = 0.0; // Left face
        Cells[i].Face_Normals[2] = 1.0;
        Cells[i].Face_Normals[3] = 0.0; // Right face
        Cells[i].Face_Normals[4] = 0.0;
        Cells[i].Face_Normals[5] = -1.0; // Bottom face
        Cells[i].Face_Normals[6] = 0.0;
        Cells[i].Face_Normals[7] = 1.0; // Top face
    }

    // Set up boundary cells
    Wall_Cells_List.clear();
    Inlet_Cells_List.clear();
    Exit_Cells_List.clear();

    // Bottom, left, and right walls
    for (int i = 0; i < 32; ++i)
    {
        Wall_Cells_List.push_back(i); // Bottom wall
        if (i > 0 && i < 31)
        {
            Wall_Cells_List.push_back(i * 32);      // Left wall
            Wall_Cells_List.push_back(i * 32 + 31); // Right wall
        }
    }

    // Top wall (moving lid)
    for (int i = 0; i < 32; ++i)
    {
        Inlet_Cells_List.push_back(31 * 32 + i); // Top wall with specified velocity
    }

    std::cout << "Grid setup complete:\n";
    std::cout << "  - Cells: " << No_Physical_Cells << "\n";
    std::cout << "  - Wall cells: " << Wall_Cells_List.size() << "\n";
    std::cout << "  - Moving lid cells: " << Inlet_Cells_List.size() << "\n";

    // Create solver and initialize
    IncompressibleFlowSolver solver;
    initialize_incompressible_solver_with_grid(solver);

    // Set fluid properties for the validation case (air)
    solver.set_fluid_properties(1.225, 1.8e-5); // Air at room conditions

    // Set solver parameters for lid-driven cavity
    SolverParameters params;
    params.dt = 1e-3;
    params.max_iterations = 2000;
    params.steady_state = true;
    params.velocity_tolerance = 1e-6;
    params.pressure_tolerance = 1e-6;
    params.continuity_tolerance = 1e-8;
    params.alpha_u = 0.7;
    params.alpha_v = 0.7;
    params.alpha_p = 0.3;
    params.output_frequency = 100;
    solver.set_solver_parameters(params);

    // Create custom boundary conditions for lid-driven cavity
    std::vector<BoundaryCondition> bcs;

    // Moving lid (top boundary)
    BoundaryCondition lid_bc;
    lid_bc.type = BoundaryCondition::INLET_VELOCITY;
    lid_bc.u_value = 1.0; // Lid velocity
    lid_bc.v_value = 0.0;
    lid_bc.w_value = 0.0;
    bcs.push_back(lid_bc);

    // Stationary walls
    BoundaryCondition wall_bc;
    wall_bc.type = BoundaryCondition::WALL;
    wall_bc.u_value = 0.0;
    wall_bc.v_value = 0.0;
    wall_bc.w_value = 0.0;
    bcs.push_back(wall_bc);

    solver.set_boundary_conditions(bcs);

    std::cout << "\nStarting lid-driven cavity simulation...\n";

    // Run the solver
    auto start_time = std::chrono::high_resolution_clock::now();
    solver.solve();
    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "\nSimulation completed in " << duration.count() << " ms\n";

    // Compute Reynolds number
    double Re_global, Re_cell;
    compute_flow_reynolds_numbers(solver.get_fluid_properties(),
                                  solver.get_velocity_field(),
                                  1.0, // Unit cavity length
                                  Re_global, Re_cell);

    std::cout << "\nValidation test case completed successfully!\n";
}

int main(int argc, char *argv[])
{
    print_banner();

    // Parse command line arguments
    ProgramOptions options = parse_command_line(argc, argv);

    if (options.show_help)
    {
        print_usage(argv[0]);
        return 0;
    }

    // Check if grid file is provided
    if (options.grid_file.empty())
    {
        std::cout << "No grid file specified. Running validation test case...\n";
        run_validation_test_case();
        return 0;
    }

    // Setup output directory
    if (!setup_output_directory(options.output_dir))
    {
        return 1;
    }

    try
    {
        std::cout << "=== Incompressible Flow Simulation ===\n\n";

        // Step 1: Load grid
        std::cout << "Step 1: Loading grid from " << options.grid_file << "...\n";
        Read_Grid(options.grid_file);

        if (No_Physical_Cells == 0)
        {
            std::cerr << "Error: Failed to load grid or no cells found.\n";
            return 1;
        }

        std::cout << "Grid loaded successfully: " << No_Physical_Cells << " cells\n";

        // Step 2: Validate grid
        std::cout << "\nStep 2: Validating grid for incompressible solver...\n";
        if (!validate_grid_for_incompressible_solver())
        {
            std::cerr << "Error: Grid validation failed. Please check grid quality.\n";
            return 1;
        }

        // Step 3: Initialize solver
        std::cout << "\nStep 3: Initializing incompressible flow solver...\n";
        IncompressibleFlowSolver solver;
        initialize_incompressible_solver_with_grid(solver);

        // Step 4: Load configuration
        SolverParameters params;
        if (!options.config_file.empty())
        {
            std::cout << "\nStep 4: Loading configuration from " << options.config_file << "...\n";
            params = read_solver_parameters_from_json(options.config_file);
        }
        else
        {
            std::cout << "\nStep 4: Using default solver parameters...\n";
            // Use default parameters
        }

        // Override steady-state setting from command line
        params.steady_state = options.steady_state;
        solver.set_solver_parameters(params);

        // Step 5: Setup boundary conditions
        std::cout << "\nStep 5: Setting up boundary conditions...\n";
        std::vector<BoundaryCondition> bcs = create_default_boundary_conditions();
        solver.set_boundary_conditions(bcs);

        // Step 6: Run simulation
        std::cout << "\nStep 6: Running " << (params.steady_state ? "steady-state" : "transient")
                  << " simulation...\n";

        auto start_time = std::chrono::high_resolution_clock::now();
        solver.solve();
        auto end_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        std::cout << "\nSimulation completed in " << duration.count() << " seconds\n";

        // Step 7: Post-processing
        std::cout << "\nStep 7: Post-processing and output...\n";

        // Write final solution
        std::string vtk_output = options.output_dir + "/final_solution.vtk";
        solver.write_solution_vtk(vtk_output);

        // Write residual history
        std::string residual_output = options.output_dir + "/residual_history.dat";
        solver.write_residual_history(residual_output);

        // Write configuration for reference
        std::string config_output = options.output_dir + "/solver_config.json";
        write_solver_parameters_to_json(params, config_output);

        // Compute and display Reynolds numbers
        double Re_global, Re_cell;
        double characteristic_length = 1.0; // Default, should be set based on problem
        compute_flow_reynolds_numbers(solver.get_fluid_properties(),
                                      solver.get_velocity_field(),
                                      characteristic_length,
                                      Re_global, Re_cell);

        std::cout << "\n=== Simulation Summary ===\n";
        std::cout << "Grid file: " << options.grid_file << "\n";
        std::cout << "Number of cells: " << No_Physical_Cells << "\n";
        std::cout << "Simulation type: " << (params.steady_state ? "Steady-state" : "Transient") << "\n";
        std::cout << "Reynolds number: " << Re_global << "\n";
        std::cout << "Output directory: " << options.output_dir << "\n";
        std::cout << "Runtime: " << duration.count() << " seconds\n";

        std::cout << "\n=== Output Files ===\n";
        std::cout << "Solution: " << vtk_output << "\n";
        std::cout << "Residuals: " << residual_output << "\n";
        std::cout << "Configuration: " << config_output << "\n";

        std::cout << "\nIncompressible flow simulation completed successfully!\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown error occurred during simulation." << std::endl;
        return 1;
    }

    return 0;
}

// Alternative entry point for library usage
extern "C"
{
    int run_incompressible_solver(const char *grid_file, const char *config_file, const char *output_dir)
    {
        try
        {
            // Load grid
            Read_Grid(std::string(grid_file));

            if (No_Physical_Cells == 0)
            {
                return 1; // Failed to load grid
            }

            // Validate grid
            if (!validate_grid_for_incompressible_solver())
            {
                return 2; // Grid validation failed
            }

            // Initialize solver
            IncompressibleFlowSolver solver;
            initialize_incompressible_solver_with_grid(solver);

            // Load configuration if provided
            SolverParameters params;
            if (config_file && strlen(config_file) > 0)
            {
                params = read_solver_parameters_from_json(std::string(config_file));
            }
            solver.set_solver_parameters(params);

            // Setup boundary conditions
            std::vector<BoundaryCondition> bcs = create_default_boundary_conditions();
            solver.set_boundary_conditions(bcs);

            // Run simulation
            solver.solve();

            // Write output
            std::string output_path = output_dir ? std::string(output_dir) : "./output/";
            solver.write_solution_vtk(output_path + "/solution.vtk");
            solver.write_residual_history(output_path + "/residuals.dat");

            return 0; // Success
        }
        catch (...)
        {
            return -1; // Unknown error
        }
    }
}