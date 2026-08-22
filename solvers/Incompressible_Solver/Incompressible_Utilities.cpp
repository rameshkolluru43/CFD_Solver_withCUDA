#include "Incompressible_Solver.h"
#include <json/json.h>
#include <fstream>
#include <iostream>
#include <cmath>

namespace IncompressibleSolver
{

    // ================================
    // Utility Functions Implementation
    // ================================

    void initialize_incompressible_solver_with_grid(IncompressibleFlowSolver &solver)
    {
        std::cout << "Initializing incompressible solver with existing grid data..." << std::endl;

        // Ensure grid is loaded
        if (No_Physical_Cells == 0 || Cells.empty())
        {
            std::cerr << "Error: No grid data found. Please load grid first." << std::endl;
            return;
        }

        // Initialize solver with grid
        solver.initialize_from_grid();

        // Set default fluid properties (water at room temperature)
        solver.set_fluid_properties(1000.0, 1e-3); // density = 1000 kg/m³, viscosity = 1e-3 Pa·s

        // Set default solver parameters
        SolverParameters params;
        params.dt = 1e-4;
        params.max_iterations = 1000;
        params.steady_state = true;
        params.velocity_tolerance = 1e-6;
        params.pressure_tolerance = 1e-6;
        params.continuity_tolerance = 1e-8;
        solver.set_solver_parameters(params);

        // Create default boundary conditions
        std::vector<BoundaryCondition> bcs = create_default_boundary_conditions();
        solver.set_boundary_conditions(bcs);

        // Initialize flow fields
        solver.initialize_fields();

        std::cout << "Incompressible solver initialization complete." << std::endl;
    }

    std::vector<BoundaryCondition> create_default_boundary_conditions()
    {
        std::vector<BoundaryCondition> boundary_conditions;

        // Create boundary conditions based on existing boundary classifications

        // Wall boundary condition
        if (!Wall_Cells_List.empty())
        {
            BoundaryCondition wall_bc;
            wall_bc.type = BoundaryCondition::WALL;
            wall_bc.u_value = 0.0;
            wall_bc.v_value = 0.0;
            wall_bc.w_value = 0.0;
            boundary_conditions.push_back(wall_bc);
            std::cout << "Created wall boundary condition for " << Wall_Cells_List.size() << " cells" << std::endl;
        }

        // Inlet boundary condition
        if (!Inlet_Cells_List.empty())
        {
            BoundaryCondition inlet_bc;
            inlet_bc.type = BoundaryCondition::INLET_VELOCITY;
            inlet_bc.u_value = 1.0; // Default inlet velocity of 1 m/s
            inlet_bc.v_value = 0.0;
            inlet_bc.w_value = 0.0;
            boundary_conditions.push_back(inlet_bc);
            std::cout << "Created inlet boundary condition for " << Inlet_Cells_List.size() << " cells" << std::endl;
        }

        // Outlet boundary condition
        if (!Exit_Cells_List.empty())
        {
            BoundaryCondition outlet_bc;
            outlet_bc.type = BoundaryCondition::OUTLET_PRESSURE;
            outlet_bc.p_value = 0.0; // Reference pressure
            boundary_conditions.push_back(outlet_bc);
            std::cout << "Created outlet boundary condition for " << Exit_Cells_List.size() << " cells" << std::endl;
        }

        // Symmetry boundary condition
        if (!Symmetry_Cells_List.empty())
        {
            BoundaryCondition symmetry_bc;
            symmetry_bc.type = BoundaryCondition::SYMMETRY;
            boundary_conditions.push_back(symmetry_bc);
            std::cout << "Created symmetry boundary condition for " << Symmetry_Cells_List.size() << " cells" << std::endl;
        }

        if (boundary_conditions.empty())
        {
            std::cout << "Warning: No boundary conditions created. Using default settings." << std::endl;
        }

        return boundary_conditions;
    }

    SolverParameters read_solver_parameters_from_json(const std::string &filename)
    {
        SolverParameters params;

        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Warning: Cannot open " << filename << ". Using default parameters." << std::endl;
            return params;
        }

        Json::Value root;
        Json::Reader reader;

        if (!reader.parse(file, root))
        {
            std::cerr << "Warning: Failed to parse " << filename << ". Using default parameters." << std::endl;
            return params;
        }

        // Read solver parameters
        if (root.isMember("solver_parameters"))
        {
            Json::Value solver_config = root["solver_parameters"];

            // Time stepping
            if (solver_config.isMember("time_step"))
            {
                params.dt = solver_config["time_step"].asDouble();
            }
            if (solver_config.isMember("max_time"))
            {
                params.max_time = solver_config["max_time"].asDouble();
            }
            if (solver_config.isMember("max_iterations"))
            {
                params.max_iterations = solver_config["max_iterations"].asInt();
            }
            if (solver_config.isMember("steady_state"))
            {
                params.steady_state = solver_config["steady_state"].asBool();
            }

            // Convergence criteria
            if (solver_config.isMember("velocity_tolerance"))
            {
                params.velocity_tolerance = solver_config["velocity_tolerance"].asDouble();
            }
            if (solver_config.isMember("pressure_tolerance"))
            {
                params.pressure_tolerance = solver_config["pressure_tolerance"].asDouble();
            }
            if (solver_config.isMember("continuity_tolerance"))
            {
                params.continuity_tolerance = solver_config["continuity_tolerance"].asDouble();
            }

            // Under-relaxation factors
            if (solver_config.isMember("alpha_u"))
            {
                params.alpha_u = solver_config["alpha_u"].asDouble();
            }
            if (solver_config.isMember("alpha_v"))
            {
                params.alpha_v = solver_config["alpha_v"].asDouble();
            }
            if (solver_config.isMember("alpha_w"))
            {
                params.alpha_w = solver_config["alpha_w"].asDouble();
            }
            if (solver_config.isMember("alpha_p"))
            {
                params.alpha_p = solver_config["alpha_p"].asDouble();
            }

            // Linear solver parameters
            if (solver_config.isMember("max_linear_iterations"))
            {
                params.max_linear_iterations = solver_config["max_linear_iterations"].asInt();
            }
            if (solver_config.isMember("linear_tolerance"))
            {
                params.linear_tolerance = solver_config["linear_tolerance"].asDouble();
            }

            // Output control
            if (solver_config.isMember("output_frequency"))
            {
                params.output_frequency = solver_config["output_frequency"].asInt();
            }
            if (solver_config.isMember("output_format"))
            {
                params.output_format = solver_config["output_format"].asString();
            }
        }

        std::cout << "Solver parameters loaded from " << filename << std::endl;
        return params;
    }

    void write_solver_parameters_to_json(const SolverParameters &params, const std::string &filename)
    {
        Json::Value root;
        Json::Value solver_config;

        // Time stepping
        solver_config["time_step"] = params.dt;
        solver_config["max_time"] = params.max_time;
        solver_config["max_iterations"] = params.max_iterations;
        solver_config["steady_state"] = params.steady_state;

        // Convergence criteria
        solver_config["velocity_tolerance"] = params.velocity_tolerance;
        solver_config["pressure_tolerance"] = params.pressure_tolerance;
        solver_config["continuity_tolerance"] = params.continuity_tolerance;

        // Under-relaxation factors
        solver_config["alpha_u"] = params.alpha_u;
        solver_config["alpha_v"] = params.alpha_v;
        solver_config["alpha_w"] = params.alpha_w;
        solver_config["alpha_p"] = params.alpha_p;

        // Linear solver parameters
        solver_config["max_linear_iterations"] = params.max_linear_iterations;
        solver_config["linear_tolerance"] = params.linear_tolerance;

        // Output control
        solver_config["output_frequency"] = params.output_frequency;
        solver_config["output_format"] = params.output_format;

        root["solver_parameters"] = solver_config;

        // Add example fluid properties
        Json::Value fluid_props;
        fluid_props["density"] = 1000.0; // kg/m³
        fluid_props["viscosity"] = 1e-3; // Pa·s
        fluid_props["description"] = "Water at room temperature";
        root["fluid_properties"] = fluid_props;

        // Add example boundary conditions
        Json::Value boundary_conditions(Json::arrayValue);

        Json::Value wall_bc;
        wall_bc["type"] = "wall";
        wall_bc["u_value"] = 0.0;
        wall_bc["v_value"] = 0.0;
        wall_bc["w_value"] = 0.0;
        wall_bc["description"] = "No-slip wall boundary";
        boundary_conditions.append(wall_bc);

        Json::Value inlet_bc;
        inlet_bc["type"] = "inlet_velocity";
        inlet_bc["u_value"] = 1.0;
        inlet_bc["v_value"] = 0.0;
        inlet_bc["w_value"] = 0.0;
        inlet_bc["description"] = "Uniform velocity inlet";
        boundary_conditions.append(inlet_bc);

        Json::Value outlet_bc;
        outlet_bc["type"] = "outlet_pressure";
        outlet_bc["p_value"] = 0.0;
        outlet_bc["description"] = "Zero gauge pressure outlet";
        boundary_conditions.append(outlet_bc);

        root["boundary_conditions"] = boundary_conditions;

        // Write to file
        std::ofstream file(filename);
        if (file.is_open())
        {
            Json::StreamWriterBuilder builder;
            builder["indentation"] = "  ";
            std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
            writer->write(root, &file);
            file.close();
            std::cout << "Solver parameters written to " << filename << std::endl;
        }
        else
        {
            std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
        }
    }

    bool validate_grid_for_incompressible_solver()
    {
        std::cout << "Validating grid for incompressible solver..." << std::endl;

        bool is_valid = true;
        std::vector<std::string> issues;

        // Check if grid data exists
        if (No_Physical_Cells == 0)
        {
            issues.push_back("No physical cells found");
            is_valid = false;
        }

        if (Cells.empty())
        {
            issues.push_back("Cells vector is empty");
            is_valid = false;
        }

        // Check cell data consistency
        if (No_Physical_Cells != static_cast<int>(Cells.size()))
        {
            issues.push_back("Number of physical cells (" + std::to_string(No_Physical_Cells) +
                             ") does not match Cells vector size (" + std::to_string(Cells.size()) + ")");
            is_valid = false;
        }

        // Check for valid cell geometry
        int cells_with_invalid_geometry = 0;
        int cells_with_missing_neighbors = 0;
        int cells_with_invalid_areas = 0;

        for (int i = 0; i < std::min(No_Physical_Cells, static_cast<int>(Cells.size())); ++i)
        {
            // Check face areas
            if (Cells[i].Face_Areas.empty())
            {
                cells_with_invalid_areas++;
            }
            else
            {
                for (double area : Cells[i].Face_Areas)
                {
                    if (area <= 0.0)
                    {
                        cells_with_invalid_areas++;
                        break;
                    }
                }
            }

            // Check neighbors
            if (Cells[i].Neighbours.empty())
            {
                cells_with_missing_neighbors++;
            }

            // Check cell volume/area
            double cell_size = (Cells[i].Vol > 0.0) ? Cells[i].Vol : Cells[i].Area;
            if (cell_size <= 0.0)
            {
                cells_with_invalid_geometry++;
            }
        }

        if (cells_with_invalid_geometry > 0)
        {
            issues.push_back(std::to_string(cells_with_invalid_geometry) + " cells have invalid geometry");
            is_valid = false;
        }

        if (cells_with_missing_neighbors > 0)
        {
            issues.push_back(std::to_string(cells_with_missing_neighbors) + " cells have no neighbors");
            is_valid = false;
        }

        if (cells_with_invalid_areas > 0)
        {
            issues.push_back(std::to_string(cells_with_invalid_areas) + " cells have invalid face areas");
            is_valid = false;
        }

        // Check boundary conditions setup
        int total_boundary_cells = Wall_Cells_List.size() + Inlet_Cells_List.size() +
                                   Exit_Cells_List.size() + Symmetry_Cells_List.size();

        if (total_boundary_cells == 0)
        {
            issues.push_back("No boundary cells identified");
            // This is a warning, not necessarily invalid
            std::cout << "Warning: No boundary cells identified. You may need to classify boundaries." << std::endl;
        }

        // Check 2D vs 3D consistency
        if (Is_2D_Flow)
        {
            std::cout << "Grid is configured for 2D flow" << std::endl;
        }
        else
        {
            std::cout << "Grid is configured for 3D flow" << std::endl;
        }

        // Print validation results
        if (is_valid)
        {
            std::cout << "Grid validation PASSED" << std::endl;
            std::cout << "Grid statistics:" << std::endl;
            std::cout << "  - Physical cells: " << No_Physical_Cells << std::endl;
            std::cout << "  - Wall cells: " << Wall_Cells_List.size() << std::endl;
            std::cout << "  - Inlet cells: " << Inlet_Cells_List.size() << std::endl;
            std::cout << "  - Exit cells: " << Exit_Cells_List.size() << std::endl;
            std::cout << "  - Symmetry cells: " << Symmetry_Cells_List.size() << std::endl;
            std::cout << "  - 2D Flow: " << (Is_2D_Flow ? "Yes" : "No") << std::endl;
        }
        else
        {
            std::cout << "Grid validation FAILED" << std::endl;
            std::cout << "Issues found:" << std::endl;
            for (const auto &issue : issues)
            {
                std::cout << "  - " << issue << std::endl;
            }
        }

        return is_valid;
    }

    void compute_flow_reynolds_numbers(const FluidProperties &props,
                                       const VelocityField &vel_field,
                                       double characteristic_length,
                                       double &Re_global, double &Re_cell)
    {

        // Compute characteristic velocity
        double u_max = 0.0, u_avg = 0.0;
        int valid_cells = 0;

        for (size_t i = 0; i < vel_field.u.size(); ++i)
        {
            if (i < vel_field.v.size())
            {
                double u_mag = std::sqrt(vel_field.u[i] * vel_field.u[i] + vel_field.v[i] * vel_field.v[i]);
                if (!vel_field.w.empty() && i < vel_field.w.size())
                {
                    u_mag = std::sqrt(u_mag * u_mag + vel_field.w[i] * vel_field.w[i]);
                }

                u_max = std::max(u_max, u_mag);
                u_avg += u_mag;
                valid_cells++;
            }
        }

        if (valid_cells > 0)
        {
            u_avg /= valid_cells;
        }

        // Global Reynolds number based on characteristic length and average velocity
        Re_global = (props.kinematic_viscosity > 0.0) ? (u_avg * characteristic_length) / props.kinematic_viscosity : 0.0;

        // Cell Reynolds number based on typical cell size
        double typical_cell_size = 0.0;
        int cells_counted = 0;

        for (int i = 0; i < std::min(No_Physical_Cells, static_cast<int>(Cells.size())); ++i)
        {
            double cell_size = (Cells[i].Vol > 0.0) ? std::pow(Cells[i].Vol, 1.0 / 3.0) : std::sqrt(Cells[i].Area);
            if (cell_size > 0.0)
            {
                typical_cell_size += cell_size;
                cells_counted++;
            }
        }

        if (cells_counted > 0)
        {
            typical_cell_size /= cells_counted;
        }

        Re_cell = (props.kinematic_viscosity > 0.0 && typical_cell_size > 0.0) ? (u_max * typical_cell_size) / props.kinematic_viscosity : 0.0;

        std::cout << "Reynolds number analysis:" << std::endl;
        std::cout << "  - Characteristic length: " << characteristic_length << " m" << std::endl;
        std::cout << "  - Typical cell size: " << typical_cell_size << " m" << std::endl;
        std::cout << "  - Average velocity: " << u_avg << " m/s" << std::endl;
        std::cout << "  - Maximum velocity: " << u_max << " m/s" << std::endl;
        std::cout << "  - Global Reynolds number: " << Re_global << std::endl;
        std::cout << "  - Cell Reynolds number: " << Re_cell << std::endl;

        // Provide guidance based on Reynolds numbers
        if (Re_global < 1.0)
        {
            std::cout << "  - Flow regime: Highly viscous (Stokes flow)" << std::endl;
        }
        else if (Re_global < 100.0)
        {
            std::cout << "  - Flow regime: Laminar, viscous-dominated" << std::endl;
        }
        else if (Re_global < 2000.0)
        {
            std::cout << "  - Flow regime: Laminar" << std::endl;
        }
        else if (Re_global < 4000.0)
        {
            std::cout << "  - Flow regime: Transitional" << std::endl;
        }
        else
        {
            std::cout << "  - Flow regime: Turbulent (consider turbulence modeling)" << std::endl;
        }

        if (Re_cell > 2.0)
        {
            std::cout << "  - Warning: Cell Reynolds number > 2. Consider refining the grid or using upwind schemes." << std::endl;
        }
    }

} // namespace IncompressibleSolver