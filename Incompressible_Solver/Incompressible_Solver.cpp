#include "Incompressible_Solver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

// Forward declarations for external variables from Grid_Computations.cpp
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
extern std::vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;
extern std::vector<int> Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
extern std::vector<double> Vertices;
extern int Total_No_Cells, No_Physical_Cells;
extern bool Is_2D_Flow;

namespace IncompressibleSolver
{

    // ================================
    // Constructor and Destructor
    // ================================

    IncompressibleFlowSolver::IncompressibleFlowSolver()
        : num_cells(0), num_faces(0), num_boundary_faces(0), is_2d_flow(true)
    {
        std::cout << "Incompressible Flow Solver initialized" << std::endl;
    }

    IncompressibleFlowSolver::~IncompressibleFlowSolver()
    {
        std::cout << "Incompressible Flow Solver destroyed" << std::endl;
    }

    // ================================
    // Initialization Methods
    // ================================

    void IncompressibleFlowSolver::initialize_from_grid()
    {
        std::cout << "Initializing incompressible solver from existing grid..." << std::endl;

        // Use existing grid data from globals
        num_cells = No_Physical_Cells;
        is_2d_flow = Is_2D_Flow;

        // Count total faces (estimate based on cell connectivity)
        num_faces = 0;
        num_boundary_faces = 0;

        for (int i = 0; i < num_cells; ++i)
        {
            if (i < Cells.size())
            {
                num_faces += Cells[i].Face_Areas.size();
                if (Cells[i].hasBoundaryface)
                {
                    num_boundary_faces += Cells[i].NoBoundaryFaces;
                }
            }
        }

        // Account for shared internal faces (each internal face counted twice)
        num_faces = (num_faces + num_boundary_faces) / 2;

        std::cout << "Grid information:" << std::endl;
        std::cout << "  Number of cells: " << num_cells << std::endl;
        std::cout << "  Number of faces: " << num_faces << std::endl;
        std::cout << "  Number of boundary faces: " << num_boundary_faces << std::endl;
        std::cout << "  2D Flow: " << (is_2d_flow ? "Yes" : "No") << std::endl;

        // Resize flow fields
        velocity_field.resize(num_cells, num_faces);
        pressure_field.resize(num_cells);

        // Initialize matrices
        momentum_matrix_u.resize(num_cells, std::vector<double>(num_cells, 0.0));
        momentum_matrix_v.resize(num_cells, std::vector<double>(num_cells, 0.0));
        if (!is_2d_flow)
        {
            momentum_matrix_w.resize(num_cells, std::vector<double>(num_cells, 0.0));
        }
        pressure_matrix.resize(num_cells, std::vector<double>(num_cells, 0.0));

        momentum_rhs_u.resize(num_cells, 0.0);
        momentum_rhs_v.resize(num_cells, 0.0);
        if (!is_2d_flow)
        {
            momentum_rhs_w.resize(num_cells, 0.0);
        }
        pressure_rhs.resize(num_cells, 0.0);

        std::cout << "Incompressible solver initialization complete." << std::endl;
    }

    void IncompressibleFlowSolver::set_fluid_properties(double density, double viscosity)
    {
        fluid_properties.density = density;
        fluid_properties.viscosity = viscosity;
        fluid_properties.update_kinematic_viscosity();

        std::cout << "Fluid properties set:" << std::endl;
        std::cout << "  Density: " << density << " kg/m³" << std::endl;
        std::cout << "  Dynamic viscosity: " << viscosity << " Pa·s" << std::endl;
        std::cout << "  Kinematic viscosity: " << fluid_properties.kinematic_viscosity << " m²/s" << std::endl;
    }

    void IncompressibleFlowSolver::set_solver_parameters(const SolverParameters &params)
    {
        solver_params = params;

        std::cout << "Solver parameters set:" << std::endl;
        std::cout << "  Time step: " << solver_params.dt << " s" << std::endl;
        std::cout << "  Maximum iterations: " << solver_params.max_iterations << std::endl;
        std::cout << "  Velocity tolerance: " << solver_params.velocity_tolerance << std::endl;
        std::cout << "  Pressure tolerance: " << solver_params.pressure_tolerance << std::endl;
        std::cout << "  Continuity tolerance: " << solver_params.continuity_tolerance << std::endl;
    }

    void IncompressibleFlowSolver::set_boundary_conditions(const std::vector<BoundaryCondition> &bcs)
    {
        boundary_conditions = bcs;
        std::cout << "Boundary conditions set: " << boundary_conditions.size() << " conditions" << std::endl;
    }

    void IncompressibleFlowSolver::initialize_fields()
    {
        std::cout << "Initializing flow fields..." << std::endl;

        // Initialize velocity field to zero (or specified initial conditions)
        std::fill(velocity_field.u.begin(), velocity_field.u.end(), 0.0);
        std::fill(velocity_field.v.begin(), velocity_field.v.end(), 0.0);
        if (!is_2d_flow)
        {
            std::fill(velocity_field.w.begin(), velocity_field.w.end(), 0.0);
        }

        // Initialize pressure field (reference pressure = 0)
        std::fill(pressure_field.p.begin(), pressure_field.p.end(), 0.0);

        // Initialize intermediate fields
        velocity_field.u_star = velocity_field.u;
        velocity_field.v_star = velocity_field.v;
        if (!is_2d_flow)
        {
            velocity_field.w_star = velocity_field.w;
        }

        std::cout << "Flow fields initialized." << std::endl;
    }

    // ================================
    // Main Solver Methods
    // ================================

    void IncompressibleFlowSolver::solve()
    {
        std::cout << "Starting incompressible flow solution..." << std::endl;

        if (solver_params.steady_state)
        {
            solve_steady_state();
        }
        else
        {
            solve_transient();
        }

        std::cout << "Incompressible flow solution complete." << std::endl;
    }

    void IncompressibleFlowSolver::solve_steady_state()
    {
        std::cout << "Solving steady-state incompressible flow using SIMPLE algorithm..." << std::endl;

        bool converged = false;
        int iteration = 0;

        // Clear residual history
        velocity_residuals.clear();
        pressure_residuals.clear();
        continuity_residuals.clear();

        while (!converged && iteration < solver_params.max_iterations)
        {
            iteration++;

            // Store old pressure field
            pressure_field.p_old = pressure_field.p;

            // SIMPLE algorithm steps
            std::cout << "\n--- Iteration " << iteration << " ---" << std::endl;

            // Step 1: Momentum predictor
            std::cout << "  Step 1: Momentum predictor..." << std::endl;
            momentum_predictor_step();

            // Step 2: Pressure correction
            std::cout << "  Step 2: Pressure correction..." << std::endl;
            pressure_correction_step();

            // Step 3: Velocity correction
            std::cout << "  Step 3: Velocity correction..." << std::endl;
            velocity_corrector_step();

            // Step 4: Update face velocities
            std::cout << "  Step 4: Update face velocities..." << std::endl;
            update_face_velocities();

            // Check convergence
            converged = check_convergence();

            // Print iteration info
            print_iteration_info(iteration);

            // Write output at specified frequency
            if (iteration % solver_params.output_frequency == 0)
            {
                std::string filename = "solution_iter_" + std::to_string(iteration) + ".vtk";
                write_solution_vtk(filename);
            }
        }

        if (converged)
        {
            std::cout << "\nSolution converged in " << iteration << " iterations!" << std::endl;
        }
        else
        {
            std::cout << "\nSolution did not converge within " << solver_params.max_iterations
                      << " iterations." << std::endl;
        }

        // Write final solution
        write_solution_vtk("final_solution.vtk");
        write_residual_history("residual_history.dat");
    }

    void IncompressibleFlowSolver::solve_transient()
    {
        std::cout << "Solving transient incompressible flow..." << std::endl;

        double current_time = 0.0;
        int time_step = 0;

        while (current_time < solver_params.max_time)
        {
            time_step++;
            current_time += solver_params.dt;

            std::cout << "\nTime step " << time_step << ", Time = " << current_time << " s" << std::endl;

            // Store old fields for time derivative
            std::vector<double> u_old = velocity_field.u;
            std::vector<double> v_old = velocity_field.v;
            std::vector<double> w_old = velocity_field.w;

            // SIMPLE iterations for this time step
            bool converged = false;
            int iteration = 0;

            while (!converged && iteration < solver_params.max_iterations)
            {
                iteration++;

                // SIMPLE algorithm with time derivative
                momentum_predictor_step();
                pressure_correction_step();
                velocity_corrector_step();
                update_face_velocities();

                converged = check_convergence();

                if (iteration % 10 == 0)
                {
                    print_iteration_info(iteration);
                }
            }

            // Write output at specified frequency
            if (time_step % solver_params.output_frequency == 0)
            {
                std::string filename = "solution_t_" + std::to_string(time_step) + ".vtk";
                write_solution_vtk(filename);
            }
        }

        std::cout << "\nTransient solution complete!" << std::endl;
    }

    // ================================
    // SIMPLE Algorithm Components
    // ================================

    void IncompressibleFlowSolver::momentum_predictor_step()
    {
        // Assemble momentum equations
        assemble_momentum_equations();

        // Apply velocity boundary conditions
        apply_velocity_boundary_conditions();

        // Solve momentum equations
        solve_linear_system_velocity('u');
        solve_linear_system_velocity('v');
        if (!is_2d_flow)
        {
            solve_linear_system_velocity('w');
        }

        // Apply under-relaxation
        for (int i = 0; i < num_cells; ++i)
        {
            velocity_field.u_star[i] = solver_params.alpha_u * velocity_field.u_star[i] +
                                       (1.0 - solver_params.alpha_u) * velocity_field.u[i];
            velocity_field.v_star[i] = solver_params.alpha_v * velocity_field.v_star[i] +
                                       (1.0 - solver_params.alpha_v) * velocity_field.v[i];
            if (!is_2d_flow)
            {
                velocity_field.w_star[i] = solver_params.alpha_w * velocity_field.w_star[i] +
                                           (1.0 - solver_params.alpha_w) * velocity_field.w[i];
            }
        }
    }

    void IncompressibleFlowSolver::pressure_correction_step()
    {
        // Assemble pressure correction equation
        assemble_pressure_equation();

        // Apply pressure boundary conditions
        apply_pressure_boundary_conditions();

        // Solve pressure correction equation
        solve_linear_system_pressure();

        // Apply under-relaxation to pressure
        for (int i = 0; i < num_cells; ++i)
        {
            pressure_field.p[i] += solver_params.alpha_p * pressure_field.p_prime[i];
        }
    }

    void IncompressibleFlowSolver::velocity_corrector_step()
    {
        // Correct velocities based on pressure correction
        for (int i = 0; i < num_cells; ++i)
        {
            if (i < Cells.size() && i < pressure_field.p_prime.size())
            {
                // Get cell volume
                double cell_volume = Cells[i].Vol;
                if (cell_volume <= 0.0)
                {
                    cell_volume = Cells[i].Area; // For 2D, use area
                }

                // Compute pressure correction gradients and correct velocities
                double dp_dx = 0.0, dp_dy = 0.0, dp_dz = 0.0;

                // Simple gradient calculation using neighboring cells
                int num_neighbors = 0;
                for (size_t face = 0; face < Cells[i].Face_Areas.size(); ++face)
                {
                    if (face < Cells[i].Neighbours.size())
                    {
                        int neighbor = Cells[i].Neighbours[face];
                        if (neighbor >= 0 && neighbor < pressure_field.p_prime.size())
                        {
                            double face_area = Cells[i].Face_Areas[face];
                            double nx = (face * 2 < Cells[i].Face_Normals.size()) ? Cells[i].Face_Normals[face * 2] : 0.0;
                            double ny = (face * 2 + 1 < Cells[i].Face_Normals.size()) ? Cells[i].Face_Normals[face * 2 + 1] : 0.0;

                            double dp = pressure_field.p_prime[neighbor] - pressure_field.p_prime[i];
                            dp_dx += dp * nx * face_area / cell_volume;
                            dp_dy += dp * ny * face_area / cell_volume;
                            num_neighbors++;
                        }
                    }
                }

                if (num_neighbors > 0)
                {
                    // Velocity correction (simplified, should use proper coefficients)
                    double dt_over_rho = solver_params.dt / fluid_properties.density;
                    velocity_field.u[i] = velocity_field.u_star[i] - dt_over_rho * dp_dx;
                    velocity_field.v[i] = velocity_field.v_star[i] - dt_over_rho * dp_dy;
                    if (!is_2d_flow)
                    {
                        velocity_field.w[i] = velocity_field.w_star[i] - dt_over_rho * dp_dz;
                    }
                }
            }
        }
    }

    void IncompressibleFlowSolver::update_face_velocities()
    {
        // Interpolate cell-centered velocities to faces
        interpolate_cell_to_face();

        // Apply mass conservation constraint at faces if needed
        // This is a simplified implementation
    }

    // ================================
    // Discretization Methods
    // ================================

    void IncompressibleFlowSolver::assemble_momentum_equations()
    {
        // Clear matrices and RHS
        for (auto &row : momentum_matrix_u)
        {
            std::fill(row.begin(), row.end(), 0.0);
        }
        for (auto &row : momentum_matrix_v)
        {
            std::fill(row.begin(), row.end(), 0.0);
        }
        if (!is_2d_flow)
        {
            for (auto &row : momentum_matrix_w)
            {
                std::fill(row.begin(), row.end(), 0.0);
            }
        }

        std::fill(momentum_rhs_u.begin(), momentum_rhs_u.end(), 0.0);
        std::fill(momentum_rhs_v.begin(), momentum_rhs_v.end(), 0.0);
        if (!is_2d_flow)
        {
            std::fill(momentum_rhs_w.begin(), momentum_rhs_w.end(), 0.0);
        }

        // Assemble for each cell
        for (int i = 0; i < num_cells; ++i)
        {
            discretize_convection_term(i, 'u');
            discretize_convection_term(i, 'v');
            if (!is_2d_flow)
            {
                discretize_convection_term(i, 'w');
            }

            discretize_diffusion_term(i, 'u');
            discretize_diffusion_term(i, 'v');
            if (!is_2d_flow)
            {
                discretize_diffusion_term(i, 'w');
            }

            discretize_pressure_gradient(i, 'u');
            discretize_pressure_gradient(i, 'v');
            if (!is_2d_flow)
            {
                discretize_pressure_gradient(i, 'w');
            }
        }
    }

    void IncompressibleFlowSolver::assemble_pressure_equation()
    {
        // Clear pressure matrix and RHS
        for (auto &row : pressure_matrix)
        {
            std::fill(row.begin(), row.end(), 0.0);
        }
        std::fill(pressure_rhs.begin(), pressure_rhs.end(), 0.0);

        // Assemble pressure correction equation based on mass conservation
        for (int i = 0; i < num_cells; ++i)
        {
            discretize_mass_conservation(i);
        }
    }

    void IncompressibleFlowSolver::discretize_convection_term(int cell_id, char component)
    {
        if (cell_id >= Cells.size())
            return;

        // Simplified convection term discretization using upwind scheme
        double cell_volume = (Cells[cell_id].Vol > 0.0) ? Cells[cell_id].Vol : Cells[cell_id].Area;

        for (size_t face = 0; face < Cells[cell_id].Face_Areas.size(); ++face)
        {
            if (face < Cells[cell_id].Neighbours.size())
            {
                int neighbor = Cells[cell_id].Neighbours[face];
                double face_area = Cells[cell_id].Face_Areas[face];

                if (neighbor >= 0 && neighbor < num_cells)
                {
                    // Mass flux through face (simplified)
                    double mass_flux = fluid_properties.density *
                                       compute_face_value_upwind(face, component) * face_area;

                    // Convection coefficients
                    if (component == 'u')
                    {
                        momentum_matrix_u[cell_id][cell_id] += std::max(mass_flux, 0.0);
                        momentum_matrix_u[cell_id][neighbor] -= std::max(-mass_flux, 0.0);
                    }
                    else if (component == 'v')
                    {
                        momentum_matrix_v[cell_id][cell_id] += std::max(mass_flux, 0.0);
                        momentum_matrix_v[cell_id][neighbor] -= std::max(-mass_flux, 0.0);
                    }
                    else if (component == 'w' && !is_2d_flow)
                    {
                        momentum_matrix_w[cell_id][cell_id] += std::max(mass_flux, 0.0);
                        momentum_matrix_w[cell_id][neighbor] -= std::max(-mass_flux, 0.0);
                    }
                }
            }
        }
    }

    void IncompressibleFlowSolver::discretize_diffusion_term(int cell_id, char component)
    {
        if (cell_id >= Cells.size())
            return;

        double cell_volume = (Cells[cell_id].Vol > 0.0) ? Cells[cell_id].Vol : Cells[cell_id].Area;

        for (size_t face = 0; face < Cells[cell_id].Face_Areas.size(); ++face)
        {
            if (face < Cells[cell_id].Neighbours.size())
            {
                int neighbor = Cells[cell_id].Neighbours[face];
                double face_area = Cells[cell_id].Face_Areas[face];

                if (neighbor >= 0 && neighbor < num_cells)
                {
                    // Distance between cell centers (simplified)
                    double distance = 1.0; // Should be computed from geometry

                    // Diffusion coefficient
                    double diff_coeff = fluid_properties.viscosity * face_area / distance;

                    if (component == 'u')
                    {
                        momentum_matrix_u[cell_id][cell_id] += diff_coeff;
                        momentum_matrix_u[cell_id][neighbor] -= diff_coeff;
                    }
                    else if (component == 'v')
                    {
                        momentum_matrix_v[cell_id][cell_id] += diff_coeff;
                        momentum_matrix_v[cell_id][neighbor] -= diff_coeff;
                    }
                    else if (component == 'w' && !is_2d_flow)
                    {
                        momentum_matrix_w[cell_id][cell_id] += diff_coeff;
                        momentum_matrix_w[cell_id][neighbor] -= diff_coeff;
                    }
                }
            }
        }
    }

    void IncompressibleFlowSolver::discretize_pressure_gradient(int cell_id, char component)
    {
        if (cell_id >= Cells.size())
            return;

        // Pressure gradient contribution to momentum equation
        for (size_t face = 0; face < Cells[cell_id].Face_Areas.size(); ++face)
        {
            if (face < Cells[cell_id].Neighbours.size())
            {
                int neighbor = Cells[cell_id].Neighbours[face];
                double face_area = Cells[cell_id].Face_Areas[face];

                if (neighbor >= 0 && neighbor < num_cells)
                {
                    double nx = (face * 2 < Cells[cell_id].Face_Normals.size()) ? Cells[cell_id].Face_Normals[face * 2] : 0.0;
                    double ny = (face * 2 + 1 < Cells[cell_id].Face_Normals.size()) ? Cells[cell_id].Face_Normals[face * 2 + 1] : 0.0;

                    double pressure_gradient_contrib = 0.0;

                    if (component == 'u')
                    {
                        pressure_gradient_contrib = nx * face_area;
                        momentum_rhs_u[cell_id] -= pressure_gradient_contrib *
                                                   (pressure_field.p[neighbor] - pressure_field.p[cell_id]);
                    }
                    else if (component == 'v')
                    {
                        pressure_gradient_contrib = ny * face_area;
                        momentum_rhs_v[cell_id] -= pressure_gradient_contrib *
                                                   (pressure_field.p[neighbor] - pressure_field.p[cell_id]);
                    }
                }
            }
        }
    }

    void IncompressibleFlowSolver::discretize_mass_conservation(int cell_id)
    {
        if (cell_id >= Cells.size())
            return;

        // Mass conservation: ∇·(ρu) = 0 becomes ∇·u = 0 for incompressible flow
        double cell_volume = (Cells[cell_id].Vol > 0.0) ? Cells[cell_id].Vol : Cells[cell_id].Area;

        for (size_t face = 0; face < Cells[cell_id].Face_Areas.size(); ++face)
        {
            if (face < Cells[cell_id].Neighbours.size())
            {
                int neighbor = Cells[cell_id].Neighbours[face];
                double face_area = Cells[cell_id].Face_Areas[face];

                if (neighbor >= 0 && neighbor < num_cells)
                {
                    // Velocity divergence from intermediate velocities
                    double u_face = 0.5 * (velocity_field.u_star[cell_id] + velocity_field.u_star[neighbor]);
                    double v_face = 0.5 * (velocity_field.v_star[cell_id] + velocity_field.v_star[neighbor]);

                    double nx = (face * 2 < Cells[cell_id].Face_Normals.size()) ? Cells[cell_id].Face_Normals[face * 2] : 0.0;
                    double ny = (face * 2 + 1 < Cells[cell_id].Face_Normals.size()) ? Cells[cell_id].Face_Normals[face * 2 + 1] : 0.0;

                    double velocity_flux = (u_face * nx + v_face * ny) * face_area;
                    pressure_rhs[cell_id] += fluid_properties.density * velocity_flux / solver_params.dt;

                    // Pressure correction matrix coefficients
                    double distance = 1.0; // Should be computed from geometry
                    double coeff = fluid_properties.density * face_area * face_area /
                                   (solver_params.dt * distance);

                    pressure_matrix[cell_id][cell_id] += coeff;
                    pressure_matrix[cell_id][neighbor] -= coeff;
                }
            }
        }
    }

    // ================================
    // Linear System Solvers
    // ================================

    void IncompressibleFlowSolver::solve_linear_system_velocity(char component)
    {
        if (component == 'u')
        {
            gauss_seidel_solver(momentum_matrix_u, velocity_field.u_star, momentum_rhs_u,
                                solver_params.linear_tolerance, solver_params.max_linear_iterations);
        }
        else if (component == 'v')
        {
            gauss_seidel_solver(momentum_matrix_v, velocity_field.v_star, momentum_rhs_v,
                                solver_params.linear_tolerance, solver_params.max_linear_iterations);
        }
        else if (component == 'w' && !is_2d_flow)
        {
            gauss_seidel_solver(momentum_matrix_w, velocity_field.w_star, momentum_rhs_w,
                                solver_params.linear_tolerance, solver_params.max_linear_iterations);
        }
    }

    void IncompressibleFlowSolver::solve_linear_system_pressure()
    {
        conjugate_gradient_solver(pressure_matrix, pressure_field.p_prime, pressure_rhs,
                                  solver_params.linear_tolerance, solver_params.max_linear_iterations);
    }

    void IncompressibleFlowSolver::gauss_seidel_solver(const std::vector<std::vector<double>> &A,
                                                       std::vector<double> &x,
                                                       const std::vector<double> &b,
                                                       double tolerance, int max_iter)
    {
        int n = x.size();
        std::vector<double> x_old = x;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            double max_change = 0.0;

            for (int i = 0; i < n; ++i)
            {
                double sum = 0.0;
                for (int j = 0; j < n; ++j)
                {
                    if (i != j)
                    {
                        sum += A[i][j] * x[j];
                    }
                }

                double x_new = (b[i] - sum) / std::max(A[i][i], 1e-16);
                max_change = std::max(max_change, std::abs(x_new - x[i]));
                x[i] = x_new;
            }

            if (max_change < tolerance)
            {
                break;
            }
        }
    }

    void IncompressibleFlowSolver::conjugate_gradient_solver(const std::vector<std::vector<double>> &A,
                                                             std::vector<double> &x,
                                                             const std::vector<double> &b,
                                                             double tolerance, int max_iter)
    {
        int n = x.size();
        std::vector<double> r(n), p(n), Ap(n);

        // Initialize residual r = b - Ax
        for (int i = 0; i < n; ++i)
        {
            double Ax_i = 0.0;
            for (int j = 0; j < n; ++j)
            {
                Ax_i += A[i][j] * x[j];
            }
            r[i] = b[i] - Ax_i;
            p[i] = r[i];
        }

        double rsold = 0.0;
        for (int i = 0; i < n; ++i)
        {
            rsold += r[i] * r[i];
        }

        for (int iter = 0; iter < max_iter; ++iter)
        {
            // Ap = A * p
            for (int i = 0; i < n; ++i)
            {
                Ap[i] = 0.0;
                for (int j = 0; j < n; ++j)
                {
                    Ap[i] += A[i][j] * p[j];
                }
            }

            // alpha = rsold / (p^T * Ap)
            double pAp = 0.0;
            for (int i = 0; i < n; ++i)
            {
                pAp += p[i] * Ap[i];
            }

            if (std::abs(pAp) < 1e-16)
                break;

            double alpha = rsold / pAp;

            // x = x + alpha * p
            // r = r - alpha * Ap
            for (int i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }

            double rsnew = 0.0;
            for (int i = 0; i < n; ++i)
            {
                rsnew += r[i] * r[i];
            }

            if (std::sqrt(rsnew) < tolerance)
            {
                break;
            }

            // beta = rsnew / rsold
            double beta = rsnew / rsold;

            // p = r + beta * p
            for (int i = 0; i < n; ++i)
            {
                p[i] = r[i] + beta * p[i];
            }

            rsold = rsnew;
        }
    }

    // ================================
    // Convergence and Monitoring
    // ================================

    bool IncompressibleFlowSolver::check_convergence()
    {
        double vel_residual = compute_velocity_residual();
        double press_residual = compute_pressure_residual();
        double cont_residual = compute_continuity_residual();

        velocity_residuals.push_back(vel_residual);
        pressure_residuals.push_back(press_residual);
        continuity_residuals.push_back(cont_residual);

        bool converged = (vel_residual < solver_params.velocity_tolerance) &&
                         (press_residual < solver_params.pressure_tolerance) &&
                         (cont_residual < solver_params.continuity_tolerance);

        return converged;
    }

    double IncompressibleFlowSolver::compute_velocity_residual()
    {
        double residual = 0.0;
        int count = 0;

        for (int i = 0; i < num_cells; ++i)
        {
            if (i < velocity_field.u.size() && i < velocity_field.u_star.size())
            {
                residual += std::pow(velocity_field.u[i] - velocity_field.u_star[i], 2);
                residual += std::pow(velocity_field.v[i] - velocity_field.v_star[i], 2);
                if (!is_2d_flow && i < velocity_field.w.size() && i < velocity_field.w_star.size())
                {
                    residual += std::pow(velocity_field.w[i] - velocity_field.w_star[i], 2);
                }
                count += is_2d_flow ? 2 : 3;
            }
        }

        return (count > 0) ? std::sqrt(residual / count) : 0.0;
    }

    double IncompressibleFlowSolver::compute_pressure_residual()
    {
        double residual = 0.0;
        int count = 0;

        for (int i = 0; i < num_cells; ++i)
        {
            if (i < pressure_field.p.size() && i < pressure_field.p_old.size())
            {
                residual += std::pow(pressure_field.p[i] - pressure_field.p_old[i], 2);
                count++;
            }
        }

        return (count > 0) ? std::sqrt(residual / count) : 0.0;
    }

    double IncompressibleFlowSolver::compute_continuity_residual()
    {
        double residual = 0.0;
        int count = 0;

        for (int i = 0; i < num_cells; ++i)
        {
            if (i < Cells.size())
            {
                double divergence = 0.0;

                for (size_t face = 0; face < Cells[i].Face_Areas.size(); ++face)
                {
                    if (face < Cells[i].Neighbours.size())
                    {
                        int neighbor = Cells[i].Neighbours[face];
                        double face_area = Cells[i].Face_Areas[face];

                        if (neighbor >= 0 && neighbor < num_cells)
                        {
                            double nx = (face * 2 < Cells[i].Face_Normals.size()) ? Cells[i].Face_Normals[face * 2] : 0.0;
                            double ny = (face * 2 + 1 < Cells[i].Face_Normals.size()) ? Cells[i].Face_Normals[face * 2 + 1] : 0.0;

                            double u_face = 0.5 * (velocity_field.u[i] + velocity_field.u[neighbor]);
                            double v_face = 0.5 * (velocity_field.v[i] + velocity_field.v[neighbor]);

                            divergence += (u_face * nx + v_face * ny) * face_area;
                        }
                    }
                }

                double cell_volume = (Cells[i].Vol > 0.0) ? Cells[i].Vol : Cells[i].Area;
                residual += std::pow(divergence / cell_volume, 2);
                count++;
            }
        }

        return (count > 0) ? std::sqrt(residual / count) : 0.0;
    }

    void IncompressibleFlowSolver::print_iteration_info(int iteration)
    {
        std::cout << "  Iteration " << std::setw(4) << iteration
                  << ":  Vel Res = " << std::scientific << std::setprecision(3) << velocity_residuals.back()
                  << ", Press Res = " << pressure_residuals.back()
                  << ", Cont Res = " << continuity_residuals.back() << std::endl;
    }

    // ================================
    // Utility Methods
    // ================================

    void IncompressibleFlowSolver::interpolate_cell_to_face()
    {
        // Simple averaging for face values
        std::fill(velocity_field.u_face.begin(), velocity_field.u_face.end(), 0.0);
        std::fill(velocity_field.v_face.begin(), velocity_field.v_face.end(), 0.0);
        if (!is_2d_flow)
        {
            std::fill(velocity_field.w_face.begin(), velocity_field.w_face.end(), 0.0);
        }

        // This is a simplified implementation - proper face interpolation would be more complex
        for (int i = 0; i < num_cells; ++i)
        {
            if (i < Cells.size())
            {
                for (size_t face = 0; face < Cells[i].Face_Areas.size() && face < velocity_field.u_face.size(); ++face)
                {
                    if (face < Cells[i].Neighbours.size())
                    {
                        int neighbor = Cells[i].Neighbours[face];
                        if (neighbor >= 0 && neighbor < num_cells)
                        {
                            velocity_field.u_face[face] = 0.5 * (velocity_field.u[i] + velocity_field.u[neighbor]);
                            velocity_field.v_face[face] = 0.5 * (velocity_field.v[i] + velocity_field.v[neighbor]);
                            if (!is_2d_flow)
                            {
                                velocity_field.w_face[face] = 0.5 * (velocity_field.w[i] + velocity_field.w[neighbor]);
                            }
                        }
                    }
                }
            }
        }
    }

    double IncompressibleFlowSolver::compute_face_value_central_difference(int face_id, char component)
    {
        // Simplified - should use proper face indexing
        return 0.0;
    }

    double IncompressibleFlowSolver::compute_face_value_upwind(int face_id, char component)
    {
        // Simplified - should use proper upwind scheme
        return 0.0;
    }

    // ================================
    // Boundary Conditions (Simplified)
    // ================================

    void IncompressibleFlowSolver::apply_velocity_boundary_conditions()
    {
        // Apply boundary conditions based on cell classification
        // This is a simplified implementation

        for (int cell_id : Wall_Cells_List)
        {
            if (cell_id < num_cells)
            {
                apply_wall_boundary_condition(cell_id);
            }
        }

        for (int cell_id : Inlet_Cells_List)
        {
            if (cell_id < num_cells)
            {
                apply_inlet_boundary_condition(cell_id);
            }
        }

        for (int cell_id : Exit_Cells_List)
        {
            if (cell_id < num_cells)
            {
                apply_outlet_boundary_condition(cell_id);
            }
        }
    }

    void IncompressibleFlowSolver::apply_pressure_boundary_conditions()
    {
        // Apply pressure boundary conditions
        // Simplified implementation
    }

    void IncompressibleFlowSolver::apply_wall_boundary_condition(int cell_id)
    {
        // No-slip wall: u = v = w = 0
        if (cell_id < velocity_field.u.size())
        {
            velocity_field.u[cell_id] = 0.0;
            velocity_field.v[cell_id] = 0.0;
            if (!is_2d_flow && cell_id < velocity_field.w.size())
            {
                velocity_field.w[cell_id] = 0.0;
            }
        }
    }

    void IncompressibleFlowSolver::apply_inlet_boundary_condition(int cell_id)
    {
        // Specified velocity inlet - would need proper BC values
        if (cell_id < velocity_field.u.size())
        {
            velocity_field.u[cell_id] = 1.0; // Example inlet velocity
            velocity_field.v[cell_id] = 0.0;
            if (!is_2d_flow && cell_id < velocity_field.w.size())
            {
                velocity_field.w[cell_id] = 0.0;
            }
        }
    }

    void IncompressibleFlowSolver::apply_outlet_boundary_condition(int cell_id)
    {
        // Zero gradient or specified pressure outlet
        // Simplified implementation
    }

    void IncompressibleFlowSolver::apply_symmetry_boundary_condition(int cell_id)
    {
        // Symmetry: normal velocity = 0, tangential velocity free
        // Simplified implementation
    }

    // ================================
    // Output Methods (Simplified)
    // ================================

    void IncompressibleFlowSolver::write_solution_vtk(const std::string &filename)
    {
        std::cout << "Writing VTK solution to " << filename << std::endl;

        std::ofstream vtk_file(filename);
        if (!vtk_file.is_open())
        {
            std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
            return;
        }

        // VTK header
        vtk_file << "# vtk DataFile Version 3.0" << std::endl;
        vtk_file << "Incompressible Flow Solution" << std::endl;
        vtk_file << "ASCII" << std::endl;
        vtk_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

        // Points (simplified - using cell centers)
        vtk_file << "POINTS " << num_cells << " float" << std::endl;
        for (int i = 0; i < num_cells && i < Cells.size(); ++i)
        {
            if (Cells[i].Cell_Center.size() >= 3)
            {
                vtk_file << Cells[i].Cell_Center[0] << " "
                         << Cells[i].Cell_Center[1] << " "
                         << Cells[i].Cell_Center[2] << std::endl;
            }
            else
            {
                vtk_file << "0.0 0.0 0.0" << std::endl;
            }
        }

        // Point data
        vtk_file << "POINT_DATA " << num_cells << std::endl;

        // Velocity
        vtk_file << "VECTORS Velocity float" << std::endl;
        for (int i = 0; i < num_cells; ++i)
        {
            double u = (i < velocity_field.u.size()) ? velocity_field.u[i] : 0.0;
            double v = (i < velocity_field.v.size()) ? velocity_field.v[i] : 0.0;
            double w = (!is_2d_flow && i < velocity_field.w.size()) ? velocity_field.w[i] : 0.0;
            vtk_file << u << " " << v << " " << w << std::endl;
        }

        // Pressure
        vtk_file << "SCALARS Pressure float 1" << std::endl;
        vtk_file << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < num_cells; ++i)
        {
            double p = (i < pressure_field.p.size()) ? pressure_field.p[i] : 0.0;
            vtk_file << p << std::endl;
        }

        vtk_file.close();
        std::cout << "VTK solution written successfully." << std::endl;
    }

    void IncompressibleFlowSolver::write_residual_history(const std::string &filename)
    {
        std::cout << "Writing residual history to " << filename << std::endl;

        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
            return;
        }

        file << "# Iteration Velocity_Residual Pressure_Residual Continuity_Residual" << std::endl;

        size_t num_iters = std::min({velocity_residuals.size(), pressure_residuals.size(), continuity_residuals.size()});

        for (size_t i = 0; i < num_iters; ++i)
        {
            file << i + 1 << " "
                 << std::scientific << std::setprecision(6)
                 << velocity_residuals[i] << " "
                 << pressure_residuals[i] << " "
                 << continuity_residuals[i] << std::endl;
        }

        file.close();
        std::cout << "Residual history written successfully." << std::endl;
    }

    void IncompressibleFlowSolver::write_solution_tecplot(const std::string &filename)
    {
        std::cout << "Tecplot output not implemented yet." << std::endl;
    }

    void IncompressibleFlowSolver::write_monitor_points(const std::string &filename)
    {
        std::cout << "Monitor points output not implemented yet." << std::endl;
    }

} // namespace IncompressibleSolver