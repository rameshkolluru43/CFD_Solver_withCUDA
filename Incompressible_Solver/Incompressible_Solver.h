#ifndef INCOMPRESSIBLE_SOLVER_H
#define INCOMPRESSIBLE_SOLVER_H

#include "definitions.h"
#include <vector>
#include <string>
#include <memory>

/**
 * @file Incompressible_Solver.h
 * @brief Header file for cell-centered staggered grid based incompressible solver using finite volume method
 *
 * This solver implements the SIMPLE algorithm with cell-centered staggered grid approach for
 * incompressible Navier-Stokes equations. It uses the same grid functions from the src folder
 * for consistency with the existing compressible solver framework.
 *
 * Key Features:
 * - Cell-centered staggered grid arrangement
 * - SIMPLE algorithm for pressure-velocity coupling
 * - Finite volume discretization
 * - Support for both 2D and 3D domains
 * - Integration with existing grid infrastructure
 * - CUDA-ready data structures
 */

namespace IncompressibleSolver
{

    // ================================
    // Data Structures
    // ================================

    /**
     * @struct VelocityField
     * @brief Stores velocity components for the incompressible solver
     */
    struct VelocityField
    {
        std::vector<double> u;      // U-velocity component at cell centers
        std::vector<double> v;      // V-velocity component at cell centers
        std::vector<double> w;      // W-velocity component at cell centers
        std::vector<double> u_star; // Intermediate U-velocity (predictor step)
        std::vector<double> v_star; // Intermediate V-velocity (predictor step)
        std::vector<double> w_star; // Intermediate W-velocity (predictor step)
        std::vector<double> u_face; // Face-centered U-velocity
        std::vector<double> v_face; // Face-centered V-velocity
        std::vector<double> w_face; // Face-centered W-velocity

        void resize(int numCells, int numFaces)
        {
            u.resize(numCells, 0.0);
            v.resize(numCells, 0.0);
            w.resize(numCells, 0.0);
            u_star.resize(numCells, 0.0);
            v_star.resize(numCells, 0.0);
            w_star.resize(numCells, 0.0);
            u_face.resize(numFaces, 0.0);
            v_face.resize(numFaces, 0.0);
            w_face.resize(numFaces, 0.0);
        }
    };

    /**
     * @struct PressureField
     * @brief Stores pressure and pressure correction fields
     */
    struct PressureField
    {
        std::vector<double> p;       // Pressure at cell centers
        std::vector<double> p_prime; // Pressure correction
        std::vector<double> p_old;   // Previous iteration pressure

        void resize(int numCells)
        {
            p.resize(numCells, 0.0);
            p_prime.resize(numCells, 0.0);
            p_old.resize(numCells, 0.0);
        }
    };

    /**
     * @struct FluidProperties
     * @brief Physical properties of the incompressible fluid
     */
    struct FluidProperties
    {
        double density;             // Fluid density (constant)
        double viscosity;           // Dynamic viscosity
        double kinematic_viscosity; // Kinematic viscosity (mu/rho)

        FluidProperties() : density(1.0), viscosity(1e-3), kinematic_viscosity(1e-3) {}

        void update_kinematic_viscosity()
        {
            kinematic_viscosity = viscosity / density;
        }
    };

    /**
     * @struct SolverParameters
     * @brief Control parameters for the incompressible solver
     */
    struct SolverParameters
    {
        // Time stepping
        double dt;          // Time step size
        double max_time;    // Maximum simulation time
        int max_iterations; // Maximum number of iterations
        bool steady_state;  // Steady state flag

        // Convergence criteria
        double velocity_tolerance;   // Velocity convergence tolerance
        double pressure_tolerance;   // Pressure convergence tolerance
        double continuity_tolerance; // Mass conservation tolerance

        // Under-relaxation factors
        double alpha_u; // U-velocity under-relaxation
        double alpha_v; // V-velocity under-relaxation
        double alpha_w; // W-velocity under-relaxation
        double alpha_p; // Pressure under-relaxation

        // Linear solver parameters
        int max_linear_iterations; // Maximum linear solver iterations
        double linear_tolerance;   // Linear solver tolerance

        // Output control
        int output_frequency;      // Output frequency
        std::string output_format; // Output format (VTK, Tecplot, etc.)

        // Default constructor with reasonable values
        SolverParameters() : dt(1e-4), max_time(1.0), max_iterations(1000), steady_state(true),
                             velocity_tolerance(1e-6), pressure_tolerance(1e-6), continuity_tolerance(1e-8),
                             alpha_u(0.7), alpha_v(0.7), alpha_w(0.7), alpha_p(0.3),
                             max_linear_iterations(1000), linear_tolerance(1e-9),
                             output_frequency(100), output_format("VTK") {}
    };

    /**
     * @struct BoundaryCondition
     * @brief Boundary condition specification for incompressible flow
     */
    struct BoundaryCondition
    {
        enum Type
        {
            WALL,             // No-slip wall
            INLET_VELOCITY,   // Specified velocity inlet
            INLET_MASS_FLOW,  // Specified mass flow inlet
            OUTLET_PRESSURE,  // Specified pressure outlet
            OUTLET_ZERO_GRAD, // Zero gradient outlet
            SYMMETRY,         // Symmetry boundary
            PERIODIC          // Periodic boundary
        };

        Type type;
        double u_value;   // Specified U-velocity
        double v_value;   // Specified V-velocity
        double w_value;   // Specified W-velocity
        double p_value;   // Specified pressure
        double mass_flow; // Specified mass flow rate

        BoundaryCondition() : type(WALL), u_value(0.0), v_value(0.0), w_value(0.0),
                              p_value(0.0), mass_flow(0.0) {}
    };

    // ================================
    // Class Definitions
    // ================================

    /**
     * @class IncompressibleSolver
     * @brief Main class for incompressible flow solver using cell-centered staggered grid
     */
    class IncompressibleFlowSolver
    {
    private:
        // Grid and geometry data (using existing infrastructure)
        int num_cells;
        int num_faces;
        int num_boundary_faces;
        bool is_2d_flow;

        // Flow fields
        VelocityField velocity_field;
        PressureField pressure_field;

        // Physical properties
        FluidProperties fluid_properties;

        // Solver parameters
        SolverParameters solver_params;

        // Boundary conditions
        std::vector<BoundaryCondition> boundary_conditions;

        // Coefficient matrices and linear system
        std::vector<std::vector<double>> momentum_matrix_u;
        std::vector<std::vector<double>> momentum_matrix_v;
        std::vector<std::vector<double>> momentum_matrix_w;
        std::vector<std::vector<double>> pressure_matrix;
        std::vector<double> momentum_rhs_u;
        std::vector<double> momentum_rhs_v;
        std::vector<double> momentum_rhs_w;
        std::vector<double> pressure_rhs;

        // Convergence monitoring
        std::vector<double> velocity_residuals;
        std::vector<double> pressure_residuals;
        std::vector<double> continuity_residuals;

    public:
        // Constructor and destructor
        IncompressibleFlowSolver();
        ~IncompressibleFlowSolver();

        // Initialization methods
        void initialize_from_grid();
        void set_fluid_properties(double density, double viscosity);
        void set_solver_parameters(const SolverParameters &params);
        void set_boundary_conditions(const std::vector<BoundaryCondition> &bcs);
        void initialize_fields();

        // Main solver methods
        void solve();
        void solve_steady_state();
        void solve_transient();

        // SIMPLE algorithm components
        void momentum_predictor_step();
        void pressure_correction_step();
        void velocity_corrector_step();
        void update_face_velocities();

        // Discretization methods
        void assemble_momentum_equations();
        void assemble_pressure_equation();
        void discretize_convection_term(int cell_id, char component);
        void discretize_diffusion_term(int cell_id, char component);
        void discretize_pressure_gradient(int cell_id, char component);
        void discretize_mass_conservation(int cell_id);

        // Boundary condition application
        void apply_velocity_boundary_conditions();
        void apply_pressure_boundary_conditions();
        void apply_wall_boundary_condition(int cell_id);
        void apply_inlet_boundary_condition(int cell_id);
        void apply_outlet_boundary_condition(int cell_id);
        void apply_symmetry_boundary_condition(int cell_id);

        // Linear system solvers
        void solve_linear_system_velocity(char component);
        void solve_linear_system_pressure();
        void gauss_seidel_solver(const std::vector<std::vector<double>> &A,
                                 std::vector<double> &x,
                                 const std::vector<double> &b,
                                 double tolerance, int max_iter);
        void conjugate_gradient_solver(const std::vector<std::vector<double>> &A,
                                       std::vector<double> &x,
                                       const std::vector<double> &b,
                                       double tolerance, int max_iter);

        // Convergence and monitoring
        bool check_convergence();
        double compute_velocity_residual();
        double compute_pressure_residual();
        double compute_continuity_residual();
        void print_iteration_info(int iteration);

        // Output methods
        void write_solution_vtk(const std::string &filename);
        void write_solution_tecplot(const std::string &filename);
        void write_residual_history(const std::string &filename);
        void write_monitor_points(const std::string &filename);

        // Utility methods
        void interpolate_cell_to_face();
        void compute_face_fluxes();
        void compute_cell_gradients();
        double compute_face_value_central_difference(int face_id, char component);
        double compute_face_value_upwind(int face_id, char component);

        // Access methods
        const VelocityField &get_velocity_field() const { return velocity_field; }
        const PressureField &get_pressure_field() const { return pressure_field; }
        const FluidProperties &get_fluid_properties() const { return fluid_properties; }
        const SolverParameters &get_solver_parameters() const { return solver_params; }
    };

    // ================================
    // Utility Functions
    // ================================

    /**
     * @brief Initialize incompressible solver with grid data from existing infrastructure
     */
    void initialize_incompressible_solver_with_grid(IncompressibleFlowSolver &solver);

    /**
     * @brief Create default boundary conditions based on existing boundary classifications
     */
    std::vector<BoundaryCondition> create_default_boundary_conditions();

    /**
     * @brief Read solver configuration from JSON file
     */
    SolverParameters read_solver_parameters_from_json(const std::string &filename);

    /**
     * @brief Write solver configuration to JSON file
     */
    void write_solver_parameters_to_json(const SolverParameters &params, const std::string &filename);

    /**
     * @brief Validate grid suitability for incompressible solver
     */
    bool validate_grid_for_incompressible_solver();

    /**
     * @brief Compute Reynolds numbers for the flow
     */
    void compute_flow_reynolds_numbers(const FluidProperties &props,
                                       const VelocityField &vel_field,
                                       double characteristic_length,
                                       double &Re_global, double &Re_cell);

} // namespace IncompressibleSolver

#endif // INCOMPRESSIBLE_SOLVER_H