#ifndef INCOMPRESSIBLE_SOLVER_STANDALONE_H
#define INCOMPRESSIBLE_SOLVER_STANDALONE_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cassert>

// Cell structure definition - standalone version
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

// Boundary condition types
enum BoundaryType
{
    BC_WALL = 1,
    BC_INLET = 2,
    BC_OUTLET = 3,
    BC_SYMMETRY = 4
};

// Velocity field structure
struct VelocityField
{
    std::vector<double> u, v, w;             // Velocity components at cell centers
    std::vector<double> u_old, v_old, w_old; // Previous time step velocities

    void resize(int numCells)
    {
        u.resize(numCells);
        v.resize(numCells);
        w.resize(numCells);
        u_old.resize(numCells);
        v_old.resize(numCells);
        w_old.resize(numCells);
    }

    void initialize(int numCells, double u_init = 0.0, double v_init = 0.0, double w_init = 0.0)
    {
        resize(numCells);
        std::fill(u.begin(), u.end(), u_init);
        std::fill(v.begin(), v.end(), v_init);
        std::fill(w.begin(), w.end(), w_init);
        u_old = u;
        v_old = v;
        w_old = w;
    }
};

// Pressure field structure
struct PressureField
{
    std::vector<double> p;            // Pressure at cell centers
    std::vector<double> p_old;        // Previous time step pressure
    std::vector<double> p_correction; // Pressure correction for SIMPLE

    void resize(int numCells)
    {
        p.resize(numCells);
        p_old.resize(numCells);
        p_correction.resize(numCells);
    }

    void initialize(int numCells, double p_init = 0.0)
    {
        resize(numCells);
        std::fill(p.begin(), p.end(), p_init);
        std::fill(p_old.begin(), p_old.end(), p_init);
        std::fill(p_correction.begin(), p_correction.end(), 0.0);
    }
};

// Fluid properties
struct FluidProperties
{
    double density;             // kg/m³
    double viscosity;           // Pa⋅s (dynamic viscosity)
    double kinematic_viscosity; // m²/s (= viscosity/density)

    FluidProperties(double rho = 1.0, double mu = 1e-3)
        : density(rho), viscosity(mu), kinematic_viscosity(mu / rho) {}

    void update()
    {
        kinematic_viscosity = viscosity / density;
    }

    double reynolds_number(double characteristic_velocity, double characteristic_length) const
    {
        return (density * characteristic_velocity * characteristic_length) / viscosity;
    }
};

// Solver parameters
struct SolverParameters
{
    double dt;                  // Time step
    int max_iterations;         // Maximum iterations per time step
    double tolerance;           // Convergence tolerance
    double pressure_relaxation; // Under-relaxation factor for pressure
    double velocity_relaxation; // Under-relaxation factor for velocity
    bool steady_state;          // Steady or unsteady flow
    int max_time_steps;         // Maximum time steps for unsteady flow
    double final_time;          // Final simulation time

    SolverParameters() : dt(0.001), max_iterations(1000), tolerance(1e-6),
                         pressure_relaxation(0.3), velocity_relaxation(0.7),
                         steady_state(true), max_time_steps(10000), final_time(1.0) {}
};

// Boundary condition data
struct BoundaryCondition
{
    BoundaryType type;
    double u_value, v_value, w_value; // Velocity values (for Dirichlet BC)
    double pressure_value;            // Pressure value

    BoundaryCondition(BoundaryType t = BC_WALL, double u = 0.0, double v = 0.0, double w = 0.0, double p = 0.0)
        : type(t), u_value(u), v_value(v), w_value(w), pressure_value(p) {}
};

// Main incompressible flow solver class
class IncompressibleFlowSolver
{
private:
    // Grid data
    std::vector<Cell> *cells_ptr;
    int num_cells;
    bool is_2d;

    // Flow fields
    VelocityField velocity;
    PressureField pressure;

    // Physical properties
    FluidProperties fluid;
    SolverParameters params;

    // Boundary conditions
    std::vector<BoundaryCondition> boundary_conditions;
    std::vector<int> wall_cells, inlet_cells, outlet_cells, symmetry_cells;

    // Linear system matrices (stored as coefficient arrays)
    std::vector<std::vector<double>> momentum_matrix_u, momentum_matrix_v, momentum_matrix_w;
    std::vector<std::vector<double>> pressure_matrix;
    std::vector<double> momentum_rhs_u, momentum_rhs_v, momentum_rhs_w;
    std::vector<double> pressure_rhs;

    // Internal methods
    void setup_boundary_cells();
    void discretize_momentum_equations();
    void discretize_pressure_correction();
    void solve_linear_system(const std::vector<std::vector<double>> &matrix,
                             const std::vector<double> &rhs, std::vector<double> &solution);
    void apply_boundary_conditions();
    void update_velocity_field();
    void update_pressure_field();
    bool check_convergence();

    // SIMPLE algorithm steps
    void solve_momentum_predictor();
    void solve_pressure_correction();
    void correct_velocity_field();
    void correct_pressure_field();

public:
    IncompressibleFlowSolver();
    ~IncompressibleFlowSolver() = default;

    // Initialization methods
    bool initialize(std::vector<Cell> &grid_cells, bool is_2d_flow = true);
    void set_fluid_properties(double density, double viscosity);
    void set_solver_parameters(const SolverParameters &solver_params);
    void set_boundary_condition(int cell_id, const BoundaryCondition &bc);
    void set_initial_conditions(double u_init = 0.0, double v_init = 0.0, double w_init = 0.0, double p_init = 0.0);

    // Solver methods
    bool solve_steady_state();
    bool solve_unsteady();
    bool solve_single_time_step();

    // Output and analysis methods
    void write_solution(const std::string &filename, int iteration = 0);
    void write_vtk_output(const std::string &filename);
    double calculate_residual();
    void print_solver_info();

    // Getter methods
    const VelocityField &get_velocity_field() const { return velocity; }
    const PressureField &get_pressure_field() const { return pressure; }
    const FluidProperties &get_fluid_properties() const { return fluid; }
    const SolverParameters &get_solver_parameters() const { return params; }
};

// External global variables (normally from Grid_Computations.cpp)
extern std::vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;
extern std::vector<int> Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
extern std::vector<double> Vertices;
extern int Total_No_Cells, No_Physical_Cells;
extern bool Is_2D_Flow;

#endif // INCOMPRESSIBLE_SOLVER_STANDALONE_H