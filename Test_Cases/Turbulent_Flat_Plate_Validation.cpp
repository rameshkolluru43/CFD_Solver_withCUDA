/**
 * @file Turbulent_Flat_Plate_Validation.cpp
 * @brief Validation test case for turbulent flow over a flat plate
 * @author Ramesh Kolluru
 * @date 2025-01-XX
 *
 * This file implements a validation test case for turbulent boundary layer
 * development over a flat plate. The case validates both K-epsilon and
 * K-omega turbulence models against theoretical and experimental data.
 *
 * Test Case Details:
 * - 2D turbulent boundary layer over flat plate
 * - Reynolds number: Re_L = 5.0e6 (based on plate length)
 * - Mach number: M = 0.3 (low-speed compressible)
 * - Plate length: L = 2.0 m
 * - Inlet turbulence intensity: Tu = 0.1%
 *
 * Expected Results:
 * - Skin friction coefficient: Cf = 0.00332 (theoretical)
 * - Boundary layer thickness: δ/L = 0.382/Re_L^(1/5)
 * - Displacement thickness: δ*/
L = 0.0463 / Re_L ^ (1 / 5) * /

#include "Turbulence_Models.h"
#include "Test_Cases.h"
#include "Solver.h"
#include "IO_Write.h"
#include <cmath>
#include <fstream>
#include <iomanip>

                        //=============================================================================
                        // VALIDATION PARAMETERS
                        //=============================================================================

                        struct FlatPlateValidationParams
{
    double plate_length = 2.0;           // Plate length (m)
    double plate_height = 0.5;           // Domain height (m)
    double inlet_velocity = 100.0;       // Inlet velocity (m/s)
    double inlet_pressure = 101325.0;    // Inlet pressure (Pa)
    double inlet_temperature = 288.15;   // Inlet temperature (K)
    double reynolds_number = 5.0e6;      // Reynolds number based on plate length
    double mach_number = 0.3;            // Mach number
    double turbulence_intensity = 0.001; // 0.1% turbulence intensity
    double viscous_length_scale = 0.1;   // Turbulent length scale ratio

    // Grid parameters
    int nx = 200;                  // Points in x-direction
    int ny = 100;                  // Points in y-direction
    double wall_spacing = 1e-5;    // First cell height at wall
    double stretching_ratio = 1.2; // Grid stretching ratio

    // Expected results (from correlations)
    double expected_cf = 0.00332;        // Skin friction coefficient
    double expected_delta_star = 0.0463; // Displacement thickness ratio
    double expected_theta = 0.0365;      // Momentum thickness ratio
};

//=============================================================================
// MAIN VALIDATION FUNCTION
//=============================================================================

/**
 * @brief Run turbulent flat plate validation test
 * @param model Turbulence model to validate
 * @return True if validation passes, false otherwise
 */
bool Run_Turbulent_Flat_Plate_Validation(TurbulenceModel model)
{
    cout << "========================================" << endl;
    cout << "TURBULENT FLAT PLATE VALIDATION TEST" << endl;
    cout << "========================================" << endl;

    FlatPlateValidationParams params;

    // Print test parameters
    Print_Validation_Parameters(params, model);

    // Generate flat plate grid
    Generate_Flat_Plate_Grid(params);

    // Initialize turbulence model
    Initialize_Turbulence_Model(model);

    // Set boundary conditions
    Set_Flat_Plate_Boundary_Conditions(params);

    // Set initial conditions
    Set_Flat_Plate_Initial_Conditions(params);

    // Run simulation
    bool converged = Run_Flat_Plate_Simulation(params);

    if (!converged)
    {
        cout << "ERROR: Simulation did not converge!" << endl;
        return false;
    }

    // Calculate and validate results
    bool validation_passed = Validate_Flat_Plate_Results(params, model);

    // Write detailed results
    Write_Flat_Plate_Results(params, model);

    cout << "========================================" << endl;
    cout << "VALIDATION " << (validation_passed ? "PASSED" : "FAILED") << endl;
    cout << "========================================" << endl;

    return validation_passed;
}

//=============================================================================
// GRID GENERATION
//=============================================================================

/**
 * @brief Generate structured grid for flat plate case
 * @param params Validation parameters
 */
void Generate_Flat_Plate_Grid(const FlatPlateValidationParams &params)
{
    cout << "Generating flat plate grid..." << endl;

    // Clear existing grid
    Cells.clear();

    // Calculate total number of cells
    Total_Cells = params.nx * params.ny;
    Cells.resize(Total_Cells);

    // Generate grid points
    vector<vector<Point>> grid_points(params.nx + 1, vector<Point>(params.ny + 1));

    // X-direction: uniform spacing
    double dx = params.plate_length / params.nx;

    // Y-direction: stretched spacing near wall
    vector<double> y_coords(params.ny + 1);
    y_coords[0] = 0.0; // Wall
    y_coords[1] = params.wall_spacing;

    for (int j = 2; j <= params.ny; j++)
    {
        y_coords[j] = y_coords[j - 1] + (y_coords[j - 1] - y_coords[j - 2]) * params.stretching_ratio;
        if (y_coords[j] > params.plate_height)
        {
            // Adjust remaining points to reach domain height
            double remaining_height = params.plate_height - y_coords[j - 1];
            int remaining_points = params.ny + 1 - j;
            double uniform_spacing = remaining_height / remaining_points;

            for (int k = j; k <= params.ny; k++)
            {
                y_coords[k] = y_coords[j - 1] + (k - j + 1) * uniform_spacing;
            }
            break;
        }
    }

    // Create grid points
    for (int i = 0; i <= params.nx; i++)
    {
        for (int j = 0; j <= params.ny; j++)
        {
            double x = i * dx;
            double y = y_coords[j];
            grid_points[i][j] = Point(x, y, 0.0);
        }
    }

    // Create cells
    int cell_index = 0;
    for (int i = 0; i < params.nx; i++)
    {
        for (int j = 0; j < params.ny; j++)
        {
            // Cell vertices (counter-clockwise)
            vector<Point> vertices(4);
            vertices[0] = grid_points[i][j];         // Bottom-left
            vertices[1] = grid_points[i + 1][j];     // Bottom-right
            vertices[2] = grid_points[i + 1][j + 1]; // Top-right
            vertices[3] = grid_points[i][j + 1];     // Top-left

            // Create cell
            Cells[cell_index].cellID = cell_index;
            Cells[cell_index].cellType = INTERIOR_CELL;
            Cells[cell_index].Dimension = 2;
            Cells[cell_index].No_of_Faces = 4;

            // Calculate cell center
            double cx = 0.25 * (vertices[0].Get_x() + vertices[1].Get_x() +
                                vertices[2].Get_x() + vertices[3].Get_x());
            double cy = 0.25 * (vertices[0].Get_y() + vertices[1].Get_y() +
                                vertices[2].Get_y() + vertices[3].Get_y());
            Cells[cell_index].Cell_Center[0] = cx;
            Cells[cell_index].Cell_Center[1] = cy;
            Cells[cell_index].Cell_Center[2] = 0.0;

            // Calculate cell area (for 2D quadrilateral)
            double area = Calculate_Quadrilateral_Area(vertices);
            Cells[cell_index].Area = area;
            Cells[cell_index].Inv_Area = 1.0 / area;
            Cells[cell_index].Volume = area; // 2D case

            // Set boundary conditions
            Cells[cell_index].hasBoundaryface = false;

            if (j == 0)
            { // Wall boundary
                Cells[cell_index].hasBoundaryface = true;
                Cells[cell_index].cellType = WALL_CELL;
            }
            else if (i == 0)
            { // Inlet boundary
                Cells[cell_index].hasBoundaryface = true;
                Cells[cell_index].cellType = INLET_CELL;
            }
            else if (i == params.nx - 1)
            { // Outlet boundary
                Cells[cell_index].hasBoundaryface = true;
                Cells[cell_index].cellType = OUTLET_CELL;
            }
            else if (j == params.ny - 1)
            { // Top boundary (symmetry)
                Cells[cell_index].hasBoundaryface = true;
                Cells[cell_index].cellType = SYMMETRY_CELL;
            }

            cell_index++;
        }
    }

    // Set up connectivity
    Setup_Flat_Plate_Connectivity(params);

    cout << "Grid generated: " << params.nx << " x " << params.ny << " cells" << endl;
    cout << "First cell height: " << params.wall_spacing << " m" << endl;
}

/**
 * @brief Calculate area of quadrilateral cell
 * @param vertices Vector of cell vertices
 * @return Cell area
 */
double Calculate_Quadrilateral_Area(const vector<Point> &vertices)
{
    // Using shoelace formula for quadrilateral
    double area = 0.0;
    int n = vertices.size();

    for (int i = 0; i < n; i++)
    {
        int j = (i + 1) % n;
        area += vertices[i].Get_x() * vertices[j].Get_y();
        area -= vertices[j].Get_x() * vertices[i].Get_y();
    }

    return abs(area) * 0.5;
}

/**
 * @brief Setup cell connectivity for flat plate grid
 * @param params Validation parameters
 */
void Setup_Flat_Plate_Connectivity(const FlatPlateValidationParams &params)
{
    cout << "Setting up cell connectivity..." << endl;

    for (int i = 0; i < params.nx; i++)
    {
        for (int j = 0; j < params.ny; j++)
        {
            int cell_index = i * params.ny + j;

            // Initialize neighbors
            Cells[cell_index].Neighbours.resize(4, -1);
            Cells[cell_index].Face_Areas.resize(4, 0.0);
            Cells[cell_index].Face_Normals.resize(4, vector<double>(3, 0.0));

            // Left neighbor (face 0)
            if (i > 0)
            {
                Cells[cell_index].Neighbours[0] = (i - 1) * params.ny + j;
            }

            // Right neighbor (face 1)
            if (i < params.nx - 1)
            {
                Cells[cell_index].Neighbours[1] = (i + 1) * params.ny + j;
            }

            // Bottom neighbor (face 2)
            if (j > 0)
            {
                Cells[cell_index].Neighbours[2] = i * params.ny + (j - 1);
            }

            // Top neighbor (face 3)
            if (j < params.ny - 1)
            {
                Cells[cell_index].Neighbours[3] = i * params.ny + (j + 1);
            }

            // Calculate face areas and normals
            Calculate_Face_Properties(cell_index);
        }
    }

    cout << "Cell connectivity established." << endl;
}

//=============================================================================
// BOUNDARY AND INITIAL CONDITIONS
//=============================================================================

/**
 * @brief Set boundary conditions for flat plate validation
 * @param params Validation parameters
 */
void Set_Flat_Plate_Boundary_Conditions(const FlatPlateValidationParams &params)
{
    cout << "Setting boundary conditions..." << endl;

    // Set inlet conditions
    Inlet_Condition.u = params.inlet_velocity;
    Inlet_Condition.v = 0.0;
    Inlet_Condition.P = params.inlet_pressure;
    Inlet_Condition.T = params.inlet_temperature;
    Inlet_Condition.Rho = params.inlet_pressure / (R_GC * params.inlet_temperature);

    // Set reference conditions
    reference_pressure = params.inlet_pressure;
    reference_temperature = params.inlet_temperature;
    reference_density = Inlet_Condition.Rho;
    reference_velocity = params.inlet_velocity;

    // Set characteristic length
    characteristic_length = params.plate_length;

    cout << "Boundary conditions set:" << endl;
    cout << "  Inlet velocity: " << params.inlet_velocity << " m/s" << endl;
    cout << "  Inlet pressure: " << params.inlet_pressure << " Pa" << endl;
    cout << "  Inlet temperature: " << params.inlet_temperature << " K" << endl;
    cout << "  Reynolds number: " << params.reynolds_number << endl;
}

//=============================================================================
// RESULTS VALIDATION
//=============================================================================

/**
 * @brief Validate flat plate results against theoretical values
 * @param params Validation parameters
 * @param model Turbulence model used
 * @return True if validation passes
 */
bool Validate_Flat_Plate_Results(const FlatPlateValidationParams &params, TurbulenceModel model)
{
    cout << "Validating results..." << endl;

    // Calculate skin friction coefficient at x = 0.8L
    double x_location = 0.8 * params.plate_length;
    double calculated_cf = Calculate_Skin_Friction_At_Location(x_location);

    // Calculate boundary layer parameters
    BoundaryLayerProperties bl_props = Calculate_Boundary_Layer_Properties(x_location);

    // Compare with theoretical values
    double cf_error = abs(calculated_cf - params.expected_cf) / params.expected_cf;
    double delta_star_error = abs(bl_props.delta_star_ratio - params.expected_delta_star) / params.expected_delta_star;

    cout << "=== VALIDATION RESULTS ===" << endl;
    cout << "Skin friction coefficient:" << endl;
    cout << "  Calculated: " << calculated_cf << endl;
    cout << "  Expected:   " << params.expected_cf << endl;
    cout << "  Error:      " << cf_error * 100.0 << "%" << endl;

    cout << "Displacement thickness ratio:" << endl;
    cout << "  Calculated: " << bl_props.delta_star_ratio << endl;
    cout << "  Expected:   " << params.expected_delta_star << endl;
    cout << "  Error:      " << delta_star_error * 100.0 << "%" << endl;

    cout << "Boundary layer thickness: " << bl_props.delta_99 << " m" << endl;
    cout << "Momentum thickness: " << bl_props.theta << " m" << endl;

    // Validation criteria (within 10% of theoretical values)
    bool cf_valid = cf_error < 0.10;
    bool delta_star_valid = delta_star_error < 0.15;

    bool validation_passed = cf_valid && delta_star_valid;

    if (!validation_passed)
    {
        cout << "VALIDATION FAILED:" << endl;
        if (!cf_valid)
            cout << "  - Skin friction coefficient error too large" << endl;
        if (!delta_star_valid)
            cout << "  - Displacement thickness error too large" << endl;
    }

    return validation_passed;
}

/**
 * @brief Calculate skin friction coefficient at specific location
 * @param x_location X-coordinate location
 * @return Skin friction coefficient
 */
double Calculate_Skin_Friction_At_Location(double x_location)
{
    // Find cells at the specified x-location near the wall
    double min_distance = 1e10;
    int target_cell = -1;

    for (int i = 0; i < Total_Cells; i++)
    {
        if (Cells[i].cellType == WALL_CELL ||
            (Cells[i].Cell_Center[1] < 0.01 && abs(Cells[i].Cell_Center[0] - x_location) < min_distance))
        {
            min_distance = abs(Cells[i].Cell_Center[0] - x_location);
            target_cell = i;
        }
    }

    if (target_cell == -1)
        return 0.0;

    // Calculate wall shear stress
    Calculate_Wall_Shear_Stress(target_cell);
    double tau_wall = turbulence_vars[target_cell].tau_wall;

    // Calculate skin friction coefficient: Cf = τw / (0.5 * ρ * U²)
    double rho_inf = reference_density;
    double U_inf = reference_velocity;
    double cf = tau_wall / (0.5 * rho_inf * U_inf * U_inf);

    return cf;
}

/**
 * @brief Structure to hold boundary layer properties
 */
struct BoundaryLayerProperties
{
    double delta_99;         // 99% boundary layer thickness
    double delta_star;       // Displacement thickness
    double theta;            // Momentum thickness
    double delta_star_ratio; // δ*/L
    double theta_ratio;      // θ/L
    double shape_factor;     // H = δ*/θ
};

/**
 * @brief Calculate boundary layer properties at specific location
 * @param x_location X-coordinate location
 * @return Boundary layer properties
 */
BoundaryLayerProperties Calculate_Boundary_Layer_Properties(double x_location)
{
    BoundaryLayerProperties props;

    // Find velocity profile at x_location
    vector<pair<double, double>> velocity_profile;

    for (int i = 0; i < Total_Cells; i++)
    {
        if (abs(Cells[i].Cell_Center[0] - x_location) < 0.01)
        {
            double y = Cells[i].Cell_Center[1];
            double rho = Cells[i].Conservative_Variables[0];
            double u = Cells[i].Conservative_Variables[1] / rho;
            velocity_profile.push_back(make_pair(y, u));
        }
    }

    // Sort by y-coordinate
    sort(velocity_profile.begin(), velocity_profile.end());

    if (velocity_profile.empty())
    {
        return props; // Return default values
    }

    double U_inf = velocity_profile.back().second; // Free stream velocity

    // Calculate boundary layer thickness (99% of free stream velocity)
    props.delta_99 = 0.0;
    for (const auto &point : velocity_profile)
    {
        if (point.second >= 0.99 * U_inf)
        {
            props.delta_99 = point.first;
            break;
        }
    }

    // Calculate displacement and momentum thickness using trapezoidal integration
    props.delta_star = 0.0;
    props.theta = 0.0;

    for (size_t i = 1; i < velocity_profile.size(); i++)
    {
        double y1 = velocity_profile[i - 1].first;
        double y2 = velocity_profile[i].first;
        double u1 = velocity_profile[i - 1].second;
        double u2 = velocity_profile[i].second;

        double dy = y2 - y1;

        // Displacement thickness integral: ∫(1 - u/U_inf)dy
        double integrand_delta_star_1 = 1.0 - u1 / U_inf;
        double integrand_delta_star_2 = 1.0 - u2 / U_inf;
        props.delta_star += 0.5 * (integrand_delta_star_1 + integrand_delta_star_2) * dy;

        // Momentum thickness integral: ∫(u/U_inf)(1 - u/U_inf)dy
        double integrand_theta_1 = (u1 / U_inf) * (1.0 - u1 / U_inf);
        double integrand_theta_2 = (u2 / U_inf) * (1.0 - u2 / U_inf);
        props.theta += 0.5 * (integrand_theta_1 + integrand_theta_2) * dy;

        if (y2 > props.delta_99)
            break;
    }

    // Calculate ratios and shape factor
    props.delta_star_ratio = props.delta_star / characteristic_length;
    props.theta_ratio = props.theta / characteristic_length;
    props.shape_factor = (props.theta > 0) ? props.delta_star / props.theta : 0.0;

    return props;
}

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Print validation test parameters
 * @param params Validation parameters
 * @param model Turbulence model
 */
void Print_Validation_Parameters(const FlatPlateValidationParams &params, TurbulenceModel model)
{
    cout << "Test Parameters:" << endl;
    cout << "  Plate length: " << params.plate_length << " m" << endl;
    cout << "  Reynolds number: " << params.reynolds_number << endl;
    cout << "  Mach number: " << params.mach_number << endl;
    cout << "  Turbulence intensity: " << params.turbulence_intensity * 100.0 << "%" << endl;

    string model_name;
    switch (model)
    {
    case TurbulenceModel::K_EPSILON:
        model_name = "K-epsilon";
        break;
    case TurbulenceModel::K_OMEGA_WILCOX:
        model_name = "K-omega (Wilcox)";
        break;
    case TurbulenceModel::K_OMEGA_SST:
        model_name = "K-omega SST";
        break;
    default:
        model_name = "Unknown";
        break;
    }
    cout << "  Turbulence model: " << model_name << endl;
    cout << "  Grid size: " << params.nx << " x " << params.ny << endl;
    cout << endl;
}

/**
 * @brief Write detailed validation results to file
 * @param params Validation parameters
 * @param model Turbulence model
 */
void Write_Flat_Plate_Results(const FlatPlateValidationParams &params, TurbulenceModel model)
{
    string model_suffix = (model == TurbulenceModel::K_EPSILON) ? "kepsilon" : "komega";

    // Write velocity profiles
    string profile_file = "flat_plate_velocity_profile_" + model_suffix + ".dat";
    Write_Velocity_Profiles(profile_file, params);

    // Write skin friction distribution
    string cf_file = "flat_plate_skin_friction_" + model_suffix + ".dat";
    Write_Skin_Friction_Distribution(cf_file, params);

    // Write turbulence variables
    string turb_file = "flat_plate_turbulence_" + model_suffix + ".dat";
    Write_Turbulence_Variables(turb_file);

    cout << "Validation results written to files." << endl;
}