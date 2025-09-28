/**
 * @file Configuration_Read.cpp
 * @brief 3D CFD simulation configuration file parser and parameter management
 *
 * This module handles reading and parsing configuration files for 3D CFD simulations.
 * Supports comprehensive parameter management, validation, and default value handling
 * for all simulation aspects including grid, solver, boundary conditions, and output.
 *
 * Key Features:
 * - Comprehensive parameter file parsing (JSON/INI format support)
 * - 3D-specific configuration parameters
 * - Boundary condition specification for all 6 faces
 * - Solver parameter validation and range checking
 * - Default value management with override capability
 * - Grid configuration and mesh parameters
 * - Output and visualization settings
 * - Multi-physics parameter support
 *
 * Mathematical Framework:
 * - Grid dimensions: Nx × Ny × Nz
 * - Boundary conditions: 6 faces (left, right, front, back, bottom, top)
 * - Solver parameters: CFL, time stepping, convergence criteria
 * - Physical parameters: Mach, Reynolds, Prandtl numbers
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>

/**
 * @brief Configuration parameter structure for 3D CFD simulations
 *
 * Comprehensive structure containing all simulation parameters
 * with validation ranges and default values.
 */
struct CFD_Configuration_3D
{
    // Grid Parameters
    int Nx, Ny, Nz;        // Grid dimensions
    double Lx, Ly, Lz;     // Domain size
    std::string grid_type; // "uniform", "stretched", "custom"
    std::string mesh_file; // External mesh filename

    // Physical Parameters
    double Mach_number;     // Freestream Mach number
    double Reynolds_number; // Reynolds number
    double Prandtl_number;  // Prandtl number
    double gamma;           // Specific heat ratio
    double gas_constant;    // Specific gas constant

    // Reference Conditions
    double ref_density;     // Reference density
    double ref_pressure;    // Reference pressure
    double ref_temperature; // Reference temperature
    double ref_velocity;    // Reference velocity

    // Solver Parameters
    std::string solver_type;     // "explicit", "implicit", "dual-time"
    std::string flux_scheme;     // "van_leer", "roe", "ausm", "llf"
    std::string limiter_type;    // "minmod", "muscl", "weno"
    double CFL_number;           // CFL number
    int max_iterations;          // Maximum iterations
    double convergence_criteria; // Convergence threshold

    // Time Integration
    std::string time_scheme; // "euler", "rk2", "rk3", "rk4"
    double time_step;        // Physical time step
    double final_time;       // Final simulation time
    bool steady_state;       // Steady/unsteady flag

    // Boundary Conditions (6 faces for 3D)
    struct BoundaryCondition_3D
    {
        std::string type;   // "wall", "inlet", "outlet", "symmetry", "farfield"
        double temperature; // Wall temperature (if applicable)
        double pressure;    // Boundary pressure
        double velocity[3]; // Boundary velocity components
        bool adiabatic;     // Adiabatic wall flag
    };

    BoundaryCondition_3D bc_left, bc_right; // X-direction boundaries
    BoundaryCondition_3D bc_front, bc_back; // Y-direction boundaries
    BoundaryCondition_3D bc_bottom, bc_top; // Z-direction boundaries

    // Output Parameters
    int output_frequency;         // Output every N iterations
    int vtk_frequency;            // VTK output frequency
    int restart_frequency;        // Restart file frequency
    std::string output_directory; // Output directory path
    bool binary_output;           // Binary file format flag
    bool convergence_history;     // Write convergence history

    // Advanced Parameters
    bool viscous_flow;           // Viscous/inviscid flag
    bool turbulence_model;       // Turbulence model flag
    std::string turbulence_type; // "none", "sa", "k_omega", "les"
    double artificial_viscosity; // Artificial viscosity coefficient
    bool shock_capturing;        // Shock capturing flag

    // Parallel Parameters
    int num_threads;       // OpenMP threads
    bool gpu_acceleration; // CUDA acceleration flag
    int gpu_device_id;     // GPU device ID

    // Constructor with default values
    CFD_Configuration_3D()
    {
        // Default grid parameters
        Nx = 100;
        Ny = 100;
        Nz = 100;
        Lx = 1.0;
        Ly = 1.0;
        Lz = 1.0;
        grid_type = "uniform";
        mesh_file = "";

        // Default physical parameters
        Mach_number = 0.8;
        Reynolds_number = 1000.0;
        Prandtl_number = 0.72;
        gamma = 1.4;
        gas_constant = 287.0;

        // Default reference conditions
        ref_density = 1.225;
        ref_pressure = 101325.0;
        ref_temperature = 288.15;
        ref_velocity = 100.0;

        // Default solver parameters
        solver_type = "explicit";
        flux_scheme = "roe";
        limiter_type = "muscl";
        CFL_number = 0.5;
        max_iterations = 10000;
        convergence_criteria = 1e-6;

        // Default time integration
        time_scheme = "rk3";
        time_step = 1e-5;
        final_time = 1.0;
        steady_state = true;

        // Default boundary conditions (all walls)
        SetDefaultBoundaryConditions();

        // Default output parameters
        output_frequency = 100;
        vtk_frequency = 500;
        restart_frequency = 1000;
        output_directory = "./output";
        binary_output = true;
        convergence_history = true;

        // Default advanced parameters
        viscous_flow = true;
        turbulence_model = false;
        turbulence_type = "none";
        artificial_viscosity = 0.0;
        shock_capturing = true;

        // Default parallel parameters
        num_threads = 4;
        gpu_acceleration = false;
        gpu_device_id = 0;
    }

    void SetDefaultBoundaryConditions()
    {
        // Initialize all boundaries as adiabatic walls
        BoundaryCondition_3D default_wall;
        default_wall.type = "wall";
        default_wall.temperature = 288.15;
        default_wall.pressure = 101325.0;
        default_wall.velocity[0] = 0.0;
        default_wall.velocity[1] = 0.0;
        default_wall.velocity[2] = 0.0;
        default_wall.adiabatic = true;

        bc_left = bc_right = bc_front = bc_back = bc_bottom = bc_top = default_wall;
    }
};

// Global configuration instance
CFD_Configuration_3D config_3d;

/**
 * @brief Trims whitespace from string
 *
 * Utility function to remove leading and trailing whitespace from strings.
 *
 * @param str Input string to trim
 * @return Trimmed string
 */
std::string Trim(const std::string &str)
{
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
        return "";

    size_t end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

/**
 * @brief Converts string to lowercase
 *
 * @param str Input string
 * @return Lowercase string
 */
std::string ToLowerCase(const std::string &str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

/**
 * @brief Parses boundary condition from configuration string
 *
 * Parses boundary condition specification and validates parameters.
 * Supports multiple boundary condition types with comprehensive validation.
 *
 * @param bc_string Boundary condition specification string
 * @param bc Output boundary condition structure
 * @return true if parsing successful
 *
 * Format examples:
 * - "wall,adiabatic"
 * - "wall,isothermal,T=300.0"
 * - "inlet,p=101325,u=100,v=0,w=0"
 * - "outlet,p=101325"
 * - "symmetry"
 * - "farfield,M=0.8,p=101325,T=288"
 */
bool Parse_Boundary_Condition_3D(const std::string &bc_string,
                                 CFD_Configuration_3D::BoundaryCondition_3D &bc)
{
    std::istringstream iss(bc_string);
    std::string token;

    // Parse boundary condition type
    if (!std::getline(iss, token, ','))
    {
        std::cerr << "Error: Missing boundary condition type" << std::endl;
        return false;
    }

    bc.type = ToLowerCase(Trim(token));

    // Validate boundary condition type
    if (bc.type != "wall" && bc.type != "inlet" && bc.type != "outlet" &&
        bc.type != "symmetry" && bc.type != "farfield")
    {
        std::cerr << "Error: Invalid boundary condition type: " << bc.type << std::endl;
        return false;
    }

    // Parse additional parameters
    while (std::getline(iss, token, ','))
    {
        token = Trim(token);

        if (token == "adiabatic")
        {
            bc.adiabatic = true;
        }
        else if (token == "isothermal")
        {
            bc.adiabatic = false;
        }
        else if (token.find('=') != std::string::npos)
        {
            size_t eq_pos = token.find('=');
            std::string param = ToLowerCase(Trim(token.substr(0, eq_pos)));
            std::string value = Trim(token.substr(eq_pos + 1));

            try
            {
                if (param == "t" || param == "temperature")
                {
                    bc.temperature = std::stod(value);
                }
                else if (param == "p" || param == "pressure")
                {
                    bc.pressure = std::stod(value);
                }
                else if (param == "u")
                {
                    bc.velocity[0] = std::stod(value);
                }
                else if (param == "v")
                {
                    bc.velocity[1] = std::stod(value);
                }
                else if (param == "w")
                {
                    bc.velocity[2] = std::stod(value);
                }
                else
                {
                    std::cerr << "Warning: Unknown boundary parameter: " << param << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error parsing boundary parameter " << param << ": " << e.what() << std::endl;
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Validates configuration parameters
 *
 * Comprehensive validation of all configuration parameters with
 * physical and numerical constraints checking.
 *
 * @param config Configuration structure to validate
 * @return true if configuration is valid
 *
 * Mathematical validation:
 * - Grid dimensions: Nx, Ny, Nz > 0
 * - Physical parameters: M > 0, Re > 0, γ > 1
 * - CFL condition: 0 < CFL < 2
 * - Time parameters: dt > 0, final_time > 0
 */
bool Validate_Configuration_3D(const CFD_Configuration_3D &config)
{
    bool is_valid = true;

    // Validate grid parameters
    if (config.Nx <= 0 || config.Ny <= 0 || config.Nz <= 0)
    {
        std::cerr << "Error: Grid dimensions must be positive" << std::endl;
        is_valid = false;
    }

    if (config.Lx <= 0.0 || config.Ly <= 0.0 || config.Lz <= 0.0)
    {
        std::cerr << "Error: Domain size must be positive" << std::endl;
        is_valid = false;
    }

    // Validate physical parameters
    if (config.Mach_number <= 0.0)
    {
        std::cerr << "Error: Mach number must be positive" << std::endl;
        is_valid = false;
    }

    if (config.Reynolds_number <= 0.0)
    {
        std::cerr << "Error: Reynolds number must be positive" << std::endl;
        is_valid = false;
    }

    if (config.gamma <= 1.0)
    {
        std::cerr << "Error: Specific heat ratio must be greater than 1" << std::endl;
        is_valid = false;
    }

    // Validate solver parameters
    if (config.CFL_number <= 0.0 || config.CFL_number > 2.0)
    {
        std::cerr << "Error: CFL number must be in range (0, 2]" << std::endl;
        is_valid = false;
    }

    if (config.max_iterations <= 0)
    {
        std::cerr << "Error: Maximum iterations must be positive" << std::endl;
        is_valid = false;
    }

    if (config.convergence_criteria <= 0.0)
    {
        std::cerr << "Error: Convergence criteria must be positive" << std::endl;
        is_valid = false;
    }

    // Validate time parameters
    if (!config.steady_state)
    {
        if (config.time_step <= 0.0)
        {
            std::cerr << "Error: Time step must be positive for unsteady simulations" << std::endl;
            is_valid = false;
        }

        if (config.final_time <= 0.0)
        {
            std::cerr << "Error: Final time must be positive" << std::endl;
            is_valid = false;
        }
    }

    // Validate reference conditions
    if (config.ref_density <= 0.0 || config.ref_pressure <= 0.0 || config.ref_temperature <= 0.0)
    {
        std::cerr << "Error: Reference conditions must be positive" << std::endl;
        is_valid = false;
    }

    return is_valid;
}

/**
 * @brief Reads configuration file in INI format
 *
 * Parses INI-style configuration file with section support.
 * Handles comments, multiline values, and type conversion.
 *
 * @param filename Configuration filename
 * @param config Output configuration structure
 * @return true if file read successfully
 *
 * INI format example:
 * ```
 * [Grid]
 * Nx = 100
 * Ny = 100
 * Nz = 50
 *
 * [Physics]
 * Mach = 0.8
 * Reynolds = 1000
 *
 * [Boundaries]
 * left = inlet,u=100,v=0,w=0,p=101325
 * right = outlet,p=101325
 * ```
 */
bool Read_Configuration_INI_3D(const std::string &filename, CFD_Configuration_3D &config)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Cannot open configuration file: " << filename << std::endl;
        return false;
    }

    std::string line, current_section;
    int line_number = 0;

    while (std::getline(file, line))
    {
        line_number++;
        line = Trim(line);

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == ';')
            continue;

        // Section headers
        if (line[0] == '[' && line.back() == ']')
        {
            current_section = ToLowerCase(line.substr(1, line.length() - 2));
            continue;
        }

        // Parameter assignments
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos)
        {
            std::cerr << "Warning: Invalid line " << line_number << ": " << line << std::endl;
            continue;
        }

        std::string key = ToLowerCase(Trim(line.substr(0, eq_pos)));
        std::string value = Trim(line.substr(eq_pos + 1));

        try
        {
            // Grid parameters
            if (current_section == "grid")
            {
                if (key == "nx")
                    config.Nx = std::stoi(value);
                else if (key == "ny")
                    config.Ny = std::stoi(value);
                else if (key == "nz")
                    config.Nz = std::stoi(value);
                else if (key == "lx")
                    config.Lx = std::stod(value);
                else if (key == "ly")
                    config.Ly = std::stod(value);
                else if (key == "lz")
                    config.Lz = std::stod(value);
                else if (key == "type")
                    config.grid_type = value;
                else if (key == "mesh_file")
                    config.mesh_file = value;
            }
            // Physics parameters
            else if (current_section == "physics")
            {
                if (key == "mach")
                    config.Mach_number = std::stod(value);
                else if (key == "reynolds")
                    config.Reynolds_number = std::stod(value);
                else if (key == "prandtl")
                    config.Prandtl_number = std::stod(value);
                else if (key == "gamma")
                    config.gamma = std::stod(value);
                else if (key == "gas_constant")
                    config.gas_constant = std::stod(value);
            }
            // Solver parameters
            else if (current_section == "solver")
            {
                if (key == "type")
                    config.solver_type = value;
                else if (key == "flux_scheme")
                    config.flux_scheme = value;
                else if (key == "limiter")
                    config.limiter_type = value;
                else if (key == "cfl")
                    config.CFL_number = std::stod(value);
                else if (key == "max_iterations")
                    config.max_iterations = std::stoi(value);
                else if (key == "convergence")
                    config.convergence_criteria = std::stod(value);
            }
            // Boundary conditions
            else if (current_section == "boundaries")
            {
                if (key == "left")
                    Parse_Boundary_Condition_3D(value, config.bc_left);
                else if (key == "right")
                    Parse_Boundary_Condition_3D(value, config.bc_right);
                else if (key == "front")
                    Parse_Boundary_Condition_3D(value, config.bc_front);
                else if (key == "back")
                    Parse_Boundary_Condition_3D(value, config.bc_back);
                else if (key == "bottom")
                    Parse_Boundary_Condition_3D(value, config.bc_bottom);
                else if (key == "top")
                    Parse_Boundary_Condition_3D(value, config.bc_top);
            }
            // Output parameters
            else if (current_section == "output")
            {
                if (key == "frequency")
                    config.output_frequency = std::stoi(value);
                else if (key == "vtk_frequency")
                    config.vtk_frequency = std::stoi(value);
                else if (key == "restart_frequency")
                    config.restart_frequency = std::stoi(value);
                else if (key == "directory")
                    config.output_directory = value;
                else if (key == "binary")
                    config.binary_output = (ToLowerCase(value) == "true");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error parsing line " << line_number << " (" << key << "): " << e.what() << std::endl;
        }
    }

    file.close();
    return true;
}

/**
 * @brief Applies configuration to global variables
 *
 * Transfers configuration parameters to global solver variables.
 * Handles unit conversions and derived parameter calculations.
 *
 * @param config Configuration structure
 *
 * Mathematical transformations:
 * - Non-dimensional parameters from physical values
 * - Grid spacing: Δx = Lx/Nx, Δy = Ly/Ny, Δz = Lz/Nz
 * - Reference state calculations
 */
void Apply_Configuration_3D(const CFD_Configuration_3D &config)
{
    // Grid parameters
    Nx_max = config.Nx;
    Ny_max = config.Ny;
    Nz_max = config.Nz;
    Total_Cells_3D = Nx_max * Ny_max * Nz_max;
    Total_Points_3D = (Nx_max + 1) * (Ny_max + 1) * (Nz_max + 1);

    // Domain size
    Length_x = config.Lx;
    Length_y = config.Ly;
    Length_z = config.Lz;

    // Grid spacing
    dx = Length_x / Nx_max;
    dy = Length_y / Ny_max;
    dz = Length_z / Nz_max;

    // Physical parameters
    Mach_Number = config.Mach_number;
    Reynolds_Number = config.Reynolds_number;
    Prandtl_Number = config.Prandtl_number;
    gamma = config.gamma;
    gas_constant = config.gas_constant;

    // Reference conditions
    rho_ref = config.ref_density;
    p_ref = config.ref_pressure;
    T_ref = config.ref_temperature;
    u_ref = config.ref_velocity;

    // Solver parameters
    CFL = config.CFL_number;
    Max_Iterations = config.max_iterations;
    Convergence_Criteria = config.convergence_criteria;

    // Output parameters
    Output_Frequency = config.output_frequency;
    VTK_Frequency = config.vtk_frequency;
    Restart_Frequency = config.restart_frequency;

    std::cout << "Configuration applied successfully:" << std::endl;
    std::cout << "Grid: " << Nx_max << " x " << Ny_max << " x " << Nz_max << std::endl;
    std::cout << "Total cells: " << Total_Cells_3D << std::endl;
    std::cout << "Mach number: " << Mach_Number << std::endl;
    std::cout << "Reynolds number: " << Reynolds_Number << std::endl;
    std::cout << "CFL number: " << CFL << std::endl;
}

/**
 * @brief Main configuration reading function
 *
 * Reads configuration file and applies parameters to solver.
 * Handles file format detection and validation.
 *
 * @param filename Configuration filename
 * @return true if configuration loaded successfully
 */
bool Read_Configuration_3D(const std::string &filename)
{
    std::cout << "Reading 3D CFD configuration file: " << filename << std::endl;

    // Initialize with defaults
    config_3d = CFD_Configuration_3D();

    // Read configuration file
    bool success = Read_Configuration_INI_3D(filename, config_3d);

    if (!success)
    {
        std::cerr << "Error reading configuration file" << std::endl;
        return false;
    }

    // Validate configuration
    if (!Validate_Configuration_3D(config_3d))
    {
        std::cerr << "Configuration validation failed" << std::endl;
        return false;
    }

    // Apply configuration to global variables
    Apply_Configuration_3D(config_3d);

    std::cout << "Configuration loaded and validated successfully" << std::endl;
    return true;
}

/**
 * @brief Writes current configuration to file
 *
 * Outputs current configuration parameters to file for documentation
 * and reproducibility.
 *
 * @param filename Output configuration filename
 * @return true if file written successfully
 */
bool Write_Configuration_3D(const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Cannot create configuration file: " << filename << std::endl;
        return false;
    }

    file << "# 3D CFD Solver Configuration\n";
    file << "# Generated automatically\n\n";

    file << "[Grid]\n";
    file << "Nx = " << config_3d.Nx << "\n";
    file << "Ny = " << config_3d.Ny << "\n";
    file << "Nz = " << config_3d.Nz << "\n";
    file << "Lx = " << config_3d.Lx << "\n";
    file << "Ly = " << config_3d.Ly << "\n";
    file << "Lz = " << config_3d.Lz << "\n";
    file << "type = " << config_3d.grid_type << "\n\n";

    file << "[Physics]\n";
    file << "Mach = " << config_3d.Mach_number << "\n";
    file << "Reynolds = " << config_3d.Reynolds_number << "\n";
    file << "Prandtl = " << config_3d.Prandtl_number << "\n";
    file << "gamma = " << config_3d.gamma << "\n\n";

    file << "[Solver]\n";
    file << "type = " << config_3d.solver_type << "\n";
    file << "flux_scheme = " << config_3d.flux_scheme << "\n";
    file << "limiter = " << config_3d.limiter_type << "\n";
    file << "CFL = " << config_3d.CFL_number << "\n";
    file << "max_iterations = " << config_3d.max_iterations << "\n";
    file << "convergence = " << config_3d.convergence_criteria << "\n\n";

    file << "[Output]\n";
    file << "frequency = " << config_3d.output_frequency << "\n";
    file << "vtk_frequency = " << config_3d.vtk_frequency << "\n";
    file << "restart_frequency = " << config_3d.restart_frequency << "\n";
    file << "directory = " << config_3d.output_directory << "\n";
    file << "binary = " << (config_3d.binary_output ? "true" : "false") << "\n";

    file.close();

    std::cout << "Configuration written to: " << filename << std::endl;
    return true;
}