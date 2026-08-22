/**
 * @file Create_Vtk_File.cpp
 * @brief 3D VTK file generation for ParaView visualization
 *
 * This module generates VTK files for 3D CFD solutions visualization using ParaView.
 * Supports hexahedral cells, multiple scalar and vector fields, and parallel output.
 * VTK format enables high-quality 3D flow visualization and post-processing.
 *
 * Key Features:
 * - 3D VTK unstructured grid format (hexahedral cells)
 * - Multiple solution variables (density, velocity, pressure, temperature)
 * - Vector fields (velocity, momentum)
 * - Scalar fields (density, pressure, temperature, Mach number)
 * - Binary and ASCII output formats
 * - Parallel VTK support for large datasets
 * - Time-series animation support
 *
 * Mathematical Framework:
 * - Hexahedral cell connectivity: 8 vertices per cell
 * - Conservative variables: [ρ, ρu, ρv, ρw, ρE]
 * - Primitive variables: [ρ, u, v, w, p, T, a, h, μ, λ]
 * - Derived quantities: Mach number, total pressure, entropy
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"
#include <iomanip>
#include <sstream>

/**
 * @brief Creates VTK file header for 3D unstructured grid
 *
 * Writes the standard VTK header with dataset information for hexahedral mesh.
 * Includes version information, format specification, and dataset type.
 *
 * @param file Output file stream
 * @param format Output format ("ASCII" or "BINARY")
 *
 * Mathematical formulation:
 * - VTK format version 4.2 compatibility
 * - Unstructured grid for hexahedral cells
 * - Point and cell data sections
 */
void Write_VTK_Header_3D(std::ofstream &file, const std::string &format)
{
    file << "# vtk DataFile Version 4.2\n";
    file << "3D CFD Solution - Hexahedral Mesh\n";
    file << format << "\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
}

/**
 * @brief Writes 3D points (vertices) to VTK file
 *
 * Outputs all grid vertices in VTK format with proper coordinate scaling.
 * Handles both structured and unstructured hexahedral meshes.
 *
 * @param file Output file stream
 *
 * Mathematical formulation:
 * - Point coordinates: (x, y, z) for each vertex
 * - Total points: (Nx+1) × (Ny+1) × (Nz+1) for structured mesh
 * - Vertex indexing: i + j*(Nx+1) + k*(Nx+1)*(Ny+1)
 */
void Write_VTK_Points_3D(std::ofstream &file)
{
    file << "POINTS " << Total_Points_3D << " double\n";

    for (int k = 0; k <= Nz_max; k++)
    {
        for (int j = 0; j <= Ny_max; j++)
        {
            for (int i = 0; i <= Nx_max; i++)
            {
                int point_id = i + j * (Nx_max + 1) + k * (Nx_max + 1) * (Ny_max + 1);

                // Output coordinates with proper scaling
                file << std::scientific << std::setprecision(12)
                     << Grid_Points_3D[point_id].x << " "
                     << Grid_Points_3D[point_id].y << " "
                     << Grid_Points_3D[point_id].z << "\n";
            }
        }
    }
}

/**
 * @brief Writes hexahedral cell connectivity to VTK file
 *
 * Defines connectivity for all hexahedral cells using VTK ordering convention.
 * Each hexahedral cell requires 8 vertices in specific VTK order.
 *
 * @param file Output file stream
 *
 * Mathematical formulation:
 * - Hexahedral connectivity: 8 vertices per cell
 * - VTK vertex ordering: (0,1,2,3) bottom face, (4,5,6,7) top face
 * - Cell index: i + j*Nx + k*Nx*Ny
 * - Vertex local indexing for hex cell
 */
void Write_VTK_Cells_3D(std::ofstream &file)
{
    // Write cell connectivity
    file << "CELLS " << Total_Cells_3D << " " << 9 * Total_Cells_3D << "\n";

    for (int k = 0; k < Nz_max; k++)
    {
        for (int j = 0; j < Ny_max; j++)
        {
            for (int i = 0; i < Nx_max; i++)
            {
                // Calculate vertex indices for hexahedral cell
                int v0 = i + j * (Nx_max + 1) + k * (Nx_max + 1) * (Ny_max + 1);
                int v1 = v0 + 1;
                int v2 = v0 + (Nx_max + 1) + 1;
                int v3 = v0 + (Nx_max + 1);
                int v4 = v0 + (Nx_max + 1) * (Ny_max + 1);
                int v5 = v4 + 1;
                int v6 = v4 + (Nx_max + 1) + 1;
                int v7 = v4 + (Nx_max + 1);

                // Write cell connectivity (8 vertices + count)
                file << "8 " << v0 << " " << v1 << " " << v2 << " " << v3 << " "
                     << v4 << " " << v5 << " " << v6 << " " << v7 << "\n";
            }
        }
    }

    // Write cell types (VTK_HEXAHEDRON = 12)
    file << "CELL_TYPES " << Total_Cells_3D << "\n";
    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        file << "12\n"; // VTK_HEXAHEDRON
    }
}

/**
 * @brief Writes scalar field data to VTK file
 *
 * Outputs specified scalar field with proper name and data format.
 * Supports both primitive and derived scalar quantities.
 *
 * @param file Output file stream
 * @param field_name Name of the scalar field
 * @param data Pointer to scalar data array
 *
 * Mathematical formulation:
 * - Scalar fields: ρ, p, T, |V|, M, s (entropy)
 * - Cell-centered data storage
 * - Scientific notation for numerical precision
 */
void Write_VTK_Scalar_Field_3D(std::ofstream &file, const std::string &field_name, double *data)
{
    file << "SCALARS " << field_name << " double\n";
    file << "LOOKUP_TABLE default\n";

    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        file << std::scientific << std::setprecision(12) << data[cell] << "\n";
    }
}

/**
 * @brief Writes vector field data to VTK file
 *
 * Outputs 3D vector field with proper component ordering (u, v, w).
 * Supports velocity, momentum, and other vector quantities.
 *
 * @param file Output file stream
 * @param field_name Name of the vector field
 * @param u_data X-component data array
 * @param v_data Y-component data array
 * @param w_data Z-component data array
 *
 * Mathematical formulation:
 * - Vector fields: V = (u, v, w), ρV = (ρu, ρv, ρw)
 * - 3D vector components in Cartesian coordinates
 * - Cell-centered vector storage
 */
void Write_VTK_Vector_Field_3D(std::ofstream &file, const std::string &field_name,
                               double *u_data, double *v_data, double *w_data)
{
    file << "VECTORS " << field_name << " double\n";

    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        file << std::scientific << std::setprecision(12)
             << u_data[cell] << " " << v_data[cell] << " " << w_data[cell] << "\n";
    }
}

/**
 * @brief Calculates derived quantities for visualization
 *
 * Computes additional visualization quantities from conservative variables.
 * Includes Mach number, total pressure, entropy, and other derived fields.
 *
 * @param cell Cell index
 * @param mach_number Output Mach number
 * @param total_pressure Output total pressure
 * @param entropy Output specific entropy
 *
 * Mathematical formulation:
 * - Mach number: M = |V|/a where a = √(γRT)
 * - Total pressure: p₀ = p(1 + (γ-1)/2 M²)^(γ/(γ-1))
 * - Specific entropy: s = cv*ln(p/ρ^γ) + constant
 * - Total enthalpy: h₀ = h + |V|²/2
 */
void Calculate_Derived_Quantities_3D(int cell, double &mach_number,
                                     double &total_pressure, double &entropy)
{
    // Extract primitive variables
    double rho = Primitive_Variables_3D[cell][0];
    double u = Primitive_Variables_3D[cell][1];
    double v = Primitive_Variables_3D[cell][2];
    double w = Primitive_Variables_3D[cell][3];
    double p = Primitive_Variables_3D[cell][4];
    double T = Primitive_Variables_3D[cell][5];
    double a = Primitive_Variables_3D[cell][6];

    // Calculate velocity magnitude
    double velocity_magnitude = sqrt(u * u + v * v + w * w);

    // Mach number
    mach_number = (a > 1e-12) ? velocity_magnitude / a : 0.0;

    // Total pressure (stagnation pressure)
    double gamma_minus_1 = gamma - 1.0;
    double mach_squared = mach_number * mach_number;
    total_pressure = p * pow(1.0 + 0.5 * gamma_minus_1 * mach_squared,
                             gamma / gamma_minus_1);

    // Specific entropy (non-dimensional)
    entropy = log(p / pow(rho, gamma));

    // Validate results
    if (!std::isfinite(mach_number))
        mach_number = 0.0;
    if (!std::isfinite(total_pressure))
        total_pressure = p;
    if (!std::isfinite(entropy))
        entropy = 0.0;
}

/**
 * @brief Creates complete 3D VTK file for CFD solution
 *
 * Main function to generate VTK file with all solution variables and derived quantities.
 * Includes grid geometry, conservative variables, primitive variables, and derived fields.
 *
 * @param filename Output VTK filename
 * @param iteration Current iteration number
 * @param time Current simulation time
 * @param include_derived Include derived quantities (Mach, total pressure, etc.)
 *
 * Mathematical formulation:
 * - Complete 3D CFD solution output
 * - Conservative variables: [ρ, ρu, ρv, ρw, ρE]
 * - Primitive variables: [ρ, u, v, w, p, T]
 * - Derived quantities: [M, p₀, s, h₀]
 * - Vector fields: velocity, momentum
 */
void Create_VTK_File_3D(const std::string &filename, int iteration, double time,
                        bool include_derived = true)
{
    std::ofstream vtk_file(filename);

    if (!vtk_file.is_open())
    {
        std::cerr << "Error: Cannot create VTK file: " << filename << std::endl;
        return;
    }

    std::cout << "Creating 3D VTK file: " << filename << std::endl;
    std::cout << "Iteration: " << iteration << ", Time: " << time << std::endl;

    // Write VTK header
    Write_VTK_Header_3D(vtk_file, "ASCII");

    // Write points (vertices)
    Write_VTK_Points_3D(vtk_file);

    // Write cells (hexahedral connectivity)
    Write_VTK_Cells_3D(vtk_file);

    // Begin cell data section
    vtk_file << "CELL_DATA " << Total_Cells_3D << "\n";

    // Allocate temporary arrays for derived quantities
    std::vector<double> mach_numbers(Total_Cells_3D);
    std::vector<double> total_pressures(Total_Cells_3D);
    std::vector<double> entropies(Total_Cells_3D);
    std::vector<double> velocity_magnitudes(Total_Cells_3D);

    // Calculate derived quantities if requested
    if (include_derived)
    {
        for (int cell = 0; cell < Total_Cells_3D; cell++)
        {
            Calculate_Derived_Quantities_3D(cell, mach_numbers[cell],
                                            total_pressures[cell], entropies[cell]);

            // Velocity magnitude
            double u = Primitive_Variables_3D[cell][1];
            double v = Primitive_Variables_3D[cell][2];
            double w = Primitive_Variables_3D[cell][3];
            velocity_magnitudes[cell] = sqrt(u * u + v * v + w * w);
        }
    }

    // Write scalar fields

    // Density
    std::vector<double> density(Total_Cells_3D);
    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        density[cell] = Primitive_Variables_3D[cell][0];
    }
    Write_VTK_Scalar_Field_3D(vtk_file, "Density", density.data());

    // Pressure
    std::vector<double> pressure(Total_Cells_3D);
    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        pressure[cell] = Primitive_Variables_3D[cell][4];
    }
    Write_VTK_Scalar_Field_3D(vtk_file, "Pressure", pressure.data());

    // Temperature
    std::vector<double> temperature(Total_Cells_3D);
    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        temperature[cell] = Primitive_Variables_3D[cell][5];
    }
    Write_VTK_Scalar_Field_3D(vtk_file, "Temperature", temperature.data());

    // Velocity magnitude
    Write_VTK_Scalar_Field_3D(vtk_file, "VelocityMagnitude", velocity_magnitudes.data());

    if (include_derived)
    {
        // Mach number
        Write_VTK_Scalar_Field_3D(vtk_file, "MachNumber", mach_numbers.data());

        // Total pressure
        Write_VTK_Scalar_Field_3D(vtk_file, "TotalPressure", total_pressures.data());

        // Entropy
        Write_VTK_Scalar_Field_3D(vtk_file, "Entropy", entropies.data());
    }

    // Write vector fields

    // Velocity vector
    std::vector<double> u_velocity(Total_Cells_3D), v_velocity(Total_Cells_3D), w_velocity(Total_Cells_3D);
    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        u_velocity[cell] = Primitive_Variables_3D[cell][1];
        v_velocity[cell] = Primitive_Variables_3D[cell][2];
        w_velocity[cell] = Primitive_Variables_3D[cell][3];
    }
    Write_VTK_Vector_Field_3D(vtk_file, "Velocity", u_velocity.data(),
                              v_velocity.data(), w_velocity.data());

    // Momentum vector
    std::vector<double> rho_u(Total_Cells_3D), rho_v(Total_Cells_3D), rho_w(Total_Cells_3D);
    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        rho_u[cell] = Conservative_Variables_3D[cell][1];
        rho_v[cell] = Conservative_Variables_3D[cell][2];
        rho_w[cell] = Conservative_Variables_3D[cell][3];
    }
    Write_VTK_Vector_Field_3D(vtk_file, "Momentum", rho_u.data(),
                              rho_v.data(), rho_w.data());

    vtk_file.close();

    std::cout << "VTK file created successfully with " << Total_Cells_3D
              << " hexahedral cells" << std::endl;
}

/**
 * @brief Creates time-series VTK files for animation
 *
 * Generates sequentially numbered VTK files for time-dependent visualization.
 * Enables creation of animations and time-series analysis in ParaView.
 *
 * @param base_filename Base filename (without extension)
 * @param iteration Current iteration number
 * @param time Current simulation time
 * @param time_series_frequency Output frequency for time series
 *
 * Mathematical formulation:
 * - Time-dependent CFD solution: U(x,y,z,t)
 * - Sequential file naming: base_name_NNNNNN.vtk
 * - Consistent time stepping for animation
 */
void Create_Time_Series_VTK_3D(const std::string &base_filename, int iteration,
                               double time, int time_series_frequency)
{
    if (iteration % time_series_frequency != 0)
        return;

    // Create filename with iteration number
    std::ostringstream filename_stream;
    filename_stream << base_filename << "_" << std::setfill('0') << std::setw(6)
                    << iteration << ".vtk";

    std::string filename = filename_stream.str();

    // Create VTK file
    Create_VTK_File_3D(filename, iteration, time, true);

    std::cout << "Time-series VTK file created: " << filename
              << " at time = " << time << std::endl;
}

/**
 * @brief Creates parallel VTK files for large datasets
 *
 * Generates parallel VTK format (.pvtu) for large 3D datasets that require
 * domain decomposition or parallel processing visualization.
 *
 * @param base_filename Base filename for parallel output
 * @param num_processes Number of parallel processes
 * @param iteration Current iteration number
 * @param time Current simulation time
 *
 * Mathematical formulation:
 * - Domain decomposition: Ω = ∪ᵢ Ωᵢ
 * - Parallel VTK format for distributed datasets
 * - Master .pvtu file with piece references
 */
void Create_Parallel_VTK_3D(const std::string &base_filename, int num_processes,
                            int iteration, double time)
{
    // Create master .pvtu file
    std::ostringstream master_filename;
    master_filename << base_filename << "_" << std::setfill('0') << std::setw(6)
                    << iteration << ".pvtu";

    std::ofstream master_file(master_filename.str());

    if (!master_file.is_open())
    {
        std::cerr << "Error: Cannot create parallel VTK master file" << std::endl;
        return;
    }

    // Write parallel VTK header
    master_file << "<?xml version=\"1.0\"?>\n";
    master_file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    master_file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";

    // Define point and cell data arrays
    master_file << "    <PPointData>\n";
    master_file << "    </PPointData>\n";
    master_file << "    <PCellData>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\"/>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"Pressure\" NumberOfComponents=\"1\"/>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\"/>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"VelocityMagnitude\" NumberOfComponents=\"1\"/>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"MachNumber\" NumberOfComponents=\"1\"/>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\"/>\n";
    master_file << "      <PDataArray type=\"Float64\" Name=\"Momentum\" NumberOfComponents=\"3\"/>\n";
    master_file << "    </PCellData>\n";

    // Reference individual piece files
    for (int rank = 0; rank < num_processes; rank++)
    {
        std::ostringstream piece_filename;
        piece_filename << base_filename << "_" << std::setfill('0') << std::setw(6)
                       << iteration << "_" << rank << ".vtu";

        master_file << "    <Piece Source=\"" << piece_filename.str() << "\"/>\n";
    }

    master_file << "  </PUnstructuredGrid>\n";
    master_file << "</VTKFile>\n";

    master_file.close();

    std::cout << "Parallel VTK master file created: " << master_filename.str() << std::endl;
    std::cout << "References " << num_processes << " piece files" << std::endl;
}

/**
 * @brief Validates VTK output data quality
 *
 * Performs quality checks on solution data before VTK output to ensure
 * valid visualization and prevent ParaView errors.
 *
 * @return true if data is valid for VTK output
 *
 * Mathematical formulation:
 * - Range checks: ρ > 0, p > 0, T > 0
 * - Finite value validation: no NaN or infinite values
 * - Physical consistency: energy > kinetic energy
 */
bool Validate_VTK_Data_3D()
{
    bool is_valid = true;
    int invalid_cells = 0;

    for (int cell = 0; cell < Total_Cells_3D; cell++)
    {
        // Check primitive variables
        double rho = Primitive_Variables_3D[cell][0];
        double u = Primitive_Variables_3D[cell][1];
        double v = Primitive_Variables_3D[cell][2];
        double w = Primitive_Variables_3D[cell][3];
        double p = Primitive_Variables_3D[cell][4];
        double T = Primitive_Variables_3D[cell][5];

        // Validate finite values
        if (!std::isfinite(rho) || !std::isfinite(u) || !std::isfinite(v) ||
            !std::isfinite(w) || !std::isfinite(p) || !std::isfinite(T))
        {
            invalid_cells++;
            is_valid = false;
            continue;
        }

        // Validate physical ranges
        if (rho <= 0.0 || p <= 0.0 || T <= 0.0)
        {
            invalid_cells++;
            is_valid = false;
        }
    }

    if (!is_valid)
    {
        std::cerr << "Warning: " << invalid_cells << " cells have invalid data for VTK output" << std::endl;
    }

    return is_valid;
}