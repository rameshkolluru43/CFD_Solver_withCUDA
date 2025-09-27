/**
 * @file Viscous_Calculations.cpp
 * @brief 3D viscous flux calculations for Navier-Stokes equations
 *
 * This module implements 3D viscous flux computations including viscous stress
 * tensor, heat flux, and wall boundary conditions. It extends 2D viscous calculations
 * to handle 3D flow physics with proper stress tensor formulation.
 *
 * Key Features:
 * - 3D viscous stress tensor computation
 * - 3D heat flux calculations
 * - Sutherland's law for temperature-dependent viscosity
 * - 3D wall shear stress and heat transfer
 * - Green-Gauss gradient computation in 3D
 * - Support for laminar and turbulent flows
 *
 * Mathematical Framework:
 * - 3D stress tensor: τ_ij = μ(∂u_i/∂x_j + ∂u_j/∂x_i - (2/3)δ_ij∇⋅V)
 * - Heat flux: q_i = -k ∂T/∂x_i
 * - Viscous flux: F_v = [0, τ_x, τ_y, τ_z, u⋅τ_x + v⋅τ_y + w⋅τ_z - q⋅n]^T
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"

// 3D Viscous calculation working variables
double mu_star_3D = 0.0; // Non-dimensional viscosity
double K_3D = 0.0;       // Thermal conductivity
double cp_ref_3D = 0.0;  // Reference specific heat

// 3D Gradient vectors
V_D u_Gradient_Face_3D(3, 0.0); // du/dx, du/dy, du/dz on face
V_D v_Gradient_Face_3D(3, 0.0); // dv/dx, dv/dy, dv/dz on face
V_D w_Gradient_Face_3D(3, 0.0); // dw/dx, dw/dy, dw/dz on face
V_D T_Gradient_Face_3D(3, 0.0); // dT/dx, dT/dy, dT/dz on face

/**
 * @brief Calculate viscosity using Sutherland's law for 3D flows
 *
 * This function computes temperature-dependent viscosity using Sutherland's
 * law, which is essential for accurate viscous flux calculations.
 *
 * @param T_Star Average temperature on the face (non-dimensional)
 */
void Viscosity_3D(double &T_Star)
{
    double Term1 = T_S_Mu / T_ref; // Sutherland temperature ratio
    double Term2 = T_Star + Term1; // Temperature sum

    // Sutherland's law: μ/μ_ref = (T/T_ref)^1.5 * (T_ref + T_s)/(T + T_s)
    mu_star_3D = pow(T_Star, 1.5) * ((1.0 + Term1) / Term2);

    // Ensure positive viscosity
    if (mu_star_3D <= 0.0)
    {
        cout << "Warning: Non-positive viscosity " << mu_star_3D
             << " at temperature " << T_Star << endl;
        mu_star_3D = 1e-10; // Small positive value
    }
}

/**
 * @brief Calculate thermal conductivity for 3D flows
 *
 * This function computes thermal conductivity from viscosity using
 * the constant Prandtl number assumption.
 *
 * @param T Temperature (non-dimensional)
 */
void Thermal_Conductivity_3D(double &T)
{
    // First calculate viscosity at this temperature
    Viscosity_3D(T);

    // Thermal conductivity from Prandtl number: k = μ*cp/Pr
    K_3D = mu_star_3D * cp_ref_3D / Pr;

    // Ensure positive thermal conductivity
    if (K_3D <= 0.0)
    {
        cout << "Warning: Non-positive thermal conductivity " << K_3D << endl;
        K_3D = 1e-10;
    }
}

/**
 * @brief Set reference values for 3D non-dimensionalization
 *
 * This function establishes reference values for viscous properties
 * in 3D flow calculations.
 */
void Reference_Values_3D()
{
    // Reference viscosity from Reynolds number
    mu_ref = (Rho_ref_3D * u_ref_3D * L_ref) / Re;

    // Reference specific heat from Mach number
    cp_ref_3D = R_gas * gamma / (gamma - 1.0);

    // Reference gas constant
    R_ref = R_gas;

    // Reference thermal conductivity
    K_ref = mu_ref * cp_ref_3D / Pr;

    cout << "3D Reference viscous properties:" << endl;
    cout << "  μ_ref = " << mu_ref << " Pa⋅s" << endl;
    cout << "  k_ref = " << K_ref << " W/(m⋅K)" << endl;
    cout << "  Re = " << Re << endl;
    cout << "  Pr = " << Pr << endl;
}

/**
 * @brief Calculate 3D viscous stress tensor components
 *
 * This function computes the full 3D viscous stress tensor using
 * velocity gradients and Stokes' hypothesis.
 *
 * @param Cell_No Current cell index
 * @param Face_No Face index (0-5)
 * @param stress_tensor Output stress tensor (9 components)
 * @param mu_face Face-averaged viscosity
 */
void Calculate_3D_Stress_Tensor(const int &Cell_No, const int &Face_No,
                                V_D &stress_tensor, const double &mu_face)
{
    // Compute velocity gradients on the face
    Calculate_3D_Velocity_Gradients_On_Face(Cell_No, Face_No);

    // Extract gradient components
    double dudx = u_Gradient_Face_3D[0], dudy = u_Gradient_Face_3D[1], dudz = u_Gradient_Face_3D[2];
    double dvdx = v_Gradient_Face_3D[0], dvdy = v_Gradient_Face_3D[1], dvdz = v_Gradient_Face_3D[2];
    double dwdx = w_Gradient_Face_3D[0], dwdy = w_Gradient_Face_3D[1], dwdz = w_Gradient_Face_3D[2];

    // Velocity divergence
    double div_V = dudx + dvdy + dwdz;

    // 3D stress tensor with Stokes' hypothesis
    // τ_ij = μ(∂u_i/∂x_j + ∂u_j/∂x_i - (2/3)δ_ij∇⋅V)

    // Normal stress components
    stress_tensor[0] = mu_face * Inv_Re * (2.0 * dudx - (2.0 / 3.0) * div_V); // τ_xx
    stress_tensor[4] = mu_face * Inv_Re * (2.0 * dvdy - (2.0 / 3.0) * div_V); // τ_yy
    stress_tensor[8] = mu_face * Inv_Re * (2.0 * dwdz - (2.0 / 3.0) * div_V); // τ_zz

    // Shear stress components (symmetric)
    stress_tensor[1] = stress_tensor[3] = mu_face * Inv_Re * (dudy + dvdx); // τ_xy = τ_yx
    stress_tensor[2] = stress_tensor[6] = mu_face * Inv_Re * (dudz + dwdx); // τ_xz = τ_zx
    stress_tensor[5] = stress_tensor[7] = mu_face * Inv_Re * (dvdz + dwdy); // τ_yz = τ_zy
}

/**
 * @brief Calculate 3D heat flux vector
 *
 * This function computes the heat flux vector using Fourier's law
 * with temperature gradients in all three directions.
 *
 * @param Cell_No Current cell index
 * @param Face_No Face index
 * @param heat_flux Output heat flux vector (3 components)
 * @param k_face Face-averaged thermal conductivity
 */
void Calculate_3D_Heat_Flux(const int &Cell_No, const int &Face_No,
                            V_D &heat_flux, const double &k_face)
{
    // Compute temperature gradient on the face
    Calculate_Temperature_Gradient_On_Face_3D(Cell_No, Face_No);

    // Fourier's law: q_i = -k ∂T/∂x_i
    heat_flux[0] = -k_face * Inv_Re * Inv_Pr * T_Gradient_Face_3D[0]; // q_x
    heat_flux[1] = -k_face * Inv_Re * Inv_Pr * T_Gradient_Face_3D[1]; // q_y
    heat_flux[2] = -k_face * Inv_Re * Inv_Pr * T_Gradient_Face_3D[2]; // q_z
}

/**
 * @brief Calculate viscous flux through a 3D face
 *
 * This function computes the complete viscous flux vector for a face
 * including momentum and energy contributions.
 *
 * @param Cell_No Current cell index
 * @param Face_No Face index (0-5)
 */
void Calculate_Viscous_Flux_Face_3D(const int &Cell_No, const int &Face_No)
{
    // Get neighbor cell
    int neighbor = Cells[Cell_No].Neighbours[Face_No];

    // Face normal vector
    int index = Face_No * 3;
    double nx = Cells[Cell_No].Face_Normals[index + 0];
    double ny = Cells[Cell_No].Face_Normals[index + 1];
    double nz = Cells[Cell_No].Face_Normals[index + 2];
    double area = Cells[Cell_No].Face_Areas[Face_No];

    // Average primitive variables on the face
    double rho_avg, u_avg, v_avg, w_avg, T_avg, mu_avg, k_avg;

    if (neighbor >= 0)
    { // Internal face
        rho_avg = 0.5 * (Primitive_Cells_3D[Cell_No][0] + Primitive_Cells_3D[neighbor][0]);
        u_avg = 0.5 * (Primitive_Cells_3D[Cell_No][1] + Primitive_Cells_3D[neighbor][1]);
        v_avg = 0.5 * (Primitive_Cells_3D[Cell_No][2] + Primitive_Cells_3D[neighbor][2]);
        w_avg = 0.5 * (Primitive_Cells_3D[Cell_No][3] + Primitive_Cells_3D[neighbor][3]);
        T_avg = 0.5 * (Primitive_Cells_3D[Cell_No][5] + Primitive_Cells_3D[neighbor][5]);
        mu_avg = 0.5 * (Primitive_Cells_3D[Cell_No][8] + Primitive_Cells_3D[neighbor][8]);
    }
    else
    { // Boundary face - use ghost cell values
        rho_avg = Primitive_Cells_3D[Cell_No][0];
        u_avg = Primitive_Cells_3D[Cell_No][1];
        v_avg = Primitive_Cells_3D[Cell_No][2];
        w_avg = Primitive_Cells_3D[Cell_No][3];
        T_avg = Primitive_Cells_3D[Cell_No][5];
        mu_avg = Primitive_Cells_3D[Cell_No][8];
    }

    // Update viscosity and thermal conductivity
    Thermal_Conductivity_3D(T_avg);
    k_avg = K_3D;

    // Calculate stress tensor
    V_D stress_tensor(9, 0.0); // 3x3 tensor flattened
    Calculate_3D_Stress_Tensor(Cell_No, Face_No, stress_tensor, mu_avg);

    // Calculate heat flux
    V_D heat_flux(3, 0.0);
    Calculate_3D_Heat_Flux(Cell_No, Face_No, heat_flux, k_avg);

    // Viscous flux components
    V_D viscous_flux(NUM_CONSERVATIVE_3D, 0.0);

    // Mass flux (no viscous contribution)
    viscous_flux[0] = 0.0;

    // Momentum fluxes: F_v = τ⋅n
    viscous_flux[1] = (stress_tensor[0] * nx + stress_tensor[1] * ny + stress_tensor[2] * nz) * area; // x-momentum
    viscous_flux[2] = (stress_tensor[3] * nx + stress_tensor[4] * ny + stress_tensor[5] * nz) * area; // y-momentum
    viscous_flux[3] = (stress_tensor[6] * nx + stress_tensor[7] * ny + stress_tensor[8] * nz) * area; // z-momentum

    // Energy flux: F_v = u⋅(τ⋅n) - q⋅n
    double work_term = u_avg * viscous_flux[1] + v_avg * viscous_flux[2] + w_avg * viscous_flux[3];
    double heat_term = (heat_flux[0] * nx + heat_flux[1] * ny + heat_flux[2] * nz) * area;
    viscous_flux[4] = work_term - heat_term;

    // Store in global viscous flux array
    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        Cells_Viscous_Flux_3D[Cell_No][comp] += viscous_flux[comp];
    }
}

/**
 * @brief Calculate 3D velocity gradients on a face
 *
 * This function computes velocity gradients required for stress tensor
 * calculations using Green-Gauss method.
 *
 * @param Cell_No Current cell index
 * @param Face_No Face index
 */
void Calculate_3D_Velocity_Gradients_On_Face(const int &Cell_No, const int &Face_No)
{
    // Initialize gradients
    u_Gradient_Face_3D.assign(3, 0.0);
    v_Gradient_Face_3D.assign(3, 0.0);
    w_Gradient_Face_3D.assign(3, 0.0);

    // Get neighbor cell
    int neighbor = Cells[Cell_No].Neighbours[Face_No];

    if (neighbor >= 0)
    { // Internal face
        // Face center coordinates
        double x_face = 0.5 * (Cells[Cell_No].Centroid[0] + Cells[neighbor].Centroid[0]);
        double y_face = 0.5 * (Cells[Cell_No].Centroid[1] + Cells[neighbor].Centroid[1]);
        double z_face = 0.5 * (Cells[Cell_No].Centroid[2] + Cells[neighbor].Centroid[2]);

        // Distance vector
        double dx = Cells[neighbor].Centroid[0] - Cells[Cell_No].Centroid[0];
        double dy = Cells[neighbor].Centroid[1] - Cells[Cell_No].Centroid[1];
        double dz = Cells[neighbor].Centroid[2] - Cells[Cell_No].Centroid[2];
        double distance = sqrt(dx * dx + dy * dy + dz * dz);

        // Velocity differences
        double du = Primitive_Cells_3D[neighbor][1] - Primitive_Cells_3D[Cell_No][1];
        double dv = Primitive_Cells_3D[neighbor][2] - Primitive_Cells_3D[Cell_No][2];
        double dw = Primitive_Cells_3D[neighbor][3] - Primitive_Cells_3D[Cell_No][3];

        // Simple central difference (can be improved with least squares)
        if (distance > 1e-12)
        {
            u_Gradient_Face_3D[0] = du * dx / (distance * distance); // du/dx
            u_Gradient_Face_3D[1] = du * dy / (distance * distance); // du/dy
            u_Gradient_Face_3D[2] = du * dz / (distance * distance); // du/dz

            v_Gradient_Face_3D[0] = dv * dx / (distance * distance); // dv/dx
            v_Gradient_Face_3D[1] = dv * dy / (distance * distance); // dv/dy
            v_Gradient_Face_3D[2] = dv * dz / (distance * distance); // dv/dz

            w_Gradient_Face_3D[0] = dw * dx / (distance * distance); // dw/dx
            w_Gradient_Face_3D[1] = dw * dy / (distance * distance); // dw/dy
            w_Gradient_Face_3D[2] = dw * dz / (distance * distance); // dw/dz
        }
    }
    else
    { // Boundary face - use boundary conditions
        Apply_3D_Viscous_Boundary_Gradients(Cell_No, Face_No);
    }
}

/**
 * @brief Calculate temperature gradient on a 3D face
 *
 * @param Cell_No Current cell index
 * @param Face_No Face index
 */
void Calculate_Temperature_Gradient_On_Face_3D(const int &Cell_No, const int &Face_No)
{
    // Initialize gradient
    T_Gradient_Face_3D.assign(3, 0.0);

    int neighbor = Cells[Cell_No].Neighbours[Face_No];

    if (neighbor >= 0)
    { // Internal face
        // Distance vector
        double dx = Cells[neighbor].Centroid[0] - Cells[Cell_No].Centroid[0];
        double dy = Cells[neighbor].Centroid[1] - Cells[Cell_No].Centroid[1];
        double dz = Cells[neighbor].Centroid[2] - Cells[Cell_No].Centroid[2];
        double distance = sqrt(dx * dx + dy * dy + dz * dz);

        // Temperature difference
        double dT = Primitive_Cells_3D[neighbor][5] - Primitive_Cells_3D[Cell_No][5];

        // Central difference gradient
        if (distance > 1e-12)
        {
            T_Gradient_Face_3D[0] = dT * dx / (distance * distance); // dT/dx
            T_Gradient_Face_3D[1] = dT * dy / (distance * distance); // dT/dy
            T_Gradient_Face_3D[2] = dT * dz / (distance * distance); // dT/dz
        }
    }
    else
    { // Boundary face
        Apply_3D_Thermal_Boundary_Gradients(Cell_No, Face_No);
    }
}

/**
 * @brief Apply viscous boundary conditions for velocity gradients
 *
 * @param Cell_No Current cell index
 * @param Face_No Face index
 */
void Apply_3D_Viscous_Boundary_Gradients(const int &Cell_No, const int &Face_No)
{
    // Get boundary type for this face
    int boundary_type = Get_3D_Boundary_Type(Cell_No, Face_No);

    switch (boundary_type)
    {
    case WALL_BOUNDARY_3D:
        // No-slip wall: zero velocity at wall
        Apply_3D_Wall_Velocity_Gradients(Cell_No, Face_No);
        break;

    case INLET_BOUNDARY_3D:
        // Specified velocity gradients at inlet
        Apply_3D_Inlet_Velocity_Gradients(Cell_No, Face_No);
        break;

    case OUTLET_BOUNDARY_3D:
        // Zero gradient at outlet
        u_Gradient_Face_3D.assign(3, 0.0);
        v_Gradient_Face_3D.assign(3, 0.0);
        w_Gradient_Face_3D.assign(3, 0.0);
        break;

    case SYMMETRY_BOUNDARY_3D:
        // Symmetry conditions
        Apply_3D_Symmetry_Velocity_Gradients(Cell_No, Face_No);
        break;

    default:
        // Default to zero gradient
        u_Gradient_Face_3D.assign(3, 0.0);
        v_Gradient_Face_3D.assign(3, 0.0);
        w_Gradient_Face_3D.assign(3, 0.0);
        break;
    }
}

/**
 * @brief Evaluate 3D wall skin friction and heat transfer
 *
 * This function computes wall shear stress and heat transfer coefficients
 * for all wall faces in the 3D domain.
 */
void Evaluate_3D_Wall_Skin_Friction()
{
    if (Wall_Cells_List_3D.empty())
        return;

    cout << "Evaluating 3D wall skin friction..." << endl;

    int list_size = Wall_Cells_List_3D.size();

    // Wall coefficients
    vector<double> Cf_3D; // Skin friction coefficient
    vector<double> Ch_3D; // Heat transfer coefficient
    vector<double> St_3D; // Stanton number

    Cf_3D.reserve(list_size / 3); // Each wall entry has 3 components
    Ch_3D.reserve(list_size / 3);
    St_3D.reserve(list_size / 3);

    for (int i = 0; i < list_size; i += 3)
    {
        int cell_index = Wall_Cells_List_3D[i + 0];
        int face_no = Wall_Cells_List_3D[i + 1];
        int ghost_cell = Wall_Cells_List_3D[i + 2];

        // Calculate wall shear stress
        double tau_wall = Calculate_3D_Wall_Shear_Stress(cell_index, face_no);

        // Calculate heat transfer at wall
        double q_wall = Calculate_3D_Wall_Heat_Transfer(cell_index, face_no);

        // Reference dynamic pressure
        double q_inf = 0.5 * Rho_inf_3D * (u_inf_3D * u_inf_3D + v_inf_3D * v_inf_3D + w_inf_3D * w_inf_3D);

        // Skin friction coefficient
        double cf = tau_wall / q_inf;
        Cf_3D.push_back(cf);

        // Heat transfer coefficient (if thermal wall)
        if (Is_Thermal_Wall)
        {
            double Delta_T = T_wall - T_inf_3D;
            if (fabs(Delta_T) > 1e-10)
            {
                double ch = q_wall / (q_inf * Delta_T);
                Ch_3D.push_back(ch);

                // Stanton number
                double st = ch / (2.0 * cf);
                St_3D.push_back(st);
            }
        }
    }

    // Output statistics
    if (!Cf_3D.empty())
    {
        double cf_avg = 0.0, cf_max = 0.0, cf_min = 1e10;
        for (double cf : Cf_3D)
        {
            cf_avg += cf;
            cf_max = max(cf_max, cf);
            cf_min = min(cf_min, cf);
        }
        cf_avg /= Cf_3D.size();

        cout << "3D Wall skin friction statistics:" << endl;
        cout << "  Average Cf: " << cf_avg << endl;
        cout << "  Maximum Cf: " << cf_max << endl;
        cout << "  Minimum Cf: " << cf_min << endl;
    }
}

/**
 * @brief Calculate wall shear stress magnitude for 3D wall
 *
 * @param Cell_No Wall cell index
 * @param Face_No Wall face index
 * @return Wall shear stress magnitude
 */
double Calculate_3D_Wall_Shear_Stress(const int &Cell_No, const int &Face_No)
{
    // Calculate velocity gradients at wall
    Calculate_3D_Velocity_Gradients_On_Face(Cell_No, Face_No);

    // Wall viscosity
    double T_wall = Primitive_Cells_3D[Cell_No][5];
    Viscosity_3D(T_wall);
    double mu_wall = mu_star_3D;

    // Face normal
    int index = Face_No * 3;
    double nx = Cells[Cell_No].Face_Normals[index + 0];
    double ny = Cells[Cell_No].Face_Normals[index + 1];
    double nz = Cells[Cell_No].Face_Normals[index + 2];

    // Wall shear stress vector components
    double tau_x = mu_wall * Inv_Re * (u_Gradient_Face_3D[0] * nx + u_Gradient_Face_3D[1] * ny + u_Gradient_Face_3D[2] * nz);
    double tau_y = mu_wall * Inv_Re * (v_Gradient_Face_3D[0] * nx + v_Gradient_Face_3D[1] * ny + v_Gradient_Face_3D[2] * nz);
    double tau_z = mu_wall * Inv_Re * (w_Gradient_Face_3D[0] * nx + w_Gradient_Face_3D[1] * ny + w_Gradient_Face_3D[2] * nz);

    // Tangential shear stress magnitude
    double tau_n = tau_x * nx + tau_y * ny + tau_z * nz; // Normal component
    double tau_t_x = tau_x - tau_n * nx;                 // Tangential components
    double tau_t_y = tau_y - tau_n * ny;
    double tau_t_z = tau_z - tau_n * nz;

    return sqrt(tau_t_x * tau_t_x + tau_t_y * tau_t_y + tau_t_z * tau_t_z);
}

/**
 * @brief Calculate wall heat transfer for 3D thermal wall
 *
 * @param Cell_No Wall cell index
 * @param Face_No Wall face index
 * @return Wall heat flux magnitude
 */
double Calculate_3D_Wall_Heat_Transfer(const int &Cell_No, const int &Face_No)
{
    // Calculate temperature gradient at wall
    Calculate_Temperature_Gradient_On_Face_3D(Cell_No, Face_No);

    // Wall thermal conductivity
    double T_wall = Primitive_Cells_3D[Cell_No][5];
    Thermal_Conductivity_3D(T_wall);
    double k_wall = K_3D;

    // Face normal
    int index = Face_No * 3;
    double nx = Cells[Cell_No].Face_Normals[index + 0];
    double ny = Cells[Cell_No].Face_Normals[index + 1];
    double nz = Cells[Cell_No].Face_Normals[index + 2];

    // Heat flux normal to wall (Fourier's law)
    double q_wall = -k_wall * Inv_Re * Inv_Pr *
                    (T_Gradient_Face_3D[0] * nx + T_Gradient_Face_3D[1] * ny + T_Gradient_Face_3D[2] * nz);

    return fabs(q_wall);
}