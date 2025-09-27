/**
 * @file Initialize.cpp
 * @brief 3D CFD solver initialization module
 *
 * This module handles the initialization of all data structures, memory allocation,
 * and initial conditions setup for 3D CFD simulations. It extends the 2D initialization
 * to handle hexahedral cells, 3D flow variables, and 6-face connectivity.
 *
 * Key Features:
 * - 3D memory allocation for hexahedral cells
 * - 5-component conservative variables [ρ, ρu, ρv, ρw, ρE]
 * - 3D primitive variables [ρ, u, v, w, p, T, a, h, μ, λ]
 * - 6-face boundary identification and setup
 * - 3D gradient computation stencils
 * - Support for multiple initialization methods
 * - 3D reference state calculations
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"

// 3D Global solution vectors
vector<V_D> U_Cells_3D, Cells_Net_Flux_3D, Cells_DelU_3D;
vector<V_D> U_Cells_RK_1_3D, U_Cells_RK_2_3D, U_Cells_RK_3_3D, U_Cells_RK_4_3D; // 3D RK storage

// 3D Flux vectors for each coordinate direction
vector<V_D> A_x_3D, A_y_3D, A_z_3D;                                     // 3D flux Jacobians
vector<V_D> A_x_L_3D, A_x_R_3D, A_y_T_3D, A_y_B_3D, A_z_F_3D, A_z_B_3D; // 3D face flux Jacobians

// 3D Boundary and mesh vectors
vector<bool> v_bool_3D(NUM_FACES_3D, false);      // 6-face boolean vector
vector<vector<bool>> Cells_Face_Boundary_Type_3D; // 6-face boundary types per cell

// 3D Reference and flow parameters
double P_ref_3D, Rho_ref_3D, u_ref_3D, v_ref_3D, w_ref_3D, T_ref_3D, a_ref_3D; // 3D reference state
double u_inf_3D, v_inf_3D, w_inf_3D, M_inf_3D, P_inf_3D, T_inf_3D, Rho_inf_3D; // 3D freestream

// 3D Grid and solver parameters
int nx_3D, ny_3D, nz_3D; // 3D grid dimensions
bool Is_3D_Simulation = true;

// 3D Viscous flow vectors
vector<V_D> Cells_Viscous_Flux_3D;                               // 3D viscous fluxes
vector<V_D> Cells_u_Grad_3D, Cells_v_Grad_3D, Cells_w_Grad_3D;   // 3D velocity gradients
vector<V_D> Cells_T_Grad_3D, Cells_Rho_Grad_3D, Cells_P_Grad_3D; // 3D scalar gradients

// 3D State vectors (Left/Right/Top/Bottom/Front/Back faces)
double Rho_L_3D, T_L_3D, P_L_3D, u_L_3D, v_L_3D, w_L_3D, M_L_3D, C_L_3D, H_L_3D;
double Rho_R_3D, T_R_3D, P_R_3D, u_R_3D, v_R_3D, w_R_3D, M_R_3D, C_R_3D, H_R_3D;
V_D Flux_L_3D(NUM_CONSERVATIVE_3D, 0.0), Flux_R_3D(NUM_CONSERVATIVE_3D, 0.0);
V_D U_L_3D(NUM_CONSERVATIVE_3D, 0.0), U_R_3D(NUM_CONSERVATIVE_3D, 0.0);

// 3D Global working vectors
V_D Global_Primitive_3D(NUM_PRIMITIVE_3D, 0.0); // 3D primitive variables
V_D Global_U_3D(NUM_CONSERVATIVE_3D, 0.0);      // 3D conservative variables
V_D Global_Flux_3D(NUM_CONSERVATIVE_3D, 0.0);   // 3D flux vector

// 3D Gradient vectors
V_D u_Gradient_3D(3, 0.0), v_Gradient_3D(3, 0.0), w_Gradient_3D(3, 0.0);   // 3D velocity gradients
V_D T_Gradient_3D(3, 0.0), Rho_Gradient_3D(3, 0.0), P_Gradient_3D(3, 0.0); // 3D scalar gradients
V_D Viscous_Flux_3D(NUM_CONSERVATIVE_3D, 0.0);                             // 3D viscous flux

// 3D Working arrays for primitive variables
vector<V_D> Primitive_Cells_3D; // 3D primitive variables for all cells

/**
 * @brief Main 3D initialization function
 *
 * This function allocates memory for all 3D data structures and initializes
 * the computational domain for 3D CFD simulations.
 *
 * @param Test_Case Test case identifier for initialization
 */
void Initialize_3D(const int &Test_Case)
{
    cout << "\n========================================" << endl;
    cout << "3D CFD SOLVER INITIALIZATION" << endl;
    cout << "========================================" << endl;

    cout << "Allocating Memory for 3D Physical Cells: " << No_Physical_Cells << endl;
    cout << "Total Number of 3D Cells (including Ghost): " << Total_No_Cells << endl;
    cout << "3D Grid Dimensions: " << nx_3D << " x " << ny_3D << " x " << nz_3D << endl;

    // Resize global working vectors for 3D
    Global_Primitive_3D.assign(NUM_PRIMITIVE_3D, 0.0);
    Global_U_3D.assign(NUM_CONSERVATIVE_3D, 0.0);
    Global_Flux_3D.assign(NUM_CONSERVATIVE_3D, 0.0);

    // 3D flux vectors initialization
    Flux_L_3D.assign(NUM_CONSERVATIVE_3D, 0.0);
    Flux_R_3D.assign(NUM_CONSERVATIVE_3D, 0.0);
    U_L_3D.assign(NUM_CONSERVATIVE_3D, 0.0);
    U_R_3D.assign(NUM_CONSERVATIVE_3D, 0.0);

    // 3D gradient vectors
    u_Gradient_3D.assign(3, 0.0);   // du/dx, du/dy, du/dz
    v_Gradient_3D.assign(3, 0.0);   // dv/dx, dv/dy, dv/dz
    w_Gradient_3D.assign(3, 0.0);   // dw/dx, dw/dy, dw/dz
    T_Gradient_3D.assign(3, 0.0);   // dT/dx, dT/dy, dT/dz
    Rho_Gradient_3D.assign(3, 0.0); // dρ/dx, dρ/dy, dρ/dz
    P_Gradient_3D.assign(3, 0.0);   // dp/dx, dp/dy, dp/dz

    // Viscous flux vector
    Viscous_Flux_3D.assign(NUM_CONSERVATIVE_3D, 0.0);

    // Boolean vector for 6 faces
    v_bool_3D.assign(NUM_FACES_3D, false);

    cout << "Initializing 3D cell-based vectors..." << endl;

    // Reserve memory for efficiency
    U_Cells_3D.reserve(Total_No_Cells);
    Cells_Net_Flux_3D.reserve(Total_No_Cells);
    Cells_DelU_3D.reserve(Total_No_Cells);
    Primitive_Cells_3D.reserve(Total_No_Cells);
    Cells_Face_Boundary_Type_3D.reserve(Total_No_Cells);

    // Runge-Kutta storage vectors
    U_Cells_RK_1_3D.reserve(Total_No_Cells);
    U_Cells_RK_2_3D.reserve(Total_No_Cells);
    U_Cells_RK_3_3D.reserve(Total_No_Cells);
    U_Cells_RK_4_3D.reserve(Total_No_Cells);

    // Initialize cell-based vectors
    for (int cell_idx = 0; cell_idx < Total_No_Cells; cell_idx++)
    {
        // Conservative variables
        U_Cells_3D.push_back(Global_U_3D);
        U_Cells_RK_1_3D.push_back(Global_U_3D);
        U_Cells_RK_2_3D.push_back(Global_U_3D);
        U_Cells_RK_3_3D.push_back(Global_U_3D);
        U_Cells_RK_4_3D.push_back(Global_U_3D);

        // Flux and solution increment vectors
        Cells_Net_Flux_3D.push_back(Global_Flux_3D);
        Cells_DelU_3D.push_back(Global_U_3D);

        // Primitive variables
        Primitive_Cells_3D.push_back(Global_Primitive_3D);

        // 6-face boundary type information
        Cells_Face_Boundary_Type_3D.push_back(v_bool_3D);
    }

    // Initialize viscous flow vectors if needed
    if (Is_Viscous_Wall)
    {
        cout << "Initializing 3D viscous flow vectors..." << endl;

        Cells_Viscous_Flux_3D.reserve(Total_No_Cells);
        Cells_u_Grad_3D.reserve(Total_No_Cells);
        Cells_v_Grad_3D.reserve(Total_No_Cells);
        Cells_w_Grad_3D.reserve(Total_No_Cells);
        Cells_T_Grad_3D.reserve(Total_No_Cells);
        Cells_Rho_Grad_3D.reserve(Total_No_Cells);
        Cells_P_Grad_3D.reserve(Total_No_Cells);

        for (int cell_idx = 0; cell_idx < Total_No_Cells; cell_idx++)
        {
            Cells_Viscous_Flux_3D.push_back(Viscous_Flux_3D);
            Cells_u_Grad_3D.push_back(u_Gradient_3D);
            Cells_v_Grad_3D.push_back(v_Gradient_3D);
            Cells_w_Grad_3D.push_back(w_Gradient_3D); // 3D w-velocity gradient
            Cells_T_Grad_3D.push_back(T_Gradient_3D);
            Cells_Rho_Grad_3D.push_back(Rho_Gradient_3D);
            Cells_P_Grad_3D.push_back(P_Gradient_3D);
        }

        // Initialize 3D gradient computation stencils
        for (int cell_idx = 0; cell_idx < No_Physical_Cells; cell_idx++)
        {
            Setup_3D_Gradient_Stencil(cell_idx);
        }
    }

    // Initialize 3D flux Jacobian matrices if needed for implicit methods
    if (Is_Implicit_Method)
    {
        cout << "Initializing 3D flux Jacobian matrices..." << endl;

        A_x_3D.reserve(Total_No_Cells);
        A_y_3D.reserve(Total_No_Cells);
        A_z_3D.reserve(Total_No_Cells); // 3D z-direction Jacobian

        // Face-based Jacobians for 6 faces
        A_x_L_3D.reserve(Total_No_Cells); // Left face (x-direction)
        A_x_R_3D.reserve(Total_No_Cells); // Right face (x-direction)
        A_y_T_3D.reserve(Total_No_Cells); // Top face (y-direction)
        A_y_B_3D.reserve(Total_No_Cells); // Bottom face (y-direction)
        A_z_F_3D.reserve(Total_No_Cells); // Front face (z-direction)
        A_z_B_3D.reserve(Total_No_Cells); // Back face (z-direction)

        V_D Jacobian_3D(NUM_CONSERVATIVE_3D * NUM_CONSERVATIVE_3D, 0.0); // 5x5 Jacobian

        for (int cell_idx = 0; cell_idx < Total_No_Cells; cell_idx++)
        {
            A_x_3D.push_back(Jacobian_3D);
            A_y_3D.push_back(Jacobian_3D);
            A_z_3D.push_back(Jacobian_3D); // 3D extension

            A_x_L_3D.push_back(Jacobian_3D);
            A_x_R_3D.push_back(Jacobian_3D);
            A_y_T_3D.push_back(Jacobian_3D);
            A_y_B_3D.push_back(Jacobian_3D);
            A_z_F_3D.push_back(Jacobian_3D); // 3D front face
            A_z_B_3D.push_back(Jacobian_3D); // 3D back face
        }
    }

    // Initialize 3D boundary conditions
    Setup_3D_Boundary_Conditions();

    // Initialize 3D reference state
    Setup_3D_Reference_State(Test_Case);

    cout << "3D Memory allocation completed successfully!" << endl;
    cout << "Total memory allocated:" << endl;
    cout << "- Conservative variables: " << Total_No_Cells * NUM_CONSERVATIVE_3D * sizeof(double) / (1024 * 1024) << " MB" << endl;
    cout << "- Primitive variables: " << Total_No_Cells * NUM_PRIMITIVE_3D * sizeof(double) / (1024 * 1024) << " MB" << endl;

    if (Is_Viscous_Wall)
    {
        cout << "- Viscous gradients: " << Total_No_Cells * 6 * 3 * sizeof(double) / (1024 * 1024) << " MB" << endl;
    }

    cout << "========================================\n"
         << endl;
}

/**
 * @brief Setup 3D gradient computation stencil for viscous terms
 *
 * This function identifies the neighbor cells required for 3D gradient
 * computations using Green-Gauss or least squares methods.
 *
 * @param Cell_No Current cell index
 */
void Setup_3D_Gradient_Stencil(const int &Cell_No)
{
    // 3D stencil requires 6 primary neighbors plus potentially more for accuracy
    vector<int> gradient_stencil;

    // Primary 6 neighbors (face neighbors)
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor = Cells[Cell_No].Neighbours[face];
        if (neighbor >= 0 && neighbor < Total_No_Cells)
        {
            gradient_stencil.push_back(neighbor);
        }
    }

    // Extended stencil for higher accuracy (edge and vertex neighbors)
    if (High_Order_Gradients)
    {
        // Add edge neighbors (12 edges for hexahedron)
        Add_3D_Edge_Neighbors(Cell_No, gradient_stencil);

        // Add vertex neighbors (8 vertices for hexahedron)
        Add_3D_Vertex_Neighbors(Cell_No, gradient_stencil);
    }

    // Store stencil information in cell structure
    Cells[Cell_No].Gradient_Stencil = gradient_stencil;
    Cells[Cell_No].Stencil_Size = gradient_stencil.size();
}

/**
 * @brief Setup 3D boundary conditions
 *
 * This function initializes boundary condition types for all 6 faces
 * of hexahedral cells and sets up boundary-specific parameters.
 */
void Setup_3D_Boundary_Conditions()
{
    cout << "Setting up 3D boundary conditions..." << endl;

    // Initialize boundary condition arrays for 6 faces
    // Face order: Left(x-), Right(x+), Bottom(y-), Top(y+), Back(z-), Front(z+)

    for (int cell = 0; cell < Total_No_Cells; cell++)
    {
        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            // Default to internal face
            Cells_Face_Boundary_Type_3D[cell][face] = false;

            // Check if this face is on a boundary
            int neighbor = Cells[cell].Neighbours[face];
            if (neighbor < 0)
            { // Boundary face (negative neighbor index)
                Cells_Face_Boundary_Type_3D[cell][face] = true;

                // Set boundary condition type based on face and location
                Set_3D_Boundary_Type(cell, face, neighbor);
            }
        }
    }

    // Setup inlet conditions (3D velocity components)
    if (Inlet_Conditions_3D.size() > 0)
    {
        for (auto &inlet : Inlet_Conditions_3D)
        {
            inlet.velocity_magnitude = sqrt(inlet.u * inlet.u + inlet.v * inlet.v + inlet.w * inlet.w);
            inlet.mach_number = inlet.velocity_magnitude / inlet.speed_of_sound;
        }
    }

    // Setup outlet conditions
    if (Outlet_Conditions_3D.size() > 0)
    {
        for (auto &outlet : Outlet_Conditions_3D)
        {
            // Calculate outlet properties based on pressure ratio
            Calculate_3D_Outlet_Properties(outlet);
        }
    }

    cout << "3D boundary conditions setup completed." << endl;
}

/**
 * @brief Setup 3D reference state for non-dimensionalization
 *
 * This function establishes reference values for 3D flow variables
 * used in non-dimensional formulations.
 *
 * @param Test_Case Test case identifier
 */
void Setup_3D_Reference_State(const int &Test_Case)
{
    cout << "Establishing 3D reference state..." << endl;

    switch (Test_Case)
    {
    case 1: // 3D Shock tube
        Setup_3D_Shock_Tube_Reference();
        break;

    case 2: // Flow over sphere
        Setup_3D_Sphere_Reference();
        break;

    case 3: // 3D Channel flow
        Setup_3D_Channel_Reference();
        break;

    case 4: // 3D Supersonic wedge
        Setup_3D_Wedge_Reference();
        break;

    default:
        // Default freestream reference state
        Setup_3D_Freestream_Reference();
        break;
    }

    // Calculate reference derived quantities
    u_ref_3D = sqrt(u_inf_3D * u_inf_3D + v_inf_3D * v_inf_3D + w_inf_3D * w_inf_3D); // 3D velocity magnitude
    a_ref_3D = sqrt(gamma * R_gas * T_ref_3D);                                        // Reference speed of sound

    // Reference Reynolds and Prandtl numbers
    if (Is_Viscous_Wall)
    {
        Re = (Rho_ref_3D * u_ref_3D * L_ref) / mu_ref;
        Pr = (mu_ref * cp_ref) / K_ref;
        Inv_Re = 1.0 / Re;
        Inv_Pr = 1.0 / Pr;
    }

    // Non-dimensional reference values
    if (Non_Dimensional_Form)
    {
        // Normalize all reference quantities
        Rho_ref_3D = 1.0;
        u_ref_3D = 1.0;
        P_ref_3D = 1.0;
        T_ref_3D = 1.0;
        L_ref = 1.0;
    }

    cout << "3D Reference state:" << endl;
    cout << "  Density: " << Rho_ref_3D << " kg/m³" << endl;
    cout << "  Velocity: (" << u_inf_3D << ", " << v_inf_3D << ", " << w_inf_3D << ") m/s" << endl;
    cout << "  Pressure: " << P_ref_3D << " Pa" << endl;
    cout << "  Temperature: " << T_ref_3D << " K" << endl;
    cout << "  Mach number: " << M_inf_3D << endl;

    if (Is_Viscous_Wall)
    {
        cout << "  Reynolds number: " << Re << endl;
        cout << "  Prandtl number: " << Pr << endl;
    }
}

/**
 * @brief Initialize solution with test case specific conditions
 *
 * This function sets initial conditions for various 3D test cases.
 *
 * @param Test_Case Test case identifier
 */
void Initialize_3D_Solution(const int &Test_Case)
{
    cout << "Initializing 3D solution for test case " << Test_Case << "..." << endl;

    for (int cell = 0; cell < No_Physical_Cells; cell++)
    {
        switch (Test_Case)
        {
        case 1:
            Initialize_3D_Shock_Tube(cell);
            break;

        case 2:
            Initialize_3D_Uniform_Flow(cell);
            break;

        case 3:
            Initialize_3D_Channel_Flow(cell);
            break;

        case 4:
            Initialize_3D_Vortex(cell);
            break;

        default:
            Initialize_3D_Freestream(cell);
            break;
        }

        // Convert primitive to conservative variables
        Convert_Primitive_to_Conservative_3D(cell);

        // Initialize time step
        Evaluate_Time_Step_3D(cell);
    }

    cout << "3D solution initialization completed." << endl;
}

/**
 * @brief Convert primitive to conservative variables for 3D
 *
 * @param cell_no Cell index
 */
void Convert_Primitive_to_Conservative_3D(const int &cell_no)
{
    double rho = Primitive_Cells_3D[cell_no][0];
    double u = Primitive_Cells_3D[cell_no][1];
    double v = Primitive_Cells_3D[cell_no][2];
    double w = Primitive_Cells_3D[cell_no][3]; // 3D w-velocity
    double p = Primitive_Cells_3D[cell_no][4];

    // Conservative variables: [ρ, ρu, ρv, ρw, ρE]
    U_Cells_3D[cell_no][0] = rho;                                 // Density
    U_Cells_3D[cell_no][1] = rho * u;                             // x-momentum
    U_Cells_3D[cell_no][2] = rho * v;                             // y-momentum
    U_Cells_3D[cell_no][3] = rho * w;                             // z-momentum (3D)
    U_Cells_3D[cell_no][4] = p / (gamma - 1.0) +                  // Total energy
                             0.5 * rho * (u * u + v * v + w * w); // 3D kinetic energy
}

/**
 * @brief Initialize uniform freestream conditions
 *
 * @param cell_no Cell index
 */
void Initialize_3D_Freestream(const int &cell_no)
{
    Primitive_Cells_3D[cell_no][0] = Rho_inf_3D; // Density
    Primitive_Cells_3D[cell_no][1] = u_inf_3D;   // u-velocity
    Primitive_Cells_3D[cell_no][2] = v_inf_3D;   // v-velocity
    Primitive_Cells_3D[cell_no][3] = w_inf_3D;   // w-velocity (3D)
    Primitive_Cells_3D[cell_no][4] = P_inf_3D;   // Pressure
    Primitive_Cells_3D[cell_no][5] = T_inf_3D;   // Temperature
    Primitive_Cells_3D[cell_no][6] = a_ref_3D;   // Speed of sound

    // Calculate derived quantities
    double velocity_mag = sqrt(u_inf_3D * u_inf_3D + v_inf_3D * v_inf_3D + w_inf_3D * w_inf_3D);
    Primitive_Cells_3D[cell_no][7] = Calculate_Enthalpy_3D(cell_no); // Enthalpy

    if (Is_Viscous_Wall)
    {
        Primitive_Cells_3D[cell_no][8] = Calculate_Viscosity_3D(T_inf_3D);            // Viscosity
        Primitive_Cells_3D[cell_no][9] = Calculate_Thermal_Conductivity_3D(T_inf_3D); // Thermal conductivity
    }
}