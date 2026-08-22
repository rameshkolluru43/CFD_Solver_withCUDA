#include "definitions.h"
#include "Boundary_Conditions.h"
#include "Globals.h"
#include "Utilities.h"

/**
 * @file Boundary_Conditions.cpp
 * @brief 3D Boundary conditions implementation for CFD solver
 *
 * This file implements various boundary conditions for 3D CFD computations including:
 * - Inlet conditions (subsonic/supersonic)
 * - Outlet conditions (subsonic/supersonic)
 * - Wall conditions (viscous/inviscid)
 * - Symmetry conditions
 * - Far-field conditions
 * - Front/back face conditions (3D extension)
 */

string BCFileName, InitCondFileName;

/**
 * @brief Main function to apply all boundary conditions in 3D
 *
 * This function coordinates the application of all boundary conditions
 * based on the flow regime and boundary types defined for the 3D domain.
 */
void Apply_Boundary_Conditions_3D()
{
    cout << "Applying 3D Boundary Conditions" << endl;

    // Apply inlet boundary conditions
    switch (Is_Inlet_SubSonic)
    {
    case true:
        Subsonic_Inlet_Condition_3D(inletCond, Inlet_Cells_List);
        break;
    case false:
        Supersonic_Inlet_Condition_3D(inletCond, Inlet_Cells_List);
        break;
    }

    // Apply exit boundary conditions
    switch (Is_Exit_SubSonic)
    {
    case true:
        Subsonic_Exit_Condition_3D(exitCond, Exit_Cells_List);
        break;
    case false:
        Supersonic_Exit_Condition_3D(exitCond, Exit_Cells_List);
        break;
    }

    // Apply wall boundary conditions
    switch (Is_Viscous_Wall)
    {
    case true:
        Viscous_Wall_Boundary_Condition_3D();
        break;
    case false:
        InViscid_Wall_Boundary_Condition_3D();
        break;
    }

    // Apply symmetry boundary conditions
    if (has_Symmetry_BC)
    {
        Symmetry_Boundary_Condition_3D();
    }

    // Apply far-field boundary conditions
    if (Enable_Far_Field)
    {
        Far_Field_Boundary_Condition_3D();
    }

    // Apply 3D-specific front/back boundary conditions
    Front_Back_Boundary_Condition_3D();

    cout << "3D Boundary Conditions Applied Successfully" << endl;
}

/**
 * @brief Apply subsonic inlet boundary condition for 3D flow
 * @param inletCond Inlet condition parameters
 * @param inlet_cells List of inlet boundary cells
 */
void Subsonic_Inlet_Condition_3D(const InletCondition &inletCond, const V_I &inlet_cells)
{
    for (size_t i = 0; i < inlet_cells.size(); i += 3)
    {
        int cell_id = inlet_cells[i];
        int face_id = inlet_cells[i + 1];
        int ghost_id = inlet_cells[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells)
        {
            // Get interior cell state
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            // Set inlet state with specified total conditions
            V_D U_inlet(NUM_CONSERVATIVE_VARS);
            V_D Prim_inlet(NUM_PRIMITIVE_VARS);

            // For subsonic inlet: specify total pressure, total temperature, and flow direction
            double gamma_m1 = GAMMA - 1.0;
            double M_inlet = inletCond.M;
            double T_total = inletCond.T;
            double P_total = inletCond.P;

            // Static properties from isentropic relations
            double T_static = T_total / (1.0 + 0.5 * gamma_m1 * M_inlet * M_inlet);
            double P_static = P_total * pow(T_static / T_total, GAMMA / gamma_m1);
            double rho_static = P_static / (R_GAS * T_static);
            double c_static = sqrt(GAMMA * P_static / rho_static);
            double V_mag = M_inlet * c_static;

            // Set velocity components (3D)
            Prim_inlet[PRIM_RHO] = rho_static;
            Prim_inlet[PRIM_U] = inletCond.u;
            Prim_inlet[PRIM_V] = inletCond.v;
            Prim_inlet[PRIM_W] = inletCond.w;
            Prim_inlet[PRIM_P] = P_static;

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_inlet, U_inlet);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_inlet;
                Primitive_Cells[ghost_id] = Prim_inlet;
            }
        }
    }
}

/**
 * @brief Apply supersonic inlet boundary condition for 3D flow
 * @param inletCond Inlet condition parameters
 * @param inlet_cells List of inlet boundary cells
 */
void Supersonic_Inlet_Condition_3D(const InletCondition &inletCond, const V_I &inlet_cells)
{
    for (size_t i = 0; i < inlet_cells.size(); i += 3)
    {
        int cell_id = inlet_cells[i];
        int face_id = inlet_cells[i + 1];
        int ghost_id = inlet_cells[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells)
        {
            // For supersonic inlet: all conditions are specified
            V_D U_inlet(NUM_CONSERVATIVE_VARS);
            V_D Prim_inlet(NUM_PRIMITIVE_VARS);

            // Set primitive variables
            Prim_inlet[PRIM_RHO] = inletCond.Rho;
            Prim_inlet[PRIM_U] = inletCond.u;
            Prim_inlet[PRIM_V] = inletCond.v;
            Prim_inlet[PRIM_W] = inletCond.w; // 3D component
            Prim_inlet[PRIM_P] = inletCond.P;

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_inlet, U_inlet);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_inlet;
                Primitive_Cells[ghost_id] = Prim_inlet;
            }
        }
    }
}

/**
 * @brief Apply subsonic exit boundary condition for 3D flow
 * @param exitCond Exit condition parameters
 * @param exit_cells List of exit boundary cells
 */
void Subsonic_Exit_Condition_3D(const ExitCondition &exitCond, const V_I &exit_cells)
{
    for (size_t i = 0; i < exit_cells.size(); i += 3)
    {
        int cell_id = exit_cells[i];
        int face_id = exit_cells[i + 1];
        int ghost_id = exit_cells[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells)
        {
            // Get interior cell state
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            // For subsonic exit: specify static pressure, extrapolate other variables
            V_D U_exit(NUM_CONSERVATIVE_VARS);
            V_D Prim_exit(NUM_PRIMITIVE_VARS);

            // Extrapolate density and velocity from interior
            Prim_exit[PRIM_RHO] = Prim_interior[PRIM_RHO];
            Prim_exit[PRIM_U] = Prim_interior[PRIM_U];
            Prim_exit[PRIM_V] = Prim_interior[PRIM_V];
            Prim_exit[PRIM_W] = Prim_interior[PRIM_W]; // 3D component

            // Specify exit pressure
            Prim_exit[PRIM_P] = exitCond.P;

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_exit, U_exit);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_exit;
                Primitive_Cells[ghost_id] = Prim_exit;
            }
        }
    }
}

/**
 * @brief Apply supersonic exit boundary condition for 3D flow
 * @param exitCond Exit condition parameters
 * @param exit_cells List of exit boundary cells
 */
void Supersonic_Exit_Condition_3D(const ExitCondition &exitCond, const V_I &exit_cells)
{
    for (size_t i = 0; i < exit_cells.size(); i += 3)
    {
        int cell_id = exit_cells[i];
        int face_id = exit_cells[i + 1];
        int ghost_id = exit_cells[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells)
        {
            // For supersonic exit: extrapolate all variables from interior
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            // Set ghost cell values equal to interior
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_interior;
                Primitive_Cells[ghost_id] = Prim_interior;
            }
        }
    }
}

/**
 * @brief Apply viscous wall boundary condition for 3D flow
 */
void Viscous_Wall_Boundary_Condition_3D()
{
    for (size_t i = 0; i < Wall_Cells_List.size(); i += 3)
    {
        int cell_id = Wall_Cells_List[i];
        int face_id = Wall_Cells_List[i + 1];
        int ghost_id = Wall_Cells_List[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells)
        {
            // Get interior cell state
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            // Apply no-slip condition: all velocity components are zero
            V_D U_wall(NUM_CONSERVATIVE_VARS);
            V_D Prim_wall(NUM_PRIMITIVE_VARS);

            // Reflect density and pressure, set velocities to zero
            Prim_wall[PRIM_RHO] = Prim_interior[PRIM_RHO];
            Prim_wall[PRIM_U] = 0.0; // No-slip condition
            Prim_wall[PRIM_V] = 0.0; // No-slip condition
            Prim_wall[PRIM_W] = 0.0; // No-slip condition (3D)
            Prim_wall[PRIM_P] = Prim_interior[PRIM_P];

            // For isothermal wall
            if (wallCond.T > 0.0)
            {
                double T_wall = wallCond.T;
                Prim_wall[PRIM_P] = Prim_wall[PRIM_RHO] * R_GAS * T_wall;
            }

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_wall, U_wall);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_wall;
                Primitive_Cells[ghost_id] = Prim_wall;
            }
        }
    }
}

/**
 * @brief Apply inviscid wall boundary condition for 3D flow
 */
void InViscid_Wall_Boundary_Condition_3D()
{
    for (size_t i = 0; i < Wall_Cells_List.size(); i += 3)
    {
        int cell_id = Wall_Cells_List[i];
        int face_id = Wall_Cells_List[i + 1];
        int ghost_id = Wall_Cells_List[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells && face_id < NUM_FACES_3D)
        {
            // Get interior cell state
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            // Get face normal
            double nx = Cells[cell_id].Face_Normals[3 * face_id];
            double ny = Cells[cell_id].Face_Normals[3 * face_id + 1];
            double nz = Cells[cell_id].Face_Normals[3 * face_id + 2];

            // Apply slip condition: reflect normal velocity component
            V_D U_wall(NUM_CONSERVATIVE_VARS);
            V_D Prim_wall(NUM_PRIMITIVE_VARS);

            double u_int = Prim_interior[PRIM_U];
            double v_int = Prim_interior[PRIM_V];
            double w_int = Prim_interior[PRIM_W];

            // Normal velocity component
            double Vn = u_int * nx + v_int * ny + w_int * nz;

            // Reflect normal component, keep tangential components
            Prim_wall[PRIM_RHO] = Prim_interior[PRIM_RHO];
            Prim_wall[PRIM_U] = u_int - 2.0 * Vn * nx;
            Prim_wall[PRIM_V] = v_int - 2.0 * Vn * ny;
            Prim_wall[PRIM_W] = w_int - 2.0 * Vn * nz; // 3D component
            Prim_wall[PRIM_P] = Prim_interior[PRIM_P];

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_wall, U_wall);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_wall;
                Primitive_Cells[ghost_id] = Prim_wall;
            }
        }
    }
}

/**
 * @brief Apply symmetry boundary condition for 3D flow
 */
void Symmetry_Boundary_Condition_3D()
{
    for (size_t i = 0; i < Symmetry_Cells_List.size(); i += 3)
    {
        int cell_id = Symmetry_Cells_List[i];
        int face_id = Symmetry_Cells_List[i + 1];
        int ghost_id = Symmetry_Cells_List[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells && face_id < NUM_FACES_3D)
        {
            // Get interior cell state
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            // Get face normal
            double nx = Cells[cell_id].Face_Normals[3 * face_id];
            double ny = Cells[cell_id].Face_Normals[3 * face_id + 1];
            double nz = Cells[cell_id].Face_Normals[3 * face_id + 2];

            // Apply symmetry condition: zero normal velocity gradient
            V_D U_symm(NUM_CONSERVATIVE_VARS);
            V_D Prim_symm(NUM_PRIMITIVE_VARS);

            double u_int = Prim_interior[PRIM_U];
            double v_int = Prim_interior[PRIM_V];
            double w_int = Prim_interior[PRIM_W];

            // Normal velocity component
            double Vn = u_int * nx + v_int * ny + w_int * nz;

            // Mirror normal velocity, preserve tangential components
            Prim_symm[PRIM_RHO] = Prim_interior[PRIM_RHO];
            Prim_symm[PRIM_U] = u_int - 2.0 * Vn * nx;
            Prim_symm[PRIM_V] = v_int - 2.0 * Vn * ny;
            Prim_symm[PRIM_W] = w_int - 2.0 * Vn * nz;
            Prim_symm[PRIM_P] = Prim_interior[PRIM_P];

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_symm, U_symm);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_symm;
                Primitive_Cells[ghost_id] = Prim_symm;
            }
        }
    }
}

/**
 * @brief Apply far-field boundary condition for 3D flow
 */
void Far_Field_Boundary_Condition_3D()
{
    for (size_t i = 0; i < Far_Field_Out_Flow_List.size(); i += 3)
    {
        int cell_id = Far_Field_Out_Flow_List[i];
        int face_id = Far_Field_Out_Flow_List[i + 1];
        int ghost_id = Far_Field_Out_Flow_List[i + 2];

        if (cell_id >= 0 && cell_id < No_Physical_Cells)
        {
            // Apply far-field conditions based on local Mach number
            V_D U_interior = U_Cells[cell_id];
            V_D Prim_interior = Primitive_Cells[cell_id];

            double M_local = Calculate_Mach_Number_3D(U_interior);

            V_D U_farfield(NUM_CONSERVATIVE_VARS);
            V_D Prim_farfield(NUM_PRIMITIVE_VARS);

            if (M_local < 1.0)
            {
                // Subsonic: use Riemann invariants
                Apply_Subsonic_Far_Field_3D(Prim_interior, Prim_farfield, face_id);
            }
            else
            {
                // Supersonic: use upwind based on flow direction
                Apply_Supersonic_Far_Field_3D(Prim_interior, Prim_farfield, face_id);
            }

            // Convert to conservative variables
            Convert_Primitive_to_Conservative_3D(Prim_farfield, U_farfield);

            // Set ghost cell values
            if (ghost_id >= 0 && ghost_id < Total_No_Cells)
            {
                U_Cells[ghost_id] = U_farfield;
                Primitive_Cells[ghost_id] = Prim_farfield;
            }
        }
    }
}

/**
 * @brief Apply front/back boundary conditions specific to 3D
 */
void Front_Back_Boundary_Condition_3D()
{
    // Front boundary (z-max)
    for (size_t i = 0; i < Front_Cells_List.size(); i += 3)
    {
        int cell_id = Front_Cells_List[i];
        int face_id = Front_Cells_List[i + 1];
        int ghost_id = Front_Cells_List[i + 2];

        // Apply appropriate boundary condition (symmetry, wall, etc.)
        Apply_Z_Direction_Boundary_3D(cell_id, face_id, ghost_id, true);
    }

    // Back boundary (z-min)
    for (size_t i = 0; i < Back_Cells_List.size(); i += 3)
    {
        int cell_id = Back_Cells_List[i];
        int face_id = Back_Cells_List[i + 1];
        int ghost_id = Back_Cells_List[i + 2];

        // Apply appropriate boundary condition (symmetry, wall, etc.)
        Apply_Z_Direction_Boundary_3D(cell_id, face_id, ghost_id, false);
    }
}

/**
 * @brief Apply boundary condition in z-direction
 * @param cell_id Interior cell ID
 * @param face_id Face ID
 * @param ghost_id Ghost cell ID
 * @param is_front True for front face (z-max), false for back face (z-min)
 */
void Apply_Z_Direction_Boundary_3D(const int &cell_id, const int &face_id,
                                   const int &ghost_id, const bool &is_front)
{
    if (cell_id >= 0 && cell_id < No_Physical_Cells)
    {
        V_D U_interior = U_Cells[cell_id];
        V_D Prim_interior = Primitive_Cells[cell_id];

        V_D U_boundary(NUM_CONSERVATIVE_VARS);
        V_D Prim_boundary(NUM_PRIMITIVE_VARS);

        // Apply symmetry condition in z-direction
        Prim_boundary[PRIM_RHO] = Prim_interior[PRIM_RHO];
        Prim_boundary[PRIM_U] = Prim_interior[PRIM_U];
        Prim_boundary[PRIM_V] = Prim_interior[PRIM_V];
        Prim_boundary[PRIM_W] = -Prim_interior[PRIM_W]; // Reflect w-velocity
        Prim_boundary[PRIM_P] = Prim_interior[PRIM_P];

        // Convert to conservative variables
        Convert_Primitive_to_Conservative_3D(Prim_boundary, U_boundary);

        // Set ghost cell values
        if (ghost_id >= 0 && ghost_id < Total_No_Cells)
        {
            U_Cells[ghost_id] = U_boundary;
            Primitive_Cells[ghost_id] = Prim_boundary;
        }
    }
}

/**
 * @brief Apply general boundary condition for a face in 3D
 * @param Cell_No Cell index
 * @param Face_No Face index
 * @param U_boundary Output boundary state
 * @param Prim_boundary Output boundary primitive variables
 */
void Apply_Boundary_Condition_3D(const int &Cell_No, const int &Face_No,
                                 V_D &U_boundary, V_D &Prim_boundary)
{
    // Default: extrapolate from interior
    U_boundary = U_Cells[Cell_No];
    Prim_boundary = Primitive_Cells[Cell_No];

    // Check if this is a specific boundary face and apply appropriate condition
    if (Cells[Cell_No].has_Inlet_Face)
    {
        // Apply inlet condition
        if (Is_Inlet_SubSonic)
        {
            // Apply subsonic inlet logic
            Prim_boundary[PRIM_RHO] = inletCond.Rho;
            Prim_boundary[PRIM_U] = inletCond.u;
            Prim_boundary[PRIM_V] = inletCond.v;
            Prim_boundary[PRIM_W] = inletCond.w;
            Prim_boundary[PRIM_P] = inletCond.P;
        }
    }
    else if (Cells[Cell_No].has_Exit_Face)
    {
        // Apply exit condition
        if (Is_Exit_SubSonic)
        {
            Prim_boundary[PRIM_P] = exitCond.P;
        }
    }
    else if (Cells[Cell_No].has_Wall_Face)
    {
        // Apply wall condition
        if (Is_Viscous_Wall)
        {
            Prim_boundary[PRIM_U] = 0.0;
            Prim_boundary[PRIM_V] = 0.0;
            Prim_boundary[PRIM_W] = 0.0;
        }
        else
        {
            // Reflect normal velocity component
            double nx = Cells[Cell_No].Face_Normals[3 * Face_No];
            double ny = Cells[Cell_No].Face_Normals[3 * Face_No + 1];
            double nz = Cells[Cell_No].Face_Normals[3 * Face_No + 2];

            double u = Prim_boundary[PRIM_U];
            double v = Prim_boundary[PRIM_V];
            double w = Prim_boundary[PRIM_W];
            double Vn = u * nx + v * ny + w * nz;

            Prim_boundary[PRIM_U] = u - 2.0 * Vn * nx;
            Prim_boundary[PRIM_V] = v - 2.0 * Vn * ny;
            Prim_boundary[PRIM_W] = w - 2.0 * Vn * nz;
        }
    }

    // Convert back to conservative variables
    Convert_Primitive_to_Conservative_3D(Prim_boundary, U_boundary);
}

/**
 * @brief Apply subsonic far-field condition
 * @param Prim_interior Interior primitive variables
 * @param Prim_farfield Output far-field primitive variables
 * @param face_id Face identifier
 */
void Apply_Subsonic_Far_Field_3D(const V_D &Prim_interior, V_D &Prim_farfield, const int &face_id)
{
    // Use Riemann invariants for subsonic far-field
    double rho_inf = Rho_inf;
    double u_inf = u_inf;
    double v_inf = v_inf;
    double w_inf = w_inf;
    double p_inf = P_inf;

    // Extrapolate based on characteristic speeds
    Prim_farfield[PRIM_RHO] = rho_inf;
    Prim_farfield[PRIM_U] = u_inf;
    Prim_farfield[PRIM_V] = v_inf;
    Prim_farfield[PRIM_W] = w_inf;
    Prim_farfield[PRIM_P] = p_inf;
}

/**
 * @brief Apply supersonic far-field condition
 * @param Prim_interior Interior primitive variables
 * @param Prim_farfield Output far-field primitive variables
 * @param face_id Face identifier
 */
void Apply_Supersonic_Far_Field_3D(const V_D &Prim_interior, V_D &Prim_farfield, const int &face_id)
{
    // For supersonic flow, use upwind based on flow direction
    // This is simplified - in practice, check normal velocity direction

    double rho_inf = Rho_inf;
    double u_inf = u_inf;
    double v_inf = v_inf;
    double w_inf = w_inf;
    double p_inf = P_inf;

    Prim_farfield[PRIM_RHO] = rho_inf;
    Prim_farfield[PRIM_U] = u_inf;
    Prim_farfield[PRIM_V] = v_inf;
    Prim_farfield[PRIM_W] = w_inf;
    Prim_farfield[PRIM_P] = p_inf;
}