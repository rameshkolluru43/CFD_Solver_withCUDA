#include "definitions.h"
#include "Globals.h"
#include "Primitive_Computational.h"
#include "Viscous_Functions.h"
#include "Utilities.h"

/**
 * @file Basic_Functions.cpp
 * @brief Basic computational functions for 3D CFD solver
 *
 * This file contains fundamental functions for converting between
 * conservative and primitive variables, calculating thermodynamic
 * properties, and handling basic vector operations in 3D.
 */

// Function to set the delta U values for each cell (extended for 3D - 5 components)
// dU: Vector containing delta U values
void Set_DelU(V_D &dU)
{
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        for (int i = 0; i < NUM_CONSERVATIVE_VARS; i++) // 5 components for 3D
        {
            Cells_DelU[Cell_No][i] = dU[NUM_CONSERVATIVE_VARS * Cell_No + i];
        }
    }
}

// Function to calculate primitive variables from conservative variables (3D version)
// Cell_No: Cell index
// U_Vect: Vector of conservative variables [rho, rho*u, rho*v, rho*w, rho*E]
void Calculate_Primitive_Variables(const int &Cell_No, V_D &U_Vect)
{
    double vmag = 0.0, inv_Density = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0;
    double Pressure = 0.0, Temperature = 0.0, C = 0.0, Rho = 0.0;

    Rho = U_Vect[RHO_INDEX];
    inv_Density = 1.0 / Rho;

    // Extract velocity components (3D)
    v1 = U_Vect[RHU_INDEX] * inv_Density; // u-velocity
    v2 = U_Vect[RHV_INDEX] * inv_Density; // v-velocity
    v3 = U_Vect[RHW_INDEX] * inv_Density; // w-velocity (3D extension)

    // Apply small velocity cutoff
    if (fabs(v1) <= EPSILON)
        v1 = 0.0;
    if (fabs(v2) <= EPSILON)
        v2 = 0.0;
    if (fabs(v3) <= EPSILON)
        v3 = 0.0;

    // Calculate velocity magnitude (3D)
    vmag = (v1 * v1 + v2 * v2 + v3 * v3);

    // Calculate pressure using 3D energy equation
    Pressure = (GAMMA - 1.0) * (U_Vect[ENERGY_INDEX] - 0.5 * Rho * vmag);

    // Calculate temperature
    if (Non_Dimensional_Form)
        Temperature = (Pressure / (Rho * R_ref)); // Non-dimensionalized EOS
    else
        Temperature = (Pressure * inv_Density / R_GAS);

    if (Pressure < 0.0)
    {
        cout << "Negative pressure detected at Cell " << Cell_No
             << ", P = " << Pressure << endl;
    }

    // Calculate speed of sound and Mach number
    C = sqrt(GAMMA * Pressure * inv_Density);
    M = sqrt(vmag) / C;

    // Calculate viscous properties if needed
    if (Is_Viscous_Wall)
    {
        Viscosity(Temperature);
        K = mu_star * C_P / PR;
    }
    else
    {
        K = 0.0;
        mu_star = 0.0;
    }

    // Store primitive variables (extended for 3D)
    Global_Primitive[PRIM_RHO] = Rho;           // Density
    Global_Primitive[PRIM_U] = v1;              // u-velocity
    Global_Primitive[PRIM_V] = v2;              // v-velocity
    Global_Primitive[PRIM_W] = v3;              // w-velocity (3D)
    Global_Primitive[PRIM_P] = Pressure;        // Pressure
    Global_Primitive[5] = Temperature;          // Temperature
    Global_Primitive[6] = C;                    // Speed of sound
    Global_Primitive[7] = U_Vect[ENERGY_INDEX]; // Total energy
    Global_Primitive[8] = M;                    // Mach number
    Global_Primitive[9] = mu_star;              // Non-dimensional viscosity
    Global_Primitive[10] = K;                   // Thermal conductivity
    Global_Primitive[11] = Pressure * pow((1.0 + 0.5 * (GAMMA - 1.0) * M * M),
                                          GAMMA / (GAMMA - 1.0)); // Total pressure

    if (isnan(C) || isnan(Pressure))
    {
        cout << "Error: Invalid thermodynamic state at Cell " << Cell_No << endl;
        cout << "Rho = " << Rho << ", P = " << Pressure << ", C = " << C << endl;
        exit(0);
    }
}

// Overloaded function to calculate primitive variables and store them in a given vector (3D)
// Cell_No: Cell index
// U_Vect: Vector of conservative variables
// G_Primitive: Vector to store primitive variables
void Calculate_Primitive_Variables(const int &Cell_No, V_D &U_Vect, V_D &G_Primitive)
{
    double vmag = 0.0, inv_Density = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0;
    double Pressure = 0.0, Temperature = 0.0, C = 0.0, Rho = 0.0;
    mu_star = 0.0;
    K = 0.0;

    Vector_Reset(G_Primitive);

    Rho = U_Vect[RHO_INDEX];
    inv_Density = 1.0 / Rho;

    // Extract velocity components (3D)
    v1 = U_Vect[RHU_INDEX] * inv_Density;
    v2 = U_Vect[RHV_INDEX] * inv_Density;
    v3 = U_Vect[RHW_INDEX] * inv_Density;

    // Apply small velocity cutoff
    if (fabs(v1) <= EPSILON)
        v1 = 0.0;
    if (fabs(v2) <= EPSILON)
        v2 = 0.0;
    if (fabs(v3) <= EPSILON)
        v3 = 0.0;

    // Calculate velocity magnitude (3D)
    vmag = (v1 * v1 + v2 * v2 + v3 * v3);

    // Calculate pressure
    Pressure = (GAMMA - 1.0) * (U_Vect[ENERGY_INDEX] - 0.5 * Rho * vmag);

    // Calculate temperature
    if (Non_Dimensional_Form)
        Temperature = (Pressure / (Rho * R_ref));
    else
        Temperature = (Pressure * inv_Density / R_GAS);

    // Calculate speed of sound and Mach number
    C = sqrt(GAMMA * Pressure * inv_Density);
    M = sqrt(vmag) / C;

    // Calculate viscous properties
    if (Is_Viscous_Wall)
    {
        Viscosity(Temperature);
        K = mu_star * C_P / PR;
    }
    else
    {
        K = 0.0;
        mu_star = 0.0;
    }

    // Store primitive variables (3D)
    G_Primitive[PRIM_RHO] = Rho;
    G_Primitive[PRIM_U] = v1;
    G_Primitive[PRIM_V] = v2;
    G_Primitive[PRIM_W] = v3; // 3D component
    G_Primitive[PRIM_P] = Pressure;
    G_Primitive[5] = Temperature;
    G_Primitive[6] = C;
    G_Primitive[7] = U_Vect[ENERGY_INDEX];
    G_Primitive[8] = M;
    G_Primitive[9] = mu_star;
    G_Primitive[10] = K;
    G_Primitive[11] = Pressure * pow((1.0 + 0.5 * (GAMMA - 1.0) * M * M),
                                     GAMMA / (GAMMA - 1.0));

    if (isnan(C) || isnan(Pressure))
    {
        cout << "Error: Invalid thermodynamic state at Cell " << Cell_No << endl;
        exit(0);
    }
}

// Function to calculate conservative variables from primitive variables (3D)
// Primitive_Vect: Vector of primitive variables
void Calculate_Computational_Variables(V_D &Primitive_Vect)
{
    V_D::iterator P_V_iter = Primitive_Vect.begin();
    Global_U[RHO_INDEX] = P_V_iter[PRIM_RHO];                    // Density
    Global_U[RHU_INDEX] = P_V_iter[PRIM_RHO] * P_V_iter[PRIM_U]; // Rho * u
    Global_U[RHV_INDEX] = P_V_iter[PRIM_RHO] * P_V_iter[PRIM_V]; // Rho * v
    Global_U[RHW_INDEX] = P_V_iter[PRIM_RHO] * P_V_iter[PRIM_W]; // Rho * w (3D)
    Global_U[ENERGY_INDEX] = P_V_iter[7];                        // Total energy
}

// Overloaded function to calculate conservative variables for a specific cell (3D)
// Cell_No: Cell index
// U: Vector to store conservative variables
void Calculate_Computational_Variables(int &Cell_No, V_D &U)
{
    Vector_Reset(U);
    U[RHO_INDEX] = Primitive_Cells[Cell_No][PRIM_RHO];                                    // Density
    U[RHU_INDEX] = Primitive_Cells[Cell_No][PRIM_RHO] * Primitive_Cells[Cell_No][PRIM_U]; // Rho * u
    U[RHV_INDEX] = Primitive_Cells[Cell_No][PRIM_RHO] * Primitive_Cells[Cell_No][PRIM_V]; // Rho * v
    U[RHW_INDEX] = Primitive_Cells[Cell_No][PRIM_RHO] * Primitive_Cells[Cell_No][PRIM_W]; // Rho * w (3D)
    U[ENERGY_INDEX] = Primitive_Cells[Cell_No][7];                                        // Total energy
}

// Function to calculate conservative variables from given inputs (3D version)
// var: Either pressure or density
// V: Velocity vector [u, v, w]
// var1: Either temperature or density
// i: Flag to identify the type of input (1: Pressure & Temperature, 2: Pressure & Density, default: Density & Temperature)
void Calculate_Computational_Variables(const double &var, const V_D &V, const double &var1, const int &i)
{
    double temp_rho = 0.0, vmag = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0, T = 0.0, P = 0.0;

    // Reset conservative variables
    Global_U[RHO_INDEX] = 0.0;
    Global_U[RHU_INDEX] = 0.0;
    Global_U[RHV_INDEX] = 0.0;
    Global_U[RHW_INDEX] = 0.0;
    Global_U[ENERGY_INDEX] = 0.0;

    switch (i)
    {
    case 1:
        P = var;
        T = var1;
        temp_rho = P / (R_GAS * T); // If Pressure and Temperature are given
        break;
    case 2:
        P = var;
        temp_rho = var1; // If Pressure and Density are given
        T = P / (temp_rho * R_GAS);
        break;
    default:
        temp_rho = var;
        T = var1;
        P = temp_rho * R_GAS * T; // If Density and Temperature are given
        break;
    }

    // Extract velocity components (3D)
    v1 = V[X_DIR];
    v2 = V[Y_DIR];
    v3 = V[Z_DIR];
    vmag = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);

    // Calculate conservative variables (3D)
    Global_U[RHO_INDEX] = temp_rho;
    Global_U[RHU_INDEX] = temp_rho * v1;
    Global_U[RHV_INDEX] = temp_rho * v2;
    Global_U[RHW_INDEX] = temp_rho * v3; // 3D component
    Global_U[ENERGY_INDEX] = (P / (GAMMA - 1.0)) + temp_rho * vmag;
}

// 3D specific utility functions

// Function to calculate velocity magnitude in 3D
double Calculate_Velocity_Magnitude_3D(const V_D &U)
{
    double rho_inv = 1.0 / U[RHO_INDEX];
    double u = U[RHU_INDEX] * rho_inv;
    double v = U[RHV_INDEX] * rho_inv;
    double w = U[RHW_INDEX] * rho_inv;
    return sqrt(u * u + v * v + w * w);
}

// Function to calculate kinetic energy in 3D
double Calculate_Kinetic_Energy_3D(const V_D &U)
{
    double rho_inv = 1.0 / U[RHO_INDEX];
    double u = U[RHU_INDEX] * rho_inv;
    double v = U[RHV_INDEX] * rho_inv;
    double w = U[RHW_INDEX] * rho_inv;
    return 0.5 * U[RHO_INDEX] * (u * u + v * v + w * w);
}

// Function to calculate internal energy in 3D
double Calculate_Internal_Energy_3D(const V_D &U)
{
    return U[ENERGY_INDEX] - Calculate_Kinetic_Energy_3D(U);
}

// Function to calculate pressure from conservative variables in 3D
double Calculate_Pressure_3D(const V_D &U)
{
    double internal_energy = Calculate_Internal_Energy_3D(U);
    return (GAMMA - 1.0) * internal_energy;
}

// Function to calculate temperature from conservative variables in 3D
double Calculate_Temperature_3D(const V_D &U)
{
    double pressure = Calculate_Pressure_3D(U);
    if (Non_Dimensional_Form)
        return pressure / (U[RHO_INDEX] * R_ref);
    else
        return pressure / (U[RHO_INDEX] * R_GAS);
}

// Function to calculate speed of sound in 3D
double Calculate_Speed_of_Sound_3D(const V_D &U)
{
    double pressure = Calculate_Pressure_3D(U);
    return sqrt(GAMMA * pressure / U[RHO_INDEX]);
}

// Function to calculate Mach number in 3D
double Calculate_Mach_Number_3D(const V_D &U)
{
    double velocity_mag = Calculate_Velocity_Magnitude_3D(U);
    double speed_of_sound = Calculate_Speed_of_Sound_3D(U);
    return velocity_mag / speed_of_sound;
}

// Function to validate thermodynamic state in 3D
bool Validate_Thermodynamic_State_3D(const V_D &U, const int &Cell_No)
{
    if (U[RHO_INDEX] <= 0.0)
    {
        cout << "Error: Negative or zero density at Cell " << Cell_No
             << ", rho = " << U[RHO_INDEX] << endl;
        return false;
    }

    double pressure = Calculate_Pressure_3D(U);
    if (pressure <= 0.0)
    {
        cout << "Error: Negative or zero pressure at Cell " << Cell_No
             << ", P = " << pressure << endl;
        return false;
    }

    double temperature = Calculate_Temperature_3D(U);
    if (temperature <= 0.0)
    {
        cout << "Error: Negative or zero temperature at Cell " << Cell_No
             << ", T = " << temperature << endl;
        return false;
    }

    return true;
}

// Function to apply entropy condition in 3D
void Apply_Entropy_Condition_3D(V_D &U, const int &Cell_No)
{
    if (!Validate_Thermodynamic_State_3D(U, Cell_No))
    {
        // Apply entropy fix by adjusting internal energy
        double kinetic_energy = Calculate_Kinetic_Energy_3D(U);
        double min_internal_energy = SMALL_NUMBER * U[RHO_INDEX];
        U[ENERGY_INDEX] = max(U[ENERGY_INDEX], kinetic_energy + min_internal_energy);
    }
}