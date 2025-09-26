#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"

/**
 * @brief Computes the Alpha parameter for the RICCA condition based on input variables.
 *
 * This function calculates the Alpha parameter, which is used in the RICCA flux computation.
 * It considers the velocity components, pressure differences, and other physical properties
 * to determine the appropriate value of Alpha.
 *
 * @param d_U Reference to the difference in conserved variables (e.g., density, momentum, energy).
 * @param d_F Reference to the difference in flux variables.
 * @param Vn_L Normal velocity on the left side of the interface.
 * @param Vn_R Normal velocity on the right side of the interface.
 * @param dP Reference to the pressure difference across the interface.
 * @param Rho_I Density at the interface.
 * @param P_I Pressure at the interface.
 * @param Alpha Reference to the Alpha parameter to be computed.
 *
 * @note If the pressure difference (dP) is below a small threshold (epsilon), it is set to zero.
 *       The function uses the maximum of the absolute values of the normal velocities (Vn_L, Vn_R)
 *       and incorporates a term proportional to the square root of the ratio of pressure to density
 *       at the interface, scaled by a constant (gamma).
 *
 * @warning Ensure that the input parameters are properly initialized before calling this function.
 *          The function assumes that gamma is defined globally or accessible within the scope.
 */
void Condition_For_RICCA(double &d_U, double &d_F, double &Vn_L, double &Vn_R, double &dP, double &Rho_I, double &P_I, double &Alpha)
{
    double epsilon = 1e-8, Max = 0.0, ita = 1.0;
    //***************************************************************************************************************
    Alpha = 0.0;
    if (dP < epsilon)
        dP = 0.0;

    Max = max(fabs(Vn_L), fabs(Vn_R));
    if ((fabs(d_F) < epsilon and fabs(d_U) < epsilon))
    {
        Alpha = 0.5 * (fabs(Vn_L) + fabs(Vn_R));
    }
    else
    {
        Alpha = Max + ita * Sign(dP) * sqrt(gamma * P_I / Rho_I);
    }
}

/**
 * @brief Computes the dissipative flux for a given cell and face using the RICCA method.
 *
 * This function calculates the dissipative flux based on the left and right state variables
 * of a cell interface. It uses the RICCA method to evaluate the flux differences and
 * applies conditions for smooth regions.
 *
 * @param Cell_No The index of the current cell.
 * @param N_Cell_No The index of the neighboring cell.
 * @param Face_No The index of the face being processed.
 */
void RICCA(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
    double P_I, dP, Rho_I, nx, ny, dl, Vdotn_L, Vdotn_R, Vmag_L, Vmag_R;
    std::fill(d_U.begin(), d_U.end(), 0.0);
    std::fill(d_F.begin(), d_F.end(), 0.0);
    std::fill(Mod_Alpha.begin(), Mod_Alpha.end(), 0.0);
    std::fill(Dissipative_Flux.begin(), Dissipative_Flux.end(), 0.0);

    int index = Face_No * 2;

    // Left state variables
    double Rho_L = Primitive_Cells[Cell_No][0];
    double P_L = Primitive_Cells[Cell_No][4];
    double u_L = Primitive_Cells[Cell_No][1];
    double v_L = Primitive_Cells[Cell_No][2];

    // Right state variables
    double Rho_R = Primitive_Cells[N_Cell_No][0];
    double P_R = Primitive_Cells[N_Cell_No][4];
    double u_R = Primitive_Cells[N_Cell_No][1];
    double v_R = Primitive_Cells[N_Cell_No][2];

    // Face geometry
    nx = Cells[Cell_No].Face_Normals[index];
    ny = Cells[Cell_No].Face_Normals[index + 1];
    dl = Cells[Cell_No].Face_Areas[Face_No];

    // Normal velocities
    Vdotn_L = u_L * nx + v_L * ny;
    Vdotn_R = u_R * nx + v_R * ny;

    // Velocity magnitudes
    Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    // Conserved variable differences
    d_U[0] = Rho_R - Rho_L;
    d_U[1] = Rho_R * u_R - Rho_L * u_L;
    d_U[2] = Rho_R * v_R - Rho_L * v_L;
    d_U[3] = ((P_R / (gamma - 1.0)) + 0.5 * Rho_R * (u_R * u_R + v_R * v_R)) -
             ((P_L / (gamma - 1.0)) + 0.5 * Rho_L * (u_L * u_L + v_L * v_L));

    // Flux differences
    d_F[0] = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
    d_F[1] = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L + P_L * nx) * dl;
    d_F[2] = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L + P_L * ny) * dl;
    d_F[3] = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R -
              (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) *
             dl;

    // Intermediate variables
    dP = fabs(P_L - P_R);
    P_I = 0.5 * (P_L + P_R);
    Rho_I = 0.5 * (Rho_L + Rho_R);

    // Condition for smooth region
    for (int i = 0; i < 4; ++i)
    {
        Condition_For_RICCA(d_U[i], d_F[i], Vdotn_L, Vdotn_R, dP, Rho_I, P_I, Mod_Alpha[i]);
        Dissipative_Flux[i] = 0.5 * Mod_Alpha[i] * dl * d_U[i];
    }
}

// This is the second order RICCA function
/**
 * @brief Computes the dissipative flux for a given cell and face using the second-order RICCA method.
 *
 * This function calculates the dissipative flux based on the left and right state variables
 * of a cell interface. It uses the second-order RICCA method to evaluate the flux differences and
 * applies conditions for smooth regions.
 *
 * @param Cell_No The index of the current cell.
 * @param N_Cell_No The index of the neighboring cell.
 * @param Face_No The index of the face being processed.
 */
void RICCA_2O(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
    double P_I, dP, Rho_I, nx, ny, dl, Vdotn_L, Vdotn_R, Vmag_L, Vmag_R;
    std::fill(d_U.begin(), d_U.end(), 0.0);
    std::fill(d_F.begin(), d_F.end(), 0.0);
    std::fill(Mod_Alpha.begin(), Mod_Alpha.end(), 0.0);
    std::fill(Dissipative_Flux.begin(), Dissipative_Flux.end(), 0.0);

    int index = Face_No * 2;

    // Left state variables
    double Rho_L = Primitive_Cells[Cell_No][0];
    double P_L = Primitive_Cells[Cell_No][4];
    double u_L = Primitive_Cells[Cell_No][1];
    double v_L = Primitive_Cells[Cell_No][2];

    // Right state variables
    double Rho_R = Primitive_Cells[N_Cell_No][0];
    double P_R = Primitive_Cells[N_Cell_No][4];
    double u_R = Primitive_Cells[N_Cell_No][1];
    double v_R = Primitive_Cells[N_Cell_No][2];

    // Face geometry
    nx = Cells[Cell_No].Face_Normals[index];
    ny = Cells[Cell_No].Face_Normals[index + 1];
    dl = Cells[Cell_No].Face_Areas[Face_No];

    // Normal velocities
    Vdotn_L = u_L * nx + v_L * ny;
    Vdotn_R = u_R * nx + v_R * ny;

    // Velocity magnitudes
    Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    // Conserved variable differences
    d_U[0] = Rho_R - Rho_L;
    d_U[1] = Rho_R * u_R - Rho_L * u_L;
    d_U[2] = Rho_R * v_R - Rho_L * v_L;
    d_U[3] = ((P_R / (gamma - 1.0)) + 0.5 * Rho_R * (u_R * u_R + v_R * v_R)) -
             ((P_L / (gamma - 1.0)) + 0.5 * Rho_L * (u_L * u_L + v_L * v_L));

    // Flux differences
    d_F[0] = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
    d_F[1] = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L + P_L * nx) * dl;
    d_F[2] = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L + P_L * ny) * dl;
    d_F[3] = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R -
              (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) *
             dl;

    // Intermediate variables
    dP = fabs(P_L - P_R);
    P_I = 0.5 * (P_L + P_R);
    Rho_I = 0.5 * (Rho_L + Rho_R);

    // Apply second-order limiter
    Second_Order_Limiter(Cell_No, Face_No, d_U);

    // Condition for smooth region
    for (int i = 0; i < 4; ++i)
    {
        Condition_For_RICCA(d_U[i], d_F[i], Vdotn_L, Vdotn_R, dP, Rho_I, P_I, Mod_Alpha[i]);
        Dissipative_Flux[i] = 0.5 * Mod_Alpha[i] * dl * d_U[i];
    }
}
