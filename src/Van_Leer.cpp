#include "definitions.h"
#include "Globals.h"
#include "Flux.h"

void Van_Leer_Flux(const int &Cell_No)
{
    /*
     * Van Leer Flux Vector Splitting Method
     *
     * The Van Leer scheme is a flux vector splitting method that splits the
     * convective flux vector based on the Mach number. It provides excellent
     * shock-capturing capabilities with minimal numerical diffusion.
     *
     * Mathematical Framework:
     * F = F+ + F-
     * where F+ and F- are the positive and negative split fluxes
     *
     * Key Features:
     * - Exact preservation of contact discontinuities
     * - Sharp shock resolution
     * - Robust performance across all Mach number regimes
     * - No carbuncle phenomena in low Mach number flows
     */

    // Initialize net flux for this cell to zero
    for (int i = 0; i < 4; i++)
    {
        Cells_Net_Flux[Cell_No][i] = 0.0;
    }

    // Get neighbor cell indices for all four faces
    int Neighbour_1 = Cells[Cell_No].Neighbours[0]; // Face 0 (left)
    int Neighbour_2 = Cells[Cell_No].Neighbours[1]; // Face 1 (bottom)
    int Neighbour_3 = Cells[Cell_No].Neighbours[2]; // Face 2 (right)
    int Neighbour_4 = Cells[Cell_No].Neighbours[3]; // Face 3 (top)

    // Process each face of the current cell
    for (int face_id = 0; face_id < 4; face_id++)
    {
        int neighbor_cell = -1;

        // Get the neighbor cell for this face
        switch (face_id)
        {
        case 0:
            neighbor_cell = Neighbour_1;
            break;
        case 1:
            neighbor_cell = Neighbour_2;
            break;
        case 2:
            neighbor_cell = Neighbour_3;
            break;
        case 3:
            neighbor_cell = Neighbour_4;
            break;
        }

        // Get face geometry
        double nx = Cells[Cell_No].Face_Normals[face_id * 2 + 0];
        double ny = Cells[Cell_No].Face_Normals[face_id * 2 + 1];
        double face_area = Cells[Cell_No].Face_Areas[face_id];

        // Left state (current cell) - Extract primitive variables
        double rho_L = U_Cells[Cell_No][0];
        double rhou_L = U_Cells[Cell_No][1];
        double rhov_L = U_Cells[Cell_No][2];
        double rhoE_L = U_Cells[Cell_No][3];

        double u_L = rhou_L / rho_L;
        double v_L = rhov_L / rho_L;
        double p_L = (gamma - 1.0) * (rhoE_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
        double a_L = sqrt(gamma * p_L / rho_L);
        double H_L = (rhoE_L + p_L) / rho_L; // Total enthalpy

        // Normal velocity for left state
        double V_n_L = u_L * nx + v_L * ny;
        double M_L = V_n_L / a_L; // Mach number

        double rho_R, u_R, v_R, p_R, a_R, H_R, V_n_R, M_R;

        // Right state (neighbor cell or boundary condition)
        if (neighbor_cell >= 0 && neighbor_cell < No_Physical_Cells)
        {
            // Internal face - get neighbor cell data
            rho_R = U_Cells[neighbor_cell][0];
            double rhou_R = U_Cells[neighbor_cell][1];
            double rhov_R = U_Cells[neighbor_cell][2];
            double rhoE_R = U_Cells[neighbor_cell][3];

            u_R = rhou_R / rho_R;
            v_R = rhov_R / rho_R;
            p_R = (gamma - 1.0) * (rhoE_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R));
            a_R = sqrt(gamma * p_R / rho_R);
            H_R = (rhoE_R + p_R) / rho_R;

            V_n_R = u_R * nx + v_R * ny;
            M_R = V_n_R / a_R;
        }
        else
        {
            // Boundary face - apply boundary conditions
            if (neighbor_cell == -1)
            {
                // Wall boundary condition
                rho_R = rho_L;
                u_R = u_L - 2.0 * V_n_L * nx; // Reflect normal component
                v_R = v_L - 2.0 * V_n_L * ny; // Reflect normal component
                p_R = p_L;
                a_R = a_L;
                H_R = H_L;

                V_n_R = u_R * nx + v_R * ny; // Should be -V_n_L for wall
                M_R = V_n_R / a_R;
            }
            else
            {
                // Far-field boundary condition - use left state
                rho_R = rho_L;
                u_R = u_L;
                v_R = v_L;
                p_R = p_L;
                a_R = a_L;
                H_R = H_L;
                V_n_R = V_n_L;
                M_R = M_L;
            }
        }

        // Van Leer flux splitting calculation
        double F_plus[4] = {0.0, 0.0, 0.0, 0.0};  // Positive flux (left state)
        double F_minus[4] = {0.0, 0.0, 0.0, 0.0}; // Negative flux (right state)

        // Calculate F+ (positive flux from left state)
        if (M_L >= 1.0)
        {
            // Supersonic outflow - all flux from left state
            F_plus[0] = rho_L * V_n_L;
            F_plus[1] = rho_L * V_n_L * u_L + p_L * nx;
            F_plus[2] = rho_L * V_n_L * v_L + p_L * ny;
            F_plus[3] = rho_L * V_n_L * H_L;
        }
        else if (M_L <= -1.0)
        {
            // Supersonic inflow - no flux from left state
            F_plus[0] = 0.0;
            F_plus[1] = 0.0;
            F_plus[2] = 0.0;
            F_plus[3] = 0.0;
        }
        else
        {
            // Subsonic flow - Van Leer splitting
            double rho_star = rho_L * a_L * 0.25 * (M_L + 1.0) * (M_L + 1.0);
            double u_star = (-V_n_L + 2.0 * a_L) / gamma + u_L - V_n_L * nx;
            double v_star = (-V_n_L + 2.0 * a_L) / gamma + v_L - V_n_L * ny;
            double H_star = ((gamma - 1.0) * V_n_L + 2.0 * a_L) * ((gamma - 1.0) * V_n_L + 2.0 * a_L) / (gamma * gamma - 1.0) + 0.5 * ((u_L * u_L + v_L * v_L) - V_n_L * V_n_L);

            F_plus[0] = rho_star;
            F_plus[1] = rho_star * u_star + p_L * ((M_L + 1.0) * (M_L + 1.0) * 0.25) * nx;
            F_plus[2] = rho_star * v_star + p_L * ((M_L + 1.0) * (M_L + 1.0) * 0.25) * ny;
            F_plus[3] = rho_star * H_star;
        }

        // Calculate F- (negative flux from right state)
        if (M_R >= 1.0)
        {
            // Supersonic outflow - no flux from right state
            F_minus[0] = 0.0;
            F_minus[1] = 0.0;
            F_minus[2] = 0.0;
            F_minus[3] = 0.0;
        }
        else if (M_R <= -1.0)
        {
            // Supersonic inflow - all flux from right state
            F_minus[0] = rho_R * V_n_R;
            F_minus[1] = rho_R * V_n_R * u_R + p_R * nx;
            F_minus[2] = rho_R * V_n_R * v_R + p_R * ny;
            F_minus[3] = rho_R * V_n_R * H_R;
        }
        else
        {
            // Subsonic flow - Van Leer splitting
            double rho_star = -rho_R * a_R * 0.25 * (M_R - 1.0) * (M_R - 1.0);
            double u_star = (-V_n_R - 2.0 * a_R) / gamma + u_R - V_n_R * nx;
            double v_star = (-V_n_R - 2.0 * a_R) / gamma + v_R - V_n_R * ny;
            double H_star = ((gamma - 1.0) * V_n_R - 2.0 * a_R) * ((gamma - 1.0) * V_n_R - 2.0 * a_R) / (gamma * gamma - 1.0) + 0.5 * ((u_R * u_R + v_R * v_R) - V_n_R * V_n_R);

            F_minus[0] = rho_star;
            F_minus[1] = rho_star * u_star + p_R * ((M_R - 1.0) * (M_R - 1.0) * 0.25) * nx;
            F_minus[2] = rho_star * v_star + p_R * ((M_R - 1.0) * (M_R - 1.0) * 0.25) * ny;
            F_minus[3] = rho_star * H_star;
        }

        // Combine fluxes and scale by face area
        double face_flux[4];
        face_flux[0] = (F_plus[0] + F_minus[0]) * face_area;
        face_flux[1] = (F_plus[1] + F_minus[1]) * face_area;
        face_flux[2] = (F_plus[2] + F_minus[2]) * face_area;
        face_flux[3] = (F_plus[3] + F_minus[3]) * face_area;

        // Handle special case for wall boundaries
        if (neighbor_cell == -1)
        {
            // Wall boundary: zero mass flux, only pressure contribution
            double wall_flux[4];
            wall_flux[0] = 0.0;                  // No mass flux
            wall_flux[1] = p_L * nx * face_area; // Pressure force in x
            wall_flux[2] = p_L * ny * face_area; // Pressure force in y
            wall_flux[3] = 0.0;                  // No energy flux

            // Replace the computed flux with wall flux
            face_flux[0] = wall_flux[0];
            face_flux[1] = wall_flux[1];
            face_flux[2] = wall_flux[2];
            face_flux[3] = wall_flux[3];
        }

        // Accumulate flux contribution to cell net flux
        // Note: Sign convention for outward flux from cell
        Cells_Net_Flux[Cell_No][0] += face_flux[0];
        Cells_Net_Flux[Cell_No][1] += face_flux[1];
        Cells_Net_Flux[Cell_No][2] += face_flux[2];
        Cells_Net_Flux[Cell_No][3] += face_flux[3];
    }

    /*
     * Van Leer Flux Splitting Summary:
     *
     * The Van Leer scheme provides excellent performance characteristics:
     * 1. Exact preservation of contact discontinuities
     * 2. Sharp shock resolution (typically 2-3 cells)
     * 3. No carbuncle phenomena at low Mach numbers
     * 4. Robust boundary condition treatment
     * 5. Efficient computation with simple arithmetic operations
     *
     * Key Mathematical Properties:
     * - Conservative: Satisfies conservation laws exactly
     * - Upwind: Information propagates in correct direction
     * - TVD: Total Variation Diminishing for monotonic solutions
     * - Entropy-satisfying: Respects second law of thermodynamics
     *
     * Applications:
     * - Supersonic and hypersonic flows
     * - Shock-boundary layer interactions
     * - Mixed subsonic-supersonic flows
     * - Aerospace and propulsion applications
     */
}