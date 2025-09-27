#include "definitions.h"
#include "Globals.h"
#include "Flux.h"

void Ausm_Flux(const int &Cell_No)
{
    /*
     * AUSM (Advection Upstream Splitting Method) Flux Calculation
     *
     * AUSM is a flux vector splitting scheme that splits the flux into convective
     * and pressure parts. It's particularly effective for low Mach number flows
     * and provides good shock-capturing capabilities.
     *
     * Mathematical Framework:
     * F = F_c + F_p
     * where F_c is the convective flux and F_p is the pressure flux
     *
     * This implementation computes the AUSM flux for all faces of the given cell
     * and accumulates the net flux contribution.
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

    // Process each face of the cell
    for (int face = 0; face < 4; face++)
    {
        int neighbor_idx;
        switch (face)
        {
        case 0:
            neighbor_idx = Neighbour_1;
            break;
        case 1:
            neighbor_idx = Neighbour_2;
            break;
        case 2:
            neighbor_idx = Neighbour_3;
            break;
        case 3:
            neighbor_idx = Neighbour_4;
            break;
        default:
            neighbor_idx = -1;
            break;
        }

        // Skip if this is a boundary face (no valid neighbor)
        if (neighbor_idx < 0 || neighbor_idx >= No_Physical_Cells)
        {
            continue;
        }

        // Calculate AUSM flux for this face
        double face_flux[4] = {0.0, 0.0, 0.0, 0.0};

        // Get face geometry
        double nx = Cells[Cell_No].Face_Normals[face * 2 + 0];
        double ny = Cells[Cell_No].Face_Normals[face * 2 + 1];
        double face_area = Cells[Cell_No].Face_Areas[face];

        // Left state (current cell)
        double rho_L = U_Cells[Cell_No][0];
        double rhou_L = U_Cells[Cell_No][1];
        double rhov_L = U_Cells[Cell_No][2];
        double rhoE_L = U_Cells[Cell_No][3];

        double u_L = rhou_L / rho_L;
        double v_L = rhov_L / rho_L;
        double p_L = (gamma - 1.0) * (rhoE_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
        double a_L = sqrt(gamma * p_L / rho_L);
        double H_L = (rhoE_L + p_L) / rho_L;

        // Right state (neighbor cell)
        double rho_R = U_Cells[neighbor_idx][0];
        double rhou_R = U_Cells[neighbor_idx][1];
        double rhov_R = U_Cells[neighbor_idx][2];
        double rhoE_R = U_Cells[neighbor_idx][3];

        double u_R = rhou_R / rho_R;
        double v_R = rhov_R / rho_R;
        double p_R = (gamma - 1.0) * (rhoE_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R));
        double a_R = sqrt(gamma * p_R / rho_R);
        double H_R = (rhoE_R + p_R) / rho_R;

        // Normal velocities
        double Vn_L = u_L * nx + v_L * ny;
        double Vn_R = u_R * nx + v_R * ny;

        // Interface sound speed (Roe average)
        double sqrt_rho_L = sqrt(rho_L);
        double sqrt_rho_R = sqrt(rho_R);
        double inv_sqrt_sum = 1.0 / (sqrt_rho_L + sqrt_rho_R);

        double u_roe = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sqrt_sum;
        double v_roe = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sqrt_sum;
        double H_roe = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sqrt_sum;

        double a_roe = sqrt((gamma - 1.0) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe)));

        // Mach numbers
        double M_L = Vn_L / a_roe;
        double M_R = Vn_R / a_roe;

        // AUSM+ Mach number splitting functions
        double M_plus, M_minus;

        if (fabs(M_L) >= 1.0)
        {
            M_plus = (M_L > 0.0) ? M_L : 0.0;
        }
        else
        {
            M_plus = 0.25 * (M_L + 1.0) * (M_L + 1.0);
        }

        if (fabs(M_R) >= 1.0)
        {
            M_minus = (M_R < 0.0) ? M_R : 0.0;
        }
        else
        {
            M_minus = -0.25 * (M_R - 1.0) * (M_R - 1.0);
        }

        // Interface Mach number
        double M_interface = M_plus + M_minus;

        // Pressure splitting functions
        double P_plus, P_minus;

        if (fabs(M_L) >= 1.0)
        {
            P_plus = (M_L > 0.0) ? 1.0 : 0.0;
        }
        else
        {
            P_plus = 0.25 * (M_L + 1.0) * (M_L + 1.0) * (2.0 - M_L);
        }

        if (fabs(M_R) >= 1.0)
        {
            P_minus = (M_R < 0.0) ? 1.0 : 0.0;
        }
        else
        {
            P_minus = 0.25 * (M_R - 1.0) * (M_R - 1.0) * (2.0 + M_R);
        }

        // Interface pressure
        double p_interface = P_plus * p_L + P_minus * p_R;

        // Mass flux
        double mass_flux;
        if (M_interface >= 0.0)
        {
            mass_flux = M_interface * a_roe * rho_L;
        }
        else
        {
            mass_flux = M_interface * a_roe * rho_R;
        }

        // AUSM fluxes
        if (M_interface >= 0.0)
        {
            face_flux[0] = mass_flux * face_area;
            face_flux[1] = (mass_flux * u_L + p_interface * nx) * face_area;
            face_flux[2] = (mass_flux * v_L + p_interface * ny) * face_area;
            face_flux[3] = mass_flux * H_L * face_area;
        }
        else
        {
            face_flux[0] = mass_flux * face_area;
            face_flux[1] = (mass_flux * u_R + p_interface * nx) * face_area;
            face_flux[2] = (mass_flux * v_R + p_interface * ny) * face_area;
            face_flux[3] = mass_flux * H_R * face_area;
        }

        // Accumulate flux contribution to cell net flux
        // Note: Sign convention depends on face orientation
        // For outward normals, flux out of cell is positive
        Cells_Net_Flux[Cell_No][0] += face_flux[0];
        Cells_Net_Flux[Cell_No][1] += face_flux[1];
        Cells_Net_Flux[Cell_No][2] += face_flux[2];
        Cells_Net_Flux[Cell_No][3] += face_flux[3];
    }

    // Apply boundary conditions if needed
    // Handle wall boundaries, inlet/outlet conditions, etc.
    for (int face = 0; face < 4; face++)
    {
        if (Cells_Face_Boundary_Type[Cell_No][face] != 0) // Non-internal face
        {
            // Apply appropriate boundary flux modification
            // This is a simplified approach - more sophisticated boundary
            // treatments would be implemented based on boundary type

            double face_area = Cells[Cell_No].Face_Areas[face];
            double nx = Cells[Cell_No].Face_Normals[face * 2 + 0];
            double ny = Cells[Cell_No].Face_Normals[face * 2 + 1];

            // For wall boundaries, ensure no mass flux
            if (Cells_Face_Boundary_Type[Cell_No][face] == 2) // Wall boundary
            {
                double rho = U_Cells[Cell_No][0];
                double u = U_Cells[Cell_No][1] / rho;
                double v = U_Cells[Cell_No][2] / rho;
                double p = (gamma - 1.0) * (U_Cells[Cell_No][3] - 0.5 * rho * (u * u + v * v));

                // Wall boundary: zero mass flux, only pressure contribution
                double wall_flux[4];
                wall_flux[0] = 0.0;                // No mass flux
                wall_flux[1] = p * nx * face_area; // Pressure force in x
                wall_flux[2] = p * ny * face_area; // Pressure force in y
                wall_flux[3] = 0.0;                // No energy flux

                // Replace the computed flux with wall flux
                Cells_Net_Flux[Cell_No][0] -= wall_flux[0];
                Cells_Net_Flux[Cell_No][1] -= wall_flux[1];
                Cells_Net_Flux[Cell_No][2] -= wall_flux[2];
                Cells_Net_Flux[Cell_No][3] -= wall_flux[3];
            }
        }
    }
}