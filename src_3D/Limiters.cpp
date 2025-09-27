#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"
#include "Limiter.h"

/**
 * @file Limiters.cpp
 * @brief 3D slope limiters and reconstruction methods for CFD solver
 *
 * This file implements various slope limiters and reconstruction techniques
 * for 3D hexahedral cells to ensure monotonicity and prevent spurious
 * oscillations in high-resolution finite volume schemes.
 */

/**
 * @brief Computes the MinMod limiter for three input values and updates the result in the provided reference.
 *
 * The MinMod function is used in numerical methods, particularly in finite volume and finite difference
 * schemes, to compute a slope limiter. It ensures that the solution remains monotonic by selecting the
 * minimum magnitude value among the inputs if they all have the same sign. If the inputs have mixed signs,
 * the limiter is set to zero.
 *
 * @param a Reference to the first input value.
 * @param b Reference to the second input value.
 * @param c Reference to the third input value.
 * @param phi Reference to the output value, which will be updated based on the MinMod calculation.
 *
 * The function uses the following logic:
 * - If all inputs have the same sign and are non-zero, the minimum magnitude value is selected.
 * - If all inputs are negative, the maximum magnitude value is selected.
 * - If the inputs have mixed signs, the output is set to zero.
 */
void MinMod(double &a, double &b, double &c, double &phi)
{
    double a1 = fabs(a), b1 = fabs(b), c1 = fabs(c);
    if (a1 > 0.0 && b1 > 0.0 && c1 > 0.0)
    {
        Minimum(a, b, c, phi);
    }
    else if (a < 0.0 && b < 0.0 && c < 0.0)
    {
        Maximum(a, b, c, phi);
    }
    else
    {
        phi = 0.0;
    }
}

/**
 * @brief Computes the MinMod limiter for two input values.
 *
 * The MinMod limiter is a slope limiter function used in numerical methods
 * to ensure stability and prevent oscillations in solutions. It takes two
 * input values, compares their magnitudes, and determines the limited value
 * based on their signs and the smaller magnitude.
 *
 * @param a Reference to the first input value.
 * @param b Reference to the second input value.
 * @param phi Reference to the output value, which will store the result of the MinMod limiter.
 */
void MinMod(double &a, double &b, double &phi)
{
    phi = 0.0;
    double a1 = fabs(a), b1 = fabs(b);
    phi = 0.5 * (Sign(a) + Sign(b)) * min(a1, b1);
}

/**
 * @brief SuperBee limiter for enhanced resolution
 * @param a First gradient estimate
 * @param b Second gradient estimate
 * @param phi Output limited value
 */
void SuperBee(double &a, double &b, double &phi)
{
    double r = (fabs(b) > 1e-12) ? a / b : 0.0;

    if (r <= 0.0)
    {
        phi = 0.0;
    }
    else if (r <= 0.5)
    {
        phi = 2.0 * r * b;
    }
    else if (r <= 1.0)
    {
        phi = b;
    }
    else if (r <= 2.0)
    {
        phi = min(2.0 * b, r * b);
    }
    else
    {
        phi = 2.0 * b;
    }
}

/**
 * @brief Van Leer limiter for smooth solutions
 * @param a First gradient estimate
 * @param b Second gradient estimate
 * @param phi Output limited value
 */
void VanLeer(double &a, double &b, double &phi)
{
    if (a * b <= 0.0)
    {
        phi = 0.0;
    }
    else
    {
        phi = 2.0 * a * b / (a + b);
    }
}

/**
 * @brief Van Albada limiter for smooth monotonic solutions
 * @param a First gradient estimate
 * @param b Second gradient estimate
 * @param phi Output limited value
 */
void VanAlbada(double &a, double &b, double &phi)
{
    if (a * b <= 0.0)
    {
        phi = 0.0;
    }
    else
    {
        phi = (a * (b * b + Limiter_Epsilon) + b * (a * a + Limiter_Epsilon)) /
              (a * a + b * b + 2.0 * Limiter_Epsilon);
    }
}

/**
 * @brief Second-order limiter for 3D hexahedral cells
 * @param Cell_Index The index of the current cell
 * @param Face_No The face number (0-5 for hexahedral cells)
 * @param d_U Vector to store the difference between reconstructed right and left states
 *
 * This function extends the 2D limiter to 3D by handling 6 faces instead of 4:
 * Face 0: Left   (i-1/2, j, k)
 * Face 1: Right  (i+1/2, j, k)
 * Face 2: Bottom (i, j-1/2, k)
 * Face 3: Top    (i, j+1/2, k)
 * Face 4: Back   (i, j, k-1/2)
 * Face 5: Front  (i, j, k+1/2)
 */
void Second_Order_Limiter_3D(const int &Cell_Index, const int &Face_No, V_D &d_U)
{
    V_D d_Var_L(NUM_CONSERVATIVE_VARS, 0.0), d_Var_R(NUM_CONSERVATIVE_VARS, 0.0);
    int Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0;
    double d1 = 0.0, d2 = 0.0, d3 = 0.0, phi = 0.0;

    enum SlopeMethod
    {
        MINMOD_THREE_ARGS = 1,
        MINMOD_TWO_ARGS = 2,
        SUPERBEE = 3,
        VAN_LEER = 4,
        VAN_ALBADA = 5
    };
    SlopeMethod Slope_Method = MINMOD_TWO_ARGS;

    auto computeSlopes = [&](int k, double &Slope1, double &Slope2, double &Slope3)
    {
        Slope1 = Limiter_Zeta * (U_Cells[Cell_Index][k] - U_Cells[Neighbour_1][k]) / d1;
        Slope2 = Limiter_Zeta1 * (U_Cells[Neighbour_1][k] - U_Cells[Neighbour_2][k]) / d2;
        Slope3 = Limiter_Zeta1 * (U_Cells[Cell_Index][k] - U_Cells[Neighbour_2][k]) / (d1 + d2);
    };

    auto applyLimiter = [&](double Slope1, double Slope2, double Slope3, double &phi)
    {
        switch (Slope_Method)
        {
        case MINMOD_THREE_ARGS:
            MinMod(Slope1, Slope2, Slope3, phi);
            break;
        case MINMOD_TWO_ARGS:
            MinMod(Slope1, Slope2, phi);
            break;
        case SUPERBEE:
            SuperBee(Slope1, Slope2, phi);
            break;
        case VAN_LEER:
            VanLeer(Slope1, Slope2, phi);
            break;
        case VAN_ALBADA:
            VanAlbada(Slope1, Slope2, phi);
            break;
        }
    };

    auto reconstructVariables = [&](int k, double phi)
    {
        d_Var_R[k] = U_Cells[Neighbour_1][k] + 0.5 * phi * d1;
        d_Var_L[k] = U_Cells[Cell_Index][k] - 0.5 * phi * d1;
        d_U[k] = d_Var_R[k] - d_Var_L[k];
    };

    // 3D face selection for hexahedral cells
    switch (Face_No)
    {
    case Face_0: // Left face (i-1/2, j, k)
        Neighbour_1 = Cells[Cell_Index].Neighbours[Face_0];
        d1 = Cells[Cell_Index].Cell_Center_Distances[Face_0];
        Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[Face_0];
        d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[Face_0];
        Neighbour_3 = Cells[Cell_Index].Neighbours[Face_1]; // Right neighbor for additional constraint
        d3 = Cells[Cell_Index].Cell_Center_Distances[Face_1];
        break;

    case Face_1: // Right face (i+1/2, j, k)
        Neighbour_1 = Cells[Cell_Index].Neighbours[Face_1];
        d1 = Cells[Cell_Index].Cell_Center_Distances[Face_1];
        Neighbour_2 = Cells[Cell_Index].Neighbours[Face_0]; // Left neighbor
        d2 = Cells[Cell_Index].Cell_Center_Distances[Face_0];
        Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[Face_1];
        d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[Face_1];
        break;

    case Face_2: // Bottom face (i, j-1/2, k)
        Neighbour_1 = Cells[Cell_Index].Neighbours[Face_2];
        d1 = Cells[Cell_Index].Cell_Center_Distances[Face_2];
        Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[Face_2];
        d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[Face_2];
        Neighbour_3 = Cells[Cell_Index].Neighbours[Face_3]; // Top neighbor
        d3 = Cells[Cell_Index].Cell_Center_Distances[Face_3];
        break;

    case Face_3: // Top face (i, j+1/2, k)
        Neighbour_1 = Cells[Cell_Index].Neighbours[Face_3];
        d1 = Cells[Cell_Index].Cell_Center_Distances[Face_3];
        Neighbour_2 = Cells[Cell_Index].Neighbours[Face_2]; // Bottom neighbor
        d2 = Cells[Cell_Index].Cell_Center_Distances[Face_2];
        Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[Face_3];
        d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[Face_3];
        break;

    case Face_4: // Back face (i, j, k-1/2)
        Neighbour_1 = Cells[Cell_Index].Neighbours[Face_4];
        d1 = Cells[Cell_Index].Cell_Center_Distances[Face_4];
        Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[Face_4];
        d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[Face_4];
        Neighbour_3 = Cells[Cell_Index].Neighbours[Face_5]; // Front neighbor
        d3 = Cells[Cell_Index].Cell_Center_Distances[Face_5];
        break;

    case Face_5: // Front face (i, j, k+1/2)
        Neighbour_1 = Cells[Cell_Index].Neighbours[Face_5];
        d1 = Cells[Cell_Index].Cell_Center_Distances[Face_5];
        Neighbour_2 = Cells[Cell_Index].Neighbours[Face_4]; // Back neighbor
        d2 = Cells[Cell_Index].Cell_Center_Distances[Face_4];
        Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[Face_5];
        d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[Face_5];
        break;

    default:
        cout << "ERROR: Invalid face number " << Face_No << " in Second_Order_Limiter_3D" << endl;
        exit(1);
    }

    // Apply limiting to all conservative variables
    for (int k = 0; k < NUM_CONSERVATIVE_VARS; k++)
    {
        double Slope1 = 0.0, Slope2 = 0.0, Slope3 = 0.0;
        computeSlopes(k, Slope1, Slope2, Slope3);
        applyLimiter(Slope1, Slope2, Slope3, phi);
        reconstructVariables(k, phi);
    }
}

/**
 * @brief Multidimensional slope limiter for 3D unstructured grids
 * @param Cell_Index Current cell index
 * @param gradients Cell gradients for each variable
 * @param limited_gradients Output limited gradients
 */
void Multidimensional_Limiter_3D(const int &Cell_Index, const VV_D &gradients, VV_D &limited_gradients)
{
    // Initialize limited gradients
    limited_gradients.resize(NUM_CONSERVATIVE_VARS);
    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        limited_gradients[var].resize(3, 0.0);
    }

    // Get cell center and value
    V_D cell_center(3);
    cell_center[0] = Cells[Cell_Index].Cell_Center[0];
    cell_center[1] = Cells[Cell_Index].Cell_Center[1];
    cell_center[2] = Cells[Cell_Index].Cell_Center[2];

    // Apply Barth-Jespersen multidimensional limiter
    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        double phi_min = 1.0;
        double cell_value = U_Cells[Cell_Index][var];

        // Find min/max among neighbors
        double U_min = cell_value;
        double U_max = cell_value;

        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            int neighbor = Cells[Cell_Index].Neighbours[face];
            if (neighbor >= 0 && neighbor < No_Physical_Cells)
            {
                double neighbor_value = U_Cells[neighbor][var];
                U_min = min(U_min, neighbor_value);
                U_max = max(U_max, neighbor_value);
            }
        }

        // Check each face for monotonicity
        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            // Get face center
            V_D face_center(3);
            face_center[0] = Cells[Cell_Index].Face_Centers[3 * face];
            face_center[1] = Cells[Cell_Index].Face_Centers[3 * face + 1];
            face_center[2] = Cells[Cell_Index].Face_Centers[3 * face + 2];

            // Calculate distance vector
            V_D dr(3);
            dr[0] = face_center[0] - cell_center[0];
            dr[1] = face_center[1] - cell_center[1];
            dr[2] = face_center[2] - cell_center[2];

            // Compute reconstructed value
            double U_face = cell_value + DOT_PRODUCT_3D(gradients[var], dr);

            // Apply Barth-Jespersen limiting
            double phi_face = 1.0;
            if (U_face > cell_value)
            {
                if (U_max > cell_value)
                {
                    phi_face = min(1.0, (U_max - cell_value) / (U_face - cell_value));
                }
                else
                {
                    phi_face = 0.0;
                }
            }
            else if (U_face < cell_value)
            {
                if (U_min < cell_value)
                {
                    phi_face = min(1.0, (U_min - cell_value) / (U_face - cell_value));
                }
                else
                {
                    phi_face = 0.0;
                }
            }

            phi_min = min(phi_min, phi_face);
        }

        // Apply limiter to gradient
        for (int dir = 0; dir < 3; dir++)
        {
            limited_gradients[var][dir] = phi_min * gradients[var][dir];
        }
    }
}

/**
 * @brief Venkatakrishnan limiter with enhanced smoothness
 * @param Cell_Index Current cell index
 * @param gradients Cell gradients
 * @param limited_gradients Output limited gradients
 */
void Venkatakrishnan_Limiter_3D(const int &Cell_Index, const VV_D &gradients, VV_D &limited_gradients)
{
    // Initialize
    limited_gradients.resize(NUM_CONSERVATIVE_VARS);
    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        limited_gradients[var].resize(3, 0.0);
    }

    double cell_volume = 1.0 / Cells[Cell_Index].Inv_Volume;
    double K_constant = 100.0; // Venkatakrishnan constant

    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        double phi_min = 1.0;
        double cell_value = U_Cells[Cell_Index][var];

        // Find extrema among neighbors
        double U_min = cell_value;
        double U_max = cell_value;

        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            int neighbor = Cells[Cell_Index].Neighbours[face];
            if (neighbor >= 0 && neighbor < No_Physical_Cells)
            {
                double neighbor_value = U_Cells[neighbor][var];
                U_min = min(U_min, neighbor_value);
                U_max = max(U_max, neighbor_value);
            }
        }

        // Venkatakrishnan parameter
        double Delta_sq = K_constant * K_constant * K_constant * cell_volume * cell_volume;

        // Check each face
        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            V_D dr(3);
            dr[0] = Cells[Cell_Index].Face_Centers[3 * face] - Cells[Cell_Index].Cell_Center[0];
            dr[1] = Cells[Cell_Index].Face_Centers[3 * face + 1] - Cells[Cell_Index].Cell_Center[1];
            dr[2] = Cells[Cell_Index].Face_Centers[3 * face + 2] - Cells[Cell_Index].Cell_Center[2];

            double Delta_plus = DOT_PRODUCT_3D(gradients[var], dr);

            double phi_face = 1.0;
            if (fabs(Delta_plus) > 1e-12)
            {
                if (Delta_plus > 0.0)
                {
                    double y = U_max - cell_value;
                    phi_face = (y * y + Delta_sq + 2.0 * y * Delta_plus) /
                               (y * y + 2.0 * Delta_plus * Delta_plus + y * Delta_plus + Delta_sq);
                }
                else
                {
                    double y = U_min - cell_value;
                    phi_face = (y * y + Delta_sq + 2.0 * y * Delta_plus) /
                               (y * y + 2.0 * Delta_plus * Delta_plus + y * Delta_plus + Delta_sq);
                }
            }

            phi_min = min(phi_min, phi_face);
        }

        // Apply limiter
        for (int dir = 0; dir < 3; dir++)
        {
            limited_gradients[var][dir] = phi_min * gradients[var][dir];
        }
    }
}

/**
 * @brief Enhanced BVD (Boundary Variation Diminishing) limiter for 3D
 * @param Cell_Index Current cell index
 * @param gradients Cell gradients
 * @param limited_gradients Output limited gradients
 */
void BVD_Limiter_3D(const int &Cell_Index, const VV_D &gradients, VV_D &limited_gradients)
{
    // Initialize
    limited_gradients.resize(NUM_CONSERVATIVE_VARS);
    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        limited_gradients[var].resize(3, 0.0);
    }

    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        double phi_final = 1.0;
        double cell_value = U_Cells[Cell_Index][var];

        // Enhanced BVD approach with 3D considerations
        V_D face_variations(NUM_FACES_3D, 0.0);
        double max_variation = 0.0;

        // Calculate face variations
        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            int neighbor = Cells[Cell_Index].Neighbours[face];
            if (neighbor >= 0 && neighbor < No_Physical_Cells)
            {
                face_variations[face] = fabs(U_Cells[neighbor][var] - cell_value);
                max_variation = max(max_variation, face_variations[face]);
            }
        }

        if (max_variation > 1e-12)
        {
            // Apply BVD criterion
            for (int face = 0; face < NUM_FACES_3D; face++)
            {
                V_D dr(3);
                dr[0] = Cells[Cell_Index].Face_Centers[3 * face] - Cells[Cell_Index].Cell_Center[0];
                dr[1] = Cells[Cell_Index].Face_Centers[3 * face + 1] - Cells[Cell_Index].Cell_Center[1];
                dr[2] = Cells[Cell_Index].Face_Centers[3 * face + 2] - Cells[Cell_Index].Cell_Center[2];

                double reconstructed_variation = fabs(DOT_PRODUCT_3D(gradients[var], dr));

                if (reconstructed_variation > face_variations[face] && face_variations[face] > 1e-12)
                {
                    double phi_face = face_variations[face] / reconstructed_variation;
                    phi_final = min(phi_final, phi_face);
                }
            }
        }

        // Apply limiter
        for (int dir = 0; dir < 3; dir++)
        {
            limited_gradients[var][dir] = phi_final * gradients[var][dir];
        }
    }
}

/**
 * @brief Compute limited gradients using specified limiter type
 * @param Cell_Index Current cell index
 * @param limiter_type Type of limiter to use
 * @param gradients Input gradients
 * @param limited_gradients Output limited gradients
 */
void Apply_Gradient_Limiter_3D(const int &Cell_Index, const LimiterType &limiter_type,
                               const VV_D &gradients, VV_D &limited_gradients)
{
    switch (limiter_type)
    {
    case LIMITER_BARTH_JESPERSEN:
        Multidimensional_Limiter_3D(Cell_Index, gradients, limited_gradients);
        break;

    case LIMITER_VENKATAKRISHNAN:
        Venkatakrishnan_Limiter_3D(Cell_Index, gradients, limited_gradients);
        break;

    case LIMITER_BVD:
        BVD_Limiter_3D(Cell_Index, gradients, limited_gradients);
        break;

    case LIMITER_NONE:
        // No limiting - copy gradients directly
        limited_gradients = gradients;
        break;

    default:
        // Default to Barth-Jespersen
        Multidimensional_Limiter_3D(Cell_Index, gradients, limited_gradients);
        break;
    }
}

/**
 * @brief Reconstruct variables at face using limited gradients
 * @param Cell_Index Current cell index
 * @param Face_Index Face index
 * @param limited_gradients Limited gradients
 * @param reconstructed_vars Output reconstructed variables
 */
void Reconstruct_Variables_at_Face_3D(const int &Cell_Index, const int &Face_Index,
                                      const VV_D &limited_gradients, V_D &reconstructed_vars)
{
    reconstructed_vars.resize(NUM_CONSERVATIVE_VARS);

    // Distance from cell center to face center
    V_D dr(3);
    dr[0] = Cells[Cell_Index].Face_Centers[3 * Face_Index] - Cells[Cell_Index].Cell_Center[0];
    dr[1] = Cells[Cell_Index].Face_Centers[3 * Face_Index + 1] - Cells[Cell_Index].Cell_Center[1];
    dr[2] = Cells[Cell_Index].Face_Centers[3 * Face_Index + 2] - Cells[Cell_Index].Cell_Center[2];

    // Reconstruct each variable
    for (int var = 0; var < NUM_CONSERVATIVE_VARS; var++)
    {
        double gradient_contribution = DOT_PRODUCT_3D(limited_gradients[var], dr);
        reconstructed_vars[var] = U_Cells[Cell_Index][var] + gradient_contribution;
    }
}