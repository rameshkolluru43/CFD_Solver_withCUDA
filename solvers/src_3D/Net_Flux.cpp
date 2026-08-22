/**
 * @file Net_Flux.cpp
 * @brief 3D Net flux calculation for hexahedral finite volume cells
 *
 * This module computes the net flux through all 6 faces of hexahedral cells
 * for 3D CFD simulations. It integrates convective and dissipative fluxes
 * from all faces to evaluate the total flux balance for each cell.
 *
 * Key Features:
 * - 6-face flux integration for hexahedral cells
 * - Support for multiple flux schemes (Van Leer, Roe, AUSM, LLF)
 * - First and second-order accurate methods
 * - 3D boundary condition handling
 * - Conservation enforcement across shared faces
 * - Volume-based flux integration
 *
 * Mathematical Framework:
 * - Volume integral: ∫∫∫ ∂U/∂t dV = -∮∮ F⋅n dS
 * - Discrete form: V * dU/dt = -Σ(F_i * A_i) for i=0,1,2,3,4,5
 * - Face ordering: Left(0), Right(1), Bottom(2), Top(3), Back(4), Front(5)
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"

// 3D Global flux working vectors
V_D Average_Convective_Flux_3D(NUM_CONSERVATIVE_3D, 0.0); // 5-component flux
V_D Dissipative_Flux_3D(NUM_CONSERVATIVE_3D, 0.0);        // 5-component dissipation

/**
 * @brief Calculate flux through all 6 faces of a 3D hexahedral cell
 *
 * This function computes the net flux for a hexahedral cell by integrating
 * convective and dissipative fluxes through all 6 faces. It handles both
 * internal faces and boundary faces with appropriate treatments.
 *
 * @param Current_Cell_No Current cell index
 * @param Dissipation_Function Pointer to dissipation scheme function
 *
 * Face numbering convention:
 * - Face_0: Left face   (x- direction, normal = [-1, 0, 0])
 * - Face_1: Right face  (x+ direction, normal = [+1, 0, 0])
 * - Face_2: Bottom face (y- direction, normal = [ 0,-1, 0])
 * - Face_3: Top face    (y+ direction, normal = [ 0,+1, 0])
 * - Face_4: Back face   (z- direction, normal = [ 0, 0,-1])
 * - Face_5: Front face  (z+ direction, normal = [ 0, 0,+1])
 */
void Calculate_Flux_For_All_Faces_3D(int &Current_Cell_No,
                                     void (*Dissipation_Function)(const int &, int &, const int &))
{
    // 6 neighbor cells for hexahedral connectivity
    int neighbors[NUM_FACES_3D];

    // Extract all 6 neighbors
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        neighbors[face] = Cells[Current_Cell_No].Neighbours[face];
    }

#ifdef DEBUG_3D
    cout << "3D Flux calculation for cell " << Current_Cell_No << endl;
    cout << "Neighbors: ";
    for (int i = 0; i < NUM_FACES_3D; i++)
    {
        cout << neighbors[i] << " ";
    }
    cout << endl;
#endif

    // Integrate flux through all 6 faces
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Calculate average convective flux for this face
        Calculate_Face_Average_Flux_3D(Current_Cell_No, neighbors[face], face,
                                       Cells_Face_Boundary_Type_3D[Current_Cell_No][face]);

        // Apply dissipation scheme for this face
        Dissipation_Function(Current_Cell_No, neighbors[face], face);

        // Accumulate net flux (convective - dissipative)
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            Cells_Net_Flux_3D[Current_Cell_No][comp] +=
                (Average_Convective_Flux_3D[comp] - Dissipative_Flux_3D[comp]);
        }

#ifdef DEBUG_3D
        if (face == 0)
        { // Debug first face only
            cout << "Face " << face << " flux components: ";
            for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
            {
                cout << (Average_Convective_Flux_3D[comp] - Dissipative_Flux_3D[comp]) << " ";
            }
            cout << endl;
        }
#endif
    }

    // Scale by inverse volume for finite volume formulation
    double inv_volume = Cells[Current_Cell_No].Inv_Volume;
    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        Cells_Net_Flux_3D[Current_Cell_No][comp] *= inv_volume;
    }
}

/**
 * @brief First-order accurate 3D net flux evaluation
 *
 * This function evaluates the net flux for all physical cells using
 * first-order accurate methods. It supports multiple dissipation schemes
 * and provides robust convergence for 3D simulations.
 */
void Evaluate_Cell_Net_Flux_1O_3D()
{
#ifdef DEBUG_3D
    cout << "Entered 3D 1st Order Flux Calculation" << endl;
    cout << "Processing " << No_Physical_Cells << " physical cells" << endl;
    cout << "Dissipation Type: " << Dissipation_Type << endl;
#endif

// Loop over all physical cells
#pragma omp parallel for if (No_Physical_Cells > 1000)
    for (int Current_Cell_No = 0; Current_Cell_No < No_Physical_Cells; Current_Cell_No++)
    {

        // Initialize net flux to zero for all 5 components
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            Cells_Net_Flux_3D[Current_Cell_No][comp] = 0.0;
        }

        // Select dissipation scheme and calculate flux through all 6 faces
        switch (Dissipation_Type)
        {
        case 1: // LLF (Local Lax-Friedrichs)
            Calculate_Flux_For_All_Faces_3D(Current_Cell_No, LLF_3D);
            break;

        case 2: // Van Leer flux splitting
            Calculate_Flux_For_All_Faces_3D(Current_Cell_No, Van_Leer_3D);
            break;

        case 3: // Roe approximate Riemann solver
            Calculate_Flux_For_All_Faces_3D(Current_Cell_No, Roe_3D);
            break;

        case 4: // AUSM scheme
            Calculate_Flux_For_All_Faces_3D(Current_Cell_No, AUSM_3D);
            break;

        case 5: // MOVERS scheme (3D extension)
            Calculate_Flux_For_All_Faces_3D(Current_Cell_No, MOVERS_3D);
            break;

        default:
            cout << "Error: Unsupported 3D dissipation type: " << Dissipation_Type << endl;
            exit(1);
        }

        // Evaluate time step for this cell
        Evaluate_Time_Step_3D(Current_Cell_No);

        // Check for NaN or invalid fluxes
        Check_Flux_Validity_3D(Current_Cell_No);
    }

#ifdef DEBUG_3D
    cout << "3D 1st Order Flux Calculation completed" << endl;
#endif
}

/**
 * @brief Second-order accurate 3D net flux evaluation
 *
 * This function evaluates the net flux using second-order accurate
 * reconstruction methods such as MUSCL or WENO schemes.
 */
void Evaluate_Cell_Net_Flux_2O_3D()
{
#ifdef DEBUG_3D
    cout << "Entered 3D 2nd Order Flux Calculation" << endl;
    cout << "Reconstruction scheme: " << (Is_WENO ? "WENO3D" : "MUSCL3D") << endl;
#endif

    // Pre-compute gradients for MUSCL reconstruction if needed
    if (!Is_WENO && Is_Second_Order)
    {
        Compute_3D_Gradients_All_Cells();
    }

// Loop over all physical cells
#pragma omp parallel for if (No_Physical_Cells > 1000)
    for (int Current_Cell_No = 0; Current_Cell_No < No_Physical_Cells; Current_Cell_No++)
    {

        // Initialize net flux to zero
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            Cells_Net_Flux_3D[Current_Cell_No][comp] = 0.0;
        }

        // Apply second-order reconstruction and flux calculation
        if (Is_WENO)
        {
            // WENO reconstruction for high-order accuracy
            Calculate_Flux_WENO_3D(Current_Cell_No);
        }
        else
        {
            // MUSCL reconstruction with limiters
            Calculate_Flux_MUSCL_3D(Current_Cell_No);
        }

        // Evaluate time step
        Evaluate_Time_Step_3D(Current_Cell_No);

        // Flux validity check
        Check_Flux_Validity_3D(Current_Cell_No);
    }
}

/**
 * @brief Calculate flux using MUSCL reconstruction for 3D
 *
 * This function applies MUSCL (Monotonic Upstream-centered Schemes for
 * Conservation Laws) reconstruction to achieve second-order accuracy.
 *
 * @param Cell_No Current cell index
 */
void Calculate_Flux_MUSCL_3D(const int &Cell_No)
{
    // MUSCL reconstruction for all 6 faces
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        int neighbor = Cells[Cell_No].Neighbours[face];

        // Left and right reconstructed states
        V_D U_Left_3D(NUM_CONSERVATIVE_3D), U_Right_3D(NUM_CONSERVATIVE_3D);

        // Apply MUSCL reconstruction
        MUSCL_Reconstruct_3D(Cell_No, neighbor, face, U_Left_3D, U_Right_3D);

        // Apply limiter
        Apply_3D_Limiter(Cell_No, neighbor, face, U_Left_3D, U_Right_3D);

        // Calculate flux using reconstructed states
        switch (Dissipation_Type)
        {
        case 1:
            LLF_Flux_3D(U_Left_3D, U_Right_3D, face, Cell_No);
            break;
        case 2:
            Van_Leer_Flux_3D(U_Left_3D, U_Right_3D, face, Cell_No);
            break;
        case 3:
            Roe_Flux_3D(U_Left_3D, U_Right_3D, face, Cell_No);
            break;
        case 4:
            AUSM_Flux_3D(U_Left_3D, U_Right_3D, face, Cell_No);
            break;
        }

        // Accumulate flux
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            Cells_Net_Flux_3D[Cell_No][comp] +=
                (Average_Convective_Flux_3D[comp] - Dissipative_Flux_3D[comp]);
        }
    }

    // Volume scaling
    double inv_volume = Cells[Cell_No].Inv_Volume;
    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        Cells_Net_Flux_3D[Cell_No][comp] *= inv_volume;
    }
}

/**
 * @brief Calculate flux using WENO reconstruction for 3D
 *
 * This function applies 5th-order WENO reconstruction for high-order
 * accuracy in smooth regions while maintaining robustness near discontinuities.
 *
 * @param Cell_No Current cell index
 */
void Calculate_Flux_WENO_3D(const int &Cell_No)
{
    // WENO reconstruction requires extended stencil
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Build WENO stencil for this face
        vector<int> weno_stencil;
        Build_WENO_Stencil_3D(Cell_No, face, weno_stencil);

        // Left and right WENO reconstructed states
        V_D U_Left_WENO(NUM_CONSERVATIVE_3D), U_Right_WENO(NUM_CONSERVATIVE_3D);

        // Apply WENO reconstruction
        WENO5_Reconstruct_3D(Cell_No, face, weno_stencil, U_Left_WENO, U_Right_WENO);

        // Calculate flux using WENO states
        switch (Dissipation_Type)
        {
        case 1:
            LLF_Flux_3D(U_Left_WENO, U_Right_WENO, face, Cell_No);
            break;
        case 2:
            Van_Leer_Flux_3D(U_Left_WENO, U_Right_WENO, face, Cell_No);
            break;
        case 3:
            Roe_Flux_3D(U_Left_WENO, U_Right_WENO, face, Cell_No);
            break;
        case 4:
            AUSM_Flux_3D(U_Left_WENO, U_Right_WENO, face, Cell_No);
            break;
        }

        // Accumulate flux
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            Cells_Net_Flux_3D[Cell_No][comp] +=
                (Average_Convective_Flux_3D[comp] - Dissipative_Flux_3D[comp]);
        }
    }

    // Volume scaling
    double inv_volume = Cells[Cell_No].Inv_Volume;
    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        Cells_Net_Flux_3D[Cell_No][comp] *= inv_volume;
    }
}

/**
 * @brief Check flux validity to prevent numerical issues
 *
 * This function validates computed fluxes to ensure they don't contain
 * NaN or infinite values that could destabilize the simulation.
 *
 * @param Cell_No Current cell index
 */
void Check_Flux_Validity_3D(const int &Cell_No)
{
    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        double flux = Cells_Net_Flux_3D[Cell_No][comp];

        if (isnan(flux) || isinf(flux))
        {
            cout << "ERROR: Invalid flux detected in 3D cell " << Cell_No << endl;
            cout << "Component " << comp << ": " << flux << endl;
            cout << "Cell primitive variables: ";
            for (int i = 0; i < NUM_PRIMITIVE_3D; i++)
            {
                cout << Primitive_Cells_3D[Cell_No][i] << " ";
            }
            cout << endl;
            cout << "Cell conservative variables: ";
            for (int i = 0; i < NUM_CONSERVATIVE_3D; i++)
            {
                cout << U_Cells_3D[Cell_No][i] << " ";
            }
            cout << endl;
            exit(1);
        }
    }
}

/**
 * @brief Compute 3D gradients for all cells using Green-Gauss method
 *
 * This function computes velocity and scalar gradients required for
 * second-order MUSCL reconstruction in 3D.
 */
void Compute_3D_Gradients_All_Cells()
{
#pragma omp parallel for if (No_Physical_Cells > 500)
    for (int cell = 0; cell < No_Physical_Cells; cell++)
    {
        // Compute velocity gradients
        Compute_3D_Velocity_Gradients(cell);

        // Compute scalar gradients
        Compute_3D_Scalar_Gradients(cell);
    }
}

/**
 * @brief Add viscous flux contributions to net flux
 *
 * This function adds viscous flux contributions for Navier-Stokes equations
 * by integrating viscous fluxes through all 6 faces.
 *
 * @param Cell_No Current cell index
 */
void Add_Viscous_Flux_3D(const int &Cell_No)
{
    if (!Is_Viscous_Wall)
        return;

    // Loop through all 6 faces
    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Calculate viscous flux for this face
        Calculate_Viscous_Flux_Face_3D(Cell_No, face);

        // Add to net flux
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            Cells_Net_Flux_3D[Cell_No][comp] += Cells_Viscous_Flux_3D[Cell_No][comp];
        }
    }
}

/**
 * @brief Main 3D flux evaluation function
 *
 * This function selects between first-order and second-order methods
 * based on the simulation configuration.
 */
void Evaluate_3D_Net_Flux()
{
    if (Is_Second_Order)
    {
        Evaluate_Cell_Net_Flux_2O_3D();
    }
    else
    {
        Evaluate_Cell_Net_Flux_1O_3D();
    }

    // Add viscous contributions if needed
    if (Is_Viscous_Wall)
    {
#pragma omp parallel for if (No_Physical_Cells > 500)
        for (int cell = 0; cell < No_Physical_Cells; cell++)
        {
            Add_Viscous_Flux_3D(cell);
        }
    }
}

/**
 * @brief Calculate flux balance error for conservation check
 *
 * This function computes the total flux imbalance across the domain
 * to verify conservation properties of the 3D finite volume scheme.
 *
 * @return Total flux imbalance for each conservative variable
 */
V_D Calculate_3D_Flux_Balance()
{
    V_D total_flux_imbalance(NUM_CONSERVATIVE_3D, 0.0);

    // Sum flux contributions from all internal faces
    for (int cell = 0; cell < No_Physical_Cells; cell++)
    {
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            total_flux_imbalance[comp] += Cells_Net_Flux_3D[cell][comp] * Cells[cell].Volume;
        }
    }

    return total_flux_imbalance;
}

/**
 * @brief Print 3D flux statistics for monitoring
 *
 * This function outputs flux statistics to monitor simulation progress
 * and detect potential numerical issues.
 */
void Print_3D_Flux_Statistics()
{
    // Find maximum flux components
    double max_flux[NUM_CONSERVATIVE_3D] = {0.0};
    int max_cell[NUM_CONSERVATIVE_3D] = {0};

    for (int cell = 0; cell < No_Physical_Cells; cell++)
    {
        for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
        {
            double abs_flux = fabs(Cells_Net_Flux_3D[cell][comp]);
            if (abs_flux > max_flux[comp])
            {
                max_flux[comp] = abs_flux;
                max_cell[comp] = cell;
            }
        }
    }

    cout << "3D Flux Statistics:" << endl;
    cout << "Maximum flux components:" << endl;
    const char *component_names[] = {"Mass", "X-Mom", "Y-Mom", "Z-Mom", "Energy"};

    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        cout << "  " << component_names[comp] << ": " << max_flux[comp]
             << " (Cell " << max_cell[comp] << ")" << endl;
    }

    // Conservation check
    V_D flux_balance = Calculate_3D_Flux_Balance();
    cout << "Conservation check (should be ~0):" << endl;
    for (int comp = 0; comp < NUM_CONSERVATIVE_3D; comp++)
    {
        cout << "  " << component_names[comp] << ": " << flux_balance[comp] << endl;
    }
}