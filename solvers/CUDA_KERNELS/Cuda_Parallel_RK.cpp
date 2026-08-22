#include "definitions.h"


// Host function to execute Runge-Kutta method with OpenMP and CUDA
void Runge_Kutta_Method() {
    int Step_Case = 0, Cell_Index;
    double inv_Area = 0.0;

    // Step 1 - First step of Runge-Kutta
    Step_Case = 1;
    if (Is_Viscous_Wall) Evaluate_Viscous_Fluxes();
    if (Is_WENO) Evaluate_Cell_Net_Flux_WENO();
    else {
        if (Is_Second_Order)
            Evaluate_Cell_Net_Flux_2O();
        else
            Evaluate_Cell_Net_Flux_1O();
    }

    Max_dt = *max_element(Cells_DelT.begin(), Cells_DelT.end());
    Min_dt = *min_element(Cells_DelT.begin(), Cells_DelT.end());

    // OpenMP parallel loop for stage 1
    #pragma omp parallel for private(Cell_Index, inv_Area)
    for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++) {
        inv_Area = Cells_Inv_Area[Cell_Index];
        for (int i = 0; i < 4; i++) {
            U_Cells_RK_1[Cell_Index][i] = U_Cells[Cell_Index][i] - Min_dt * inv_Area * 
                (Cells_Net_Flux[Cell_Index][i] - Cells_Viscous_Flux[Cell_Index][i]);
        }
    }

    // Apply boundary conditions, update fluxes
    Step_Case = 2;
    Apply_Boundary_Conditions();
    if (Is_Viscous_Wall) Evaluate_Viscous_Fluxes();
    if (Is_WENO) Evaluate_Cell_Net_Flux_WENO();
    else {
        if (Is_Second_Order)
            Evaluate_Cell_Net_Flux_2O();
        else
            Evaluate_Cell_Net_Flux_1O();
    }

    // OpenMP parallel loop for stage 2
    #pragma omp parallel for private(Cell_Index, inv_Area)
    for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++) {
        inv_Area = Cells_Inv_Area[Cell_Index];
        for (int i = 0; i < 4; i++) {
            U_Cells_RK_2[Cell_Index][i] = 0.25 * (U_Cells_RK_1[Cell_Index][i] + 
                3.0 * U_Cells[Cell_Index][i] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][i] - Cells_Viscous_Flux[Cell_Index][i]));
        }
    }

    // Final step (Stage 3)
    Step_Case = 3;
    Apply_Boundary_Conditions();
    if (Is_Viscous_Wall) Evaluate_Viscous_Fluxes();
    if (Is_WENO) Evaluate_Cell_Net_Flux_WENO();
    else {
        if (Is_Second_Order)
            Evaluate_Cell_Net_Flux_2O();
        else
            Evaluate_Cell_Net_Flux_1O();
    }

    #pragma omp parallel for private(Cell_Index, inv_Area)
    for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++) {
        inv_Area = Cells_Inv_Area[Cell_Index];
        for (int i = 0; i < 4; i++) {
            Cells_DelU[Cell_Index][i] = (2.0 / 3.0) * (U_Cells_RK_2[Cell_Index][i] - 
                U_Cells[Cell_Index][i] - Min_dt * inv_Area * (Cells_Net_Flux[Cell_Index][i] - Cells_Viscous_Flux[Cell_Index][i]));
        }
    }
}
