#include "definitions.h"

// CUDA Kernel for different Runge-Kutta stages
__global__ void runge_kutta_stage_kernel(
    double* d_U_Cells, double* d_Net_Flux, double* d_Viscous_Flux, double* d_U_Cells_RK, 
    double* d_Inv_Area, double Min_dt, int No_Physical_Cells, int stage) {
    
    int Cell_Index = blockIdx.x * blockDim.x + threadIdx.x;
    if (Cell_Index >= No_Physical_Cells) return;

    double inv_Area = d_Inv_Area[Cell_Index];
    
    // Stage-specific logic
    if (stage == 1) {
        // Stage 1 updates
        for (int i = 0; i < 4; ++i) {
            d_U_Cells_RK[Cell_Index * 4 + i] = d_U_Cells[Cell_Index * 4 + i] - 
                Min_dt * inv_Area * (d_Net_Flux[Cell_Index * 4 + i] - d_Viscous_Flux[Cell_Index * 4 + i]);
        }
    } else if (stage == 2) {
        // Stage 2 updates
        for (int i = 0; i < 4; ++i) {
            d_U_Cells_RK[Cell_Index * 4 + i] = 0.25 * (d_U_Cells_RK[Cell_Index * 4 + i] + 
                3.0 * d_U_Cells[Cell_Index * 4 + i] - Min_dt * inv_Area * (d_Net_Flux[Cell_Index * 4 + i] - d_Viscous_Flux[Cell_Index * 4 + i]));
        }
    } else if (stage == 3) {
        // Stage 3 updates
        for (int i = 0; i < 4; ++i) {
            d_U_Cells_RK[Cell_Index * 4 + i] = (2.0 / 3.0) * (d_U_Cells_RK[Cell_Index * 4 + i] - 
                d_U_Cells[Cell_Index * 4 + i] - Min_dt * inv_Area * (d_Net_Flux[Cell_Index * 4 + i] - d_Viscous_Flux[Cell_Index * 4 + i]));
        }
    }
}
