#ifndef MOVERS_RICCA_FLUX_CUDA_H
#define MOVERS_RICCA_FLUX_CUDA_H

/**
 * GPU net-flux / resident explicit step for Dissipation_Type 2/4/5 (and WENO+RICCA)
 * with Jiang–Shu WENO5 on conservatives when Is_WENO.
 */
bool Evaluate_Cell_Net_Flux_CUDA_Movers_Ricca(bool second_order);

bool Resident_GPU_Explicit_Available();
bool Resident_GPU_Explicit_Init();
bool Resident_GPU_Explicit_Step(double &min_dt_out, double err_out[4]);
bool Resident_GPU_Explicit_Download_Host();
void Resident_GPU_Explicit_Shutdown();

#endif
