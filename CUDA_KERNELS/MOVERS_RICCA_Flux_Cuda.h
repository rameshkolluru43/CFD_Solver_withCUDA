#ifndef MOVERS_RICCA_FLUX_CUDA_H
#define MOVERS_RICCA_FLUX_CUDA_H

/**
 * GPU net-flux evaluation for MOVERS (Dissipation_Type 2), RICCA (4), MOVERS_NWSC (5).
 * Returns true if the CUDA path ran; false to use the CPU implementation (e.g. 2nd order MUSCL).
 */
bool Evaluate_Cell_Net_Flux_CUDA_Movers_Ricca(bool second_order);

#endif
