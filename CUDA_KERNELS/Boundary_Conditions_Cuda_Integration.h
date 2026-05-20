#ifndef BOUNDARY_CONDITIONS_CUDA_INTEGRATION_H
#define BOUNDARY_CONDITIONS_CUDA_INTEGRATION_H

/** Viscous no-slip wall on ghost cells (Wall_Cells_List triplets). Returns true if CUDA ran. */
bool Apply_Viscous_Wall_Boundary_Condition_CUDA();

#endif
