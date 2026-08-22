#ifndef FLUX_JACOBIAN_H
#define FLUX_JACOBIAN_H

#include "definitions.h"
#include "Globals.h"

vector<V_D> Compute_Flux_Jacobian(int &, vector<V_D> &Ac, int &);
vector<V_D> ComputeGhostCell_Flux_Jacobian(int &, int &, vector<V_D> &, int &);
void CheckMatrixForErrors(vector<V_D> &);

#endif // FLUX_JACOBIAN_H