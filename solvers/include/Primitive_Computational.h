#ifndef PRIMITIVE_COMPUTATIONAL_H
#define PRIMITIVE_COMPUTATIONAL_H
#include "definitions.h"
#include "Globals.h"
/*---------------------------------------------------------------------------------*/

// Function for calculating Primitive Variables
void Calculate_Primitive_Variables(const int &, V_D &);
void Calculate_Primitive_Variables(const int &, V_D &, V_D &);
void Reconstruct_Primitive_Variables(const int &, int &, const int &);
// Function for calculaiting Computational Varialbes
void Calculate_Computational_Variables(V_D &);
void Calculate_Computational_Variables(const double &, const V_D &, const double &, const int &);
void Calculate_Computational_Variables(int &, V_D &);

#endif // #ifndef PRIMITIVE_COMPUTATIONAL_H