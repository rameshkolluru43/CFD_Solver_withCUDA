// Project : 2D Compressible Navier-Stokes Solver
// File    : Viscous_Functions.h
// Author  : Ramesh Kolluru
// Date    : 2025-04-13
#ifndef VISCOUS_TERMS_H
#define VISCOUS_TERMS_H
#include "definitions.h"
#include "Globals.h"
/*-------Functions to evaluate viscous terms using Simple Green Gauss Theorem---------------*/

void Viscosity(double &);
void Thermal_Conductivity(double &);
void Calculate_Gradient_On_Face(const int &, const int &, const int &);
void Viscous_Flux_on_Face(const int &, const int &);
void Evaluate_Viscous_Fluxes();
void Evaluate_Wall_Skin_Friction();
void Reference_Values();
void Identify_Neighbours_For_Second_Gradients(int &);
void Calculate_Gradient(V_D &, V_D &, V_D &, V_D &, double &, V_D &);
void Calculate_Gradients_At_Cell_Centers();
void Calculate_Gradient_At_Cell_Center(int &, int &, V_D &);
void Calculate_Vertex_Average(const V_D &, const V_D &, V_D &);

#endif // #ifndef VISCOUS_TERMS_H
// --------------------------------------------------------------
