// File: Gradients.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef GRADIENTS_H
#define GRADIENTS_H
#include "definitions.h"
#include "Globals.h"
void Calculate_Gradient(V_D &, V_D &, V_D &, V_D &, double &, V_D &);
void Calculate_Gradients_At_Cell_Centers();
void Calculate_Gradient_At_Cell_Center(int &, int &, V_D &);
void Calculate_Vertex_Average(const V_D &, const V_D &, V_D &);
#endif // GRADIENTS_H
       // --------------------------------------------------------------