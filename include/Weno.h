// File: Weno.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef WENO_H
#define WENO_H
#include "definitions.h"
#include "Globals.h"
/*------------------------------------------------------------------------------------*/
// Weno Functions
void Calculate_Face_WENO_Flux(int &, int &, const int &, bool);
void WENO_Reconstruction_X(int &, const int &, V_D &, V_D &);
void Evaluate_Cell_Net_Flux_WENO();
void WENO_Reconstruction(double &, double &, double &, double &, double &, int &, double &);
void Get_LR(int &, int &, const int &, V_D &, V_D &);
void MatVecMul(V_D &, V_D &, V_D &);
void Get_Reconstructed_U(int &, const int &, int &, int &, int &, int &, int &, int &, V_D &);
#endif // #ifndef WENO_H
       //------------------------------------------------------------------------------------

// #endif