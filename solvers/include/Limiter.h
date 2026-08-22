// File: Limiter.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef LIMITER_H
#define LIMITER_H

#include "definitions.h"
#include "Globals.h"

// ------------------ Limiter Functions -------------------------
/** Venkatakrishnan (1995) face limiter factor for MUSCL-2D; same formula as Venkatakrishnan_Limiter_3D. */
double Venkatakrishnan_Phi_MUSCL2D(int cell_index, int var, double delta_plus);
void Second_Order_Limiter(const int &, const int &, double &, double &, int &);
void Second_Order_Limiter(const int &, const int &, vector<double> &);
void Phi1(double &, double &, int &);
void MinMod(double &, double &, double &, double &);
void MinMod(double &, double &, double &);

/* Pair-wise slope limiters for MUSCL (same signature as two-arg MinMod: sets phi from slopes a,b). */
void VanLeerSlope2(double a, double b, double &phi);
void SuperbeeSlope2(double a, double b, double &phi);
void VanAlbadaSlope2(double a, double b, double &phi);
#endif // #ifndef LIMITER_H
       // --------------------------------------------------------------