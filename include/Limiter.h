// File: Limiter.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef LIMITER_H
#define LIMITER_H

#include "definitions.h"
#include "Globals.h"

// ------------------ Limiter Functions -------------------------
void Venkatkrishnan_Limiter(const int &, const int &, double &, double &, int &);
void Second_Order_Limiter(const int &, const int &, double &, double &, int &);
void Second_Order_Limiter(const int &, const int &, vector<double> &);
void Phi1(double &, double &, int &);
void MinMod(double &, double &, double &, double &);
void MinMod(double &, double &, double &);
#endif // #ifndef LIMITER_H
       // --------------------------------------------------------------