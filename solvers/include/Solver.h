// File: Solver.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef Solver_h
#define Solver_h
#include "definitions.h"
#include "Globals.h"
/*-----------Functions required for Solver Part---------------------------*/
bool Inviscid_Solver(string &, string &);
bool Viscous_Solver(string &, string &);
void Explicit_Method();
void Runge_Kutta_Method();
void Implicit_Method();
bool runSolver();
void Lax_Fedrichs();

#endif // Solver_h
       // --------------------------------------------------------------