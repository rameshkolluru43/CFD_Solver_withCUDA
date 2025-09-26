// File: Solver.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef Solver_h
#define Solver_h
#include "definitions.h"
#include "Globals.h"
/*-----------Functions required for Solver Part---------------------------*/
void Inviscid_Solver(string &, string &);
void Viscous_Solver(string &, string &);
void Explicit_Method();
void Runge_Kutta_Method();
void Implicit_Method();
void runSolver();
void Lax_Fedrichs();

#endif // Solver_h
       // --------------------------------------------------------------