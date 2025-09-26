// File: Error_Update.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef ERROR_UPDATE_H
#define ERROR_UPDATE_H
#include "definitions.h"
#include "Globals.h"


void Estimate_Error();
void Update();
void Update(int &, int &);
void Update(const int &);
void Update(const int &, V_D &);

#endif // ERROR_UPDATE_H
       // --------------------------------------------------------------