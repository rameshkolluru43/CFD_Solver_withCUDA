// File: Utilities.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef UTILITIES_H
#define UTILITIES_H

#include "definitions.h"
#include "Globals.h"

void Print(V_I &);
void Print(V_D &);
void Print(vector<V_D> &);
void Print(vector<bool> &);
void Vector_Reset(V_D &);
void Print(Cell &);
void Vector_Reset(vector<V_D> &);
void printCellEntry(const Cell &);
void Print(Cell &);

//---------------------------------------------------------------
void Maximum(double &, double &, double &, double &);
void Maximum(double &, double &, double &);
void Minimum(double &, double &, double &, double &);
void Minimum(double &, double &, double &);
#endif // #ifndef UTILITIES_H
       //--------------------------------------------------------------