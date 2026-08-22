// File: Assemble.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef ASSEMBLE_MATRIX_H
#define ASSEMBLE_MATRIX_H
#include "definitions.h"
#include "Globals.h"

vector<V_D> Assemble_A(vector<V_D> &, double &);
V_D Assemble_b(V_D &);
void Set_DelU(V_D &);
void Assemble_A1(double &);
vector<int> get_row_indices();
vector<int> get_col_indices();
V_D get_Values();

#endif // ASSEMBLE_MATRIX_H