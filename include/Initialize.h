// File: Initialize.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef INITIALIZE_H
#define INITIALIZE_H
#include "definitions.h"
#include "Globals.h"
// Function used for initialization of computational vectors from the exit conditions
void Initialize(const int &);
// Function used for initialization from a given file at a given time step
void Initialize(const string &);
#endif // INITIALIZE_H
       // --------------------------------------------------------------