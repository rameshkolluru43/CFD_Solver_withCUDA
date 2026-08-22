// File: Time_Step.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef TIMESTEP_H
#define TIMESTEP_H
#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"

void Evaluate_Time_Step(int &);
void Viscous_Time_Step_1(int &);
void Viscous_Time_Step_3(int &);
void Viscous_Time_Step_2(int &);
void Invisicd_Time_Step(int &);

double get_Min_dt();
double get_Max_dt();

#endif // #ifndef TIMESTEP_H
       // --------------------------------------------------------------