// File: Boundary_Conditions.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
/*----------------Functions for Boundary Conditions----------*/
#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H
#include "definitions.h"
#include "Globals.h"

// This function groups all the other boundary conditions
void Apply_Boundary_Conditions();
void Apply_Boundary_Conditions(const int &, const bool &, const int &);
void Viscous_Wall_Boundary_Condition();
void Identify_Wall_Boundary_Faces(const int &);
void InViscid_Wall_Boundary_Condition();
void Symmetry_Boundary_Condition(unsigned int &);
void Symmetry_Boundary_Condition();
void Subsonic_Inlet_Boundary_Condition();
void Supersonic_Inlet_Boundary_Condition();
void Supersonic_Inlet_Boundary_Condition(double &, double &, double &, double &, double &);
void Subsonic_Exit_Boundary_Condition(unsigned int &);
void Supersonic_Exit_Boundary_Condition(unsigned int &);
void Subsonic_Exit_Condition(ExitCondition &, V_I &);
void Subsonic_Inlet_Condition(InletCondition &, V_I &);
void Supersonic_Inlet_Condition(InletCondition &, V_I &);
void Supersonic_Exit_Condition(ExitCondition &, V_I &);
void Wall_Boundary_Condition(const double &);
void Far_Field_Boundary_Condition();
void Super_Sonic_Inflow(int &, int &, int &);
void Super_Sonic_outflow();
void Sub_Sonic_Inflow();
void Sub_Sonic_Outflow();
void Identify_Boundary_Cells(const int &, const bool &);
#endif // #ifndef BOUNDARY_CONDITIONS_H
       // --------------------------------------------------------------