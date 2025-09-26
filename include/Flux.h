// File: Flux.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#ifndef FLUX_H
#define FLUX_H
#include "definitions.h"
#include "Globals.h"

/*------------------Functions for Flux and Dissipation Evaluations----------------*/
//---------------------Functions for Dissipation Schemes-----------------------
void Evaluate_Cell_Net_Flux_1O();
void Evaluate_Cell_Net_Flux_2O();
void Calculate_Face_Average_Flux(const int &, const int &, const int &, bool);
void Calculate_Flux_For_All_Faces(int &, void (*Dissipation_Function)(const int &, int &, const int &));

// --------------- Functions which evaluate the individual Dissipative Fluxes -------------------------
void LLF(const int &, int &, const int &);
void LLF_2O(const int &, int &, const int &);
void ROE(const int &, int &, const int &);
void ROE_2O(const int &, int &, const int &);
void MOVERS(const int &, int &, const int &);
void MOVERS_NWSC(const int &, int &, const int &);
void MOVERS_2O(const int &, int &, const int &);
void MOVERS_NWSC_2O(const int &, int &, const int &);
void RICCA(const int &, int &, const int &);
void RICCA_2O(const int &, int &, const int &);
void PVU(const int &, int &, const int &);
void Entropy_Fix(double &, double &);
void Condition_For_MOVERS(double &, double &, double &, double &, double &);
void Condition_For_MOVERS_NWSC(double &, double &, double &);
void Condition_For_RICCA(double &, double &, double &, double &, double &, double &, double &, double &);
void Van_Leer_Flux(const int &);
void Van_Leer_Flux_2O(const int &);
void AUSM_Flux(const int &);
//---------------------------------------------------------------
#endif // #ifndef FLUX_H