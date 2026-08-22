/**
 * @file Viscous_Functions.h
 * @brief Viscous (Navier–Stokes) helpers: Green–Gauss co-volume gradients, wall fluxes, corner diagnostics.
 *
 * Structured TXT grids use LBRT face indexing aligned with Neighbours/BC lists.
 * See docs/MESH_AND_GRID.md and docs/HALF_CYLINDER_VALIDATION.md.
 */
// Project : 2D Compressible Navier-Stokes Solver
// Author  : Ramesh Kolluru
#ifndef VISCOUS_TERMS_H
#define VISCOUS_TERMS_H
#include "definitions.h"
#include "Globals.h"
/*-------Functions to evaluate viscous terms using Simple Green Gauss Theorem---------------*/

void Viscosity(double &);
void Thermal_Conductivity(double &);
void Calculate_Gradient_On_Face(const int &, const int &, const int &);
void Viscous_Flux_on_Face(const int &, const int &);
/** @brief Assemble viscous face fluxes into Cells_Viscous_Flux (host co-volume path; CUDA integration optional). */
void Evaluate_Viscous_Fluxes();
void Evaluate_Wall_Skin_Friction();
void Evaluate_Wall_Heat_Flux();
void Reference_Values();
/**
 * @brief Build secondary (diagonal) neighbours for co-volume face gradients.
 * Corner SW/SE/NE/NW logic must not collapse to cell 0 via an invalid ghost walk.
 */
void Identify_Neighbours_For_Second_Gradients(int &);
/**
 * @brief Write multi-BC corner cell face/normal/neighbour dump (viscous diagnostics).
 * @param opfile Output path (typically under the solution directory or plots/).
 */
void Dump_Viscous_Corner_Diagnostics(const string &opfile);
void Calculate_Gradient(V_D &, V_D &, V_D &, V_D &, double &, V_D &);
void Calculate_Gradients_At_Cell_Centers();
void Calculate_Gradient_At_Cell_Center(int &, int &, V_D &);
void Calculate_Vertex_Average(const V_D &, const V_D &, V_D &);

#endif // #ifndef VISCOUS_TERMS_H
// --------------------------------------------------------------
