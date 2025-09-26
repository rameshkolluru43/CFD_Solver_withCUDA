
#include "Flux.h"
#include "definitions.h"
#include "Globals.h"

// Function to calculate the average flux at the interface between two cells in 2D
/**
 * @brief Calculates the average convective flux at a given face of a cell.
 *
 * This function computes the average flux at the interface between two cells
 * using the left and right state variables. It also accounts for wall boundary
 * conditions if the face is a wall face. The flux is calculated based on the
 * velocity, pressure, and density of the cells, as well as the face geometry.
 *
 * @param Cell_No Index of the current cell (left cell).
 * @param N_Cell_No Index of the neighboring cell (right cell).
 * @param Face_No Index of the face for which the flux is being calculated.
 * @param Is_Wall_Face Boolean flag indicating whether the face is a wall face.
 */
void Calculate_Face_Average_Flux(const int &Cell_No, const int &N_Cell_No, const int &Face_No, bool Is_Wall_Face)
{

	// Left state variables
	const auto &Left_Cell = Primitive_Cells[Cell_No];
	double Rho_L = Left_Cell[0];
	double u_L = Left_Cell[1];
	double v_L = Left_Cell[2];
	double P_L = Left_Cell[4];

	// Right state variables
	const auto &Right_Cell = Primitive_Cells[N_Cell_No];
	double Rho_R = Right_Cell[0];
	double u_R = Right_Cell[1];
	double v_R = Right_Cell[2];
	double P_R = Right_Cell[4];

	// cout << "Cell_No\t" << Cell_No << "\tN_Cell_No\t" << N_Cell_No << "\tFace_No\t" << Face_No << "\tIs_Wall_Face\t" << Is_Wall_Face << endl;
	// cout << "Rho_L\t" << Rho_L << "\tu_L\t" << u_L << "\tv_L\t" << v_L << "\tP_L\t" << P_L << endl;
	// cout << "Rho_R\t" << Rho_R << "\tu_R\t" << u_R << "\tv_R\t" << v_R << "\tP_R\t" << P_R << endl;

	// Face properties
	const auto &Cell = Cells[Cell_No];
	double nx = Cell.Face_Normals[Face_No * 2];		// nx = dy/dl
	double ny = Cell.Face_Normals[Face_No * 2 + 1]; // ny = -dx/dl
	double dl = Cell.Face_Areas[Face_No];			// Length of the face

	// Velocity dot normal
	double Vdotn_L = (u_L * nx + v_L * ny);
	double Vdotn_R = (u_R * nx + v_R * ny);

	// Velocity magnitudes
	double Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	double Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

	// Handle wall face condition
	if (Is_Wall_Face)
	{
		P_R = P_L;
	}

	// Precompute common terms
	double Vdotn_L_dl = Vdotn_L * dl;
	double Vdotn_R_dl = Vdotn_R * dl;

	// Left state flux
	Flux_L[0] = Rho_L * Vdotn_L_dl;
	Flux_L[1] = Rho_L * u_L * Vdotn_L_dl + P_L * nx * dl;
	Flux_L[2] = Rho_L * v_L * Vdotn_L_dl + P_L * ny * dl;
	Flux_L[3] = (gamma1 * P_L + Rho_L * Vmag_L) * Vdotn_L_dl;

	// Right state flux
	Flux_R[0] = Rho_R * Vdotn_R_dl;
	Flux_R[1] = Rho_R * u_R * Vdotn_R_dl + P_R * nx * dl;
	Flux_R[2] = Rho_R * v_R * Vdotn_R_dl + P_R * ny * dl;
	Flux_R[3] = (gamma1 * P_R + Rho_R * Vmag_R) * Vdotn_R_dl;

	// Average flux
	for (int i = 0; i < NUM_FLUX_COMPONENTS; ++i)
	{
		Average_Convective_Flux[i] = 0.5 * (Flux_L[i] + Flux_R[i]);
	}
}
