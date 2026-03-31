#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"
#include "Limiter.h"
#include <algorithm>

/**
 * @brief Computes the MinMod limiter for three input values and updates the result in the provided reference.
 *
 * The MinMod function is used in numerical methods, particularly in finite volume and finite difference
 * schemes, to compute a slope limiter. It ensures that the solution remains monotonic by selecting the
 * minimum magnitude value among the inputs if they all have the same sign. If the inputs have mixed signs,
 * the limiter is set to zero.
 *
 * @param a Reference to the first input value.
 * @param b Reference to the second input value.
 * @param c Reference to the third input value.
 * @param phi Reference to the output value, which will be updated based on the MinMod calculation.
 *
 * The function uses the following logic:
 * - If all inputs have the same sign and are non-zero, the minimum magnitude value is selected.
 * - If all inputs are negative, the maximum magnitude value is selected.
 * - If the inputs have mixed signs, the output is set to zero.
 */
void MinMod(double &a, double &b, double &c, double &phi)
{
	double a1 = fabs(a), b1 = fabs(b), c1 = fabs(c);
	if (a1 > 0.0 && b1 > 0.0 && c1 > 0.0)
	{
		Minimum(a, b, c, phi);
	}
	else if (a < 0.0 && b < 0.0 && c < 0.0)
	{
		Maximum(a, b, c, phi);
	}
	else
	{
		phi = 0.0;
	}
}

/**
 * @brief Computes the MinMod limiter for two input values.
 *
 * The MinMod limiter is a slope limiter function used in numerical methods
 * to ensure stability and prevent oscillations in solutions. It takes two
 * input values, compares their magnitudes, and determines the limited value
 * based on their signs and the smaller magnitude.
 *
 * @param a Reference to the first input value.
 * @param b Reference to the second input value.
 * @param phi Reference to the output value, which will store the result of the MinMod limiter.
 *
 * @note This implementation assumes the existence of a `Sign` function that
 *       returns the sign of a number (-1 for negative, 1 for positive, and 0 for zero).
 *       It also uses the `fabs` function to compute the absolute value and the `min`
 *       function to find the smaller of two values.
 */
void MinMod(double &a, double &b, double &phi)
{
	phi = 0.0;
	double a1 = fabs(a), b1 = fabs(b);
	phi = 0.5 * (Sign(a) + Sign(b)) * min(a1, b1);
}

/* Van Leer (1974) harmonic mean: more compressive than MinMod, smooth ψ. */
void VanLeerSlope2(double a, double b, double &phi)
{
	const double ab = a * b;
	if (ab <= 0.0)
	{
		phi = 0.0;
		return;
	}
	const double den = a + b;
	if (fabs(den) < 1e-30)
		phi = 0.0;
	else
		phi = 2.0 * ab / den;
}

/* Superbee: most compressive second-order TVD limiter (steeper shocks, can add noise). */
void SuperbeeSlope2(double a, double b, double &phi)
{
	if (a * b <= 0.0)
	{
		phi = 0.0;
		return;
	}
	const double s = (a > 0.0) ? 1.0 : -1.0;
	const double aa = fabs(a), bb = fabs(b);
	phi = s * max(min(2.0 * aa, bb), min(aa, 2.0 * bb));
}

/* van Albada (1982): smooth, between MinMod and Superbee — good accuracy on smooth flows. */
void VanAlbadaSlope2(double a, double b, double &phi)
{
	if (a * b <= 0.0)
	{
		phi = 0.0;
		return;
	}
	const double a2 = a * a, b2 = b * b;
	const double den = a2 + b2;
	if (den < 1e-30)
		phi = 0.0;
	else
		phi = (a * b2 + b * a2) / den;
}

/**
 * Unlimited change from cell center to the face along the segment to Neighbour_1 (half the center-to-center jump).
 * Uses local extrema among the cell and its face neighbors (physical cells only), matching Venkatakrishnan_Limiter_3D.
 */
double Venkatakrishnan_Phi_MUSCL2D(int cell_index, int var, double Delta_plus)
{
	if (fabs(Delta_plus) < 1e-14)
		return 1.0;

	const double u_c = U_Cells[cell_index][var];
	double u_min = u_c, u_max = u_c;
	const int nf = Cells[cell_index].numFaces;
	for (int f = 0; f < nf; ++f)
	{
		const int nb = Cells[cell_index].Neighbours[f];
		if (nb >= 0 && nb < No_Physical_Cells)
		{
			const double uv = U_Cells[nb][var];
			u_min = std::min(u_min, uv);
			u_max = std::max(u_max, uv);
		}
	}

	double A = Cells[cell_index].Area;
	if (A <= 0.0)
	{
		const double L = Cell_Minimum_Length;
		A = L * L;
	}
	if (A <= 0.0)
		A = 1.0;

	const double Delta_sq = Venkat_K * Venkat_K * Venkat_K * A * A;
	double y = 0.0;
	if (Delta_plus > 0.0)
		y = u_max - u_c;
	else
		y = u_min - u_c;

	const double num = y * y + Delta_sq + 2.0 * y * Delta_plus;
	const double den = y * y + 2.0 * Delta_plus * Delta_plus + y * Delta_plus + Delta_sq;
	if (fabs(den) < 1e-30)
		return 1.0;

	double phi = num / den;
	return std::max(0.0, std::min(1.0, phi));
}

/** One-sided MUSCL state at face `Face_No`: value on the side of `Cell_Index` facing that face. */
static void MUSCL_ReconstructOwnerSide(const int Cell_Index, const int Face_No, vector<double> &U_face)
{
	U_face.assign(4, 0.0);
	int Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0;
	double d1 = 0.0, d2 = 0.0, d3 = 0.0, phi = 0.0;

	const double soScale = Second_Order_Limiter_Scale;
	auto computeSlopes = [&](int k, double &Slope1, double &Slope2, double &Slope3)
	{
		Slope1 = soScale * Limiter_Zeta * (U_Cells[Cell_Index][k] - U_Cells[Neighbour_1][k]) / d1;
		Slope2 = soScale * Limiter_Zeta1 * (U_Cells[Neighbour_1][k] - U_Cells[Neighbour_2][k]) / d2;
		Slope3 = soScale * Limiter_Zeta1 * (U_Cells[Cell_Index][k] - U_Cells[Neighbour_2][k]) / (d1 + d2);
	};

	auto applyLimiter = [&](double Slope1, double Slope2, double Slope3, double &phiOut)
	{
		switch (Limiter_Case)
		{
		case 2:
			MinMod(Slope1, Slope2, Slope3, phiOut);
			break;
		case 1:
			SuperbeeSlope2(Slope1, Slope2, phiOut);
			break;
		case 3:
			VanLeerSlope2(Slope1, Slope2, phiOut);
			break;
		case 4:
			VanAlbadaSlope2(Slope1, Slope2, phiOut);
			break;
		default:
			MinMod(Slope1, Slope2, phiOut);
			break;
		}
	};

	const int nF = Cells[Cell_Index].numFaces;
	const double recon_half_sign = -1.0;

	if (nF != 4)
	{
		if (Face_No < static_cast<int>(Cells[Cell_Index].Neighbours.size()))
			Neighbour_1 = Cells[Cell_Index].Neighbours[Face_No];
		if (Face_No < static_cast<int>(Cells[Cell_Index].Cell_Center_Distances.size()))
			d1 = Cells[Cell_Index].Cell_Center_Distances[Face_No];
		if (d1 <= 0.0)
			d1 = 1.0;
		Neighbour_2 = Neighbour_1;
		d2 = d1;
		Neighbour_3 = Neighbour_1;
		d3 = d1;
	}
	else
		switch (Face_No)
		{
		case 0:
			Neighbour_1 = Cells[Cell_Index].Neighbours[0];
			d1 = Cells[Cell_Index].Cell_Center_Distances[0];
			Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[0];
			d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[0];
			Neighbour_3 = Cells[Cell_Index].Neighbours[2];
			d3 = Cells[Cell_Index].Cell_Center_Distances[2];
			break;
		case 1:
			Neighbour_1 = Cells[Cell_Index].Neighbours[1];
			d1 = Cells[Cell_Index].Cell_Center_Distances[1];
			Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[1];
			d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[1];
			Neighbour_3 = Cells[Cell_Index].Neighbours[3];
			d3 = Cells[Cell_Index].Cell_Center_Distances[3];
			break;
		case 2:
			Neighbour_1 = Cells[Cell_Index].Neighbours[2];
			d1 = Cells[Cell_Index].Cell_Center_Distances[2];
			Neighbour_2 = Cells[Cell_Index].Neighbours[0];
			d2 = Cells[Cell_Index].Cell_Center_Distances[0];
			Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[2];
			d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[2];
			break;
		case 3:
			Neighbour_1 = Cells[Cell_Index].Neighbours[3];
			d1 = Cells[Cell_Index].Cell_Center_Distances[3];
			Neighbour_2 = Cells[Cell_Index].Neighbours[1];
			d2 = Cells[Cell_Index].Cell_Center_Distances[1];
			Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[3];
			d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[3];
			break;
		default:
			break;
		}

	for (int k = 0; k < 4; ++k)
	{
		double Slope1 = 0.0, Slope2 = 0.0, Slope3 = 0.0;
		computeSlopes(k, Slope1, Slope2, Slope3);
		if (Limiter_Case == 5)
		{
			const double d_unlim = 0.5 * (U_Cells[Neighbour_1][k] - U_Cells[Cell_Index][k]);
			phi = Venkatakrishnan_Phi_MUSCL2D(Cell_Index, k, d_unlim);
		}
		else
			applyLimiter(Slope1, Slope2, Slope3, phi);
		phi *= Second_Order_Phi_Blend;
		U_face[k] = U_Cells[Cell_Index][k] + 0.5 * recon_half_sign * phi * d1;
	}
}

// This Function reconstructs the Conservative Variable on each face for a given cell
/**
 * @brief Computes the second-order limited slopes and reconstructs variables for a given cell and face.
 *
 * This function calculates the limited slopes using the MinMod limiter and reconstructs the left and right
 * states of the variables at the cell interface. It is used in finite volume methods to ensure stability
 * and accuracy of the solution.
 *
 * @param Cell_Index The index of the current cell.
 * @param Face_No The face number (0, 1, 2, or 3) corresponding to the interface direction.
 * @param d_U A vector to store the difference between reconstructed right and left states of variables.
 *
 * The function performs the following steps:
 * 1. Identifies neighboring cells and computes distances between cell centers.
 * 2. Computes slopes and applies the limiter from Limiter_Case (two-point TVD limiters or three-point MinMod for MCD).
 * 3. Applies the limiter to the computed slopes.
 * 4. Reconstructs the left and right states of variables at the interface.
 *
 * The function uses the following helper lambdas:
 * - computeSlopes: Computes the slopes for the current variable.
 * - applyLimiter: Applies the slope limiter selected by Limiter_Case.
 * - reconstructVariables: Reconstructs the left and right states of the variable at the interface.
 *
 * Face indices match Cells[*].Neighbours and Grid_Computations (structured quads):
 * - Face_No = 0: left   (neighbor Neighbours[0], -i)
 * - Face_No = 1: bottom (neighbor Neighbours[1], -j)
 * - Face_No = 2: right  (neighbor Neighbours[2], +i)
 * - Face_No = 3: top    (neighbor Neighbours[3], +j)
 *
 * Reconstruction uses one sign for all faces given Slope1 = (U_cell - U_neighbour) / d1; see body.
 *
 * @note `Limiter_Case` (JSON / Globals) selects the slope limiter; names match `output_files.cpp`:
 *       0 MinMod (two-point), 1 Superbee, 2 MCD (three-point MinMod), 3 van Leer,
 *       4 van Albada (filenames use legacy _Log_Limiter), 5 Venkatakrishnan (`Venkat_K` in LimiterCoefficients),
 *       6 and other values fall back to MinMod.
 */
void Second_Order_Limiter(const int &Cell_Index, const int &Face_No, vector<double> &d_U)
{
	vector<double> U_L(4, 0.0), U_R(4, 0.0);
	MUSCL_ReconstructOwnerSide(Cell_Index, Face_No, U_L);

	const int nF = Cells[Cell_Index].numFaces;
	int nb = -1;
	if (Face_No < static_cast<int>(Cells[Cell_Index].Neighbours.size()))
		nb = Cells[Cell_Index].Neighbours[Face_No];

	if (nb >= 0 && nF == 4)
	{
		const int oppFace = (Face_No + 2) % 4;
		MUSCL_ReconstructOwnerSide(nb, oppFace, U_R);
	}
	else if (nb >= 0)
	{
		for (int k = 0; k < 4; ++k)
			U_R[k] = U_Cells[nb][k];
	}
	else
	{
		for (int k = 0; k < 4; ++k)
			U_R[k] = U_Cells[Cell_Index][k];
	}

	const double djBlend = Second_Order_Dissipation_Blend;
	d_U.resize(4);
	for (int k = 0; k < 4; ++k)
	{
		const double d_muscl = U_R[k] - U_L[k];
		const double d_first = (nb >= 0) ? (U_Cells[nb][k] - U_Cells[Cell_Index][k]) : 0.0;
		d_U[k] = djBlend * d_muscl + (1.0 - djBlend) * d_first;
	}

	MUSCL_Face_U_L.resize(4);
	MUSCL_Face_U_R.resize(4);
	MUSCL_d_U.resize(4);
	for (int k = 0; k < 4; ++k)
	{
		MUSCL_Face_U_L[k] = U_L[k];
		MUSCL_Face_U_R[k] = U_R[k];
		MUSCL_d_U[k] = d_U[k];
	}
}
