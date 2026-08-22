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

namespace
{
constexpr double kSmall = 1.0e-14;

bool valid_cell_index(const int cell_index)
{
	return cell_index >= 0 && cell_index < static_cast<int>(Cells.size()) &&
		   cell_index < static_cast<int>(U_Cells.size());
}

void center_to_neighbor_vector(const int cell_index, const int face_no, double &dx, double &dy)
{
	dx = 0.0;
	dy = 0.0;
	const Cell &cell = Cells[cell_index];
	if (3 * face_no + 1 < static_cast<int>(cell.Cell_Center_Vector.size()))
	{
		dx = cell.Cell_Center_Vector[3 * face_no + 0];
		dy = cell.Cell_Center_Vector[3 * face_no + 1];
		return;
	}

	if (face_no < static_cast<int>(cell.Neighbours.size()) && valid_cell_index(cell.Neighbours[face_no]))
	{
		const Cell &nb = Cells[cell.Neighbours[face_no]];
		if (cell.Cell_Center.size() >= 2 && nb.Cell_Center.size() >= 2)
		{
			dx = nb.Cell_Center[0] - cell.Cell_Center[0];
			dy = nb.Cell_Center[1] - cell.Cell_Center[1];
		}
	}
}

bool least_squares_gradient(const int cell_index, const int var, double &gx, double &gy)
{
	gx = 0.0;
	gy = 0.0;
	if (!valid_cell_index(cell_index))
		return false;

	double a00 = 0.0, a01 = 0.0, a11 = 0.0;
	double b0 = 0.0, b1 = 0.0;
	const Cell &cell = Cells[cell_index];
	const double u0 = U_Cells[cell_index][var];
	const int n_faces = std::min(static_cast<int>(cell.Neighbours.size()), cell.numFaces > 0 ? cell.numFaces : static_cast<int>(cell.Neighbours.size()));

	for (int f = 0; f < n_faces; ++f)
	{
		const int nb = cell.Neighbours[f];
		if (!valid_cell_index(nb))
			continue;

		double dx = 0.0, dy = 0.0;
		center_to_neighbor_vector(cell_index, f, dx, dy);
		const double r2 = dx * dx + dy * dy;
		if (r2 < kSmall)
			continue;

		const double w = 1.0 / r2;
		const double du = U_Cells[nb][var] - u0;
		a00 += w * dx * dx;
		a01 += w * dx * dy;
		a11 += w * dy * dy;
		b0 += w * dx * du;
		b1 += w * dy * du;
	}

	const double det = a00 * a11 - a01 * a01;
	if (fabs(det) < 1.0e-30)
		return false;

	gx = (b0 * a11 - b1 * a01) / det;
	gy = (a00 * b1 - a01 * b0) / det;
	return true;
}

void neighbor_bounds(const int cell_index, const int var, double &u_min, double &u_max)
{
	const double u0 = U_Cells[cell_index][var];
	u_min = u0;
	u_max = u0;

	const Cell &cell = Cells[cell_index];
	const int n_faces = std::min(static_cast<int>(cell.Neighbours.size()), cell.numFaces > 0 ? cell.numFaces : static_cast<int>(cell.Neighbours.size()));
	for (int f = 0; f < n_faces; ++f)
	{
		const int nb = cell.Neighbours[f];
		if (valid_cell_index(nb))
		{
			u_min = std::min(u_min, U_Cells[nb][var]);
			u_max = std::max(u_max, U_Cells[nb][var]);
		}
	}
}

double bounded_tvd_phi(const int cell_index, const int var, const double delta)
{
	if (fabs(delta) < kSmall)
		return 1.0;

	double u_min = 0.0, u_max = 0.0;
	neighbor_bounds(cell_index, var, u_min, u_max);
	const double u0 = U_Cells[cell_index][var];
	const double bound = (delta > 0.0) ? (u_max - u0) : (u_min - u0);
	const double r = std::max(0.0, bound / delta);
	double phi = std::min(1.0, r);

	switch (Limiter_Case)
	{
	case 1: // Superbee
		phi = std::max(std::min(2.0 * r, 1.0), std::min(r, 2.0));
		break;
	case 3: // van Leer
		phi = (r > 0.0) ? (2.0 * r / (1.0 + r)) : 0.0;
		break;
	case 4: // van Albada
		phi = (r > 0.0) ? ((r * r + r) / (r * r + 1.0)) : 0.0;
		break;
	default: // MinMod / MCD / fallback
		break;
	}

	// Preserve the local extrema bound even for compressive limiter curves.
	return std::max(0.0, std::min(phi, std::min(1.0, r)));
}

} // namespace

/**
 * Venkatakrishnan limiter factor for the unlimited change from a cell center to a face.
 */
double Venkatakrishnan_Phi_MUSCL2D(int cell_index, int var, double Delta_plus)
{
	if (fabs(Delta_plus) < kSmall)
		return 1.0;

	double u_min = 0.0, u_max = 0.0;
	neighbor_bounds(cell_index, var, u_min, u_max);
	const double u_c = U_Cells[cell_index][var];
	const double y = (Delta_plus > 0.0) ? (u_max - u_c) : (u_min - u_c);
	double h2 = Cells[cell_index].Area;
	if (h2 <= 0.0)
		h2 = Cell_Minimum_Length * Cell_Minimum_Length;
	if (h2 <= 0.0)
		h2 = 1.0;

	const double eps2 = Venkat_K * Venkat_K * Venkat_K * h2 * h2;
	const double num = y * y + eps2 + 2.0 * y * Delta_plus;
	const double den = y * y + 2.0 * Delta_plus * Delta_plus + y * Delta_plus + eps2;
	if (fabs(den) < 1.0e-30)
		return 1.0;
	return std::max(0.0, std::min(1.0, num / den));
}

/** One-sided least-squares MUSCL state at face `Face_No`. */
static void MUSCL_ReconstructOwnerSide(const int Cell_Index, const int Face_No, vector<double> &U_face)
{
	U_face.assign(4, 0.0);
	if (!valid_cell_index(Cell_Index))
		return;

	double dx = 0.0, dy = 0.0;
	center_to_neighbor_vector(Cell_Index, Face_No, dx, dy);
	const double rx = 0.5 * dx;
	const double ry = 0.5 * dy;

	for (int k = 0; k < 4; ++k)
	{
		double gx = 0.0, gy = 0.0;
		if (!least_squares_gradient(Cell_Index, k, gx, gy))
		{
			U_face[k] = U_Cells[Cell_Index][k];
			continue;
		}

		const double delta = Second_Order_Limiter_Scale * (gx * rx + gy * ry);
		double phi = (Limiter_Case == 5)
						 ? Venkatakrishnan_Phi_MUSCL2D(Cell_Index, k, delta)
						 : bounded_tvd_phi(Cell_Index, k, delta);
		phi *= std::max(0.0, std::min(1.0, Second_Order_Phi_Blend));
		U_face[k] = U_Cells[Cell_Index][k] + phi * delta;
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
 * 1. Computes a least-squares gradient from all face neighbors.
 * 2. Limits the face extrapolation with Barth-Jespersen or Venkatakrishnan limiting.
 * 3. Reconstructs left and right conservative states at the interface.
 *
 * The function uses the following helper lambdas:
 * - computeSlopes: Computes the slopes for the current variable.
 * - applyLimiter: Applies the slope limiter selected by Limiter_Case.
 * - reconstructVariables: Reconstructs the left and right states of the variable at the interface.
 *
 * @note `Limiter_Case` (JSON / Globals) selects the slope limiter; names match `output_files.cpp`:
 *       1/3/4 select bounded Superbee / van Leer / van Albada curves, 5 selects
 *       Venkatakrishnan (`Venkat_K` in LimiterCoefficients), and other values use MinMod-style limiting.
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
