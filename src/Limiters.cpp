#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"
#include "Limiter.h"

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
 * 2. Computes slopes using the specified slope method (MinMod with two or three arguments).
 * 3. Applies the limiter to the computed slopes.
 * 4. Reconstructs the left and right states of variables at the interface.
 *
 * The function uses the following helper lambdas:
 * - computeSlopes: Computes the slopes for the current variable.
 * - applyLimiter: Applies the MinMod limiter to the computed slopes.
 * - reconstructVariables: Reconstructs the left and right states of the variable at the interface.
 *
 * The function supports four face directions:
 * - Face_No = 0: i-1/2, j interface
 * - Face_No = 1: i, j-1/2 interface
 * - Face_No = 2: i+1/2, j interface
 * - Face_No = 3: i, j+1/2 interface
 *
 * @note The function assumes that global variables such as `U_Cells`, `Cells`, `Limiter_Zeta`, and `Limiter_Zeta1`
 *       are defined and accessible. It also assumes that the `MinMod` function is implemented elsewhere.
 */
void Second_Order_Limiter(const int &Cell_Index, const int &Face_No, vector<double> &d_U)
{
	vector<double> d_Var_L(4, 0.0), d_Var_R(4, 0.0);
	int Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0;
	double d1 = 0.0, d2 = 0.0, d3 = 0.0, phi = 0.0;
	enum SlopeMethod
	{
		MINMOD_THREE_ARGS = 1,
		MINMOD_TWO_ARGS = 2
	};
	SlopeMethod Slope_Method = MINMOD_TWO_ARGS;

	auto computeSlopes = [&](int k, double &Slope1, double &Slope2, double &Slope3)
	{
		Slope1 = Limiter_Zeta * (U_Cells[Cell_Index][k] - U_Cells[Neighbour_1][k]) / d1;
		Slope2 = Limiter_Zeta1 * (U_Cells[Neighbour_1][k] - U_Cells[Neighbour_2][k]) / d2;
		Slope3 = Limiter_Zeta1 * (U_Cells[Cell_Index][k] - U_Cells[Neighbour_2][k]) / (d1 + d2);
	};

	auto applyLimiter = [&](double Slope1, double Slope2, double Slope3, double &phi)
	{
		switch (Slope_Method)
		{
		case MINMOD_THREE_ARGS:
			MinMod(Slope1, Slope2, Slope3, phi);
			break;
		case MINMOD_TWO_ARGS:
			MinMod(Slope1, Slope2, phi);
			break;
		}
	};

	auto reconstructVariables = [&](int k, double phi)
	{
		d_Var_R[k] = U_Cells[Neighbour_1][k] + 0.5 * phi * d1;
		d_Var_L[k] = U_Cells[Cell_Index][k] - 0.5 * phi * d1;
		d_U[k] = d_Var_R[k] - d_Var_L[k];
	};

	switch (Face_No)
	{
	case 0: // i-1/2, j interface
		Neighbour_1 = Cells[Cell_Index].Neighbours[0];
		d1 = Cells[Cell_Index].Cell_Center_Distances[0];
		Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[0];
		d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[0];
		Neighbour_3 = Cells[Cell_Index].Neighbours[2];
		d3 = Cells[Cell_Index].Cell_Center_Distances[2];
		break;

	case 1: // i, j-1/2 interface
		Neighbour_1 = Cells[Cell_Index].Neighbours[1];
		d1 = Cells[Cell_Index].Cell_Center_Distances[1];
		Neighbour_2 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[1];
		d2 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[1];
		Neighbour_3 = Cells[Cell_Index].Neighbours[3];
		d3 = Cells[Cell_Index].Cell_Center_Distances[3];
		break;

	case 2: // i+1/2, j interface
		Neighbour_1 = Cells[Cell_Index].Neighbours[2];
		d1 = Cells[Cell_Index].Cell_Center_Distances[2];
		Neighbour_2 = Cells[Cell_Index].Neighbours[0];
		d2 = Cells[Cell_Index].Cell_Center_Distances[0];
		Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[2];
		d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[2];
		break;

	case 3: // i, j+1/2 interface
		Neighbour_1 = Cells[Cell_Index].Neighbours[3];
		d1 = Cells[Cell_Index].Cell_Center_Distances[3];
		Neighbour_2 = Cells[Cell_Index].Neighbours[1];
		d2 = Cells[Cell_Index].Cell_Center_Distances[1];
		Neighbour_3 = (Neighbour_1 >= No_Physical_Cells) ? Neighbour_1 : Cells[Neighbour_1].Neighbours[3];
		d3 = (Neighbour_1 >= No_Physical_Cells) ? d1 : Cells[Neighbour_1].Cell_Center_Distances[3];
		break;
	}

	for (int k = 0; k < static_cast<int>(d_Var_L.size()); k++)
	{
		double Slope1 = 0.0, Slope2 = 0.0, Slope3 = 0.0;
		computeSlopes(k, Slope1, Slope2, Slope3);
		applyLimiter(Slope1, Slope2, Slope3, phi);
		reconstructVariables(k, phi);
	}
}
