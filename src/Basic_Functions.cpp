#include "definitions.h"
#include "Globals.h"
#include "Primitive_Computational.h"
#include "Viscous_Functions.h"
#include "Utilities.h"

// Function to set the delta U values for each cell
// dU: Vector containing delta U values
void Set_DelU(V_D &dU)
{
	for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
	{
		for (int i = 0; i < 4; i++)
		{
			Cells_DelU[Cell_No][i] = dU[4 * Cell_No + i];
		}
	}
}

// Function to calculate primitive variables from conservative variables
// Cell_No: Cell index
// U_Vect: Vector of conservative variables
void Calculate_Primitive_Variables(const int &Cell_No, V_D &U_Vect)
{
	double vmag = 0.0, inv_Density = 0.0, v1 = 0.0, v2 = 0.0, Pressure = 0.0, Temperature = 0.0, C = 0.0, Rho = 0.0;

	Rho = U_Vect[0];
	inv_Density = 1.0 / Rho;
	v1 = U_Vect[1] * inv_Density;
	v2 = U_Vect[2] * inv_Density;

	if (fabs(v1) <= 1e-10)
		v1 = 0.0;
	if (fabs(v2) <= 1e-10)
		v2 = 0.0;

	vmag = (v1 * v1 + v2 * v2);
	Pressure = gamma_M_1 * (U_Vect[3] - 0.5 * Rho * vmag);

	if (Non_Dimensional_Form)
		Temperature = (Pressure / (Rho * R_ref)); // Non-dimensionalized EOS
	else
		Temperature = (Pressure * inv_Density / R_GC);

	if (Pressure < 0.0)
		cout << "Negative pressure detected at Cell " << Cell_No << endl;

	C = sqrt(gamma * Pressure * inv_Density); // Speed of sound
	M = sqrt(vmag) / C;						  // Mach number

	if (Is_Viscous_Wall)
	{
		Viscosity(Temperature);
		K = mu_star * cp_ref / Pr;
	}
	else
	{
		K = 0.0;
		mu_star = 0.0;
	}

	Global_Primitive[0] = Rho;		   // Density
	Global_Primitive[1] = v1;		   // Velocity u
	Global_Primitive[2] = v2;		   // Velocity v
	Global_Primitive[3] = Temperature; // Temperature
	Global_Primitive[4] = Pressure;	   // Pressure
	Global_Primitive[5] = C;		   // Speed of sound
	Global_Primitive[6] = U_Vect[3];   // Total energy
	Global_Primitive[7] = M;		   // Mach number
	Global_Primitive[8] = mu_star;	   // Non-dimensional viscosity
	Global_Primitive[9] = K;		   // Thermal conductivity
	Global_Primitive[10] = Pressure * pow((1.0 + 0.5 * gamma_M_1 * M * M), gamma / gamma_M_1);

	if (isnan(C))
	{
		cout << "Error: Negative pressure occurred at Cell " << Cell_No << endl;
		exit(0);
	}
}

// Overloaded function to calculate primitive variables and store them in a given vector
// Cell_No: Cell index
// U_Vect: Vector of conservative variables
// G_Primitive: Vector to store primitive variables
void Calculate_Primitive_Variables(const int &Cell_No, V_D &U_Vect, V_D &G_Primitive)
{
	double vmag = 0.0, inv_Density = 0.0, v1 = 0.0, v2 = 0.0, Pressure = 0.0, Temperature = 0.0, C = 0.0, Rho = 0.0;
	mu_star = 0.0;
	K = 0.0;

	Vector_Reset(G_Primitive);

	Rho = U_Vect[0];
	inv_Density = 1.0 / Rho;
	v1 = U_Vect[1] * inv_Density;
	v2 = U_Vect[2] * inv_Density;

	if (fabs(v1) <= 1e-10)
		v1 = 0.0;
	if (fabs(v2) <= 1e-10)
		v2 = 0.0;

	vmag = (v1 * v1 + v2 * v2);
	Pressure = gamma_M_1 * (U_Vect[3] - 0.5 * Rho * vmag);

	if (Non_Dimensional_Form)
		Temperature = (Pressure / (Rho * R_ref)); // Non-dimensionalized EOS
	else
		Temperature = (Pressure * inv_Density / R_GC);

	C = sqrt(gamma * Pressure * inv_Density); // Speed of sound
	M = sqrt(vmag) / C;						  // Mach number

	if (Is_Viscous_Wall)
	{
		Viscosity(Temperature);
		K = mu_star * cp_ref / Pr;
	}
	else
	{
		K = 0.0;
		mu_star = 0.0;
	}

	G_Primitive[0] = Rho;		  // Density
	G_Primitive[1] = v1;		  // Velocity u
	G_Primitive[2] = v2;		  // Velocity v
	G_Primitive[3] = Temperature; // Temperature
	G_Primitive[4] = Pressure;	  // Pressure
	G_Primitive[5] = C;			  // Speed of sound
	G_Primitive[6] = U_Vect[3];	  // Total energy
	G_Primitive[7] = M;			  // Mach number
	G_Primitive[8] = mu_star;	  // Non-dimensional viscosity
	G_Primitive[9] = K;			  // Thermal conductivity
	G_Primitive[10] = Pressure * pow((1.0 + 0.5 * gamma_M_1 * M * M), gamma / gamma_M_1);

	if (isnan(C))
	{
		cout << "Error: Negative pressure occurred at Cell " << Cell_No << endl;
		exit(0);
	}
}

// Function to calculate conservative variables from primitive variables
// Primitive_Vect: Vector of primitive variables
void Calculate_Computational_Variables(V_D &Primitive_Vect)
{
	V_D::iterator P_V_iter = Primitive_Vect.begin();
	Global_U[0] = P_V_iter[0];				 // Density
	Global_U[1] = P_V_iter[0] * P_V_iter[1]; // Rho * u
	Global_U[2] = P_V_iter[0] * P_V_iter[2]; // Rho * v
	Global_U[3] = P_V_iter[6];				 // Total energy
}

// Overloaded function to calculate conservative variables for a specific cell
// Cell_No: Cell index
// U: Vector to store conservative variables
void Calculate_Computational_Variables(int &Cell_No, V_D &U)
{
	Vector_Reset(U);
	U[0] = Primitive_Cells[Cell_No][0];								  // Density
	U[1] = Primitive_Cells[Cell_No][0] * Primitive_Cells[Cell_No][1]; // Rho * u
	U[2] = Primitive_Cells[Cell_No][0] * Primitive_Cells[Cell_No][2]; // Rho * v
	U[3] = Primitive_Cells[Cell_No][6];								  // Total energy
}

// Function to calculate conservative variables from given inputs
// var: Either pressure or density
// V: Velocity vector
// var1: Either temperature or density
// i: Flag to identify the type of input (1: Pressure & Temperature, 2: Pressure & Density, default: Density & Temperature)
void Calculate_Computational_Variables(const double &var, const V_D &V, const double &var1, const int &i)
{
	double temp_rho = 0.0, vmag = 0.0, v1 = 0.0, v2 = 0.0, T = 0.0, P = 0.0;

	Global_U[0] = 0.0;
	Global_U[1] = 0.0;
	Global_U[2] = 0.0;
	Global_U[3] = 0.0;

	switch (i)
	{
	case 1:
		T = var1;
		temp_rho = var / (R_GC * T); // If Pressure and Temperature are given
		break;
	case 2:
		P = var;
		temp_rho = var1; // If Pressure and Density are given
		break;
	default:
		temp_rho = var;
		T = var1;
		P = temp_rho * R_GC * T; // If Density and Temperature are given
		break;
	}

	v1 = V[0];
	v2 = V[1];
	vmag = 0.5 * (v1 * v1 + v2 * v2);

	Global_U[0] = temp_rho;
	Global_U[1] = temp_rho * v1;
	Global_U[2] = temp_rho * v2;
	Global_U[3] = ((P / gamma_M_1) + temp_rho * vmag);
}
