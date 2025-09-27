#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"

void ROE(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
	/*
	 * First-Order Roe Scheme Implementation
	 *
	 * This function implements the classical first-order Roe approximate Riemann solver
	 * for computing dissipative fluxes. The Roe scheme provides excellent shock-capturing
	 * capabilities by linearizing the Euler equations around Roe-averaged states.
	 *
	 * Mathematical Framework:
	 * 1. Compute Roe averages using density-weighted averaging
	 * 2. Calculate eigenvalues and eigenvectors of the linearized flux Jacobian
	 * 3. Decompose state differences into characteristic variables
	 * 4. Apply upwind dissipation based on wave speeds
	 *
	 * Features:
	 * - Exact resolution of contact discontinuities in 1D
	 * - Sharp shock capturing without spurious oscillations
	 * - Entropy fix for sonic points to prevent rarefaction shocks
	 * - Robust performance across all Mach number regimes
	 */

	// Variables for eigenvalues (wave speeds)
	double Lambda1 = 0.0, Lambda2 = 0.0, Lambda3 = 0.0, Lambda4 = 0.0;
	double Term1 = 0.0, Term2 = 0.0;

	// Right eigenvector components (4 eigenvectors × 4 components)
	double R11 = 0.0, R12 = 0.0, R13 = 0.0, R14 = 0.0; // Eigenvector 1: λ = u - a
	double R21 = 0.0, R22 = 0.0, R23 = 0.0, R24 = 0.0; // Eigenvector 2: λ = u
	double R31 = 0.0, R32 = 0.0, R33 = 0.0, R34 = 0.0; // Eigenvector 3: λ = u
	double R41 = 0.0, R42 = 0.0, R43 = 0.0, R44 = 0.0; // Eigenvector 4: λ = u + a

	// Roe averaging variables
	double sqrt_Rho_L = 0.0, sqrt_Rho_R = 0.0;

	// Roe-averaged quantities
	double Roe_u = 0.0, Roe_v = 0.0, Roe_Rho = 0.0, Roe_H = 0.0, Roe_Vmag = 0.0;
	double Un = 0.0, Ut = 0.0, Roe_a = 0.0;

	// Left and right state primitive variables
	double Rho_L = 0.0, Rho_R = 0.0, u_L = 0.0, u_R = 0.0, v_L = 0.0, v_R = 0.0;
	double P_L = 0.0, P_R = 0.0, H_L = 0.0, H_R = 0.0, Vmag_L = 0.0, Vmag_R = 0.0;
	double C_L = 0.0, C_R = 0.0;

	// Face geometry
	double nx = 0.0, ny = 0.0, dl = 0.0;

	// Wave strength coefficients (characteristic variables)
	double alpha_1 = 0.0, alpha_2 = 0.0, alpha_3 = 0.0, alpha_4 = 0.0;

	// State differences for characteristic decomposition
	double du = 0.0, dv = 0.0, dP = 0.0, dRho = 0.0, dUn = 0.0, dUt = 0.0;

	// Initialize dissipative flux to zero
	Dissipative_Flux[0] = 0.0;
	Dissipative_Flux[1] = 0.0;
	Dissipative_Flux[2] = 0.0;
	Dissipative_Flux[3] = 0.0;

	// Extract face geometry (normal vectors and area)
	nx = Cells[Cell_No].Face_Normals[Face_No * 2 + 0];
	ny = Cells[Cell_No].Face_Normals[Face_No * 2 + 1];
	dl = Cells[Cell_No].Face_Areas[Face_No];

	// Validate face geometry
	double normal_magnitude = sqrt(nx * nx + ny * ny);
	if (normal_magnitude < 1e-12)
	{
		// Degenerate face - skip flux calculation
		return;
	}

	// Check for boundary faces and handle appropriately
	if (N_Cell_No < 0 || N_Cell_No >= No_Physical_Cells)
	{
		// Boundary face - may need special treatment
		// For now, use the current cell state for both sides (zero gradient)
		N_Cell_No = Cell_No;
	}

	// Extract left state variables (current cell)
	Rho_L = Primitive_Cells[Cell_No][0];
	u_L = Primitive_Cells[Cell_No][1];
	v_L = Primitive_Cells[Cell_No][2];
	P_L = Primitive_Cells[Cell_No][4];
	C_L = Primitive_Cells[Cell_No][5];

	// Validate left state
	if (Rho_L <= 0.0 || P_L <= 0.0 || C_L <= 0.0)
	{
		// Invalid state - skip flux calculation or apply correction
		return;
	}

	Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	H_L = ((C_L * C_L) / (gamma - 1.0)) + Vmag_L;

	// Extract right state variables (neighbor cell)
	Rho_R = Primitive_Cells[N_Cell_No][0];
	u_R = Primitive_Cells[N_Cell_No][1];
	v_R = Primitive_Cells[N_Cell_No][2];
	P_R = Primitive_Cells[N_Cell_No][4];
	C_R = Primitive_Cells[N_Cell_No][5];

	// Validate right state
	if (Rho_R <= 0.0 || P_R <= 0.0 || C_R <= 0.0)
	{
		// Invalid state - skip flux calculation or apply correction
		return;
	}

	Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);
	H_R = ((C_R * C_R) / (gamma - 1.0)) + Vmag_R;

	// Compute Roe averages using density-weighted averaging
	sqrt_Rho_L = sqrt(Rho_L);
	sqrt_Rho_R = sqrt(Rho_R);
	Term1 = 1.0 / (sqrt_Rho_L + sqrt_Rho_R);

	// Roe-averaged density (geometric mean)
	Roe_Rho = sqrt(Rho_L * Rho_R);

	// Roe-averaged velocities (density-weighted arithmetic mean)
	Roe_u = Term1 * (u_R * sqrt_Rho_R + u_L * sqrt_Rho_L);
	Roe_v = Term1 * (v_R * sqrt_Rho_R + v_L * sqrt_Rho_L);

	// Roe-averaged total enthalpy
	Roe_H = Term1 * (H_R * sqrt_Rho_R + H_L * sqrt_Rho_L);

	// Velocity magnitude and sound speed
	Roe_Vmag = 0.5 * (Roe_u * Roe_u + Roe_v * Roe_v);
	Roe_a = sqrt((gamma - 1.0) * (Roe_H - Roe_Vmag));

	// Validate Roe-averaged sound speed
	if (Roe_a <= 0.0)
	{
		// Invalid sound speed - use fallback
		Roe_a = 0.5 * (C_L + C_R);
		if (Roe_a <= 0.0)
		{
			return; // Cannot proceed with invalid sound speed
		}
	}

	// Normal and tangential velocity components
	Un = nx * Roe_u + ny * Roe_v;  // Normal velocity
	Ut = -ny * Roe_u + nx * Roe_v; // Tangential velocity

	// Calculate state differences for characteristic decomposition
	du = (u_R - u_L);
	dv = (v_R - v_L);
	dP = (P_R - P_L);
	dRho = (Rho_R - Rho_L);

	// Normal and tangential velocity differences
	dUn = nx * du + ny * dv;  // Change in normal velocity
	dUt = -ny * du + nx * dv; // Change in tangential velocity

	// Right eigenvectors of the Roe matrix A = ∂F/∂U
	// These correspond to the four characteristic directions

	// Eigenvector 1: λ₁ = u - a (left-running acoustic wave)
	R11 = 1.0;
	R12 = Roe_u - nx * Roe_a;
	R13 = Roe_v - ny * Roe_a;
	R14 = Roe_H - Roe_a * Un;

	// Eigenvector 2: λ₂ = u (entropy wave)
	R21 = 1.0;
	R22 = Roe_u;
	R23 = Roe_v;
	R24 = Roe_Vmag;

	// Eigenvector 3: λ₃ = u (shear wave)
	R31 = 0.0;
	R32 = -ny;
	R33 = nx;
	R34 = Ut;

	// Eigenvector 4: λ₄ = u + a (right-running acoustic wave)
	R41 = 1.0;
	R42 = Roe_u + nx * Roe_a;
	R43 = Roe_v + ny * Roe_a;
	R44 = Roe_H + Roe_a * Un;

	// Wave strengths (characteristic variable coefficients)
	// These are obtained by projecting the state difference onto the left eigenvectors
	Term2 = 1.0 / (Roe_a * Roe_a);

	alpha_1 = 0.5 * Term2 * (dP - Roe_Rho * Roe_a * dUn); // Left-running acoustic wave
	alpha_2 = dRho - Term2 * dP;						  // Entropy wave
	alpha_3 = Roe_Rho * dUt;							  // Shear wave
	alpha_4 = 0.5 * Term2 * (dP + Roe_Rho * Roe_a * dUn); // Right-running acoustic wave

	// Eigenvalues (wave speeds) - these determine the direction of information propagation
	Lambda1 = fabs(Un - Roe_a); // Left-running acoustic wave speed
	Lambda2 = fabs(Un);			// Entropy wave speed
	Lambda3 = fabs(Un);			// Shear wave speed
	Lambda4 = fabs(Un + Roe_a); // Right-running acoustic wave speed

	// Apply entropy fix to prevent expansion shocks (rarefaction shocks)
	// This is crucial for physical consistency near sonic points
	double entropy_fix = 0.1 * Roe_a; // Entropy fix parameter

	if (Lambda1 < entropy_fix)
	{
		Lambda1 = 0.5 * (Lambda1 * Lambda1 / entropy_fix + entropy_fix);
	}
	if (Lambda4 < entropy_fix)
	{
		Lambda4 = 0.5 * (Lambda4 * Lambda4 / entropy_fix + entropy_fix);
	}
	// Compute dissipative flux using Roe's approximate Riemann solver
	// D = 0.5 * Σₖ |λₖ| αₖ Rₖ
	// This represents the upwind dissipation that stabilizes the central difference

	Dissipative_Flux[0] = 0.5 * (Lambda1 * alpha_1 * R11 + Lambda2 * alpha_2 * R21 + Lambda3 * alpha_3 * R31 + Lambda4 * alpha_4 * R41) * dl;

	Dissipative_Flux[1] = 0.5 * (Lambda1 * alpha_1 * R12 + Lambda2 * alpha_2 * R22 + Lambda3 * alpha_3 * R32 + Lambda4 * alpha_4 * R42) * dl;

	Dissipative_Flux[2] = 0.5 * (Lambda1 * alpha_1 * R13 + Lambda2 * alpha_2 * R23 + Lambda3 * alpha_3 * R33 + Lambda4 * alpha_4 * R43) * dl;

	Dissipative_Flux[3] = 0.5 * (Lambda1 * alpha_1 * R14 + Lambda2 * alpha_2 * R24 + Lambda3 * alpha_3 * R34 + Lambda4 * alpha_4 * R44) * dl;

	/*
	 * First-Order Roe Scheme Summary:
	 *
	 * Mathematical Foundation:
	 * - Linearizes the Euler equations around Roe-averaged states
	 * - Decomposes flux differences into characteristic waves
	 * - Applies upwind dissipation based on wave propagation direction
	 *
	 * Key Properties:
	 * - Exact resolution of contact discontinuities in 1D
	 * - Sharp shock capturing without spurious oscillations
	 * - Entropy-satisfying through entropy fix
	 * - Conservative and consistent finite volume method
	 *
	 * Wave Structure:
	 * - λ₁ = u - a: Left-running acoustic wave (compression/expansion)
	 * - λ₂ = u:     Entropy wave (density/temperature changes)
	 * - λ₃ = u:     Shear wave (tangential velocity changes)
	 * - λ₄ = u + a: Right-running acoustic wave (compression/expansion)
	 *
	 * Applications:
	 * - Shock-dominated flows with exact contact preservation
	 * - High-speed aerodynamics and gas dynamics
	 * - Foundation for higher-order schemes (MUSCL-Roe, etc.)
	 * - Robust baseline scheme for complex geometries
	 */
}
void ROE_2O(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
	/*
	 * Second-Order Roe Scheme Implementation
	 *
	 * This function implements the second-order accurate Roe scheme for computing
	 * dissipative fluxes. The key difference from the first-order scheme is the
	 * use of higher-order reconstruction through slope limiters.
	 *
	 * Mathematical Framework:
	 * 1. Compute Roe averages from left and right states
	 * 2. Calculate eigenvalues and eigenvectors of the Roe matrix
	 * 3. Apply second-order reconstruction with slope limiting
	 * 4. Compute wave strengths and dissipative flux
	 *
	 * Features:
	 * - Higher-order spatial accuracy (second-order)
	 * - TVD property through slope limiting
	 * - Robust shock capturing with minimal numerical diffusion
	 * - Entropy fix for sonic points
	 */

	// Variables required for Calculating Roe Averages
	double Lambda1 = 0.0, Lambda2 = 0.0, Lambda3 = 0.0, Lambda4 = 0.0, Term1 = 0.0, Term2 = 0.0;
	double R11 = 0.0, R12 = 0.0, R13 = 0.0, R14 = 0.0;
	double R21 = 0.0, R22 = 0.0, R23 = 0.0, R24 = 0.0;
	double R31 = 0.0, R32 = 0.0, R33 = 0.0, R34 = 0.0;
	double R41 = 0.0, R42 = 0.0, R43 = 0.0, R44 = 0.0;
	double sqrt_Rho_L = 0.0, sqrt_Rho_R = 0.0;

	// Variables for Roe Averages
	double Roe_u = 0.0, Roe_v = 0.0, Roe_Rho = 0.0, Roe_H = 0.0, Roe_Vmag = 0.0;
	double Un = 0.0, Ut = 0.0, Roe_a = 0.0;

	// Primitive variables
	double Rho_L = 0.0, Rho_R = 0.0, u_L = 0.0, u_R = 0.0, v_L = 0.0, v_R = 0.0;
	double P_L = 0.0, P_R = 0.0, H_L = 0.0, H_R = 0.0, Vmag_L = 0.0, Vmag_R = 0.0;
	double max_eigen_value = 0.0, nx = 0.0, ny = 0.0, Vdotn_L = 0.0, Vdotn_R = 0.0, dl = 0.0;
	double C_L = 0.0, C_R = 0.0, mev_L = 0.0, mev_R = 0.0;

	// Wave strength variables
	double alpha_1 = 0.0, alpha_2 = 0.0, alpha_3 = 0.0, alpha_4 = 0.0;
	double du = 0.0, dv = 0.0, dP = 0.0, dRho = 0.0, dUn = 0.0, dUt = 0.0;

	// Second-order reconstruction variables
	vector<double> d_U(4, 0.0); // Limited difference vector for 2nd order accuracy

	// Initialize dissipative flux to zero
	Dissipative_Flux[0] = 0.0;
	Dissipative_Flux[1] = 0.0;
	Dissipative_Flux[2] = 0.0;
	Dissipative_Flux[3] = 0.0;

	// Face geometry
	nx = Cells[Cell_No].Face_Normals[Face_No * 2 + 0];
	ny = Cells[Cell_No].Face_Normals[Face_No * 2 + 1];
	dl = Cells[Cell_No].Face_Areas[Face_No];

	// Extract primitive variables for left and right states
	// Left State Variables (current cell)
	Rho_L = Primitive_Cells[Cell_No][0];
	u_L = Primitive_Cells[Cell_No][1];
	v_L = Primitive_Cells[Cell_No][2];
	P_L = Primitive_Cells[Cell_No][4];
	C_L = Primitive_Cells[Cell_No][5];
	Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	H_L = ((C_L * C_L) / (gamma - 1.0)) + Vmag_L;

	// Right state Variables (neighbor cell)
	Rho_R = Primitive_Cells[N_Cell_No][0];
	u_R = Primitive_Cells[N_Cell_No][1];
	v_R = Primitive_Cells[N_Cell_No][2];
	P_R = Primitive_Cells[N_Cell_No][4];
	C_R = Primitive_Cells[N_Cell_No][5];
	Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);
	H_R = ((C_R * C_R) / (gamma - 1.0)) + Vmag_R;

	// Calculate Roe averages
	sqrt_Rho_R = sqrt(Rho_R);
	sqrt_Rho_L = sqrt(Rho_L);
	Term1 = 1.0 / (sqrt_Rho_L + sqrt_Rho_R);

	// Roe-averaged quantities
	Roe_Rho = sqrt(Rho_L * Rho_R);
	Roe_u = Term1 * (u_R * sqrt_Rho_R + u_L * sqrt_Rho_L);
	Roe_v = Term1 * (v_R * sqrt_Rho_R + v_L * sqrt_Rho_L);
	Roe_H = Term1 * (H_R * sqrt_Rho_R + H_L * sqrt_Rho_L);
	Roe_Vmag = 0.5 * (Roe_u * Roe_u + Roe_v * Roe_v);

	// Normal and tangential velocities
	Un = nx * Roe_u + ny * Roe_v;					  // Normal velocity
	Ut = -ny * Roe_u + nx * Roe_v;					  // Tangential velocity
	Roe_a = sqrt((gamma - 1.0) * (Roe_H - Roe_Vmag)); // Roe-averaged sound speed

	// Apply second-order reconstruction with slope limiting
	// This computes the limited differences d_U for higher-order accuracy
	Second_Order_Limiter(Cell_No, Face_No, d_U);

	// Calculate differences for wave strength computation
	// These are based on the limited reconstructed states
	du = d_U[1] / Roe_Rho - d_U[0] * Roe_u / (Roe_Rho * Roe_Rho);
	dv = d_U[2] / Roe_Rho - d_U[0] * Roe_v / (Roe_Rho * Roe_Rho);

	// Calculate pressure difference using limited energy difference
	double dE = d_U[3] - 0.5 * Roe_Rho * (2.0 * (Roe_u * du + Roe_v * dv) + du * du + dv * dv) - 0.5 * d_U[0] * Roe_Vmag;
	dP = (gamma - 1.0) * dE;
	dRho = d_U[0];

	// Normal and tangential velocity differences
	dUn = nx * du + ny * dv;
	dUt = -ny * du + nx * dv;

	// Eigenvalues (wave speeds)
	Lambda1 = fabs(Un - Roe_a);
	Lambda2 = fabs(Un);
	Lambda3 = fabs(Un);
	Lambda4 = fabs(Un + Roe_a);

	// Apply entropy fix for sonic points
	double entropy_fix = 0.1 * Roe_a;
	if (Lambda1 < entropy_fix)
		Lambda1 = 0.5 * (Lambda1 * Lambda1 / entropy_fix + entropy_fix);
	if (Lambda4 < entropy_fix)
		Lambda4 = 0.5 * (Lambda4 * Lambda4 / entropy_fix + entropy_fix);

	// Right eigenvectors of the Roe matrix
	// Eigenvector 1: u - a
	R11 = 1.0;
	R12 = Roe_u - nx * Roe_a;
	R13 = Roe_v - ny * Roe_a;
	R14 = Roe_H - Roe_a * Un;

	// Eigenvector 2: u
	R21 = 1.0;
	R22 = Roe_u;
	R23 = Roe_v;
	R24 = Roe_Vmag;

	// Eigenvector 3: u (tangential component)
	R31 = 0.0;
	R32 = -ny;
	R33 = nx;
	R34 = Ut;

	// Eigenvector 4: u + a
	R41 = 1.0;
	R42 = Roe_u + nx * Roe_a;
	R43 = Roe_v + ny * Roe_a;
	R44 = Roe_H + Roe_a * Un;

	// Wave strengths (coefficients of characteristic variables)
	Term2 = 1.0 / (Roe_a * Roe_a);

	alpha_1 = 0.5 * Term2 * (dP - Roe_Rho * Roe_a * dUn);
	alpha_2 = dRho - Term2 * dP;
	alpha_3 = Roe_Rho * dUt;
	alpha_4 = 0.5 * Term2 * (dP + Roe_Rho * Roe_a * dUn);

	// Compute dissipative flux using Roe's approximate Riemann solver
	// D = 0.5 * Σ |λₖ| αₖ Rₖ
	Dissipative_Flux[0] = 0.5 * (Lambda1 * alpha_1 * R11 + Lambda2 * alpha_2 * R21 + Lambda3 * alpha_3 * R31 + Lambda4 * alpha_4 * R41) * dl;

	Dissipative_Flux[1] = 0.5 * (Lambda1 * alpha_1 * R12 + Lambda2 * alpha_2 * R22 + Lambda3 * alpha_3 * R32 + Lambda4 * alpha_4 * R42) * dl;

	Dissipative_Flux[2] = 0.5 * (Lambda1 * alpha_1 * R13 + Lambda2 * alpha_2 * R23 + Lambda3 * alpha_3 * R33 + Lambda4 * alpha_4 * R43) * dl;

	Dissipative_Flux[3] = 0.5 * (Lambda1 * alpha_1 * R14 + Lambda2 * alpha_2 * R24 + Lambda3 * alpha_3 * R34 + Lambda4 * alpha_4 * R44) * dl;

	/*
	 * Second-Order Roe Scheme Summary:
	 *
	 * Key Improvements over First-Order:
	 * 1. Higher spatial accuracy through slope limiting
	 * 2. Reduced numerical diffusion while maintaining TVD property
	 * 3. Better resolution of contact discontinuities and shear layers
	 * 4. Enhanced capturing of fine-scale flow features
	 *
	 * Mathematical Properties:
	 * - Total Variation Diminishing (TVD) through limiters
	 * - Entropy-satisfying with entropy fix
	 * - Conservative and consistent
	 * - Second-order accurate in smooth regions
	 * - First-order accurate near discontinuities (maintains stability)
	 *
	 * Applications:
	 * - High-resolution shock capturing
	 * - Boundary layer computations
	 * - Turbulent flow simulations
	 * - Complex wave interaction problems
	 */
}
