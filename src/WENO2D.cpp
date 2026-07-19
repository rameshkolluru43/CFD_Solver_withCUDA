#include "definitions.h"
#include "Globals.h"
#include "Weno.h"
#include "Flux.h"
#include "Primitive_Computational.h"
#include "Utilities.h"
#include "Timestep.h"
#include <cmath>	// For std::isfinite and other math functions
#include <iostream> // For error messages
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
	double Pressure_From_U(const V_D &U)
	{
		const double rho = U[0];
		if (rho <= 1e-14 || !std::isfinite(rho))
			return -1.0;
		const double u = U[1] / rho;
		const double v = U[2] / rho;
		const double p = gamma_M_1 * (U[3] - 0.5 * rho * (u * u + v * v));
		return p;
	}

	/* Thread-safe primitive recovery for WENO face states (no global M / mu / K). */
	void WenoPrimitiveFromU(const V_D &U_Vect, V_D &G_Primitive)
	{
		Vector_Reset(G_Primitive);
		const double Rho = U_Vect[0];
		CFD_RequireFinitePositive(Rho, "WenoPrimitiveFromU density");
		const double inv_Density = 1.0 / Rho;
		double v1 = U_Vect[1] * inv_Density;
		double v2 = U_Vect[2] * inv_Density;
		if (fabs(v1) <= 1e-10)
			v1 = 0.0;
		if (fabs(v2) <= 1e-10)
			v2 = 0.0;
		const double vmag = (v1 * v1 + v2 * v2);
		const double Pressure = gamma_M_1 * (U_Vect[3] - 0.5 * Rho * vmag);
		CFD_RequireFinitePositive(Pressure, "WenoPrimitiveFromU pressure");
		double Temperature;
		if (Non_Dimensional_Form)
			Temperature = (Pressure / (Rho * R_ref));
		else
			Temperature = (Pressure * inv_Density / R_GC);
		const double C = sqrt(gamma * Pressure * inv_Density);
		CFD_RequireFinitePositive(C, "WenoPrimitiveFromU sound speed");
		const double Mloc = sqrt(vmag) / C;

		G_Primitive[0] = Rho;
		G_Primitive[1] = v1;
		G_Primitive[2] = v2;
		G_Primitive[3] = Temperature;
		G_Primitive[4] = Pressure;
		G_Primitive[5] = C;
		G_Primitive[6] = U_Vect[3];
		G_Primitive[7] = Mloc;
		G_Primitive[8] = 0.0;
		G_Primitive[9] = 0.0;
		G_Primitive[10] = Pressure * pow((1.0 + 0.5 * gamma_M_1 * Mloc * Mloc), gamma / gamma_M_1);

		if (std::isnan(C))
		{
			throw runtime_error("WenoPrimitiveFromU: non-finite sound speed");
		}
	}
}
void WENO_Reconstruction(double &a, double &b, double &c, double &d, double &e, int &shift, double &U)
{
	// Check for NaN or infinite values in input
	if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d) || !std::isfinite(e))
	{
		std::cout << "Warning: Non-finite values in WENO reconstruction input" << std::endl;
		U = c; // Fall back to central value
		return;
	}

	double d0 = 0.0, d1 = 0.0, d2 = 0.0, b0 = 0.0, b1 = 0.0, b2 = 0.0, a0 = 0.0, a1 = 0.0, a2 = 0.0;
	double v0 = 0.0, v1 = 0.0, v2 = 0.0, w0 = 0.0, w1 = 0.0, w2 = 0.0, sum = 0.0, epsilon = 1e-6;

	switch (shift)
	{
	case 0: // for Left Values
		d0 = 3.0 / 10.0;
		d1 = 6.0 / 10.0;
		d2 = 1.0 / 10.0;

		v2 = (2.0 * a - 7.0 * b + 11.0 * c) / 6.0;
		v1 = (-1.0 * b + 5.0 * c + 2.0 * d) / 6.0;
		v0 = (2.0 * c + 5.0 * d - 1.0 * e) / 6.0;

		break;
	case 1: // for Right Values
		d0 = 1.0 / 10.0;
		d1 = 6.0 / 10.0;
		d2 = 3.0 / 10.0;

		v2 = (-1.0 * a + 5.0 * b + 2.0 * c) / 6.0;
		v1 = (2.0 * b + 5.0 * c - 1.0 * d) / 6.0;
		v0 = (11.0 * c - 7.0 * d + 2.0 * e) / 6.0;

		break;
	}

	{
		const double t2 = (a - 2.0 * b + c);
		const double t2b = (a - 4.0 * b + 3.0 * c);
		b2 = (13.0 / 12.0) * (t2 * t2) + (1.0 / 4.0) * (t2b * t2b);
	}
	{
		const double t1 = (b - 2.0 * c + d);
		const double t1b = (b - d);
		b1 = (13.0 / 12.0) * (t1 * t1) + (1.0 / 4.0) * (t1b * t1b);
	}
	{
		const double t0 = (c - 2.0 * d + e);
		const double t0b = (3.0 * c - 4.0 * d + e);
		b0 = (13.0 / 12.0) * (t0 * t0) + (1.0 / 4.0) * (t0b * t0b);
	}

	{
		const double s0 = epsilon + b0;
		const double s1 = epsilon + b1;
		const double s2 = epsilon + b2;
		a0 = d0 / (s0 * s0);
		a1 = d1 / (s1 * s1);
		a2 = d2 / (s2 * s2);
	}

	sum = a0 + a1 + a2;

	// Check for degenerate case where sum is too small
	if (sum < epsilon)
	{
		// Fall back to simple average of the three candidate values
		U = (v0 + v1 + v2) / 3.0;
	}
	else
	{
		w0 = a0 / sum;
		w1 = a1 / sum;
		w2 = a2 / sum;
		U = w0 * v0 + w1 * v1 + w2 * v2;
	}

	// Final validation
	if (!std::isfinite(U))
	{
		std::cout << "Warning: Non-finite result in WENO reconstruction, using central value" << std::endl;
		U = c;
	}
}

void Get_LR(int &CellNo, int &Neighbour, const int &Face_No, V_D &L, V_D &IL)
{

	int index = Face_No * 2;

	double dL = 0.0, uL = 0.0, vL = 0.0, aL = 0.0, dR = 0.0, uR = 0.0, vR = 0.0, aR = 0.0, h = 0.0;
	double ek = 0.0, vn = 0.0, nx = 0.0, ny = 0.0;
	double a_RL = 0.0, u_RL = 0.0, v_RL = 0.0, t2, t1;

	//	cout<<"in Get LR"<<endl;
	//	cout<<CellNo<<"\t"<<Neighbour<<"\t"<<Face_No<<"\t"<<index<<endl;

	dL = Primitive_Cells[CellNo][0];
	uL = Primitive_Cells[CellNo][1];
	vL = Primitive_Cells[CellNo][2];
	aL = Primitive_Cells[CellNo][5];

	dR = Primitive_Cells[Neighbour][0];
	uR = Primitive_Cells[Neighbour][1];
	vR = Primitive_Cells[Neighbour][2];
	aR = Primitive_Cells[Neighbour][5];

	double sqrt_dL = sqrt(fmax(dL, 1e-14));
	double sqrt_dR = sqrt(fmax(dR, 1e-14));
	double denom = sqrt_dL + sqrt_dR;

	if (denom < 1e-14)
	{
		// Handle degenerate case - use simple average
		u_RL = 0.5 * (uL + uR);
		v_RL = 0.5 * (vL + vR);
		a_RL = 0.5 * (aL + aR);
	}
	else
	{
		u_RL = (uL * sqrt_dL + uR * sqrt_dR) / denom;
		v_RL = (vL * sqrt_dL + vR * sqrt_dR) / denom;
		a_RL = (aL * sqrt_dL + aR * sqrt_dR) / denom;
	}

	nx = Cells[CellNo].Face_Normals[index + 0]; //----------------- nx = dy/dl
	ny = Cells[CellNo].Face_Normals[index + 1]; //----------------- ny = -dx/dl

	//	cout<<" Cell Normals\t"<<nx<<"\t"<<ny<<endl;

	vn = (u_RL * nx + v_RL * ny); // V dot n = (u dy - v dx) / dl
	ek = 0.5 * (u_RL * u_RL + v_RL * v_RL);
	h = (a_RL * a_RL / gamma_M_1) + ek;

	// Add safety check for speed of sound
	if (a_RL < 1e-14)
	{
		std::cout << "Warning: Very small speed of sound in WENO, a_RL = " << a_RL << std::endl;
		a_RL = 1e-14;
	}

	t1 = 0.5 / (a_RL * a_RL); // Fixed: added parentheses for correct division
	t2 = gamma_M_1 * t1;

	// Matrix to Conserved to Characteristic

	L[0] = 1.0 - 2.0 * t2 * ek;
	L[4] = 2.0 * t2 * u_RL;
	L[8] = 2.0 * t2 * v_RL;
	L[12] = -2.0 * t2;
	L[1] = v_RL * nx - u_RL * ny;
	L[5] = ny;
	L[9] = -nx;
	L[13] = 0.0;
	L[2] = t2 * ek - a_RL * vn * t1;
	L[6] = -t2 * u_RL + a_RL * t1 * nx;
	L[10] = -t2 * v_RL + a_RL * t1 * ny;
	L[14] = t2;
	L[3] = t2 * ek + a_RL * vn * t1;
	L[7] = -t2 * u_RL - a_RL * t1 * nx;
	L[11] = -t2 * v_RL - a_RL * t1 * ny;
	L[15] = t2;

	//	Print(L);
	for (int i = 0; i < 16; i++)
	{
		if (isnan(L[i]))
		{
			throw runtime_error("Get_LR: non-finite WENO characteristic matrix L");
		}
	}

	// Matrix to convert characteristic to conserved
	IL[0] = 1.0;
	IL[4] = 0.0;
	IL[8] = 1.0;
	IL[12] = 1.0;
	IL[1] = u_RL;
	IL[5] = ny;
	IL[9] = u_RL + a_RL * nx;
	IL[13] = u_RL - a_RL * nx;
	IL[2] = v_RL;
	IL[6] = -nx;
	IL[10] = v_RL + a_RL * ny;
	IL[14] = v_RL - a_RL * ny;
	//	IL[1] = vn;		IL[5] = ny;					IL[9]  = vn + a_RL*nx;		IL[13] = vn - a_RL*nx;
	//	IL[2] = vn;		IL[6] = -nx;				IL[10] = vn + a_RL*ny;		IL[14] = vn - a_RL*ny;
	IL[3] = ek;
	IL[7] = u_RL * ny - v_RL * nx;
	IL[11] = h + a_RL * vn;
	IL[15] = h - a_RL * vn;

	//		Print(IL);
	for (int i = 0; i < 16; i++)
	{
		if (isnan(IL[i]))
		{
			throw runtime_error("Get_LR: non-finite WENO inverse characteristic matrix IL");
		}
	}
}

void MatVecMul(V_D &A, V_D &x, V_D &B)
{
	int size = B.size();
	for (int i = 0; i < size; i++)
	{
		B[i] = 0.0;
		for (int j = 0; j < size; j++)
		{
			//			cout<<i<<"\t"<<i + size * j<<"\t"<<j<<endl;
			B[i] += A[i + size * j] * x[j];
		}
	}
}

void Get_Reconstructed_U(int &Cell_No, const int &Face_No, int &i1, int &i2, int &i3, int &i4, int &i5, int &LR, V_D &U)
{

	int Neighbour = 0;
	switch (Face_No)
	{
	case 0:
		Neighbour = Cells[Cell_No].Neighbours[0];
		break;
	case 1:
		Neighbour = Cells[Cell_No].Neighbours[1];
		break;
	case 2:
		Neighbour = Cells[Cell_No].Neighbours[2];
		break;
	case 3:
		Neighbour = Cells[Cell_No].Neighbours[3];
		break;
	}

	thread_local V_D U1(4, 0.0), U2(4, 0.0), U3(4, 0.0), U4(4, 0.0), U5(4, 0.0);
	thread_local V_D W1(4, 0.0), W2(4, 0.0), W3(4, 0.0), W4(4, 0.0), W5(4, 0.0), W(4, 0.0);
	thread_local V_D L(16, 0.0), InvL(16, 0.0);

	//	cout<<"Starting reconstruction"<<endl;
	//	cout<<i1<<"\t"<<i2<<"\t"<<i3<<"\t"<<i4<<"\t"<<i5<<endl;
	Calculate_Computational_Variables(i1, U1);
	Calculate_Computational_Variables(i2, U2);
	Calculate_Computational_Variables(i3, U3);
	Calculate_Computational_Variables(i4, U4);
	Calculate_Computational_Variables(i5, U5);

	if (Is_Char)
	{
		Get_LR(Cell_No, Neighbour, Face_No, L, InvL);
		MatVecMul(L, U1, W1);
		MatVecMul(L, U2, W2);
		MatVecMul(L, U3, W3);
		MatVecMul(L, U4, W4);
		MatVecMul(L, U5, W5);
		for (int i = 0; i < 4; i++)
		{
			WENO_Reconstruction(W1[i], W2[i], W3[i], W4[i], W5[i], LR, W[i]);
		}
		MatVecMul(InvL, W, U);
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			WENO_Reconstruction(U1[i], U2[i], U3[i], U4[i], U5[i], LR, U[i]);
		}
	}
}

void WENO_Reconstruction_X(int &Cell_No, const int &Face_No, V_D &U_L, V_D &U_R)
{
	int im3 = 0, im2 = 0, im1 = 0, ip2 = 0, ip1 = 0, ip3 = 0, jm3 = 0, jm2 = 0, jm1 = 0, jp2 = 0, jp1 = 0, jp3 = 0, LR = 0;

	// Obtaining the indicies of neighbour cells

	//	cout<<"In weno construction x\t"<< Cell_No<<"\t"<<Face_No<<endl;

	/* Neighbours[0..3] = left, bottom, right, top (same as Grid_Computations / flux). */
	im1 = Cells[Cell_No].Neighbours[0];
	ip1 = Cells[Cell_No].Neighbours[2];
	jm1 = Cells[Cell_No].Neighbours[1];
	jp1 = Cells[Cell_No].Neighbours[3];

	if (im1 >= No_Physical_Cells)
	{
		// left most boundary im1 is ghost cell then im2 and im3 = im1
		im2 = im1;
		im3 = im1;
	}
	else
	{
		im2 = Cells[im1].Neighbours[0];
		if (im2 >= No_Physical_Cells || im2 < 0)
		{
			im3 = im2;
		}
		else
		{
			im3 = Cells[im2].Neighbours[0];
			if (im3 < 0)
				im3 = im2; // Additional safety check
		}
	}

	//	cout<<"left side stencil fixed"<<endl;
	if (ip1 >= No_Physical_Cells)
	{
		ip2 = ip1;
		ip3 = ip1;
	}
	else
	{
		ip2 = Cells[ip1].Neighbours[2];
		//		cout<<"value of ip2\t"<<ip2<<endl;
		if (ip2 >= No_Physical_Cells || ip2 < 0)
		{
			ip3 = ip2;
		}
		else
		{
			ip3 = Cells[ip2].Neighbours[2];
			if (ip3 < 0)
				ip3 = ip2; // Additional safety check
		}
		//		cout<<"value of ip3\t"<<ip3<<endl;
	}

	// For 	y direction neighbours
	if (jm1 >= No_Physical_Cells)
	{
		// left most boundary im1 is ghost cell then im2  and im3 = im1
		jm2 = jm1;
		jm3 = jm1;
	}
	else
	{
		jm2 = Cells[jm1].Neighbours[1];
		if (jm2 >= No_Physical_Cells)
		{
			jm3 = jm2;
		}
		else
		{
			jm3 = Cells[jm2].Neighbours[1];
		}
	}

	if (jp1 >= No_Physical_Cells)
	{
		jp2 = jp1;
		jp3 = jp1;
	}
	else
	{
		jp2 = Cells[jp1].Neighbours[3];
		if (jp2 >= No_Physical_Cells)
		{
			jp3 = jp2;
		}
		else
		{
			jp3 = Cells[jp2].Neighbours[3];
		}
	}
	//	cout<<"In weno construction x\t"<< Cell_No<<"\t"<<Face_No<<endl;
	//	cout<<im3<<"\t"<<im2<<"\t"<<im1<<"\t"<<ip1<<"\t"<<ip2<<"\t"<<ip3<<endl;
	switch (Face_No)
	{
	case 0:
		LR = 0;
		Get_Reconstructed_U(Cell_No, Face_No, im3, im2, im1, Cell_No, ip1, LR, U_L);
		LR = 1;
		Get_Reconstructed_U(Cell_No, Face_No, im2, im1, Cell_No, ip1, ip2, LR, U_R);
		break;
	case 1:
		LR = 0;
		Get_Reconstructed_U(Cell_No, Face_No, jm3, jm2, jm1, Cell_No, jp1, LR, U_L);
		LR = 1;
		Get_Reconstructed_U(Cell_No, Face_No, jm2, jm1, Cell_No, jp1, jp2, LR, U_R);
		break;
	case 2:
		LR = 0;
		Get_Reconstructed_U(Cell_No, Face_No, im2, im1, Cell_No, ip1, ip2, LR, U_L);
		LR = 1;
		Get_Reconstructed_U(Cell_No, Face_No, im1, Cell_No, ip1, ip2, ip3, LR, U_R);
		break;
	case 3:
		LR = 0;
		Get_Reconstructed_U(Cell_No, Face_No, jm2, jm1, Cell_No, jp1, jp2, LR, U_L);
		LR = 1;
		Get_Reconstructed_U(Cell_No, Face_No, jm1, Cell_No, jp1, jp2, jp3, LR, U_R);
		break;
	}
}

void Evaluate_Cell_Net_Flux_WENO()
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (No_Physical_Cells > 512)
#endif
	for (int Current_Cell_No = 0; Current_Cell_No < No_Physical_Cells; Current_Cell_No++)
	{
		if (Cells[Current_Cell_No].numFaces != 4)
		{
			for (int i = 0; i < NUM_FLUX_COMPONENTS; i++)
				Cells_Net_Flux[Current_Cell_No][i] = 0.0;
			Evaluate_Time_Step(Current_Cell_No);
			continue;
		}

		const int N_1 = Cells[Current_Cell_No].Neighbours[0];
		const int N_2 = Cells[Current_Cell_No].Neighbours[1];
		const int N_3 = Cells[Current_Cell_No].Neighbours[2];
		const int N_4 = Cells[Current_Cell_No].Neighbours[3];

		for (int i = 0; i < NUM_FLUX_COMPONENTS; i++)
			Cells_Net_Flux[Current_Cell_No][i] = 0.0;

		V_D fc_avg(4, 0.0), fc_diss(4, 0.0);

		Calculate_Face_WENO_Flux(Current_Cell_No, N_1, Face_0, Cells_Face_Boundary_Type[Current_Cell_No][Face_0], fc_avg, fc_diss);
		for (int i = 0; i < 4; i++)
			Cells_Net_Flux[Current_Cell_No][i] += fc_avg[i] - fc_diss[i];

		Calculate_Face_WENO_Flux(Current_Cell_No, N_2, Face_1, Cells_Face_Boundary_Type[Current_Cell_No][Face_1], fc_avg, fc_diss);
		for (int i = 0; i < 4; i++)
			Cells_Net_Flux[Current_Cell_No][i] += fc_avg[i] - fc_diss[i];

		Calculate_Face_WENO_Flux(Current_Cell_No, N_3, Face_2, Cells_Face_Boundary_Type[Current_Cell_No][Face_2], fc_avg, fc_diss);
		for (int i = 0; i < 4; i++)
			Cells_Net_Flux[Current_Cell_No][i] += fc_avg[i] - fc_diss[i];

		Calculate_Face_WENO_Flux(Current_Cell_No, N_4, Face_3, Cells_Face_Boundary_Type[Current_Cell_No][Face_3], fc_avg, fc_diss);
		for (int i = 0; i < 4; i++)
			Cells_Net_Flux[Current_Cell_No][i] += fc_avg[i] - fc_diss[i];

		Evaluate_Time_Step(Current_Cell_No);
	}
}

void Calculate_Face_WENO_Flux(int Cell_No, int N_Cell_No, const int &Face_No, bool Is_Wall_Face, V_D &out_avg, V_D &out_diss)
{
	V_D U_L(4, 0.0), U_R(4, 0.0), Flux_L(4, 0.0), Flux_R(4, 0.0);
	V_D d_U(4, 0.0), d_F(4, 0.0), Mod_Alpha(4, 0.0);
	V_D ricca_diss(4, 0.0), llf_diss(4, 0.0);
	V_D primL(11, 0.0), primR(11, 0.0);

	double Rho_L = 0.0, P_L = 0.0, u_L = 0.0, v_L = 0.0, Vdotn_L = 0.0, C_L = 0.0, Vmag_L = 0.0;
	double Rho_R = 0.0, P_R = 0.0, u_R = 0.0, v_R = 0.0, Vdotn_R = 0.0, C_R = 0.0, Vmag_R = 0.0;
	double nx = 0.0, ny = 0.0, dl = 0.0;

	const int index = Face_No * 2;

	/* Wall faces: WENO5 needs a 5-point stencil; host ghosts collapse to one cell
	   and frequently produce negative pressure at M=6 (e.g. cell 5504). Use
	   robust 1O interior/ghost conservatives from the BC-filled neighbour. */
	if (Is_Wall_Face || N_Cell_No < 0 || N_Cell_No >= static_cast<int>(U_Cells.size()))
	{
		U_L = U_Cells[Cell_No];
		if (N_Cell_No >= 0 && N_Cell_No < static_cast<int>(U_Cells.size()))
			U_R = U_Cells[N_Cell_No];
		else
			U_R = U_L;
	}
	else
	{
		WENO_Reconstruction_X(Cell_No, Face_No, U_L, U_R);

		if (Face_No == 0 || Face_No == 1)
			std::swap(U_L, U_R);

		const double pL = Pressure_From_U(U_L);
		const double pR = Pressure_From_U(U_R);
		if (pL < 1e-12 || pR < 1e-12 || !std::isfinite(pL) || !std::isfinite(pR))
		{
			U_L = U_Cells[Cell_No];
			U_R = U_Cells[N_Cell_No];
		}
	}

	WenoPrimitiveFromU(U_L, primL);
	Rho_L = primL[0];
	P_L = primL[4];
	u_L = primL[1];
	v_L = primL[2];
	C_L = primL[5];
	WenoPrimitiveFromU(U_R, primR);
	Rho_R = primR[0];
	P_R = primR[4];
	u_R = primR[1];
	v_R = primR[2];
	C_R = primR[5];

	nx = Cells[Cell_No].Face_Normals[index + 0];
	ny = Cells[Cell_No].Face_Normals[index + 1];
	dl = Cells[Cell_No].Face_Areas[Face_No];

	Vdotn_L = (u_L * nx + v_L * ny);
	Vdotn_R = (u_R * nx + v_R * ny);

	Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

	if (Is_Wall_Face)
		P_R = P_L;

	Flux_L[0] = Rho_L * Vdotn_L * dl;
	Flux_L[1] = Rho_L * u_L * Vdotn_L * dl + P_L * nx * dl;
	Flux_L[2] = Rho_L * v_L * Vdotn_L * dl + P_L * ny * dl;
	Flux_L[3] = (gamma1 * P_L + Rho_L * Vmag_L) * Vdotn_L * dl;

	Flux_R[0] = Rho_R * Vdotn_R * dl;
	Flux_R[1] = Rho_R * u_R * Vdotn_R * dl + P_R * nx * dl;
	Flux_R[2] = Rho_R * v_R * Vdotn_R * dl + P_R * ny * dl;
	Flux_R[3] = (gamma1 * P_R + Rho_R * Vmag_R) * Vdotn_R * dl;

	if (C_L < 1e-14)
		C_L = 1e-14;
	if (C_R < 1e-14)
		C_R = 1e-14;

	double S0 = fabs(Vdotn_L - C_L) * dl;
	double S1 = fabs(Vdotn_L + C_L) * dl;
	double S2 = fabs(Vdotn_L) * dl;
	double S3 = fabs(Vdotn_R - C_R) * dl;
	double S4 = fabs(Vdotn_R + C_R) * dl;
	double S5 = fabs(Vdotn_R) * dl;

	double Max1 = 0.0, Max2 = 0.0, max_eigen_value = 0.0;
	Maximum(S0, S1, S2, Max1);
	Maximum(S3, S4, S5, Max2);
	Maximum(Max1, Max2, max_eigen_value);

	for (int i = 0; i < 4; i++)
		out_avg[i] = 0.5 * (Flux_L[i] + Flux_R[i]);

	for (int i = 0; i < 4; i++)
		llf_diss[i] = 0.5 * max_eigen_value * (U_R[i] - U_L[i]);

	if (Dissipation_Type == 4 || Dissipation_Type == 6)
	{
		d_U[0] = U_R[0] - U_L[0];
		d_U[1] = U_R[1] - U_L[1];
		d_U[2] = U_R[2] - U_L[2];
		d_U[3] = U_R[3] - U_L[3];

		d_F[0] = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
		d_F[1] = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L + P_L * nx) * dl;
		d_F[2] = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L + P_L * ny) * dl;
		d_F[3] = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R -
				  (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) *
				 dl;

		double dP = fabs(P_L - P_R);
		double P_I = 0.5 * (P_L + P_R);
		double Rho_I = 0.5 * (Rho_L + Rho_R);

		for (int i = 0; i < 4; ++i)
		{
			Condition_For_RICCA(d_U[i], d_F[i], Vdotn_L, Vdotn_R, dP, Rho_I, P_I, Mod_Alpha[i]);
			ricca_diss[i] = 0.5 * Mod_Alpha[i] * dl * d_U[i];
		}

		if (Dissipation_Type == 6)
		{
			const double denom = fabs(P_L) + fabs(P_R) + 1e-14;
			const double pressure_sensor = std::min(1.0, 10.0 * fabs(P_R - P_L) / denom);
			for (int i = 0; i < 4; ++i)
			{
				if (!std::isfinite(ricca_diss[i]))
					out_diss[i] = llf_diss[i];
				else
					out_diss[i] = (1.0 - pressure_sensor) * ricca_diss[i] + pressure_sensor * llf_diss[i];
			}
		}
		else
		{
			for (int i = 0; i < 4; ++i)
				out_diss[i] = std::isfinite(ricca_diss[i]) ? ricca_diss[i] : llf_diss[i];
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
			out_diss[i] = llf_diss[i];
	}
}