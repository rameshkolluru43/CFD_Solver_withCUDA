#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"

void ROE(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
	// 	cout<<Cell_No<<"\t"<<N_Cell_No<<"\t"<<Face_No<<endl;
	// Variables required for Calculating Roe Averages
	//	Variables for Eigenvalues
	double Lambda1 = 0.0, Lambda2 = 0.0, Lambda3 = 0.0, Lambda4 = 0.0, Term1 = 0.0, Term2 = 0.0;
	double R11 = 0.0, R12 = 0.0, R13 = 0.0, R14 = 0.0, R21 = 0.0, R22 = 0.0, R23 = 0.0, R24 = 0.0, R31 = 0.0, R32 = 0.0, R33 = 0.0, R34 = 0.0, R41 = 0.0, R42 = 0.0, R43 = 0.0, R44 = 0.0, sqrt_Rho_L = 0.0, sqrt_Rho_R = 0.0;

	//	Variables for Roe Averages
	Roe_u = 0.0, Roe_v = 0.0, Roe_Rho = 0.0, Roe_H = 0.0, Roe_Vmag = 0.0, Un = 0.0, Ut = 0.0, Roe_a = 0.0;

	Rho_L = 0.0, Rho_R = 0.0, u_L = 0.0, u_R = 0.0, v_L = 0.0, v_R = 0.0, P_L = 0.0, P_R = 0.0, H_L = 0.0, H_R = 0.0, Vmag_L = 0.0, Vmag_R = 0.0;
	max_eigen_value = 0.0, nx = 0.0, ny = 0.0, Vdotn_L = 0.0, Vdotn_R = 0.0, dl = 0.0, C_L = 0.0, C_R = 0.0, mev_L = 0.0, mev_R = 0.0;

	double alpha_1 = 0.0, alpha_2 = 0.0, alpha_3 = 0.0, alpha_4 = 0.0, du = 0.0, dv = 0.0, dP = 0.0, dRho = 0.0, dUn = 0.0, dUt = 0.0;

	//      Initializing the Cell dissipation to be zero to avoid any type of over writing
	Dissipative_Flux[0] = 0.0;
	Dissipative_Flux[1] = 0.0;
	Dissipative_Flux[2] = 0.0;
	Dissipative_Flux[3] = 0.0;

	// Face Normals
	/*nx = Cell_Face_Normals[Cell_No][Face_No*2+0];
	ny = Cell_Face_Normals[Cell_No][Face_No*2+1];
	dl = Cell_Face_Areas[Cell_No][Face_No];*/
	nx = Cells[Cell_No].Face_Normals[Face_No * 2 + 0];
	ny = Cells[Cell_No].Face_Normals[Face_No * 2 + 1];
	dl = Cells[Cell_No].Face_Areas[Face_No];

	//	Right state Variables
	Rho_R = Primitive_Cells[N_Cell_No][0];
	u_R = Primitive_Cells[N_Cell_No][1];
	v_R = Primitive_Cells[N_Cell_No][2];
	P_R = Primitive_Cells[N_Cell_No][4];
	C_R = Primitive_Cells[N_Cell_No][5];
	Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);
	H_R = ((C_R * C_R) / (gamma - 1.0)) + Vmag_R;
	//	Left State Variables
	Rho_L = Primitive_Cells[Cell_No][0];
	u_L = Primitive_Cells[Cell_No][1];
	v_L = Primitive_Cells[Cell_No][2];
	P_L = Primitive_Cells[Cell_No][4];
	Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	H_L = ((C_L * C_L) / (gamma - 1.0)) + Vmag_L;

	sqrt_Rho_R = sqrt(Rho_R);
	sqrt_Rho_L = sqrt(Rho_L);
	Term1 = 1.0 / (sqrt_Rho_L + sqrt_Rho_R);

	// 		Evaluating the Roe averages
	Roe_Rho = sqrt(Rho_L * Rho_R);
	Roe_u = Term1 * (u_R * sqrt_Rho_R + u_L * sqrt_Rho_L);
	Roe_v = Term1 * (v_R * sqrt_Rho_R + v_L * sqrt_Rho_L);
	Roe_H = Term1 * (H_R * sqrt_Rho_R + H_L * sqrt_Rho_L);
	Roe_Vmag = 0.5 * (Roe_u * Roe_u + Roe_v * Roe_v);

	Un = nx * Roe_u + ny * Roe_v;
	Ut = -ny * Roe_u + nx * Roe_v;
	Roe_a = sqrt((gamma - 1.0) * (Roe_H - Roe_Vmag));

	du = (u_R - u_L);
	dv = (v_R - v_L);
	dP = (P_R - P_L);
	dRho = (Rho_R - Rho_L);

	dUn = nx * du + ny * dv;
	dUt = -ny * du + nx * dv;

	// 	Eigen Vector Evaluation
	R11 = 1.0;
	R12 = Roe_u - nx * Roe_a;
	R13 = Roe_v - ny * Roe_a;
	R14 = Roe_H - Roe_a * Un; // Eigen Vector corresponding to u -a

	R21 = 1.0;
	R22 = Roe_u;
	R23 = Roe_v;
	R24 = Roe_Vmag; // Eigen Vector corresponding to u

	R31 = 0.0;
	R32 = -ny;
	R33 = nx;
	R34 = Ut; // Eigen Vector corresponding to u

	R41 = 1.0;
	R42 = Roe_u + nx * Roe_a;
	R43 = Roe_v + ny * Roe_a;
	R44 = Roe_H + Roe_a * Un; // Eigen Vector corresponding to u + a

	//	Evaluating Wave Strengths
	Term2 = 1.0 / (Roe_a * Roe_a);

	alpha_1 = 0.5 * Term2 * (dP - Roe_Rho * Roe_a * dUn);

	alpha_2 = dRho - Term2 * dP;

	alpha_3 = Roe_Rho * dUt;

	alpha_4 = 0.5 * Term2 * (dP + Roe_Rho * Roe_a * dUn);
	//	Evaluating Eigen Values
	Lambda1 = fabs(Un - Roe_a);
	Lambda2 = fabs(Un);
	Lambda3 = fabs(Un);
	Lambda4 = fabs(Un + Roe_a);
	//      Evaluating the Eigen Values
	Dissipative_Flux[0] = 0.5 * (Lambda1 * alpha_1 * R11 + Lambda2 * alpha_2 * R21 + Lambda3 * alpha_3 * R31 + Lambda4 * alpha_4 * R41) * dl;
	Dissipative_Flux[1] = 0.5 * (Lambda1 * alpha_1 * R12 + Lambda2 * alpha_2 * R22 + Lambda3 * alpha_3 * R32 + Lambda4 * alpha_4 * R42) * dl;
	Dissipative_Flux[2] = 0.5 * (Lambda1 * alpha_1 * R13 + Lambda2 * alpha_2 * R23 + Lambda3 * alpha_3 * R33 + Lambda4 * alpha_4 * R43) * dl;
	Dissipative_Flux[3] = 0.5 * (Lambda1 * alpha_1 * R14 + Lambda2 * alpha_2 * R24 + Lambda3 * alpha_3 * R34 + Lambda4 * alpha_4 * R44) * dl;
}
void ROE_2O(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
}
