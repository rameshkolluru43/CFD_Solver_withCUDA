#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"

void Condition_For_MOVERS_NWSC(double &d_U, double &d_F, double &Alpha)
{
	double epsilon = 1e-8;
	//***************************************************************************************************************
	//	cout<<"NWSC funtion\t"<<d_U<<"\t"<<d_F<<"\t"<<L_Min<<Alpha<<endl;
	Alpha = 0.0;
	// Condition for discontinuity.... across which Flux is zero
	if (fabs(d_F) <= epsilon and fabs(d_U) >= epsilon)
	{
		Alpha = 0.0;
	}
	else
	{
		Alpha = Sign(d_U) * fabs(d_F);
	}
}

void MOVERS_NWSC(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{

	Rho_L = 0.0, P_L = 0.0, C_L = 0.0, u_L = 0.0, v_L = 0.0, Vdotn_L = 0.0, nx = 0.0, ny = 0.0, Vmag_L = 0.0;
	Rho_R = 0.0, P_R = 0.0, C_R = 0.0, u_R = 0.0, v_R = 0.0, Vdotn_R = 0.0, Vmag_R = 0.0, dl = 0.0;

	Mod_Alpha0 = 0.0, Mod_Alpha1 = 0.0, Mod_Alpha2 = 0.0, Mod_Alpha3 = 0.0, d_F_0 = 0.0, d_F_1 = 0.0, d_F_2 = 0.0;
	d_F_3 = 0.0, d_U_0 = 0.0, d_U_1 = 0.0, d_U_2 = 0.0, d_U_3 = 0.0;

	double Alpha_P = 0.0, P_I = 0.0, d_P = 0.0, Sensor = 0.0, beta = 0.21, a = 0.5;
	// vector<double> d_U(4,0.0);
	int index = Face_No * 2, Which_Var = 0;

	Dissipative_Flux[0] = 0.0;
	Dissipative_Flux[1] = 0.0;
	Dissipative_Flux[2] = 0.0;
	Dissipative_Flux[3] = 0.0;

	//      Left state Variables, Density, Pressure, u, v, speed of cound C

	Rho_L = Primitive_Cells[Cell_No][0];
	P_L = Primitive_Cells[Cell_No][4];
	u_L = Primitive_Cells[Cell_No][1];
	v_L = Primitive_Cells[Cell_No][2];

	//      Right state Variables, Density, Pressure, u, v, speed of cound C
	Rho_R = Primitive_Cells[N_Cell_No][0];
	P_R = Primitive_Cells[N_Cell_No][4];
	u_R = Primitive_Cells[N_Cell_No][1];
	v_R = Primitive_Cells[N_Cell_No][2];

	/*        	if(P_R<0.0  or P_L<0.0)
				{
					cout<<Cell_No<<"\t"<<P_R<<"\t"<<P_L<<"\t"<<Primitive_Cells[N_Cell_No][7]<<"\t"<<Primitive_Cells[Cell_No][7]<<endl;
				}*/

	//      Magnitudes of Velocity on Left side and Right side of an interface
	Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

	//	cout<<u_L<<"\t"<<u_R<<"\t"<<v_L<<"\t"<<v_R<<endl;
	//      nx and ny are normals to the face and dl - face length
	/*   nx = Cell_Face_Normals[Cell_No][index + 0];
	   ny = Cell_Face_Normals[Cell_No][index + 1];
	   dl = Cell_Face_Areas[Cell_No][Face_No];*/
	nx = Cells[Cell_No].Face_Normals[index + 0];
	ny = Cells[Cell_No].Face_Normals[index + 1];
	dl = Cells[Cell_No].Face_Areas[Face_No];

	//	cout<<nx<<"\t"<<ny<<"\t"<<dl<<endl;

	//      Normal Velocity taking Left state Velocity and Right State Velocity of an interface
	Vdotn_L = (u_L * nx + v_L * ny);
	Vdotn_R = (u_R * nx + v_R * ny);

	//          Evaluating Conserved variable Difference between Left and Right faces

	for (int i = 0; i < 4; i++)
	{
		d_U[i] = U_Cells[N_Cell_No][i] - U_Cells[Cell_No][i];
	}

	//          Evaluating Flux Difference between Left and Right Faces
	d_F_0 = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
	d_F_1 = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L + P_L * nx) * dl;
	d_F_2 = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L + P_L * ny) * dl;
	d_F_3 = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R - (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) * dl;
	Condition_For_MOVERS_NWSC(d_U[0], d_F_0, Mod_Alpha0);
	Condition_For_MOVERS_NWSC(d_U[1], d_F_1, Mod_Alpha1);
	Condition_For_MOVERS_NWSC(d_U[2], d_F_2, Mod_Alpha2);
	Condition_For_MOVERS_NWSC(d_U[3], d_F_3, Mod_Alpha3);

	d_P = P_R - P_L;
	P_I = 0.5 * (P_R + P_L);
	Alpha_P = 0.5 * (fabs(Vdotn_R * dl) + fabs(Vdotn_L * dl)); // + (1-3.5*beta)*0.5*(Primitive_Cells[Cell_No][5] +Primitive_Cells[N_Cell_No][5])*dl;
	Sensor = beta * fabs(d_P / (2.0 * P_I));

	//			      cout<<Sensor<<"\t"<< Mod_Alpha0<<"\t"<< Mod_Alpha1<<"\t"<<Mod_Alpha2<<"\t"<<Mod_Alpha3<<"\t"<<Alpha_P<<endl;
	Dissipative_Flux[0] = 0.5 * (Sensor * Mod_Alpha0 + Alpha_P * d_U[0]);
	Dissipative_Flux[1] = 0.5 * (Sensor * Mod_Alpha1 + Alpha_P * d_U[1]);
	Dissipative_Flux[2] = 0.5 * (Sensor * Mod_Alpha2 + Alpha_P * d_U[2]);
	Dissipative_Flux[3] = 0.5 * (Sensor * Mod_Alpha3 + Alpha_P * d_U[3]);
}

void MOVERS_NWSC_2O(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{

	Rho_L = 0.0, P_L = 0.0, C_L = 0.0, u_L = 0.0, v_L = 0.0, Vdotn_L = 0.0, nx = 0.0, ny = 0.0, Vmag_L = 0.0;
	Rho_R = 0.0, P_R = 0.0, C_R = 0.0, u_R = 0.0, v_R = 0.0, Vdotn_R = 0.0, Vmag_R = 0.0, dl = 0.0;

	Mod_Alpha0 = 0.0, Mod_Alpha1 = 0.0, Mod_Alpha2 = 0.0, Mod_Alpha3 = 0.0, d_F_0 = 0.0, d_F_1 = 0.0, d_F_2 = 0.0;
	d_F_3 = 0.0, d_U_0 = 0.0, d_U_1 = 0.0, d_U_2 = 0.0, d_U_3 = 0.0;
	Lambda_Max = 0.0, Lambda_Min = 0.0, Max1 = 0.0, Max2 = 0.0, Min1 = 0.0, Min2 = 0.0, d_Var = 0.0;

	double Alpha_P = 0.0, P_I = 0.0, d_P = 0.0, Sensor = 0.0, beta = 0.21, a = 0.5;
	// vector<double> d_U(4,0.0);
	int index = Face_No * 2, Which_Var = 0;

	mev_L = 0.0, mev_R = 0.0, max_eigen_value = 0.0;

	V_D S(6, 0.0);
	U_L[0] = 0.0;
	U_L[1] = 0.0;
	U_L[2] = 0.0;
	U_L[3] = 0.0;
	U_R[0] = 0.0;
	U_R[1] = 0.0;
	U_R[2] = 0.0;
	U_R[3] = 0.0;

	Dissipative_Flux[0] = 0.0;
	Dissipative_Flux[1] = 0.0;
	Dissipative_Flux[2] = 0.0;
	Dissipative_Flux[3] = 0.0;

	//      Left state Variables, Density, Pressure, u, v, speed of cound C

	Rho_L = Primitive_Cells[Cell_No][0];
	P_L = Primitive_Cells[Cell_No][4];
	u_L = Primitive_Cells[Cell_No][1];
	v_L = Primitive_Cells[Cell_No][2];
	C_L = Primitive_Cells[Cell_No][5];

	//      Right state Variables, Density, Pressure, u, v, speed of cound C
	Rho_R = Primitive_Cells[N_Cell_No][0];
	P_R = Primitive_Cells[N_Cell_No][4];
	u_R = Primitive_Cells[N_Cell_No][1];
	v_R = Primitive_Cells[N_Cell_No][2];
	C_R = Primitive_Cells[N_Cell_No][5];

	//      Magnitudes of Velocity on Left side and Right side of an interface
	Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
	Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

	//	cout<<u_L<<"\t"<<u_R<<"\t"<<v_L<<"\t"<<v_R<<endl;
	//      nx and ny are normals to the face and dl - face length
	/*nx = Cell_Face_Normals[Cell_No][index + 0];
	ny = Cell_Face_Normals[Cell_No][index + 1];
	dl = Cell_Face_Areas[Cell_No][Face_No];*/
	nx = Cells[Cell_No].Face_Normals[index + 0];
	ny = Cells[Cell_No].Face_Normals[index + 1];
	dl = Cells[Cell_No].Face_Areas[Face_No];

	//	cout<<nx<<"\t"<<ny<<"\t"<<dl<<endl;

	//      Normal Velocity taking Left state Velocity and Right State Velocity of an interface
	Vdotn_L = (u_L * nx + v_L * ny);
	Vdotn_R = (u_R * nx + v_R * ny);

	//          Evaluating Conserved variable Difference between Left and Right faces

	d_U[0] = 0.0;
	d_U[1] = 0.0;
	d_U[2] = 0.0;
	d_U[3] = 0.0;

	//          Evaluating Flux Difference between Left and Right Faces
	d_F_0 = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
	d_F_1 = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L - P_L * nx) * dl;
	d_F_2 = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L - P_L * ny) * dl;
	d_F_3 = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R - (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) * dl;

	Second_Order_Limiter(Cell_No, Face_No, d_U);

	Condition_For_MOVERS_NWSC(d_U[0], d_F_0, Mod_Alpha0);
	Condition_For_MOVERS_NWSC(d_U[1], d_F_1, Mod_Alpha1);
	Condition_For_MOVERS_NWSC(d_U[2], d_F_2, Mod_Alpha2);
	Condition_For_MOVERS_NWSC(d_U[3], d_F_3, Mod_Alpha3);

	d_P = P_R - P_L;
	P_I = 0.5 * (P_R + P_L);
	Alpha_P = 0.5 * (fabs(Vdotn_R * dl) + fabs(Vdotn_L * dl));
	Sensor = beta * fabs(d_P / (2.0 * P_I));

	//			      cout<<Sensor<<"\t"<< Mod_Alpha0<<"\t"<< Mod_Alpha1<<"\t"<<Mod_Alpha2<<"\t"<<Mod_Alpha3<<"\t"<<Alpha_P<<endl;
	Dissipative_Flux[0] = 0.5 * (Sensor * Mod_Alpha0 + Alpha_P * d_U[0]);
	Dissipative_Flux[1] = 0.5 * (Sensor * Mod_Alpha1 + Alpha_P * d_U[1]);
	Dissipative_Flux[2] = 0.5 * (Sensor * Mod_Alpha2 + Alpha_P * d_U[2]);
	Dissipative_Flux[3] = 0.5 * (Sensor * Mod_Alpha3 + Alpha_P * d_U[3]);
}
