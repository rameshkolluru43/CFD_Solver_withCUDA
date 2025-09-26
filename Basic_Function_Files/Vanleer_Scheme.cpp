#include"../Basic_Function_Files/headers.hpp"

//VanLeer Flux Splitting Scheme
void Face::Flux_From_Face_VL(Cell* CC,Cell * NC)
{
	vector<double> Q1(5,0.0),Q2(5,0.0),Q_Avg(5,0.0),Flux_Plus(5,0.0),Flux_Minus(5,0.0);
	double M_L_Plus,M_R_minus,Rho_Plus,Rho_Minus=0.0,M_L=0.0,M_R=0.0,C_L=0.0,C_R=0.0,V_L=0.0,V_R=0.0,Temp=0.0;
	Q1 = CC->Get_QatCellCenter();
	Q2 = NC->Get_QatCellCenter();
	Cal_Primitive(Q1);
	Face_a = sqrt(gamma*R*Face_T);
	V_L = (Face_Velocity*Face_Normal);
	Face_M =V_L/Face_a;
	M_L = Face_M;
	C_L = Face_a;
	M_L_Plus = Evaluate_Mach(Face_M,1);
// 	cout<<M_L<<"\t"<<C_L<<"\t"<<V_L<<endl;
	Rho_Plus = Q1[0]*C_L*0.25*((M_L+1.0)*(M_L+1.0));

	Cal_Primitive(Q2);
	Face_a = sqrt(gamma*R*Face_T);
	V_R = (Face_Velocity*Face_Normal);
	Face_M = V_R/Face_a;
	C_R = Face_a;
	M_R = Face_M;
	M_R_minus = Evaluate_Mach(Face_M,2);
// 	cout<<M_R<<"\t"<<C_R<<"\t"<<V_R<<endl;
	Rho_Minus = -Q2[0]*C_R*0.25*((M_R-1.0)*(M_R-1.0));
	Q_Avg[0]=(Q1[0] + Q2[0])*0.5;
	Q_Avg[1]=(Q1[1] + Q2[1])*0.5;
	Q_Avg[2]=(Q1[2] + Q2[2])*0.5;
	Q_Avg[3]=(Q1[3] + Q2[3])*0.5;
	Q_Avg[4]=(Q1[4] + Q2[4])*0.5;
	Cal_Primitive(Q_Avg);

//Right state Flux denoted by Plus
	Flux_Plus[0] = Rho_Plus;
	Flux_Plus[1] = Rho_Plus*(Face_Normal(1)*((-V_L+2.0*C_L)/gamma) + v1);
	Flux_Plus[2] = Rho_Plus*(Face_Normal(2)*((-V_L+2.0*C_L)/gamma) + v2);
	Flux_Plus[3] = Rho_Plus*(Face_Normal(3)*((-V_L+2.0*C_L)/gamma) + v3);
	Temp = (gamma-1.0)*V_L+2*Face_a;
	Flux_Plus[4] = 0.5*Rho_Plus*((Temp*Temp)/(gamma*gamma-1)+(((Face_Velocity*Face_Velocity)-V_L*V_L)));
// 	cout<<Flux_Plus[0]<<"\t\t"<<Flux_Plus[1]<<"\t\t"<<Flux_Plus[2]<<"\t\t"<<Flux_Plus[3]<<"\t\t"<<Flux_Plus[4]<<"\n";

//Left state Flux denoted by Minus
	Flux_Minus[0] = Rho_Minus;
	Flux_Minus[1] = Rho_Minus*(Face_Normal(1)*((-V_R-2*C_R)/gamma) + v1);
	Flux_Minus[2] = Rho_Minus*(Face_Normal(2)*((-V_R-2*C_R)/gamma) + v2);
	Flux_Minus[3] = Rho_Minus*(Face_Normal(3)*((-V_R-2*C_R)/gamma) + v3);
	Temp = (gamma-1.0)*V_R-2.0*Face_a;
	Flux_Minus[4] =( 0.5*Rho_Minus)*((Temp*Temp)/(gamma*gamma-1)+(((Face_Velocity*Face_Velocity)-V_R*V_R)));
// 	cout<<Flux_Minus[0]<<"\t\t"<<Flux_Minus[1]<<"\t\t"<<Flux_Minus[2]<<"\t\t"<<Flux_Minus[3]<<"\t\t"<<Flux_Minus[4]<<"\n";

	Flux[0] = (Flux_Minus [0] + Flux_Plus[0])*Face_area;
	Flux[1] = (Flux_Minus [1] + Flux_Plus[1])*Face_area;
	Flux[2] = (Flux_Minus [2] + Flux_Plus[2])*Face_area;
	Flux[3] = (Flux_Minus [3] + Flux_Plus[3])*Face_area;
	Flux[4] = (Flux_Minus [4] + Flux_Plus[4])*Face_area;
	Find_Gradients();
// 	cout<<Flux[0]<<"\t\t"<<Flux[1]<<"\t\t"<<Flux[2]<<"\t\t"<<Flux[3]<<"\t\t"<<Flux[4]<<"\n";
// 	 cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

double Face::Evaluate_Mach(const double & M,const int & i)
{
	switch(i)
	{
		case 1:
			if(M>=1.0)
				return (M);
			else if ((M>-1.0 )||( M<1.0))
				return (0.0);
			else
				return (0.25*(M+1.0)*(M+1.0));
		case 2:
			if(M>=1.0)
				return 0.0;
			else if ((M>-1.0 )||( M<1.0))
				return (M);
			else
				return (0.25*(M-1.0)*(M-1.0));
		default:
				cout<<"Evaluating Left and Right states of Mach number so please check value of i, it should be 1 for Left state,2 for Right state\n";
				return (0.0);
	}
}