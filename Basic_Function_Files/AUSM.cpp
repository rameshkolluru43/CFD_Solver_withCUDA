#include"../Basic_Function_Files/headers.hpp"

const vector<double> & Face::Flux_From_Face(Cell* CC,Cell * NC,const int & i)
{
	switch(i)
	{
		case 1:
				Flux_From_Face_AUSM(CC, NC);
				return Flux;
		case 2:
				Flux_From_Face_VL(CC, NC);
				return Flux;
		default:
				Flux_From_Face_Central(CC, NC);
				return Flux;
	}
}

void Face::Flux_From_Face_AUSM(Cell* CC,Cell * NC)
{
// 	cout<<"In Face AUSM Function\n";
	vector<double> Q1(5,0.0),Q2(5,0.0),Q_Avg(5,0.0),Flux_Plus(5,0.0),Flux_Minus(5,0.0);
	double M_L_Plus=0.0,M_R_Minus=0.0,V_L=0.0,V_R=0.0,P_Plus=0.0,P_Minus=0.0;

	Q1 = CC->Get_QatCellCenter();
	Q2 = NC->Get_QatCellCenter();
	Cal_Primitive(Q1);
	Face_a = sqrt(gamma_R*Face_T);
	V_L = (Face_Velocity*Face_Normal);
	Face_M = V_L/Face_a;
	M_L_Plus = Evaluate_Mach(Face_M,1);
	P_Plus = Evaluate_Pressure(Face_P,Face_M,1);
// 	cout<<Face_a<<"\t"<<V_L<<"\t"<<Face_M<<"\t"<<M_L_Plus<<"\t"<<P_Plus<<endl;
	Flux_Plus[0] = Face_Rho*Face_a;
	Flux_Plus[1] = Face_Rho*Face_a*v1;
	Flux_Plus[2] = Face_Rho*Face_a*v2;
	Flux_Plus[3] = Face_Rho*Face_a*v3;
	Flux_Plus[4] = (Face_Rho*Face_ET+Face_P)*Face_a;
// 	cout<<"Plus\t"<<Flux_Plus[0]<<"\t\t"<<Flux_Plus[1]<<"\t\t"<<Flux_Plus[2]<<"\t\t"<<Flux_Plus[3]<<"\t\t"<<Flux_Plus[4]<<"\n";

	Cal_Primitive(Q2);
	Face_a = sqrt(gamma_R*Face_T);
	V_R = (Face_Velocity*Face_Normal);
	Face_M = V_R/Face_a;
	M_R_Minus = Evaluate_Mach(Face_M,2);
	P_Minus = Evaluate_Pressure(Face_P,Face_M,2);
//  	cout<<Face_a<<"\t"<<V_R<<"\t"<<Face_M<<"\t"<<M_R_Minus<<"\t"<<P_Minus<<endl;
	Flux_Minus[0] = Face_Rho*Face_a;
	Flux_Minus[1] = Face_Rho*Face_a*v1;
	Flux_Minus[2] = Face_Rho*Face_a*v2;
	Flux_Minus[3] = Face_Rho*Face_a*v3;
	Flux_Minus[4] = (Face_Rho*Face_ET+Face_P)*Face_a;
// 	cout<<"Minus\t"<<Flux_Minus[0]<<"\t\t"<<Flux_Minus[1]<<"\t\t"<<Flux_Minus[2]<<"\t\t"<<Flux_Minus[3]<<"\t\t"<<Flux_Minus[4]<<"\n"; 
	
	Flux[0]=0.0;Flux[1]=0.0;Flux[2]=0.0;Flux[3]=0.0;Flux[4]=0.0;
	Flux[0] = (Face_area*0.5)*((M_L_Plus+M_R_Minus)*(Flux_Minus [0] + Flux_Plus[0])
						- fabs((M_L_Plus+M_R_Minus))*(-Flux_Plus [0] + Flux_Minus[0]));
						
	Flux[1] = (Face_area*0.5)*((M_L_Plus+M_R_Minus)*(Flux_Minus [1]+ Flux_Plus[1])
						- fabs((M_L_Plus+M_R_Minus))*(-Flux_Plus[1]+Flux_Minus [1]))+(P_Plus+P_Minus)*Area_Comp1;
						
	Flux[2] = (Face_area*0.5)*((M_L_Plus+M_R_Minus)*(Flux_Minus [2] + Flux_Plus[2])
						- fabs((M_L_Plus+M_R_Minus))*(-Flux_Plus[2]+Flux_Minus [2]))+(P_Plus+P_Minus)*Area_Comp2;
						
	Flux[3] = (Face_area*0.5)*((M_L_Plus+M_R_Minus)*(Flux_Minus [3] + Flux_Plus[3])
						- fabs((M_L_Plus+M_R_Minus))*(-Flux_Plus[3]+Flux_Minus [3]))+(P_Plus+P_Minus)*Area_Comp3;
						 
	Flux[4] = (Face_area*0.5)*((M_L_Plus+M_R_Minus)*(Flux_Minus [4] + Flux_Plus[4])
							- fabs((M_L_Plus+M_R_Minus))*(- Flux_Plus[4]+Flux_Minus [4]));
// 	cout<<"Net flux\t"<<Flux[0]<<"\t\t"<<Flux[1]<<"\t\t"<<Flux[2]<<"\t\t"<<Flux[3]<<"\t\t"<<Flux[4]<<"\n";
// 	cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

	Q_Avg[0]=(Q1[0] + Q2[0])*0.5;
	Q_Avg[1]=(Q1[1] + Q2[1])*0.5;
	Q_Avg[2]=(Q1[2] + Q2[2])*0.5;
	Q_Avg[3]=(Q1[3] + Q2[3])*0.5;
	Q_Avg[4]=(Q1[4] + Q2[4])*0.5;
	Cal_Primitive(Q_Avg);
	Find_Gradients();
}

double Face::Evaluate_Pressure(const double & P,const double & M,const int & i)
{
	switch(i)
	{
		case 1:
			//Evaluating P_Plus
			if(M>=1.0)
				return P;
			else if ((M>-1.0 )||( M<1.0))
				return ((P*0.25)*(M+1.0)*(M+1.0)*(2.0-M));
			else
				return 0.0;
		case 2:
			//Evaluating P_Minus
			if(M>=1.0)
				return 0.0;
			else if ((M>-1.0 )||( M<1.0))
				return ((P*0.25)*(M-1.0)*(M-1.0)*(2.0+M));
			else
				return P;
		default:
			cout<<"Evaluating Left and Right states of Mach number so please check value of i, it should be 1 for Left state,2 for Right state\n";
			break;
	}
}