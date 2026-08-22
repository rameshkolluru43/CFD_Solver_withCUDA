#include "headers.hpp"

const vector<double>& Face::Add_Dissipation(Cell * Cell_i, Cell * Cell_i_plus_1,Cell * Cell_i_minus_1,Cell *Cell_i_plus_2)
{
	vector<double> Q_i(5,0.0),Q_i_plus_1(5,0.0),Q_i_plus_2(5,0.0),Q_i_minus_1(5,0.0),Q_i_minus_2(5,0.0);
	bool flag=false;
	// alpha2 = 1/4; alpha4 = 1/256
	double epsilon4=0.0,epsilon2=0.0;
	Flux[0]=0.0;Flux[1]=0.0;Flux[2]=0.0;Flux[3]=0.0;Flux[4]=0.0;
	flag = Cell_i_plus_1->Is_Ghost_Cell();
	Q_i = Cell_i->Get_QatCellCenter();
	Q_i_plus_1 = Cell_i_plus_1->Get_QatCellCenter();
	Q_i_minus_1= Cell_i_minus_1->Get_QatCellCenter();
	Q_i_plus_2 = Cell_i_plus_2->Get_QatCellCenter();
	epsilon2 = 0.5*max_eigen_value*Face_area;
	epsilon4 = 0.0625*max_eigen_value*Face_area*Face_area*Face_area;
	switch(flag)
	{
		case true:
			Flux[0] =
			(epsilon2*(Q_i_plus_1[0]-Q_i[0])-epsilon4*(-Q_i_plus_2[0]+Q_i_plus_1[0]-3*Q_i[0]-3*Q_i_minus_1[0]));
			Flux[1] =
			(epsilon2*(Q_i_plus_1[1]-Q_i[1])-epsilon4*(-Q_i_plus_2[1]+Q_i_plus_1[1]-3*Q_i[1]-3*Q_i_minus_1[1]));
			Flux[2] =
			(epsilon2*(Q_i_plus_1[2]-Q_i[2])-epsilon4*(-Q_i_plus_2[2]+Q_i_plus_1[2]-3*Q_i[2]-3*Q_i_minus_1[2]));
			Flux[3] =
			(epsilon2*(Q_i_plus_1[3]-Q_i[3])-epsilon4*(-Q_i_plus_2[3]+Q_i_plus_1[3]-3*Q_i[3]-3*Q_i_minus_1[3]));
			Flux[4] =
			(epsilon2*(Q_i_plus_1[4]-Q_i[4])-epsilon4*(-Q_i_plus_2[4]+Q_i_plus_1[4]-3*Q_i[4]-3*Q_i_minus_1[4]));
			flag==false;
			break;
		default:
			Flux[0] =
			(epsilon2*(Q_i_plus_1[0]-Q_i[0])-epsilon4*(Q_i_plus_2[0]-3*Q_i_plus_1[0]+3*Q_i[0]-Q_i_minus_1[0]));
			Flux[1] =
			(epsilon2*(Q_i_plus_1[1]-Q_i[1])-epsilon4*(Q_i_plus_2[1]-3*Q_i_plus_1[1]+3*Q_i[1]-Q_i_minus_1[1]));
			Flux[2] =
			(epsilon2*(Q_i_plus_1[2]-Q_i[2])-epsilon4*(Q_i_plus_2[2]-3*Q_i_plus_1[2]+3*Q_i[2]-Q_i_minus_1[2]));
			Flux[3] =
			(epsilon2*(Q_i_plus_1[3]-Q_i[3])-epsilon4*(Q_i_plus_2[3]-3*Q_i_plus_1[3]+3*Q_i[3]-Q_i_minus_1[3]));
			Flux[4] =
			(epsilon2*(Q_i_plus_1[4]-Q_i[4])-epsilon4*(Q_i_plus_2[4]-3*Q_i_plus_1[4]+3*Q_i[4]-Q_i_minus_1[4]));
			flag==false;
			break;
	}
return Flux;
}

// Cell_ip2 => cell[i+2], cell_im1=>Cell[i-1], Neighbour_Cell = Cell[i+1], Current_Cell cell[i]
void Face::Flux_From_Face_Central(Cell* Current_Cell, Cell * Neighbour_Cell)
{
	vector<double> Q1(5,0.0),Q2(5,0.0),Q_Avg(5,0.0);
	double vnds=0.0;
	Flux[0]=0.0;Flux[1]=0.0;Flux[2]=0.0;Flux[3]=0.0;Flux[4]=0.0;
	Q1 = Current_Cell->Get_QatCellCenter();	Q2 = Neighbour_Cell->Get_QatCellCenter();
	
	Q_Avg[0]= (Q1[0]+Q2[0])*0.5;
	Q_Avg[1]= (Q1[1]+Q2[1])*0.5;
	Q_Avg[2]= (Q1[2]+Q2[2])*0.5;
	Q_Avg[3]= (Q1[3]+Q2[3])*0.5;
	Q_Avg[4]= (Q1[4]+Q2[4])*0.5;

	Cal_Primitive(Q_Avg);
/*	Face_Velocity.Print();	Area_Vector.Print();	cout<<"v.nds\t"<<vnds<<endl;*/
	if((Neighbour_Cell->Is_Ghost_Cell())&&(Face_No==4))
	{
		Face_P = Current_Cell->Get_Cell_Pressure();
		Face_T = Current_Cell->Get_Cell_Temparature();
	}
	vnds=Face_Velocity*Area_Vector;

//	Flux Determines the Convective Fluxes on the face
	Flux[0] = Face_Rho*vnds;
	Flux[1] = Face_Rho*v1*vnds + Face_P*Area_Comp1;
	Flux[2] = Face_Rho*v2*vnds + Face_P*Area_Comp2;
	Flux[3] = Face_Rho*v3*vnds + Face_P*Area_Comp3;
	Flux[4] = (Face_Rho*Face_ET + Face_P)*vnds;
//   	cout<<Flux[0]<<"\t\t"<<Flux[1]<<"\t\t"<<Flux[2]<<"\t\t"<<Flux[3]<<"\t\t"<<Flux[4]<<"\n";
	
	undS=Scalar_dot_ndS(v1);
	vndS=Scalar_dot_ndS(v2);
	wndS=Scalar_dot_ndS(v3);
	TndS=Scalar_dot_ndS(Face_T);
	
	rho_ndS = Scalar_dot_ndS(Face_Rho);
	rhou_ndS = Scalar_dot_ndS(Face_Rho*v1);
	rhov_ndS = Scalar_dot_ndS(Face_Rho*v2);
	rhow_ndS = Scalar_dot_ndS(Face_Rho*v3);
	rhoET_ndS = Scalar_dot_ndS(Face_Rho*Face_ET);
	
// 	rho_ndS.Print();rhou_ndS.Print();rhov_ndS.Print();rhow_ndS.Print();rhoET_ndS.Print();
// 	cout<<"---------------------------------------------------------------------\n";
	Viscosity();
	Thermal_Conductivity();
}


const Vector& Face::Scalar_dot_ndS(const double & phi)
{
	return (Area_Vector*phi);
}


double Face::Gradient_dot_ndS(const Vector & grad_phi)
{
	return (Area_Vector*grad_phi);
}


const Vector& Face::Get_QnDS(const int & i)
{
	switch(i)
	{
		case 0:
			return rho_ndS;
		case 1:
			return rhou_ndS;
		case 2:
			return rhov_ndS;
		case 3:
			return rhow_ndS;
		case 4:
			return rhoET_ndS;
		default :
			cout<<"Enter 0 -4 for fetching Q Gradients, default returning rho Gradient\n";
			return rho_ndS;
	}
}
