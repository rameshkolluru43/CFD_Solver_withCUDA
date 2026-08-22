#include "../Basic_Files/definitions.h"
//Given Pressure, Temparature and Velocity Vector calculates Q for a given cell
void Cell::Q_atCellCenter(const double & P,const double & T,const Vector& V)
{
	Cell_Pressure = P;
	Cell_Temparature = T;
	Cell_Rho = P/(R*T);
	Cell_Q[0]=Cell_Rho;
	Cell_Q[1]=Cell_Rho*V(1);
	Cell_Q[2]=Cell_Rho*V(2);
	Cell_Q[3]=Cell_Rho*V(3);
	Cell_Et = (cv*T +  0.5*(V*V));
	Cell_Velocity=V;
	Cell_Q[4]=Cell_Rho*Cell_Et;
}
//Function for updating Q 
void Cell::operator()(const vector<double>& Q_temp)
{
	Cell_Q=Q_temp;
	Cal_Primitive(Q_temp);
}

//calculates p,T,v1,v2,v3,rho,a,M when Q is given
void Cell::Cal_Primitive(const vector<double>& Qt)
{
//  	cout<<"In cal Primitive function\n";
//  	cout<<Qt[0]<<"\t"<<Qt[1]<<"\t"<<Qt[2]<<"\t"<<Qt[3]<<"\t"<<Qt[4]<<endl;
	double temp=0.0;
	Cell_Rho=Qt[0];
	temp=1.0/Cell_Rho;
	v1=Qt[1]*temp;
	v2=Qt[2]*temp;
	v3=Qt[3]*temp;
	Cell_Velocity(v1,v2,v3);
	Cell_Et=Qt[4]*temp;
	Cell_Temparature = (Cell_Et-(0.5*(v1*v1+v2*v2+v3*v3)))/cv;
	Cell_Pressure=Cell_Rho*R*Cell_Temparature;
}

//Function called for summation of fluxes from all the faces 0 for convective fluxes and 1 for viscous and Dissipation quantities
void Cell::Add_Fluxes(const vector<double>& T_F,const int & i)
{
	switch(i)
	{
		case 0:
			Cell_Flux[0] += T_F[0];
			Cell_Flux[1] += T_F[1];
			Cell_Flux[2] += T_F[2];
			Cell_Flux[3] += T_F[3];
			Cell_Flux[4] += T_F[4];
			break;
		case 1:
			Cell_Flux[0] -= T_F[0];
			Cell_Flux[1] -= T_F[1];
			Cell_Flux[2] -= T_F[2];
			Cell_Flux[3] -= T_F[3];
			Cell_Flux[4] -= T_F[4];
			break;
		default:
			cout<<"pass either 0 or 1............Please check\n";
			break;
	}	
//	cout<<"Final\t"<<Cell_Flux[0]<<"\t"<<Cell_Flux[1]<<"\t"<<Cell_Flux[2]<<"\t"<<Cell_Flux[3]<<"\t"<<Cell_Flux[4]<<endl;
}


// calculates flux on all the faces
void Cell::Evaluate_FluxonFaces(const int & i)
{
// 	cout<<"in Cell Get_FluxonFaces function\n";
	vector<double> T_F(5,0.0);
	double P1=0.0,P2=0.0,P3=0.0,P4=0.0,P5=0.0,P6=0.0;
	P1 = Front_Cell->Get_Cell_Pressure();
	P2 = Back_Cell->Get_Cell_Pressure();
	P3 = Top_Cell->Get_Cell_Pressure();
	P5 = Left_Cell->Get_Cell_Pressure();
	P4 = Right_Cell->Get_Cell_Pressure();
	
	lambda = fabs((P1+P2+P3+P4+P5 - 5.0*Cell_Pressure)/(P1+P2+P3+P4+P5 + 5.0*Cell_Pressure));
	
	Cell_Flux[0]=0.0;Cell_Flux[1]=0.0;Cell_Flux[2]=0.0;Cell_Flux[3]=0.0;Cell_Flux[4]=0.0;

	T_F=Front_Face.Flux_From_Face(this,Front_Cell,i);
	Add_Fluxes(T_F,0);
	T_F=Back_Face.Flux_From_Face(this,Back_Cell,i);
	Add_Fluxes(T_F,0);
	T_F=Left_Face.Flux_From_Face(this, Left_Cell,i);
	Add_Fluxes(T_F,0);
	T_F=Right_Face.Flux_From_Face(this,Right_Cell,i);
	Add_Fluxes(T_F,0);
	T_F=Top_Face.Flux_From_Face(this, Top_Cell,i);
	Add_Fluxes(T_F,0);
	if(No_of_Faces==6)
	{
		T_F=Bottom_Face.Flux_From_Face(this, Bottom_Cell,i);
		P6 = Bottom_Cell->Get_Cell_Pressure();
		lambda = fabs((P1+P2+P3+P4+P5+P6 - 6.0*Cell_Pressure)/(P1+P2+P3+P4+P5+P6 + 5.0*Cell_Pressure));
		Add_Fluxes(T_F,0);
	}
	Velocity_GradatCenter();
	Cal_Q_Gradient();
/*
	cout<<self_index<<"\t"<<Cell_Flux[0]<<"\t"<<Cell_Flux[1]<<"\t"<<Cell_Flux[2]<<"\t"<<Cell_Flux[3]<<"\t"<<Cell_Flux[4]<<endl;
	cout<<"-----------------------------------------------------------------\n";*/
}

const vector<double>& Cell::Get_FluxofCell() const
{
	return Cell_Flux;
}

void  Cell::Set_FluxofCell(vector<double>& Flux) 
{
	Cell_Flux=Flux;
}


vector<double>& Cell::Get_QatCellCenter()
{
// 	cout<<"Cell index\t"<<self_index<<endl;
	return Cell_Q;
}

const Vector& Cell::Get_Cell_Velocity() const
{
	return Cell_Velocity;
}

void Cell::Set_Cell_Velocity(const Vector & V)
{
	Cell_Velocity=V;
}

const double& Cell::Get_Cell_Pressure() const
{
	return Cell_Pressure;
}

const double& Cell::Get_Cell_Density() const
{
	return Cell_Rho;
}

const double& Cell::Get_Cell_Temparature() const
{
	return Cell_Temparature;
}

void Cell::UpdateQ()
{
// 	cout<<"Updating Cell Q\n";
	Cell_Q = Q_new;
	Cal_Primitive(Q_new);
}


void Cell::Test_Fluxes(const int & j)
{
	vector<double> Temp_Flux(5,0.0),Temp_Flux1(5,0.0);
	Face Temp_Face;
	ofstream fluxfile("../Euler_Cube_Results/Fluxfile.txt",ios::out|ios::app);
	if(fluxfile.is_open())
	{
		fluxfile<<"Current Cell Number\t"<<self_index<<endl;
		// 		cout<<"Current Cell Number\t"<<self_index<<endl;
		fluxfile<<"Cell_Flux Back face\t";
		Temp_Flux=Back_Face.Flux_From_Face(this, Back_Cell,j);
		Temp_Face =Back_Cell->Get_Face(0);
		Temp_Flux1=Temp_Face.Flux_From_Face(Back_Cell, this,j);
		for(int i=0;i<5;i++)
		{
			fluxfile<<Temp_Flux.at(i)+Temp_Flux1.at(i)<<"\t";
		}
		fluxfile<<endl;
		Temp_Flux=Front_Face.Flux_From_Face(this, Front_Cell,j);
		Temp_Face=Front_Cell->Get_Face(1);
		Temp_Flux1=Temp_Face.Flux_From_Face(Front_Cell, this,j);
		fluxfile<<"Cell_Flux Front face\t";
		for(int i=0;i<5;i++)
		{
			fluxfile<<Temp_Flux.at(i)+Temp_Flux1.at(i)<<"\t";
		}
		fluxfile<<endl;
		Temp_Flux=Left_Face.Flux_From_Face(this, Left_Cell,j);
		Temp_Face =Left_Cell->Get_Face(3);
		Temp_Flux1=Temp_Face.Flux_From_Face(Left_Cell,this,j);
		fluxfile<<"Cell_Flux Left\t";
		for(int i=0;i<5;i++)
		{
			fluxfile<<Temp_Flux.at(i)+Temp_Flux1.at(i)<<"\t";
		}
		fluxfile<<endl;
		Temp_Flux=Right_Face.Flux_From_Face(this, Right_Cell,j);
		Temp_Face =Right_Cell->Get_Face(2);
		Temp_Flux1=Temp_Face.Flux_From_Face(Right_Cell,this,j);
		fluxfile<<"Cell_Flux Right Face\t";
		for(int i=0;i<5;i++)
		{
			fluxfile<<Temp_Flux.at(i)+Temp_Flux1.at(i)<<"\t";
		}
		fluxfile<<endl;
		Temp_Flux=Top_Face.Flux_From_Face(this, Top_Cell,j);
		Temp_Face =Top_Cell->Get_Face(5);
		Temp_Flux1=Temp_Face.Flux_From_Face(Top_Cell,this,j);
		fluxfile<<"Cell_Flux Top face\t";
		for(int i=0;i<5;i++)
			fluxfile<<Temp_Flux.at(i)+Temp_Flux1.at(i)<<"\t";
		fluxfile<<endl;
		if(No_of_Faces==6)
		{
			Temp_Flux=Bottom_Face.Flux_From_Face(this, Bottom_Cell,j);
			Temp_Face =Bottom_Cell->Get_Face(4);
			Temp_Flux1=Temp_Face.Flux_From_Face(Bottom_Cell,this,j);
			fluxfile<<"Cell_Flux Bottom face\t";
			for(int i=0;i<5;i++)
				fluxfile<<Temp_Flux.at(i)+Temp_Flux1.at(i)<<"\t";
			fluxfile<<endl;
		}
		fluxfile<<"-----------------------------------------------------------\n";
	}
	else
	cout<<"Unable to open the file\n";
}
