#include "../Basic_Function_Files/headers.hpp"
// Implementing Range-Kutta 4th order method for a given cell
vector<double>& Cell::Rk4(const double& time_step,const int & i)
{
	vector<double> Flux0(5,0.0),Flux1(5,0.0),Flux2(5,0.0),Temp_Q(5,0.0),D(5,0.0);
	double deltvol=0.0,eps1=0.0,eps2=0.0,dx=0.0,dy=0.0,dz=0.0;
	if(i==3)
	{
		eps1 = 0.25*lambda; eps2 =1.0/256.0;
		Diagonal_Vector.Print();
		dx = Diagonal_Vector(1); dy = Diagonal_Vector(2); dz = Diagonal_Vector(3);
		dx = dy=dz = Cell_Avg_Length;
 		cout<<self_index<<"\t"<<lambda<<"\t"<<dx<<"\t"<<dy<<"\t"<<dz<<"\t"<<time_step<<endl;
		D[0]		=	time_step*(eps1*(dx*dx*Del2_Q_1[0] + dy*dy*Del2_Q_2[0]+dz*dz*Del2_Q_3[0])	-
					(eps2-eps1)*(dx*dx*dx*dx*Del4_Q_1[0] + dy*dy*dy*dy*Del4_Q_2[0]+dz*dz*dz*dz*Del4_Q_3[0]));
		D[1]		=	time_step*(eps1*(dx*dx*Del2_Q_1[1] + dy*dy*Del2_Q_2[1]+dz*dz*Del2_Q_3[1])	-
					(eps2-eps1)*(dx*dx*dx*dx*Del4_Q_1[1] + dy*dy*dy*dy*Del4_Q_2[1]+dz*dz*dz*dz*Del4_Q_3[1]));
		D[2]		=	time_step*(eps1*(dx*dx*Del2_Q_1[2] + dy*dy*Del2_Q_2[2]+dz*dz*Del2_Q_3[2])	-
					(eps2-eps1)*(dx*dx*dx*dx*Del4_Q_1[2] + dy*dy*dy*dy*Del4_Q_2[2]+dz*dz*dz*dz*Del4_Q_3[2]));
		D[3]		=	time_step*(eps1*(dx*dx*Del2_Q_1[3] + dy*dy*Del2_Q_2[3]+dz*dz*Del2_Q_3[3])	-
					(eps2-eps1)*(dx*dx*dx*dx*Del4_Q_1[3] + dy*dy*dy*dy*Del4_Q_2[3]+dz*dz*dz*dz*Del4_Q_3[3]));
		D[4]		=	time_step*(eps1*(dx*dx*Del2_Q_1[4] + dy*dy*Del2_Q_2[4]+dz*dz*Del2_Q_3[4])	-
					(eps2-eps1)*(dx*dx*dx*dx*Del4_Q_1[4] + dy*dy*dy*dy*Del4_Q_2[4]+dz*dz*dz*dz*Del4_Q_3[4]));
	}
	else
	{eps1=0.0;eps2=0.0;}

	Temp_Q = Cell_Q;
	Flux0 = Cell_Flux;			//k1
	deltvol	=	0.5*time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * Cell_Flux[0]-D[0];
	Q_new[1]= Temp_Q[1] - deltvol * Cell_Flux[1]-D[1];
	Q_new[2]= Temp_Q[2] - deltvol * Cell_Flux[2]-D[2];
	Q_new[3]= Temp_Q[3] - deltvol * Cell_Flux[3]-D[3];
	Q_new[4]= Temp_Q[4] - deltvol * Cell_Flux[4]-D[4];
	UpdateQ();
	Evaluate_FluxonFaces(i);		// Cell_Flux consists of flux of updated Q-->Q2
//  	Cal_Stresses();
	Flux1 = Cell_Flux;			//k1
	deltvol	=	0.5*time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * Cell_Flux[0]-D[0];
	Q_new[1]= Temp_Q[1] - deltvol * Cell_Flux[1]-D[1];
	Q_new[2]= Temp_Q[2] - deltvol * Cell_Flux[2]-D[2];
	Q_new[3]= Temp_Q[3] - deltvol * Cell_Flux[3]-D[3];
	Q_new[4]= Temp_Q[4] - deltvol * Cell_Flux[4]-D[4];
	UpdateQ();
	Evaluate_FluxonFaces(i);		// Cell_Flux consists of flux of updated Q-->Q3
//  	Cal_Stresses();
	Flux2 = Cell_Flux;			//k1
	deltvol	=	time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * Cell_Flux[0]-D[0];
	Q_new[1]= Temp_Q[1] - deltvol * Cell_Flux[1]-D[1];
	Q_new[2]= Temp_Q[2] - deltvol * Cell_Flux[2]-D[2];
	Q_new[3]= Temp_Q[3] - deltvol * Cell_Flux[3]-D[3];
	Q_new[4]= Temp_Q[4] - deltvol * Cell_Flux[4]-D[4];
	UpdateQ();
	Evaluate_FluxonFaces(i);			// Cell_Flux consists of flux of updated Q-->Q4
// 	Cal_Stresses();
	deltvol=(1.0/6.0)*time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * (Flux0[0]+2.0*Flux1[0]+2.0*Flux2[0]+Cell_Flux[0])-D[0];
	Q_new[1]= Temp_Q[1] - deltvol * (Flux0[1]+2.0*Flux1[1]+2.0*Flux2[1]+Cell_Flux[1])-D[1];
	Q_new[2]= Temp_Q[2] - deltvol * (Flux0[2]+2.0*Flux1[2]+2.0*Flux2[2]+Cell_Flux[2])-D[2];
	Q_new[3]= Temp_Q[3] - deltvol * (Flux0[3]+2.0*Flux1[3]+2.0*Flux2[3]+Cell_Flux[3])-D[3];
	Q_new[4]= Temp_Q[4] - deltvol * (Flux0[4]+2.0*Flux1[4]+2.0*Flux2[4]+Cell_Flux[4])-D[4];
	return Q_new;
}

void Cell_Property_Manager::Rk4(const double& time_step,const int &index)
{
	vector<double> Flux0(5,0.0),Flux1(5,0.0),Flux2(5,0.0),Flux3(5,0.0),Temp_Q(5,0.0),Fdisp(5,0.0),Fdisp_1(5,0.0),Fdisp_2(5,0.0),Fdisp_3(5,0.0);
	double deltvol=0.0,inv_vol=0.0,dx=0.0,dy=0.0,dz=0.0,mu6 = -0.0078125,inv_deltvol=0.0;
	Vector Diag_Vector = Cell_List_h.at(index)->Get_Diagonal_Vector();
	Temp_Q = Cell_List_h.at(index)->Get_QatCellCenter();
	Flux0 = Cell_List_h.at(index)->Get_FluxofCell();			//k1
	Fdisp_1 = Cell_List_h.at(index)->Get_DelnQ(3,1);
	Fdisp_2 = Cell_List_h.at(index)->Get_DelnQ(3,2);
	Fdisp_3 = Cell_List_h.at(index)->Get_DelnQ(3,3);
	dx = Diag_Vector(1); dy = Diag_Vector(2); dz = Diag_Vector(3);
	for(int i=0;i<5;i++)
		Fdisp[i] = (dx*dx*dx*dx*dx*dx)*Fdisp_1[i]+(dy*dy*dy*dy*dy*dy)*Fdisp_2[i]+(dz*dz*dz*dz*dz*dz)*Fdisp_3[i];
	inv_vol = Cell_List_h.at(index)->Get_Volume();
	deltvol=time_step*inv_vol;
	inv_deltvol = mu6 / deltvol;
	Q_new[0]= Temp_Q[0] - deltvol * Flux0[0] + inv_deltvol*Fdisp[0];
	Q_new[1]= Temp_Q[1] - deltvol * Flux0[1] + inv_deltvol*Fdisp[1];
	Q_new[2]= Temp_Q[2] - deltvol * Flux0[2] + inv_deltvol*Fdisp[2];
	Q_new[3]= Temp_Q[3] - deltvol * Flux0[3] + inv_deltvol*Fdisp[3];
	Q_new[4]= Temp_Q[4] - deltvol * Flux0[4] + inv_deltvol*Fdisp[4];
	(*Cell_List_h.at(index))(Q_new);
	(*Cell_List_3h.at(index))(Q_new);
	Cal_QAverageonFaces(index);
	Flux1 = Cell_List_h.at(index)->Get_FluxofCell();			//k2
	deltvol=0.5*time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * Flux1[0] + inv_deltvol*Fdisp[0];
	Q_new[1]= Temp_Q[1] - deltvol * Flux1[1] + inv_deltvol* Fdisp[1];
	Q_new[2]= Temp_Q[2] - deltvol * Flux1[2] + inv_deltvol*Fdisp[2];
	Q_new[3]= Temp_Q[3] - deltvol * Flux1[3] + inv_deltvol*Fdisp[3];
	Q_new[4]= Temp_Q[4] - deltvol * Flux1[4] + inv_deltvol*Fdisp[4];
	(*Cell_List_h.at(index))(Q_new);
	(*Cell_List_3h.at(index))(Q_new);
	Cal_QAverageonFaces(index);
	Flux2 = Cell_List_h.at(index)->Get_FluxofCell();			//k3
	deltvol=0.5*time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * Flux2[0] + inv_deltvol*Fdisp[0];
	Q_new[1]= Temp_Q[1] - deltvol * Flux2[1] + inv_deltvol*Fdisp[1];
	Q_new[2]= Temp_Q[2] - deltvol * Flux2[2] + inv_deltvol*Fdisp[2];
	Q_new[3]= Temp_Q[3] - deltvol * Flux2[3] + inv_deltvol*Fdisp[3];
	Q_new[4]= Temp_Q[4] - deltvol * Flux2[4] + inv_deltvol*Fdisp[4];
	(*Cell_List_h.at(index))(Q_new);
	(*Cell_List_3h.at(index))(Q_new);
	Cal_QAverageonFaces(index);
	Flux3 = Cell_List_h.at(index)->Get_FluxofCell();			//k3
	deltvol=(1.0/6.0)*time_step*inv_vol;
	Q_new[0]= Temp_Q[0] - deltvol * (Flux0[0]+2.0*Flux1[0]+2.0*Flux2[0]+Flux3[0]) + inv_deltvol*Fdisp[0];
	Q_new[1]= Temp_Q[1] - deltvol * (Flux0[1]+2.0*Flux1[1]+2.0*Flux2[1]+Flux3[1]) + inv_deltvol*Fdisp[1];
	Q_new[2]= Temp_Q[2] - deltvol * (Flux0[2]+2.0*Flux1[2]+2.0*Flux2[2]+Flux3[2]) + inv_deltvol*Fdisp[2];
	Q_new[3]= Temp_Q[3] - deltvol * (Flux0[3]+2.0*Flux1[3]+2.0*Flux2[3]+Flux3[3]) + inv_deltvol*Fdisp[3];
	Q_new[4]= Temp_Q[4] - deltvol * (Flux0[4]+2.0*Flux1[4]+2.0*Flux2[4]+Flux3[4]) + inv_deltvol*Fdisp[4];
}
