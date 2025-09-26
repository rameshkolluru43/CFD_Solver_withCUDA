#include "../Basic_Files/definitions.h"
void Cell::Cal_Q_Gradient()
{
	rho_Gradient.Clear();rhou_Gradient.Clear();rhov_Gradient.Clear();rhow_Gradient.Clear();rhoET_Gradient.Clear();
	
	rho_Gradient =(Front_Face.Get_QnDS(0) + Back_Face.Get_QnDS(0) + Top_Face.Get_QnDS(0)
			+ Left_Face.Get_QnDS(0) + Right_Face.Get_QnDS(0));

	rhou_Gradient =(Front_Face.Get_QnDS(1) + Back_Face.Get_QnDS(1) + Top_Face.Get_QnDS(1)
			+ Left_Face.Get_QnDS(1) + Right_Face.Get_QnDS(1));
	
	rhov_Gradient =(Front_Face.Get_QnDS(2) + Back_Face.Get_QnDS(2) + Top_Face.Get_QnDS(2)
			+ Left_Face.Get_QnDS(2) + Right_Face.Get_QnDS(2));
	
	rhow_Gradient =(Front_Face.Get_QnDS(3) + Back_Face.Get_QnDS(3) + Top_Face.Get_QnDS(3)
			+ Left_Face.Get_QnDS(3) + Right_Face.Get_QnDS(3));

	rhoET_Gradient =(Front_Face.Get_QnDS(4) + Back_Face.Get_QnDS(4) + Top_Face.Get_QnDS(4)
			+ Left_Face.Get_QnDS(4) + Right_Face.Get_QnDS(4));
	if(No_of_Faces==6)
	{
		rho_Gradient += Bottom_Face.Get_QnDS(0);
		rhou_Gradient += Bottom_Face.Get_QnDS(1);
		rhov_Gradient += Bottom_Face.Get_QnDS(2);
		rhow_Gradient += Bottom_Face.Get_QnDS(3);
		rhoET_Gradient += Bottom_Face.Get_QnDS(4);
	}
	rho_Gradient *= inv_vol;
	rhou_Gradient *= inv_vol;
	rhov_Gradient *= inv_vol;
	rhow_Gradient *= inv_vol;
	rhoET_Gradient *= inv_vol;
}

const Vector & Cell::Get_Q_Grad(const int & i)
{
	switch(i)
	{
		case 0:
			return rho_Gradient;
		case 1:
			return rhou_Gradient;
		case 2:
			return rhov_Gradient;
		case 3:
			return rhow_Gradient;
		case 4:
			return rhoET_Gradient;
		default :
			cout<<"Enter 0 -4 for fetching Q Gradients, default returning rho Gradient\n";
			return rho_Gradient;
	}
}

const vector<double>& Cell::Get_DelnQ(const int & i, const int & j)
{
	switch(i)
	{
		case 1:
			switch(i)
			{
				case 1:
					return Del2_Q_1;
				case 2:
					return Del2_Q_2;
				case 3:
					return Del2_Q_3;
				default:
					cout<<"For parameter j, enter 1,2 or 3 to get Del2_Q_1,Del2_Q_2 and Del2_Q_3 \n";
					return Del2_Q_1;
			}
		case 2:
			switch(i)
			{
				case 1:
					return Del4_Q_1;
				case 2:
					return Del4_Q_2;
				case 3:
					return Del4_Q_3;
				default:
					cout<<"For parameter j, enter 1,2 or 3 to get Del4_Q_1,Del4_Q_2 and Del4_Q_3 \n";
					return Del4_Q_1;
			}
		case 3:
			switch(i)
			{
				case 1:
					return Del6_Q_1;
				case 2:
					return Del6_Q_2;
				case 3:
					return Del6_Q_3;
				default:
					cout<<"For parameter j, enter 1,2 or 3 to get Del6_Q_1,Del6_Q_2 and Del6_Q_3 \n";
					return Del6_Q_1;
			}
		default:
			cout<<"For parameter i, enter  1,2 or 3 to get Del2_Q,Del4_Q and Del6_Q \n";
			break;
	}
}
