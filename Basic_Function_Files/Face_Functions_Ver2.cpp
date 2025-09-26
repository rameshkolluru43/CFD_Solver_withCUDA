#include "../Basic_Function_Files/headers.hpp"

//Gets the first gradient and calculates Del2_Q_1,Del2_Q_2,Del2_Q_3
void Face::Cal_Del2Q(Cell * Current_Cell,Cell* Neighbour_Cell)
{
	Del_Q_1[0]=0.0;Del_Q_1[1]=0.0;Del_Q_1[2]=0.0;Del_Q_1[3]=0.0;Del_Q_1[4]=0.0;
	Del_Q_2[0]=0.0;Del_Q_2[1]=0.0;Del_Q_2[2]=0.0;Del_Q_2[3]=0.0;Del_Q_2[4]=0.0;
	Del_Q_3[0]=0.0;Del_Q_3[1]=0.0;Del_Q_3[2]=0.0;Del_Q_3[3]=0.0;Del_Q_3[4]=0.0;

	if(Neighbour_Cell->Is_Ghost_Cell())
	{
		rho_ndS = Current_Cell->Get_Q_Grad(0);
		rhou_ndS = Current_Cell->Get_Q_Grad(1);
		rhov_ndS = Current_Cell->Get_Q_Grad(2);
		rhow_ndS = Current_Cell->Get_Q_Grad(3);
		rhoET_ndS = Current_Cell->Get_Q_Grad(4);
	}
	else
	{
		rho_ndS = (Current_Cell->Get_Q_Grad(0) + Neighbour_Cell->Get_Q_Grad(0))*0.5;
		rhou_ndS = (Current_Cell->Get_Q_Grad(1) + Neighbour_Cell->Get_Q_Grad(1))*0.5;
		rhov_ndS = (Current_Cell->Get_Q_Grad(2) + Neighbour_Cell->Get_Q_Grad(2))*0.5;
		rhow_ndS = (Current_Cell->Get_Q_Grad(3) + Neighbour_Cell->Get_Q_Grad(3))*0.5;
		rhoET_ndS = (Current_Cell->Get_Q_Grad(4) + Neighbour_Cell->Get_Q_Grad(4))*0.5;
	}
	
	Del_Q_1[0] = rho_ndS(1)*Area_Comp1;
	Del_Q_2[0] = rho_ndS(2)*Area_Comp2;
	Del_Q_3[0] = rho_ndS(3)*Area_Comp3;

	Del_Q_1[1] = rhou_ndS(1)*Area_Comp1;
	Del_Q_2[1] = rhou_ndS(2)*Area_Comp2;
	Del_Q_3[1] = rhou_ndS(3)*Area_Comp3;

	Del_Q_1[2] = rhov_ndS(1)*Area_Comp1;
	Del_Q_2[2] = rhov_ndS(2)*Area_Comp2;
	Del_Q_3[2] = rhov_ndS(3)*Area_Comp3;

	Del_Q_1[3] = rhow_ndS(1)*Area_Comp1;
	Del_Q_2[3] = rhow_ndS(2)*Area_Comp2;
	Del_Q_3[3] = rhow_ndS(3)*Area_Comp3;

	Del_Q_1[4] = rhoET_ndS(1)*Area_Comp1;
	Del_Q_2[4] = rhoET_ndS(2)*Area_Comp2;
	Del_Q_3[4] = rhoET_ndS(3)*Area_Comp3;
}


// this function returns a  stl vector which calculates the nth gradient with ith component
const vector<double>& Face::Calculate_Grad_n_Q(const Cell * const CC, const Cell * const NC, const int &n, const int & i)
{
	Del_Q_1[0]=0.0;Del_Q_1[1]=0.0;Del_Q_1[2]=0.0;Del_Q_1[3]=0.0;Del_Q_1[4]=0.0;
	Del_Q_2[0]=0.0;Del_Q_2[1]=0.0;Del_Q_2[2]=0.0;Del_Q_2[3]=0.0;Del_Q_2[4]=0.0;
	Del_Q_3[0]=0.0;Del_Q_3[1]=0.0;Del_Q_3[2]=0.0;Del_Q_3[3]=0.0;Del_Q_3[4]=0.0;

	vector<double> gradq1(5,0.0),grad_avg(5,0.0),gradq2(5,0.0);
	// n refers for the nth gradient and i refers to the direction
	gradq1 = CC->Get_DelQ_fromCell(n,i);
	if(NC->Is_Ghost_Cell())
		gradq2 = gradq1;
	else
		gradq2 = NC->Get_DelQ_fromCell(n,i);
	
	for(int k=0;k<5;k++)
		grad_avg[k] = 0.5*(gradq1[k]+gradq2[k]);
	
	switch(i)
	{
		case 1:
				for(int p=0;p<5;p++)
					Del_Q_1[p] = grad_avg[p]*Area_Comp1;
				return Del_Q_1;
		case 2:
				for(int p=0;p<5;p++)
					Del_Q_2[p] = grad_avg[p]*Area_Comp2;
				return Del_Q_2;
		case 3:
				for(int p=0;p<5;p++)
					Del_Q_3[p] = grad_avg[p]*Area_Comp3;
				return Del_Q_3;
		default:
				cout<<"Check  of i , should be 1,2,3 to get the corresponding gradients\n";
				break;
	}
}

const vector<double>& Face::Get_DelQ_onFace(const int & i)
{
	switch(i)
	{
		case 1:
			return Del_Q_1;
		case 2:
			return Del_Q_2;
		case 3:
			return Del_Q_3;
		default:
			cout<<"Enter 1,2 or 3 for getting DelnQ\n";
			break;
	}
}
// n refers to  nth gradient and k refers to kth direction of the gradient
const vector<double>& Cell::Get_DelQ_fromCell(const int &n,const int &k) const
{
	switch(n)
	{
		case 2:
				switch(k)
				{
					case 1:
						return Del2_Q_1;
					case 2:
						return Del2_Q_2;
					case 3:
						return Del2_Q_3;
					default:
						cout<<"Trying to fetch Del2Q, should be 1,2,3 please check the number k in passint to this function\n";	
				}
		case 3:
			switch(k)
			{
				case 1:
					return Del3_Q_1;
				case 2:
					return Del3_Q_2;
				case 3:
					return Del3_Q_3;
				default:
					cout<<"Trying to fetch Del2Q, should be 1,2,3 please check the number k in passint to this function\n";
			}
		case 4:
			switch(k)
			{
				case 1:
					return Del4_Q_1;
				case 2:
					return Del4_Q_2;
				case 3:
					return Del4_Q_3;
				default:
					cout<<"Trying to fetch Del2Q, should be 1,2,3 please check the number k in passint to this function\n";
			}
		case 5:
			switch(k)
			{
				case 1:
					return Del5_Q_1;
				case 2:
					return Del5_Q_2;
				case 3:
					return Del5_Q_3;
				default:
					cout<<"Trying to fetch Del2Q, should be 1,2,3 please check the number k in passint to this function\n";
			}
		case 6:
			switch(k)
			{
				case 1:
					return Del6_Q_1;
				case 2:
					return Del6_Q_2;
				case 3:
					return Del6_Q_3;
				default:
					cout<<"Trying to fetch Del2Q, should be 1,2,3 please check the number k in passint to this function\n";
			}
		default:
			cout<<"Trying to fetch the nth gradient so please check the nuber n passing, it should be 2,3,4,5,6\n";
	}
}

void Cell::Cal_Del2_Q()
{
	Del2_Q_1[0]=0.0;Del2_Q_1[1]=0.0;Del2_Q_1[2]=0.0;Del2_Q_1[3]=0.0;Del2_Q_1[4]=0.0;
	Del2_Q_2[0]=0.0;Del2_Q_2[1]=0.0;Del2_Q_2[2]=0.0;Del2_Q_2[3]=0.0;Del2_Q_2[4]=0.0;
	Del2_Q_3[0]=0.0;Del2_Q_3[1]=0.0;Del2_Q_3[2]=0.0;Del2_Q_3[3]=0.0;Del2_Q_3[4]=0.0;

	Front_Face.Cal_Del2Q(this,Front_Cell);
	Back_Face.Cal_Del2Q(this,Back_Cell);
	Left_Face.Cal_Del2Q(this,Left_Cell);
	Right_Face.Cal_Del2Q(this,Right_Cell);
	Top_Face.Cal_Del2Q(this,Top_Cell);
	
	for(int i=0;i<5;i++)
	{
		Del2_Q_1[i] += Front_Face.Get_DelQ_onFace(1)[i];
		Del2_Q_1[i] += Back_Face.Get_DelQ_onFace(1)[i];
		Del2_Q_1[i] += Left_Face.Get_DelQ_onFace(1)[i];
		Del2_Q_1[i] += Right_Face.Get_DelQ_onFace(1)[i];
		Del2_Q_1[i] += Top_Face.Get_DelQ_onFace(1)[i];

		Del2_Q_2[i] += Front_Face.Get_DelQ_onFace(2)[i];
		Del2_Q_2[i] += Back_Face.Get_DelQ_onFace(2)[i];
		Del2_Q_2[i] += Left_Face.Get_DelQ_onFace(2)[i];
		Del2_Q_2[i] += Right_Face.Get_DelQ_onFace(2)[i];
		Del2_Q_2[i] += Top_Face.Get_DelQ_onFace(2)[i];

		Del2_Q_3[i] += Front_Face.Get_DelQ_onFace(3)[i];
		Del2_Q_3[i] += Back_Face.Get_DelQ_onFace(3)[i];
		Del2_Q_3[i] += Left_Face.Get_DelQ_onFace(3)[i];
		Del2_Q_3[i] += Right_Face.Get_DelQ_onFace(3)[i];
		Del2_Q_3[i] += Top_Face.Get_DelQ_onFace(3)[i];
	}
	if(No_of_Faces==6)
	{
		for(int i=0;i<5;i++)
		{
			Del2_Q_1[i] += Bottom_Face.Get_DelQ_onFace(1)[i];
			Del2_Q_2[i] += Bottom_Face.Get_DelQ_onFace(2)[i];
			Del2_Q_3[i] += Bottom_Face.Get_DelQ_onFace(3)[i];
		}
	}
	for(int i=0;i<5;i++)
	{
		Del2_Q_1[i]  *= inv_vol;
		Del2_Q_2[i]  *= inv_vol;
		Del2_Q_3[i]  *= inv_vol;
	}
/*cout<<self_index<<"\tDel2Q_1\t"<<Del2_Q_1[0]<<"\t"<<Del2_Q_1[1]<<"\t"<<Del2_Q_1[2]<<"\t"<<Del2_Q_1[3]<<"\t"<<Del2_Q_1[4]<<"\n";
cout<<self_index<<"\tDel2Q_2\t"<<Del2_Q_2[0]<<"\t"<<Del2_Q_2[1]<<"\t"<<Del2_Q_2[2]<<"\t"<<Del2_Q_2[3]<<"\t"<<Del2_Q_2[4]<<"\n";
cout<<self_index<<"\tDel2Q_3\t"<<Del2_Q_3[0]<<"\t"<<Del2_Q_3[1]<<"\t"<<Del2_Q_3[2]<<"\t"<<Del2_Q_3[3]<<"\t"<<Del2_Q_3[4]<<"\n";*/
}

void Cell:: Cal_Del3_Q()
{
	Del3_Q_1[0]=0.0;Del3_Q_1[1]=0.0;Del3_Q_1[2]=0.0;Del3_Q_1[3]=0.0;Del3_Q_1[4]=0.0;
	Del3_Q_2[0]=0.0;Del3_Q_2[1]=0.0;Del3_Q_2[2]=0.0;Del3_Q_2[3]=0.0;Del3_Q_2[4]=0.0;
	Del3_Q_3[0]=0.0;Del3_Q_3[1]=0.0;Del3_Q_3[2]=0.0;Del3_Q_3[3]=0.0;Del3_Q_3[4]=0.0;

	for(int i=0;i<5;i++)
	{
		Del3_Q_1[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,2,1)[i];
		Del3_Q_1[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,2,1)[i];
		Del3_Q_1[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,2,1)[i];
		Del3_Q_1[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,2,1)[i];
		Del3_Q_1[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,2,1)[i];
		
		Del3_Q_2[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,2,2)[i];
		Del3_Q_2[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,2,2)[i];
		Del3_Q_2[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,2,2)[i];
		Del3_Q_2[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,2,2)[i];
		Del3_Q_2[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,2,2)[i];

		Del3_Q_3[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,2,3)[i];
		Del3_Q_3[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,2,3)[i];
		Del3_Q_3[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,2,3)[i];
		Del3_Q_3[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,2,3)[i];
		Del3_Q_3[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,2,3)[i];
	}
	if(No_of_Faces ==6)
	{
		for(int i=0;i<5;i++)
		{
			Del3_Q_1[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,2,1)[i];
			Del3_Q_2[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,2,2)[i];
			Del3_Q_3[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,2,3)[i];
		}
	}
	for(int i=0;i<5;i++)
	{
		Del3_Q_1[i]  *= inv_vol;
		Del3_Q_2[i]  *= inv_vol;
		Del3_Q_3[i]  *= inv_vol;
	}
/*cout<<self_index<<"\tDel3Q_1\t"<<Del3_Q_1[0]<<"\t"<<Del3_Q_1[1]<<"\t"<<Del3_Q_1[2]<<"\t"<<Del3_Q_1[3]<<"\t"<<Del3_Q_1[4]<<"\n";
cout<<self_index<<"\tDel3Q_2\t"<<Del3_Q_2[0]<<"\t"<<Del3_Q_2[1]<<"\t"<<Del3_Q_2[2]<<"\t"<<Del3_Q_2[3]<<"\t"<<Del3_Q_2[4]<<"\n";
cout<<self_index<<"\tDel3Q_3\t"<<Del3_Q_3[0]<<"\t"<<Del3_Q_3[1]<<"\t"<<Del3_Q_3[2]<<"\t"<<Del3_Q_3[3]<<"\t"<<Del3_Q_3[4]<<"\n";*/
		
}

void Cell:: Cal_Del4_Q()
{
	Del4_Q_1[0]=0.0;Del4_Q_1[1]=0.0;Del4_Q_1[2]=0.0;Del4_Q_1[3]=0.0;Del4_Q_1[4]=0.0;
	Del4_Q_2[0]=0.0;Del4_Q_2[1]=0.0;Del4_Q_2[2]=0.0;Del4_Q_2[3]=0.0;Del4_Q_2[4]=0.0;
	Del4_Q_3[0]=0.0;Del4_Q_3[1]=0.0;Del4_Q_3[2]=0.0;Del4_Q_3[3]=0.0;Del4_Q_3[4]=0.0;

	for(int i=0;i<5;i++)
	{
		Del4_Q_1[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,3,1)[i];
		Del4_Q_1[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,3,1)[i];
		Del4_Q_1[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,3,1)[i];
		Del4_Q_1[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,3,1)[i];
		Del4_Q_1[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,3,1)[i];
		
		Del4_Q_2[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,3,2)[i];
		Del4_Q_2[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,3,2)[i];
		Del4_Q_2[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,3,2)[i];
		Del4_Q_2[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,3,2)[i];
		Del4_Q_2[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,3,2)[i];

		Del4_Q_3[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,3,3)[i];
		Del4_Q_3[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,3,3)[i];
		Del4_Q_3[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,3,3)[i];
		Del4_Q_3[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,3,3)[i];
		Del4_Q_3[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,3,3)[i];
	}
	if(No_of_Faces ==6)
	{
		for(int i=0;i<5;i++)
		{
			Del4_Q_1[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,3,1)[i];
			Del4_Q_2[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,3,2)[i];
			Del4_Q_3[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,3,3)[i];
		}
	}
	for(int i=0;i<5;i++)
	{
		Del4_Q_1[i]  *= inv_vol;
		Del4_Q_2[i]  *= inv_vol;
		Del4_Q_3[i]  *= inv_vol;
	}
/*cout<<self_index<<"\tDel4Q_1\t"<<Del4_Q_1[0]<<"\t"<<Del4_Q_1[1]<<"\t"<<Del4_Q_1[2]<<"\t"<<Del4_Q_1[3]<<"\t"<<Del4_Q_1[4]<<"\n";
cout<<self_index<<"\tDel4Q_2\t"<<Del4_Q_2[0]<<"\t"<<Del4_Q_2[1]<<"\t"<<Del4_Q_2[2]<<"\t"<<Del4_Q_2[3]<<"\t"<<Del4_Q_2[4]<<"\n";
cout<<self_index<<"\tDel4Q_3\t"<<Del4_Q_3[0]<<"\t"<<Del4_Q_3[1]<<"\t"<<Del4_Q_3[2]<<"\t"<<Del4_Q_3[3]<<"\t"<<Del4_Q_3[4]<<"\n";	*/
}


void Cell:: Cal_Del5_Q()
{
	Del5_Q_1[0]=0.0;Del5_Q_1[1]=0.0;Del5_Q_1[2]=0.0;Del5_Q_1[3]=0.0;Del5_Q_1[4]=0.0;
	Del5_Q_2[0]=0.0;Del5_Q_2[1]=0.0;Del5_Q_2[2]=0.0;Del5_Q_2[3]=0.0;Del5_Q_2[4]=0.0;
	Del5_Q_3[0]=0.0;Del5_Q_3[1]=0.0;Del5_Q_3[2]=0.0;Del5_Q_3[3]=0.0;Del5_Q_3[4]=0.0;

	for(int i=0;i<5;i++)
	{
		Del5_Q_1[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,4,1)[i];
		Del5_Q_1[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,4,1)[i];
		Del5_Q_1[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,4,1)[i];
		Del5_Q_1[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,4,1)[i];
		Del5_Q_1[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,4,1)[i];
		
		Del5_Q_2[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,4,2)[i];
		Del5_Q_2[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,4,2)[i];
		Del5_Q_2[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,4,2)[i];
		Del5_Q_2[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,4,2)[i];
		Del5_Q_2[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,4,2)[i];

		Del5_Q_3[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,4,3)[i];
		Del5_Q_3[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,4,3)[i];
		Del5_Q_3[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,4,3)[i];
		Del5_Q_3[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,4,3)[i];
		Del5_Q_3[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,4,3)[i];
	}
	if(No_of_Faces ==6)
	{
		for(int i=0;i<5;i++)
		{
			Del5_Q_1[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,4,1)[i];
			Del5_Q_2[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,4,2)[i];
			Del5_Q_3[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,4,3)[i];
		}
	}
	for(int i=0;i<5;i++)
	{
		Del5_Q_1[i]  *= inv_vol;
		Del5_Q_2[i]  *= inv_vol;
		Del5_Q_3[i]  *= inv_vol;
	}
}

void Cell:: Cal_Del6_Q()
{
	Del6_Q_1[0]=0.0;Del6_Q_1[1]=0.0;Del6_Q_1[2]=0.0;Del6_Q_1[3]=0.0;Del6_Q_1[4]=0.0;
	Del6_Q_2[0]=0.0;Del6_Q_2[1]=0.0;Del6_Q_2[2]=0.0;Del6_Q_2[3]=0.0;Del6_Q_2[4]=0.0;
	Del6_Q_3[0]=0.0;Del6_Q_3[1]=0.0;Del6_Q_3[2]=0.0;Del6_Q_3[3]=0.0;Del6_Q_3[4]=0.0;

	for(int i=0;i<5;i++)
	{
		Del6_Q_1[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,5,1)[i];
		Del6_Q_1[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,5,1)[i];
		Del6_Q_1[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,5,1)[i];
		Del6_Q_1[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,5,1)[i];
		Del6_Q_1[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,5,1)[i];
		
		Del6_Q_2[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,5,2)[i];
		Del6_Q_2[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,5,2)[i];
		Del6_Q_2[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,5,2)[i];
		Del6_Q_2[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,5,2)[i];
		Del6_Q_2[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,5,2)[i];

		Del6_Q_3[i] += Front_Face.Calculate_Grad_n_Q(this,Front_Cell,5,3)[i];
		Del6_Q_3[i] += Back_Face.Calculate_Grad_n_Q(this,Back_Cell,5,3)[i];
		Del6_Q_3[i] += Left_Face.Calculate_Grad_n_Q(this,Left_Cell,5,3)[i];
		Del6_Q_3[i] += Right_Face.Calculate_Grad_n_Q(this,Right_Cell,5,3)[i];
		Del6_Q_3[i] += Top_Face.Calculate_Grad_n_Q(this,Top_Cell,5,3)[i];
	}
	if(No_of_Faces ==6)
	{
		for(int i=0;i<5;i++)
		{
			Del6_Q_1[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,5,1)[i];
			Del6_Q_2[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,5,2)[i];
			Del6_Q_3[i] += Bottom_Face.Calculate_Grad_n_Q(this,Bottom_Cell,5,3)[i];
		}
	}
	for(int i=0;i<5;i++)
	{
		Del6_Q_1[i]  *= inv_vol;
		Del6_Q_2[i]  *= inv_vol;
		Del6_Q_3[i]  *= inv_vol;
	}
}
