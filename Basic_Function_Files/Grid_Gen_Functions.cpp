#include "Geometry_Header.h"

void Grid::operator()(vector<Point> & List1, vector<Point> & List2,vector<Point> & List3, vector<Point> & List4 )
{
	Point P1,P2,P3,P4,P5,P6,P7,P8;
	int No_P_Curve1, No_P_Curve2,No_P_Curve3,No_P_Curve4,i,j,Total_Points;
	No_P_Curve1 = List1.size();
	No_P_Curve2 = List2.size();
	No_P_Curve3 = List3.size();
	No_P_Curve4 = List4.size();
	
	if((No_P_Curve1 == No_P_Curve3)&&(No_P_Curve2==No_P_Curve4))
	{	
 		List4[No_P_Curve4-1].Print();
 		List1[0].Print();
 		List1[No_P_Curve1-1].Print();
		List2[0].Print();
		List2[No_P_Curve2-1].Print();
		List3[0].Print();
 		List3[No_P_Curve3-1].Print();
 		List4[0].Print();
		cout<<"Checking out for Closed Curve of Points List"<<endl;
		if(		(List4[No_P_Curve4-1]==List1[0])
			&&	(List1[No_P_Curve1-1]==List2[0])
			&&	(List3[No_P_Curve3-1]==List4[0])
			&&	(List2[No_P_Curve2-1]==List3[0])
		)
		{
			cout<<"The Curves form a closed loop.....\n Proceeding for Grid Generation"<<endl;
			nx = No_P_Curve1;
			ny = No_P_Curve2;
			Total_Points = nx*ny;
		// 	cout<<"=====================\n";
			x.resize(nx,vector<double>(ny,0.0));
			y.resize (nx,vector<double>(ny,0.0));
			xtemp.resize (nx,vector<double>(ny,0.0));
			ytemp.resize (nx,vector<double>(ny,0.0));
			erx.resize (nx,vector<double>(ny,0.0));
			ery.resize (nx,vector<double>(ny,0.0));

			for ( i=0;i<nx;i++)
			{
				j=0;
//  		 		cout<<i<<"\t"<<j<<"\t"<<i+j<<endl;
				x[i][j] = List1[i].Get_x();
				y[i][j] = List1[i].Get_y();
// 		 		cout<<x[i][j]<<"\t"<<y[i][j]<<endl;
			}
// 		 cout<<"=====================\n";
			for (j=0;j<ny;j++)
			{
				i=nx-1;
// 		 		cout<<i<<"\t"<<j<<"\t"<<i+j<<endl;
				x[i][j]=List2[j].Get_x();
				y[i][j]=List2[j].Get_y();
// 		 		cout<<x[i][j]<<"\t"<<y[i][j]<<endl;
			}
// 		 cout<<"=====================\n";
			for(i=0;i<List3.size();i++)
			{
				j=ny-1;
// 		 		cout<<i<<"\t"<<j<<"\t"<<( i+j)<<endl;
				x[nx-1-i][j]=List3[i].Get_x();
				y[nx-1-i][j]=List3[i].Get_y();
// 		 		cout<<x[nx-1-i][j]<<"\t"<<y[nx-1-i][j]<<endl;
			}
// 		cout<<"=====================\n";
			for(j=0;j<List4.size();j++)
			{
				i=0;
// 				cout<<i<<"\t"<<j<<"\t"<<(i+j)<<endl;
				x[i][ny-1-j] = List4[j].Get_x();
				y[i][ny-1-j] = List4[j].Get_y();
// 		 		cout<<x[i][ny-1-j]<<"\t"<<y[i][ny-1-j]<<endl;
			}
// 		cout<<"=====================\n";
			
		}
		else
		{
			cout<<"The Points List doesn't Form a Closed Loop, Please Check them"<<endl;
			exit(0);
		}
	}
	else
	{
		cout<<"Miss-Match of number of Points on the curves.........\n Can not proceed to generate a grid"<<endl;
		exit(0);
	}
	
}