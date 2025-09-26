#include "../Basic_Function_Files/Geometry_Header.h"


void Generate_Arc(double & radius,int & N_Points,vector<Point> & Arc,int & sign)
{
	Point P;
	if(sign)
	{
		for(double theta=180;theta>=0;theta-=(180.0/(N_Points-1)))
		{
			P(radius,theta,1);
			Arc.push_back(P);
		}
	}
	else
	{
		for(double theta=0.0;theta<=180.0;theta+=(180.0/(N_Points-1)))
		{
			P(radius,theta,1);
			Arc.push_back(P);
		}
	}
	
}

int main()
{
	Grid Grid_Half_Cylinder;
	string Solver_Input_File,Grid_View_File;
	
	int sign =1;
	
	Solver_Input_File="../Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_61_21.txt";
	Grid_View_File="../Grid_Files/Half_Cylinder_Files/Flow_Over_Cylinder_61_21.vtk";
	
	double radius1 =1.0,radius2=8;
	int noP_x,noP_y;
	bool With_Periodic_Boundary =false,EnableElliptic=false;
	int iterations = 10000;
	
	noP_x=61; 
	noP_y=21;
	
	Line Line1,Line2,Line3,Line4;
	
	vector<Point> Line_List1,Line_List2,Line_List3,Line_List4,Stacked_Grid_List,Grid_Plane_List;

	Point P1,P2,P3,P4,P;
	P1(-radius2,0.0,0.0);
	P2(-radius1,0.0,0.0);
	P3(radius1,0.0,0.0);
	P4(radius2,0.0,0.0);

	Line1.generate(P1,P2,noP_y);
	Line3.generate(P3,P4,noP_y);


	Line_List4 = Line1.Get_Point_list();
	Line_List2 = Line3.Get_Point_list();
	Generate_Arc(radius1,noP_x,Line_List1,sign);
	sign = 0;
	Generate_Arc(radius2,noP_x,Line_List3,sign);

	
	cout<<Line_List1.size()<<"\t"<<Line_List2.size()<<"\t"<<Line_List3.size()<<"\t"<<Line_List4.size()<<endl;
	Grid_Half_Cylinder(Line_List1,Line_List2,Line_List3,Line_List4);
//	Grid_Half_Cylinder.Generate_Grid(Is_Periodic_Boundary);
	Grid_Half_Cylinder.Generate_Grid(With_Periodic_Boundary,EnableElliptic,iterations);
	Grid_Plane_List = Grid_Half_Cylinder.Get_Grid_List();
	
	Identify_Cells(noP_x,noP_y);
	Identify_Neighbours(noP_x,noP_y);
//	Identify_Cells_and_Assign_Neighbours(noP_x,noP_y);	
	write_VTK(noP_x,noP_y,Grid_Plane_List,Grid_View_File);
	write_inputfile(noP_x,noP_y,Grid_Plane_List,Solver_Input_File);
return 0;
}

