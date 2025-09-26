#include "../Basic_Function_Files/Geometry_Header.h"


void Generate_Arc(double & radius,int & N_Points,vector<Point> & Arc,int & sign,Point & P_Ref)
{
	Point P;
	if(sign)
	{
		for(double theta=180;theta>=0;theta-=(180.0/(N_Points-1)))
		{
			P(radius,theta,1);
			Arc.push_back(P + P_Ref);
		}
	}
	else
	{
		for(double theta=0.0;theta<=180.0;theta+=(180.0/(N_Points-1)))
		{
			P(radius,theta,1);
//			P+=P_Ref;
			P.Print();
			Arc.push_back(P);
		}
	}
}

void Generate_Arc(int & N_Points,vector<Point> & Arc,Point & Initial_Point, Point & Final_Point)
{
	Point P;
	double dx=0.0,dy=0.0,y=0.0,a,b,c;
	dx = (Final_Point.Get_x() - Initial_Point.Get_x())/(N_Points-1);
	dy = (Final_Point.Get_y() - Initial_Point.Get_y())/(N_Points-1);
	cout<<N_Points<<"\t"<<dx<<"\t"<<dy<<endl;
//	Arc.push_back(Initial_Point);
	for(int i=0;i<N_Points;i++)
	{
		
//		if(i*dx<=0.1)
//			y =0.29;
//		else
		a = 0.29;
		b = -2.0*a + tan(25*radian);
		c = a - tan(25*radian);
		y = a + b*i*dx + c*pow(i*dx,2);
		P(i*dx,y,0.0);
//		P.Print();
		Arc.push_back(P); 
	}
//		Arc.push_back(Final_Point);
}

void Generate_Arc(Point & P1, Point & P2,int & N_Points,vector<Point> & Arc, double & Percent_Chord)
{
	Point P;
	double Radius =0.0,dx=0.0,dy=0.0,Chord_Length=0.0,Y_Max=0.0,Tan_Theta;
	
	P2.Print();
	P1.Print();
	dx = P2.Get_x() - P1.Get_x();
	dy = P2.Get_y() - P1.Get_y();
	cout<<dx<<"\t"<<dy<<endl;
	Chord_Length = sqrt(dx*dx + dy*dy);
	cout<<Chord_Length<<endl;
	Y_Max = (Percent_Chord/100.0)*Chord_Length;
	Tan_Theta = Y_Max/(0.5*Chord_Length);
	
//	cout<<Y_Max<<"\t"<<Tan_Theta<<endl;
	
	dx = Chord_Length/(N_Points-1);
	
	for(int i=0;i<N_Points;i++)
	{
		dy=Percent_Chord*sin((M_PI)*(i*dx/Chord_Length));   	
		P(P1.Get_x()+i*dx,P1.Get_y()+dy,0.0);
//		cout<<i<<"\t";
//		P.Print();
		Arc.push_back(P);
	}
	
	
	
}

void Append_List(vector<Point> & List1, vector<Point> & List2)
{
	cout<<"Appending Lists \n";
	List1[List1.size()-1].Print();
	List2[0].Print();
	if(List1[List1.size()-1]==List2[0])
	{
		for(unsigned int i=1;i<List2.size();i++)
		{
			List1.push_back(List2[i]);
		}
	}
	else
	{
		for(unsigned int i=0;i<List2.size();i++)
		{
			List1.push_back(List2[i]);
		}
	}
}

void Reverse_List(vector<Point> & List)
{
	cout<<"Reversing the list of points"<<endl;
	int size = List.size();
	vector<Point> Temp_List;
	for(int i = size-1;i>=0;i--)
	{
		Temp_List.push_back(List[i]);
	}
	List = Temp_List;
}

int main()
{
	Grid Grid_Flow_Over_Cylinder;
	string Solver_Input_File,Grid_View_File;


	Solver_Input_File="../Grid_Files/Flow_Over_Cylinder_Files/Double_Bump_301_301.txt";
	Grid_View_File="../Grid_Files/Flow_Over_Cylinder_Files/Double_301_301.vtk";

//	Solver_Input_File="../Grid_Files/Expansion_Ramp_Files/Expansion_Ramp_1001_1001.txt";
//	Grid_View_File="../Grid_Files/Expansion_Ramp_Files/Expansion_Ramp_1001_1001.vtk";

	
	int Nop_Circle,Nop_X,Nop_X1,Nop_Y,sign=0,Total_Points_X,Test_Case;
	double Radius1,Chord_Length,Percent_Chord=0.285,Percent_Chord1=-0.285;
	Line L1,L2,L3,L4,L5,L6;
	Point P1,P2,P3,P4,P5,P6,P7,P8;
	vector<Point> Line_List1,Line_List2,Line_List3,Line_List4,Grid_Plane_List,Arc_Point_List1,Arc_Point_List2,Temp_List,Arc_Point_List;
	bool Is_Periodic_Boundary = false;

	Chord_Length = 1.0;
	Radius1 = 0.06;

	Nop_X = 101;
	Nop_X1 = 101;
	Nop_Circle = 101;
	Nop_Y = 101;
//	Total_Points_X = 2*Nop_X + Nop_Circle - 2;
	Total_Points_X = Nop_X1 + Nop_X + Nop_Circle - 2;
	Test_Case = 2;
	switch(Test_Case)
	{
		case 1: // for flow over a circular bump
		P1(0.0*Chord_Length,0.0,0.0);
		P2 (0.2*Chord_Length,0.0,0.0);
		P3(0.6*Chord_Length,0.0,0.0);
		P4(1.0*Chord_Length,0.0,0.0);
		P5(1.0*Chord_Length,1.0*Chord_Length,0.0);
		P6(0.0*Chord_Length,1.0*Chord_Length,0.0);
	
		L1.generate(P1,P2,Nop_X);	
	//	Generate_Arc(,Nop_Circle,Arc_Point_List,sign,P2);
//		Generate_Arc(P2,P3,Nop_Circle,Arc_Point_List,Percent_Chord);
		Generate_Arc(P2,P3,Nop_Circle,Arc_Point_List,Percent_Chord);		
		L2.generate(P3,P4,Nop_X1);
	
		L3.generate(P4,P5,Nop_Y);
		L4.generate(P5,P6,Total_Points_X);
		L5.generate(P6,P1,Nop_Y);
	
		Line_List1 = L1.Get_Point_list();
		Append_List(Line_List1,Arc_Point_List);
		Temp_List = L2.Get_Point_list();
		Append_List(Line_List1,Temp_List);
		Line_List2 = L3.Get_Point_list();
		Line_List3 = L4.Get_Point_list();
		Line_List4 = L5.Get_Point_list();
		cout<<Line_List1.size()<<"\t"<<Line_List2.size()<<"\t"<<Line_List3.size()<<"\t"<<Line_List4.size()<<endl;
		Grid_Flow_Over_Cylinder(Line_List1,Line_List2,Line_List3,Line_List4);
		Grid_Flow_Over_Cylinder.Generate_Grid(Is_Periodic_Boundary);
		Grid_Plane_List = Grid_Flow_Over_Cylinder.Get_Grid_List();
	
	
	 	Identify_Cells(Total_Points_X,Nop_Y);
	 	Identify_Neighbours(Total_Points_X,Nop_Y);
	 	write_VTK(Total_Points_X,Nop_Y,Grid_Plane_List,Grid_View_File);
	 	write_inputfile(Total_Points_X,Nop_Y,Grid_Plane_List,Solver_Input_File);	
		
		case 2: // for flow over a Double circular bump
		
		P1(0.0*Chord_Length,0.0,0.0);
		P2(0.5*Chord_Length,0.0,0.0);
		P3(1.5*Chord_Length,0.0,0.0);
		P4(2.5*Chord_Length,0.0,0.0);
		
		P5(0.0,1.0*Chord_Length,0.0);
		P6(1.0*Chord_Length,1.0*Chord_Length,0.0);
		P7(2.0*Chord_Length,1.0*Chord_Length,0.0);
		P8(2.5*Chord_Length,1.0*Chord_Length,0.0);
		
	
		L1.generate(P1,P2,Nop_X);	
		Generate_Arc(P2,P3,Nop_Circle,Arc_Point_List1,Percent_Chord);
		L2.generate(P3,P4,Nop_X1);
	
		
		L3.generate(P4,P8,Nop_Y);
		
		L4.generate(P5,P6,Nop_X);	
		Generate_Arc(P6,P7,Nop_Circle,Arc_Point_List2,Percent_Chord1);
		L5.generate(P7,P8,Nop_X1);
		

		L6.generate(P5,P1,Nop_Y);
	
		Line_List1 = L1.Get_Point_list();
		Append_List(Line_List1,Arc_Point_List1);
		Temp_List = L2.Get_Point_list();
		Append_List(Line_List1,Temp_List);
		
		Line_List2 = L3.Get_Point_list();
		
		Line_List3 = L4.Get_Point_list();
		Append_List(Line_List3,Arc_Point_List2);
		Temp_List = L5.Get_Point_list();
		Append_List(Line_List3,Temp_List);
		
		Reverse_List(Line_List3);		

		
		Line_List4 = L6.Get_Point_list();
		cout<<Line_List1.size()<<"\t"<<Line_List2.size()<<"\t"<<Line_List3.size()<<"\t"<<Line_List4.size()<<endl;
		Grid_Flow_Over_Cylinder(Line_List1,Line_List2,Line_List3,Line_List4);
		Grid_Flow_Over_Cylinder.Generate_Grid(Is_Periodic_Boundary);
		Grid_Plane_List = Grid_Flow_Over_Cylinder.Get_Grid_List();
	
	
	 	Identify_Cells(Total_Points_X,Nop_Y);
	 	Identify_Neighbours(Total_Points_X,Nop_Y);
	 	write_VTK(Total_Points_X,Nop_Y,Grid_Plane_List,Grid_View_File);
	 	write_inputfile(Total_Points_X,Nop_Y,Grid_Plane_List,Solver_Input_File);	
		
		break;
		case 3: // For flow over a expansion Ramp
		P1(0.0,0.29,0.0);	
		P2 (1.0,0.0,0.0);
		P3(1.0,1.0,0.0);
		P4(0.0,1.0,0.0);
		Nop_X = 1001;
		Nop_Y = 1001;				
	//	Generate_Arc(,Nop_Circle,Arc_Point_List,sign,P2);
//		Generate_Arc(P1,P2,Nop_X,Line_List1,Percent_Chord);
		Generate_Arc(Nop_X,Line_List1,P1,P2);
		P2 = Line_List1[Nop_X-1];
		P2.Print();
		L2.generate(P2,P3,Nop_Y);
		L3.generate(P3,P4,Nop_X);
		L4.generate(P4,P1,Nop_Y);
		cout<<"Size of Arc points list\t"<<Line_List1.size()<<endl;
//		for(unsigned int i=0;i<Line_List1.size();i++)
//			Line_List1[i].Print();
		cout<<"Last point on the arc\t";
			Line_List1[Nop_X-1].Print(); cout<<endl;

		Line_List2 = L2.Get_Point_list();
				Line_List2[0].Print();
		Line_List3 = L3.Get_Point_list();
		Line_List4 = L4.Get_Point_list();
		cout<<Line_List1.size()<<"\t"<<Line_List2.size()<<"\t"<<Line_List3.size()<<"\t"<<Line_List4.size()<<endl;
		Grid_Flow_Over_Cylinder(Line_List1,Line_List2,Line_List3,Line_List4);
		Grid_Flow_Over_Cylinder.Generate_Grid(Is_Periodic_Boundary);
		Grid_Plane_List = Grid_Flow_Over_Cylinder.Get_Grid_List();
	
	
	 	Identify_Cells(Nop_X,Nop_Y);
	 	Identify_Neighbours(Nop_X,Nop_Y);
	 	write_VTK(Nop_X,Nop_Y,Grid_Plane_List,Grid_View_File);
	 	write_inputfile(Nop_X,Nop_Y,Grid_Plane_List,Solver_Input_File);	
		
		break;
	}
	
		
return 0;
}

