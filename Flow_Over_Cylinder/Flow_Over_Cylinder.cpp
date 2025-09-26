#include "../Basic_Function_Files/Geometry_Header.h"


void Generate_Arc(double & radius,int & N_Points,vector<Point> & Arc,int & sign,Point & P_Ref)
{
	cout<<"Generating Semi Circle with start point as P_Ref"<<endl;
	P_Ref.Print();
	Point P;
	if(sign)
	{
		for(double theta=180;theta>=0;theta-=(180.0/(N_Points-1)))
		{
			P(radius,theta,1);
			P+=P_Ref;
			P.Print();
			Arc.push_back(P);
		}
	}
	else
	{
		for(double theta=0.0;theta<=180.0;theta+=(180.0/(N_Points-1)))
		{
			P(radius,theta,1);
			P.Print();
			if(theta<=90)
				P-=P_Ref;
			else 
				P+=P_Ref;
			Arc.push_back(P);
		}
	}
}


void Generate_Arc(int &N_Points, vector<Point> &Arc, Point &Initial_Point, Point &Final_Point)
{
    Point P, Mid_P;
    double dx, dy, r;

    // Calculate the midpoint between Initial_Point and Final_Point
    Mid_P = (Initial_Point + Final_Point) * 0.5;

    // Calculate dx, dy, and radius
    dx = Final_Point.Get_x() - Initial_Point.Get_x();
    dy = Final_Point.Get_y() - Initial_Point.Get_y();
    r = 0.5 * sqrt(dx * dx + dy * dy);

    // Calculate the center of the circle (perpendicular bisector pointing downward)
    double cx = Mid_P.Get_x();
    double cy = Mid_P.Get_y() - sqrt(r * r - (dx * dx) / 4);

    Point Center(cx, cy, 0.0);

    // Calculate the angle step
    double start_angle = atan2(Initial_Point.Get_y() - cy, Initial_Point.Get_x() - cx);
    double end_angle = atan2(Final_Point.Get_y() - cy, Final_Point.Get_x() - cx);

    // Ensure that the arc is generated in the correct direction
    if (start_angle < end_angle) {
        end_angle -= 2 * M_PI;
    }

    double angle_step = (end_angle - start_angle) / (N_Points - 1);

    // Include the Initial_Point as the first point in the arc
    Arc.push_back(Initial_Point);

    // Generate points along the arc
    for (int i = 1; i < N_Points - 1; ++i)
    {
        double angle = start_angle + i * angle_step;
        double x = cx + r * cos(angle);
        double y = cy + r * sin(angle);
        P(x, y, 0.0);
        Arc.push_back(P);
    }

    // Include the Final_Point as the last point in the arc
    Arc.push_back(Final_Point);

    // Print generated arc points
    for (const auto &point : Arc)
    {
        point.Print();
    }
}
/*void Generate_Arc(int &N_Points, vector<Point> &Arc, Point &Initial_Point, Point &Final_Point)
{
    Point P, Mid_P;
    double dx, dy, r;

    // Calculate the midpoint between Initial_Point and Final_Point
    Mid_P = (Initial_Point + Final_Point) * 0.5;

    // Calculate dx, dy, and radius
    dx = Final_Point.Get_x() - Initial_Point.Get_x();
    dy = Final_Point.Get_y() - Initial_Point.Get_y();
    r = 0.5 * sqrt(dx * dx + dy * dy);

    // Calculate the center of the circle (perpendicular bisector pointing upward)
    double cx = Mid_P.Get_x();
    double cy = Mid_P.Get_y() + sqrt(r * r - (dx * dx) / 4);

    Point Center(cx, cy, 0.0);

    // Calculate the angle step
    double start_angle = atan2(Initial_Point.Get_y() - cy, Initial_Point.Get_x() - cx);
    double end_angle = atan2(Final_Point.Get_y() - cy, Final_Point.Get_x() - cx);

    // Ensure that the arc is generated in the correct direction
    if (start_angle > end_angle) {
        end_angle += 2 * M_PI;
    }

    double angle_step = (end_angle - start_angle) / (N_Points - 1);

    // Include the Initial_Point as the first point in the arc
    Arc.push_back(Initial_Point);

    // Generate points along the arc
    for (int i = 1; i < N_Points - 1; ++i)
    {
        double angle = start_angle + i * angle_step;
        double x = cx + r * cos(angle);
        double y = cy + r * sin(angle);
        P(x, y, 0.0);
        Arc.push_back(P);
    }

    // Include the Final_Point as the last point in the arc
    Arc.push_back(Final_Point);

    // Print generated arc points
    for (const auto &point : Arc)
    {
        point.Print();
    }
}

void Generate_Arc(int &N_Points, vector<Point> &Arc, Point &Initial_Point, Point &Final_Point)
{
    Point P;
    double dx, dy, r;

    // Calculate the chord length and radius
    dx = Final_Point.Get_x() - Initial_Point.Get_x();
    dy = Final_Point.Get_y() - Initial_Point.Get_y();
    r = 0.5 * sqrt(dx * dx + dy * dy);

    // Calculate the midpoint
    Point Mid_P = (Initial_Point + Final_Point) * 0.5;

    // Calculate the center of the circle
    double cx = Mid_P.Get_x();
    double cy = Mid_P.Get_y() + sqrt(r * r - (dx * dx) / 4);
    Point Center(cx, cy, 0.0);

    // Generate points using cosine clustering
    for (int i = 0; i < N_Points; ++i)
    {
        double theta = M_PI * (1 - cos(M_PI * i / (N_Points - 1))) / 2;
        double x = cx + r * cos(theta);
        double y = cy + r * sin(theta);
        P(x, y, 0.0);
        Arc.push_back(P);
    }

    // Ensure the arc includes the initial and final points
    Arc[0] = Initial_Point;
    Arc[N_Points - 1] = Final_Point;
}*/

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


	Solver_Input_File="../Grid_Files/Flow_Over_Bump/Single_Bump_181_61.txt";
	Grid_View_File="../Grid_Files/Flow_Over_Bump/Single_Bump_181_61.vtk";

//	Solver_Input_File="../Grid_Files/Expansion_Ramp_Files/Expansion_Ramp_1001_1001.txt";
//	Grid_View_File="../Grid_Files/Expansion_Ramp_Files/Expansion_Ramp_1001_1001.vtk";

	
	int Nop_Circle,Nop_X,Nop_X1,Nop_Y,sign=0,Total_Points_X,Test_Case,iterations = 50000;
	double Radius1,Chord_Length,Percent_Chord=0.5,Percent_Chord1=-1.0;
	Line L1,L2,L3,L4,L5,L6;
	Point P1,P2,P3,P4,P5,P6,P7,P8;
	vector<Point> Line_List1,Line_List2,Line_List3,Line_List4,Grid_Plane_List,Arc_Point_List1,Arc_Point_List2,Temp_List,Arc_Point_List;
	bool Is_Periodic_Boundary = false, Enable_Elliptic_Solver = true;

	Chord_Length = 1.0;
	Radius1 = 1.0;

	Nop_X = 61;
	Nop_X1 = 61;
	Nop_Circle = 61;
	Nop_Y = 61;
//	Total_Points_X = 2*Nop_X + Nop_Circle - 2;
	Total_Points_X = Nop_X1 + Nop_X + Nop_Circle - 2;
	Test_Case = 1;
	switch(Test_Case)
	{
		case 1: // for flow over a circular bump
		P1(-2.0*Chord_Length,0.0,0.0);
		P2 (1.5*Chord_Length,0.0,0.0);
		P3(2.5*Chord_Length,0.0,0.0);
		P4(10.0*Chord_Length,0.0,0.0);
		P5(10.0*Chord_Length,5.0*Chord_Length,0.0);
		P6(-2.0*Chord_Length,5.0*Chord_Length,0.0);
	
		L1.generate(P1,P2,Nop_X);	
//		Generate_Arc(Radius1,Nop_Circle,Arc_Point_List,sign,P2);
//		Generate_Arc(P2,P3,Nop_Circle,Arc_Point_List,Percent_Chord);	
		Generate_Arc(Nop_Circle,Arc_Point_List,P2,P3);	
		
		//P3 = Arc_Point_List[Nop_Circle-1];	
		L2.generate(P3,P4,Nop_X1);
	
		L3.generate(P4,P5,Nop_Y);
		L4.generate(P5,P6,Total_Points_X);
		L5.generate(P6,P1,Nop_Y);
	
		Line_List1 = L1.Get_Point_list();
		Append_List(Line_List1,Arc_Point_List);
		Temp_List = L2.Get_Point_list();
		Append_List(Line_List1,Temp_List);
		for(unsigned int i=0;i<Line_List1.size();i++)
		{
			Line_List1[i].Print();
		}
		Line_List2 = L3.Get_Point_list();
		Line_List3 = L4.Get_Point_list();
		Line_List4 = L5.Get_Point_list();
		cout<<Line_List1.size()<<"\t"<<Line_List2.size()<<"\t"<<Line_List3.size()<<"\t"<<Line_List4.size()<<endl;
		Grid_Flow_Over_Cylinder(Line_List1,Line_List2,Line_List3,Line_List4);
		Grid_Flow_Over_Cylinder.Generate_Grid(Is_Periodic_Boundary,Enable_Elliptic_Solver,iterations);
		Grid_Plane_List = Grid_Flow_Over_Cylinder.Get_Grid_List();
	
	
	 	Identify_Cells(Total_Points_X,Nop_Y);
//	 	Identify_Neighbours(Total_Points_X,Nop_Y);
		Identify_Cells_and_Assign_Neighbours(Nop_X,Nop_Y);			
	 	write_VTK(Total_Points_X,Nop_Y,Grid_Plane_List,Grid_View_File);
	 	write_inputfile(Total_Points_X,Nop_Y,Grid_Plane_List,Solver_Input_File);	

		break;
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
		Grid_Flow_Over_Cylinder.Generate_Grid(Is_Periodic_Boundary,Enable_Elliptic_Solver,iterations);
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
//	 	Identify_Neighbours(Nop_X,Nop_Y);
		Identify_Cells_and_Assign_Neighbours_Custom(Nop_X,Nop_Y);
	 	write_VTK(Nop_X,Nop_Y,Grid_Plane_List,Grid_View_File);
	 	write_inputfile(Nop_X,Nop_Y,Grid_Plane_List,Solver_Input_File);	
		
		break;
	}
	
		
return 0;
}

