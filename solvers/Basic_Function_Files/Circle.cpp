#include "Geometry_Header.h"




/*Constructor function which creates a circle object of radius =0 and angle =0.0*/
Circle::Circle()
{
        radius=0.0;
        theta = 0.0;
}

// Destructor function
Circle::~Circle()
{
}

// Function creates a circle with points at radius 
void Circle::generate(double & radius)
{
        Point temp_Point;
        no_of_points=0.0;
        for(theta=0.0;theta<=360;theta++)
        {
                 temp_Point(radius,theta,1);
                 point_list.push_back(temp_Point);
                 no_of_points+=1;
        }
}

void Circle::operator() (double &radius)
{
        Point temp_Point;
        no_of_points=0.0;
        for(theta=0.0;theta<=360;theta++)
        {
                temp_Point(radius,theta,1);
                point_list.push_back(temp_Point);
                no_of_points+=1;
        }
}

void Circle::generate(double &radius,int &n)
{
        Point temp_Point;
        no_of_points=n;
        cout<<n<<endl;
        for(theta=0.0;theta<=360.0;theta+=(360.0/(n-1)))
        {
                temp_Point(radius,theta,1);
                point_list.push_back(temp_Point);
                no_of_points+=1;
        }
        cout<<point_list.size()<<endl;
}

void Circle::operator()(double &radius,int &n)
{
        Point temp_Point;
        no_of_points=0.0;
         cout<<n<<endl;
        for(theta=0.0;theta<=360.0;theta+=(360.0/(n-1)))
        {
                temp_Point(radius,theta,1);
//                 cout<<theta<<"\t"; temp_Point.Print();
                point_list.push_back(temp_Point);
                no_of_points+=1;
        }
//         cout<<point_list.size()<<endl;
}

void Circle::generate(double &radius1, double &radius2,int &n)
{
        generate(radius1,n);
        generate(radius2,n);
}

void Circle::operator() (double &radius1, double &radius2,int &n)
{
        generate(radius1,n);
        generate(radius2,n);
}

Point & Circle::get_Point(int &i)
{
        return point_list[i];
}

int & Circle::get_Nop()
{
        return no_of_points;
}

void Circle::write_output(string &filename)
{
        Point temp_Point;
        //cout<<point_list.size()<<endl;
        double nx,ny;
        nx=ny=point_list.size();
        ofstream myfileout(filename.c_str(),ios::app|ios::binary);
        if(myfileout.is_open())
        {
	//	myfileout<<"POINTS "<<point_list.size()<<endl;
		for(int i=0;i<point_list.size();i++)
		{
			temp_Point=get_Point(i);
			myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y();
			myfileout<<"\t"<<0.0<<endl;
		}
        }
        else
        {
                cout<<" unble to open output file\n"; 
        }
        myfileout.close();
}

vector<Point>& Circle::Get_Point_list()
{
        return point_list;
}

void Circle::Reverse_Points()
{
	int size = point_list.size();
	vector<Point> Reverselist;
	for(int i = size-1;i>=0;i--)
		Reverselist.push_back(point_list[i]);
	point_list = Reverselist;
}

void Circle::Print()
{
	for(unsigned int i=0;i<point_list.size();i++)
	{
		point_list[i].Print();
		cout<<endl;
	}
}
