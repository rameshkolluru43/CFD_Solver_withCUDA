#include "headers.hpp"

Line::Line()
{
        Start_Point(0.0,0.0);
        End_Point(0.0,0.0);
}

Point Line::get_Point(int i)
{
        return point_list[i];
}

void Line::generate(Point sp,Point ep,int n)
{
        Point Temp_point;
        sp.Print();
        ep.Print();
        delx=(sp.Get_x()-ep.Get_x())/(n-1);
        dely=(sp.Get_y()-ep.Get_y())/(n-1);
        for(int i=0;i<n;i++)
        {
                Temp_point(sp.Get_x()+i*delx,sp.Get_y()+i*dely);
                point_list.push_back(Temp_point);
        }
        cout<<point_list.size()<<endl;
}

void Line::write_output()
{
        Point temp_Point;
        ofstream myfileout("line.txt",ios::app|ios::binary); 
        if(myfileout.is_open())
        {
  //              myfileout<<"POINTS "<<point_list.size()<<endl;
                for(int i=0;i<point_list.size();i++)
                {
                        temp_Point=get_Point(i);
                        temp_Point.Print();
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

vector<Point > Line::Get_Point_list()
{
        return point_list;
}



Circle::Circle()
{
        origin(0.0,0.0);
        radius=0.0;
        theta = 0.0;
}

/*Circle::~Circle(){}*/

void Circle::generate(double radius)
{
        Point temp_Point;
        no_of_points=0.0;
        for(theta=0.0;theta<360;theta++)
        {
                 temp_Point(radius,theta,1);
                 point_list.push_back(temp_Point);
                 no_of_points+=1;
        }
}

void Circle::operator() (double radius)
{
        Point temp_Point;
        no_of_points=0.0;
        for(theta=0.0;theta<360;theta++)
        {
                temp_Point(radius,theta,1);
                point_list.push_back(temp_Point);
                no_of_points+=1;
        }
}

void Circle::generate(double radius,int n)
{
        Point temp_Point;
        no_of_points=0.0;
        cout<<n<<endl;
        for(theta=0.0;theta<360.0;theta+=(360.0/n))
        {
                temp_Point(radius,theta,1);
                point_list.push_back(temp_Point);
                no_of_points+=1;
        }
        cout<<point_list.size()<<endl;
}

void Circle::operator()(double radius,int n)
{
        Point temp_Point;
        no_of_points=0.0;
        cout<<n<<endl;
        for(theta=0.0;theta<360.0;theta+=(360.0/n))
        {
                temp_Point(radius,theta,1);
                //cout<<theta+i<<"\t"; temp_Point.Print();
                point_list.push_back(temp_Point);
                no_of_points+=1;
        }
        cout<<point_list.size()<<endl;
}

void Circle::generate(double radius1, double radius2,int n)
{
        generate(radius1,n);
        generate(radius2,n);
}

void Circle::operator() (double radius1, double radius2,int n)
{
        generate(radius1,n);
        generate(radius2,n);
}

Point Circle::get_Point(int i)
{
        return point_list[i];
}

int Circle::get_Nop()
{
        return no_of_points;
}

void Circle::write_output(string filename)
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

vector<Point> Circle::Get_Point_list()
{
        return point_list;
}


Object::Object()
{
        np_c=10;
        np_l=np_c;
        radius =1.0;
}

void Object::generate(Circle cir1,Circle cir2,string filename)
{
        
        c1=cir1;
        cplist1=c1.Get_Point_list();
        c2=cir2;
        cplist2=c2.Get_Point_list();
        write_output(cplist1,cplist2,filename);
        cout<<"object generated"<<endl;
}

void Object::operator() (Circle cir1,Circle cir2,string filename)
{
        c1=cir1;
        cplist1=c1.Get_Point_list();
        c2=cir2;
        cplist2=c2.Get_Point_list();
        write_output(cplist1,cplist2,filename);
        cout<<"object generated"<<endl;
}

void Object::operator() (Line line1,Line line2,string filename)
{
	l1=line1;
	l2=line2;
	lplist1=l1.Get_Point_list();
	lplist2=l2.Get_Point_list();
	cout<<"am here\n";
	//lplist1.size();
	//lplist2.size();
	write_output(lplist1,lplist2,filename);
        cout<<"object generated"<<endl;
}

void Object::generate()
{
}
void Object::generate(Circle cir1,string filename)
{
        c1=cir1;
        p1=c1.get_Point(0);
        cplist1=c1.Get_Point_list();
        write_output(cplist1,filename);
        cout<<"Object generated and written"<<endl;
}

void Object::write_output(vector<Point> list1,string filename)
{
        ofstream myfileout;
        myfileout.open(filename.c_str(),ios::app|ios::binary); 
        Point temp_Point;
        if(myfileout.is_open())
        {
                myfileout<<1<<endl;
                myfileout<<list1.size()/4<<"\t"<<list1.size()/4<<endl;
                for(int i=0;i<list1.size();i++)
                {
                        temp_Point=c1.get_Point(i);
                        myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<endl;
                }
        }
        else
        {
                cout<<" unble to open output file\n"; 
        }
        myfileout.close();
}

void Object::write_output(vector<Point> list1,vector<Point> list2,string filename)
{
        ofstream myfileout;
        myfileout.open(filename.c_str(),ios::app|ios::binary); 
        Point temp_Point;
        if(myfileout.is_open())
        {
                myfileout<<2<<endl;
                myfileout<<list1.size()<<"\t"<<list2.size()<<endl;
                for(int i=0;i<list1.size();i++)
                {
                        temp_Point=l1.get_Point(i);
                        myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<endl;
                       // cout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<endl;
                }
                for(int i=0;i<list2.size();i++)
                {
                        temp_Point=l2.get_Point(i);
                        myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<endl;
                        //cout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<endl;
                }
        }
        else
        {
                cout<<" unble to open output file\n"; 
        }
        myfileout.close();
}

Cylinder::Cylinder()
{

}

void Cylinder::operator()(double diameter,int length,int nop,int nopz)
{
        C1(diameter/2.0,nop); // creates a circle with number of points nop
	plist1=C1.Get_Point_list();
	Cylinder_Length=length;
	nr=ntheta=nop;
	nz=nopz;
}

void Cylinder::operator()(double d1,double d2,double length,int nx,int ny,int nl,string filename)
{
	nr=nx;
	ntheta=ny;
	nz=nl;
	C1(d1/2.0,nr); // creates a circle with number of points nop
	C2(d2/2.0,nr);
	Cylinder_Length = length;
	plist1=C1.Get_Point_list();
	plist2=C2.Get_Point_list();
	write_output(plist1,plist2,filename);
}

void Cylinder::Grid_Cylinder(string ipfile,string grid_plane_file,string opfile)
{
	grid.Generate_Grid(ipfile);
	Grid_Plane = grid.Get_Grid_List();
	cout<<"Size of Grid Plane\t"<<Grid_Plane.size()<<endl;
	write_output(Grid_Plane,grid_plane_file);
	cout<<Cylinder_Length<<"\t"<<nz<<endl;
	Cylinder_Grid=grid.Stack_Grid(Grid_Plane,Cylinder_Length,nz);
	cout<<"Size of Grid Cylinder \t"<<Cylinder_Grid.size()<<endl;
	write_output(Cylinder_Grid,opfile);
}

void Cylinder::write_output(vector<Point> list,string filename)
{
        ofstream myfileout;
        myfileout.open(filename.c_str(),ios::out);
        Point temp_Point;
	cout<<"Writing Grid Plane vtk file\n ";
	cout<<nr<<"\t"<<ntheta<<"\t"<<nz<<endl;
	myfileout<<nr<<"\t"<<ntheta<<"\t"<<nz<<endl;
	cout<<"Grid plane size\t"<<Grid_Plane.size()<<endl;
	myfileout<<Grid_Plane.size()<<endl;
	cout<<"List size\t"<<list.size()<<endl;
	myfileout<<list.size()<<endl;
        if(myfileout.is_open())
        {
                for(int i=0;i<list.size();i++)
                {
                        temp_Point=list[i];
                        myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<"\t"<<temp_Point.Get_z()<<endl;
                }
        }
        else
        {
                cout<<" unble to open output file\n"; 
        }
        myfileout.close();
}

void Cylinder::write_output(vector<Point> list1,vector<Point> list2,string filename)
{
        ofstream myfileout;
        myfileout.open(filename.c_str(),ios::out);
        Point temp_Point;
        if(myfileout.is_open())
        {
		myfileout<<2<<endl;
		myfileout<<nr<<"\t"<<ntheta<<endl;
                for(int i=0;i<list1.size();i++)
                {
                        temp_Point=list1[i];
                        myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<"\t"<<temp_Point.Get_z()<<endl;
                }
                for(int i=0;i<list2.size();i++)
                {
                        temp_Point=list2[i];
                        myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<"\t"<<temp_Point.Get_z()<<endl;
                }
        }
        else
        {
                cout<<" unble to open output file\n";
        }
        myfileout.close();
}
