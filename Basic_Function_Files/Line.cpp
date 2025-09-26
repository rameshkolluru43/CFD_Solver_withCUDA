#include "Geometry_Header.h"

Line::Line()
{
        Start_Point(0.0,0.0,0.0);
        End_Point(0.0,0.0,0.0);
		ds_x=0.0;
		ds_y=0.0;
		ds_z=0.0;
		Length = 0.0;
}

const Point & Line::get_Point(int i) const
{
        return point_list[i];
}

void Line::generate(Point & sp,Point & ep,int & n)
{
        Point Temp_point;
		no_of_points = n;
        sp.Print();
        ep.Print();
        
		ds_x=(ep.Get_x()-sp.Get_x());
        ds_y=(ep.Get_y()-sp.Get_y());
		ds_z=(ep.Get_z()-sp.Get_z());
		
		Length = sqrt(ds_x*ds_x + ds_y*ds_y + ds_z*ds_z);
		
        delx=(ep.Get_x()-sp.Get_x())/(n-1);
        dely=(ep.Get_y()-sp.Get_y())/(n-1);
		delz=(ep.Get_z()-sp.Get_z())/(n-1);
		
        for(int i=0;i<n;i++)
        {
                Temp_point(sp.Get_x()+i*delx,sp.Get_y()+i*dely,sp.Get_z()+i*delz);
                point_list.push_back(Temp_point);
        }
        cout<<"Created Line\t"<<point_list.size()<<endl;
}

void Line::generate(Point & sp,Point & ep,int & n,bool & both_sides)
{
		double beta=0.0,alpha=0.0,a=0.0,b=0.0,y=0.0,eta=0.0;
        Point Temp_point;
		no_of_points = n;
		cout<<"Starting Point\t";
        sp.Print();
		cout<<"Ending Point\t";
        ep.Print();
		
		ds_x=(ep.Get_x()-sp.Get_x());
        ds_y=(ep.Get_y()-sp.Get_y());
		ds_z=(ep.Get_z()-sp.Get_z());
		
		Length = sqrt(ds_x*ds_x + ds_y*ds_y + ds_z*ds_z);
		
		cout<<"Length of the line\t"<<Length<<endl;
        delx=(ep.Get_x()-sp.Get_x())/(n-1);
        dely=(ep.Get_y()-sp.Get_y())/(n-1);
		delz=(ep.Get_z()-sp.Get_z())/(n-1);
		
		
		del_s= Length/(n-1);
		
		for(int j=0;j<n;j++)
		{	
			both_sides = true;
			
			if(both_sides)
			{
				beta = 1.01;
				alpha = 0.5;
				eta = j*del_s/Length;
				a = (eta - alpha)/(1.0 - alpha);
				b = (beta + 1.0)/(beta - 1.0);
				y = Length*((((2*alpha + beta)*pow(b,a)) + 2*alpha - beta)/((2*alpha + 1.0)*(pow(b,a)+1.0)));
			}
			else
			{

				beta = 1.01;
				alpha = 0.01;
				eta = j*del_s/Length;
				a = (1.0-eta);
				b = (beta + 1.0)/(beta - 1.0);
				y = Length*((beta + 1.0) - (beta - 1.0)*pow(b,a))/(pow(b,a) + 1.0);
//				cout<<"y=\t"<<y<<endl;
			}
			if(ds_x == 0.0 and ds_z == 0.0)
				Temp_point(sp.Get_x()+j*ds_x,sp.Get_y()+y ,sp.Get_z()+j*ds_z);
			else if(ds_x == 0.0 and ds_y == 0.0)
				Temp_point(sp.Get_x()+j*ds_x,sp.Get_y()+j*ds_y,sp.Get_z()+y);
			else if(ds_z == 0.0 and ds_y == 0.0)
				Temp_point(sp.Get_x()+y,sp.Get_y()+j*ds_y ,sp.Get_z()+j*ds_z);
			else if(ds_x !=0.0 and ds_y !=0.0 and ds_z ==0.0)
				Temp_point(sp.Get_x() + y,sp.Get_y() + y,sp.Get_z()+j*ds_z);
			else if(ds_x !=0.0 and ds_z !=0.0 and ds_y ==0.0)
				Temp_point(sp.Get_x()+y,sp.Get_y()+j*ds_y,sp.Get_z()+y);
			else if(ds_y !=0.0 and ds_z !=0.0 and ds_x ==0.0)
				Temp_point(sp.Get_x()+j*ds_x,sp.Get_y()+y,sp.Get_z()+y);
			else 
				Temp_point(sp.Get_x() + y,sp.Get_y()+y,sp.Get_z()+y);
	
//				Temp_point.Print();
                point_list.push_back(Temp_point);
		}
}

void Line::Reverse_Points()
{
	int size = point_list.size();
	vector<Point> Reverselist;
	for(int i = size-1;i>=0;i--)
		Reverselist.push_back(point_list[i]);
	point_list = Reverselist;
}

void Line::Merge(Line & Line1)
{
	cout<<"Merging Two lInes\n";
	vector<Point> Temp_Point_List;
	Temp_Point_List = Line1.Get_Point_list();
 	cout<<point_list.size()<<"\t"<<Temp_Point_List.size()<<endl;
//  	point_list[no_of_points-1].Print();
//  	Temp_Point_List[0].Print();
	if(point_list[no_of_points-1]==Temp_Point_List[0])
	{
		for(unsigned int i=1;i<Temp_Point_List.size();i++)
			point_list.push_back(Temp_Point_List[i]);
	}
	else
	{
		for(unsigned int i=0;i<Temp_Point_List.size();i++)
			point_list.push_back(Temp_Point_List[i]);
	}
}

void Line::Merge(vector<Point> & List1)
{
	cout<<"Merging List of Points to current Line\n";
 	cout<<point_list.size()<<"\t"<<List1.size()<<endl;
//  	point_list[no_of_points-1].Print();
//  	Temp_Point_List[0].Print();
	if(point_list[no_of_points-1]==List1[0])
	{
		for(unsigned int i=1;i<List1.size();i++)
			point_list.push_back(List1[i]);
	}
	else
	{
		for(unsigned int i=0;i<List1.size();i++)
			point_list.push_back(List1[i]);
	}
}

void Line::write_output()
{
        Point temp_Point;
        ofstream myfileout("line.txt",ios::app|ios::binary); 
        if(myfileout.is_open())
        {
  //              myfileout<<"POINTS "<<point_list.size()<<endl;
                for(unsigned int i=0;i<point_list.size();i++)
                {
                        temp_Point=get_Point(i);
                        temp_Point.Print();
			myfileout<<temp_Point.Get_x()<<"\t"<<temp_Point.Get_y()<<"\t"<<temp_Point.Get_z();
                        myfileout<<"\t"<<0.0<<endl;
                }
        }
        else
        {
                cout<<" unble to open output file\n"; 
        }
        myfileout.close();
}

const vector<Point> & Line::Get_Point_list() const
{
        return point_list;
}

void Line::Print()
{
	for(unsigned int i=0;i<point_list.size();i++)
	{
		point_list[i].Print();
	}
	cout<<"---------------------------------\n";
}

double  Line::Size()
{
	return point_list.size();
}