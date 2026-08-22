#include "headers.hpp"

Edge::Edge()
{
}

Edge ::~Edge()
{
}

Edge:: Edge(const Point& sp,const Point& ep)
{
	Point1=sp;
	Point2=ep;
	Edge_Vector(Point2.Get_x()-Point1.Get_x(),Point2.Get_y()-Point1.Get_y(),Point2.Get_z()-Point1.Get_z());
	length = Edge_Vector.magnitude();
	Center_Point=(Point1+Point2)*0.5;
}

void Edge::operator() (const Point& sp,const Point& ep)
{
	Point1=sp;
	Point2=ep;
	Edge_Vector(Point2.Get_x()-Point1.Get_x(),Point2.Get_y()-Point1.Get_y(),Point2.Get_z()-Point1.Get_z());
	length = Edge_Vector.magnitude();
	Center_Point=(Point1+Point2)*0.5;
}

const Vector& Edge::Get_EdgeVector() const
{
	return Edge_Vector;
}

const Point& Edge::Get_Point(const int& i) const
{
	switch(i)
	{
		case 1:
			return Point1;
			break;
		case 2:
			return Point2;
			break;
		default:
			cout<<"Edge has only two points check and enter either 1 or 2, returning default Point Point1"<<endl;
			return Point1;
			break;
	}
}

const double& Edge::Edge_Length() const
{
	return length;
}

const Point& Edge::Edge_Center() const
{
	return Center_Point;
}

const Point& Edge::have_Common_Point(const Edge& e1) const
{
	if((Point1==e1.Point2)||(Point1==e1.Point1))
	{
		return Point1;
	}
	else if((Point2==e1.Point1)||(Point2==e1.Point2))
	{
		return Point2;
	}
	else
		cout<<"No Common Points for edges\n";
}