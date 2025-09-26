#include "headers.hpp"

Point::Point()
{
	x=0.0; y=0.0; z=0.0;
	r = 0.0; theta = 0.0; phi = 0.0;
}

Point::Point(const double& a , const double & b , const double & c)
{
        x=a;
        y=b;
        z=c;
// 	Convert_Cylindrical();
//      Convert_Spherical();
}

void Point::operator () (const double & a , const double & b)
{
	x=a;
	y=b;
	z=0.0;
	Convert_Polar();
}

void Point::operator () (const double & a , const double & b ,const int& c) //
{
	r=a;
	theta =b;
	Convert_Cartesian();
}

void Point::operator () (const double & a , const double & b , const double & c)
{
	x=a;
	y=b;
	z=c;
// 	Convert_Spherical();
}

void Point::operator () (const double & a , const double & b ,const double & c ,const int& s)
{
	r=a;
	theta =b;
	phi=c;
	Convert_Cartesian(1);
}

const double& Point::Get_x() const
{
	return x;
}

const double& Point::Get_y() const
{
	return y;
}

const double& Point::Get_z() const 
{
	return z;
}

const double& Point::Get_r() const 
{
	return r;
}

double& Point::Get_theta()
{
	theta= radian*theta;
	return theta;
}

double& Point::Get_phi()
{
	phi = (radian*phi);
	return phi;
}

void Point::Convert_Polar()
{
	r=sqrt(x*x + y*y);
	theta = atan(y/x);
}

void Point::Convert_Cartesian()
{
// 	cout<<r<<"\t"<<theta<<"\t"<<radian<<"\t"<<radian*theta<<endl;
	double rtheta=radian*theta;
	x = r*cos(rtheta);
	y = r*sin(rtheta);
}

void Point::Convert_Cartesian(const int& i)
{
	double rtheta,rphi;
	rtheta=radian*theta;
	rphi=radian*phi;
	x = r*sin(rtheta)*cos(rphi);
	y = r*sin(rtheta)*sin(rphi);
	z = r*cos(rtheta);
}

void Point::Convert_Spherical()
{
	r = sqrt(x*x + y*y + z*z);
	phi = atan(y/x);
	double r1=sqrt(x*x + y*y); 
	theta = atan(r1/z);
}

void Point::Print() const
{
	cout<<"("<<x<<","<<y<<","<<z<<")\n";
	//cout<<"("<<r<<","<<theta<<","<<phi<<")"<<endl;
}

Point::~Point()
{
// 	cout<<"Point Destructed \n";
}

Point& Point::operator= (const Point& p1)
{
	if(this == &p1)
		return *this;
        x=p1.x;
        y=p1.y;
        z=p1.z;
        return *this;
}

const Point& Point::operator-=(const Point& p1)
{
        x-=p1.x;
        y-=p1.y;
        z-=p1.z;
        return *this;
}

Point Point::operator-(const Point& p1) const
{
	Point Temp;
	Temp =*this;
	Temp-=p1;
	return Temp;
}

const Point& Point::operator +=(const Point& p1)
{
	x+=p1.x;
	y+=p1.y;
	z+=p1.z;
	return *this;
}

Point Point::operator +(const Point& p1) const
{
	Point Temp;
	Temp =*this;
	Temp+=p1;
	return Temp;
}

const Point& Point::operator /=(const double& alpha)
{
	x/=alpha;
	y/=alpha;
	z/=alpha;
	return *this;
}

Point Point::operator /(const double& alpha) const
{
	Point Temp;
	Temp=*this;
	Temp/=alpha;
	return Temp;
}
const Point& Point::operator *=(const double& alpha)
{
	x*=alpha;
	y*=alpha;
	z*=alpha;
	return *this;
}

Point Point::operator *(const double& alpha) const
{
	Point Temp;
	Temp=*this;
	Temp*=alpha;
	return Temp;
}


bool operator==(const Point& p1,const Point& p2)
{
//  	p1.Print();
//  	p2.Print();
	cout<<"Printing the differences in points\t"<<p1.x-p2.x<<"\t"<<p1.y-p2.y<<"\t"<<p1.z-p2.z<<endl;
        if(	((p1.x-p2.x)<1e-10)		&&		((p1.y-p2.y)<1e-10)		&&		((p1.z-p2.z)<1e-10)		)
                return true;
        else
                return false;
}

bool operator!= (const Point& p1,const Point& p2)
{
        return !(p1==p2);
}