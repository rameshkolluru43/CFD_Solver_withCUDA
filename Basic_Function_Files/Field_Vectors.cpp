#include "headers.hpp"

Vector Vector::Vect;

Vector::Vector()
{
	no_of_components = 3;
	comp1=0.0;comp2=0.0;comp3=0.0;
	mag=0.0;
}

Vector::~Vector(){}

void Vector::operator()(const double& a, const double& b, const double& c)
{
	comp1=a;comp2=b;comp3=c;
	has_Vector_Changed=true;
}

void Vector::operator()(const double& a, const double& b)
{
	comp1=a;comp2=b;comp3=0.0;
	has_Vector_Changed = true;
}

void Vector::operator()(const Point &p1, const Point & p2 )
{
             comp1 = (p2.Get_x() - p1.Get_x());
	     comp2 = (p2.Get_y() - p1.Get_y());
	     comp3 = (p2.Get_z() - p1.Get_z());
	     has_Vector_Changed = true;
}

void Vector::operator()(const Point &p1)
{
	comp1 = p1.Get_x();
	comp2 = p1.Get_y();
	comp3 = p1.Get_z();
	has_Vector_Changed = true;
}

const Vector& Vector::operator = (const Vector& V)
{
	if(this == &V)
		return *this;
	comp1=V.comp1;
	comp2=V.comp2;
	comp3=V.comp3;
	mag=V.mag;
	has_Vector_Changed = true;
	return *this;
}

const Vector& Vector::operator+=(const Vector& temp)
{
	comp1 += temp.comp1;
	comp2 += temp.comp2;
	comp3 += temp.comp3;
	has_Vector_Changed = true;
	return *this;
}


const Vector& Vector::operator+ (const Vector& V_temp) const
{

	Vect = *this;
	Vect+=V_temp;
	return Vect;
}

const Vector& Vector::operator-=(const Vector& V_temp)
{
	comp1-=V_temp.comp1;
	comp2-=V_temp.comp2;
	comp3-=V_temp.comp3;
	has_Vector_Changed = true;
	return *this;
}

const Vector& Vector::operator- (const Vector& V_temp) const
{
	Vect = *this;
	Vect-=V_temp;
	return Vect;
}

const Vector& Vector::operator*= (const double& alpha)
{
	comp1*=alpha;
	comp2*=alpha;
	comp3*=alpha;
	has_Vector_Changed = true;
	return *this;
}

//operator overloading performing product of a vector with a scalar ---- scaling
const Vector& Vector::operator* (const double& alpha) const
{
// 	cout<<"Sclar value\t"<<alpha<<endl;
	Vect = *this;
	Vect*=alpha;
	return Vect;
}
/*
//operator overloading for dot product of vector
double Vector::operator*(const Vector& v) const
{
	return (comp1*v.comp1+comp2*v.comp2+comp3*v.comp3);
}

*/

const Vector& Vector::Cross(const Vector& V) const
{
	double c1,c2,c3;
	c1= comp2*V.comp3-comp3*V.comp2;
	c2= -(comp1*V.comp3-comp3*V.comp1);
	c3= comp1*V.comp2-comp2*V.comp1;
	Vect(c1,c2,c3);
	return Vect;
}

const Vector& Vector::operator/= (const double& alpha)
{
	comp1=comp1/alpha;
	comp2=comp2/alpha;
	comp3=comp3/alpha;
	has_Vector_Changed = true;
//	cout<<"comp1\t"<<comp1<<"\tcomp2\t"<<comp2<<"\tcomp3\t"<<comp3<<endl;
	return *this;
}



const double& Vector::operator()(const int& i) const
{
	switch(i)
	{
		case 1:
			return comp1;
		case 2:
			return comp2;
		case 3:
			return comp3;
		default:
			cout<<"Takes 1,2 or 3 to the components please check the number, returning default component 1\n";
			return comp1;
	}
}

double& Vector::magnitude()
{
	if(has_Vector_Changed)
	{
		mag = sqrt(comp1*comp1+comp2*comp2+comp3*comp3);
		return mag;
	}
	else
		return mag;
}

//This function evaluates the unit vector of a given vector  vect a/ |a| , but you have to evaluate the magnitude of the vector before calling this function other wise this function will return inf 
const Vector& Vector::Find_UnitVector() const
{
// 	cout<<"In unit vector function\n";
	double c1,c2,c3;
	c1=comp1/mag;
	c2=comp2/mag;
	c3=comp3/mag;
	Vect(c1,c2,c3);
	return Vect;
}

void Vector::Print() const
{
	cout<<comp1<<"\t"<<comp2<<"\t"<<comp3<<endl;
}

void Vector::Clear()
{	
	comp1=0.0;
	comp2=0.0;
	comp3=0.0;
}

const double& Vector::Get_Component(const int& i) const
{
	switch(i)
	{
		case 1:
			return comp1;
		case 2:
			return comp2;
		case 3:
			return comp3;
		default:
			cout<<"Takes 1,2 or 3 to the components please check the number, returning default comp1\n";
			return comp1;
	}
}
