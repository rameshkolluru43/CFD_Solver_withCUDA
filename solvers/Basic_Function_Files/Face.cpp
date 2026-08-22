//******************************************************************************************************
//	Implementation of Face Class
//	Face has two objectives
//		1. Geometrical Feature, has Face Area, Face Center, Face Normal,
//		2. Property, It calculates the flux from the face called by the cell
//******************************************************************************************************

#include "headers.hpp"

Face::Face()
{
	Flux.resize(5,0.0);Del_Q_1.resize(5,0.0);Del_Q_2.resize(5,0.0);Del_Q_3.resize(5,0.0);
	Area_Comp1=0.0;Area_Comp2=0.0;Area_Comp3=0.0;
}

Face::~Face()
{
}


//overloaded ()-operator, takes 4 points and creates a face, calculates the face area and face normal, which are the cross product of the diagonal vectors formed by the diagonal points.
void Face::operator()(const Point& a,const Point& b,const Point& c,const Point& d,const int & face_no)
{
//  	cout<<"Face with 4 points\n";
	Face_No = face_no;
	face_nop=4;
	Point1=a;Point2=b;Point3=c;Point4=d;
// 	a.Print();b.Print();c.Print();d.Print();
// 	Point1.Print();Point2.Print();Point3.Print();Point4.Print();
	Vector Av1,Av2;
	Av1(Point1,Point3);Av2(Point2,Point4);
	Area_Vector=(Av1.Cross(Av2))*0.5;
	Face_area = Area_Vector.magnitude();
	Area_Comp1=Area_Vector(1);
	Area_Comp2=Area_Vector(2);
	Area_Comp3=Area_Vector(3);
	Face_Normal=Area_Vector.Find_UnitVector();
// 	Area_Vector->Print();
// 	Face_Normal->Print();
	Face_Center();
}

void Face::operator()(const Point& a,const Point& b,const Point& c,const int & face_no)
{
//  	cout<<"Face with 3 points\n";
	Face_No = face_no;
	face_nop=3;
        Point1=a;Point2=b;Point3=c;
	Vector V1,V2;
	V1(Point1,Point2);V2(Point1,Point3);
	Area_Vector=(V1.Cross(V2)*0.5);
	Face_area = Area_Vector.magnitude();
	Area_Comp1=Area_Vector.Get_Component(1);Area_Comp2=Area_Vector.Get_Component(2);Area_Comp3=Area_Vector.Get_Component(3);
	Face_Normal=Area_Vector.Find_UnitVector();
// 	Face_Normal.Print();
	Face_Center();
}

void Face::Face_Center()
{
        if(face_nop==4)
        {
		Center_Point=(Point1+Point2+Point3+Point4)*0.25;
        }
        else
        {
		Center_Point=(Point1+Point2+Point3)/3.0;
        }
	R_mid(Center_Point);
}

void Face :: Form_Centroid_Vector(Point & P1,Point & P2)
{
	double cv_mag=0.0;
	Centroid_Vector(P1,P2);
	cv_mag=Centroid_Vector.magnitude();
	Centroid_Normal = Centroid_Vector.Find_UnitVector();
}

const Face& Face::operator=(const Face& face1)
{
	if(this == &face1)
		return *this;
	Point1=face1.Point1;Point2=face1.Point2;Point3=face1.Point3;Point4=face1.Point4;
	face_nop=face1.face_nop;
	Area_Vector=face1.Area_Vector;
	Area_Comp1=Area_Vector.Get_Component(1);Area_Comp2=Area_Vector.Get_Component(2);Area_Comp3=Area_Vector.Get_Component(3);
	Face_Normal=face1.Face_Normal;
	Center_Point=face1.Center_Point;
	R_mid=face1.R_mid;
	return *this;
}

const Vector& Face::Get_FaceNormal() const
{
	return Face_Normal;
}

const Vector& Face::Get_FaceArea_Vector() const
{
	return Area_Vector;
}

const Vector& Face::Get_Rmid() const
{
	return R_mid;
}

const double& Face::Get_Face_Area() const
{
	return Face_area;
}

const int& Face::Get_FaceNop() const 
{
	return face_nop;
}


double Face::Get_RmiddotA() const
{
	return R_mid*Area_Vector;
}


//calculates p,T,v1,v2,v3,rho,a,M when Q is given
void Face::Cal_Primitive(const vector<double>& Qt)
{
/*	cout<<"In Face cal Primitive function\n";
	cout<<Qt[0]<<"\t"<<Qt[1]<<"\t"<<Qt[2]<<"\t"<<Qt[3]<<"\t"<<Qt[4]<<endl;*/
	double temp=0.0;
	Face_Rho = Qt[0];
	temp = 1.0/Face_Rho;
	v1 = Qt[1]*temp;
	v2 = Qt[2]*temp;
	v3 = Qt[3]*temp;
	Face_Velocity(v1,v2,v3);
	Face_ET = Qt[4]*temp;
	Face_T = (Face_ET-((v1*v1+v2*v2+v3*v3)*0.5))/cv;
	Face_P = Face_Rho*R*Face_T;
/*	cout<<Face_P<<"\t"<<Face_Rho<<"\t"<<Face_T<<"\t"<<v1<<"\t"<<v2<<"\t"<<v3<<endl;
	Face_Velocity.Print();*/
}


const vector<double>& Face::Flux_From_Face_Central(const vector<double> & Q1,const vector<double> &  Q2)
{
	double vnds=0.0;
	vector<double> Q(5,0.0);
	Flux[0]=0.0;Flux[1]=0.0;Flux[2]=0.0;Flux[3]=0.0;Flux[4]=0.0;
/*	cout<<Q1[0]<<"\t"<<Q1[1]<<"\t"<<Q1[2]<<"\t"<<Q1[3]<<"\t"<<Q1[4]<<"\n";
	cout<<Q2[0]<<"\t"<<Q2[1]<<"\t"<<Q2[2]<<"\t"<<Q2[3]<<"\t"<<Q2[4]<<"\n";*/
// averge of Q1 and Q2 , current and Neighbouring cells
	Q[0]=(Q1[0] + Q2[0])*0.5;
	Q[1]=(Q1[1] + Q2[1])*0.5;
	Q[2]=(Q1[2] + Q2[2])*0.5;
	Q[3]=(Q1[3] + Q2[3])*0.5;
	Q[4]=(Q1[4] + Q2[4])*0.5;
/*	cout<<Q[0]<<"\t"<<Q[1]<<"\t"<<Q[2]<<"\t"<<Q[3]<<"\t"<<Q[4]<<"\n";
	cout<<"Area Vector\t";Area_Vector.Print();*/
	Cal_Primitive(Q);
	vnds=v1*Area_Comp1 + v2*Area_Comp2 + v3*Area_Comp3 ;//Face_Velocity*Area_Vector;
//	Flux Determines the Convective Fluxes on the face
	Flux[0]=Face_Rho*vnds ;
	Flux[1]=Face_Rho*v1*vnds + Face_P*Area_Comp1;
	Flux[2]=Face_Rho*v2*vnds + Face_P*Area_Comp2;
	Flux[3]=Face_Rho*v3*vnds + Face_P*Area_Comp3;
	Flux[4]=(Face_Rho*Face_ET + Face_P)*vnds;
// 	cout<<Flux[0]<<"\t\t"<<Flux[1]<<"\t\t"<<Flux[2]<<"\t\t"<<Flux[3]<<"\t\t"<<Flux[4]<<"\n";
//	 cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	return Flux;
}

const vector<double>& Face::Flux_From_Face_Central(const double &P,const double &T,const Vector& V)
{
	double vnds=0.0;
	Flux[0]=0.0;Flux[1]=0.0;Flux[2]=0.0;Flux[3]=0.0;Flux[4]=0.0;
	Face_Rho = P/(R*T);
	Face_P = P;
	Face_T = T;
	Face_Velocity = V;
	v1=V.Get_Component(1);v2=V.Get_Component(2);v3=V.Get_Component(3);
// 	cout<<Face_P<<"\t"<<Face_T<<"\t"<<Face_Rho<<"\t";
// 	Face_Velocity.Print();	Area_Vector.Print();
	Face_a = sqrt(gamma_R*Face_T);
	max_eigen_value = (fabs(Face_Velocity*Area_Vector) + fabs(Face_a * Face_area));
	vnds = v1*Area_Comp1 + v2*Area_Comp2 + v3*Area_Comp3 ;//Face_Velocity*Area_Vector;
	Flux[0]= Face_Rho*vnds;
	Flux[1]= Face_Rho*v1*vnds  + Face_P*Area_Comp1;
	Flux[2]= Face_Rho*v2*vnds  + Face_P*Area_Comp2;
	Flux[3]= Face_Rho*v3*vnds  + Face_P*Area_Comp3;
	Face_ET = (cv*Face_T + (Face_Velocity*Face_Velocity)*0.5);
	Flux[4]= (Face_Rho*Face_ET+ Face_P)*vnds;
// 	cout<<Flux[0]<<"\t\t"<<Flux[1]<<"\t\t"<<Flux[2]<<"\t\t"<<Flux[3]<<"\t\t"<<Flux[4]<<"\n";
//   cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	return Flux;
}

void Face::Find_Gradients()
{
	undS=Scalar_dot_ndS(v1);
	vndS=Scalar_dot_ndS(v2);
	wndS=Scalar_dot_ndS(v3);
	TndS=Scalar_dot_ndS(Face_T);
	rho_ndS = Scalar_dot_ndS(Face_Rho);
	rhou_ndS = Scalar_dot_ndS(Face_Rho*v1);
	rhov_ndS = Scalar_dot_ndS(Face_Rho*v2);
	rhow_ndS = Scalar_dot_ndS(Face_Rho*v3);
	rhoET_ndS = Scalar_dot_ndS(Face_Rho*Face_ET);
	Viscosity();
	Thermal_Conductivity();
}

const  Vector& Face::Get_Face_undS() const
{
	return undS;
}
const Vector& Face::Get_Face_vndS() const 
{
	return vndS;
}
const Vector& Face::Get_Face_wndS() const 
{
	return wndS;
}

const Vector& Face::Get_Face_TndS() const
{
	return TndS;
}
const double& Face::Get_Face_Pressure() const
{
	return Face_P;
}
void Face::Set_Face_Pressure(const double & p)
{
	Face_P=p;
}
const double& Face::Get_Face_Temparature() const 
{
	return Face_T;
}
void Face::Set_Face_Temparature(const double & T)
{
	Face_T = T;
}
const Vector& Face::Get_Face_Velocity() const
{
	return Face_Velocity;
}

void Face::Viscosity()
{
	mu  = (V_C1*pow(Face_T,1.5))/(Face_T + V_S_T);
}

void Face::Thermal_Conductivity()
{
	K   = (T_C1*pow(Face_T,1.5))/(Face_T + T_S_T);
}
