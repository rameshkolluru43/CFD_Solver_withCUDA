#ifndef _Geometry_Header_H
#define _Geometry_Header_H


#include "headers.hpp"

class Shape
{
	public:
		Shape(){ };
		virtual void generate(){};
		virtual void write_output(){};
//              virtual  Get_Point_list(){};
	private:
		Point p;
};

class Line : public Shape
{
	public:
		Line();
		using Shape::generate;
		using Shape::write_output;
		void generate(Point &,Point &,int &);
		void generate(Point &,Point &,int &,bool &);
		void Merge(Line &);
		void Merge(vector<Point> &);
		void write_output();
		const Point & get_Point(int) const;
		const vector<Point> & Get_Point_list() const;
		void Print();
		double Size();
		void Reverse_Points();
	private:
		vector<Point> point_list;
		int no_of_points;
		double delx,dely,delz,Length,ds_x,ds_y,ds_z,del_s;
		Point Start_Point,End_Point;
};

class Circle : public Shape
{
	public:
		Circle();
		~Circle();
		using Shape::generate;
		using Shape::write_output;
		void operator()(double &);
		void operator()(double &,int &);
		void operator()(double &,double &,int &);
		void operator()(double &,double &,double &,int &);
		void generate(double &);
		void generate(double &,int &);
		void generate(double &, double &,int &);
		void write_output(string &);
		Point & get_Point(int  &);
		vector<Point> & Get_Point_list();
		int & get_Nop();
		void Print();
		void Reverse_Points();
	private:
		Point origin;
		vector<Point> point_list;
		double radius,theta;
		int no_of_points,nloops;
		vector<double > x,y;
};

class Object : public Shape
{
	public:
		Object();
		using Shape::generate;
		using Shape::write_output;
		void generate();
		void generate(Circle,Circle,string);
		void operator() (Circle,Circle,string);
		void operator() (Line,Line,string);
		void generate(Circle,string);
		void write_output(vector<Point>,string);
		void write_output(vector<Point>,vector<Point>,string);
		//                vector<Point> Get_Point_list();
	private:
		Circle c1,c2;
		Line l1,l2;
		int np_c,np_l;
		double radius;
		Point p1,p2;
		vector<Point> object_point_list,cplist1,cplist2,lplist1,lplist2;
};

class Grid
{
	public:
		Grid();
		void operator()(Circle & );
		void operator()(vector<Point> &, vector<Point> &,vector<Point> &, vector<Point> &);
		void Read_data(string);
		void Generate_Grid_List();
		void Generate_Grid(bool &);
		void Generate_Grid(bool & ,bool & ,int &);
		void TFI();
		void Elliptic();
		void Boundary_Condition();
		vector<Point>& Get_Grid_List();
		void Stack_Grid(double & ,int &,vector<Point> &);
		void Estimate_Error_and_Update();
		void Periodic_Boundary_Condition();
	private:
		vector< vector <double>  > x,y,xtemp,ytemp,erx,ery;
		int nx,ny,nloops,noit,np1,np2,np3,np4;
		double eta,zi,t,s,xzi,xeta,yzi,yeta,alpha,beta,gamma1,A,B,C,D,inv_eta,inv_zi,ertot,ermax;
		vector<Point> Grid_List;
};

class Cylinder:public Shape
{
	public:
		Cylinder();
                ~Cylinder();
		using Shape::write_output;
		void operator()(double&,int&,int&,int&);//diameter,length
		void operator()(double&,double&,double&,int&,int&,int&,bool &,bool &);//dia1,dia2,length,nz,nr
		void Write_Grid_Plane(const string&,const bool &);
		void write_output(string&);
		void Grid_Cylinder(bool &);
		void Grid_Cylinder_Plane(bool &);
		void Cluster_Grid_Cylinder_Plane(bool &);
		void Stack_Grid(vector<Point>&,double&,int&,bool&);
		void write_VTK(const string &,const bool &);
		void write_inputfile(const string & ,const bool & );
		void Identify_Cells_Polar();
		void Identify_Cells_Cartesian();
		void Identify_Neighbours_Polar(const bool &);
		void Get_BPoint_Details(vector<int> & );
		void Identify_Neighbours_Cartesian();
		void Identify_Boundaries(const bool & ,const bool &);
		void Append_Boundary_Information(const string &,const bool &);
	private:
		vector<Point> plist1,plist2,Grid_Plane_Cartesian,Grid_Plane_Polar,Point_List,Core_Grid_List,Polar_Grid_List;
		Circle C1,C2;
		string filename;
		int nr,ntheta,nlength,no_of_Cells_Cartesian,no_of_Cells_Polar,TGCFP,TGCBaP,TGCWall;
		int nx_c,ny_c,nz_c,Total_No_Cells,cells_in_Plane_Cartesian,cells_in_Plane_Polar,Cold_Layers,Hot_Layers,I_N_L;
// 		number of points in r,theta and length wise
		double Cylinder_Length,delx,dia1,dia2,Stretch_Ratio,a;
		vector<int> C_N_DC,C_N_DP,C_P_DC,C_P_DP,IB_List,WB_List,EB_List,Indicies;
		Grid Elliptic_Core;
};


extern vector<Point> Line_List,Grid_list,grid_cube_list;
extern vector<int> Cell_Point_Data,Cell_Neighbour_Data,Inlet_Boundary_List,Wall_Boundary_List,Exit_Boundary_List,Symmetry_Symmetry_List;


void Append_Boundary_Information(const string & );
void Stack_Grid(vector<Point> &,double &,int &,vector<Point>& );

void Identify_Cells(int &,int &,int &);
void Identify_Neighbours(int &,int &,int &);

void Identify_Cells(int &,int &,int &);
void Identify_Neighbours(int &,int &,int &);

void write_inputfile(int & ,int & ,int & ,vector<Point> & ,string & );
void write_VTK(int &,int &,int &,vector<Point> &,string & );

void Identify_Cells(int &,int &);
void Identify_Neighbours(int &,int &);

void Identify_Cells(int &,int &,int &,int &);
void Identify_Neighbours(int &,int &,int &,int &);
void Identify_Cells_and_Assign_Neighbours(int &, int &);
void Identify_Cells_and_Assign_Neighbours_Custom(int & Nx, int & Ny);

void write_inputfile(int & ,int & ,vector<Point> & ,string & );
void write_VTK(int &,int &,vector<Point> &,string & );

void write_VTK(int &,int &,int &,int &,vector<Point> &,string & );
void write_inputfile(int & ,int & ,int &,int &, vector<Point> & ,string & );
#endif // _Geometry_Header_H
