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
		void generate_stretched(Point &sp, Point &ep, int &n,
								double beta_param, double alpha_param,
								bool both_sides);
		void generate_boundary_layer(Point &sp, Point &ep, int &n,
									 double first_cell_height,
									 double growth_rate);
		void generate_tanh(Point &sp, Point &ep, int &n,
						   double stretching_factor);
		void Merge(Line &);
		void Merge(vector<Point> &);
		void write_output();
		const Point & get_Point(int) const;
		const vector<Point> & Get_Point_list() const;
		void Print();
		double Size();
		void Reverse_Points();
		void Clear();
		double GetLength() const { return Length; }
	private:
		vector<Point> point_list;
		int no_of_points;
		double delx,dely,delz,Length,ds_x,ds_y,ds_z,del_s;
		Point Start_Point,End_Point;
		void distribute_points(Point &sp, const vector<double> &s_distribution);
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

struct BoundaryLayerParams
{
	double first_cell_height;
	double growth_rate;
	int    num_layers;
	double max_stretching;
	bool   enabled;
	BoundaryLayerParams() : first_cell_height(1e-4), growth_rate(1.2),
		num_layers(20), max_stretching(1.5), enabled(false) {}
};

struct SourceTermParams
{
	bool   enabled;
	double amplitude;
	double decay_rate;
	bool   use_thomas_middlecoff;
	bool   use_boundary_attraction;
	double attraction_strength;
	SourceTermParams() : enabled(false), amplitude(1.0), decay_rate(1.0),
		use_thomas_middlecoff(true), use_boundary_attraction(false),
		attraction_strength(5.0) {}
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
		void Generate_Grid(bool &periodic, bool &enableElliptic, int &iterations,
						   double omega, double convergence_tol);
		void TFI();
		void Elliptic();
		void Elliptic_Poisson();
		void Boundary_Condition();
		vector<Point>& Get_Grid_List();
		void Stack_Grid(double & ,int &,vector<Point> &);
		void Stack_Grid(double &length, int &nop, vector<Point> &out,
						double growth_ratio);
		void Estimate_Error_and_Update();
		void Estimate_Error_and_Update_SOR(double omega);
		void Periodic_Boundary_Condition();

		void SetSourceTermParams(const SourceTermParams &params);
		void SetBoundaryLayerParams(int boundary_id, const BoundaryLayerParams &params);
		void ComputeThomasMiddlecoffSources();
		void ComputeBoundaryAttractionSources();

		int    Get_nx() const { return nx; }
		int    Get_ny() const { return ny; }
		double Get_x(int i, int j) const { return x[i][j]; }
		double Get_y(int i, int j) const { return y[i][j]; }
		double Get_Convergence() const { return ertot; }
		int    Get_Iterations() const { return noit; }
		void   ComputeGridQuality(double &min_jac, double &max_ar,
								  double &max_skew, double &min_orth) const;

	private:
		vector< vector <double>  > x,y,xtemp,ytemp,erx,ery;
		vector< vector <double>  > P_source, Q_source;
		int nx,ny,nloops,noit,np1,np2,np3,np4;
		double eta,zi,t,s,xzi,xeta,yzi,yeta,alpha,beta,gamma1,A,B,C,D,inv_eta,inv_zi,ertot,ermax;
		vector<Point> Grid_List;

		SourceTermParams source_params;
		BoundaryLayerParams bl_south, bl_north, bl_east, bl_west;
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
