#ifndef _Header_H
#define _Header_H

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<typeinfo>
#include<iomanip>
//#include <omp.h>
#include<cstdlib>

//# define NUM_THREADS 16
#define radian M_PI/180.0
#define R 287.5
#define gamma 1.4
#define cp 1005.0
#define cv 717.5
#define	V_S_T  110.4						//	mu_ref=1.716e-5; T_ref=273.15;
#define	V_C1   1.45793265452e-06				// 	C1=(mu_ref/(pow(T_ref,1.5)))*(T_ref+S_T);
#define	T_S_T  194.4						//	K_ref=0.02414;	T_ref=273.15;
#define	T_C1   0.0025001353447				// 	C1=(K_ref/(pow(T_ref,1.5)))*(T_ref + S_T);
#define gamma_R 402.5
	
using namespace std;

//forward decleration of classes
class Face;
class Cell;
class Point;


class Vector
{
	public:
		Vector();
		~Vector();
		void operator()(const double &,const double &,const double &);
		void operator()(const double &, const double &);
                void operator()(const Point &,const Point &);
		void operator()(const Point &);
		const double& Get_Component(const int &) const;
		const double & operator()(const int& ) const ; // returns the component with the given index i
		double & magnitude();
		double dot(const Vector&) const;
		double  operator* (const Vector& v) const{return (comp1*v.comp1+comp2*v.comp2+comp3*v.comp3);};
		//dot product of two vectors
		const Vector& Cross(const Vector&)const;
		const Vector& operator+=(const Vector&);
		const Vector& operator+(const Vector&) const;
		const Vector& operator-=(const Vector&);
		const Vector& operator-(const Vector&) const;
		const Vector& operator=(const Vector&);
		const Vector& operator/=(const double&);
		const Vector& operator*= (const double&);
		const Vector& operator* (const double &) const;
		const Vector& Find_UnitVector() const;
		void Print() const;
		void Clear();
	private:
		int no_of_components;
		double comp1,comp2,comp3,mag;
		bool has_Vector_Changed;
		static Vector Vect;
};


class Point
{
        public:
                Point();
		~Point();
                Point(const double&,const double&,const double&);
                void operator() (const double&,const double&);
                void operator() (const double &,const double &,const int&);
                void operator() (const double &,const double &,const double &);
                void operator() (const double &,const double &,const double &,const int&);
		Point& operator = (const Point& );
		const Point& operator -= (const Point& );
		Point operator - (const Point& ) const;
		const Point& operator += (const Point& );
		Point operator + (const Point& ) const;
		const Point& operator /= (const double &);
		Point operator / (const double &) const;
		const Point& operator *= (const double &);
		Point operator * (const double &) const;
                friend bool operator== (const Point&,const Point&);
                friend bool operator!= (const Point&,const Point&);
                const double& Get_x() const;                 //return x
		const double& Get_y() const;                 //return y
		const double& Get_z() const;                 //return z
		const double& Get_r() const;                 //return r
		double& Get_theta();             //return theta
		double& Get_phi();               //return phi
                void Print() const;   		//writes point data
                void Convert_Polar();
                void Convert_Cartesian(); 	// for 2d
                void Convert_Cartesian(const int&); 	// for 3d
                void Convert_Spherical();
        private:
                double x,y,z,r,theta,phi;
};


class Face
{
        public:
                Face();
		~Face();
		void operator()(const Point&,const Point&,const Point&,const int &);
		void operator()(const Point&,const Point&,const Point&,const Point&,const int &);
		void Cal_Primitive(const vector<double>& );
		void Face_Center();
		void Set_Face_Pressure(const double &);
		void Set_Face_Temparature(const double &);
		void Viscosity();
		void Thermal_Conductivity();
		void Form_Centroid_Vector(Point & P1,Point & P2);
		void Cal_GradnQ(Cell *,Cell*,const int &);
		void Cal_DelnQ(Cell *,Cell*,const int &);
		const Face& operator=(const Face&);
		const int& Get_FaceNop() const;
		double Get_RmiddotA() const;
		const vector<double>& Flux_From_Face_Central(const vector<double>&, const vector<double>&);
		const vector<double>& Flux_From_Face_Central(const double &,const double &,const Vector &);
		void Flux_From_Face_VL(Cell* ,Cell * );
		void Flux_From_Face_AUSM(Cell* ,Cell * );
		void Flux_From_Face_Central(Cell* ,Cell* );
		const vector<double>& Flux_From_Face(Cell* ,Cell* ,const int &);
		const vector<double>& Add_Dissipation(Cell *,Cell *,Cell* ,Cell* );
		const vector<double>& Cal_Viscous_Stress_On_Face(const Cell* const ,const Cell* const) ;
		const Vector& Get_FaceNormal() const;
		const Vector& Get_FaceArea_Vector() const;
		const Vector& Get_Rmid() const;
		const Vector& Get_Face_undS() const ;
		const Vector& Get_Face_vndS() const ;
		const Vector& Get_Face_wndS() const ;
		const Vector& Get_Face_TndS() const ;
		const Vector& Get_Face_Velocity() const;
		const Vector & Scalar_dot_ndS(const double &);
		const double& Get_Face_Temparature() const;
		const double& Get_Face_Pressure() const;
		double Gradient_dot_ndS(const Vector &);
		const double& Get_Face_Area() const;
		void Cal_Del2Q(Cell *,Cell*);
		const vector<double>& Calculate_Grad_n_Q(const Cell * const , const Cell * const, const int &, const int & );
		const vector<double>& Get_DelQ_onFace(const int & );
		const Vector & Get_QnDS(const int &);
		double  Evaluate_Mach(const double & ,const int & );
		double  Evaluate_Pressure(const double & ,const double &,const int & );
		void Find_Gradients();
        private:
                Point Point1,Point2,Point3,Point4,Center_Point;
		double Face_P,Face_T,Face_Rho,Face_ET,Face_a,Area_Comp1,Area_Comp2,Area_Comp3,max_eigen_value,Face_area,mu,K,v1,v2,v3,Face_M,Heat_Flux;
		Vector Area_Vector,Face_Normal,R_mid,Face_Velocity,undS,vndS,wndS,TndS,rho_ndS,rhou_ndS,rhov_ndS,rhow_ndS,rhoET_ndS;
		int face_nop,Face_No;
		Vector Avg_U_GradonFace,Avg_V_GradonFace,Avg_W_GradonFace, Avg_T_GradonFace,Centroid_Vector,Centroid_Normal;
		vector<double> Flux,Del_Q_1,Del_Q_2,Del_Q_3;
};

class Cell
{
        public:
                Cell();
               ~Cell();
	        void operator()(const Point&,const Point&,const Point&,const Point&,const Point&,const Point&,const Point&,const Point&);
	        void operator()(const Point&,const Point&,const Point&,const Point&,const Point&,const Point&);
		void operator()(const vector<double>& );// assigns Q for cell
		void operator()(const int&,const int &,const int &,const int &,const int &,const int &,const int &);
		void Cell_Center(const Point&,const Point&,const Point&,const Point&,const Point&,const Point&,const Point&,const Point&);
		void Volume();
		void Evaluate_FluxonFaces(const int &);
		void Q_atCellCenter(const double&,const double&, const Vector &);
		void Set_NeighbourCells(const vector<Cell*>* const);
		void Print_Neighbours() const;
		void UpdateQ();
		void Cal_Primitive(const vector<double>& );
		void Set_Cell_Velocity(const Vector&);
		void Check_Cell() const;
// 		void Add_Dissipation();
		void Write_Cell_Info(const string &) const;
		void Test_Fluxes(const int &);
		void Add_Fluxes(const vector<double>&, const int & );
		void Cal_Stresses();
		void Set_FluxofCell(vector<double>& );
		void Velocity_GradatCenter();
		void Cal_Q_Gradient();
		void Cal_Del2_Q();
		void Cal_Del3_Q();
		void Cal_Del4_Q();
		void Cal_Del5_Q();
		void Cal_Del6_Q();
		const Point& Get_Cell_Center() const;
		const Face& Get_Face(const int&) const;
		const int& Get_Neighbours(const int&) const;// return neighbours based on face		
		const int& no_of_Faces() const;
		const double& Get_Volume() const;
		const double& Get_Celldx() const;
		const double& Cell_Face_Area(const int&) const;
		const double& Get_Cell_Temparature() const;
		const double& Get_Cell_Pressure() const;
		const double& Get_Cell_Density() const;
		vector<double>& Get_QatCellCenter();
		const vector<double>& Get_DelnQ(const int & i, const int & j);
		const vector<double>& Get_DelQ_fromCell(const int &,const int &) const;
		vector<double>& Rk4(const double&,const int &);
		const vector<double>& Get_FluxofCell() const;
		const Vector& Get_Cell_Velocity() const;
		const Vector& Get_Diagonal_Vector() const;
		const Vector& Get_Grad_UatCenter() const;
		const Vector& Get_Grad_VatCenter() const;
		const Vector& Get_Grad_WatCenter() const;
		const Vector& Get_Grad_TatCenter() const;
		const Vector & Get_Q_Grad(const int & );
		const bool& Is_Ghost_Cell() const;
		Cell* Get_Cell(const int &);
        private:
			 Face Top_Face,Bottom_Face,Left_Face,Right_Face,Front_Face,Back_Face;
			Point Cell_MidPoint;
			Cell *Front_Cell,*Back_Cell,*Top_Cell,*Bottom_Cell,*Left_Cell,*Right_Cell;
			int cell_nop,No_of_Faces, self_index,index1,index2,index3,index4,index5,index6,scheme_type;
			vector<double> Cell_Q,Cell_Flux,Q_new;
			vector<double> Del2_Q_1,Del2_Q_2,Del2_Q_3,Del3_Q_1,Del3_Q_2,Del3_Q_3,Del4_Q_1,Del4_Q_2,Del4_Q_3,Del5_Q_1,Del5_Q_2,Del5_Q_3,Del6_Q_1,Del6_Q_2,Del6_Q_3,Viscous_Flux;
			double Cell_Pressure,Cell_Temparature,Cell_Rho,Cell_Et,v1,v2,v3,M,a,Average_total_area,inv_vol,Cell_Avg_Length,volume,max_eigen_value;
			double FF_Area,BaF_Area,TF_Area,LF_Area,BoF_Area,RF_Area,lambda;
			Vector Cell_Velocity,u_GradatCenter,v_GradatCenter,w_GradatCenter,T_GradatCenter,Diagonal_Vector;
			Vector rho_Gradient,rhou_Gradient,rhov_Gradient,rhow_Gradient,rhoET_Gradient;
			bool is_Ghost_Cell;
};


class Cell_Property_Manager
{
	public:
		Cell_Property_Manager();
		~Cell_Property_Manager();
		void Cal_Dissipation();
		void Cal_QAverageonFaces();
		void Cal_QAverageonFaces(const int &);
		void Cal_Primitive(vector<double>&);
		void Cal_Computational_Variables(const double&,const double&,const Vector&);
		void Cal_Computational_Variables(const double&,const double&,const double&,const double&,const double&,const double&);
		void Initialize();
		void Initialize(const string&);
		void Apply_Boundary_Conditions();
		void Wall_Bc();
		void Left_WallBc();
		void Right_WallBc();
		void Top_WallBc();
		void Bottom_WallBc();
		void Inlet_Bc();
		void Exit_Bc();
		void Cal_StressOnFace();
		void Find_Relative_Change();
		void Solve();
		void Cal_QatCellCenter(const int&);
		void write_var(const int&);
		void write_residue(const int&,const double&,const double&,const double&,const double&,const double&) const;
		void Form_Cells();
		void Read_Grid(const string&);
		void Rk4(const double&,const int & );
		void Flux_Richardson_Extrapolation(const int &);
		void Calculate_Dissipation();
	private:
		vector<double> Q,Q1,Q_new,temp_Q;
		Vector V;
		double Rho,P_Total_in,T_Total_in,P_Static_in,P_Static_out,v1,v2,v3,M,Et,a,Pressure,dist,Temparature;
		double total_error,error,delx,mindelx,CFL,time_step,cell_volume,max_volume;
		double rho_error,rhou_error,rhov_error,rhow_error,rhoet_error;
		string filename,filename_3h;
		vector<Cell*> Cell_List,Cell_List_h,Cell_List_3h;
		int nx,ny,nz,No_of_Cells,No_of_Ghost_Cells,Cells_in_Plane,Total_Cells,no_of_timesteps,scheme_type;
};

#endif  //#ifndef _Header_H
