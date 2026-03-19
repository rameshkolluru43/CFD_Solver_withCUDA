#include "definitions.h"
#include "Globals.h"
#include "Viscous_Functions.h"
#include "Primitive_Computational.h"

// string pwd = std::filesystem::current_path();
// pwd = pwd  + "/QPDE_FVM_CFD_Solver";

vector<V_D> U_Cells, Cells_Net_Flux, Cells_DelU, A_x, A_y, A_x_L, A_x_R, A_y_T, A_y_B, A;
vector<V_D> U_Cells_RK_1, U_Cells_RK_2;

vector<bool> v_bool;
vector<vector<bool>> Cells_Face_Boundary_Type;
double Pressure_Total_Inlet, Temperature_Total_Inlet, Pressure_Static_Exit, Re, Pr, Inv_Re, Inv_Pr;
double Max_dt, Min_dt, u_ref, t_ref, L_ref, M_ref, mu_ref, P_ref, Rho_ref, K1, mu_star, R_ref, cp_ref;
double Pressure_Static_Inlet, Rho_Static_Inlet, Temperature_Static_Inlet, Inlet_Mach_No, V_1, V_2, V_Magnitude, V_Magnitude_Inlet, C_Acoustic, M;
int nx_1, nx_2, ny_1, ny_2;
int Dissipation_Type, Flux_Type, Test_Case, Initialize_Type, iterations, Grid_Size, Area_Weighted_Average, Total_Iterations, Limiter_Case;
bool Is_Second_Order, Is_MOVERS_1, Enable_Entropy_Fix, Is_Time_Dependent, Time_Accurate, Local_Time_Stepping, Non_Dimensional_Form, Is_WENO, Is_Char, Is_Implicit_Method;
bool Enable_AMR = false;
int AMR_Period = 100;
double AMR_Gradient_Threshold = 0.1, AMR_Max_Fraction = 0.3;
vector<double> Gradient_Refinement_Indicator;
string Grid_File, Initial_Solution_File, Solution_File, Error_File, Limiter_File, Final_Solution_File, Grid_Vtk_File, CF_File;
// variables to read the face normals and face lenghts and face area
double nx, ny, dl, dA;

double Mod_Alpha0, Mod_Alpha1, Mod_Alpha2, Mod_Alpha3, d_F_0, d_F_1, d_F_2, d_F_3, d_U_0, d_U_1, d_U_2, d_U_3;
double Lambda_Max, Lambda_Min, Max1, Max2, Min1, Min2, d_Var;
vector<V_D> Cells_Viscous_Flux, Cells_T_Grad, Cells_u_Grad, Cells_v_Grad, Cells_Rho_Grad;
double mu = 0.0, K = 0.0;
int Viscous_Time_Case;
// Left State Variables
double Rho_L, T_L, P_L, u_L, v_L, Vdotn_L, Vmag_L, M_L, C_L, H_L;
V_D Flux_L, U_L, CF, b;
int NUM_FLUX_COMPONENTS;
// Right State Varialbes
double Rho_R, T_R, P_R, u_R, v_R, Vdotn_R, Vmag_R, M_R, C_R, H_R, Qx, Qy;
V_D Flux_R, U_R;

int Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0, Neighbour_4 = 0, Neighbour_5 = 0, Neighbour_6 = 0, Neighbour_7 = 0, Neighbour_8 = 0;
double mev_L, mev_R, max_eigen_value, Terminating_Time, Total_Time;
double Roe_u, Roe_v, Roe_Rho, Roe_P, Roe_e, Roe_Vmag, U_avg, Roe_a, Roe_H, Un, Ut, cs;
double P_inf, Rho_inf, R_inf, M_inf, u_inf, T_inf, q_inf; // Free Stream Varialbles
V_D d_U(4, 0.0), d_F(4, 0.0), Mod_Alpha(4, 0.0);
V_D Global_Primitive, Global_U, Average_Convective_Flux, Dissipative_Flux;
// Global_Primitive represents global Vector for storing Primitive Variables
V_D Cells_DelT;
vector<V_D> Primitive_Cells;

V_D u_Gradient, v_Gradient, T_Gradient, Rho_Gradient, P_Gradient, Viscous_Flux, Grad_Var, S;

vector<int> row_indices;
vector<int> col_indices;
vector<double> values;

// Function used for initialization of computational vectors from the exit conditions
void Initialize(const int &Test_Case)
{
	cout << "Allocating Memory for Physical Cells\t" << No_Physical_Cells << endl;
	cout << "Total Number of Cells including Ghost Cells\t" << Total_No_Cells << endl;
	Average_Convective_Flux.resize(NUM_FLUX_COMPONENTS, 0.0);
	Dissipative_Flux.resize(NUM_FLUX_COMPONENTS, 0.0);
	Flux_L.resize(NUM_FLUX_COMPONENTS, 0.0);
	Flux_R.resize(NUM_FLUX_COMPONENTS, 0.0);
	Global_U.resize(NUM_FLUX_COMPONENTS, 0.0);
	U_L.resize(NUM_FLUX_COMPONENTS, 0.0);
	U_R.resize(NUM_FLUX_COMPONENTS, 0.0);
	v_bool.resize(NUM_FLUX_COMPONENTS, false);
	Global_Primitive.resize(11, 0.0);
	u_Gradient.resize(2, 0.0);
	v_Gradient.resize(2, 0.0);
	T_Gradient.resize(2, 0.0);
	Viscous_Flux.resize(NUM_FLUX_COMPONENTS, 0.0);
	Grad_Var.resize(2, 0.0);
	V_I Second_Neighbours(4, 0.0);
	S.resize(6, 0.0);
	//	cout<<"am here"<<endl;
	for (unsigned int i = 0; i < Wall_Cells_List.size(); i += 3)
	{
		CF.push_back(0.0);
	}
	cout << CF.size() << endl;
	for (int index = 0; index < No_Physical_Cells; index++)
		Cells_DelT.push_back(0.0);
	for (int index = 0; index < Total_No_Cells; index++)
	{
		Cells_Net_Flux.push_back(Average_Convective_Flux);
		Cells_DelU.push_back(Global_U);
		int nf = 4;
		if (index < No_Physical_Cells && Cells[index].numFaces > 0)
			nf = Cells[index].numFaces;
		Cells_Face_Boundary_Type.push_back(vector<bool>(nf, false));
		U_Cells.push_back(Global_U);
		U_Cells_RK_1.push_back(Global_U);
		U_Cells_RK_2.push_back(Global_U);
		Primitive_Cells.push_back(Global_Primitive);
		Cells_Viscous_Flux.push_back(Viscous_Flux);
	}
	// Identify_Wall_Boundary_Faces(Grid_Type);
	//		cout<<"am here"<<endl;
	if (Is_Viscous_Wall)
	{
		for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
			Identify_Neighbours_For_Second_Gradients(Cell_Index);
	}
	cout << "Memory Allocation Done for all Cells\n";

	if (Is_Implicit_Method)
	{
		cout << "Memory for creating Jacobian and A and b matrices started " << endl;
		A_x.resize(4, V_D(4, 0.0));
		A_y.resize(4, V_D(4, 0.0));
		A_x_L.resize(4, V_D(4, 0.0));
		A_y_B.resize(4, V_D(4, 0.0));
		A_x_R.resize(4, V_D(4, 0.0));
		A_y_T.resize(4, V_D(4, 0.0));

		A.resize(4 * Total_No_Cells, V_D(4 * Total_No_Cells, 0.0));
		b.resize(4 * Total_No_Cells, 0.0);
		cout << "Memory for creating Jacobian and A and b matrices ended " << endl;
	}
}

// Function used for initialization from a given file at a given time step
void Initialize(const string &file_name)
{
	ifstream ipfile(file_name.c_str(), ios::in);
	int n1, n2, n3, N_Cells, CellIndex;
	V_D V(2, 0.0);
	double P = 0.0, Rho = 0.0, T = 0.0, v1 = 0.0, v2 = 0.0, dt = 0.0, P0 = 0.0;
	Initialize(Test_Case);
	cout << file_name << endl;

	if (ipfile.is_open())
	{
		cout << "File opened for initialization\t" << file_name << endl;
		ipfile >> n1 >> n2 >> n3;
		cout << n1 << "\t" << n2 << "\t" << n3 << endl;
		ipfile >> N_Cells;
		cout << "Total Number of Cells\t" << N_Cells << endl;
		ipfile >> iterations;
		cout << "Solution Initialized from Iteration Number \t" << iterations << endl;
		if (N_Cells == No_Physical_Cells)
		{
			cout << "Number of Physical Cells are equal\n";
			for (int index = 0; index < N_Cells; index++)
			{
				ipfile >> CellIndex >> dt >> Rho >> P >> T >> v1 >> v2 >> M >> P0;
				//  				cout<<index<<"\t"<<Rho<<"\t"<<T<<"\t"<<P<<"\t"<<v1<<"\t"<<v2<<endl;
				V[0] = v1;
				V[1] = v2;
				Calculate_Computational_Variables(P, V, Rho, 2);
				for (unsigned int j = 0; j < Global_U.size(); j++)
					U_Cells[index][j] = Global_U[j];
				Calculate_Primitive_Variables(index, U_Cells[index]);
				for (unsigned int j = 0; j < Global_Primitive.size(); j++)
					Primitive_Cells[index][j] = Global_Primitive[j];
			}
		}
		else
		{
			cout << "Mismatch in Solution file.........Please verify that correct solution file is given for input\n";
		}
		cout << "Initialization done from file\t" << file_name << endl;
	}
	else
	{
		cout << "Could not Open Inputfile for reading......... Please check  the file name\n";
		cout << file_name << endl;
		exit(0);
	}
	//	Print(U_Cells);
}

// This Function identifies the Wall boundary faces give a grid for channel flow for both Rectangular and Circular cross sections
// Wall_Cells_List contains the list of Wall Cells with Cell Index, Face Number and Ghost Cell Index
// Cells_Face_Boundary_Type is a 2D vector which stores the boundary type of the face of the cell
// For all the cells in the Wall_Cells_List the face boundary type is set to true, This function is used in wall boundary condition implementation

void Identify_Wall_Boundary_Faces(const int &Grid_Type)
{
	cout << "Identify_Boundary_Types \t" << Grid_Type << endl;
	for (unsigned int i = 0; i < Wall_Cells_List.size(); i += 3)
	{
		int cell = Wall_Cells_List[i];
		int face = Wall_Cells_List[i + 1];
		if (cell < static_cast<int>(Cells_Face_Boundary_Type.size()) && face < static_cast<int>(Cells_Face_Boundary_Type[cell].size()))
			Cells_Face_Boundary_Type[cell][face] = true;
	}
	for (unsigned int i = 0; i < Inlet_Cells_List.size(); i += 3)
	{
		int cell = Inlet_Cells_List[i];
		int face = Inlet_Cells_List[i + 1];
		if (cell < static_cast<int>(Cells_Face_Boundary_Type.size()) && face < static_cast<int>(Cells_Face_Boundary_Type[cell].size()))
			Cells_Face_Boundary_Type[cell][face] = true;
	}
	for (unsigned int i = 0; i < Exit_Cells_List.size(); i += 3)
	{
		int cell = Exit_Cells_List[i];
		int face = Exit_Cells_List[i + 1];
		if (cell < static_cast<int>(Cells_Face_Boundary_Type.size()) && face < static_cast<int>(Cells_Face_Boundary_Type[cell].size()))
			Cells_Face_Boundary_Type[cell][face] = true;
	}
	for (unsigned int i = 0; i < Symmetry_Cells_List.size(); i += 3)
	{
		int cell = Symmetry_Cells_List[i];
		int face = Symmetry_Cells_List[i + 1];
		if (cell < static_cast<int>(Cells_Face_Boundary_Type.size()) && face < static_cast<int>(Cells_Face_Boundary_Type[cell].size()))
			Cells_Face_Boundary_Type[cell][face] = true;
	}
}
