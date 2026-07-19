#include "definitions.h"
#include "Globals.h"
#include "Boundary_Conditions.h"
#include "Primitive_Computational.h"
#include "Utilities.h"
#include <cmath>

// Inviscid Wall Boundary Condition
void Symmetry_Boundary_Condition()
{
	double v1_Ghost_Cell = 0.0, v2_Ghost_Cell = 0.0, v1_Interior_Cell = 0.0, v2_Interior_Cell = 0.0, n1 = 0.0, n2 = 0.0, vmag = 0.0, vdotn = 0.0;
	int Face_Index = 0, Cell_Index, Ghost_Cell_Index, Face_No = 0;
	//  	cout<<"Size of Wall_Cells_List\t"<<Wall_Cells_List.size()<<endl;

	for (unsigned int i = 0; i < Symmetry_Cells_List.size(); i += 3)
	{

		Cell_Index = Symmetry_Cells_List[i + 0];
		Face_No = Symmetry_Cells_List[i + 1];
		Ghost_Cell_Index = Symmetry_Cells_List[i + 2];

		v1_Interior_Cell = Primitive_Cells[Cell_Index][1];
		v2_Interior_Cell = Primitive_Cells[Cell_Index][2];

		Face_Index = Face_No * 2; // Face_Index = Face_No*2

		/*n1 = Cell_Face_Normals[Cell_Index][Face_Index+0];		//Face_No*2+0
		n2 = Cell_Face_Normals[Cell_Index][Face_Index+1];		//Face_No*2+1*/
		n1 = Cells[Cell_Index].Face_Normals[Face_Index + 0];
		n2 = Cells[Cell_Index].Face_Normals[Face_Index + 1];

		//   		cout<<Cell_Index<<"\t"<<Ghost_Cell_Index<<"\t"<<Face_No<<"\t"<<n1<<"\t"<<n2<<endl;

		vdotn = (v1_Interior_Cell * n1 + v2_Interior_Cell * n2); // V.n
																 //  		cout<<"Vdotn\t"<<vdotn<<endl;
		v1_Ghost_Cell = v1_Interior_Cell - 2.0 * vdotn * n1;
		v2_Ghost_Cell = v2_Interior_Cell - 2.0 * vdotn * n2; // V = V - 2*(V.n) component wise

		vmag = 0.5 * (v1_Ghost_Cell * v1_Ghost_Cell + v2_Ghost_Cell * v2_Ghost_Cell);
		//                 cout<<v1_Interior_Cell<<"\t"<<v2_Interior_Cell<<"\t"<<v1_Ghost_Cell<<"\t"<<v2_Ghost_Cell<<endl;

		U_Cells[Ghost_Cell_Index][0] = Primitive_Cells[Cell_Index][0];
		U_Cells[Ghost_Cell_Index][1] = Primitive_Cells[Cell_Index][0] * v1_Ghost_Cell;
		U_Cells[Ghost_Cell_Index][2] = Primitive_Cells[Cell_Index][0] * v2_Ghost_Cell;
		U_Cells[Ghost_Cell_Index][3] = ((Primitive_Cells[Cell_Index][4] / (gamma - 1.0)) + Primitive_Cells[Cell_Index][0] * vmag);

		Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);

		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		for (unsigned int i = 0; i < Global_Primitive.size(); i++)
			Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
	}
}

// Ghost cells for wall starts after No_Physical_Cells + Front plane gc Cells_in_Plane
void Viscous_Wall_Boundary_Condition()
{
	int Cell_Index = 0, Ghost_Cell_Index = 0, Face_No = 0;
	for (unsigned int i = 0; i < Wall_Cells_List.size(); i += 3)
	{

		Cell_Index = Wall_Cells_List[i + 0];
		Face_No = Wall_Cells_List[i + 1];
		Ghost_Cell_Index = Wall_Cells_List[i + 2];

		/* No-slip: reflect momentum so face velocity is zero. */
		const double rho_i = Primitive_Cells[Cell_Index][0];
		const double u_i = Primitive_Cells[Cell_Index][1];
		const double v_i = Primitive_Cells[Cell_Index][2];
		const double T_i = Primitive_Cells[Cell_Index][3];
		const double p_i = Primitive_Cells[Cell_Index][4];

		double rho_g = rho_i;
		double p_g = p_i;
		if (wallCond.T > 0.0)
		{
			/* Isothermal wall: Dirichlet T at face via mirrored ghost temperature. */
			double T_g = 2.0 * wallCond.T - T_i;
			if (!(T_g > 1.0e-12) || !std::isfinite(T_g))
				T_g = wallCond.T;
			p_g = p_i;
			const double Rgas = Non_Dimensional_Form ? R_ref : R_GC;
			rho_g = p_g / (Rgas * T_g);
			if (!(rho_g > 1.0e-14) || !std::isfinite(rho_g))
				rho_g = rho_i;
		}

		U_Cells[Ghost_Cell_Index][0] = rho_g;
		U_Cells[Ghost_Cell_Index][1] = -rho_g * u_i;
		U_Cells[Ghost_Cell_Index][2] = -rho_g * v_i;
		U_Cells[Ghost_Cell_Index][3] = p_g / (gamma - 1.0) + 0.5 * rho_g * (u_i * u_i + v_i * v_i);

		Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);
		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		for (unsigned int k = 0; k < Global_Primitive.size(); k++)
			Primitive_Cells[Ghost_Cell_Index][k] = Global_Primitive[k];
	}
}

// Inviscid Wall Boundary Condition
void InViscid_Wall_Boundary_Condition()
{
	double v1_Ghost_Cell = 0.0, v2_Ghost_Cell = 0.0, v1_Interior_Cell = 0.0, v2_Interior_Cell = 0.0, n1 = 0.0, n2 = 0.0, vmag = 0.0, vdotn = 0.0;
	int Face_Index = 0, Cell_Index, Ghost_Cell_Index, Face_No = 0;
	//  	cout<<"Size of Wall_Cells_List\t"<<Wall_Cells_List.size()<<endl;

	for (unsigned int i = 0; i < Wall_Cells_List.size(); i += 3)
	{

		Cell_Index = Wall_Cells_List[i + 0];
		Face_No = Wall_Cells_List[i + 1];
		Ghost_Cell_Index = Wall_Cells_List[i + 2];

		v1_Interior_Cell = Primitive_Cells[Cell_Index][1];
		v2_Interior_Cell = Primitive_Cells[Cell_Index][2];

		Face_Index = Face_No * 2; // Face_Index = Face_No*2

		/*n1 = Cell_Face_Normals[Cell_Index][Face_Index+0];		//Face_No*2+0
		n2 = Cell_Face_Normals[Cell_Index][Face_Index+1];		//Face_No*2+1*/
		n1 = Cells[Cell_Index].Face_Normals[Face_Index + 0];
		n2 = Cells[Cell_Index].Face_Normals[Face_Index + 1];

		//   		cout<<Cell_Index<<"\t"<<Ghost_Cell_Index<<"\t"<<Face_No<<"\t"<<n1<<"\t"<<n2<<endl;

		vdotn = (v1_Interior_Cell * n1 + v2_Interior_Cell * n2); // V.n
																 //  		cout<<"Vdotn\t"<<vdotn<<endl;
		v1_Ghost_Cell = v1_Interior_Cell - 2.0 * vdotn * n1;
		v2_Ghost_Cell = v2_Interior_Cell - 2.0 * vdotn * n2; // V = V - 2*(V.n) component wise

		vmag = 0.5 * (v1_Ghost_Cell * v1_Ghost_Cell + v2_Ghost_Cell * v2_Ghost_Cell);

		//                 cout<<v1_Interior_Cell<<"\t"<<v2_Interior_Cell<<"\t"<<v1_Ghost_Cell<<"\t"<<v2_Ghost_Cell<<endl;

		U_Cells[Ghost_Cell_Index][0] = Primitive_Cells[Cell_Index][0];
		U_Cells[Ghost_Cell_Index][1] = Primitive_Cells[Cell_Index][0] * v1_Ghost_Cell;
		U_Cells[Ghost_Cell_Index][2] = Primitive_Cells[Cell_Index][0] * v2_Ghost_Cell;
		U_Cells[Ghost_Cell_Index][3] = ((Primitive_Cells[Cell_Index][4] / gamma_M_1) + Primitive_Cells[Cell_Index][0] * vmag);

		Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);

		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		for (unsigned int i = 0; i < Global_Primitive.size(); i++)
			Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
	}
}
