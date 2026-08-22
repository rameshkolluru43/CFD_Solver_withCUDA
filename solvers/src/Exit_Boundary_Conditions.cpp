#include "definitions.h"
#include "Globals.h"
#include "Boundary_Conditions.h"
#include "Utilities.h"
#include "Primitive_Computational.h"

// Case_Type Indicates whether applicable for Supersonic case or subsnoic case
void Subsonic_Exit_Boundary_Condition(unsigned int &i)
{
	int Cell_Index = 0, Ghost_Cell_Index = 0, index, Face_No;
	double v1 = 0.0, v2 = 0.0, Rho = 0.0, n1 = 0.0, n2 = 0.0;
	double rho_f = 0.0, u_f = 0.0, v_f = 0.0, v_mag_f = 0.0, P_f = 0.0, P = 0.0, temp = 0.0, rho_Et;
	double Rho_o = 1.0, C_o = 1.0; // Representing Free Stream Conditions
	V_D U_Face(5, 0.0);
	// values from the interior cell
	Cell_Index = Exit_Cells_List[i + 0];
	Ghost_Cell_Index = Exit_Cells_List[i + 2];

	// 		cout<<"\ncell indexat exit\t"<<Cell_Index<<endl;
	// 		cout<<"Ghost_Cell_Index exit\t"<<Ghost_Cell_Index<<endl;
	Rho = Primitive_Cells[Cell_Index][0];
	v1 = Primitive_Cells[Cell_Index][1];
	v2 = Primitive_Cells[Cell_Index][2];
	P = Primitive_Cells[Cell_Index][4];
	rho_Et = Primitive_Cells[Cell_Index][6];

	// 		cout<<"\nRho v1 v2 P C rho_Et at exit\n"<<Rho<<"\t"<<v1<<"\t"<<v2<<"\t"<<P<<"\t"<<C<<"\t"<<rho_Et<<endl;
	Face_No = Exit_Cells_List[i + 1];
	// 		cout<<"Face no on exit cell\t"<<Face_No<<"\n"<<endl;
	index = Face_No * 2;
	// 		cout<<"index\t"<<index<<endl;

	/*n1=Cell_Face_Normals[Cell_Index][index+0];
	n2=Cell_Face_Normals[Cell_Index][index+1];*/
	n1 = Cells[Cell_Index].Face_Normals[index + 0];
	n2 = Cells[Cell_Index].Face_Normals[index + 1];
	// 		cout<<"normals\t"<<n1<<"\t"<<n2<<"\t"<<endl;

	// Pressure Value to be prescribed at the exit boundary
	P_f = 0.9 / 1.4; // Pressure_Static_Exit;
	//		cout<<Pressure_Static_Exit<<endl;
	// 		cout<<Rho_o<<"\t"<<C_o<<endl;
	temp = (P - P_f) / (Rho_o * C_o * C_o);
	rho_f = Rho - temp * Rho_o;
	u_f = v1 + n1 * temp * C_o;
	v_f = v2 + n2 * temp * C_o;
	v_mag_f = 0.5 * (u_f * u_f + v_f * v_f);

	// 		cout<<"Values at exit boundary\t"<<rho_f<<"\t"<<u_f<<"\t"<<v_f<<"\t"<<T_f<<"\t"<<P_f<<endl;
	U_Face[0] = rho_f;
	U_Face[1] = rho_f * u_f;
	U_Face[2] = rho_f * v_f;
	U_Face[3] = ((P_f / (gamma - 1.0)) + rho_f * v_mag_f);

	U_Cells[Ghost_Cell_Index][0] = 2.0 * U_Face[0] - Rho;
	U_Cells[Ghost_Cell_Index][1] = 2.0 * U_Face[1] - Rho * v1;
	U_Cells[Ghost_Cell_Index][2] = 2.0 * U_Face[2] - Rho * v2;
	U_Cells[Ghost_Cell_Index][3] = 2.0 * U_Face[3] - rho_Et;

	Calculate_Primitive_Variables(Cell_Index, U_Cells[Ghost_Cell_Index]);
	Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);

	for (unsigned int i = 0; i < Global_Primitive.size(); i++)
		Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
	// 	cout<<"*******************Exit Boundary Conditions***********************************\n";
}

void Supersonic_Exit_Boundary_Condition(unsigned int &i)
{
	//     cout<<"Applying Supersonic Exit Boundary Condition"<<endl;
	int Cell_Index = 0, Ghost_Cell_Index = 0;
	Cell_Index = Exit_Cells_List[i + 0];
	Ghost_Cell_Index = Exit_Cells_List[i + 2];
	for (unsigned int i = 0; i < U_Cells[Cell_Index].size(); i++)
		U_Cells[Ghost_Cell_Index][i] = U_Cells[Cell_Index][i];
	Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);
	Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
	for (unsigned int i = 0; i < Global_Primitive.size(); i++)
		Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
	// 	cout<<"*******************Exit Boundary Conditions***********************************\n";
}

void Subsonic_Exit_Condition(ExitCondition &exitCond, V_I &Exit_Cells_List)
{
}

void Supersonic_Exit_Condition(ExitCondition &exitCond, V_I &Exit_Cells_List)
{
	int Cell_Index = 0, Ghost_Cell_Index = 0;
	for (unsigned int i = 0; i < Exit_Cells_List.size(); i += 3)
	{
		Supersonic_Exit_Boundary_Condition(i);
	}
}