#include "definitions.h"
#include "Globals.h"
#include "Boundary_Conditions.h"
#include "Utilities.h"
#include "Primitive_Computational.h"

void Far_Field_Boundary_Condition(int &Case_Type)
{
}

void Super_Sonic_Inflow(int &Cell_No, int &Face_No, int &Ghost_Cell_Index)
{
  //   cout<<"Incoming Flow bc\n"<<Cell_No<<"\t"<<Face_No<<"\t"<<Ghost_Cell_Index<<endl;
  double Rho = 0.0, P = 0.0, u = 0.0, v = 0.0, vmag = 0.0, T = 0.0;

  v = -161.67;
  u = 915.0;
  P = 213947.0;
  T = 377.0;
  Rho = P / (R_GC * T);
  vmag = 0.5 * (u * u + v * v);
  U_Cells[Ghost_Cell_Index][0] = Rho;
  U_Cells[Ghost_Cell_Index][1] = Rho * u;
  U_Cells[Ghost_Cell_Index][2] = Rho * v;
  U_Cells[Ghost_Cell_Index][3] = Rho * (cv * T + vmag);
  Calculate_Primitive_Variables(Cell_No, U_Cells[Ghost_Cell_Index]);
  Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
  for (unsigned int i = 0; i < Global_Primitive.size(); i++)
    Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
  // 		Print(Primitive_Cells[Ghost_Cell_Index]);
}

void SuperSonic_outflow()
{
}

void SubSonic_Inflow()
{
}

void SubSonic_Outflow()
{
}