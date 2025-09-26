#include "definitions.h"
#include "Globals.h"
#include "Viscous_Functions.h"

// Function to calculate viscosity using Sutherland's law
// T_Star: Average temperature on the face (non-dimensional)
void Viscosity(double &T_Star)
{
    double Term1 = 0, Term2 = 0;
    Term1 = T_S_Mu / T_ref;
    Term2 = T_Star + Term1;
    mu_star = pow(T_Star, 1.5) * ((1.0 + Term1) / Term2); // Non-dimensional viscosity
}

// Function to calculate thermal conductivity using Sutherland's law
// T: Temperature (non-dimensional)
void Thermal_Conductivity(double &T)
{
    K = mu_star * cp_ref / Pr; // Non-dimensional thermal conductivity
}

// Function to set reference values for non-dimensionalization
void Reference_Values()
{
    mu_ref = L_ref / Re;
    cp_ref = 1.0 / (gamma_M_1 * M_ref * M_ref);
    R_ref = 1.0 / (gamma * M_ref * M_ref);
    Rho_inf = 1.0;
    T_inf = 1.0;
    P_inf = Rho_inf * R_ref * T_inf;
    q_inf = 1.0;
}

void Evaluate_Wall_Skin_Friction()
{
    int Cell_Index = 0, Ghost_Cell_Index = 0, Face_No = 0, Grad_Type, index, ListSize, StartGridPoint;
    double u11, u12, u21, u22, dudn, dvdn, dVdn, dVxdn, dVydn, sgn, dVdna, vmag;
    double v1 = 0.0, v2 = 0.0, mu = 0.0, tx = 0.0, ty = 0.0, tw = 0.0, T11 = 0.0, T12 = 0.0, T21 = 0.0, T22 = 0.0;
    ListSize = Wall_Cells_List.size();
    StartGridPoint = Wall_Cells_List[0];
    //	cout<<"Wall List size"<<ListSize<<"\t"<<StartGridPoint<<endl;

    for (unsigned int i = 0; i < ListSize; i += 3)
    {

        Cell_Index = Wall_Cells_List[i + 0];
        Face_No = Wall_Cells_List[i + 1];
        Ghost_Cell_Index = Wall_Cells_List[i + 2];
        index = Face_No * 2;

        /*Methodology to find Skin Friction on the wall*/
        /*nx = Cell_Face_Normals[Cell_Index][index+0];
        ny = Cell_Face_Normals[Cell_Index][index+1];
        dl = Cell_Face_Areas[Cell_Index][Face_No];*/
        nx = Cells[Cell_Index].Face_Normals[index + 0];
        ny = Cells[Cell_Index].Face_Normals[index + 1];
        dl = Cells[Cell_Index].Face_Areas[Face_No];

        tx = ny;
        ty = -nx;

        v1 = 0.5 * (Primitive_Cells[Cell_Index][1] + Primitive_Cells[Ghost_Cell_Index][1]);
        v2 = 0.5 * (Primitive_Cells[Cell_Index][2] + Primitive_Cells[Ghost_Cell_Index][2]);
        mu = 0.5 * (Primitive_Cells[Cell_Index][8] + Primitive_Cells[Ghost_Cell_Index][8]);

        //			cout<<nx<<"\t"<<ny<<"\t"<<dl<<endl;

        //			cout<<v1<<"\t"<<v2<<"\t"<<mu<<endl;

        Grad_Type = 1; // U velocity Gradient		// Finds the u Velocity gradient on the face required
        Calculate_Gradient_On_Face(Cell_Index, Grad_Type, Face_No);
        u11 = u_Gradient[0];
        u12 = u_Gradient[1];
        Grad_Type = 2; // V Velocity Gradient
        Calculate_Gradient_On_Face(Cell_Index, Grad_Type, Face_No);
        u21 = v_Gradient[0];
        u22 = v_Gradient[1];

        T11 = (2.0 / 3.0) * mu * Inv_Re * (2.0 * u11 - u22);
        T12 = mu * Inv_Re * (u12 + u21);
        T21 = T12;
        T22 = (2.0 / 3.0) * mu * Inv_Re * (2.0 * u22 - u11);

        Pressure_Static_Inlet = P_ref;
        Rho_Static_Inlet = Rho_ref;
        int u, v;
        u = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
        v = 0.0;

        vmag = u * u + v * v;

        tw = (T11 * nx + T12 * ny) * tx + (T21 * nx + T22 * ny) * ty;
        // cout<<"Value of tw\t"<<tw<<endl;
        CF[Cell_Index - StartGridPoint] = 2.0 * tw / (Rho_inf * q_inf * q_inf);
        //	cout<<"Skin friction\t"<<Cell_Index-StartGridPoint<<"\t"<<CF[Cell_Index-StartGridPoint]<<endl;
    }
}