#include "definitions.h"
#include "Globals.h"
#include "Timestep.h"

void Viscous_Time_Step_1(int &Cell_No)
{
    double nx_0 = 0.0, ny_0 = 0.0, dl_0 = 0.0, nx_1 = 0.0, ny_1 = 0.0, dl_1 = 0.0, nx_2 = 0.0, ny_2 = 0.0, dl_2 = 0.0, nx_3 = 0.0, ny_3 = 0.0, dl_3 = 0.0;
    double u0 = 0.0, v0 = 0.0, C0 = 0.0, Avg_nx_i = 0.0, Avg_ny_i = 0.0, Avg_nx_j = 0.0, Avg_ny_j = 0.0, Avg_dl_j = 0.0, Avg_dl_i = 0.0;
    double Inv_Area = 0.0, C11 = 0.0;
    int Face_No, index;
    Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0, Neighbour_4 = 0;
    double Rho0 = 0.0, Rho1 = 0.0, Rho2 = 0.0, Rho3 = 0.0, Rho4 = 0.0, u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0, v4 = 0.0;
    double Lambda_C = 0.0, Lambda_V = 0.0, C1 = 0.0, C2 = 0.0, C3 = 0.0, C4 = 0.0;
    double mu0 = 0.0, mu1 = 0.0, mu2 = 0.0, mu3 = 0.0, mu4 = 0.0;

    // Following Neighbouring information is required for gradient computation
    /*
      (i-1,j+1)				(i,j+1)				(i+1,j+1)
      (i-1,j)					(i,j)				(i+1,j)
      (i-1,j-1)				(i,j-1)				(i+1,j-1)
   */
    // Obtaining The neighbours for a given cell
    //			(i-1,j)												(i,j-1)
    //	      		Neighbour_1 = Cell_Neighbours[Cell_No][1];	Neighbour_2 = Cell_Neighbours[Cell_No][2];
    //			(i+1,j)												(i,j+1)
    //		Neighbour_3 = Cell_Neighbours[Cell_No][3];	Neighbour_4 = Cell_Neighbours[Cell_No][4];

    Neighbour_1 = Cells[Cell_No].Neighbours[0];
    Neighbour_2 = Cells[Cell_No].Neighbours[1];
    Neighbour_3 = Cells[Cell_No].Neighbours[2];
    Neighbour_4 = Cells[Cell_No].Neighbours[3];
    //		cout<<Neighbour_1<<"\t"<<Neighbour_2<<"\t"<<Neighbour_3<<"\t"<<Neighbour_4<<endl;

    Face_No = 0;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_0 = Cells[Cell_No].Face_Normals[index + 0];
    ny_0 = Cells[Cell_No].Face_Normals[index + 1];
    dl_0 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 1;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_1 = Cells[Cell_No].Face_Normals[index + 0];
    ny_1 = Cells[Cell_No].Face_Normals[index + 1];
    dl_1 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 2;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_2 = Cells[Cell_No].Face_Normals[index + 0];
    ny_2 = Cells[Cell_No].Face_Normals[index + 1];
    dl_2 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 3;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    // Delta_Sm as per Blazek's Notation
    nx_3 = Cells[Cell_No].Face_Normals[index + 0];
    ny_3 = Cells[Cell_No].Face_Normals[index + 1];
    dl_3 = Cells[Cell_No].Face_Areas[Face_No];

    Avg_nx_i = 0.5 * (nx_2 - nx_0);
    Avg_ny_i = 0.5 * (ny_2 - ny_0);
    Avg_dl_i = 0.5 * (dl_0 + dl_2);

    Avg_nx_j = 0.5 * (nx_3 - nx_1);
    Avg_ny_j = 0.5 * (ny_3 - ny_1);
    Avg_dl_j = 0.5 * (dl_3 + dl_1);

    // Velocity and Acoustic Components of the Cell
    u0 = Primitive_Cells[Cell_No][1];
    v0 = Primitive_Cells[Cell_No][2];
    C0 = Primitive_Cells[Cell_No][5];
    Rho0 = Primitive_Cells[Cell_No][0];
    mu0 = Primitive_Cells[Cell_No][8];

    u1 = Primitive_Cells[Neighbour_1][1];
    v1 = Primitive_Cells[Neighbour_1][2];
    C1 = Primitive_Cells[Neighbour_1][5];
    Rho1 = Primitive_Cells[Neighbour_1][0];
    mu1 = Primitive_Cells[Neighbour_1][8];

    u2 = Primitive_Cells[Neighbour_2][1];
    v2 = Primitive_Cells[Neighbour_2][2];
    C2 = Primitive_Cells[Neighbour_2][5];
    Rho2 = Primitive_Cells[Neighbour_2][0];
    mu2 = Primitive_Cells[Neighbour_2][8];

    u3 = Primitive_Cells[Neighbour_3][1];
    v3 = Primitive_Cells[Neighbour_3][2];
    C3 = Primitive_Cells[Neighbour_3][5];
    Rho3 = Primitive_Cells[Neighbour_3][0];
    mu3 = Primitive_Cells[Neighbour_3][8];

    u4 = Primitive_Cells[Neighbour_4][1];
    v4 = Primitive_Cells[Neighbour_4][2];
    C4 = Primitive_Cells[Neighbour_4][5];
    Rho4 = Primitive_Cells[Neighbour_4][0];
    mu4 = Primitive_Cells[Neighbour_4][8];

    Inv_Area = Cells[Cell_No].Inv_Area;

    //			cout<<mu0<<"\t"<<mu1<<"\t"<<mu2<<"\t"<<mu3<<"\t"<<mu4<<"\t"<<mu<<endl;
    mu = (mu0 + mu1 + mu2 + mu3 + mu4) / 5; // Primitive_Cells[Cell_No][8];

    C11 = 4.0;
    Lambda_C = (fabs(0.5 * (u0 + u1) * nx_0 + 0.5 * (v0 + v1) * ny_0) + 0.5 * (C0 + C1)) * dl_0 + (fabs(0.5 * (u0 + u2) * nx_1 + 0.5 * (v0 + v2) * ny_1) + 0.5 * (C0 + C2)) * dl_1 + (fabs(0.5 * (u0 + u3) * nx_2 + 0.5 * (v0 + v3) * ny_2) + 0.5 * (C0 + C3)) * dl_2 + (fabs(0.5 * (u0 + u4) * nx_3 + 0.5 * (v0 + v4) * ny_3) + 0.5 * (C0 + C4)) * dl_3;

    Lambda_V = Inv_Area * (max((4.0 / (3.0 * (0.5 * (Rho0 + Rho1)))), (gamma / (0.5 * (Rho0 + Rho1)))) * mu * Inv_Pr * dl_0 * dl_0 + max((4.0 / (3.0 * (0.5 * (Rho0 + Rho2)))), (gamma / (0.5 * (Rho0 + Rho2)))) * mu * Inv_Pr * dl_1 * dl_1 + max((4.0 / (3.0 * (0.5 * (Rho0 + Rho3)))), (gamma / (0.5 * (Rho0 + Rho3)))) * mu * Inv_Pr * dl_2 * dl_2 + max((4.0 / (3.0 * (0.5 * (Rho0 + Rho4)))), (gamma / (0.5 * (Rho0 + Rho4)))) * mu * Inv_Pr * dl_3 * dl_3);

    Cells[Cell_No].del_t = CFL * Cells[Cell_No].Area / (Lambda_C + C11 * Lambda_V);
    //			cout<< dl_0*dl_0 <<"\t"<<dl_1*dl_1<<"\t"<<dl_2*dl_2<<"\t"<<dl_3*dl_3<<"\t"<<1.0/Inv_Area<<endl;
    //			cout<< (dl_0*dl_0  + dl_1*dl_1 + dl_2*dl_2 + dl_3*dl_3)*Inv_Area<<endl;
    //			cout<<Cell_No<<"\t"<<mu<<"\t"<<Lambda_C<<"\t"<<Lambda_V<<"\t"<<CFL<<"\t"<<Cells_Area[Cell_No]<<"\t"<<Cells[Cell_No].del_t<<endl;
    //			cout<<"----------------------------------------------------------------\n";
}

void Viscous_Time_Step_2(int &Cell_No)
{
    double nx_0 = 0.0, ny_0 = 0.0, dl_0 = 0.0, nx_1 = 0.0, ny_1 = 0.0, dl_1 = 0.0, nx_2 = 0.0, ny_2 = 0.0, dl_2 = 0.0, nx_3 = 0.0, ny_3 = 0.0, dl_3 = 0.0;
    double u0 = 0.0, v0 = 0.0, C0 = 0.0, Avg_dl_j = 0.0, Avg_dl_i = 0.0;
    double Inv_Area = 0.0, C11 = 0.0;
    int Face_No, index;
    Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0, Neighbour_4 = 0;
    double Rho0 = 0.0, P0, Term1, Term2, Term3;
    double Lambda_C = 0.0, Lambda_V = 0.0;
    double mu0 = 0.0;

    // Following Neighbouring information is required for gradient computation
    /*
        (i-1,j+1)				(i,j+1)				(i+1,j+1)
        (i-1,j)					(i,j)				(i+1,j)
        (i-1,j-1)				(i,j-1)				(i+1,j-1)
      */
    // Obtaining The neighbours for a given cell
    //			(i-1,j)												(i,j-1)
    // Neighbour_1 = Cell_Neighbours[Cell_No][1];
    Neighbour_1 = Cells[Cell_No].Neighbours[0];
    Neighbour_2 = Cells[Cell_No].Neighbours[1];
    Neighbour_3 = Cells[Cell_No].Neighbours[2];
    Neighbour_4 = Cells[Cell_No].Neighbours[3];
    //		cout<<Neighbour_1<<"\t"<<Neighbour_2<<"\t"<<Neighbour_3<<"\t"<<Neighbour_4<<endl;

    Face_No = 0;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_0 = Cells[Cell_No].Face_Normals[index + 0];
    ny_0 = Cells[Cell_No].Face_Normals[index + 1];
    dl_0 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 1;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_1 = Cells[Cell_No].Face_Normals[index + 0];
    ny_1 = Cells[Cell_No].Face_Normals[index + 1];
    dl_1 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 2;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_2 = Cells[Cell_No].Face_Normals[index + 0];
    ny_2 = Cells[Cell_No].Face_Normals[index + 1];
    dl_2 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 3;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_3 = Cells[Cell_No].Face_Normals[index + 0];
    ny_3 = Cells[Cell_No].Face_Normals[index + 1];
    dl_3 = Cells[Cell_No].Face_Areas[Face_No];

    // Delta_Sm as per Blazek's Notation

    Avg_dl_i = 0.5 * (dl_0 + dl_2); // Average dS in perpendicular to x direction .. Magnitude of Face Area

    Avg_dl_j = 0.5 * (dl_1 + dl_3); // Average dS in perpendicular to y direction .. Magnitude of Face Area
                                    //  cout<<Avg_dl_i<<"\t"<<Avg_dl_j<<endl;

    // Velocity and Acoustic Components of the Cell
    u0 = Primitive_Cells[Cell_No][1];
    v0 = Primitive_Cells[Cell_No][2];
    C0 = Primitive_Cells[Cell_No][5];
    Rho0 = Primitive_Cells[Cell_No][0];
    P0 = Primitive_Cells[Cell_No][4];
    mu0 = Primitive_Cells[Cell_No][8];
    Inv_Area = Cells[Cell_No].Inv_Area;

    // cout<<u0<<"\t"<<v0<<"\t"<<C0<<"\t"<<Rho0<<"\t"<<Inv_Area<<endl;

    C11 = 2.0;

    Term1 = (fabs(u0) + C0) * (0.25 * (dl_0 + dl_1 + dl_2 + dl_3));
    Term2 = (fabs(v0) + C0) * (0.25 * (dl_0 + dl_1 + dl_2 + dl_3));

    Lambda_C = Term1 + Term2;

    Term3 = (dl_0 * dl_0 + dl_1 * dl_1 + dl_2 * dl_2 + dl_3 * dl_3); /// Cells_Area[Cell_No];
    //		Term3 = (Avg_dl_i*Avg_dl_i + Avg_dl_j*Avg_dl_j)/Cells_Area[Cell_No];

    //		cout<<dl_0<<"\t"<<dl_1<<"\t"<<dl_2<<"\t"<<dl_3<<"\t"<<Cells_Area[Cell_No]<<endl;
    //		cout<<dl_0*dl_0 + dl_1*dl_1 + dl_2*dl_2 + dl_3*dl_3<<"\t"<<Cells_Area[Cell_No]<<endl;
    //		cout<<(dl_0*dl_0 + dl_1*dl_1 + dl_2*dl_2 + dl_3*dl_3)/Cells_Area[Cell_No]<<endl;
    // cout<<Term3<<endl;
    Lambda_V = max((4.0 / (3.0 * Rho0)), gamma / Rho0) * (mu0 / Pr) * Term3;
    Cells[Cell_No].del_t = CFL * Cells[Cell_No].Area / (Lambda_C + C11 * Lambda_V);
    //		cout<<Lambda_C<<"\t"<<Lambda_V<<endl;
    //				cout<<Cell_No<<"\t"<<mu0<<"\t"<<Pr<<"\t"<<Lambda_C<<"\t"<<Lambda_V<<"\t"<<Cells_DelT[Cell_No]<<endl;
    //		cout<<"------------------------------------------------\n";
}

void Viscous_Time_Step_3(int &Cell_No)
{
    double nx_0 = 0.0, ny_0 = 0.0, dl_0 = 0.0, nx_1 = 0.0, ny_1 = 0.0, dl_1 = 0.0, nx_2 = 0.0, ny_2 = 0.0, dl_2 = 0.0, nx_3 = 0.0, ny_3 = 0.0, dl_3 = 0.0;
    double u0 = 0.0, v0 = 0.0, C0 = 0.0, Avg_nx_i = 0.0, Avg_ny_i = 0.0, Avg_nx_j = 0.0, Avg_ny_j = 0.0, Lambda_x = 0.0, Lambda_y = 0.0, Avg_dl_j = 0.0, Avg_dl_i = 0.0;
    double Inv_Area = 0.0, tau = 0.0, Re_d = 0.0, Term1 = 0.0, Term2 = 0.0, Rho0, Rex, Rey;
    int Face_No, index;
    Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0, Neighbour_4 = 0;

    // Following Neighbouring information is required for gradient computation
    /*
      (i-1,j+1)				(i,j+1)				(i+1,j+1)
      (i-1,j)					(i,j)				(i+1,j)
      (i-1,j-1)				(i,j-1)				(i+1,j-1)
    */
    // Obtaining The neighbours for a given cell
    //			(i-1,j)												(i,j-1)
    /*Neighbour_1 = Cell_Neighbours[Cell_No][1];
    Neighbour_2 = Cell_Neighbours[Cell_No][2];
    //			(i+1,j)												(i,j+1)
    Neighbour_3 = Cell_Neighbours[Cell_No][3];
    Neighbour_4 = Cell_Neighbours[Cell_No][4];*/
    Neighbour_1 = Cells[Cell_No].Neighbours[0];
    Neighbour_2 = Cells[Cell_No].Neighbours[1];
    Neighbour_3 = Cells[Cell_No].Neighbours[2];
    Neighbour_4 = Cells[Cell_No].Neighbours[3];
    //		cout<<Neighbour_1<<"\t"<<Neighbour_2<<"\t"<<Neighbour_3<<"\t"<<Neighbour_4<<endl;

    Face_No = 0;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_0 = Cells[Cell_No].Face_Normals[index + 0];
    ny_0 = Cells[Cell_No].Face_Normals[index + 1];
    dl_0 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 1;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_1 = Cells[Cell_No].Face_Normals[index + 0];
    ny_1 = Cells[Cell_No].Face_Normals[index + 1];
    dl_1 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 2;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_2 = Cells[Cell_No].Face_Normals[index + 0];
    ny_2 = Cells[Cell_No].Face_Normals[index + 1];
    dl_2 = Cells[Cell_No].Face_Areas[Face_No];
    Face_No = 3;
    index = Face_No * 2;
    //      nx and ny are normals to the face and dl - face length
    nx_3 = Cells[Cell_No].Face_Normals[index + 0];
    ny_3 = Cells[Cell_No].Face_Normals[index + 1];
    dl_3 = Cells[Cell_No].Face_Areas[Face_No];

    Avg_nx_i = 0.5 * (nx_2 - nx_0);
    Avg_ny_i = 0.5 * (ny_2 - ny_0);
    Avg_dl_i = 0.5 * (dl_0 + dl_2);

    Avg_nx_j = 0.5 * (nx_3 - nx_1);
    Avg_ny_j = 0.5 * (ny_3 - ny_1);
    Avg_dl_j = 0.5 * (dl_3 + dl_1);

    // Velocity and Acoustic Components of the Cell
    u0 = Primitive_Cells[Cell_No][1];
    v0 = Primitive_Cells[Cell_No][2];
    C0 = Primitive_Cells[Cell_No][5];
    Rho0 = Primitive_Cells[Cell_No][0];
    mu = Primitive_Cells[Cell_No][8];

    // Inv_Area = Cells_Inv_Area[Cell_No];
    Inv_Area = Cells[Cell_No].Inv_Area;

    // cout<<Cell_No<<"\t"<<u<<"\t"<<v<<"\t"<<C<<"\t"<<"\t";

    // Maximum Eigen Value in x and y directions respectively as per BLAZEK 6.14 Formula
    Lambda_x = ((fabs(u0 * Avg_nx_i + v0 * Avg_ny_i) + C0) * Avg_dl_i);
    Lambda_y = ((fabs(u0 * Avg_nx_j + v0 * Avg_ny_j) + C0) * Avg_dl_j);

    tau = 1.0;
    Rex = (Rho0 * fabs(u0) + C0 * Avg_nx_i * Avg_dl_i) / mu;
    Rey = (Rho0 * fabs(v0) + C0 * Avg_ny_j * Avg_dl_j) / mu;
    Re_d = min(Rex, Rey);

    Term1 = 1.0 + 2.0 / Re_d;
    Term2 = (Lambda_x + Lambda_y) * Inv_Area;
    Cells[Cell_No].del_t = tau / (Term1 * Term2);
    //		cout<<Cell_No<<"\t"<<mu<<"\t"<<Re_d<<"\t"<<Term1<<"\t"<<Term2<<"\t"<<Re_d<<"\t"<<Cells_DelT[Cell_No]<<endl;
}

void Inviscid_Time_Step(int &Cell_No)
{
    // Generic inviscid time step for arbitrary 2D polygon cells:
    // dt = CFL * Area / sum_faces( (|V·n| + a) * dl )
    const int nFaces = (Cells[Cell_No].numFaces > 0)
                           ? Cells[Cell_No].numFaces
                           : static_cast<int>(Cells[Cell_No].Face_Areas.size());
    if (nFaces <= 0)
    {
        Cells[Cell_No].del_t = 0.0;
        return;
    }

    const double u0 = Primitive_Cells[Cell_No][1];
    const double v0 = Primitive_Cells[Cell_No][2];
    const double a0 = Primitive_Cells[Cell_No][5];

    double denom = 0.0;
    for (int face = 0; face < nFaces; face++)
    {
        const int idx = face * 2;
        const double nx_f = Cells[Cell_No].Face_Normals[idx + 0];
        const double ny_f = Cells[Cell_No].Face_Normals[idx + 1];
        const double dl_f = Cells[Cell_No].Face_Areas[face];

        const double vn = u0 * nx_f + v0 * ny_f;
        denom += (fabs(vn) + a0) * dl_f;
    }

    if (denom <= 0.0 || !std::isfinite(denom))
    {
        Cells[Cell_No].del_t = 0.0;
        return;
    }

    Cells[Cell_No].del_t = CFL * Cells[Cell_No].Area / denom;
}

void Evaluate_Time_Step(int &Cell_No)
{
    if (Is_Viscous_Wall)
    {
        // The existing viscous time-step formulations assume structured-like 4-face topology.
        // For non-quad cells, fall back to generic inviscid estimate to remain stable.
        if (Cells[Cell_No].numFaces != 4)
        {
            Inviscid_Time_Step(Cell_No);
            return;
        }
        switch (Viscous_Time_Case)
        {
        case 1:
            Viscous_Time_Step_1(Cell_No);
            break;
        case 2:
            Viscous_Time_Step_2(Cell_No);
            break;
        case 3:
            Viscous_Time_Step_3(Cell_No);
            break;
        }
    }
    else
    {
        Inviscid_Time_Step(Cell_No);
    }
}
