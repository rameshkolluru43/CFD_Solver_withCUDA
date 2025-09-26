#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"
#include "Utilities.h"

void Entropy_Fix(double &Alpha, double &L_Max)
{
    double k = 1.0, delta = k * L_Max;

    if (fabs(Alpha) < delta)
        Alpha = (Alpha * Alpha + delta * delta) / (2.0 * delta);
    else
        Alpha = L_Max;

    //   cout<<Alpha<<"\t"<<L_Max<<endl;
}

void Condition_For_MOVERS(double &d_U, double &d_F, double &L_Max, double &L_Min, double &Alpha)
{
    double epsilon = 1e-10;
    //***************************************************************************************************************

    if (fabs(d_F) < epsilon and fabs(d_U) > epsilon) // Condition for discontinuity.... across which Flux is zero
    {
        Alpha = 0.0;
    }
    else if (fabs(d_F) < epsilon and fabs(d_U) < epsilon)
    {
        Alpha = L_Min;
    }
    else if (fabs(d_F) > epsilon and fabs(d_U) > epsilon)
    {
        Alpha = fabs(d_F / d_U);
        //      Wave speed correction limiting the L_min<=|alpha|<=L_Max
        if (Alpha >= L_Max)
            Alpha = Sign(Alpha) * L_Max;
        else if (Alpha <= L_Min)
            Alpha = Sign(Alpha) * L_Min;
    }
    else
        Alpha = L_Min;
}

void MOVERS(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{

    Rho_L = 0.0, P_L = 0.0, C_L = 0.0, u_L = 0.0, v_L = 0.0, Vdotn_L = 0.0, nx = 0.0, ny = 0.0, Vmag_L = 0.0;
    Rho_R = 0.0, P_R = 0.0, C_R = 0.0, u_R = 0.0, v_R = 0.0, Vdotn_R = 0.0, Vmag_R = 0.0, dl = 0.0;

    Mod_Alpha0 = 0.0, Mod_Alpha1 = 0.0, Mod_Alpha2 = 0.0, Mod_Alpha3 = 0.0, d_F_0 = 0.0, d_F_1 = 0.0, d_F_2 = 0.0, d_F_3 = 0.0, d_U_0 = 0.0, d_U_1 = 0.0, d_U_2 = 0.0, d_U_3 = 0.0;
    Lambda_Max = 0.0, Lambda_Min = 0.0, Max1 = 0.0, Max2 = 0.0, Min1 = 0.0, Min2 = 0.0, d_Var = 0.0;
    int index = Face_No * 2;

    mev_L = 0.0, mev_R = 0.0, max_eigen_value = 0.0;

    V_D S(6, 0.0);
    U_L[0] = 0.0;
    U_L[1] = 0.0;
    U_L[2] = 0.0;
    U_L[3] = 0.0;
    U_R[0] = 0.0;
    U_R[1] = 0.0;
    U_R[2] = 0.0;
    U_R[3] = 0.0;

    Dissipative_Flux[0] = 0.0;
    Dissipative_Flux[1] = 0.0;
    Dissipative_Flux[2] = 0.0;
    Dissipative_Flux[3] = 0.0;

    // 	cout<<Cell_No<<"\t"<<N_Cell_No<<"\t";

    //      Left state Variables, Density, Pressure, u, v, speed of cound C

    Rho_L = Primitive_Cells[Cell_No][0];
    P_L = Primitive_Cells[Cell_No][4];
    u_L = Primitive_Cells[Cell_No][1];
    v_L = Primitive_Cells[Cell_No][2];
    C_L = Primitive_Cells[Cell_No][5];

    //      Right state Variables, Density, Pressure, u, v, speed of cound C
    Rho_R = Primitive_Cells[N_Cell_No][0];
    P_R = Primitive_Cells[N_Cell_No][4];
    u_R = Primitive_Cells[N_Cell_No][1];
    v_R = Primitive_Cells[N_Cell_No][2];
    C_R = Primitive_Cells[N_Cell_No][5];
    //		Conserved Variables from Left and Right States
    U_L = U_Cells[Cell_No];
    U_R = U_Cells[N_Cell_No];

    //          Evaluating Conserved variable Difference between Left and Right faces

    d_U_0 = (U_R[0] - U_L[0]);
    d_U_1 = (U_R[1] - U_L[1]);
    d_U_2 = (U_R[2] - U_L[2]);
    d_U_3 = (U_R[3] - U_L[3]);

    //      nx and ny are normals to the face and dl - face length
    /*nx = Cell_Face_Normals[Cell_No][index + 0];
    ny = Cell_Face_Normals[Cell_No][index + 1];
    dl = Cell_Face_Areas[Cell_No][Face_No];*/

    nx = Cells[Cell_No].Face_Normals[index + 0];
    ny = Cells[Cell_No].Face_Normals[index + 1];
    dl = Cells[Cell_No].Face_Areas[Face_No];

    //      Normal Velocity taking Left state Velocity and Right State Velocity of an interface
    Vdotn_L = (u_L * nx + v_L * ny);
    Vdotn_R = (u_R * nx + v_R * ny);

    //      Magnitudes of Velocity on Left side and Right side of an interface
    Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    //          Evaluating Flux Difference between Left and Right Faces
    d_F_0 = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
    d_F_1 = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L + P_L * nx) * dl;
    d_F_2 = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L + P_L * ny) * dl;
    d_F_3 = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R - (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) * dl;

    //          Wave Speed evaluation
    S[0] = fabs(Vdotn_L - C_L) * dl;
    S[1] = fabs(Vdotn_L + C_L) * dl;
    S[2] = fabs(Vdotn_L) * dl;
    S[3] = fabs(Vdotn_R - C_R) * dl;
    S[4] = fabs(Vdotn_R + C_R) * dl;
    S[5] = fabs(Vdotn_R) * dl;

    //              Finding Minimum and Maximum Wave speeds from neighbouring cells
    Maximum(S[0], S[1], S[2], Max1); // Maximum of Left state Eigen values
    Maximum(S[3], S[4], S[5], Max2); // Maximum of Right state Eigen values

    Maximum(Max1, Max2, Lambda_Max); // Maximum of Left state and Right state Eigen values

    Minimum(S[0], S[1], S[2], Min1); // Minimum of Left state Eigen values
    Minimum(S[3], S[4], S[5], Min2); // Minimum of Right state Eigen values

    Maximum(Min1, Min2, Lambda_Min); // Minimum of Left state and Right state Eigen values

    //      cout<<Lambda_Max<<"\t"<<Lambda_Min<<endl;
    Condition_For_MOVERS(d_U_3, d_F_3, Lambda_Max, Lambda_Min, Mod_Alpha3);
    Mod_Alpha0 = Mod_Alpha3;
    Mod_Alpha1 = Mod_Alpha3;
    Mod_Alpha2 = Mod_Alpha3;

    if (!Is_MOVERS_1)
    {
        Condition_For_MOVERS(d_U_0, d_F_0, Lambda_Max, Lambda_Min, Mod_Alpha0);
        Condition_For_MOVERS(d_U_1, d_F_1, Lambda_Max, Lambda_Min, Mod_Alpha1);
        Condition_For_MOVERS(d_U_2, d_F_2, Lambda_Max, Lambda_Min, Mod_Alpha2);
    }
    if (Enable_Entropy_Fix)
    {
        //         cout<<"Entropy_Fix Enabled"<<endl;
        Entropy_Fix(Mod_Alpha0, Lambda_Max);
        Entropy_Fix(Mod_Alpha1, Lambda_Max);
        Entropy_Fix(Mod_Alpha2, Lambda_Max);
        Entropy_Fix(Mod_Alpha3, Lambda_Max);
    }

    //      cout<< Mod_Alpha0<<"\t"<< Mod_Alpha1<<"\t"<<Mod_Alpha2<<"\t"<<Mod_Alpha3<<endl;
    Dissipative_Flux[0] = 0.5 * Mod_Alpha0 * d_U_0;
    Dissipative_Flux[1] = 0.5 * Mod_Alpha1 * d_U_1;
    Dissipative_Flux[2] = 0.5 * Mod_Alpha2 * d_U_2;
    Dissipative_Flux[3] = 0.5 * Mod_Alpha3 * d_U_3;
}

void MOVERS_2O(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{

    Rho_L = 0.0, P_L = 0.0, C_L = 0.0, u_L = 0.0, v_L = 0.0, Vdotn_L = 0.0, nx = 0.0, ny = 0.0, Vmag_L = 0.0;
    Rho_R = 0.0, P_R = 0.0, C_R = 0.0, u_R = 0.0, v_R = 0.0, Vdotn_R = 0.0, Vmag_R = 0.0, dl = 0.0;

    Mod_Alpha0 = 0.0, Mod_Alpha1 = 0.0, Mod_Alpha2 = 0.0, Mod_Alpha3 = 0.0, d_F_0 = 0.0, d_F_1 = 0.0, d_F_2 = 0.0, d_F_3 = 0.0, d_U_0 = 0.0, d_U_1 = 0.0, d_U_2 = 0.0, d_U_3 = 0.0;
    Lambda_Max = 0.0, Lambda_Min = 0.0, Max1 = 0.0, Max2 = 0.0, Min1 = 0.0, Min2 = 0.0;
    double d_Var_L = 0.0, d_Var_R = 0.0;
    int index = Face_No * 2, Which_Var = 0;

    mev_L = 0.0, mev_R = 0.0, max_eigen_value = 0.0;

    V_D S(6, 0.0);

    Dissipative_Flux[0] = 0.0;
    Dissipative_Flux[1] = 0.0;
    Dissipative_Flux[2] = 0.0;
    Dissipative_Flux[3] = 0.0;

    // 	cout<<Cell_No<<"\t"<<N_Cell_No<<"\t";

    //      Left state Variables, Density, Pressure, u, v, speed of cound C

    Rho_L = Primitive_Cells[Cell_No][0];
    P_L = Primitive_Cells[Cell_No][4];
    u_L = Primitive_Cells[Cell_No][1];
    v_L = Primitive_Cells[Cell_No][2];

    //      Right state Variables, Density, Pressure, u, v, speed of cound C
    Rho_R = Primitive_Cells[N_Cell_No][0];
    P_R = Primitive_Cells[N_Cell_No][4];
    u_R = Primitive_Cells[N_Cell_No][1];
    v_R = Primitive_Cells[N_Cell_No][2];

    d_U_0 = Rho_R - Rho_L;
    d_U_1 = Rho_R * u_R - Rho_L * u_L;
    d_U_2 = Rho_R * v_R - Rho_L * v_L;
    d_U_3 = ((P_R / (gamma - 1.0)) + 0.5 * Rho_R * (u_R * u_R + v_R * v_R)) - ((P_L / (gamma - 1.0)) + 0.5 * Rho_L * (u_L * u_L + v_L * v_L));

    //      nx and ny are normals to the face and dl - face length
    /*nx = Cell_Face_Normals[Cell_No][index + 0];
    ny = Cell_Face_Normals[Cell_No][index + 1];
    dl = Cell_Face_Areas[Cell_No][Face_No];*/
    nx = Cells[Cell_No].Face_Normals[index + 0];
    ny = Cells[Cell_No].Face_Normals[index + 1];
    dl = Cells[Cell_No].Face_Areas[Face_No];

    //      Normal Velocity taking Left state Velocity and Right State Velocity of an interface
    Vdotn_L = (u_L * nx + v_L * ny);
    Vdotn_R = (u_R * nx + v_R * ny);

    //      Magnitudes of Velocity on Left side and Right side of an interface
    Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
    Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

    //          Evaluating Flux Difference between Left and Right Faces
    d_F_0 = (Rho_R * Vdotn_R - Rho_L * Vdotn_L) * dl;
    d_F_1 = (Rho_R * u_R * Vdotn_R + P_R * nx - Rho_L * u_L * Vdotn_L + P_L * nx) * dl;
    d_F_2 = (Rho_R * v_R * Vdotn_R + P_R * ny - Rho_L * v_L * Vdotn_L + P_L * ny) * dl;
    d_F_3 = ((((P_R / (gamma - 1.0)) + Rho_R * Vmag_R) + P_R) * Vdotn_R - (((P_L / (gamma - 1.0)) + Rho_L * Vmag_L) + P_L) * Vdotn_L) * dl;

    //          Wave Speed evaluation
    S[0] = fabs(Vdotn_L - C_L) * dl;
    S[1] = fabs(Vdotn_L + C_L) * dl;
    S[2] = fabs(Vdotn_L) * dl;
    S[3] = fabs(Vdotn_R - C_R) * dl;
    S[4] = fabs(Vdotn_R + C_R) * dl;
    S[5] = fabs(Vdotn_R) * dl;

    //              Finding Minimum and Maximum Wave speeds from neighbouring cells
    Maximum(S[0], S[1], S[2], Max1); // Maximum of Left state Eigen values
    Maximum(S[3], S[4], S[5], Max2); // Maximum of Right state Eigen values

    Maximum(Max1, Max2, Lambda_Max); // Maximum of Left state and Right state Eigen values

    Minimum(S[0], S[1], S[2], Min1); // Minimum of Left state Eigen values
    Minimum(S[3], S[4], S[5], Min2); // Minimum of Right state Eigen values

    Maximum(Min1, Min2, Lambda_Min); // Minimum of Left state and Right state Eigen values

    //      cout<<Lambda_Max<<"\t"<<Lambda_Min<<endl;
    Condition_For_MOVERS(d_U_3, d_F_3, Lambda_Max, Lambda_Min, Mod_Alpha3);
    Mod_Alpha0 = Mod_Alpha3;
    Mod_Alpha1 = Mod_Alpha3;
    Mod_Alpha2 = Mod_Alpha3;

    if (!Is_MOVERS_1)
    {
        Condition_For_MOVERS(d_U_0, d_F_0, Lambda_Max, Lambda_Min, Mod_Alpha0);
        Condition_For_MOVERS(d_U_1, d_F_1, Lambda_Max, Lambda_Min, Mod_Alpha1);
        Condition_For_MOVERS(d_U_2, d_F_2, Lambda_Max, Lambda_Min, Mod_Alpha2);
    }
    if (Enable_Entropy_Fix)
    {
        //         cout<<"Entropy_Fix Enabled"<<endl;
        //		      		  Entropy_Fix(Mod_Alpha0,Lambda_Max);
        //					  Entropy_Fix(Mod_Alpha1,Lambda_Max);
        //					  Entropy_Fix(Mod_Alpha2,Lambda_Max);
        Entropy_Fix(Mod_Alpha3, Lambda_Max);
    }

    //      cout<< Mod_Alpha0<<"\t"<< Mod_Alpha1<<"\t"<<Mod_Alpha2<<"\t"<<Mod_Alpha3<<endl;
    Dissipative_Flux[0] = 0.5 * Mod_Alpha0 * d_U_0;
    Dissipative_Flux[1] = 0.5 * Mod_Alpha1 * d_U_1;
    Dissipative_Flux[2] = 0.5 * Mod_Alpha2 * d_U_2;
    Dissipative_Flux[3] = 0.5 * Mod_Alpha3 * d_U_3;
}
