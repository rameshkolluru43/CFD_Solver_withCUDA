#include "definitions.h"
#include "Globals.h"
#include "Flux.h"
#include "Limiter.h"
#include "Utilities.h"

// This Function evaluates the amount of dissipation that is to be added on to the flux at a givne face, the face which is common to the current cell and the neighbour cell .
void LLF(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
  cout << "Evaluating LLF Dissipation" << endl;
  max_eigen_value = 0.0, u_L = 0.0, v_L = 0.0, u_R = 0.0, v_R = 0.0, nx = 0.0, ny = 0.0, Vdotn_L = 0.0, Vdotn_R = 0.0, dl = 0.0, C_L = 0.0, C_R = 0.0, mev_L = 0.0, mev_R = 0.0;

  //        V_D S(6,0.0);
  double Max1 = 0.0, Max2 = 0.0;

  //      Initializing the Left sate and Right State Conservative variables to be zero
  U_L[0] = 0.0;
  U_R[0] = 0.0;
  U_L[1] = 0.0;
  U_R[1] = 0.0;
  U_L[2] = 0.0;
  U_R[2] = 0.0;
  U_L[3] = 0.0;
  U_R[3] = 0.0;

  //      Initializing the Cell dissipation to be zero to avoid any type of over writing
  Dissipative_Flux[0] = 0.0;
  Dissipative_Flux[1] = 0.0;
  Dissipative_Flux[2] = 0.0;
  Dissipative_Flux[3] = 0.0;

  //      Left state Varialbes for each face in a given cell
  u_L = Primitive_Cells[Cell_No][1];
  v_L = Primitive_Cells[Cell_No][2];
  C_L = Primitive_Cells[Cell_No][5];

  U_L[0] = U_Cells[Cell_No][0];
  U_L[1] = U_Cells[Cell_No][1];
  U_L[2] = U_Cells[Cell_No][2];
  U_L[3] = U_Cells[Cell_No][3];

  //       cout<<"Left State\t"<<U_L[0]<<"\t"<<U_L[1]<<"\t"<<U_L[2]<<"\t"<<U_L[3]<<"\n";
  // Right State Variables
  U_R[0] = U_Cells[N_Cell_No][0];
  U_R[1] = U_Cells[N_Cell_No][1];
  U_R[2] = U_Cells[N_Cell_No][2];
  U_R[3] = U_Cells[N_Cell_No][3];

  //     cout<<"Right State\t"<<U_R[0]<<"\t"<<U_R[1]<<"\t"<<U_R[2]<<"\t"<<U_R[3]<<"\n";
  u_R = Primitive_Cells[N_Cell_No][1];
  v_R = Primitive_Cells[N_Cell_No][2];
  C_R = Primitive_Cells[N_Cell_No][5];

  nx = Cells[Cell_No].Face_Normals[Face_No * 2 + 0];
  ny = Cells[Cell_No].Face_Normals[Face_No * 2 + 1];
  dl = Cells[Cell_No].Face_Areas[Face_No];

  // Normal Velocity of left and Right state
  Vdotn_L = (u_L * nx + v_L * ny);
  Vdotn_R = (u_R * nx + v_R * ny);

  //              Wave Speed evaluation
  S[0] = fabs(Vdotn_L - C_L) * dl;
  S[1] = fabs(Vdotn_L + C_L) * dl;
  S[2] = fabs(Vdotn_L) * dl;
  S[3] = fabs(Vdotn_R - C_R) * dl;
  S[4] = fabs(Vdotn_R + C_R) * dl;
  S[5] = fabs(Vdotn_R) * dl;

  //              Finding Minimum and Maximum Wave speeds from neighbouring cells
  Maximum(S[0], S[1], S[2], Max1);      // Maximum of Left state Eigen values
  Maximum(S[3], S[4], S[5], Max2);      // Maximum of Right state Eigen values
  Maximum(Max1, Max2, max_eigen_value); // Maximum of Left state and Right state Eigen values

  // cout<<Max1<<"\t"<<Max2<<"\t"<<max_eigen_value<<endl;

  // Evaluation of Dissipation on Face

  Dissipative_Flux[0] = 0.5 * max_eigen_value * (U_R[0] - U_L[0]);
  Dissipative_Flux[1] = 0.5 * max_eigen_value * (U_R[1] - U_L[1]);
  Dissipative_Flux[2] = 0.5 * max_eigen_value * (U_R[2] - U_L[2]);
  Dissipative_Flux[3] = 0.5 * max_eigen_value * (U_R[3] - U_L[3]);
}

void LLF_2O(const int &Cell_No, int &N_Cell_No, const int &Face_No)
{
  //     cout<<"Evaluating LLF Dissipation"<<endl;
  max_eigen_value = 0.0, u_L = 0.0, v_L = 0.0, u_R = 0.0, v_R = 0.0, nx = 0.0, ny = 0.0, Vdotn_L = 0.0, Vdotn_R = 0.0, dl = 0.0, C_L = 0.0, C_R = 0.0, mev_L = 0.0, mev_R = 0.0;

  V_D S(6, 0.0);
  double Max1 = 0.0, Max2 = 0.0;
  int index = Face_No * 2;
  vector<double> d_U(4, 0.0);

  //      Initializing the Cell dissipation to be zero to avoid any type of over writing
  Dissipative_Flux[0] = 0.0;
  Dissipative_Flux[1] = 0.0;
  Dissipative_Flux[2] = 0.0;
  Dissipative_Flux[3] = 0.0;

  C_L = sqrt(gamma * P_L / Rho_L);
  C_R = sqrt(gamma * P_R / Rho_R);

  Vmag_L = 0.5 * (u_L * u_L + v_L * v_L);
  Vmag_R = 0.5 * (u_R * u_R + v_R * v_R);

  /*  nx = Cell_Face_Normals[Cell_No][Face_No*2+0];
    ny = Cell_Face_Normals[Cell_No][Face_No*2+1];
    dl = Cell_Face_Areas[Cell_No][Face_No];*/
  nx = Cells[Cell_No].Face_Normals[Face_No * 2 + 0];
  ny = Cells[Cell_No].Face_Normals[Face_No * 2 + 1];
  dl = Cells[Cell_No].Face_Areas[Face_No];

  // Normal Velocity of left and Right state
  Vdotn_L = (u_L * nx + v_L * ny);
  Vdotn_R = (u_R * nx + v_R * ny);

  //              Wave Speed evaluation
  S[0] = fabs(Vdotn_L - C_L) * dl;
  S[1] = fabs(Vdotn_L + C_L) * dl;
  S[2] = fabs(Vdotn_L) * dl;
  S[3] = fabs(Vdotn_R - C_R) * dl;
  S[4] = fabs(Vdotn_R + C_R) * dl;
  S[5] = fabs(Vdotn_R) * dl;

  //              Finding Minimum and Maximum Wave speeds from neighbouring cells
  Maximum(S[0], S[1], S[2], Max1); // Maximum of Left state Eigen values
  Maximum(S[3], S[4], S[5], Max2); // Maximum of Right state Eigen values

  Maximum(Max1, Max2, max_eigen_value); // Maximum of Left state and Right state Eigen values

  Second_Order_Limiter(Cell_No, Face_No, d_U);

  //              cout<<Max1<<"\t"<<Max2<<"\t"<<max_eigen_value<<endl;

  // Evaluation of Dissipation on Face

  Dissipative_Flux[0] = 0.5 * max_eigen_value * d_U[0];
  Dissipative_Flux[1] = 0.5 * max_eigen_value * d_U[1];
  Dissipative_Flux[2] = 0.5 * max_eigen_value * d_U[2];
  Dissipative_Flux[3] = 0.5 * max_eigen_value * d_U[3];
}
