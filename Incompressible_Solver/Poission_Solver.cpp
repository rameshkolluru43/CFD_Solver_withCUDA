#include "definitions.h"
#include "Poisson_Solver.h"
#include "Globals.h"
void solvePressurePoissonExplicit(int &No_Physical_Cells, int &Variable)
{
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0; // Coefficients for the Poisson equation
    V_D df(3, 0.0);                            // Placeholder for distance vector

    // Assemble the coefficient matrix and RHS vector
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        V_D Coeff(Cells[Cell_No].numFaces, 0.0), Phi(Cells[Cell_No].numFaces, 0.0); // Coefficients for the Poisson equation
        for (int Face_No = 0; Face_No < Cells[Cell_No].numFaces; Face_No++)
        {
            int neighbor = Cells[Cell_No].Neighbours[Face_No];          // Neighbor cell index corresponding to the face
            int index = Face_No * 2;                                    // Index for accessing face normals
            double nx = Cells[Cell_No].Face_Normals[index];             // x-component of face normal
            double ny = Cells[Cell_No].Face_Normals[index + 1];         // y-component of face normal
            double A_f = Cells[Cell_No].Face_Areas[Face_No];            // Face area
            double d_f = Cells[Cell_No].Cell_Center_Distances[Face_No]; // Placeholder: Distance between cell centers (compute appropriately)
            df[0] = Cells[Cell_No].Cell_Center_Vector[index];           // cell center vector df_x
            df[1] = Cells[Cell_No].Cell_Center_Vector[index + 1];       // cell center vector df_y
            // Compute the coefficients for the Poisson equation
            // Coefficients are based on the face normals and the normalized cell center vector
            Coeff[Face_No] = ((nx * df[0] + ny * df[1]) / d_f) * A_f; // Coefficient for the Poisson equation
            Phi[Face_No] = Primitive_Cells[neighbor][Variable];       // Placeholder: Neighbor cell pressure
        }
        double Area = Cells[Cell_No].Area; // Area of the cell

        // Source is the source term for the Poisson equation
        double Source = 0.0;

        double sum = 1.0 / Coeff[0] + Coeff[1] + Coeff[2] + Coeff[3];
        double Phi_New = sum * (Coeff[0] * Phi[0] + Coeff[1] * Phi[1] + Coeff[2] * Phi[2] + Coeff[3] * Phi[3] + Source * Area);
    }
}