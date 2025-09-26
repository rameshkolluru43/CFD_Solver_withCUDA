// File: Assemble_Matrix.cpp
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13
#include "definitions.h"
#include "Globals.h"
#include "Flux_Jacobian.h"
#include "Assemble_Matrix.h"

// Function to assemble the matrix A and vector b in terms of total number of cells

// V_D b;

vector<V_D> Assemble_A(vector<V_D> &A, double &dt)
{
    //	Ac.resize(4*No_Physical_Cells,V_D(4*No_Physical_Cells,0.0));
    // Loop over each cell by flattening the 2D index into a 1D index
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        //            int cell_index = (i - 1) * (ny - 1) + (j - 1); // Compute the 1D index for the current cell
        //	    cout<<"Assembling the matrix for Cell Number\t"<<Cell_No<<endl;
        // Get the state variables for the current cell
        double Omega = Cells[Cell_No].Area;
        int Face_No;

        Neighbour_1 = Cells[Cell_No].Neighbours[0]; //(i-1,j,k)
        Neighbour_2 = Cells[Cell_No].Neighbours[1]; //(i,j-1,k)
        Neighbour_3 = Cells[Cell_No].Neighbours[2]; //(i+1,j,k)
        Neighbour_4 = Cells[Cell_No].Neighbours[3]; //(i,j+1,k)

        // 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

        // Scale the Jacobians by the face lengths for the specific cell
        double dx_right = Cells[Cell_No].Face_Areas[2];  // Right face length in x-direction
        double dx_left = Cells[Cell_No].Face_Areas[0];   // Left face length in x-direction
        double dy_top = Cells[Cell_No].Face_Areas[3];    // Top face length in y-direction
        double dy_bottom = Cells[Cell_No].Face_Areas[1]; // Bottom face length in y-direction

        // Compute the Jacobian for the current cell
        // Flux Jacobian in x direction
        Face_No = 0;
        A_x_L = Compute_Flux_Jacobian(Cell_No, A_x_L, Face_No);
        Face_No = 2;
        A_x_R = Compute_Flux_Jacobian(Cell_No, A_x_R, Face_No);

        // Flux Jacobian in y direction
        Face_No = 1;
        A_y_T = Compute_Flux_Jacobian(Cell_No, A_y_T, Face_No);
        Face_No = 3;
        A_y_B = Compute_Flux_Jacobian(Cell_No, A_y_B, Face_No);

        // Scale the Jacobians by the face lengths
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                A_x_L[row][col] *= dy_bottom; // Scale by face length in y-direction for x-fluxes
                A_y_T[row][col] *= dx_right;  // Scale by face length in x-direction for y-fluxes
                A_x_R[row][col] *= dy_top;    // Scale by face length in y-direction for x-fluxes
                A_y_B[row][col] *= dx_left;   // Scale by face length in x-direction for y-fluxes
            }
        }

        // Add the time derivative term (Omega / dt * I) to the diagonal
        for (int d = 0; d < 4; ++d)
        {
            A[4 * Cell_No + d][4 * Cell_No + d] += Omega / dt; // Identity matrix scaling
        }
        // Loop through the conserved variables (4 variables per cell: rho, rho*u, rho*v, E)
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                // Self contribution (current cell)
                A[4 * Cell_No + row][4 * Cell_No + col] +=
                    (dt / (2.0 * (dx_left + dx_right))) * (A_x_L[row][col] + A_x_R[row][col]) + (dt / (2.0 * (dy_bottom + dy_top))) * (A_y_T[row][col] + A_y_B[row][col]);
            }
        }

        // Right neighbor (i+1, j)
        Face_No = 2;
        A_x = ComputeGhostCell_Flux_Jacobian(Neighbour_3, Cell_No, A_x, Face_No);
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                A_x[row][col] *= dx_right;
                A[4 * Cell_No + row][4 * Neighbour_3 + col] +=
                    (dt / (2.0 * dx_right)) * A_x[row][col];
            }
        }

        // Left neighbor (i-1, j)
        Face_No = 0;
        A_x = ComputeGhostCell_Flux_Jacobian(Neighbour_1, Cell_No, A_x, Face_No);
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                A_x[row][col] *= dx_left;
                A[4 * Cell_No + row][4 * Neighbour_1 + col] -=
                    (dt / (2.0 * dx_left)) * A_x[row][col];
            }
        }

        Face_No = 3;
        A_y = ComputeGhostCell_Flux_Jacobian(Neighbour_4, Cell_No, A_y, Face_No);
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                // Upper neighbor (i, j+1)
                A_y[row][col] *= dy_top;
                A[4 * Cell_No + row][4 * Neighbour_4 + col] +=
                    (dt / (2.0 * dy_top)) * A_y[row][col];
            }
        }

        // Lower neighbor (i, j-1)
        Face_No = 1;
        A_y = ComputeGhostCell_Flux_Jacobian(Neighbour_2, Cell_No, A_y, Face_No);
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                A_y[row][col] *= dy_bottom;
                A[4 * Cell_No + row][4 * Neighbour_2 + col] -=
                    (dt / (2.0 * dy_bottom)) * A_y[row][col];
            }
        }
    }
    return A;
}

// Function to assemble the vector b for implicit methods
// Ensure that these functions are called before assembling  Evaluate_Cell_Net_Flux_2O(); Evaluate_Cell_Net_Flux_1O();

V_D Assemble_b(V_D &b)
{
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        for (int d = 0; d < 4; ++d)
        {
            b[4 * Cell_No + d] = -Cells_Net_Flux[Cell_No][d];
        }
    }
    return b;
}

void Assemble_A1(double &dt)
{
    row_indices.clear();
    col_indices.clear();
    values.clear();
    int Face_No = 0;
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        double Omega = Cells[Cell_No].Area;
        int Neighbour_1 = Cells[Cell_No].Neighbours[0]; // Left neighbor (i-1, j)
        int Neighbour_2 = Cells[Cell_No].Neighbours[1]; // Bottom neighbor (i, j-1)
        int Neighbour_3 = Cells[Cell_No].Neighbours[2]; // Right neighbor (i+1, j)
        int Neighbour_4 = Cells[Cell_No].Neighbours[3]; // Top neighbor (i, j+1)

        double dx_right = Cells[Cell_No].Face_Areas[2];  // Right face length in x-direction
        double dx_left = Cells[Cell_No].Face_Areas[0];   // Left face length in x-direction
        double dy_top = Cells[Cell_No].Face_Areas[3];    // Top face length in y-direction
        double dy_bottom = Cells[Cell_No].Face_Areas[1]; // Bottom face length in y-direction

        // Compute the Jacobians for the current cell
        Face_No = 0;
        A_x_L = Compute_Flux_Jacobian(Cell_No, A_x_L, Face_No); // Left face (i-1)
        Face_No = 2;
        A_x_R = Compute_Flux_Jacobian(Cell_No, A_x_R, Face_No); // Right face (i+1)
        Face_No = 1;
        A_y_B = Compute_Flux_Jacobian(Cell_No, A_y_B, Face_No); // Bottom face (i, j-1)
        Face_No = 3;
        A_y_T = Compute_Flux_Jacobian(Cell_No, A_y_T, Face_No); // Top face (i, j+1)

        // Update the center cell and add time derivative term to diagonal elements
        for (int row = 0; row < 4; ++row)
        {

            for (int col = 0; col < 4; ++col)
            {

                row_indices.push_back(4 * Cell_No + row);
                col_indices.push_back(4 * Cell_No + col);

                // Calculate the flux contributions
                double flux_contrib = (dt / (2.0 * (dy_bottom + dy_top))) * (A_x_L[row][col] + A_x_R[row][col]) +
                                      (dt / (2.0 * (dx_left + dx_right))) * (A_y_T[row][col] + A_y_B[row][col]);

                // Add time derivative term to diagonal elements
                if (row == col)
                {
                    values.push_back(flux_contrib + Omega / dt); // Add diagonal contribution
                }
                else
                {
                    values.push_back(flux_contrib); // Non-diagonal entries, only flux contribution
                }
            }
        }

        // Right neighbor (i+1, j)
        if (Neighbour_3 >= No_Physical_Cells)
        {
            // cout<<Neighbour_3<<"\t"<<Cell_No<<endl;
            Face_No = 2;
            A_x = ComputeGhostCell_Flux_Jacobian(Neighbour_3, Cell_No, A_x, Face_No);
            for (int row = 0; row < 4; ++row)
            {
                for (int col = 0; col < 4; ++col)
                {
                    row_indices.push_back(4 * Cell_No + row);
                    col_indices.push_back(4 * Neighbour_3 + col);
                    values.push_back((dt / (2.0 * dx_right)) * A_x[row][col]);
                }
            }
        }

        // Left neighbor (i-1, j)
        if (Neighbour_1 >= No_Physical_Cells)
        {
            // cout<<Neighbour_1<<"\t"<<Cell_No<<endl;
            Face_No = 0;
            A_x = ComputeGhostCell_Flux_Jacobian(Neighbour_1, Cell_No, A_x, Face_No);
            for (int row = 0; row < 4; ++row)
            {
                for (int col = 0; col < 4; ++col)
                {
                    row_indices.push_back(4 * Cell_No + row);
                    col_indices.push_back(4 * Neighbour_1 + col);
                    values.push_back(-(dt / (2.0 * dx_left)) * A_x[row][col]);
                }
            }
        }

        // Top neighbor (i, j+1)
        if (Neighbour_4 >= No_Physical_Cells)
        {
            // cout<<Neighbour_4<<"\t"<<Cell_No<<endl;
            Face_No = 3;
            A_y = ComputeGhostCell_Flux_Jacobian(Neighbour_4, Cell_No, A_y, Face_No);
            for (int row = 0; row < 4; ++row)
            {
                for (int col = 0; col < 4; ++col)
                {
                    row_indices.push_back(4 * Cell_No + row);
                    col_indices.push_back(4 * Neighbour_4 + col);
                    values.push_back((dt / (2.0 * dy_top)) * A_y[row][col]);
                }
            }
        }

        // Bottom neighbor (i, j-1)
        if (Neighbour_2 >= No_Physical_Cells)
        {
            // cout<<Neighbour_2<<"\t"<<Cell_No<<endl;
            Face_No = 1;
            A_y = ComputeGhostCell_Flux_Jacobian(Neighbour_2, Cell_No, A_y, Face_No);
            for (int row = 0; row < 4; ++row)
            {
                for (int col = 0; col < 4; ++col)
                {
                    row_indices.push_back(4 * Cell_No + row);
                    col_indices.push_back(4 * Neighbour_2 + col);
                    values.push_back(-(dt / (2.0 * dy_bottom)) * A_y[row][col]);
                }
            }
        }
    }

    // cout<<row_indices.size()<<"\t"<<col_indices.size()<<"\t"<<values.size()<<endl;
}

// Functions to get the stored sparse matrix data
vector<int> get_row_indices()
{
    return row_indices;
}

vector<int> get_col_indices()
{
    return col_indices;
}

vector<double> get_Values()
{
    return values;
}
