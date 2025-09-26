void Velocity_Predictor_Step()
{
    double convective_u = 0.0, convective_v = 0.0, diffusive_u = 0.0, diffusive_v = 0.0;
    double pressure_grad_u = 0.0, pressure_grad_v = 0.0;
    double uf = 0.0, vf = 0.0;
    double u_diff = 0.0, v_diff = 0.0;
    double inv_Area = 0.0;
    double grad_u = 0.0, grad_v = 0.0;
    const int Grad_Type =

        /*	Neighbour_1 = Cell_Neighbours[Cell_No][1];		//(i-1,j,k) Left Neighbour
            Neighbour_2 = Cell_Neighbours[Cell_No][2];		//(i,j-1,k) bottom Neighbour
            Neighbour_3 = Cell_Neighbours[Cell_No][3];		//(i+1,j,k) right Neighbour
            Neighbour_4 = Cell_Neighbours[Cell_No][4];		//(i,j+1,k) top Neighbour
        */
        for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        inv_Area = Cells_Inv_Area[Cell_No]; // 2D cell area or volume of the cell
        convective_u = 0.0;
        convective_v = 0.0;
        for (int Face_No = 0; Face_No < 4; Face_No++)
        {
            index = Face_No * 2;
            nx = Cell_Face_Normals[Cell_No][index + 0]; // Normal in x direction
            ny = Cell_Face_Normals[Cell_No][index + 1]; // normal in y direction
            A_f = Cell_Face_Areas[Cell_No][Face_No];    // area of the face in 2D the length of the face

            Neighbour = Cell_Neighbours[Cell_No][Face_No + 1]; // Get the neighbour corresponding to the face of a 2D cell

            // Convective terms
            uf = 0.5 * (Primitive_Cells[Cell_No][1] + Primitive_Cells[Neighbor][1]); // uf - average u velocity on face
            vf = 0.5 * (Primitive_Cells[Cell_No][2] + Primitive_Cells[Neighbor][2]); // vf - avarage v velocity on face

            convective_u += (uf * nx * A_f) * uf; // uf*(uf.n)dA
            convective_v += (vf * ny * A_f) * vf; // vf*(vf.n)dA

            Calculate_Gradient_On_Face(Cell_No, 1, Face_No);
            diffusive_u += mu * (u_gradient[0] * nx + u_gradient[1] * ny) * A_f;
            Calculate_Gradient_On_Face(Cell_No, 2, Face_No);
            diffusive_v += mu * (v_gradient[0] * nx + v_gradient[1] * ny) * A_f;
        }

        rho = Primitive_Cells[Cell_no][0];
        // Update intermediate velocities
        V_Star[Cell_No][0] = Primitive_Cells[Cell_No][1] + dt / rho * (-convective_u + diffusive_u + fb[Cell_No]);
        V_Star[Cell_No][1] = Primitive_Cells[Cell_No][2] + dt / rho * (-convective_v + diffusive_v + fb[Cell_No]);
        // Pressure
    }
}

void Velocity_Corrector_Step()
{
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No)
    {
        for (int Face_no = 0; Face_no < 4; Face_no++)
        {
            int index = Face_No * 2;                           // Face normal indices
            double nx = Cell_Face_Normals[Cell_No][index];     // x-component of face normal
            double ny = Cell_Face_Normals[Cell_No][index + 1]; // y-component of face normal
            double A_f = Cell_Face_Areas[Cell_No][Face_No];    // Length of the face

            // Velocity correction from this face
            Calculate_Gradient_On_Face(Cell_No, 4, Face_No);
            correction_x += P_Gradient * nx * A_f;
            correction_y += P_Gradient * ny * A_f;
        }
        // Normalize by cell volume and apply to update the velocity
        Primitive_Cells[Cell_No][1] = V_Star[Cell_No][0] - dt * correction_x / (rho * Cells_Area[Cell_No]);
        Primitive_Cells[Cell_No][2] = V_Star[Cell_No][1] - dt * correction_y / (rho * Cells_Area[Cell_No]);
    }
}

void solvePressurePoisson(
    int numCells,
    int numFacesPerCell,
    const std::vector<std::vector<int>> &cellNeighbors,  // Neighbor indices [cell][face]
    const std::vector<std::vector<double>> &faceNormals, // Face normals [cell][face]
    const std::vector<std::vector<double>> &faceAreas,   // Face areas [cell][face]
    const std::vector<double> &cellVolumes,              // Volumes of the cells
    const std::vector<std::vector<double>> &uStar,       // Intermediate velocity [cell][components]
    double rho,                                          // Density
    double dt,                                           // Time step
    std::vector<double> &pressure,                       // Output: Pressure at cell centers
    SparseMatrix &A,                                     // Matrix for linear system
    std::vector<double> &b                               // RHS vector
)
{
    // Assemble the coefficient matrix and RHS vector
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        double diag = 0.0; // Diagonal entry for the current cell
        double rhs = 0.0;  // RHS entry for the current cell

        for (int Face_No = 0; Face_No < NumFacesPerCell[Cell_No]; Face_No++)
        {
            int neighbor = cellNeighbors[Cell_No][face];       // Neighbor cell index
            double nx = Cell_Face_Normals[Cell_No][index];     // x-component of face normal
            double ny = Cell_Face_Normals[Cell_No][index + 1]; // y-component of face normal
            double A_f = Cell_Face_Areas[Cell_No][Face_No];    // Face area
            double d_f = Distance_Bw_Cell_Centers[Face_No];    // Placeholder: Distance between cell centers (compute appropriately)

            // Contribution to the RHS from velocity divergence
            double u_star_f = (uStar[cell][0] * nx + uStar[cell][1] * ny);
            rhs += u_star_f * A_f;

            // Contribution to the matrix
            if (neighbor >= 0)
            { // Internal face
                double coeff = A_f / d_f;
                diag += coeff;                      // Accumulate diagonal term
                A.addEntry(cell, neighbor, -coeff); // Off-diagonal term
            }
            else
            {
                // Boundary face: Apply Neumann or Dirichlet conditions
                // Example: Neumann (zero normal pressure gradient)
                // No contribution to RHS or matrix
            }
        }

        // Add diagonal entry and normalize RHS
        A.addEntry(cell, cell, diag);
        b[cell] = (rho / dt) * rhs / cellVolumes[cell];
    }

    // Solve the linear system A * p = b using an iterative solver
    solveLinearSystem(A, pressure, b);
}

void solvePressurePoisson(int &maxIters)
{
    const int Grad_Type = 4;
    double div_u = 0.0 for (int iter = 0; iter < maxIters; iter++)
    {
        for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
        {
            double rhs = 0.0, lhs = 0.0;
            // Loop over faces
            for (int Face_No = 0; Face_No < NumFacesPerCell[Cell_No]; Face_No++)
            {
                index = Face_No * 2;
                nx = Cell_Face_Normals[Cell_No][index + 0]; // Normal in x direction
                ny = Cell_Face_Normals[Cell_No][index + 1]; // normal in y direction
                dl = Cell_Face_Areas[Cell_No][Face_No];     // area of the face in 2D the length of the face

                Neighbour = Cell_Neighbours[Cell_No][Face_No + 1]; // Get the neighbour corresponding to the face of a 2D cell
                df = Distance_Bw_Cell_Centers[Face_No];            // distance between the cell center and the neighbouring cell center

                    div_u += (0.5*(u_star[Cell_No] + u_star[Neighbour]) * nx * dl;
					div_v += (0.5*(v_star[Cell_No] + v_star[Neighbour]) * ny * dl;
		    
					rhs += (div_u + div_v) / (dt * rho);

					double grad_p = (p_Nf[i][j][f] - p_P[i][j]) / df;
					lhs += grad_p * faceAreas[f];
            }
            // Update pressure
            P_new[Cell_No] = (rhs - lhs) / Cells_Area[Cell_No];
        }
    }
}

void computeCellCenteredVelocities()
{
    for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
    {
        double u_x = 0.0, u_y = 0.0; // Cell-centered velocity components

        for (int Face_No = 0; Face_no < NumFacesPerCell[Cell_No]; Face_No++)
        {
            int index = Face_No * 2;
            nx = Cell_Face_Normals[Cell_No][index + 0]; // Normal in x direction
            ny = Cell_Face_Normals[Cell_No][index + 1]; // normal in y direction
            A_f = Cell_Face_Areas[Cell_No][Face_No];    // area of the face in 2D the length of the face

            double u_f_x = FaceVelocities[Cell_No][index];     // x-component of velocity on face
            double u_f_y = FaceVelocities[Cell_No][index + 1]; // y-component of velocity on face

            // Accumulate contributions from face velocities
            u_x += (u_f_x * nx + u_f_y * ny) * A_f;
            u_y += (u_f_y * ny + u_f_x * nx) * A_f;
        }

        // Normalize by cell volume
        double V_P = Cells_Inv_Area[Cell_No];    // 1/ Area of the Cell
        Primitive_Cells[Cell_No][1] = u_x * V_P; // x-component of cell-centered velocity
        Primitive_Cells[Cell_No][2] = u_y * V_P; // y-component of cell-centered velocity
    }
}
