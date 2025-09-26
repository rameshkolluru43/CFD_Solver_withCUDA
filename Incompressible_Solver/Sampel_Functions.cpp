void solveMomentumImplicit(
    double** u, double** v, double** u_old, double** v_old,
    double** p, int nx, int ny, double dx, double dy, double dt, double rho, double mu
) {
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            // Implicit terms: diffusion and pressure
            double diffusion_u = mu * (u_old[i+1][j] - 2 * u_old[i][j] + u_old[i-1][j]) / (dx * dx)
                               + mu * (u_old[i][j+1] - 2 * u_old[i][j] + u_old[i][j-1]) / (dy * dy);

            double pressure_u = (p[i][j] - p[i-1][j]) / dx;

            // Update u-velocity using implicit formulation
            u[i][j] = (rho / dt * u_old[i][j] + diffusion_u - pressure_u) / (rho / dt);
        }
    }

    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            // Implicit terms: diffusion and pressure for v
            double diffusion_v = mu * (v_old[i+1][j] - 2 * v_old[i][j] + v_old[i-1][j]) / (dx * dx)
                               + mu * (v_old[i][j+1] - 2 * v_old[i][j] + v_old[i][j-1]) / (dy * dy);

            double pressure_v = (p[i][j] - p[i][j-1]) / dy;

            // Update v-velocity using implicit formulation
            v[i][j] = (rho / dt * v_old[i][j] + diffusion_v - pressure_v) / (rho / dt);
        }
    }
}

void computeIntermediateVelocity(
    int numCells,
    const std::vector<std::vector<double>>& Cell_Face_Normals,
    const std::vector<double>& Cell_Face_Areas,
    const std::vector<int>& Cell_Neighbours,
    const std::vector<double>& Distance_Bw_Cell_Centers,
    const std::vector<std::vector<double>>& primitive_cells,
    std::vector<std::vector<double>>& u_star, // Output: Intermediate velocity (u, v)
    double nu, // Kinematic viscosity
    double dt, // Time step
    const std::vector<double>& body_forces // External body forces, e.g., gravity
) {
    for (int cell_no = 0; cell_no < numCells; ++cell_no) {
        double conv_u_x = 0.0, conv_u_y = 0.0; // Convective terms
        double diff_u_x = 0.0, diff_u_y = 0.0; // Diffusive terms

        for (int face_no = 0; face_no < 4; ++face_no) {
            int index = face_no * 2; // Face normal indices
            double nx = Cell_Face_Normals[cell_no][index];     // x-component of face normal
            double ny = Cell_Face_Normals[cell_no][index + 1]; // y-component of face normal
            double dl = Cell_Face_Areas[cell_no][face_no];     // Length of the face
            int neighbor = Cell_Neighbours[cell_no][face_no + 1]; // Neighbor corresponding to the face
            double df = Distance_Bw_Cell_Centers[face_no];     // Distance between cell centroids

            // Velocities at the current cell and neighboring cell
            double u_P_x = primitive_cells[cell_no][1]; // u_x in current cell
            double u_P_y = primitive_cells[cell_no][2]; // u_y in current cell
            double u_N_x = (neighbor >= 0) ? primitive_cells[neighbor][1] : u_P_x; // u_x in neighbor
            double u_N_y = (neighbor >= 0) ? primitive_cells[neighbor][2] : u_P_y; // u_y in neighbor

            // Convective fluxes
            double u_f_x = 0.5 * (u_P_x + u_N_x); // Interpolated u_x at face
            double u_f_y = 0.5 * (u_P_y + u_N_y); // Interpolated u_y at face
            conv_u_x += (u_f_x * nx + u_f_y * ny) * u_f_x * dl; // Convective term for u_x
            conv_u_y += (u_f_x * nx + u_f_y * ny) * u_f_y * dl; // Convective term for u_y

            // Diffusive fluxes
            diff_u_x += nu * ((u_N_x - u_P_x) / df) * dl; // Diffusive term for u_x
            diff_u_y += nu * ((u_N_y - u_P_y) / df) * dl; // Diffusive term for u_y
        }

        // Update intermediate velocity
        double volume = Cell_Face_Areas[cell_no][4]; // Assume area (volume in 2D) stored at index 4
        u_star[cell_no][0] = primitive_cells[cell_no][1] + dt * (-conv_u_x + diff_u_x + body_forces[0] * volume);
        u_star[cell_no][1] = primitive_cells[cell_no][2] + dt * (-conv_u_y + diff_u_y + body_forces[1] * volume);
    }
}


void solvePressurePoisson(
    int numCells,
    const std::vector<std::vector<double>>& Cell_Face_Normals,
    const std::vector<double>& Cell_Face_Areas,
    const std::vector<int>& Cell_Neighbours,
    const std::vector<double>& Distance_Bw_Cell_Centers,
    const std::vector<std::vector<double>>& u_star,
    std::vector<std::vector<double>>& primitive_cells, // Input/Output: Contains pressure values
    SparseMatrix& A, // Sparse matrix for the pressure Poisson system
    std::vector<double>& b, // Right-hand side of the system
    double rho, // Fluid density
    double dt // Time step
) {
    // Assemble the system
    for (int cell_no = 0; cell_no < numCells; ++cell_no) {
        double diag = 0.0;
        double rhs = 0.0;

        for (int face_no = 0; face_no < 4; ++face_no) {
            int index = face_no * 2; // Face normal indices
            double nx = Cell_Face_Normals[cell_no][index];     // x-component of face normal
            double ny = Cell_Face_Normals[cell_no][index + 1]; // y-component of face normal
            double dl = Cell_Face_Areas[cell_no][face_no];     // Length of the face
            int neighbor = Cell_Neighbours[cell_no][face_no + 1]; // Neighbor corresponding to the face
            double df = Distance_Bw_Cell_Centers[face_no];     // Distance between cell centroids

            // Intermediate velocity flux
            double u_f = u_star[cell_no][0] * nx + u_star[cell_no][1] * ny; // Flux normal to face
            rhs += (u_f * dl) / dt;

            // Pressure gradient contribution
            if (neighbor >= 0) {
                double coeff = dl / df;
                diag += coeff;
                A.addEntry(cell_no, neighbor, -coeff);
            }
        }

        // Add diagonal and RHS
        A.addEntry(cell_no, cell_no, diag);
        b[cell_no] = rhs * rho;
    }

    // Solve the linear system (e.g., Conjugate Gradient)
    iterativeSolver(A, primitive_cells, b);
}



// Pressure Poisson solver
void SolvePressurePoisson()
{
    double residual = std::numeric_limits<double>::max();
    int iteration = 0;

    while (residual > tolerance && iteration < maxIterations)
    {
        residual = 0.0;
        // Loop over cells
		for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
		 {
			 double sumFlux = 0.0;
			 double normalization = 0.0;

			 // Loop over faces connected to the cell
			 for (int Face_No =0;Face_No<4;Face_No++)
			 {
				index = Face_No*2;
				nx = Cell_Face_Normals[Cell_No][index + 0]; // Normal in x direction 
				ny = Cell_Face_Normals[Cell_No][index + 1]; // normal in y direction 
				dl = Cell_Face_Areas[Cell_No][Face_No]; // area of the face in 2D the length of the face 

				Neighbour = Cell_Neighbours[Cell_No][Face_No + 1]; // Get the neighbour corresponding to the face of a 2D cell 
				dPN = Distance_Bw_Cell_Centers[Face_No]; // distance between the cell center and the neighbouring cell center 
				P_N = Primitive_Cells[Neighbour][4]; // Pressure at the Neighbour Cell
				
				SumFlux += (dl / dPN) * P_N;
				Normalization += dl/dPN;
				
				div_u += (0.5*(u_star[Cell_No] + u_star[Neighbour]) * nx * dl;
				div_v += (0.5*(v_star[Cell_No] + v_star[Neighbour]) * ny * dl;
	    
				RHS += (div_u + div_v) / (dt * rho); 
               	 	} 
            	}

            // Update pressure
            double P_Old = Primitive_Cells[Neighbour][4]; ;
            P_New = (SumFlux - cells[i].f * cells[i].volume) / normalization;

            // Compute residual
            residual = std::max(residual, std::fabs(cells[i].p - oldPressure));
        }

        iteration++;
    }
