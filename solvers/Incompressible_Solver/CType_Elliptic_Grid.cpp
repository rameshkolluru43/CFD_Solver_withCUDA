#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Parameters for the grid
const int NX = 100;    // Grid points in ξ direction (around the cylinder)
const int NY = 100;    // Grid points in η direction (from cylinder to far field)
const double Lx = 10.0; // Extended far-field distance for smoother downstream grid
const double cylinderRadius = 0.1; // Cylinder radius
const double tolerance = 1e-7; // Tighter convergence tolerance
const int maxIterations = 5000; // Maximum number of iterations for relaxation

// Computational grid
std::vector<std::vector<double> > x(NX, std::vector<double>(NY, 0.0));
std::vector<std::vector<double> > y(NX, std::vector<double>(NY, 0.0));

// Control functions for Poisson-type smoothing
double P(int i, int j) {
    double xi = (double)i / (NX - 1); // Normalized ξ direction (0 to 1)
    double eta = (double)j / (NY - 1); // Normalized η direction (0 to 1)

    // Control function for downstream stretching
    // Increasing with η (more influence in downstream)
    double downstream_influence = 5.0 * eta * eta; // Quadratic growth in downstream

    // Control function for near-cylinder refinement
    // Higher near cylinder (closer to ξ=0.5)
    double near_cylinder_influence = 1.0 / (1.0 + 100.0 * (xi - 0.5) * (xi - 0.5));

    return downstream_influence + near_cylinder_influence;
}

double Q(int i, int j) {
    double xi = (double)i / (NX - 1); // Normalized ξ direction (0 to 1)
    double eta = (double)j / (NY - 1); // Normalized η direction (0 to 1)

    // Control function for downstream stretching
    double downstream_influence = 5.0 * eta * eta; // Quadratic growth in downstream

    // Control function for near-cylinder refinement
    double near_cylinder_influence = 1.0 / (1.0 + 100.0 * (xi - 0.5) * (xi - 0.5));

    return downstream_influence + near_cylinder_influence;
}

// Initialize computational grid and boundary conditions
void initializeGrid() {
    double dTheta = M_PI / (NX - 1);      // Half-circle around the cylinder
    double dR = (Lx - cylinderRadius) / (NY - 1); // Radial spacing

    // Wrap around the cylinder from θ = -π to θ = π
    for (int i = 0; i < NX; ++i) {
        double theta = -M_PI + i * dTheta;
        x[i][0] = cylinderRadius * cos(theta); // Inner boundary at cylinder surface
        y[i][0] = cylinderRadius * sin(theta);

        // Far-field boundary (more spaced outward)
        x[i][NY - 1] = Lx * cos(theta);
        y[i][NY - 1] = Lx * sin(theta);
    }

    // Initialize downstream boundary with smoother spacing
    for (int j = 0; j < NY; ++j) {
        double r = cylinderRadius + j * dR;
        x[0][j] = -Lx; y[0][j] = r;       // Left side of the grid, downstream
        x[NX - 1][j] = Lx; y[NX - 1][j] = r; // Right side of the grid, downstream
    }
}

// Poisson smoothing solver with control functions
void solveEllipticGrid() {
    double maxError;
    int iteration = 0;

    do {
        maxError = 0.0;

        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                double x_old = x[i][j];
                double y_old = y[i][j];

                // Poisson equation with control functions for smoother grid distribution
                x[i][j] = 0.25 * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1] - P(i, j));
                y[i][j] = 0.25 * (y[i + 1][j] + y[i - 1][j] + y[i][j + 1] + y[i][j - 1] - Q(i, j));

                double error = std::max(std::abs(x[i][j] - x_old), std::abs(y[i][j] - y_old));
                maxError = std::max(maxError, error);
            }
        }

        // Increment iteration count and print progress every 100 iterations
        iteration++;
        if (iteration % 100 == 0) {
            std::cout << "Iteration: " << iteration << ", Max Error: " << maxError << std::endl;
        }

    } while (maxError > tolerance && iteration < maxIterations);

    std::cout << "Final iteration: " << iteration << ", Final Max Error: " << maxError << std::endl;
}

// Output grid to VTK file for visualization
void writeVTK(const std::string &filename) {
    std::ofstream vtkFile(filename);
    vtkFile << "# vtk DataFile Version 2.0\n";
    vtkFile << "Improved Elliptic C-type grid with control functions\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << NX << " " << NY << " 1\n";
    vtkFile << "POINTS " << NX * NY << " float\n";

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            vtkFile << x[i][j] << " " << y[i][j] << " 0\n";
        }
    }

    vtkFile.close();
}

int main() {
    initializeGrid();
    solveEllipticGrid();
    writeVTK("elliptic_c_grid.vtk");
    std::cout << "Improved C-type grid generated with control functions. Output written to elliptic_c_grid_control_functions.vtk." << std::endl;
    return 0;
}