#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Parameters for the grid
const int NX = 100;    // Grid points in ξ direction
const int NY = 100;    // Grid points in η direction
const double Lx = 2.0; // Outer boundary distance from the cylinder center
const double cylinderRadius = 0.1; // Cylinder radius
const double tolerance = 1e-6; // Convergence tolerance for iterative solver

// Computational grid
std::vector<std::vector<double> > x(NX, std::vector<double>(NY, 0.0));
std::vector<std::vector<double> > y(NX, std::vector<double>(NY, 0.0));

// Initialize computational grid and boundary conditions
void initializeGrid() {
    double dTheta = 2.0 * M_PI / (NX - 1);
    double dR = (Lx - cylinderRadius) / (NY - 1);

    // Boundary conditions
    for (int i = 0; i < NX; ++i) {
        double theta = i * dTheta;
        // Cylinder surface
        x[i][0] = cylinderRadius * cos(theta);
        y[i][0] = cylinderRadius * sin(theta);
        // Outer boundary
        x[i][NY - 1] = Lx * cos(theta);
        y[i][NY - 1] = Lx * sin(theta);
    }

    for (int j = 0; j < NY; ++j) {
        double r = cylinderRadius + j * dR;
        // Interpolate between the inner and outer boundary
        x[0][j] = r;          y[0][j] = 0.0; // Right side
        x[NX - 1][j] = -r;     y[NX - 1][j] = 0.0; // Left side
    }
}

// Laplace smoothing solver
void solveEllipticGrid() {
    double maxError = 0.0;

    do {
        maxError = 0.0;

        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                double x_old = x[i][j];
                double y_old = y[i][j];

                // Laplace equation for smoothing
                x[i][j] = 0.25 * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1]);
                y[i][j] = 0.25 * (y[i + 1][j] + y[i - 1][j] + y[i][j + 1] + y[i][j - 1]);

                double error = std::max(std::abs(x[i][j] - x_old), std::abs(y[i][j] - y_old));
                maxError = std::max(maxError, error);
            }
        }

        // Print progress
        std::cout << "Current max error: " << maxError << std::endl;

    } while (maxError > tolerance);
}

// Output grid to VTK file for visualization
void writeVTK(const std::string &filename) {
    std::ofstream vtkFile(filename);
    vtkFile << "# vtk DataFile Version 2.0\n";
    vtkFile << "Elliptic grid for flow over a cylinder\n";
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
    writeVTK("elliptic_grid.vtk");
    std::cout << "Elliptic grid generated and output written to elliptic_grid.vtk." << std::endl;
    return 0;
}