#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Parameters for the grid and domain
const int NX = 300;        // Grid points in X direction
const int NY = 100;        // Grid points in Y direction
const double Lx = 3.0;     // Half-width of the rectangular domain (X direction)
const double Ly = 1.0;     // Half-height of the rectangular domain (Y direction)
const double cylinderRadius = 0.1;  // Radius of the embedded cylinder
const double cylinderCenterX = 0.0; // X coordinate of the cylinder center
const double cylinderCenterY = 0.0; // Y coordinate of the cylinder center
const double tolerance = 1e-6;      // Convergence tolerance
const int maxIterations = 10000;    // Maximum number of iterations

// Computational grid
std::vector<std::vector<double> > x(NX, std::vector<double>(NY, 0.0));
std::vector<std::vector<double> > y(NX, std::vector<double>(NY, 0.0));

void initializeGrid() {
    // Boundary conditions for the outer rectangular boundary
    for (int i = 0; i < NX; ++i) {
        double xi = (double)i / (NX - 1);
        x[i][0] = -Lx + 2 * Lx * xi;    // Bottom boundary (y = -Ly)
        y[i][0] = -Ly;

        x[i][NY - 1] = -Lx + 2 * Lx * xi; // Top boundary (y = +Ly)
        y[i][NY - 1] = Ly;
    }

    for (int j = 0; j < NY; ++j) {
        double eta = (double)j / (NY - 1);
        x[0][j] = -Lx;              // Left boundary (x = -Lx)
        y[0][j] = -Ly + 2 * Ly * eta;

        x[NX - 1][j] = Lx;          // Right boundary (x = +Lx)
        y[NX - 1][j] = -Ly + 2 * Ly * eta;
    }

    // Boundary conditions for the cylinder
    int numCylinderPoints = NX - 2; // Number of points to place around the cylinder
    for (int i = 0; i < numCylinderPoints; ++i) {
        double theta = 2.0 * M_PI * i / numCylinderPoints; // Theta from 0 to 2π
        x[i][1] = cylinderCenterX + cylinderRadius * cos(theta); // Circle around cylinder
        y[i][1] = cylinderCenterY + cylinderRadius * sin(theta);
    }

    // Ensure closure of the cylinder by duplicating the first point at the end
    x[numCylinderPoints][1] = x[0][1];
    y[numCylinderPoints][1] = y[0][1];

    // Initial guess for internal points (linear interpolation)
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            double xi = (double)i / (NX - 1);
            double eta = (double)j / (NY - 1);
            x[i][j] = (1 - eta) * x[i][0] + eta * x[i][NY - 1];  // Interpolate vertically
            y[i][j] = (1 - eta) * y[i][0] + eta * y[i][NY - 1];
        }
    }
}

// Laplace smoothing solver to adjust interior points
void solveEllipticGrid() {
    double maxError;
    int iteration = 0;

    do {
        maxError = 0.0;

        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 2; j < NY - 2; ++j) { // Skip the cylinder boundary at j = 1
                double x_old = x[i][j];
                double y_old = y[i][j];

                // Laplace equation for smoothing (interior points only)
                x[i][j] = 0.25 * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1]);
                y[i][j] = 0.25 * (y[i + 1][j] + y[i - 1][j] + y[i][j + 1] + y[i][j - 1]);

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
    vtkFile << "Elliptic grid for rectangular domain with embedded cylinder\n";
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
    writeVTK("Initial_rectangular_cylinder.vtk");    
    solveEllipticGrid();
    writeVTK("elliptic_rectangular_cylinder.vtk");
    std::cout << "Elliptic grid for rectangular domain with embedded cylinder generated and output written to elliptic_rectangular_cylinder.vtk." << std::endl;
    return 0;
}