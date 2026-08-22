// Quick debug to understand the issue
#include <iostream>
#include <vector>
using namespace std;

int main() {
    // Simulate the gradient calculation
    // For uniform u=2 everywhere except ghost cell u=-2
    
    // Interior cell and its 8 neighbors - if all have u=2 except ghost
    // the gradient should be non-zero
    
    double cell_avg_interior = 2.0;  // u = 2
    double cell_avg_ghost = -2.0;    // u = -2 (reflected)
    
    // Simple 2D gradient on a face:
    // du/dn = (u_neighbor - u_current) / distance
    double distance = 0.025;  // typical grid spacing
    double du_dn = (cell_avg_ghost - cell_avg_interior) / distance;
    
    cout << "Simple gradient: du/dn = " << du_dn << endl;
    
    // With viscosity and Re
    double mu = 2.5;  // approx from earlier calculation
    double Re = 100000.0;
    double tau = mu / Re * du_dn;
    
    cout << "Shear stress tau = mu/Re * du/dn = " << tau << endl;
    
    // Now the issue: if the gradient calculation averages over 8 neighbors
    // and only 1 is the ghost cell, the average will be dominated by
    // the other 7 cells which all have u=2
    
    // For Face 0 gradient, using Calculate_Vertex_Average:
    // o = average of {current, neighbor1, neighbor2, neighbor5}
    // If all have u=2, o = 2
    // c = average of {current, neighbor1, neighbor4, neighbor8}
    // If all have u=2, c = 2
    
    // Then av[0] = 0.5*(current + o) = 0.5*(2 + 2) = 2
    // av[1] = 0.5*(current + c) = 0.5*(2 + 2) = 2
    // av[2] = 0.5*(neighbor1 + c) = 0.5*(2 + 2) = 2
    // av[3] = 0.5*(neighbor1 + o) = 0.5*(2 + 2) = 2
    
    // The gradient would be:
    // grad_x = sum(av[i] * nx[i] * dl[i]) / area
    // grad_y = sum(av[i] * ny[i] * dl[i]) / area
    
    // If all av[i] = 2, and the co-volume is symmetric,
    // the gradient components would cancel out to ~0!
    
    cout << "\nThe issue: With Green-Gauss gradient on a co-volume," << endl;
    cout << "if the 'av' values are nearly equal, the gradient is near zero." << endl;
    cout << "The wall BC sets ghost cell primitive, but the gradient calculation" << endl;
    cout << "may not be using the ghost cell directly in 'neighbors' array." << endl;
    
    return 0;
}
