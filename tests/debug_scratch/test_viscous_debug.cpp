#include <iostream>
#include <cmath>
using namespace std;

double T_S_Mu = 110.4;  // Sutherland constant for air
double T_ref = 288.0;   // Reference temperature
double Pr = 0.72;
double cp_ref = 1004.5;  // Specific heat for air

int main() {
    // Test with non-dimensional temperature from the simulation
    double T_Star = 3.9984;  // From viscous simulation
    
    double Term1 = T_S_Mu / T_ref;
    double Term2 = T_Star + Term1;
    double mu_star = pow(T_Star, 1.5) * ((1.0 + Term1) / Term2);
    
    cout << "Non-dimensional temperature: " << T_Star << endl;
    cout << "Term1 (T_S_Mu/T_ref): " << Term1 << endl;
    cout << "Term2 (T_Star + Term1): " << Term2 << endl;
    cout << "mu_star: " << mu_star << endl;
    
    // Check with Re and gradients
    double Re = 100000.0;
    double Inv_Re = 1.0 / Re;
    double du_dy = 4.0 / 0.02;  // Approximate wall gradient: delta_u / delta_y
    
    double tau = mu_star * Inv_Re * du_dy;
    cout << "\nInv_Re: " << Inv_Re << endl;
    cout << "Approx du/dy at wall: " << du_dy << endl;
    cout << "Approx shear stress: " << tau << endl;
    
    return 0;
}
