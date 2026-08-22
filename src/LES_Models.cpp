/**
 * @file LES_Models.cpp
 * @brief Implementation of Large Eddy Simulation (LES) subgrid-scale models
 * @author CFD Solver Team
 * @date 2026
 *
 * Implements SGS models: Smagorinsky, Dynamic Smagorinsky, WALE, Vreman, Sigma
 * for compressible flow simulations.
 */

#include "LES_Models.h"
#include "Viscous_Functions.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

//=============================================================================
// GLOBAL LES VARIABLES DEFINITIONS
// Note: Is_LES, Use_Wall_Damping, Use_GPU_LES, Cs_Smagorinsky, Cw_WALE,
// Cv_Vreman, Csigma, Pr_turb_LES are defined in Initialize.cpp
//=============================================================================

LES_SGS_Model current_les_model = LES_SGS_Model::NONE;
vector<LES_Variables> les_vars;

//=============================================================================
// INITIALIZATION FUNCTIONS
//=============================================================================

void Initialize_LES_Model(LES_SGS_Model model)
{
    current_les_model = model;
    Is_LES = (model != LES_SGS_Model::NONE);

    cout << "=== LES MODEL INITIALIZATION ===" << endl;

    switch (model)
    {
    case LES_SGS_Model::NONE:
        cout << "LES disabled (DNS or RANS mode)" << endl;
        return;

    case LES_SGS_Model::SMAGORINSKY:
        cout << "Using Smagorinsky SGS Model (Cs = " << Cs_Smagorinsky << ")" << endl;
        break;

    case LES_SGS_Model::DYNAMIC_SMAGORINSKY:
        cout << "Using Dynamic Smagorinsky SGS Model" << endl;
        break;

    case LES_SGS_Model::WALE:
        cout << "Using WALE SGS Model (Cw = " << Cw_WALE << ")" << endl;
        break;

    case LES_SGS_Model::VREMAN:
        cout << "Using Vreman SGS Model (Cv = " << Cv_Vreman << ")" << endl;
        break;

    case LES_SGS_Model::SIGMA:
        cout << "Using Sigma SGS Model (Csigma = " << Csigma << ")" << endl;
        break;

    default:
        cout << "Warning: Unknown LES model. Disabling LES." << endl;
        current_les_model = LES_SGS_Model::NONE;
        Is_LES = false;
        return;
    }

    Initialize_LES_Variables();
    Set_LES_Filter_Width();

    if (Use_Wall_Damping)
    {
        Calculate_Wall_Distances_LES_All();
        cout << "Wall damping enabled for near-wall treatment" << endl;
    }

    cout << "=== LES MODEL INITIALIZED ===" << endl;
}

void Initialize_LES_Variables()
{
    les_vars.clear();
    les_vars.resize(Total_No_Cells);

    for (int i = 0; i < Total_No_Cells; ++i)
    {
        les_vars[i] = LES_Variables();
    }

    cout << "LES variables initialized for " << Total_No_Cells << " cells" << endl;
}

void Set_LES_Filter_Width()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_LES_Filter_Width(i);
        }
    }
    cout << "LES filter widths computed for all cells" << endl;
}

void Calculate_LES_Filter_Width(int cell_index)
{
    les_vars[cell_index].delta = Calculate_Filter_Width_Cube_Root(cell_index);
}

//=============================================================================
// FILTER WIDTH CALCULATIONS
//=============================================================================

double Calculate_Filter_Width_Cube_Root(int cell_index)
{
    double area = Cells[cell_index].Area;
    return sqrt(area);
}

double Calculate_Filter_Width_Max_Edge(int cell_index)
{
    double max_edge = 0.0;
    int num_faces = static_cast<int>(Cells[cell_index].Face_Areas.size());
    
    for (int f = 0; f < num_faces; ++f)
    {
        double edge_length = Cells[cell_index].Face_Areas[f];
        max_edge = max(max_edge, edge_length);
    }
    return max_edge;
}

double Calculate_Filter_Width_Face_Area(int cell_index)
{
    double total_area = 0.0;
    int num_faces = static_cast<int>(Cells[cell_index].Face_Areas.size());
    
    for (int f = 0; f < num_faces; ++f)
    {
        total_area += Cells[cell_index].Face_Areas[f];
    }
    return total_area / num_faces;
}

//=============================================================================
// VELOCITY GRADIENT TENSOR
//=============================================================================

void Calculate_Velocity_Gradient_Tensor(int cell_index, double G[3][3])
{
    double dudx = 0.0, dudy = 0.0, dvdx = 0.0, dvdy = 0.0;

    int num_neighbors = static_cast<int>(Cells[cell_index].Neighbours.size());
    if (num_neighbors < 2)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                G[i][j] = 0.0;
        return;
    }

    double sum_w = 0.0;
    double u0 = Primitive_Cells[cell_index][1];
    double v0 = Primitive_Cells[cell_index][2];
    double x0 = Cells[cell_index].Cell_Center[0];
    double y0 = Cells[cell_index].Cell_Center[1];

    double sum_wx = 0.0, sum_wy = 0.0;
    double sum_wxx = 0.0, sum_wyy = 0.0, sum_wxy = 0.0;
    double sum_wudx = 0.0, sum_wudy = 0.0;
    double sum_wvdx = 0.0, sum_wvdy = 0.0;

    for (int n = 0; n < num_neighbors; ++n)
    {
        int neigh = Cells[cell_index].Neighbours[n];
        if (neigh < 0 || neigh >= Total_No_Cells)
            continue;

        double xn = Cells[neigh].Cell_Center[0];
        double yn = Cells[neigh].Cell_Center[1];
        double dx = xn - x0;
        double dy = yn - y0;
        double dist = sqrt(dx * dx + dy * dy);
        if (dist < 1e-14)
            continue;

        double w = 1.0 / (dist * dist);
        double un = Primitive_Cells[neigh][1];
        double vn = Primitive_Cells[neigh][2];
        double du = un - u0;
        double dv = vn - v0;

        sum_w += w;
        sum_wx += w * dx;
        sum_wy += w * dy;
        sum_wxx += w * dx * dx;
        sum_wyy += w * dy * dy;
        sum_wxy += w * dx * dy;
        sum_wudx += w * du * dx;
        sum_wudy += w * du * dy;
        sum_wvdx += w * dv * dx;
        sum_wvdy += w * dv * dy;
    }

    double det = sum_wxx * sum_wyy - sum_wxy * sum_wxy;
    if (fabs(det) > 1e-14)
    {
        dudx = (sum_wyy * sum_wudx - sum_wxy * sum_wudy) / det;
        dudy = (sum_wxx * sum_wudy - sum_wxy * sum_wudx) / det;
        dvdx = (sum_wyy * sum_wvdx - sum_wxy * sum_wvdy) / det;
        dvdy = (sum_wxx * sum_wvdy - sum_wxy * sum_wvdx) / det;
    }

    G[0][0] = dudx;
    G[0][1] = dudy;
    G[0][2] = 0.0;
    G[1][0] = dvdx;
    G[1][1] = dvdy;
    G[1][2] = 0.0;
    G[2][0] = 0.0;
    G[2][1] = 0.0;
    G[2][2] = 0.0;
}

void Calculate_Strain_Rate_Tensor_LES(int cell_index, double S[3][3])
{
    double G[3][3];
    Calculate_Velocity_Gradient_Tensor(cell_index, G);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            S[i][j] = 0.5 * (G[i][j] + G[j][i]);
        }
    }
}

void Calculate_Rotation_Rate_Tensor(int cell_index, double Omega[3][3])
{
    double G[3][3];
    Calculate_Velocity_Gradient_Tensor(cell_index, G);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Omega[i][j] = 0.5 * (G[i][j] - G[j][i]);
        }
    }
}

double Calculate_Strain_Rate_Magnitude_LES(const double S[3][3])
{
    double S_sq = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            S_sq += S[i][j] * S[i][j];
        }
    }
    return sqrt(2.0 * S_sq);
}

double Calculate_Vorticity_Magnitude_LES(const double Omega[3][3])
{
    double Omega_sq = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Omega_sq += Omega[i][j] * Omega[i][j];
        }
    }
    return sqrt(2.0 * Omega_sq);
}

//=============================================================================
// SMAGORINSKY MODEL
//=============================================================================

double Smagorinsky_Nu_SGS(double delta, double S_mag, double Cs)
{
    return Cs * Cs * delta * delta * S_mag;
}

void Calculate_Smagorinsky_Viscosity(int cell_index)
{
    double S[3][3];
    Calculate_Strain_Rate_Tensor_LES(cell_index, S);

    double S_mag = Calculate_Strain_Rate_Magnitude_LES(S);
    les_vars[cell_index].S_mag = S_mag;

    double delta = les_vars[cell_index].delta;
    double Cs = Cs_Smagorinsky;

    if (Use_Wall_Damping && les_vars[cell_index].y_plus < 100.0)
    {
        Cs *= les_vars[cell_index].f_damping;
    }

    double nu_sgs = Smagorinsky_Nu_SGS(delta, S_mag, Cs);
    les_vars[cell_index].nu_sgs = nu_sgs;

    double rho = Primitive_Cells[cell_index][0];
    if (rho < 1e-14)
        rho = 1e-14;
    les_vars[cell_index].mu_sgs = rho * nu_sgs;
}

void Calculate_Smagorinsky_Viscosity_All()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_Smagorinsky_Viscosity(i);
        }
    }
}

//=============================================================================
// DYNAMIC SMAGORINSKY MODEL
//=============================================================================

void Calculate_Leonard_Stress(int cell_index)
{
    int num_neighbors = static_cast<int>(Cells[cell_index].Neighbours.size());
    if (num_neighbors < 1)
        return;

    double u0 = Primitive_Cells[cell_index][1];
    double v0 = Primitive_Cells[cell_index][2];

    double u_bar = u0, v_bar = v0;
    double uu_bar = u0 * u0, vv_bar = v0 * v0, uv_bar = u0 * v0;
    double u_sq_bar = u0, v_sq_bar = v0;
    double weight_sum = 1.0;

    for (int n = 0; n < num_neighbors; ++n)
    {
        int neigh = Cells[cell_index].Neighbours[n];
        if (neigh < 0 || neigh >= Total_No_Cells)
            continue;

        double un = Primitive_Cells[neigh][1];
        double vn = Primitive_Cells[neigh][2];
        double w = 1.0;

        u_bar += w * un;
        v_bar += w * vn;
        uu_bar += w * un * un;
        vv_bar += w * vn * vn;
        uv_bar += w * un * vn;
        weight_sum += w;
    }

    u_bar /= weight_sum;
    v_bar /= weight_sum;
    uu_bar /= weight_sum;
    vv_bar /= weight_sum;
    uv_bar /= weight_sum;

    les_vars[cell_index].Lij[0][0] = uu_bar - u_bar * u_bar;
    les_vars[cell_index].Lij[0][1] = uv_bar - u_bar * v_bar;
    les_vars[cell_index].Lij[1][0] = les_vars[cell_index].Lij[0][1];
    les_vars[cell_index].Lij[1][1] = vv_bar - v_bar * v_bar;
}

void Calculate_Resolved_Stress(int cell_index)
{
    double S[3][3];
    Calculate_Strain_Rate_Tensor_LES(cell_index, S);
    double S_mag = Calculate_Strain_Rate_Magnitude_LES(S);

    double delta = les_vars[cell_index].delta;
    double delta_test = LES_Constants::test_filter_ratio * delta;

    double alpha_sq = delta * delta;
    double beta_sq = delta_test * delta_test - delta * delta;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            les_vars[cell_index].Mij[i][j] = 2.0 * beta_sq * S_mag * S[i][j];
        }
    }
}

double Calculate_Dynamic_Cs(int cell_index)
{
    Calculate_Leonard_Stress(cell_index);
    Calculate_Resolved_Stress(cell_index);

    double LijMij = 0.0, MijMij = 0.0;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            LijMij += les_vars[cell_index].Lij[i][j] * les_vars[cell_index].Mij[i][j];
            MijMij += les_vars[cell_index].Mij[i][j] * les_vars[cell_index].Mij[i][j];
        }
    }

    double Cs_sq = 0.0;
    if (fabs(MijMij) > 1e-14)
    {
        Cs_sq = LijMij / MijMij;
    }

    Cs_sq = max(LES_Constants::Cs_min * LES_Constants::Cs_min, 
                min(Cs_sq, LES_Constants::Cs_max * LES_Constants::Cs_max));

    return sqrt(max(0.0, Cs_sq));
}

void Calculate_Dynamic_Smagorinsky_Viscosity(int cell_index)
{
    double S[3][3];
    Calculate_Strain_Rate_Tensor_LES(cell_index, S);
    double S_mag = Calculate_Strain_Rate_Magnitude_LES(S);
    les_vars[cell_index].S_mag = S_mag;

    double Cs = Calculate_Dynamic_Cs(cell_index);
    les_vars[cell_index].Cs_dynamic = Cs;

    double delta = les_vars[cell_index].delta;

    if (Use_Wall_Damping && les_vars[cell_index].y_plus < 100.0)
    {
        Cs *= les_vars[cell_index].f_damping;
    }

    double nu_sgs = Cs * Cs * delta * delta * S_mag;
    les_vars[cell_index].nu_sgs = nu_sgs;

    double rho = Primitive_Cells[cell_index][0];
    if (rho < 1e-14)
        rho = 1e-14;
    les_vars[cell_index].mu_sgs = rho * nu_sgs;
}

void Calculate_Dynamic_Smagorinsky_All()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_Dynamic_Smagorinsky_Viscosity(i);
        }
    }
}

//=============================================================================
// WALE MODEL
//=============================================================================

void Calculate_Sd_Tensor(int cell_index)
{
    double G[3][3];
    Calculate_Velocity_Gradient_Tensor(cell_index, G);

    double g2[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                g2[i][j] += G[i][k] * G[k][j];
            }
        }
    }

    double trace_g2 = g2[0][0] + g2[1][1] + g2[2][2];

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            les_vars[cell_index].Sd_ij[i][j] = 0.5 * (g2[i][j] + g2[j][i]);
            if (i == j)
            {
                les_vars[cell_index].Sd_ij[i][j] -= trace_g2 / 3.0;
            }
        }
    }

    double Sd_sq = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Sd_sq += les_vars[cell_index].Sd_ij[i][j] * les_vars[cell_index].Sd_ij[i][j];
        }
    }
    les_vars[cell_index].Sd_mag = sqrt(Sd_sq);
}

double WALE_Nu_SGS(double delta, double S_mag, double Sd_mag, double Cw)
{
    double Sd_52 = pow(Sd_mag, 2.5);
    double S_52 = pow(S_mag, 2.5);
    double denom = pow(S_mag, 5.0) + Sd_52 * Sd_52;

    if (denom < 1e-20)
        return 0.0;

    return Cw * Cw * delta * delta * Sd_52 * Sd_52 / sqrt(denom);
}

void Calculate_WALE_Viscosity(int cell_index)
{
    double S[3][3];
    Calculate_Strain_Rate_Tensor_LES(cell_index, S);
    double S_mag = Calculate_Strain_Rate_Magnitude_LES(S);
    les_vars[cell_index].S_mag = S_mag;

    Calculate_Sd_Tensor(cell_index);
    double Sd_mag = les_vars[cell_index].Sd_mag;

    double delta = les_vars[cell_index].delta;
    double nu_sgs = WALE_Nu_SGS(delta, S_mag, Sd_mag, Cw_WALE);

    les_vars[cell_index].nu_sgs = nu_sgs;

    double rho = Primitive_Cells[cell_index][0];
    if (rho < 1e-14)
        rho = 1e-14;
    les_vars[cell_index].mu_sgs = rho * nu_sgs;
}

void Calculate_WALE_Viscosity_All()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_WALE_Viscosity(i);
        }
    }
}

//=============================================================================
// VREMAN MODEL
//=============================================================================

void Calculate_Vreman_B(int cell_index)
{
    double G[3][3];
    Calculate_Velocity_Gradient_Tensor(cell_index, G);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            les_vars[cell_index].alpha_ij[i][j] = G[i][j];

    double delta = les_vars[cell_index].delta;

    double beta[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                beta[i][j] += delta * delta * G[k][i] * G[k][j];
            }
        }
    }

    double B_beta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
                    beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
                    beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];

    les_vars[cell_index].B_beta = max(0.0, B_beta);
}

double Vreman_Nu_SGS(double delta, double B_beta, double alpha_sq, double Cv)
{
    if (alpha_sq < 1e-14)
        return 0.0;

    return Cv * sqrt(B_beta / alpha_sq);
}

void Calculate_Vreman_Viscosity(int cell_index)
{
    Calculate_Vreman_B(cell_index);

    double alpha_sq = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            alpha_sq += les_vars[cell_index].alpha_ij[i][j] * les_vars[cell_index].alpha_ij[i][j];
        }
    }

    double delta = les_vars[cell_index].delta;
    double nu_sgs = Vreman_Nu_SGS(delta, les_vars[cell_index].B_beta, alpha_sq, Cv_Vreman);

    les_vars[cell_index].nu_sgs = nu_sgs;

    double rho = Primitive_Cells[cell_index][0];
    if (rho < 1e-14)
        rho = 1e-14;
    les_vars[cell_index].mu_sgs = rho * nu_sgs;
}

void Calculate_Vreman_Viscosity_All()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_Vreman_Viscosity(i);
        }
    }
}

//=============================================================================
// SIGMA MODEL
//=============================================================================

void Calculate_Singular_Values(int cell_index)
{
    double G[3][3];
    Calculate_Velocity_Gradient_Tensor(cell_index, G);

    double GtG[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                GtG[i][j] += G[k][i] * G[k][j];
            }
        }
    }

    double I1 = GtG[0][0] + GtG[1][1] + GtG[2][2];
    double I2 = GtG[0][0] * GtG[1][1] + GtG[0][0] * GtG[2][2] + GtG[1][1] * GtG[2][2] -
                GtG[0][1] * GtG[0][1] - GtG[0][2] * GtG[0][2] - GtG[1][2] * GtG[1][2];
    double I3 = GtG[0][0] * (GtG[1][1] * GtG[2][2] - GtG[1][2] * GtG[1][2]) -
                GtG[0][1] * (GtG[0][1] * GtG[2][2] - GtG[1][2] * GtG[0][2]) +
                GtG[0][2] * (GtG[0][1] * GtG[1][2] - GtG[1][1] * GtG[0][2]);

    double p = I1 / 3.0;
    double q = (p * p * p - I1 * p * p / 2.0 + I2 * p / 2.0 - I3 / 2.0);
    double r = I2 / 3.0 - p * p;

    double sigma_sq[3] = {0, 0, 0};

    if (fabs(r) < 1e-14)
    {
        sigma_sq[0] = sigma_sq[1] = sigma_sq[2] = p;
    }
    else
    {
        double phi = acos(max(-1.0, min(1.0, q / sqrt(-r * r * r)))) / 3.0;
        double sqrtr = sqrt(-r);

        sigma_sq[0] = p + 2.0 * sqrtr * cos(phi);
        sigma_sq[1] = p + 2.0 * sqrtr * cos(phi + 2.0 * M_PI / 3.0);
        sigma_sq[2] = p + 2.0 * sqrtr * cos(phi + 4.0 * M_PI / 3.0);
    }

    for (int i = 0; i < 3; ++i)
        sigma_sq[i] = max(0.0, sigma_sq[i]);

    sort(sigma_sq, sigma_sq + 3, greater<double>());

    les_vars[cell_index].sigma_1 = sqrt(sigma_sq[0]);
    les_vars[cell_index].sigma_2 = sqrt(sigma_sq[1]);
    les_vars[cell_index].sigma_3 = sqrt(sigma_sq[2]);
}

double Sigma_Nu_SGS(double delta, double sigma1, double sigma2, double sigma3, double Csigma)
{
    if (sigma1 < 1e-14)
        return 0.0;

    double D_sigma = sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3) / (sigma1 * sigma1);
    return Csigma * Csigma * delta * delta * D_sigma;
}

void Calculate_Sigma_Viscosity(int cell_index)
{
    Calculate_Singular_Values(cell_index);

    double delta = les_vars[cell_index].delta;
    double nu_sgs = Sigma_Nu_SGS(delta,
                                  les_vars[cell_index].sigma_1,
                                  les_vars[cell_index].sigma_2,
                                  les_vars[cell_index].sigma_3,
                                  Csigma);

    les_vars[cell_index].nu_sgs = max(0.0, nu_sgs);

    double rho = Primitive_Cells[cell_index][0];
    if (rho < 1e-14)
        rho = 1e-14;
    les_vars[cell_index].mu_sgs = rho * les_vars[cell_index].nu_sgs;
}

void Calculate_Sigma_Viscosity_All()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_Sigma_Viscosity(i);
        }
    }
}

//=============================================================================
// WALL DAMPING FUNCTIONS
//=============================================================================

void Calculate_Wall_Distance_LES(int cell_index)
{
    les_vars[cell_index].y_plus = 1000.0;

    if (!Cells[cell_index].has_Wall_Face)
        return;

    double min_dist = 1e30;
    int num_faces = static_cast<int>(Cells[cell_index].Neighbours.size());
    int num_boundary_kinds = static_cast<int>(Cells[cell_index].Face_Boundary_Kind.size());

    double cx = Cells[cell_index].Cell_Center[0];
    double cy = Cells[cell_index].Cell_Center[1];

    for (int f = 0; f < num_faces; ++f)
    {
        bool is_wall = false;
        if (f < num_boundary_kinds)
            is_wall = (Cells[cell_index].Face_Boundary_Kind[f] == 5);

        if (!is_wall)
            continue;

        int neigh = Cells[cell_index].Neighbours[f];
        if (neigh >= 0 && neigh < Total_No_Cells)
        {
            double wx = Cells[neigh].Cell_Center[0];
            double wy = Cells[neigh].Cell_Center[1];
            double dist = sqrt((cx - wx) * (cx - wx) + (cy - wy) * (cy - wy));
            min_dist = min(min_dist, dist);
        }
    }

    if (min_dist < 1e20)
    {
        double rho = Primitive_Cells[cell_index][0];
        double mu = Primitive_Cells[cell_index][8];
        if (rho < 1e-14)
            rho = 1e-14;
        if (mu < 1e-14)
            mu = 1e-14;

        double u_tau = 0.1;
        les_vars[cell_index].y_plus = min_dist * rho * u_tau / mu;
    }
}

void Calculate_Wall_Distances_LES_All()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_Wall_Distance_LES(i);
            Calculate_Van_Driest_Damping(i);
        }
    }
}

double Van_Driest_Function(double y_plus)
{
    return 1.0 - exp(-y_plus / LES_Constants::A_plus);
}

void Calculate_Van_Driest_Damping(int cell_index)
{
    double y_plus = les_vars[cell_index].y_plus;
    les_vars[cell_index].f_damping = Van_Driest_Function(y_plus);
}

//=============================================================================
// EFFECTIVE VISCOSITY UPDATE
//=============================================================================

void Update_LES_Effective_Viscosity()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            double mu_lam = Primitive_Cells[i][8];
            if (mu_lam < 1e-14)
                mu_lam = 1e-14;

            double mu_sgs = les_vars[i].mu_sgs;
            les_vars[i].mu_effective = mu_lam + mu_sgs;
        }
    }
}

void Update_LES_Effective_Conductivity()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            double mu_lam = Primitive_Cells[i][8];
            double mu_sgs = les_vars[i].mu_sgs;

            double k_lam = mu_lam * cp_ref / Pr;
            double k_sgs = mu_sgs * cp_ref / Pr_turb_LES;

            les_vars[i].k_effective = k_lam + k_sgs;
        }
    }
}

double Get_Effective_Viscosity_LES(int cell_index)
{
    if (!Is_LES || cell_index < 0 || cell_index >= static_cast<int>(les_vars.size()))
    {
        return Primitive_Cells[cell_index][8];
    }
    return les_vars[cell_index].mu_effective;
}

double Get_Effective_Conductivity_LES(int cell_index)
{
    if (!Is_LES || cell_index < 0 || cell_index >= static_cast<int>(les_vars.size()))
    {
        double mu = Primitive_Cells[cell_index][8];
        return mu * cp_ref / Pr;
    }
    return les_vars[cell_index].k_effective;
}

//=============================================================================
// MAIN SOLVER INTERFACE
//=============================================================================

void Update_LES_Variables(double dt)
{
    if (!Is_LES || current_les_model == LES_SGS_Model::NONE)
        return;

    switch (current_les_model)
    {
    case LES_SGS_Model::SMAGORINSKY:
        Calculate_Smagorinsky_Viscosity_All();
        break;

    case LES_SGS_Model::DYNAMIC_SMAGORINSKY:
        Calculate_Dynamic_Smagorinsky_All();
        break;

    case LES_SGS_Model::WALE:
        Calculate_WALE_Viscosity_All();
        break;

    case LES_SGS_Model::VREMAN:
        Calculate_Vreman_Viscosity_All();
        break;

    case LES_SGS_Model::SIGMA:
        Calculate_Sigma_Viscosity_All();
        break;

    default:
        break;
    }

    Update_LES_Effective_Viscosity();
    Update_LES_Effective_Conductivity();
}

void Solve_LES_Step(double dt)
{
    Update_LES_Variables(dt);
}

bool LES_Solver_Available()
{
    return Is_LES && (current_les_model != LES_SGS_Model::NONE);
}

void Setup_LES_From_Config(const Json::Value &config)
{
    if (!config.isMember("LES"))
    {
        Is_LES = false;
        current_les_model = LES_SGS_Model::NONE;
        return;
    }

    const Json::Value &lesConfig = config["LES"];

    Is_LES = lesConfig.get("Enabled", false).asBool();
    if (!Is_LES)
    {
        current_les_model = LES_SGS_Model::NONE;
        return;
    }

    string model_name = lesConfig.get("Model", "Smagorinsky").asString();

    if (model_name == "Smagorinsky" || model_name == "SMAGORINSKY")
        current_les_model = LES_SGS_Model::SMAGORINSKY;
    else if (model_name == "Dynamic" || model_name == "DynamicSmagorinsky" || model_name == "DYNAMIC_SMAGORINSKY")
        current_les_model = LES_SGS_Model::DYNAMIC_SMAGORINSKY;
    else if (model_name == "WALE" || model_name == "wale")
        current_les_model = LES_SGS_Model::WALE;
    else if (model_name == "Vreman" || model_name == "VREMAN")
        current_les_model = LES_SGS_Model::VREMAN;
    else if (model_name == "Sigma" || model_name == "SIGMA")
        current_les_model = LES_SGS_Model::SIGMA;
    else
    {
        cout << "Warning: Unknown LES model '" << model_name << "'. Using Smagorinsky." << endl;
        current_les_model = LES_SGS_Model::SMAGORINSKY;
    }

    Cs_Smagorinsky = lesConfig.get("Cs", LES_Constants::Cs_default).asDouble();
    Cw_WALE = lesConfig.get("Cw", LES_Constants::Cw_default).asDouble();
    Cv_Vreman = lesConfig.get("Cv", LES_Constants::Cv_default).asDouble();
    Csigma = lesConfig.get("Csigma", LES_Constants::Csigma_default).asDouble();
    Pr_turb_LES = lesConfig.get("Pr_t", LES_Constants::Pr_t).asDouble();
    Use_Wall_Damping = lesConfig.get("WallDamping", true).asBool();
    Use_GPU_LES = lesConfig.get("UseGPU", true).asBool();

    Initialize_LES_Model(current_les_model);
}

//=============================================================================
// VORTEX IDENTIFICATION
//=============================================================================

void Calculate_Q_Criterion(int cell_index)
{
    double S[3][3], Omega[3][3];
    Calculate_Strain_Rate_Tensor_LES(cell_index, S);
    Calculate_Rotation_Rate_Tensor(cell_index, Omega);

    double S_sq = 0.0, Omega_sq = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            S_sq += S[i][j] * S[i][j];
            Omega_sq += Omega[i][j] * Omega[i][j];
        }
    }

    les_vars[cell_index].Q_criterion = 0.5 * (Omega_sq - S_sq);
}

void Identify_Vortex_Structures()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            Calculate_Q_Criterion(i);
        }
    }
}

//=============================================================================
// OUTPUT AND STATISTICS
//=============================================================================

void Write_LES_Variables(const string &filename)
{
    ofstream outfile(filename);
    if (!outfile.is_open())
    {
        cerr << "Error: Cannot open file " << filename << " for writing." << endl;
        return;
    }

    outfile << "# LES Variables Output" << endl;
    outfile << "# Cell_ID, nu_sgs, mu_sgs, delta, S_mag, Q_criterion" << endl;

    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType != -1)
        {
            outfile << i << ", "
                    << les_vars[i].nu_sgs << ", "
                    << les_vars[i].mu_sgs << ", "
                    << les_vars[i].delta << ", "
                    << les_vars[i].S_mag << ", "
                    << les_vars[i].Q_criterion << endl;
        }
    }

    outfile.close();
    cout << "LES variables written to " << filename << endl;
}

//=============================================================================
// BOUNDARY CONDITIONS
//=============================================================================

void Apply_LES_Wall_BC(int cell_index)
{
    les_vars[cell_index].nu_sgs = 0.0;
    les_vars[cell_index].mu_sgs = 0.0;
}

void Apply_LES_Inlet_BC(int cell_index)
{
    // Keep SGS viscosity from calculation
}

void Apply_LES_Outlet_BC(int cell_index)
{
    // Zero gradient for SGS viscosity
    int interior_neighbor = -1;
    int num_neighbors = static_cast<int>(Cells[cell_index].Neighbours.size());
    
    for (int n = 0; n < num_neighbors; ++n)
    {
        int neigh = Cells[cell_index].Neighbours[n];
        if (neigh >= 0 && neigh < Total_No_Cells && Cells[neigh].cellType != -1)
        {
            interior_neighbor = neigh;
            break;
        }
    }

    if (interior_neighbor >= 0)
    {
        les_vars[cell_index].nu_sgs = les_vars[interior_neighbor].nu_sgs;
        les_vars[cell_index].mu_sgs = les_vars[interior_neighbor].mu_sgs;
    }
}

void Apply_All_LES_Boundary_Conditions()
{
    for (int i = 0; i < Total_No_Cells; ++i)
    {
        if (Cells[i].cellType == -1)
            continue;

        if (Cells[i].has_Wall_Face)
        {
            Apply_LES_Wall_BC(i);
        }
    }
}
