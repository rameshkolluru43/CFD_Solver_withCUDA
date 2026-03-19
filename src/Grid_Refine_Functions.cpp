#include "definitions.h"
#include "Globals.h"
#include "Viscous_Functions.h"
#include "Utilities.h"

// Compute per-cell gradient-based refinement indicator (density + pressure), scale-invariant.
// indicator = |grad(rho)|*sqrt(Area) + (|grad(P)|/max(P,eps))*sqrt(Area). Stored in Gradient_Refinement_Indicator.
void Compute_Gradient_Refinement_Indicator()
{
    if (Gradient_Refinement_Indicator.size() != static_cast<size_t>(No_Physical_Cells))
        Gradient_Refinement_Indicator.resize(No_Physical_Cells, 0.0);

    int Grad_Type_Rho = 0, Grad_Type_P = 4;
    V_D grad(2, 0.0);
    const double eps = 1e-14;

    for (int i = 0; i < No_Physical_Cells; i++)
    {
        double ind = 0.0;
        const double h = sqrt(Cells[i].Area);

        Calculate_Gradient_At_Cell_Center(i, Grad_Type_Rho, grad);
        ind += sqrt(grad[0] * grad[0] + grad[1] * grad[1]) * h;

        Calculate_Gradient_At_Cell_Center(i, Grad_Type_P, grad);
        double P = Primitive_Cells[i][4];
        if (P > eps)
            ind += (sqrt(grad[0] * grad[0] + grad[1] * grad[1]) / P) * h;

        Gradient_Refinement_Indicator[i] = ind;
    }
}

// Tag cells for refinement: indicator > threshold, optionally cap by fraction (refine top N% by indicator).
void TagRefinableCells(vector<Cell> &cells, double &threshold)
{
    Compute_Gradient_Refinement_Indicator();

    for (int i = 0; i < No_Physical_Cells; i++)
        cells[i].Is_Splittable = false;

    if (AMR_Max_Fraction > 0.0 && AMR_Max_Fraction < 1.0)
    {
        // Refine top AMR_Max_Fraction of cells by indicator value
        vector<pair<double, int>> order(No_Physical_Cells);
        for (int i = 0; i < No_Physical_Cells; i++)
            order[i] = {Gradient_Refinement_Indicator[i], i};
        sort(order.begin(), order.end(), [](const pair<double, int> &a, const pair<double, int> &b)
             { return a.first > b.first; });
        int nRefine = static_cast<int>(AMR_Max_Fraction * No_Physical_Cells);
        if (nRefine < 1)
            nRefine = 1;
        for (int k = 0; k < nRefine && k < No_Physical_Cells; k++)
            if (order[k].first >= threshold)
                cells[order[k].second].Is_Splittable = true;
    }
    else
    {
        for (int i = 0; i < No_Physical_Cells; i++)
            if (Gradient_Refinement_Indicator[i] > threshold)
                cells[i].Is_Splittable = true;
    }
}

// Helper function to calculate vertex averages
double Calculate_Vertex_Average(const V_D &weights, const V_D &cell_averages)
{
    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        weighted_sum += weights[i] * cell_averages[i];
        weight_sum += weights[i];
    }
    return weighted_sum / weight_sum;
}

// Green-Gauss gradient: grad = (1/Area) * sum_f (phi_face * n_f * dl_f). Works for any polygon.
void Calculate_Gradient(V_D &av, V_D &nx, V_D &ny, V_D &dl, double &inv_area, V_D &grad)
{
    if (grad.empty())
        grad.resize(2, 0.0);
    else
    {
        grad[0] = 0.0;
        grad[1] = 0.0;
    }
    const int n = static_cast<int>(av.size());
    for (int i = 0; i < n; ++i)
    {
        grad[0] += av[i] * nx[i] * dl[i];
        grad[1] += av[i] * ny[i] * dl[i];
    }
    grad[0] *= inv_area;
    grad[1] *= inv_area;
}
// Function to identify secondary neighbors for second gradients
/**
 * @brief Identifies the secondary neighbors of a given cell to evaluate second gradients.
 *
 * This function determines the secondary neighbors of the cell specified by `Current_Cell_Index`.
 * Secondary neighbors are used for extended stencil calculations in gradient evaluation.
 *
 * @param Current_Cell_Index Reference to the index of the current cell for which neighbors are identified.
 */

void Identify_Neighbours_For_Second_Gradients(int &Current_Cell_Index)
{
    // Initialize secondary neighbors to zero
    if (Cells[Current_Cell_Index].Secondary_Neighbours.empty())
        Cells[Current_Cell_Index].Secondary_Neighbours.resize(4, 0.0);
    else
    {
        for (int i = 0; i < 4; ++i)
            Cells[Current_Cell_Index].Secondary_Neighbours[i] = 0.0;
    }

    // Fetch primary neighbors
    const auto &Neighbours = Cells[Current_Cell_Index].Neighbours;
    int Neighbour_1 = Neighbours[0];
    int Neighbour_2 = Neighbours[1];
    int Neighbour_3 = Neighbours[2];
    int Neighbour_4 = Neighbours[3];

    // Lambda to fetch a neighbor safely
    auto getNeighbour = [](int cellIndex, int direction) -> int
    {
        return (cellIndex >= 0 && cellIndex < No_Physical_Cells) ? Cells[cellIndex].Neighbours[direction] : 0;
    };

    // Determine secondary neighbors based on boundary conditions
    if (Neighbour_1 >= No_Physical_Cells && Neighbour_2 >= No_Physical_Cells) // Bottom-left corner
    {
        Cells[Current_Cell_Index].Secondary_Neighbours[0] = Neighbour_2;
        Cells[Current_Cell_Index].Secondary_Neighbours[1] = getNeighbour(Neighbour_3, 2);
        Cells[Current_Cell_Index].Secondary_Neighbours[2] = getNeighbour(Neighbour_4, 1);
        Cells[Current_Cell_Index].Secondary_Neighbours[3] = getNeighbour(Neighbour_1, 0);
    }
    else if (Neighbour_2 >= No_Physical_Cells && Neighbour_3 >= No_Physical_Cells) // Bottom-right corner
    {
        Cells[Current_Cell_Index].Secondary_Neighbours[0] = getNeighbour(Neighbour_1, 1);
        Cells[Current_Cell_Index].Secondary_Neighbours[1] = Neighbour_2;
        Cells[Current_Cell_Index].Secondary_Neighbours[2] = getNeighbour(Neighbour_4, 2);
        Cells[Current_Cell_Index].Secondary_Neighbours[3] = getNeighbour(Neighbour_1, 3);
    }
    else if (Neighbour_3 >= No_Physical_Cells && Neighbour_4 >= No_Physical_Cells) // Top-right corner
    {
        Cells[Current_Cell_Index].Secondary_Neighbours[0] = getNeighbour(Neighbour_1, 1);
        Cells[Current_Cell_Index].Secondary_Neighbours[1] = getNeighbour(Neighbour_2, 0);
        Cells[Current_Cell_Index].Secondary_Neighbours[2] = Neighbour_4;
        Cells[Current_Cell_Index].Secondary_Neighbours[3] = getNeighbour(Neighbour_1, 3);
    }
    else if (Neighbour_1 >= No_Physical_Cells && Neighbour_4 >= No_Physical_Cells) // Top-left corner
    {
        Cells[Current_Cell_Index].Secondary_Neighbours[0] = getNeighbour(Neighbour_2, 0);
        Cells[Current_Cell_Index].Secondary_Neighbours[1] = getNeighbour(Neighbour_3, 2);
        Cells[Current_Cell_Index].Secondary_Neighbours[2] = getNeighbour(Neighbour_4, 3);
        Cells[Current_Cell_Index].Secondary_Neighbours[3] = Neighbour_1;
    }
    else // Interior or boundary cells
    {
        Cells[Current_Cell_Index].Secondary_Neighbours[0] = getNeighbour(Neighbour_1, 1);
        Cells[Current_Cell_Index].Secondary_Neighbours[1] = getNeighbour(Neighbour_2, 2);
        Cells[Current_Cell_Index].Secondary_Neighbours[2] = getNeighbour(Neighbour_3, 3);
        Cells[Current_Cell_Index].Secondary_Neighbours[3] = getNeighbour(Neighbour_4, 0);
    }
}

// Evaluate gradients at cell centers

void Calculate_Gradients_At_Cell_Centers()
{
    // Loop over all physical cells
    V_D gradient(2, 0.0);
    for (int Current_Cell_Index = 0; Current_Cell_Index < No_Physical_Cells; ++Current_Cell_Index)
    {
        // Calculate gradients for each variable
        for (int Grad_Type = 0; Grad_Type < 5; ++Grad_Type)
        {
            // Calculate gradient at cell center
            Calculate_Gradient_At_Cell_Center(Current_Cell_Index, Grad_Type, gradient);
            // Grad_type = 0, density
            // Grad_type = 1, u velocity
            // Grad_type = 2, v velocity
            // Grad_type = 3, temperature
            // Gradent = 4, pressure
            // Store gradient in the appropriate cell
            // without if else store the appropriate gradients in rho_Gradient, u_Gradient, v_Gradient, T_Gradient,P_Gradient

            if (Grad_Type == 0)
            {
                Rho_Gradient = gradient;
            }
            else if (Grad_Type == 1)
            {
                u_Gradient = gradient;
            }
            else if (Grad_Type == 2)
            {
                v_Gradient = gradient;
            }
            else if (Grad_Type == 3)
            {
                T_Gradient = gradient;
            }
            else if (Grad_Type == 4)
            {
                P_Gradient = gradient;
            }
            else
            {
                cout << "Invalid Grad_Type: " << Grad_Type << endl;
                exit(0);
            }
        }
    }
}

// Green-Gauss gradient at cell center: phi_face = 0.5*(phi_c + phi_neighbor), grad = (1/Area)*sum_f (phi_face * n_f * dl_f).
// Works for any polygon (tri/quad); uses this cell's Face_Normals and Face_Areas.
void Calculate_Gradient_At_Cell_Center(int &Current_Cell_Index, int &Grad_Type, V_D &grad)
{
    grad.assign(2, 0.0);
    const Cell &cell = Cells[Current_Cell_Index];
    const int nF = (cell.numFaces > 0) ? cell.numFaces : static_cast<int>(cell.Face_Areas.size());
    if (nF <= 0 || cell.Area <= 0.0)
        return;

    const double phi_c = Primitive_Cells[Current_Cell_Index][Grad_Type];
    for (int f = 0; f < nF; f++)
    {
        const int idx = f * 2;
        const double nx_f = cell.Face_Normals[idx + 0];
        const double ny_f = cell.Face_Normals[idx + 1];
        const double dl_f = cell.Face_Areas[f];
        int neigh = (f < static_cast<int>(cell.Neighbours.size())) ? cell.Neighbours[f] : -1;
        double phi_n = phi_c;
        if (neigh >= 0 && neigh < No_Physical_Cells)
            phi_n = Primitive_Cells[neigh][Grad_Type];
        const double phi_face = 0.5 * (phi_c + phi_n);
        grad[0] += phi_face * nx_f * dl_f;
        grad[1] += phi_face * ny_f * dl_f;
    }
    grad[0] *= cell.Inv_Area;
    grad[1] *= cell.Inv_Area;
}

void Calculate_Vertex_Average(const V_D &weights, const V_D &cell_averages, V_D &av)
{
    // Calculate vertex averages based on weights and cell averages
    for (int i = 0; i < 4; ++i)
    {
        av[i] = Calculate_Vertex_Average(weights, cell_averages);
    }
}

/**
 * @brief Calculate the gradient on a face of a cell.
 *
 * This function calculates the gradient on a face of a cell using the provided cell index, gradient type, and face number.
 *
 * @param Current_Cell_Index Reference to the index of the current cell for which the gradient is calculated.
 * @param Grad_Type Reference to the type of gradient to calculate (e.g., velocity, temperature, etc.).
 * @param Face_No Reference to the face number (0 to 3) on which to calculate the gradient.
 */
void Calculate_Gradient_On_Face(const int &Current_Cell_Index, const int &Grad_Type, const int &Face_No)
{
    // Retrieve neighbors
    V_I neighbors(8, 0);
    neighbors = {
        Cells[Current_Cell_Index].Neighbours[0],
        Cells[Current_Cell_Index].Neighbours[1],
        Cells[Current_Cell_Index].Neighbours[2],
        Cells[Current_Cell_Index].Neighbours[3],
        Cells[Current_Cell_Index].Secondary_Neighbours[0],
        Cells[Current_Cell_Index].Secondary_Neighbours[1],
        Cells[Current_Cell_Index].Secondary_Neighbours[2],
        Cells[Current_Cell_Index].Secondary_Neighbours[3]};

    // Retrieve cell averages
    V_D cell_averages(8, 0.0);
    cell_averages = {
        Primitive_Cells[Current_Cell_Index][Grad_Type], // Current Cell
        Primitive_Cells[neighbors[0]][Grad_Type],       // Neighbour 1
        Primitive_Cells[neighbors[1]][Grad_Type],       // Neighbour 2
        Primitive_Cells[neighbors[2]][Grad_Type],       // Neighbour 3
        Primitive_Cells[neighbors[3]][Grad_Type],       // Neighbour 4
        Primitive_Cells[neighbors[4]][Grad_Type],       // Neighbour 5
        Primitive_Cells[neighbors[5]][Grad_Type],       // Neighbour 6
        Primitive_Cells[neighbors[6]][Grad_Type],       // Neighbour 7
        Primitive_Cells[neighbors[7]][Grad_Type]};      // Neighbour 8

    // Calculate weights
    V_D weights(9, 1.0); // Default to equal weights
    if (Area_Weighted_Average == 3)
    { // Area Weighted Average
        for (int i = 0; i < 9; ++i)
        {
            // weights[i] = Cells_Area[i == 0 ? Current_Cell_Index : neighbors[i - 1]];
            weights[i] = (i == 0 ? Cells[Current_Cell_Index].Area : Cells[neighbors[i - 1]].Area);
        }
    }
    else if (Area_Weighted_Average == 2)
    { // Inverse Area Weighted Average
        for (int i = 0; i < 9; ++i)
        {
            // weights[i] = 1.0 / Cells_Area[i == 0 ? Current_Cell_Index : neighbors[i - 1]];
            weights[i] = (i == 0 ? Cells[Current_Cell_Index].Inv_Area : Cells[neighbors[i - 1]].Inv_Area);
        }
    }

    // Retrieve face normals and areas obtained from co-volume cells
    V_D nx(4, 0.0), ny(4, 0.0), dl(4, 0.0);
    for (int i = 0; i < 4; ++i)
    {
        nx[i] = Co_Volume_Cells[Current_Cell_Index].Face_Normals[Face_No * 8 + i * 2];
        ny[i] = Co_Volume_Cells[Current_Cell_Index].Face_Normals[Face_No * 8 + i * 2 + 1];
        dl[i] = Co_Volume_Cells[Current_Cell_Index].Face_Areas[Face_No * 4 + i];
    }
    double inv_area = Co_Volume_Cells[Current_Cell_Index].Inv_Area;

    // Calculate vertex averages based on Face_No
    V_D av(4, 0.0);
    // Values to store vertex averages of variables
    double o = 0.0, a = 0.0, b = 0.0, c = 0.0;
    // Calculate vertex averages based on weights and cell averages

    switch (Face_No)
    {
    case 0: // For Face_0
        //		cout<<"Face 0\t"<<endl;
        // Vertices of the face are o and c from the primary cell, with cell center as current cell
        // Calculate the vertex averages from cell averages using the weights
        // current cell, Neighbour 1, Neighbour 2, Neighbour 5
        o = Calculate_Vertex_Average({weights[0], weights[1], weights[2], weights[5]}, {cell_averages[0], cell_averages[1], cell_averages[2], cell_averages[5]});
        // current cell, Neighbour 1, Neighbour 4, Neighbour 8
        c = Calculate_Vertex_Average({weights[0], weights[1], weights[4], weights[8]}, {cell_averages[0], cell_averages[1], cell_averages[4], cell_averages[8]});

        // Calculate face averages from the verticies of co volume cells

        av[0] = 0.5 * (cell_averages[0] + o); // Average for Face 0
        av[1] = 0.5 * (cell_averages[0] + c); // Average for Face 1
        av[2] = 0.5 * (cell_averages[1] + c); // Average for Face 2
        av[3] = 0.5 * (cell_averages[1] + o); // Average for Face 3
        break;
    case 1: // For Face_1
        //				cout<<"Face 1\t"<<endl;
        // Vertices of the face are o and a from the primary cell, with cell center as current cell
        // Current cell, Neighbour 1, Neighbour 2, Neighbour 5
        o = Calculate_Vertex_Average({weights[0], weights[1], weights[2], weights[5]}, {cell_averages[0], cell_averages[1], cell_averages[2], cell_averages[5]});
        // Current cell, Neighbour 2, Neighbour 6, Neighbour 3
        a = Calculate_Vertex_Average({weights[0], weights[2], weights[5], weights[3]}, {cell_averages[0], cell_averages[2], cell_averages[5], cell_averages[3]});
        // Calculate face averages from the verticies of co volume cells
        av[0] = 0.5 * (cell_averages[0] + a); // Average for Face 0
        av[1] = 0.5 * (cell_averages[0] + o); // Average for Face 1
        av[2] = 0.5 * (cell_averages[2] + o); // Average for Face 2
        av[3] = 0.5 * (cell_averages[2] + a); // Average for Face 3
        break;
    case 2: // For Face_2
        //				cout<<"Face 2\t"<<endl;
        // Vertices of the face are a and b from the primary cell, with cell center as current cell
        // Current cell, Neighbour 2, Neighbour 3, Neighbour 6
        a = Calculate_Vertex_Average({weights[0], weights[2], weights[3], weights[6]}, {cell_averages[0], cell_averages[2], cell_averages[3], cell_averages[6]});
        // Current cell, Neighbour 3, Neighbour 4, Neighbour 7
        b = Calculate_Vertex_Average({weights[0], weights[3], weights[4], weights[7]}, {cell_averages[0], cell_averages[3], cell_averages[4], cell_averages[7]});
        // Calculate face averages from the verticies of co volume cells
        av[0] = 0.5 * (cell_averages[0] + a); // Average for Face 0
        av[1] = 0.5 * (cell_averages[0] + b); // Average for Face 1
        av[2] = 0.5 * (cell_averages[3] + b); // Average for Face 2
        av[3] = 0.5 * (cell_averages[3] + a); // Average for Face 3

        break;
    case 3: // For Face_3
        //				cout<<"Face 3\t"<<endl;
        // Vertices of the face are b and c from the primary cell, with cell center as current cell
        // Current cell, Neighbour 3, Neighbour 4, Neighbour 7
        b = Calculate_Vertex_Average({weights[0], weights[3], weights[4], weights[7]}, {cell_averages[0], cell_averages[3], cell_averages[4], cell_averages[7]});
        // Current cell, Neighbour 1, Neighbour 4  , Neighbour 8
        c = Calculate_Vertex_Average({weights[0], weights[1], weights[4], weights[8]}, {cell_averages[0], cell_averages[1], cell_averages[4], cell_averages[8]});
        // Calculate face averages from the verticies of co volume cells
        av[0] = 0.5 * (cell_averages[0] + b); // Average for Face 0
        av[1] = 0.5 * (cell_averages[0] + c); // Average for Face 1
        av[2] = 0.5 * (cell_averages[4] + c); // Average for Face 2
        av[3] = 0.5 * (cell_averages[4] + b); // Average for Face 3
        break;
    }

    V_D gradient(2, 0.0);
    // Calculate gradients
    Calculate_Gradient(av, nx, ny, dl, inv_area, gradient);

    // Assign gradients based on Grad_Type
    switch (Grad_Type)
    {
    case 0:
        Rho_Gradient[0] = gradient[0];
        Rho_Gradient[1] = gradient[1];
        break;
    case 1:
        u_Gradient[0] = gradient[0];
        u_Gradient[1] = gradient[1];
        break;
    case 2:
        v_Gradient[0] = gradient[0];
        v_Gradient[1] = gradient[1];
        break;
    case 3:
        T_Gradient[0] = gradient[0];
        T_Gradient[1] = gradient[1];
        break;
    case 4:
        P_Gradient[0] = gradient[0];
        P_Gradient[1] = gradient[1];
        break;
    }
}

void Viscous_Flux_on_Face(const int &Cell_No, const int &Face_No)
{
    //	cout<<"Calculating Viscous flux for Cell No\t"<<Cell_No<<"\t on face \t"<<Face_No<<endl;
    double u11, u12, u21, u22, T11, T12, T21, T22;
    double v1 = 0.0, v2 = 0.0, mu = 0.0;
    int index = Face_No * 2, Grad_Type = 0, N_Cell_No = 0;

    Qx = 0.0;
    Qy = 0.0;

    // Normals and Face Area (in 2D Length of the face)

    nx = Cells[Cell_No].Face_Normals[index];
    ny = Cells[Cell_No].Face_Normals[index + 1];
    dl = Cells[Cell_No].Face_Areas[Face_No];

    // Fetching the Corresponding Neighbours for the given cell and given face
    N_Cell_No = Cells[Cell_No].Neighbours[Face_No];

    // Find the arthematic average of the velocities and Temperature on the Face

    v1 = 0.5 * (Primitive_Cells[Cell_No][1] + Primitive_Cells[N_Cell_No][1]);
    v2 = 0.5 * (Primitive_Cells[Cell_No][2] + Primitive_Cells[N_Cell_No][2]);
    mu = 0.5 * (Primitive_Cells[Cell_No][8] + Primitive_Cells[N_Cell_No][8]);

    //	cout<<"Velocities on Face\t"<<v1<<"\t"<<v2<<endl;

    Grad_Type = 1; // U velocity Gradient		// Finds the u Velocity gradient on the face required
    Calculate_Gradient_On_Face(Cell_No, Grad_Type, Face_No);
    u11 = u_Gradient[0];
    u12 = u_Gradient[1];
    Grad_Type = 2; // V Velocity Gradient
    Calculate_Gradient_On_Face(Cell_No, Grad_Type, Face_No);
    u21 = v_Gradient[0];
    u22 = v_Gradient[1];
    Grad_Type = 3; // Temperature Gradient
    Calculate_Gradient_On_Face(Cell_No, Grad_Type, Face_No);
    //	cout<<"Gradients\t"<<Cell_No<<"\t"<<u11<<"\t"<<u12<<"\t"<<u21<<"\t"<<u22<<"\t"<<T_Gradient[0]<<"\t"<<T_Gradient[1]<<endl;
    //	cout<<T<<"\t"<<mu<<"\t"<<Re<<"\t"<<Pr<<endl;
    T11 = (2.0 / 3.0) * mu * Inv_Re * (2.0 * u11 - u22);
    T12 = mu * Inv_Re * (u12 + u21);
    T21 = T12;
    T22 = (2.0 / 3.0) * mu * Inv_Re * (2.0 * u22 - u11);
    //	cout<<T11<<"\t"<<T12<<"\t"<<T21<<"\t"<<T22<<endl;
    //	K1 = mu/(gamma_M_1*M_ref*M_ref*Re*Pr);
    Qx = K1 * T_Gradient[0];
    Qy = K1 * T_Gradient[1];

    Cells_Viscous_Flux[Cell_No][0] += 0.0;
    Cells_Viscous_Flux[Cell_No][1] += (T11 * nx + T21 * ny) * dl; // Ti1Ai
    Cells_Viscous_Flux[Cell_No][2] += (T12 * nx + T22 * ny) * dl; // Ti2Ai
    Cells_Viscous_Flux[Cell_No][3] += ((T11 * v1 + T12 * v2 + Qx) * nx * dl + (T21 * v1 + T22 * v2 + Qy) * ny * dl);

    if (isnan(T11) or isnan(T21) or isnan(T22) or isnan(Qx) or isnan(Qy))
    {

        cout << Cell_No << "\t" << Face_No << "\t" << T11 << "\t" << T12 << "\t" << T21 << "\t" << T22 << "\t" << Qx << "\t" << Qy << endl;
        cout << mu << "\t" << Inv_Re << "\t" << K1 << endl;
        cout << u11 << "\t" << u12 << "\t" << u21 << "\t" << u22 << endl;
        exit(0);
    }
}

//  This calculates Viscous flux on each face by averaging corresponding gradients evaluated at cell centers.
void Evaluate_Viscous_Fluxes()
{

    //	cout<< "Evaluating Viscous Fluxes\n";
    for (int Current_Cell_Index = 0; Current_Cell_Index < No_Physical_Cells; Current_Cell_Index++)
    {
        // Resetting all viscous flux to zero
        for (int i = 0; i < Cells_Viscous_Flux[Current_Cell_Index].size(); i++)
        {
            Cells_Viscous_Flux[Current_Cell_Index][i] = 0.0;
        }
        // Evaluating Viscous flux which includes heat flux on each face for a given cell
        Viscous_Flux_on_Face(Current_Cell_Index, Face_0);
        Viscous_Flux_on_Face(Current_Cell_Index, Face_1);
        Viscous_Flux_on_Face(Current_Cell_Index, Face_2);
        Viscous_Flux_on_Face(Current_Cell_Index, Face_3);
    }
    //	 	cout<<"Evaluating Viscous Fluxes Done"<<endl;
}

// Apply gradient-based adaptive refinement: compute indicator, tag cells. Returns true if mesh was changed (future: actual refinement).
bool Apply_Adaptive_Refinement()
{
    TagRefinableCells(Cells, AMR_Gradient_Threshold);
    int nTagged = 0;
    for (int i = 0; i < No_Physical_Cells; i++)
        if (Cells[i].Is_Splittable)
            nTagged++;
    if (nTagged > 0)
        cout << "AMR: " << nTagged << " cells tagged for refinement (gradient threshold=" << AMR_Gradient_Threshold << ")" << endl;
    return false; // no mesh change in this implementation (tagging only)
}
