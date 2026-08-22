#include "definitions.h"
#include "Globals.h"
#include "Grid.h"
#include "Initialize.h"
#include "Primitive_Computational.h"
#include "Viscous_Functions.h"
#include "LES_Models.h"

#ifdef USE_CUDA
#include "definitions.h"
#include "Globals.h"
#include "Viscous_Functions.h"
#include "LES_Models.h"
#include "../CUDA_KERNELS/Viscous_Flux_Cuda_Integration.h"
#endif
#include "Utilities.h"
#include <filesystem>
#include <limits>
#include <regex>

bool Is_Leaf_Cell(const int &cellIndex)
{
    if (cellIndex < 0 || cellIndex >= static_cast<int>(Cells.size()))
        return false;
    if (cellIndex >= static_cast<int>(U_Cells.size()))
        return false;
    return Cells[cellIndex].AMR_IsLeaf;
}

void Build_Leaf_Cell_List(V_I &leafCells)
{
    leafCells.clear();
    leafCells.reserve(U_Cells.size());
    const int n = static_cast<int>(std::min(Cells.size(), U_Cells.size()));
    for (int i = 0; i < n; ++i)
    {
        if (Is_Leaf_Cell(i))
            leafCells.push_back(i);
    }
}

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
    /* Secondary_Neighbours = {SW, SE, NE, NW} diagonal stencil for co-volume
       face gradients. Face map: 0=left, 1=bottom, 2=right, 3=top. */
    if (Cells[Current_Cell_Index].Secondary_Neighbours.empty())
        Cells[Current_Cell_Index].Secondary_Neighbours.resize(4, -1);
    else
    {
        for (int i = 0; i < 4; ++i)
            Cells[Current_Cell_Index].Secondary_Neighbours[i] = -1;
    }

    const auto &Neighbours = Cells[Current_Cell_Index].Neighbours;
    if (Neighbours.size() < 4)
        return;

    const int N_L = Neighbours[0];
    const int N_B = Neighbours[1];
    const int N_R = Neighbours[2];
    const int N_T = Neighbours[3];

    const int n_all = static_cast<int>(Cells.size());
    auto is_valid = [n_all](int c) -> bool { return c >= 0 && c < n_all; };
    auto is_phys = [](int c) -> bool { return c >= 0 && c < No_Physical_Cells; };

    /* Neighbour-of-neighbour; never return Current_Cell_Index or -1 silently as 0. */
    auto getNeighbour = [&](int cellIndex, int direction) -> int {
        if (!is_valid(cellIndex) || direction < 0 || direction >= 4)
            return -1;
        if (static_cast<int>(Cells[cellIndex].Neighbours.size()) <= direction)
            return -1;
        const int n = Cells[cellIndex].Neighbours[direction];
        if (!is_valid(n) || n == Current_Cell_Index)
            return -1;
        return n;
    };

    auto pick = [&](int a, int b) -> int {
        if (is_valid(a) && a != Current_Cell_Index)
            return a;
        if (is_valid(b) && b != Current_Cell_Index)
            return b;
        return Current_Cell_Index; /* last resort: collapse stencil */
    };

    const bool ghost_L = !is_phys(N_L);
    const bool ghost_B = !is_phys(N_B);
    const bool ghost_R = !is_phys(N_R);
    const bool ghost_T = !is_phys(N_T);

    int SW = -1, SE = -1, NE = -1, NW = -1;

    if (ghost_L && ghost_B)
    {
        /* Wall/exit bottom-left corner (e.g. half-cylinder cell 0). */
        SW = pick(N_B, N_L);
        SE = pick(getNeighbour(N_R, 1), N_B);
        NE = pick(getNeighbour(N_R, 3), getNeighbour(N_T, 2));
        NW = pick(getNeighbour(N_T, 0), N_L);
    }
    else if (ghost_B && ghost_R)
    {
        /* Wall/exit bottom-right corner (e.g. cell 479). */
        SW = pick(getNeighbour(N_L, 1), N_B);
        SE = pick(N_B, N_R);
        NE = pick(getNeighbour(N_T, 2), N_R);
        NW = pick(getNeighbour(N_L, 3), getNeighbour(N_T, 0));
    }
    else if (ghost_R && ghost_T)
    {
        /* Inlet/exit top-right corner. */
        SW = pick(getNeighbour(N_L, 1), getNeighbour(N_B, 0));
        SE = pick(getNeighbour(N_B, 2), N_R);
        NE = pick(N_T, N_R);
        NW = pick(getNeighbour(N_L, 3), N_T);
    }
    else if (ghost_L && ghost_T)
    {
        /* Inlet/exit top-left corner. */
        SW = pick(getNeighbour(N_B, 0), N_L);
        SE = pick(getNeighbour(N_R, 1), getNeighbour(N_B, 2));
        NE = pick(getNeighbour(N_R, 3), N_T);
        NW = pick(N_L, N_T);
    }
    else
    {
        /* Interior / single-face boundary: walk around primary faces. */
        SW = pick(getNeighbour(N_L, 1), getNeighbour(N_B, 0));
        SE = pick(getNeighbour(N_B, 2), getNeighbour(N_R, 1));
        NE = pick(getNeighbour(N_R, 3), getNeighbour(N_T, 2));
        NW = pick(getNeighbour(N_T, 0), getNeighbour(N_L, 3));
    }

    Cells[Current_Cell_Index].Secondary_Neighbours[0] = SW;
    Cells[Current_Cell_Index].Secondary_Neighbours[1] = SE;
    Cells[Current_Cell_Index].Secondary_Neighbours[2] = NE;
    Cells[Current_Cell_Index].Secondary_Neighbours[3] = NW;
}

void Dump_Viscous_Corner_Diagnostics(const string &opfile)
{
    ofstream out(opfile);
    if (!out.is_open())
    {
        cerr << "Dump_Viscous_Corner_Diagnostics: cannot open " << opfile << endl;
        return;
    }

    /* Multi-BC physical cells (corners). */
    map<int, vector<string>> bc_by_cell;
    auto add_list = [&](const char *tag, const V_I &list) {
        for (size_t i = 0; i + 2 < list.size(); i += 3)
            bc_by_cell[list[i]].push_back(string(tag) + ":face=" + to_string(list[i + 1]) +
                                         ":ghost=" + to_string(list[i + 2]));
    };
    add_list("Inlet", Inlet_Cells_List);
    add_list("Exit", Exit_Cells_List);
    add_list("Wall", Wall_Cells_List);
    add_list("Symmetry", Symmetry_Cells_List);

    out << "=== Multi-BC corner cells ===\n";
    for (const auto &kv : bc_by_cell)
    {
        if (kv.second.size() < 2)
            continue;
        const int c = kv.first;
        if (c < 0 || c >= No_Physical_Cells || c >= static_cast<int>(Cells.size()))
            continue;
        const Cell &cell = Cells[c];
        out << "\nCell " << c;
        if (cell.Cell_Center.size() >= 2)
            out << " center=(" << cell.Cell_Center[0] << "," << cell.Cell_Center[1] << ")";
        out << "\n  BC:";
        for (const auto &s : kv.second)
            out << " " << s;
        out << "\n  Neighbours LBRT:";
        for (int f = 0; f < 4 && f < static_cast<int>(cell.Neighbours.size()); ++f)
            out << " " << cell.Neighbours[f];
        out << "\n  Secondary SW,SE,NE,NW:";
        for (int f = 0; f < 4 && f < static_cast<int>(cell.Secondary_Neighbours.size()); ++f)
            out << " " << cell.Secondary_Neighbours[f];
        out << "\n";
        for (int f = 0; f < 4; ++f)
        {
            const double nx = (2 * f + 1 < static_cast<int>(cell.Face_Normals.size()))
                                  ? cell.Face_Normals[2 * f]
                                  : 0.0;
            const double ny = (2 * f + 1 < static_cast<int>(cell.Face_Normals.size()))
                                  ? cell.Face_Normals[2 * f + 1]
                                  : 0.0;
            const double dl = (f < static_cast<int>(cell.Face_Areas.size())) ? cell.Face_Areas[f] : 0.0;
            out << "  Face " << f << " n=(" << nx << "," << ny << ") L=" << dl;
            if (c < static_cast<int>(Co_Volume_Cells.size()) &&
                !Co_Volume_Cells[c].Face_Normals.empty())
            {
                out << "  co-vol n0=(" << Co_Volume_Cells[c].Face_Normals[f * 8 + 0] << ","
                    << Co_Volume_Cells[c].Face_Normals[f * 8 + 1] << ")";
            }
            out << "\n";
        }
    }

    /* Wall face histogram + first/last entries. */
    out << "\n=== Wall face-number histogram ===\n";
    map<int, int> hist;
    for (size_t i = 0; i + 2 < Wall_Cells_List.size(); i += 3)
        hist[Wall_Cells_List[i + 1]]++;
    for (const auto &kv : hist)
        out << "  face " << kv.first << " count=" << kv.second << "\n";

    out << "\n=== BC list samples (first 3 + last 3 of each) ===\n";
    auto sample = [&](const char *name, const V_I &list) {
        out << name << " total=" << (list.size() / 3) << "\n";
        const size_t n = list.size() / 3;
        for (size_t k = 0; k < n && k < 3; ++k)
            out << "  [" << k << "] cell=" << list[3 * k] << " face=" << list[3 * k + 1]
                << " ghost=" << list[3 * k + 2] << "\n";
        for (size_t k = (n > 3 ? n - 3 : 3); k < n; ++k)
            out << "  [" << k << "] cell=" << list[3 * k] << " face=" << list[3 * k + 1]
                << " ghost=" << list[3 * k + 2] << "\n";
    };
    sample("Inlet", Inlet_Cells_List);
    sample("Exit", Exit_Cells_List);
    sample("Wall", Wall_Cells_List);

    out.close();
    cout << "Viscous corner diagnostics written to " << opfile << endl;
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
                throw runtime_error("Calculate_Green_Gauss_Gradient: invalid Grad_Type " + to_string(Grad_Type));
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
    const int n_all = static_cast<int>(Cells.size());
    auto safe_prim = [&](int c) -> double {
        if (c < 0 || c >= n_all || c >= static_cast<int>(Primitive_Cells.size()))
            return Primitive_Cells[Current_Cell_Index][Grad_Type];
        if (Grad_Type < 0 || Grad_Type >= static_cast<int>(Primitive_Cells[c].size()))
            return Primitive_Cells[Current_Cell_Index][Grad_Type];
        return Primitive_Cells[c][Grad_Type];
    };

    /* Primary (0..3) then Secondary SW,SE,NE,NW (4..7). */
    V_I neighbors(8, Current_Cell_Index);
    for (int i = 0; i < 4; ++i)
    {
        if (i < static_cast<int>(Cells[Current_Cell_Index].Neighbours.size()))
            neighbors[i] = Cells[Current_Cell_Index].Neighbours[i];
        if (i < static_cast<int>(Cells[Current_Cell_Index].Secondary_Neighbours.size()))
            neighbors[4 + i] = Cells[Current_Cell_Index].Secondary_Neighbours[i];
        if (neighbors[i] < 0 || neighbors[i] >= n_all)
            neighbors[i] = Current_Cell_Index;
        if (neighbors[4 + i] < 0 || neighbors[4 + i] >= n_all)
            neighbors[4 + i] = Current_Cell_Index;
    }

    /* 9 entries: current + 8 neighbours (was incorrectly sized to 8). */
    V_D cell_averages(9, 0.0);
    cell_averages[0] = safe_prim(Current_Cell_Index);
    for (int i = 0; i < 8; ++i)
        cell_averages[i + 1] = safe_prim(neighbors[i]);

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
        throw runtime_error("Viscous flux calculation produced non-finite values");
    }
}

//  This calculates Viscous flux on each face by averaging corresponding gradients evaluated at cell centers.
void Evaluate_Viscous_Fluxes()
{
    /* Refresh SGS viscosity into Primitive_Cells[*][8] before host/CUDA viscous flux. */
    if (Is_LES)
        Solve_LES_Step(Min_dt);

#ifdef USE_CUDA
    if (Evaluate_Viscous_Fluxes_CUDA())
        return;
#endif

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
    if (nTagged <= 0)
        return false;

    cout << "AMR: " << nTagged << " cells tagged for refinement (gradient threshold=" << AMR_Gradient_Threshold << ")" << endl;

    // Refine by switching to the next available structured grid level (Nx,Ny)->(2Nx-1,2Ny-1)
    // and transferring U by nearest-cell interpolation.
    static int amr_refinements_done = 0;
    const int AMR_MAX_REFINEMENTS = 2;
    if (amr_refinements_done >= AMR_MAX_REFINEMENTS)
    {
        cout << "AMR: max refinement levels reached (" << AMR_MAX_REFINEMENTS << ")." << endl;
        return false;
    }

    const int old_nx = meshParams.nx;
    const int old_ny = meshParams.ny;
    const int new_nx = 2 * old_nx - 1;
    const int new_ny = 2 * old_ny - 1;

    std::regex re("_(\\d+)_(\\d+)\\.(txt|vtk)$");
    std::smatch m;
    std::string candidate = Grid_File;
    if (std::regex_search(Grid_File, m, re))
    {
        candidate = std::regex_replace(Grid_File, re, "_" + std::to_string(new_nx) + "_" + std::to_string(new_ny) + "." + m[3].str());
    }
    else
    {
        cout << "AMR: unable to infer next grid filename from " << Grid_File << endl;
        return false;
    }

    if (!std::filesystem::exists(candidate))
    {
        cout << "AMR: no finer grid file found: " << candidate << endl;
        return false;
    }

    // Backup old solution on physical cells.
    const int old_n_phys = No_Physical_Cells;
    vector<V_D> old_centers(old_n_phys, V_D(2, 0.0));
    vector<V_D> old_u(old_n_phys, V_D(4, 0.0));
    for (int i = 0; i < old_n_phys; ++i)
    {
        old_centers[i] = Cells[i].Cell_Center;
        old_u[i] = U_Cells[i];
    }

    cout << "AMR: refining grid " << old_nx << "x" << old_ny << " -> " << new_nx << "x" << new_ny << endl;

    // Rebuild mesh on finer grid.
    Grid_File = candidate;
    meshParams.nx = new_nx;
    meshParams.ny = new_ny;
    Cells.clear();
    Boundary_Cells.clear();
    Co_Volume_Cells.clear();
    Total_No_Cells = 0;
    No_Physical_Cells = 0;
    No_Ghost_Cells = 0;
    Inlet_Cells_List.clear();
    Exit_Cells_List.clear();
    Wall_Cells_List.clear();
    Symmetry_Cells_List.clear();
    if (!Form_Cells(Grid_File))
    {
        cout << "AMR: Form_Cells failed on refined grid, keeping current grid." << endl;
        return false;
    }

    // Reinitialize solver storage for new grid.
    Cells_Net_Flux.clear();
    Cells_DelU.clear();
    Cells_Face_Boundary_Type.clear();
    U_Cells.clear();
    U_Cells_RK_1.clear();
    U_Cells_RK_2.clear();
    Primitive_Cells.clear();
    Cells_Viscous_Flux.clear();
    CF.clear();
    Initialize(Test_Case);

    // Nearest-cell interpolation of conservative variables onto refined physical cells.
    for (int i = 0; i < No_Physical_Cells; ++i)
    {
        const double cx = Cells[i].Cell_Center[0];
        const double cy = Cells[i].Cell_Center[1];
        int best = 0;
        double best_d2 = std::numeric_limits<double>::max();
        for (int j = 0; j < old_n_phys; ++j)
        {
            const double dx = cx - old_centers[j][0];
            const double dy = cy - old_centers[j][1];
            const double d2 = dx * dx + dy * dy;
            if (d2 < best_d2)
            {
                best_d2 = d2;
                best = j;
            }
        }
        U_Cells[i] = old_u[best];
        Calculate_Primitive_Variables(i, U_Cells[i], Primitive_Cells[i]);
    }

    amr_refinements_done++;
    cout << "AMR: refinement applied using grid file " << Grid_File << endl;
    return true;
}
