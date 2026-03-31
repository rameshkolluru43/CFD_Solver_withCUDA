#include "AMR.hpp"
#include "Primitive_Computational.h"
#include "Viscous_Functions.h"
#include "Utilities.h"
#include <algorithm>
#include <limits>

namespace
{
double vec_norm2(double x, double y)
{
    return sqrt(x * x + y * y);
}

bool is_quad_leaf_cell(const int c)
{
    if (c < 0 || c >= No_Physical_Cells)
        return false;
    if (!Cells[c].AMR_IsLeaf)
        return false;
    return Cells[c].numFaces == 4;
}

void append_child_storage(const V_D &u_parent, const V_D &prim_parent, bool hasBoundary)
{
    const int pos = No_Physical_Cells;
    U_Cells.insert(U_Cells.begin() + pos, u_parent);
    U_Cells_RK_1.insert(U_Cells_RK_1.begin() + pos, u_parent);
    U_Cells_RK_2.insert(U_Cells_RK_2.begin() + pos, u_parent);
    Primitive_Cells.insert(Primitive_Cells.begin() + pos, prim_parent);
    Cells_Net_Flux.insert(Cells_Net_Flux.begin() + pos, V_D(NUM_FLUX_COMPONENTS, 0.0));
    Cells_DelU.insert(Cells_DelU.begin() + pos, V_D(NUM_FLUX_COMPONENTS, 0.0));
    Cells_Viscous_Flux.insert(Cells_Viscous_Flux.begin() + pos, V_D(NUM_FLUX_COMPONENTS, 0.0));
    Cells_Face_Boundary_Type.insert(Cells_Face_Boundary_Type.begin() + pos, vector<bool>(4, hasBoundary));
    Gradient_Refinement_Indicator.insert(Gradient_Refinement_Indicator.begin() + pos, 0.0);
}

void remap_neighbor_reference(int neigh, int oldParent, int newChild)
{
    if (neigh < 0 || neigh >= No_Physical_Cells)
        return;
    auto &N = Cells[neigh].Neighbours;
    for (int f = 0; f < static_cast<int>(N.size()); ++f)
    {
        if (N[f] == oldParent)
            N[f] = newChild;
    }
}

void remap_neighbor_reference_all(int neigh, const vector<int> &oldCells, int newCell)
{
    if (neigh < 0 || neigh >= No_Physical_Cells)
        return;
    auto &N = Cells[neigh].Neighbours;
    for (int f = 0; f < static_cast<int>(N.size()); ++f)
    {
        for (const int oldc : oldCells)
        {
            if (N[f] == oldc)
            {
                N[f] = newCell;
                break;
            }
        }
    }
}

int pick_first_valid_neighbor(int a, int b)
{
    if (a >= 0 && a < No_Physical_Cells)
        return a;
    if (b >= 0 && b < No_Physical_Cells)
        return b;
    return -1;
}

void shift_index_vectors_after_insert(const int pos)
{
    auto shift_vec = [&](vector<int> &v)
    {
        for (int &x : v)
            if (x >= pos)
                x++;
    };

    for (Cell &c : Cells)
    {
        shift_vec(c.Neighbours);
        shift_vec(c.Secondary_Neighbours);
        shift_vec(c.AMR_Children);
        if (c.AMR_Parent >= pos)
            c.AMR_Parent++;
    }

    for (int &x : Inlet_Cells_List)
        if (x >= pos)
            x++;
    for (int &x : Exit_Cells_List)
        if (x >= pos)
            x++;
    for (int &x : Wall_Cells_List)
        if (x >= pos)
            x++;
    for (int &x : Symmetry_Cells_List)
        if (x >= pos)
            x++;
}

Cell make_child_from_parent(const Cell &p, const int childType, const int parentIdx, const int childID)
{
    // childType: 0=SW, 1=SE, 2=NE, 3=NW
    Cell ch = p;
    ch.cellID = childID;
    ch.AMR_Parent = parentIdx;
    ch.AMR_Level = p.AMR_Level + 1;
    ch.AMR_IsLeaf = true;
    ch.AMR_Children.clear();
    ch.ParentCellID = parentIdx;
    ch.Is_Splittable = false;
    ch.numFaces = 4;
    ch.numNodes = 4;
    ch.Neighbours.assign(4, -1);
    ch.Secondary_Neighbours.assign(4, -1);
    ch.faceID.assign(4, -1);
    ch.Face_Boundary_Kind.assign(4, 0);
    ch.Face_Neighbour_List.assign(4, vector<int>());
    ch.Cell_Center_Distances.assign(4, 1.0);
    ch.Face_Areas.assign(4, 1.0);
    ch.Face_Normals.assign(8, 0.0);
    ch.Area = p.Area * 0.25;
    ch.Inv_Area = (ch.Area > 0.0) ? (1.0 / ch.Area) : 0.0;
    ch.del_t = p.del_t;

    // Map parent face normals: 0-left,1-bottom,2-right,3-top.
    for (int f = 0; f < 4; ++f)
    {
        ch.Face_Normals[2 * f + 0] = p.Face_Normals[2 * f + 0];
        ch.Face_Normals[2 * f + 1] = p.Face_Normals[2 * f + 1];
    }

    const double leftLen = (p.Face_Areas.size() > 0) ? p.Face_Areas[0] : 1.0;
    const double botLen = (p.Face_Areas.size() > 1) ? p.Face_Areas[1] : 1.0;
    const double rightLen = (p.Face_Areas.size() > 2) ? p.Face_Areas[2] : leftLen;
    const double topLen = (p.Face_Areas.size() > 3) ? p.Face_Areas[3] : botLen;
    const double vHalf = 0.25 * (leftLen + rightLen);  // half of avg vertical span
    const double hHalf = 0.25 * (botLen + topLen);     // half of avg horizontal span
    const double leftHalf = 0.5 * leftLen;
    const double rightHalf = 0.5 * rightLen;
    const double botHalf = 0.5 * botLen;
    const double topHalf = 0.5 * topLen;

    // Cell center offset using parent face normals (structured quad convention).
    const double nxL = p.Face_Normals[0], nyL = p.Face_Normals[1];
    const double nxB = p.Face_Normals[2], nyB = p.Face_Normals[3];
    const double nxR = p.Face_Normals[4], nyR = p.Face_Normals[5];
    const double nxT = p.Face_Normals[6], nyT = p.Face_Normals[7];
    const double tx = 0.5 * (nxR - nxL);
    const double ty = 0.5 * (nyR - nyL);
    const double sx = 0.5 * (nxT - nxB);
    const double sy = 0.5 * (nyT - nyB);
    const double scale = 0.25 * sqrt(std::max(p.Area, 1e-12));
    double ox = 0.0, oy = 0.0;
    if (childType == 0)      // SW
    {
        ox = -(tx + sx) * scale;
        oy = -(ty + sy) * scale;
    }
    else if (childType == 1) // SE
    {
        ox = (tx - sx) * scale;
        oy = (ty - sy) * scale;
    }
    else if (childType == 2) // NE
    {
        ox = (tx + sx) * scale;
        oy = (ty + sy) * scale;
    }
    else                     // NW
    {
        ox = (-tx + sx) * scale;
        oy = (-ty + sy) * scale;
    }
    if (ch.Cell_Center.size() < 2)
        ch.Cell_Center.assign(2, 0.0);
    ch.Cell_Center[0] = p.Cell_Center[0] + ox;
    ch.Cell_Center[1] = p.Cell_Center[1] + oy;

    // Initialize face lengths & boundary kinds from parent topology.
    if (childType == 0) // SW
    {
        ch.Face_Areas[0] = leftHalf;
        ch.Face_Areas[1] = botHalf;
        ch.Face_Areas[2] = vHalf;
        ch.Face_Areas[3] = hHalf;
        ch.Face_Boundary_Kind[0] = p.Face_Boundary_Kind.size() > 0 ? p.Face_Boundary_Kind[0] : 0;
        ch.Face_Boundary_Kind[1] = p.Face_Boundary_Kind.size() > 1 ? p.Face_Boundary_Kind[1] : 0;
    }
    else if (childType == 1) // SE
    {
        ch.Face_Areas[0] = vHalf;
        ch.Face_Areas[1] = botHalf;
        ch.Face_Areas[2] = rightHalf;
        ch.Face_Areas[3] = hHalf;
        ch.Face_Boundary_Kind[1] = p.Face_Boundary_Kind.size() > 1 ? p.Face_Boundary_Kind[1] : 0;
        ch.Face_Boundary_Kind[2] = p.Face_Boundary_Kind.size() > 2 ? p.Face_Boundary_Kind[2] : 0;
    }
    else if (childType == 2) // NE
    {
        ch.Face_Areas[0] = vHalf;
        ch.Face_Areas[1] = hHalf;
        ch.Face_Areas[2] = rightHalf;
        ch.Face_Areas[3] = topHalf;
        ch.Face_Boundary_Kind[2] = p.Face_Boundary_Kind.size() > 2 ? p.Face_Boundary_Kind[2] : 0;
        ch.Face_Boundary_Kind[3] = p.Face_Boundary_Kind.size() > 3 ? p.Face_Boundary_Kind[3] : 0;
    }
    else // NW
    {
        ch.Face_Areas[0] = leftHalf;
        ch.Face_Areas[1] = hHalf;
        ch.Face_Areas[2] = vHalf;
        ch.Face_Areas[3] = topHalf;
        ch.Face_Boundary_Kind[0] = p.Face_Boundary_Kind.size() > 0 ? p.Face_Boundary_Kind[0] : 0;
        ch.Face_Boundary_Kind[3] = p.Face_Boundary_Kind.size() > 3 ? p.Face_Boundary_Kind[3] : 0;
    }

    return ch;
}
} // namespace

void AMR_Compute_Gradient_Indicator()
{
    if (Gradient_Refinement_Indicator.size() != static_cast<size_t>(No_Physical_Cells))
        Gradient_Refinement_Indicator.assign(No_Physical_Cells, 0.0);

    V_D grad(2, 0.0);
    const double eps = 1e-14;
    int gradRho = 0;

    for (int i = 0; i < No_Physical_Cells; ++i)
    {
        if (!Cells[i].AMR_IsLeaf)
        {
            Gradient_Refinement_Indicator[i] = 0.0;
            continue;
        }
        const double h = sqrt(std::max(Cells[i].Area, eps));
        double ind = 0.0;
        Calculate_Gradient_At_Cell_Center(i, gradRho, grad);
        ind += vec_norm2(grad[0], grad[1]) * h;
        Gradient_Refinement_Indicator[i] = ind;
    }
}

void AMR_Tag_Cells_For_Split()
{
    AMR_Compute_Gradient_Indicator();
    for (int i = 0; i < No_Physical_Cells; ++i)
        Cells[i].Is_Splittable = false;

    vector<pair<double, int>> order;
    order.reserve(No_Physical_Cells);
    for (int i = 0; i < No_Physical_Cells; ++i)
    {
        if (!is_quad_leaf_cell(i))
            continue;
        if (Cells[i].has_Wall_Face || Cells[i].has_Inlet_Face || Cells[i].has_Exit_Face || Cells[i].has_Symmetry_Face)
            continue;
        if (Cells[i].AMR_Level >= 8) // hard safety cap
            continue;
        order.push_back({Gradient_Refinement_Indicator[i], i});
    }
    if (order.empty())
        return;

    sort(order.begin(), order.end(), [](const auto &a, const auto &b) { return a.first > b.first; });

    int nRefine = static_cast<int>(AMR_Max_Fraction * static_cast<double>(order.size()));
    if (nRefine < 1)
        nRefine = 1;
    if (nRefine > static_cast<int>(order.size()))
        nRefine = static_cast<int>(order.size());

    for (int k = 0; k < nRefine; ++k)
    {
        if (order[k].first >= AMR_Gradient_Threshold)
            Cells[order[k].second].Is_Splittable = true;
    }
}

bool AMR_Split_Tagged_Cells(int maxLevels)
{
    vector<int> targets;
    for (int i = 0; i < No_Physical_Cells; ++i)
    {
        if (!Cells[i].AMR_IsLeaf || !Cells[i].Is_Splittable)
            continue;
        if (Cells[i].numFaces != 4)
            continue;
        if (Cells[i].AMR_Level >= maxLevels)
            continue;
        targets.push_back(i);
    }
    if (targets.empty())
        return false;

    // Build children.
    for (int pidx : targets)
    {
        if (!Cells[pidx].AMR_IsLeaf)
            continue;
        const Cell p = Cells[pidx]; // snapshot: Cells vector may reallocate during insert
        vector<int> kids = Cells[pidx].AMR_Children;
        const bool can_reactivate = (kids.size() == 4 &&
                                     kids[0] >= 0 && kids[0] < No_Physical_Cells &&
                                     kids[1] >= 0 && kids[1] < No_Physical_Cells &&
                                     kids[2] >= 0 && kids[2] < No_Physical_Cells &&
                                     kids[3] >= 0 && kids[3] < No_Physical_Cells);

        if (!can_reactivate)
        {
            kids.assign(4, -1);
            for (int k = 0; k < 4; ++k)
            {
                const int cid = No_Physical_Cells;
                Cell ch = make_child_from_parent(p, k, pidx, cid);
                shift_index_vectors_after_insert(cid);
                ch.cellID = cid;
                Cells.insert(Cells.begin() + cid, ch);
                bool hasBoundary = false;
                for (int f = 0; f < 4; ++f)
                    hasBoundary = hasBoundary || (ch.Face_Boundary_Kind[f] != 0);
                append_child_storage(U_Cells[pidx], Primitive_Cells[pidx], hasBoundary);
                kids[k] = cid;
                No_Physical_Cells++;
                Total_No_Cells++;
            }
        }
        else
        {
            for (int k = 0; k < 4; ++k)
            {
                Cells[kids[k]].AMR_IsLeaf = true;
                Cells[kids[k]].Is_Splittable = false;
                U_Cells[kids[k]] = U_Cells[pidx];
                Primitive_Cells[kids[k]] = Primitive_Cells[pidx];
                Cells[kids[k]].del_t = Cells[pidx].del_t;
            }
        }

        // Sibling connectivity (0=SW,1=SE,2=NE,3=NW)
        // face order: 0-left,1-bottom,2-right,3-top
        Cells[kids[0]].Neighbours[2] = kids[1];
        Cells[kids[0]].Neighbours[3] = kids[3];
        Cells[kids[1]].Neighbours[0] = kids[0];
        Cells[kids[1]].Neighbours[3] = kids[2];
        Cells[kids[2]].Neighbours[0] = kids[3];
        Cells[kids[2]].Neighbours[1] = kids[1];
        Cells[kids[3]].Neighbours[1] = kids[0];
        Cells[kids[3]].Neighbours[2] = kids[2];

        // External coarse neighbors (kept coarse-fine for now).
        const int nL = p.Neighbours.size() > 0 ? p.Neighbours[0] : -1;
        const int nB = p.Neighbours.size() > 1 ? p.Neighbours[1] : -1;
        const int nR = p.Neighbours.size() > 2 ? p.Neighbours[2] : -1;
        const int nT = p.Neighbours.size() > 3 ? p.Neighbours[3] : -1;
        Cells[kids[0]].Neighbours[0] = nL;
        Cells[kids[3]].Neighbours[0] = nL;
        Cells[kids[0]].Neighbours[1] = nB;
        Cells[kids[1]].Neighbours[1] = nB;
        Cells[kids[1]].Neighbours[2] = nR;
        Cells[kids[2]].Neighbours[2] = nR;
        Cells[kids[2]].Neighbours[3] = nT;
        Cells[kids[3]].Neighbours[3] = nT;

        // Redirect neighboring leaf cells that previously pointed to parent.
        if (nL >= 0 && nL < No_Physical_Cells && Cells[nL].AMR_IsLeaf)
        {
            const int useChild = (Cells[nL].Cell_Center[1] > p.Cell_Center[1]) ? kids[3] : kids[0];
            remap_neighbor_reference(nL, pidx, useChild);
        }
        if (nR >= 0 && nR < No_Physical_Cells && Cells[nR].AMR_IsLeaf)
        {
            const int useChild = (Cells[nR].Cell_Center[1] > p.Cell_Center[1]) ? kids[2] : kids[1];
            remap_neighbor_reference(nR, pidx, useChild);
        }
        if (nB >= 0 && nB < No_Physical_Cells && Cells[nB].AMR_IsLeaf)
        {
            const int useChild = (Cells[nB].Cell_Center[0] > p.Cell_Center[0]) ? kids[1] : kids[0];
            remap_neighbor_reference(nB, pidx, useChild);
        }
        if (nT >= 0 && nT < No_Physical_Cells && Cells[nT].AMR_IsLeaf)
        {
            const int useChild = (Cells[nT].Cell_Center[0] > p.Cell_Center[0]) ? kids[2] : kids[3];
            remap_neighbor_reference(nT, pidx, useChild);
        }

        Cells[pidx].AMR_Children = kids;
        Cells[pidx].AMR_IsLeaf = false;
        Cells[pidx].Is_Splittable = false;
        Cells[pidx].hasBoundaryface = false;
    }
    return true;
}

void AMR_Tag_Cells_For_Merge()
{
    for (int i = 0; i < No_Physical_Cells; ++i)
        Cells[i].Is_Splittable = false;

    if (AMR_Coarsen_Threshold <= 0.0)
        return;

    for (int p = 0; p < No_Physical_Cells; ++p)
    {
        if (Cells[p].AMR_IsLeaf)
            continue;
        if (Cells[p].AMR_Level <= 0)
            continue;
        if (Cells[p].AMR_Children.size() != 4)
            continue;

        const vector<int> &kids = Cells[p].AMR_Children;
        bool all_leaf = true;
        bool all_valid = true;
        double max_ind = 0.0;
        for (const int k : kids)
        {
            if (k < 0 || k >= No_Physical_Cells)
            {
                all_valid = false;
                break;
            }
            if (!Cells[k].AMR_IsLeaf)
                all_leaf = false;
            if (k < static_cast<int>(Gradient_Refinement_Indicator.size()))
                max_ind = std::max(max_ind, Gradient_Refinement_Indicator[k]);
        }
        if (!all_valid || !all_leaf)
            continue;

        if (max_ind <= AMR_Coarsen_Threshold)
            Cells[p].Is_Splittable = true; // reuse tag as "merge parent"
    }
}

bool AMR_Merge_Tagged_Cells()
{
    bool changed = false;
    for (int p = 0; p < No_Physical_Cells; ++p)
    {
        if (!Cells[p].Is_Splittable)
            continue;
        if (Cells[p].AMR_IsLeaf || Cells[p].AMR_Children.size() != 4)
            continue;

        vector<int> kids = Cells[p].AMR_Children;
        bool valid = true;
        for (const int k : kids)
        {
            if (k < 0 || k >= No_Physical_Cells || !Cells[k].AMR_IsLeaf)
            {
                valid = false;
                break;
            }
        }
        if (!valid)
            continue;

        // Redirect external neighbors from children back to parent.
        for (const int k : kids)
        {
            const auto &N = Cells[k].Neighbours;
            for (int f = 0; f < static_cast<int>(N.size()); ++f)
            {
                const int neigh = N[f];
                bool is_sibling = false;
                for (const int s : kids)
                    if (neigh == s)
                        is_sibling = true;
                if (!is_sibling && neigh >= 0 && neigh < No_Physical_Cells)
                    remap_neighbor_reference_all(neigh, kids, p);
            }
        }

        // Restore coarse neighbors on parent from children boundary faces.
        Cells[p].Neighbours[0] = pick_first_valid_neighbor(Cells[kids[0]].Neighbours[0], Cells[kids[3]].Neighbours[0]);
        Cells[p].Neighbours[1] = pick_first_valid_neighbor(Cells[kids[0]].Neighbours[1], Cells[kids[1]].Neighbours[1]);
        Cells[p].Neighbours[2] = pick_first_valid_neighbor(Cells[kids[1]].Neighbours[2], Cells[kids[2]].Neighbours[2]);
        Cells[p].Neighbours[3] = pick_first_valid_neighbor(Cells[kids[3]].Neighbours[3], Cells[kids[2]].Neighbours[3]);

        // Restrict child states to parent by simple average.
        for (int c = 0; c < NUM_FLUX_COMPONENTS; ++c)
            U_Cells[p][c] = 0.25 * (U_Cells[kids[0]][c] + U_Cells[kids[1]][c] + U_Cells[kids[2]][c] + U_Cells[kids[3]][c]);
        Primitive_Cells[p] = Primitive_Cells[kids[0]];

        Cells[p].AMR_IsLeaf = true;
        Cells[p].Is_Splittable = false;
        for (const int k : kids)
        {
            Cells[k].AMR_IsLeaf = false;
            Cells[k].Is_Splittable = false;
        }
        changed = true;
    }
    return changed;
}

bool AMR_Adaptive_Step()
{
    AMR_Compute_Gradient_Indicator();

    AMR_Tag_Cells_For_Merge();
    int tagged_merge = 0;
    for (int i = 0; i < No_Physical_Cells; ++i)
        if (Cells[i].Is_Splittable && !Cells[i].AMR_IsLeaf)
            tagged_merge++;

    bool changed_merge = false;
    if (tagged_merge > 0)
    {
        cout << "AMR: " << tagged_merge << " parents tagged for 4->1 merge (threshold=" << AMR_Coarsen_Threshold << ")" << endl;
        changed_merge = AMR_Merge_Tagged_Cells();
        if (changed_merge)
            cout << "AMR: merge applied and parent cells reactivated." << endl;
    }

    AMR_Tag_Cells_For_Split();
    int tagged_split = 0;
    for (int i = 0; i < No_Physical_Cells; ++i)
        if (Cells[i].Is_Splittable && Cells[i].AMR_IsLeaf)
            tagged_split++;

    bool changed_split = false;
    if (tagged_split > 0)
    {
        cout << "AMR: " << tagged_split << " leaf cells tagged for 1->4 split (threshold=" << AMR_Gradient_Threshold << ")" << endl;
        changed_split = AMR_Split_Tagged_Cells(/*maxLevels*/ 2);
        if (changed_split)
            cout << "AMR: split applied and children activated. Active leaf set expanded." << endl;
    }
    return changed_merge || changed_split;
}
