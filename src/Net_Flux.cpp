#include "definitions.h"
#include "Globals.h"
#include "Timestep.h"
#include "Flux.h"
#include "Limiter.h"
#include "Grid.h"

#ifdef USE_CUDA
#include "../CUDA_KERNELS/MOVERS_RICCA_Flux_Cuda.h"
#endif

namespace
{
bool Is_AMR_Transition_Face(const int owner, const int neigh)
{
    if (!Enable_AMR)
        return false;
    if (owner < 0 || neigh < 0 || owner >= static_cast<int>(Cells.size()) || neigh >= static_cast<int>(Cells.size()))
        return false;

    if (Cells[owner].AMR_Level != Cells[neigh].AMR_Level)
        return true;
    if (Cells[owner].AMR_Parent == neigh || Cells[neigh].AMR_Parent == owner)
        return true;
    return false;
}

void Dissipation_Fallback_1O(const int &Current_Cell_No, int &neighbour, const int &face)
{
    switch (Dissipation_Type)
    {
    case 1:
        LLF(Current_Cell_No, neighbour, face);
        break;
    case 2:
        MOVERS(Current_Cell_No, neighbour, face);
        break;
    case 3:
        ROE(Current_Cell_No, neighbour, face);
        break;
    case 4:
        RICCA(Current_Cell_No, neighbour, face);
        break;
    case 5:
        MOVERS_NWSC(Current_Cell_No, neighbour, face);
        break;
    default:
        LLF(Current_Cell_No, neighbour, face);
        break;
    }
}
} // namespace

// Cell Net Flux evaluates both the average flux and Dissipative flux of a given cell
void Calculate_Flux_For_All_Faces(int &Current_Cell_No, void (*Dissipation_Function)(const int &, int &, const int &))
{
#ifdef DEBUG
    std::cerr << "Entered Flux Calculation (generic faces)" << std::endl;
#endif

    const int nFaces = (Cells[Current_Cell_No].numFaces > 0)
                           ? Cells[Current_Cell_No].numFaces
                           : static_cast<int>(Cells[Current_Cell_No].Face_Areas.size());
    const int nNeigh = static_cast<int>(Cells[Current_Cell_No].Neighbours.size());

    for (int face = 0; face < nFaces; face++)
    {
        if (face >= nNeigh)
        {
            // No neighbor recorded for this face; skip
            continue;
        }

        int neighbour = Cells[Current_Cell_No].Neighbours[face];
        if (neighbour < 0)
            continue;

        bool isBoundaryFace = false;
        if (Current_Cell_No < static_cast<int>(Cells_Face_Boundary_Type.size()) &&
            face < static_cast<int>(Cells_Face_Boundary_Type[Current_Cell_No].size()))
        {
            isBoundaryFace = Cells_Face_Boundary_Type[Current_Cell_No][face];
        }

        const bool amr_transition_face = Is_AMR_Transition_Face(Current_Cell_No, neighbour);
        if (Is_Second_Order && !amr_transition_face)
        {
            V_D muscl_dU_scratch(4, 0.0);
            Second_Order_Limiter(Current_Cell_No, face, muscl_dU_scratch);
            Calculate_Face_Average_Flux_MUSCL(Current_Cell_No, neighbour, face, isBoundaryFace);
            Dissipation_Function(Current_Cell_No, neighbour, face);
        }
        else
        {
            Calculate_Face_Average_Flux(Current_Cell_No, neighbour, face, isBoundaryFace);
            if (Is_Second_Order && amr_transition_face)
                Dissipation_Fallback_1O(Current_Cell_No, neighbour, face);
            else
                Dissipation_Function(Current_Cell_No, neighbour, face);
        }
        for (int v = 0; v < NUM_FLUX_COMPONENTS; v++)
        {
            Cells_Net_Flux[Current_Cell_No][v] += Average_Convective_Flux[v] - Dissipative_Flux[v];
        }
    }
}

void Evaluate_Cell_Net_Flux_1O()
{

#ifdef DEBUG
    std::cerr << "Entered 1st Order Flux Calculation" << std::endl;
#endif
#ifdef USE_CUDA
    if (Evaluate_Cell_Net_Flux_CUDA_Movers_Ricca(false))
        return;
#endif
    // Parallelize the loop with prgma
    // Loop over all physical cells
    // This function iterates over all physical cells and calculates the net flux for each cell
    // based on the specified dissipation type. It initializes the net flux for each cell to zero,
    // applies the appropriate flux calculation method for all faces of the cell, and evaluates
    // the time step for the current cell.
    // The dissipation type determines which flux calculation method is used:
    // 1: LLF (Local Lax-Friedrichs)
    // 2: MOVERS
    // 3: ROE (Roe's method)
    // 4: RICCA (Ricca's method)
    // 5: MOVERS_NWSC (MOVERS with no wall shear correction)
    // The function also evaluates the time step for each cell after computing the flux.
    // The function supports multiple dissipation types, each corresponding to a specific
    // flux calculation method.
    V_I leafCells;
    Build_Leaf_Cell_List(leafCells);
    for (int idx = 0; idx < static_cast<int>(leafCells.size()); idx++)
    {
        int Current_Cell_No = leafCells[idx];
        for (int i = 0; i < NUM_FLUX_COMPONENTS; i++)
        {
            Cells_Net_Flux[Current_Cell_No][i] = 0.0;
        }

        switch (Dissipation_Type)
        {
        case 1:
            Calculate_Flux_For_All_Faces(Current_Cell_No, LLF);
            break;
        case 2:
            Calculate_Flux_For_All_Faces(Current_Cell_No, MOVERS);
            break;
        case 3:
            Calculate_Flux_For_All_Faces(Current_Cell_No, ROE);
            break;
        case 4:
            Calculate_Flux_For_All_Faces(Current_Cell_No, RICCA);
            break;
        case 5:
            Calculate_Flux_For_All_Faces(Current_Cell_No, MOVERS_NWSC);
            break;
        }

        Evaluate_Time_Step(Current_Cell_No);
    }
}

void Evaluate_Cell_Net_Flux_2O()
{
#ifdef DEBUG
    std::cerr << "Entered 2nd Order Flux Calculation" << std::endl;
#endif
#ifdef USE_CUDA
    /* 2nd order uses MUSCL reconstruction on CPU; CUDA path is 1st order only. */
    if (Evaluate_Cell_Net_Flux_CUDA_Movers_Ricca(true))
        return;
#endif
    V_I leafCells;
    Build_Leaf_Cell_List(leafCells);
    for (int idx = 0; idx < static_cast<int>(leafCells.size()); idx++)
    {
        int Current_Cell_No = leafCells[idx];
        for (int k = 0; k < NUM_FLUX_COMPONENTS; k++)
        {
            Cells_Net_Flux[Current_Cell_No][k] = 0.0;
        }

        switch (Dissipation_Type)
        {
        case 1:
            Calculate_Flux_For_All_Faces(Current_Cell_No, LLF_2O);
            break;
        case 2:
            Calculate_Flux_For_All_Faces(Current_Cell_No, MOVERS_2O);
            break;
        case 3:
            Calculate_Flux_For_All_Faces(Current_Cell_No, ROE_2O);
            break;
        case 4:
            Calculate_Flux_For_All_Faces(Current_Cell_No, RICCA_2O);
            break;
        case 5:
            Calculate_Flux_For_All_Faces(Current_Cell_No, MOVERS_NWSC_2O);
            break;
        }
        Evaluate_Time_Step(Current_Cell_No);
    }
}