#include "definitions.h"
#include "Globals.h"
#include "Timestep.h"
#include "Flux.h"

// Cell Net Flux evaluates both the average flux and Dissipative flux of a given cell
void Calculate_Flux_For_All_Faces(int &Current_Cell_No, void (*Dissipation_Function)(const int &, int &, const int &))
{
    int Neighbour_1 = 0, Neighbour_2 = 0, Neighbour_3 = 0, Neighbour_4 = 0; // Indicates the numbers to neighbours of the cell
#ifdef DEBUG
    std::cerr << "Entered 2nd Order Flux Calculation" << std::endl;
#endif
    Neighbour_1 = Cells[Current_Cell_No].Neighbours[0]; //(i-1,j,k)
    Neighbour_2 = Cells[Current_Cell_No].Neighbours[1]; //(i,j-1,k)
    Neighbour_3 = Cells[Current_Cell_No].Neighbours[2]; //(i+1,j,k)
    Neighbour_4 = Cells[Current_Cell_No].Neighbours[3]; //(i,j+1,k)

    Calculate_Face_Average_Flux(Current_Cell_No, Neighbour_1, Face_0, Cells_Face_Boundary_Type[Current_Cell_No][Face_0]);
    Dissipation_Function(Current_Cell_No, Neighbour_1, Face_0);
    Cells_Net_Flux[Current_Cell_No][0] += Average_Convective_Flux[0] - Dissipative_Flux[0];
    Cells_Net_Flux[Current_Cell_No][1] += Average_Convective_Flux[1] - Dissipative_Flux[1];
    Cells_Net_Flux[Current_Cell_No][2] += Average_Convective_Flux[2] - Dissipative_Flux[2];
    Cells_Net_Flux[Current_Cell_No][3] += Average_Convective_Flux[3] - Dissipative_Flux[3];

    Calculate_Face_Average_Flux(Current_Cell_No, Neighbour_2, Face_1, Cells_Face_Boundary_Type[Current_Cell_No][Face_1]);
    Dissipation_Function(Current_Cell_No, Neighbour_2, Face_1);
    Cells_Net_Flux[Current_Cell_No][0] += Average_Convective_Flux[0] - Dissipative_Flux[0];
    Cells_Net_Flux[Current_Cell_No][1] += Average_Convective_Flux[1] - Dissipative_Flux[1];
    Cells_Net_Flux[Current_Cell_No][2] += Average_Convective_Flux[2] - Dissipative_Flux[2];
    Cells_Net_Flux[Current_Cell_No][3] += Average_Convective_Flux[3] - Dissipative_Flux[3];

    Calculate_Face_Average_Flux(Current_Cell_No, Neighbour_3, Face_2, Cells_Face_Boundary_Type[Current_Cell_No][Face_2]);
    Dissipation_Function(Current_Cell_No, Neighbour_3, Face_2);
    Cells_Net_Flux[Current_Cell_No][0] += Average_Convective_Flux[0] - Dissipative_Flux[0];
    Cells_Net_Flux[Current_Cell_No][1] += Average_Convective_Flux[1] - Dissipative_Flux[1];
    Cells_Net_Flux[Current_Cell_No][2] += Average_Convective_Flux[2] - Dissipative_Flux[2];
    Cells_Net_Flux[Current_Cell_No][3] += Average_Convective_Flux[3] - Dissipative_Flux[3];

    Calculate_Face_Average_Flux(Current_Cell_No, Neighbour_4, Face_3, Cells_Face_Boundary_Type[Current_Cell_No][Face_3]);
    Dissipation_Function(Current_Cell_No, Neighbour_4, Face_3);
    Cells_Net_Flux[Current_Cell_No][0] += Average_Convective_Flux[0] - Dissipative_Flux[0];
    Cells_Net_Flux[Current_Cell_No][1] += Average_Convective_Flux[1] - Dissipative_Flux[1];
    Cells_Net_Flux[Current_Cell_No][2] += Average_Convective_Flux[2] - Dissipative_Flux[2];
    Cells_Net_Flux[Current_Cell_No][3] += Average_Convective_Flux[3] - Dissipative_Flux[3];
}

void Evaluate_Cell_Net_Flux_1O()
{

#ifdef DEBUG
    std::cerr << "Entered 1st Order Flux Calculation" << std::endl;
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
    for (int Current_Cell_No = 0; Current_Cell_No < No_Physical_Cells; Current_Cell_No++)
    {

        for (int i = 0; i < 4; i++)
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
    for (int Current_Cell_No = 0; Current_Cell_No < No_Physical_Cells; Current_Cell_No++)
    {
        for (int k = 0; k < 4; k++)
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