
#include "../Basic_Files/definitions.h"
/**
 * @brief Evaluates flux on a given face.
 * This function calculates the flux for a given face by considering the current cell number,
 * its neighboring cell corresponding to the face, the face number, and whether the face is a boundary face.
 *
 * @param Current_Cell_No The number of the current cell.
 * @param Neighbour_Cell_No The number of the neighboring cell corresponding to the face.
 * @param Face_No The face number for which the flux is being calculated.
 * @param Is_Boundary_Face A boolean indicating whether the face is a boundary face.
 */

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

// This function calculates the flux for a given face by considering the current cell number, its neighboring cell corresponding to the face, the face number, and whether the face is a boundary face.

/// @brief Evaluates the first-order net flux for all cells in the computational domain.

__global__ void Kernel_Evaluate_Cell_Net_Flux_1O()
{
#ifdef DEBUG
    std::cerr << "Entered 1st Order Flux Calculation" << std::endl;
#endif
}

void Evaluate_Cell_Net_Flux_1O()
{

#ifdef DEBUG
    std::cerr << "Entered 1st Order Flux Calculation" << std::endl;
#endif
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
    std::cerr << "Entered 1st Order Flux Calculation" << std::endl;
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
