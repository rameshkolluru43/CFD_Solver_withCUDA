#include "definitions.h"
#include "Globals.h"
#include "Flux.h"

void Smooth_Residuals()
{
	int i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i6 = 0;
	double omega = 0.5;
	/*for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
	{
		i1 = Cell_Neighbours[Cell_No][1];
		i2 = Cell_Neighbours[Cell_No][2];
		i3 = Cell_Neighbours[Cell_No][3];
		i4 = Cell_Neighbours[Cell_No][4];
		i5 = Cell_Neighbours[Cell_No][5];
		i6 = Cell_Neighbours[Cell_No][6];
		for (int i = 0; i < 5; i++)
			Smoothed_Cells_Net_Flux[Cell_No][i] = ((1.0 - omega) / 6.0) * (Cells_Net_Flux[i1][i] + Cells_Net_Flux[i2][i] + Cells_Net_Flux[i3][i] + Cells_Net_Flux[i4][i] + Cells_Net_Flux[i5][i] + Cells_Net_Flux[i6][i]) +
												  omega * Cells_Net_Flux[Cell_No][i];
	}
	for (int Cell_No = 0; Cell_No < No_Physical_Cells; Cell_No++)
	{
		for (int i = 0; i < 5; i++)
			Cells_Net_Flux[Cell_No][i] = Smoothed_Cells_Net_Flux[Cell_No][i];
	}*/
}
