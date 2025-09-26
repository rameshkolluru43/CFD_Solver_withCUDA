#include "definitions.h"
#include "Globals.h"
#include "Error_Update.h"
#include "Utilities.h"
#include "Primitive_Computational.h"

V_D Error(4, 0.0);

// This Function estimates the error from the iteratons
void Estimate_Error()
{
	V_D temp(4, 0.0);
	Vector_Reset(Error);
	for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		for (int i = 0; i < NUM_FLUX_COMPONENTS; i++)
		{
			Vector_Reset(temp);
			if (Cells_DelU[Cell_Index][i] < 1e-9)
				temp[i] = fabs(Cells_DelU[Cell_Index][i]);
			else
				temp[i] = (fabs(Cells_DelU[Cell_Index][i]) / U_Cells[Cell_Index][i]);
			Error[i] += temp[i] * temp[i];
		}
	}
	for (int i = 0; i < 4; i++)
		Error[i] = sqrt(Error[i]);
}

// This function updates all the cells values after finding the error.
void Update()
{
	int Cell_Index;
	//   cout<<"Updating the Function"<<endl;
	for (Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		U_Cells[Cell_Index][0] += Cells_DelU[Cell_Index][0];
		U_Cells[Cell_Index][1] += Cells_DelU[Cell_Index][1];
		U_Cells[Cell_Index][2] += Cells_DelU[Cell_Index][2];
		U_Cells[Cell_Index][3] += Cells_DelU[Cell_Index][3];
		Calculate_Primitive_Variables(Cell_Index, U_Cells[Cell_Index], Primitive_Cells[Cell_Index]);
	}
}

// This function updates a given cell
void Update(int &Cell_Index, int &Update_Type)
{
	switch (Update_Type)
	{
	case 1:
		Calculate_Primitive_Variables(Cell_Index, U_Cells_RK_1[Cell_Index]);
		Vector_Reset(Primitive_Cells[Cell_Index]);
		Primitive_Cells[Cell_Index] = Global_Primitive;
		break;
	case 2:
		Calculate_Primitive_Variables(Cell_Index, U_Cells_RK_2[Cell_Index]);
		Vector_Reset(Primitive_Cells[Cell_Index]);
		Primitive_Cells[Cell_Index] = Global_Primitive;
		break;
	}
}

void Update(const int &Cell_Index, V_D &U)
{
	// 	Print(U);
	Calculate_Primitive_Variables(Cell_Index, U);
	Vector_Reset(Primitive_Cells[Cell_Index]);
	Primitive_Cells[Cell_Index] = Global_Primitive;
}

void Update(const int &Cell_Index)
{
	// 	Print(U);
	Calculate_Primitive_Variables(Cell_Index, U_Cells[Cell_Index]);
	Vector_Reset(Primitive_Cells[Cell_Index]);
	Primitive_Cells[Cell_Index] = Global_Primitive;
}
