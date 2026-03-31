#include "definitions.h"
#include "Globals.h"
#include "Error_Update.h"
#include "Grid.h"
#include "Utilities.h"
#include "Primitive_Computational.h"
#ifdef _OPENMP
#include <omp.h>
#endif

V_D Error(4, 0.0);

// This Function estimates the error from the iteratons
void Estimate_Error()
{
	double e0 = 0.0, e1 = 0.0, e2 = 0.0, e3 = 0.0;
	V_I leafCells;
	Build_Leaf_Cell_List(leafCells);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : e0, e1, e2, e3) if (leafCells.size() > 1024)
#endif
	for (int li = 0; li < static_cast<int>(leafCells.size()); li++)
	{
		int Cell_Index = leafCells[li];
		for (int i = 0; i < NUM_FLUX_COMPONENTS; i++)
		{
			double t;
			if (Cells_DelU[Cell_Index][i] < 1e-9)
				t = fabs(Cells_DelU[Cell_Index][i]);
			else
				t = (fabs(Cells_DelU[Cell_Index][i]) / U_Cells[Cell_Index][i]);
			const double ts = t * t;
			if (i == 0)
				e0 += ts;
			else if (i == 1)
				e1 += ts;
			else if (i == 2)
				e2 += ts;
			else
				e3 += ts;
		}
	}
	Error[0] = sqrt(e0);
	Error[1] = sqrt(e1);
	Error[2] = sqrt(e2);
	Error[3] = sqrt(e3);
}

// This function updates all the cells values after finding the error.
void Update()
{
	int Cell_Index;
	V_I leafCells;
	Build_Leaf_Cell_List(leafCells);
	//   cout<<"Updating the Function"<<endl;
	for (int li = 0; li < static_cast<int>(leafCells.size()); li++)
	{
		Cell_Index = leafCells[li];
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
