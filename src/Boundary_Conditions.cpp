#include "definitions.h"
#include "Boundary_Conditions.h"
#include "Globals.h"

string BCFileName, InitCondFileName;

// InletCondition inletCond;
// ExitCondition exitCond;
// InitialCondition initCond;

// This function groups all the other boundary conditions
void Apply_Boundary_Conditions()
{
	// 	cout<<"Applying Boundary Conditions\n";
	switch (Is_Inlet_SubSonic)
	{
	case true:
		// Subsonic_Inlet_Boundary_Condition();
		Subsonic_Inlet_Condition(inletCond, Inlet_Cells_List);
		break;
	case false:
		// Supersonic_Inlet_Boundary_Condition();
		Supersonic_Inlet_Condition(inletCond, Inlet_Cells_List);
		break;
	}
	switch (Is_Exit_SubSonic)
	{
	case true:
		for (unsigned int i = 0; i < Exit_Cells_List.size(); i += 3)
		{
			// Subsonic_Exit_Boundary_Condition(i);
			Subsonic_Exit_Condition(exitCond, Exit_Cells_List);
		}
		break;
	case false:
		for (unsigned int i = 0; i < Exit_Cells_List.size(); i += 3)
		{
			// Supersonic_Exit_Boundary_Condition(i);
			Supersonic_Exit_Condition(exitCond, Exit_Cells_List);
		}
		break;
	}
	// 	cout<<"Inlet and Exit applied"<<endl;
	switch (Is_Viscous_Wall)
	{
	case true:
		//		cout<<"Viscous wall enabled"<<endl;
		Viscous_Wall_Boundary_Condition();
		break;
	case false:
		InViscid_Wall_Boundary_Condition();
		break;
	}
	// 	cout<<"Boundary Conditions Properly applied"<<endl;
	switch (has_Symmetry_BC)
	{
	case true:
		Symmetry_Boundary_Condition();
		break;
	}
}

