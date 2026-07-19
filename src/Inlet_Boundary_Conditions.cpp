#include "definitions.h"
#include "Globals.h"
#include "Boundary_Conditions.h"
#include "Primitive_Computational.h"
#include "Utilities.h"
#include <cmath>
#include <iostream>

void Subsonic_Inlet_Condition(InletCondition &inletCond, V_I &Inlet_Cells_List)
{
	// The implementation is as described in the Blazek Text book Page 283
	// Pressure and Density are prescribed and the velocity is extrapolated from the interior cells
	double V_Magnitude = 0.0;
	int Cell_Index = 0, Ghost_Cell_Index = 0;
	if (Inlet_Cells_List.size() % 3 != 0)
	{
		std::cerr << "Error: Inlet_Cells_List size is not a multiple of 3." << std::endl;
		return;
	}
	for (unsigned int i = 0; i < Inlet_Cells_List.size(); i += 3)
	{
		Cell_Index = Inlet_Cells_List[i + 0]; // fetch the current cell index

		Ghost_Cell_Index = Inlet_Cells_List[i + 2]; // fetch the ghost cell index

		inletCond.u = Primitive_Cells[Cell_Index][1]; // extrapolate the velocity from the interior cells u
		inletCond.v = Primitive_Cells[Cell_Index][2]; // extrapolate the velocity from the interior cells v

		V_Magnitude = 0.5 * (inletCond.u * inletCond.u + inletCond.v * inletCond.v); // 0.5*V^2

		U_Cells[Ghost_Cell_Index][0] = inletCond.Rho;														  // rho
		U_Cells[Ghost_Cell_Index][1] = inletCond.Rho * inletCond.u;											  // rho*u
		U_Cells[Ghost_Cell_Index][2] = inletCond.Rho * inletCond.v;											  // rho*v
		U_Cells[Ghost_Cell_Index][3] = (inletCond.P / gamma_M_1) + inletCond.Rho * V_Magnitude * V_Magnitude; // P/(gamma-1) + 0.5*rho*V^2

		Calculate_Primitive_Variables(Cell_Index, U_Cells[Ghost_Cell_Index]);
		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		std::copy(Global_Primitive.begin(), Global_Primitive.end(), Primitive_Cells[Ghost_Cell_Index].begin());
	}
	// 	cout<<"******************* Subsonic Inlet Boundary Conditions***********************************\n";
}

void Supersonic_Inlet_Condition(InletCondition &inlet, V_I &Inlet_Cells_List)
{
	/* Blazek §8.4: all characteristics enter — prescribe freestream (do NOT
	   extrapolate velocity or mutate the global inletCond freestream state). */
	const double u_inf = inlet.u;
	const double v_inf = inlet.v;
	const double rho_inf = inlet.Rho;
	const double p_inf = inlet.P;
	const double ke = 0.5 * (u_inf * u_inf + v_inf * v_inf);

	if (!(rho_inf > 0.0) || !(p_inf > 0.0) || !std::isfinite(ke))
	{
		std::cerr << "Error: Supersonic_Inlet_Condition invalid freestream "
				  << "rho=" << rho_inf << " p=" << p_inf
				  << " u=" << u_inf << " v=" << v_inf << std::endl;
		return;
	}
	if (ke < 1.0e-24)
	{
		std::cerr << "Warning: Supersonic_Inlet_Condition freestream velocity is ~0 "
				  << "(u=" << u_inf << ", v=" << v_inf << "). "
				  << "Check Flow_Conditions.inlet_conditions in the case JSON.\n";
	}

	if (Inlet_Cells_List.size() % 3 != 0)
	{
		std::cerr << "Error: Inlet_Cells_List size is not a multiple of 3." << std::endl;
		return;
	}
	for (unsigned int i = 0; i < Inlet_Cells_List.size(); i += 3)
	{
		const int Ghost_Cell_Index = Inlet_Cells_List[i + 2];

		U_Cells[Ghost_Cell_Index][0] = rho_inf;
		U_Cells[Ghost_Cell_Index][1] = rho_inf * u_inf;
		U_Cells[Ghost_Cell_Index][2] = rho_inf * v_inf;
		U_Cells[Ghost_Cell_Index][3] = (p_inf / gamma_M_1) + rho_inf * ke;

		Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);
		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		for (unsigned int k = 0; k < Global_Primitive.size(); k++)
			Primitive_Cells[Ghost_Cell_Index][k] = Global_Primitive[k];
	}
}

// Super Sonic Conditions Prescribe all the conditions in the Ghost Cells
/*void Supersonic_Inlet_Boundary_Condition()
{
	int Cell_Index = 0, Ghost_Cell_Index = 0, Face_No = 0; // Face_Index=0,n1=0,n2=0;

	//	cout<<"Applying Supersonic Inlet Boundary Condition\t"<<endl;
	//	cout<<"Number of Inlet Cells\t"<<Inlet_Cells_List.size()/3<<endl;
	for (unsigned int i = 0; i < Inlet_Cells_List.size(); i += 3)
	{
		Cell_Index = Inlet_Cells_List[i + 0];
		Face_No = Inlet_Cells_List[i + 1];
		Ghost_Cell_Index = Inlet_Cells_List[i + 2];

		switch (Test_Case)
		{

		case 1:
			Pressure_Static_Inlet = 1.0;
			Inlet_Mach_No = 6.0;
			Rho_Static_Inlet = 1.4;

			V_1 = 0.0;
			V_2 = -6.0;

			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);

			break;
		case 4:
			Pressure_Static_Inlet = 1.0;
			Inlet_Mach_No = 3.0;
			Rho_Static_Inlet = 1.4;

			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;

			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);

			break;
		case 3:
			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 10:
			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;

		case 2:
			if (Face_No == 3)
			{
				if (Is_Viscous_Wall)
				{
					Pressure_Static_Inlet = 1.2473 * P_ref;
					Rho_Static_Inlet = 1.1706 * Rho_ref;
					V_1 = 2.0077 * cos(radian * 3.8131) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
					V_2 = -2.0077 * sin(radian * 3.8131) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
				}
				else
				{
					Pressure_Static_Inlet = 1.52819;
					Rho_Static_Inlet = 1.69997;
					V_1 = 2.61934;
					V_2 = -0.50633;
				}
			}
			else if (Face_No == 0)
			{
				if (Is_Viscous_Wall)
				{
					// Test Conditions for Visocus Case
					Pressure_Static_Inlet = P_ref;
					Rho_Static_Inlet = Rho_ref;
					V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
					V_2 = 0.0;
					if (Cell_Index >= 14000)
					{
						Pressure_Static_Inlet = 1.24729 * P_ref;
						Rho_Static_Inlet = 1.17069013 * Rho_ref;
						V_1 = 2.0077 * cos(radian * 3.8131) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
						V_2 = -2.0077 * sin(radian * 3.8131) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
					}
				}
				else
				{ // Test Conditions for Invisicd Case
					Pressure_Static_Inlet = (1.0 / 1.4);
					Rho_Static_Inlet = 1.0;
					V_1 = 2.9;
					V_2 = 0.0;
				}
			}
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 12:
			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			V_2 = 0.0;
			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			//					cout<<"Inlet Velocity at Boundary\t"<<V_1<<endl;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;

		case 11:
			Pressure_Static_Inlet = 1.0;
			Rho_Static_Inlet = 1.4;
			V_2 = 0.0;
			if (Cell_Index < 0.5 * No_Physical_Cells)
				V_1 = 2.0;
			else
				V_1 = 3.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 13:
			Pressure_Static_Inlet = 35.125;
			Rho_Static_Inlet = 7.208510638;
			V_1 = 4.43182;
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;

		case 14:
			// 			Conditions used by Maruti in his code Ic_Unsteady_Flow(id1, id2, jd1, jd2, 1, 22, 7.041132896, 1.4, 4.0779, 0.0, 0.0, 0.0, 30.05945, 1.0, x, y, cv);

			Pressure_Static_Inlet = 30.05945;
			Rho_Static_Inlet = 7.041132896;
			V_1 = 4.0779;
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 8:
			Pressure_Static_Inlet = P_ref;
			Rho_Static_Inlet = Rho_ref;
			V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;

		case 16:
			Pressure_Static_Inlet = 1.0 / 1.4;
			Rho_Static_Inlet = 1.25;
			V_1 = 1.22 * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 17:
			if (Is_Viscous_Wall)
			{
				Pressure_Static_Inlet = 2.12 * P_ref;
				Rho_Static_Inlet = (2.12 / 1.25469387) * Rho_ref;
				//						Pressure_Static_Inlet = 1.01761659*P_ref;
				//						Rho_Static_Inlet=1.68946061*Rho_ref;

				V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
				V_2 = 0.0;
			}
			else
			{
				Pressure_Static_Inlet = P_ref;
				Rho_Static_Inlet = Rho_ref;
				V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
				V_2 = 0.0;
			}
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 18: // Ic_Unsteady_Flow(id1, id2, jd1, jd2, 1, 22, 7.37561, 1.4, 4.8611, 0.0, 0.0, 0.0, 41.8333, 1.0, x, y, cv); //oed from Maruthis Code for OED
			Pressure_Static_Inlet = 41.8333;
			Rho_Static_Inlet = 7.37561;
			V_1 = 4.8611;
			V_2 = 0.0;
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 15:
			if (Face_No == 3)
			{
				if (Is_Viscous_Wall)
				{
					Pressure_Static_Inlet = 1.24729080 * P_ref;
					Rho_Static_Inlet = 1.17060913 * Rho_ref;
					V_1 = 2.00771568 * cos(radian * 3.81303835) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
					V_2 = -2.00771568 * sin(radian * 3.81303835) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
				}
				else
				{
					Pressure_Static_Inlet = 1.52819;
					Rho_Static_Inlet = 1.69997;
					V_1 = 2.61934;
					V_2 = -0.50633;
				}
			}
			else if (Face_No == 0)
			{
				if (Is_Viscous_Wall)
				{
					// Test Conditions for Visocus Case
					Pressure_Static_Inlet = P_ref;
					Rho_Static_Inlet = Rho_ref;
					V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
					V_2 = 0.0;
					if (Cell_Index >= 400 * 160)
					{
						Pressure_Static_Inlet = 1.24729080 * P_ref;
						Rho_Static_Inlet = 1.17060913 * Rho_ref;
						V_1 = 2.00771568 * cos(radian * 3.8131) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
						V_2 = -2.00771568 * sin(radian * 3.8131) * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
					}
				}
				else
				{ // Test Conditions for Invisicd Case
					Pressure_Static_Inlet = (1.0 / 1.4);
					Rho_Static_Inlet = 1.0;
					V_1 = 2.9;
					V_2 = 0.0;
				}
			}
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		case 19:
			Inlet_Mach_No = 3.0;
			Pressure_Static_Inlet = 1.0 / (gamma * Inlet_Mach_No * Inlet_Mach_No);
			Rho_Static_Inlet = 1.0;

			if (Cell_Index == ((nx_c - 1) * 0.5 * (ny_c - 1)))
			{
				V_1 = 0.0;
				V_2 = 0.0;
			}
			else
			{
				V_1 = Inlet_Mach_No * sqrt(gamma * Pressure_Static_Inlet / Rho_Static_Inlet);
				V_2 = 0.0;
			}
			V_Magnitude_Inlet = sqrt(V_1 * V_1 + V_2 * V_2);
			break;
		}
		//   		cout<<Cell_Index<<"\t"<<Pressure_Static_Inlet<<"\t"<<Rho_Static_Inlet<<"\t"<<Temperature_Static_Inlet<<"\t"<<V_1<<"\t"<<V_2<<"\t"<<V_Magnitude<<endl;
		U_Cells[Ghost_Cell_Index][0] = Rho_Static_Inlet;
		U_Cells[Ghost_Cell_Index][1] = Rho_Static_Inlet * V_1;
		U_Cells[Ghost_Cell_Index][2] = Rho_Static_Inlet * V_2;
		U_Cells[Ghost_Cell_Index][3] = ((Pressure_Static_Inlet / gamma_M_1) + Rho_Static_Inlet * 0.5 * V_Magnitude_Inlet * V_Magnitude_Inlet); // (P/(gamma-1) + 0.5*rho*|V|^2)

		Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);
		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		for (unsigned int i = 0; i < Global_Primitive.size(); i++)
			Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
	}
}

// Super Sonic Conditions Prescribe all the conditions in the Ghost Cells
void Supersonic_Inlet_Boundary_Condition(double &Pressure_Static_Inlet, double &Rho_Static_Inlet, double &u_Inlet, double &v_Inlet, double &Inlet_Mach_No)
{
	int Cell_Index = 0, Ghost_Cell_Index = 0, Face_No = 0; // Face_Index=0,n1=0,n2=0;
	double VMagnitude = 0.0;

	//	cout<<"Applying Supersonic Inlet Boundary Condition\t"<<endl;
	//	cout<<"Number of Inlet Cells\t"<<Inlet_Cells_List.size()/3<<endl;
	VMagnitude = sqrt(u_Inlet * u_Inlet + v_Inlet * v_Inlet);
	//   		cout<<Cell_Index<<"\t"<<Pressure_Static_Inlet<<"\t"<<Rho_Static_Inlet<<"\t"<<Temperature_Static_Inlet<<"\t"<<V_1<<"\t"<<V_2<<"\t"<<V_Magnitude<<endl;
	U_Cells[Ghost_Cell_Index][0] = Rho_Static_Inlet;
	U_Cells[Ghost_Cell_Index][1] = Rho_Static_Inlet * u_Inlet;
	U_Cells[Ghost_Cell_Index][2] = Rho_Static_Inlet * v_Inlet;
	U_Cells[Ghost_Cell_Index][3] = ((Pressure_Static_Inlet / gamma_M_1) + Rho_Static_Inlet * 0.5 * VMagnitude * VMagnitude); // (P/(gamma-1) + 0.5*rho*|V|^2)

	Calculate_Primitive_Variables(Ghost_Cell_Index, U_Cells[Ghost_Cell_Index]);
	Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
	for (unsigned int i = 0; i < Global_Primitive.size(); i++)
		Primitive_Cells[Ghost_Cell_Index][i] = Global_Primitive[i];
}

/*void Supersonic_Inlet_Boundary_Condition(double & Pressure_Static_Inlet,double & Rho_Static_Inlet,double & u_Inlet,double & v_Inlet,double & Inlet_Mach_No)
{
	int Cell_Index=0,Ghost_Cell_Index=0,Face_No=0;//Face_Index=0,n1=0,n2=0;

//  	  cout<<"Applying Supersonic Inlet Boundary Condition\t"<<endl;
//           cout<<"Number of Inlet Cells\t"<<Inlet_Cells_List.size()/3<<endl;
	for(unsigned int i=0;i<Inlet_Cells_List.size();i+=3)
	{
		Cell_Index = Inlet_Cells_List[i+0];
		Face_No = Inlet_Cells_List[i+1];
		Ghost_Cell_Index=Inlet_Cells_List[i+2];

		V_Magnitude_Inlet = sqrt(u_Inlet*u_Inlet + v_Inlet*v_Inlet);
//   		cout<<Cell_Index<<"\t"<<Pressure_Static_Inlet<<"\t"<<Rho_Static_Inlet<<"\t"<<Temperature_Static_Inlet<<"\t"<<V_1<<"\t"<<V_2<<"\t"<<V_Magnitude<<endl;
		U_Cells[Ghost_Cell_Index][0]	=	Rho_Static_Inlet;
		U_Cells[Ghost_Cell_Index][1]	=	Rho_Static_Inlet*V_1;
		U_Cells[Ghost_Cell_Index][2]	=	Rho_Static_Inlet*V_2;
		U_Cells[Ghost_Cell_Index][3]	=	((Pressure_Static_Inlet/(gamma-1.0)) + Rho_Static_Inlet*0.5*V_Magnitude_Inlet*V_Magnitude_Inlet); // (P/(gamma-1) + 0.5*rho*|V|^2)

		Calculate_Primitive_Variables(Ghost_Cell_Index,U_Cells[Ghost_Cell_Index]);
		Vector_Reset(Primitive_Cells[Ghost_Cell_Index]);
		for(unsigned int i=0;i<Global_Primitive.size();i++)
			Primitive_Cells[Ghost_Cell_Index][i]=Global_Primitive[i];
	}
}*/