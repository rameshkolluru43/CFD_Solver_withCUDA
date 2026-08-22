
#include "../Basic_Function_Files/headers.hpp"

void Cell::Velocity_GradatCenter()
{
	u_GradatCenter.Clear();
	v_GradatCenter.Clear();
	w_GradatCenter.Clear();
	T_GradatCenter.Clear();
	u_GradatCenter = Front_Face.Get_Face_undS() + Back_Face.Get_Face_undS() + Left_Face.Get_Face_undS() +
					 Right_Face.Get_Face_undS() + Top_Face.Get_Face_undS();

	v_GradatCenter = Front_Face.Get_Face_vndS() + Back_Face.Get_Face_vndS() + Left_Face.Get_Face_vndS() +
					 Right_Face.Get_Face_vndS() + Top_Face.Get_Face_vndS();

	w_GradatCenter = Front_Face.Get_Face_wndS() + Back_Face.Get_Face_wndS() + Left_Face.Get_Face_wndS() +
					 Right_Face.Get_Face_wndS() + Top_Face.Get_Face_wndS();

	T_GradatCenter = Front_Face.Get_Face_TndS() + Back_Face.Get_Face_TndS() + Left_Face.Get_Face_TndS() +
					 Right_Face.Get_Face_TndS() + Top_Face.Get_Face_TndS();
	if (No_of_Faces == 6)
	{
		u_GradatCenter += Bottom_Face.Get_Face_undS();
		v_GradatCenter += Bottom_Face.Get_Face_vndS();
		w_GradatCenter += Bottom_Face.Get_Face_wndS();
		T_GradatCenter += Bottom_Face.Get_Face_TndS();
	}
	// calculating velocity gradients  and Temperature for determining stresses on faces
	u_GradatCenter *= inv_vol;
	v_GradatCenter *= inv_vol;
	w_GradatCenter *= inv_vol;
	T_GradatCenter *= inv_vol;
}

void Cell::Cal_Stresses()
{
	// 	cout<<"In Stresses Function\t"<<self_index<<endl;
	Viscous_Flux.at(0) = 0.0;
	Viscous_Flux.at(1) = 0.0;
	Viscous_Flux.at(2) = 0.0;
	Viscous_Flux.at(3) = 0.0;
	Viscous_Flux.at(4) = 0.0;
	Viscous_Flux = Front_Face.Cal_Viscous_Stress_On_Face(this, Front_Cell);
	Add_Fluxes(Viscous_Flux, 1);
	Viscous_Flux = Back_Face.Cal_Viscous_Stress_On_Face(this, Back_Cell);
	Add_Fluxes(Viscous_Flux, 1);
	Viscous_Flux = Left_Face.Cal_Viscous_Stress_On_Face(this, Left_Cell);
	Add_Fluxes(Viscous_Flux, 1);
	Viscous_Flux = Right_Face.Cal_Viscous_Stress_On_Face(this, Right_Cell);
	Add_Fluxes(Viscous_Flux, 1);
	Viscous_Flux = Top_Face.Cal_Viscous_Stress_On_Face(this, Top_Cell);
	Add_Fluxes(Viscous_Flux, 1);
	if (No_of_Faces == 6)
	{
		Viscous_Flux = Bottom_Face.Cal_Viscous_Stress_On_Face(this, Bottom_Cell);
		Add_Fluxes(Viscous_Flux, 1);
	}
}

/**
 * @file Cell_Viscous_functions.cpp
 * @brief Contains the implementation of viscous stress calculations on a face.
 *
 *
 * @brief Calculates the viscous stress on a face based on the properties of the current cell and its neighbor.
 *
 * This function computes the viscous stress tensor components and fluxes on a face using the velocity and temperature
 * gradients of the current and neighboring cells. It also accounts for heat flux contributions.
 *
 * @param Current_Cell Pointer to the current cell.
 * @param Neighbour_Cell Pointer to the neighboring cell. This can be a ghost cell.
 * @return A constant reference to a vector containing the calculated flux values:
 *         - Flux[0]: Placeholder (always 0.0).
 *         - Flux[1]: Viscous stress flux in the x-direction.
 *         - Flux[2]: Viscous stress flux in the y-direction.
 *         - Flux[3]: Viscous stress flux in the z-direction.
 *         - Flux[4]: Total energy flux including heat flux contributions.
 *
 * @note If the neighboring cell is a ghost cell, the gradients from the current cell are used directly.
 *       Otherwise, the gradients are averaged between the current and neighboring cells.
 *
 * @details The function performs the following steps:
 *          1. Initializes the flux vector to zero.
 *          2. Determines the average velocity and temperature gradients on the face.
 *          3. Computes the components of the viscous stress tensor (Tij).
 *          4. Calculates the flux contributions due to viscous stresses and heat flux.
 *
 * @warning Ensure that the input cells are properly initialized and that the neighbor cell is correctly identified
 *          as either a ghost cell or a real cell.
 *
 * const std::vector<double> &Face::Cal_Viscous_Stress_On_Face(const Cell *const Current_Cell, const Cell *const Neighbour_Cell);*/

const vector<double> &Face::Cal_Viscous_Stress_On_Face(const Cell *const Current_Cell, const Cell *const Neighbour_Cell)
{
	double temp, T11, T12, T13, T21, T22, T23, T31, T32, T33, u11, u12, u13, u21, u22, u23, u31, u32, u33;
	Flux.at(0) = 0.0;
	Flux.at(1) = 0.0;
	Flux.at(2) = 0.0;
	Flux.at(3) = 0.0;
	Flux.at(4) = 0.0;
	if (Neighbour_Cell->Is_Ghost_Cell())
	{
		Avg_U_GradonFace = Current_Cell->Get_Grad_UatCenter();
		Avg_V_GradonFace = Current_Cell->Get_Grad_VatCenter();
		Avg_W_GradonFace = Current_Cell->Get_Grad_WatCenter();
		Avg_T_GradonFace = Current_Cell->Get_Grad_TatCenter();
	}
	else
	{
		Avg_U_GradonFace = (Current_Cell->Get_Grad_UatCenter() + Neighbour_Cell->Get_Grad_UatCenter()) * 0.5; // ui,j
		Avg_V_GradonFace = (Current_Cell->Get_Grad_VatCenter() + Neighbour_Cell->Get_Grad_VatCenter()) * 0.5; // vi,j
		Avg_W_GradonFace = (Current_Cell->Get_Grad_WatCenter() + Neighbour_Cell->Get_Grad_WatCenter()) * 0.5; // wi,j
		Avg_T_GradonFace = (Current_Cell->Get_Grad_TatCenter() + Neighbour_Cell->Get_Grad_TatCenter()) * 0.5;
	}
	// 	Avg_U_GradonFace.Print();Avg_V_GradonFace.Print();Avg_W_GradonFace.Print();
	u11 = Avg_U_GradonFace.Get_Component(1);
	u12 = Avg_U_GradonFace.Get_Component(2);
	u13 = Avg_U_GradonFace.Get_Component(3);
	u21 = Avg_V_GradonFace.Get_Component(1);
	u22 = Avg_V_GradonFace.Get_Component(2);
	u23 = Avg_V_GradonFace.Get_Component(3);
	u31 = Avg_W_GradonFace.Get_Component(1);
	u32 = Avg_W_GradonFace.Get_Component(2);
	u33 = Avg_W_GradonFace.Get_Component(3);

	// 	cout<<"Viscosity\t"<<mu<<endl;
	temp = (2.0 / 3.0) * (u11 + u22 + u23);

	T11 = mu * (2.0 * u11 - temp);
	T12 = mu * (u12 + u21);
	T13 = mu * (u13 + u31);
	T21 = T12;
	T22 = mu * (2.0 * u22 - temp);
	T23 = mu * (u23 + u32);
	T31 = T13;
	T32 = T23;
	T33 = mu * (2.0 * u33 - temp);
	Flux[0] = 0.0;
	Flux[1] = (T11 * Area_Comp1 + T21 * Area_Comp2 + T31 * Area_Comp3);																							 // Ti1Ai
	Flux[2] = (T12 * Area_Comp1 + T22 * Area_Comp2 + T32 * Area_Comp3);																							 // Ti2Ai
	Flux[3] = (T13 * Area_Comp1 + T23 * Area_Comp2 + T33 * Area_Comp3);																							 // Ti3Ai
	Flux[4] = (((T11 * v1 + T12 * v2 + T13 * v3) * Area_Comp1) + (T21 * v1 + T22 * v2 + T23 * v3) * Area_Comp2 + (T31 * v1 + T32 * v2 + T33 * v3) * Area_Comp3); // TijviAj
																																								 // Addition of Heat flux on the face via k*grad(T)
	Flux[4] += K * (Avg_T_GradonFace * Area_Vector);
	// 	cout<<Flux[0]<<"\t\t"<<Flux[1]<<"\t\t"<<Flux[2]<<"\t\t"<<Flux[3]<<"\t\t"<<Flux[4]<<"\n";
	// 	cout<<"-------------------------------------------------\n";
	return Flux;
}

const Vector &Cell::Get_Grad_UatCenter() const
{
	return u_GradatCenter;
}

const Vector &Cell::Get_Grad_VatCenter() const
{
	return v_GradatCenter;
}

const Vector &Cell::Get_Grad_WatCenter() const
{
	return w_GradatCenter;
}

const Vector &Cell::Get_Grad_TatCenter() const
{
	return T_GradatCenter;
}
