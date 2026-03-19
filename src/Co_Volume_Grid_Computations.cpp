#include "definitions.h"
#include "Globals.h"
#include "Grid.h"

void Construct_Co_Volumes(int &Current_Cell_No)
{
	if (Cells[Current_Cell_No].numFaces != 4)
		return;
	Cell Grid_Cells = {};
	Grid_Cells.cellID = Current_Cell_No;
	Grid_Cells.Dimension = 2;
	Grid_Cells.Cell_Vertices.resize(4, 0.0);
	V_D o(3, 0.0), a(3, 0.0), b(3, 0.0), c(3, 0.0), MP1(3, 0.0), MP2(3, 0.0), MP3(3, 0.0), MP4(3, 0.0), MP5(3, 0.0), Diagonal_Vector1(2, 0.0), Diagonal_Vector2(2, 0.0);
	V_D Vertex_Mid_Point(3, 0.0);
	V_D Co_Area(4, 0.0);
	V_D Face_Area_Components(16, 0.0);
	V_D Face_Normal_Components(32, 0.0);
	Neighbour_1 = 0;
	Neighbour_2 = 0;
	Neighbour_3 = 0;
	Neighbour_4 = 0;
	double Area = 0.0, Mid_Point_Distance = 0.0;

	Neighbour_1 = Cells[Current_Cell_No].Neighbours[0];
	Neighbour_2 = Cells[Current_Cell_No].Neighbours[1];
	Neighbour_3 = Cells[Current_Cell_No].Neighbours[2];
	Neighbour_4 = Cells[Current_Cell_No].Neighbours[3];

	//	cout<<Cells_Vertices.size()<<endl;
	//	cout<<Cells_Vertices[Current_Cell_No].size()<<endl;
	o = Cells[Current_Cell_No].Cell_Vertices;
	a = Cells[Current_Cell_No].Cell_Vertices;
	b = Cells[Current_Cell_No].Cell_Vertices;
	c = Cells[Current_Cell_No].Cell_Vertices;

	/*	Print(o);
		Print(a);
		Print(b);
		Print(c);
	*/

	MP1 = Cells[Current_Cell_No].Cell_Center;
	MP2 = Cells[Neighbour_1].Cell_Center;
	MP3 = Cells[Neighbour_2].Cell_Center;
	MP4 = Cells[Neighbour_3].Cell_Center;
	MP5 = Cells[Neighbour_4].Cell_Center;

	/*	Print(MP1);
	Print(MP2);
	Print(MP3);
	Print(MP4);
	Print(MP5);*/

	// exit(0);
	/* ------------- Co Volumes for faces  Face zero is formed by Face-0 ---- Cell ceter of i-1,j , vertex o, cell_center of i,j and vertex c   */
	Construct_Face(MP2, o, MP1, c);

	Diagonal_Vector1[0] = o[0] - c[0];
	Diagonal_Vector1[1] = o[1] - c[1];

	Diagonal_Vector2[0] = MP2[0] - MP1[0];
	Diagonal_Vector2[1] = MP2[1] - MP1[1];
	Evaluate_Cross_Product(Diagonal_Vector2, Diagonal_Vector1, Area);
	Co_Area[0] = Area;
	//	cout<<"On Face 0 done"<<endl;
	/**Face-1 ---- vertex o, Cell ceter of i,j-1 ,  vertex a, cell_center of i,j   */
	Construct_Face(o, MP3, a, MP1);
	Diagonal_Vector1[0] = MP3[0] - MP1[0];
	Diagonal_Vector1[1] = MP3[1] - MP1[1];

	Diagonal_Vector2[0] = a[0] - o[0];
	Diagonal_Vector2[1] = a[1] - o[1];
	//		cout<<"On Face 1"<<endl;
	Evaluate_Cross_Product(Diagonal_Vector1, Diagonal_Vector2, Area);
	Co_Area[1] = Area;
	//	cout<<"On Face 1 done"<<endl;
	/*Face-2 ---- Cell ceter of i,j , vertex a, cell_center of i+1,j and vertex b */
	Construct_Face(MP1, a, MP4, b);
	Diagonal_Vector1[0] = b[0] - a[0];
	Diagonal_Vector1[1] = b[1] - a[1];

	Diagonal_Vector2[0] = MP1[0] - MP4[0];
	Diagonal_Vector2[1] = MP1[1] - MP4[1];
	Evaluate_Cross_Product(Diagonal_Vector1, Diagonal_Vector2, Area);
	Co_Area[2] = Area;
	/*Face-3 ---- vertex b, Cell ceter of i,j+1 , vertex c, cell_center of i,j   */
	Construct_Face(c, MP1, b, MP5);
	Diagonal_Vector2[0] = c[0] - b[0];
	Diagonal_Vector2[1] = c[1] - b[1];

	Diagonal_Vector1[0] = MP5[0] - MP1[0];
	Diagonal_Vector1[1] = MP5[1] - MP1[1];
	Evaluate_Cross_Product(Diagonal_Vector1, Diagonal_Vector2, Area);
	Co_Area[3] = Area;
	Co_Volume_Cells.push_back(Grid_Cells);
	// Initialize and populate Face_Area_Components

	if (Co_Volume_Cells[Current_Cell_No].Face_Areas.empty())
	{
		Co_Volume_Cells[Current_Cell_No].Face_Areas.resize(16, 0.0);
	}
	else
	{
		Co_Volume_Cells[Current_Cell_No].Face_Areas = Face_Area_Components;
	}
	if (Co_Volume_Cells[Current_Cell_No].Face_Normals.empty())
	{
		Co_Volume_Cells[Current_Cell_No].Face_Normals.resize(32, 0.0);
	}
	else
	{
		Co_Volume_Cells[Current_Cell_No].Face_Normals = Face_Normal_Components;
	}
	if (Co_Volume_Cells[Current_Cell_No].Cell_Areas.empty())
	{
		Co_Volume_Cells[Current_Cell_No].Cell_Areas.resize(4, 0.0);
	}
	else
	{
		Co_Volume_Cells[Current_Cell_No].Cell_Areas = Co_Area;
	}
	Face_Area_Components.clear();
	Face_Normal_Components.clear();
	// cout<<"Construction of Co Volume done"<<endl;
}
