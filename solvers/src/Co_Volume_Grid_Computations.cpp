#include "definitions.h"
#include "Globals.h"
#include "Grid.h"
#include <cmath>

namespace
{
void Append_Co_Face_Geometry(const V_D &p0, const V_D &p1, const V_D &p2, const V_D &p3,
							 V_D &Face_Area_Components, V_D &Face_Normal_Components)
{
	const double edges[4][2] = {
		{p0[0] - p3[0], p0[1] - p3[1]},
		{p1[0] - p0[0], p1[1] - p0[1]},
		{p2[0] - p1[0], p2[1] - p1[1]},
		{p3[0] - p2[0], p3[1] - p2[1]},
	};
	for (int e = 0; e < 4; ++e)
	{
		const double dx = edges[e][0];
		const double dy = edges[e][1];
		double L = std::sqrt(dx * dx + dy * dy);
		if (!(L > 1e-16) || !std::isfinite(L))
			L = 1e-16;
		Face_Area_Components.push_back(L);
		Face_Normal_Components.push_back(dy / L);
		Face_Normal_Components.push_back(-dx / L);
	}
}

double Safe_Cross_Area(const V_D &A, const V_D &B)
{
	const double Area = 0.5 * (B[1] * A[0] - B[0] * A[1]);
	if (!std::isfinite(Area))
		return 1e-16;
	return (std::fabs(Area) < 1e-16) ? 1e-16 : std::fabs(Area);
}
} // namespace

void Construct_Co_Volumes(int &Current_Cell_No)
{
	if (Cells[Current_Cell_No].numFaces != 4)
		return;
	if (Cells[Current_Cell_No].Cell_Vertices.size() < 12)
		return;

	Cell Grid_Cells = {};
	Grid_Cells.cellID = Current_Cell_No;
	Grid_Cells.Dimension = 2;

	V_D o(3, 0.0), a(3, 0.0), b(3, 0.0), c(3, 0.0);
	V_D MP1(3, 0.0), MP2(3, 0.0), MP3(3, 0.0), MP4(3, 0.0), MP5(3, 0.0);
	V_D Diagonal_Vector1(2, 0.0), Diagonal_Vector2(2, 0.0);
	V_D Co_Area(4, 0.0);
	V_D Face_Area_Components;
	V_D Face_Normal_Components;
	Face_Area_Components.reserve(16);
	Face_Normal_Components.reserve(32);

	const int N1 = Cells[Current_Cell_No].Neighbours[0];
	const int N2 = Cells[Current_Cell_No].Neighbours[1];
	const int N3 = Cells[Current_Cell_No].Neighbours[2];
	const int N4 = Cells[Current_Cell_No].Neighbours[3];
	if (N1 < 0 || N2 < 0 || N3 < 0 || N4 < 0)
		return;
	if (N1 >= static_cast<int>(Cells.size()) || N2 >= static_cast<int>(Cells.size()) ||
		N3 >= static_cast<int>(Cells.size()) || N4 >= static_cast<int>(Cells.size()))
		return;

	const auto &verts = Cells[Current_Cell_No].Cell_Vertices;
	/* Quad vertices stored as [v0,v1,v2,v3] × (x,y,z). */
	o[0] = verts[0];
	o[1] = verts[1];
	o[2] = verts[2];
	a[0] = verts[3];
	a[1] = verts[4];
	a[2] = verts[5];
	b[0] = verts[6];
	b[1] = verts[7];
	b[2] = verts[8];
	c[0] = verts[9];
	c[1] = verts[10];
	c[2] = verts[11];

	MP1 = Cells[Current_Cell_No].Cell_Center;
	MP2 = Cells[N1].Cell_Center;
	MP3 = Cells[N2].Cell_Center;
	MP4 = Cells[N3].Cell_Center;
	MP5 = Cells[N4].Cell_Center;

	/* Face-0: left neighbour midpoints + vertices o,c */
	Append_Co_Face_Geometry(MP2, o, MP1, c, Face_Area_Components, Face_Normal_Components);
	Diagonal_Vector1[0] = o[0] - c[0];
	Diagonal_Vector1[1] = o[1] - c[1];
	Diagonal_Vector2[0] = MP2[0] - MP1[0];
	Diagonal_Vector2[1] = MP2[1] - MP1[1];
	Co_Area[0] = Safe_Cross_Area(Diagonal_Vector2, Diagonal_Vector1);

	/* Face-1: bottom neighbour */
	Append_Co_Face_Geometry(o, MP3, a, MP1, Face_Area_Components, Face_Normal_Components);
	Diagonal_Vector1[0] = MP3[0] - MP1[0];
	Diagonal_Vector1[1] = MP3[1] - MP1[1];
	Diagonal_Vector2[0] = a[0] - o[0];
	Diagonal_Vector2[1] = a[1] - o[1];
	Co_Area[1] = Safe_Cross_Area(Diagonal_Vector1, Diagonal_Vector2);

	/* Face-2: right neighbour */
	Append_Co_Face_Geometry(MP1, a, MP4, b, Face_Area_Components, Face_Normal_Components);
	Diagonal_Vector1[0] = b[0] - a[0];
	Diagonal_Vector1[1] = b[1] - a[1];
	Diagonal_Vector2[0] = MP1[0] - MP4[0];
	Diagonal_Vector2[1] = MP1[1] - MP4[1];
	Co_Area[2] = Safe_Cross_Area(Diagonal_Vector1, Diagonal_Vector2);

	/* Face-3: top neighbour */
	Append_Co_Face_Geometry(c, MP1, b, MP5, Face_Area_Components, Face_Normal_Components);
	Diagonal_Vector2[0] = c[0] - b[0];
	Diagonal_Vector2[1] = c[1] - b[1];
	Diagonal_Vector1[0] = MP5[0] - MP1[0];
	Diagonal_Vector1[1] = MP5[1] - MP1[1];
	Co_Area[3] = Safe_Cross_Area(Diagonal_Vector1, Diagonal_Vector2);

	const double area_sum = Co_Area[0] + Co_Area[1] + Co_Area[2] + Co_Area[3];
	Grid_Cells.Area = area_sum;
	Grid_Cells.Inv_Area = (area_sum > 1e-16) ? (1.0 / area_sum) : 0.0;
	Grid_Cells.Face_Areas = Face_Area_Components;
	Grid_Cells.Face_Normals = Face_Normal_Components;
	Grid_Cells.Cell_Areas = Co_Area;
	Grid_Cells.numFaces = 4;

	if (Current_Cell_No == static_cast<int>(Co_Volume_Cells.size()))
		Co_Volume_Cells.push_back(Grid_Cells);
	else if (Current_Cell_No < static_cast<int>(Co_Volume_Cells.size()))
		Co_Volume_Cells[Current_Cell_No] = Grid_Cells;
	else
	{
		Co_Volume_Cells.resize(Current_Cell_No + 1);
		Co_Volume_Cells[Current_Cell_No] = Grid_Cells;
	}
}
