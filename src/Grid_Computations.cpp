#include "definitions.h"
#include "Grid.h"
#include "Globals.h"
#include "Utilities.h"

vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;

V_I Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List, Far_Field_Out_Flow_List, Far_Field_InFlow_List;
V_D Vertices;
int Total_No_Cells, No_Cartesian_Cells, No_Polar_Cells, No_Physical_Cells, No_Ghost_Cells, Cells_in_Plane, nx_c, ny_c, nz_c, nx_p, ny_p, nz_p, Grid_Type;

double global_temp, R_Mid_dot_A, Cell_Minimum_Length;
// int nx_1,nx_2,ny_1,ny_2;

bool Is_Viscous_Wall, Is_2D_Flow, Is_Inlet_SubSonic, Is_Exit_SubSonic, Enable_Far_Field, has_Symmetry_BC;
double CFL;
int numNodes, nodeIndex, PointCellType, LineCellType, TriangleCellType, QuadrilateralCellType, HexahedronCellType, TetrahedronCellType, WedgeCellType;

int get_NoPhysical_Cells()
{
	return Total_No_Cells;
}

/**
 * @brief Reads the input grid file and constructs physical cells, evaluates face areas, normals, and cell volumes.
 *
 * @param ipfile The path to the input grid file.
 *
 * Function for reading input grid file. This file has particular format written by me for reading the grid data.
 * The grid file contains the information of the vertices of the cell followed by the cell and its neighbours information in .txt file.
 * It takes a file containing data of grid and constructs Physical Cells, evaluates its face areas, normals and cell volumes.
 * Creates the storage vectors for the above said variables.
 **/

void Read_Grid(const string &ipfile)
{

	V_D P1(3, 0.0), P2(3, 0.0), P3(3, 0.0), P4(3, 0.0);
	Cell Grid_Cell = {};
	double x = 0.0, y = 0.0, z = 0.0;
	int cellID, neighbor1, neighbor2, neighbor3, neighbor4, Conversion_Type;
	Vertices.resize(12, 0.0);
	cout << "Reading \t" << ipfile.c_str() << endl;
	ifstream Grid_File(ipfile.c_str(), ios::in);
	if (Grid_File.is_open())
	{
		cout << "File opened for reading \n";
		Grid_File >> Grid_Type;
		Grid_File >> Conversion_Type;
		if (Conversion_Type == 1)
			cout << "Grid is created in mm\n Caution convert point data to meters\n";
		else
			cout << "Grid is Created in meters\n";
		switch (Grid_Type)
		{
		case 0: // For Cartesian Mesh
			cout << "Reading Grid of type\t" << Grid_Type << "\t Cartesian Mesh\n";
			Grid_File >> nx_c >> ny_c;		// Represents Points on Cartesian Mesh
			Grid_File >> No_Physical_Cells; // Represents Points on combination of Cartesian cout<<"nx\t"<<nx_c<<"\tny\t"<<ny_c<<"\tnz\t"<<nz_c<<"\t"<<No_Physical_Cells<<endl;
			Cells_in_Plane = (nx_c - 1) * (ny_c - 1);
			No_Ghost_Cells = 2 * (nx_c - 1) + 2 * (ny_c - 1);
			cout << "Total Number of Ghost Cells to be constructed \t" << No_Ghost_Cells << endl;
			Total_No_Cells = No_Physical_Cells + No_Ghost_Cells;
			cout << "Number of Physical Cells\t" << No_Physical_Cells << endl;
			cout << "Number of Ghost Cells\t" << No_Ghost_Cells << endl;
			cout << "Total Number of Cells\t" << Total_No_Cells << endl;
			break;
		case 1: // For Shock Diffraction and Forward Step cases
			cout << "Reading Grid of type\t" << Grid_Type << "\t Dual block mesh\n";
			Grid_File >> nx_1 >> ny_1 >> nx_2 >> ny_2; // Mesh with Two blocks
			Grid_File >> No_Physical_Cells;			   //
			cout << "nx1\t" << nx_1 << "\tny1\t" << ny_1 << "\tnx2\t" << nx_2 << "\tny2\t" << ny_2 << "\t" << No_Physical_Cells << endl;
			Cells_in_Plane = (nx_1 - 1) * (ny_1 - 1) + (nx_2 - 1) * (ny_2 - 1);
			No_Ghost_Cells = (nx_1 - 1) + 2 * (ny_1 - 1) + (nx_2 - 1) + 2 * (ny_2 - 1) + (nx_2 - 1) - (nx_1 - 1);
			cout << "Total Number of Ghost Cells to be constructed \t" << No_Ghost_Cells << endl;
			if (No_Physical_Cells == Cells_in_Plane)
			{
				cout << "Dual Block Grid is accurate" << endl;
			}
			else
			{
				cout << "No of Physical Cells\t" << No_Physical_Cells << endl;
				cout << "No of Cells in plane\t" << Cells_in_Plane << endl;
				cout << "Check the Grid data" << endl;
				exit(0);
			}
			Total_No_Cells = No_Physical_Cells + No_Ghost_Cells;
			break;
		}
		// Reading the Physical Cells Information which are the vertices of the cell
		// followed by the cell and its neighbours information
		for (int i = 0; i < No_Physical_Cells; i++)
		{
			//			cout<<"Physical Cell Number "<<i<<endl;
			x = 0.0;
			y = 0.0;
			z = 0.0;
			Grid_File >> x >> y >> z;
			P1[0] = x;
			P1[1] = y;
			P1[2] = z;
			Grid_File >> x >> y >> z;
			P2[0] = x;
			P2[1] = y;
			P2[2] = z;
			Grid_File >> x >> y >> z;
			P3[0] = x;
			P3[1] = y;
			P3[2] = z;
			Grid_File >> x >> y >> z;
			P4[0] = x;
			P4[1] = y;
			P4[2] = z;
			switch (Conversion_Type)
			{
			case 1:
				Conversion_Factor(P1);
				Conversion_Factor(P2);
				Conversion_Factor(P3);
				Conversion_Factor(P4);
				break;
			default:
				// 					cout<<"Grid is created in m, No conversion requited"<<endl;
				break;
			}
			// --------------------Vertices of Current cell -----------------------------------
			Vertices[0] = P1[0];
			Vertices[1] = P1[1];
			Vertices[2] = P1[2];
			Vertices[3] = P2[0];
			Vertices[4] = P2[1];
			Vertices[5] = P2[2];
			Vertices[6] = P3[0];
			Vertices[7] = P3[1];
			Vertices[8] = P3[2];
			Vertices[9] = P4[0];
			Vertices[10] = P4[1];
			Vertices[11] = P4[2];
			//-----------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------
			Grid_Cell.cellID = i;
			Grid_Cell.Dimension = 2;

			if (Grid_Cell.Cell_Vertices.empty())
				Grid_Cell.Cell_Vertices.resize(12, 0.0);
			Grid_Cell.Cell_Vertices = Vertices;
			if (Grid_Cell.Neighbours.empty())
				Grid_Cell.Neighbours.resize(4, 0);
			// Read the Neighbouring Cell Information
			Grid_File >> cellID >> neighbor1 >> neighbor2 >> neighbor3 >> neighbor4;
			Grid_Cell.Neighbours = {neighbor1, neighbor2, neighbor3, neighbor4};
			Grid_Cell.faceID.resize(4, 0);
			Grid_Cell.nodeIndices.resize(4, 0);
			Grid_Cell.Area = 0.0;
			Grid_Cell.Inv_Area = 0.0;
			// Enforce anti-clockwise vertex ordering and node indices for geometric consistency
			Grid_Cell.nodeIndices = {0, 1, 2, 3};
			Sort_Points_AntiClockWise(Grid_Cell.Cell_Vertices, Grid_Cell.nodeIndices);
			Construct_Cell(Grid_Cell);
			Cells.push_back(Grid_Cell);
			// Delete all the data in the Grid_Cell
			Grid_Cell.Cell_Vertices.clear();
			Grid_Cell.Neighbours.clear();
			Grid_Cell.faceID.clear();
			Grid_Cell.nodeIndices.clear();
			Grid_Cell.Face_Areas.clear();
			Grid_Cell.Face_Normals.clear();
			//-----------------------------------------------------------------------------------
		}
		cout << "Construction of Physical Cells Completed\n";
		Check_Cells();
		cout << "Checking Cells done" << endl;
		cout << "Reading Boundary Cells Information for Boundary Conditions\n";
		string Str;
		int Boundary_Type = 0, No_Cells = 0, a = 0, b = 0, c = 0;
		Grid_File >> Str >> Boundary_Type >> No_Cells;
		cout << Str << "\t" << Boundary_Type << "\t" << No_Cells << endl;
		if (Boundary_Type == 0)
		{
			for (int i = 0; i < No_Cells; i++)
			{
				Grid_File >> a >> b >> c;
				// cout << a << "\t" << b << "\t" << c << endl;
				Inlet_Cells_List.push_back(a);
				Inlet_Cells_List.push_back(b);
				Inlet_Cells_List.push_back(c);
			}
		}
		Grid_File >> Str >> Boundary_Type >> No_Cells;
		cout << Str << "\t" << Boundary_Type << "\t" << No_Cells << endl;
		if (Boundary_Type == 1)
		{
			for (int i = 0; i < No_Cells; i++)
			{
				Grid_File >> a >> b >> c;
				Exit_Cells_List.push_back(a);
				Exit_Cells_List.push_back(b);
				Exit_Cells_List.push_back(c);
			}
		}
		Grid_File >> Str >> Boundary_Type >> No_Cells;
		cout << Str << "\t" << Boundary_Type << "\t" << No_Cells << endl;
		if (Boundary_Type == 2)
		{
			for (int i = 0; i < No_Cells; i++)
			{
				Grid_File >> a >> b >> c;
				// 				cout<<a<<"\t"<<b<<"\t"<<c<<endl;
				Wall_Cells_List.push_back(a);
				Wall_Cells_List.push_back(b);
				Wall_Cells_List.push_back(c);
			}
		}
		Grid_File >> Str >> Boundary_Type >> No_Cells;
		cout << Str << "\t" << Boundary_Type << "\t" << No_Cells << endl;
		if (Boundary_Type == 3)
		{
			for (int i = 0; i < No_Cells; i++)
			{
				Grid_File >> a >> b >> c;
				// 				cout<<a<<"\t"<<b<<"\t"<<c<<endl;
				Symmetry_Cells_List.push_back(a);
				Symmetry_Cells_List.push_back(b);
				Symmetry_Cells_List.push_back(c);
			}
		}

		/*		Grid_File>>Str>>Boundary_Type>>No_Cells;
				cout<<Str<<"\t"<<Boundary_Type<<"\t"<<No_Cells<<endl;
				if(Boundary_Type==4)
				{
					for(int i=0;i<No_Cells;i++)
					{
						Grid_File>>a>>b>>c;
		// 				cout<<a<<"\t"<<b<<"\t"<<c<<endl;
						Far_Field_InFlow_List.push_back(a);
						Far_Field_InFlow_List.push_back(b);
						Far_Field_InFlow_List.push_back(c);
					}
				}

				cout<<Str<<"\t"<<Boundary_Type<<"\t"<<No_Cells<<endl;
				cout<<Str<<"\t"<<No_Cells<<endl;
				if(Boundary_Type==5)
				{
					for(int i=0;i<No_Cells;i++)
					{
						Grid_File>>a>>b>>c;
		// 				cout<<a<<"\t"<<b<<"\t"<<c<<endl;
						Far_Field_Out_Flow_List.push_back(a);
						Far_Field_Out_Flow_List.push_back(b);
						Far_Field_Out_Flow_List.push_back(c);
					}
				}
				*/
		cout << "Reading and assigning Boundary Cells Done\n";
		Grid_File.close();
	}
	else
	{
		cout << "Could not open the file....Please Check the Grid File Name and path of the file\n";
		cout << ipfile.c_str() << endl;
		exit(0);
	}
}

void Conversion_Factor(V_D &Data)
{
	for (unsigned int i = 0; i < Data.size(); i++)
		Data[i] *= 0.001;
}

/*		Function for Constructing Cells from the data read from grid file				*/
void Form_Cells(const string &ipfile)
{
	cout << "Reading \t" << ipfile.c_str() << endl;
	// Auto-detect loader (JSON: mesh.xnodes/mesh.ynodes or mesh.vtk/mesh.txt, VTK, CSV, or legacy TXT)
	Load_Mesh(ipfile);
	cout << "Constructing Physical Cells done" << endl;
	cout << Cells.size() << endl;
	cout << Inlet_Cells_List.size() << "\t" << Wall_Cells_List.size() << "\t" << Exit_Cells_List.size() << "\t" << Symmetry_Cells_List.size() << endl;
	/*for (int i = 0; i < No_Ghost_Cells; i++)
	{
		Cell Temp;
		Cells.push_back(Temp);
	}*/
	cout << Cells.size() << endl;
	Construct_Ghost_Cells();
	cout << "Constructing Ghost Cells done" << endl;
	Calculate_Cell_Center_Distances();
	cout << "Evaluating Cell Center Distances done" << endl;
	cout << No_Physical_Cells << "\t" << No_Ghost_Cells << "\t" << Total_No_Cells << endl;
	Check_Cells();
	cout << "Checking Cells done" << endl;
	cout << "Construction Co Volumes for NS Solver" << endl;
	if (Is_Viscous_Wall)
	{
		for (int i = 0; i < No_Physical_Cells; i++)
		{
			//			cout<<i<<endl;
			Construct_Co_Volumes(i);
			//
		}
		cout << "Construction of Co-Volumes Completed\n";
	}

	//	exit(0);
}

void Construct_Ghost_Cells()
{
	int Cell_Index = 0, Ghost_Cell_Index = 0, Face_No = 0;
	Construct_Cell();
	cout << "Constructing Ghost Cells" << endl;
	cout << "Number of Cells With Inlet Boundary\t" << Inlet_Cells_List.size() / 3 << endl;
	for (unsigned int i = 0; i < Inlet_Cells_List.size(); i += 3)
	{
		Cell_Index = Inlet_Cells_List[i];
		Face_No = Inlet_Cells_List[i + 1];
		Ghost_Cell_Index = Inlet_Cells_List[i + 2];
		// cout << Cell_Index << "\t" << Face_No << "\t" << Ghost_Cell_Index << endl;
		Construct_Cell(Cell_Index, Face_No, Ghost_Cell_Index);
		//		cout<<Cells_Volume.size()<<"\t"<<Cells_Inv_Volume.size()<<endl;
		// 	Print(Cells_Inv_Volume);
	}

	cout << "Number of Cells With Wall Boundary\t" << Wall_Cells_List.size() / 3 << endl;
	for (unsigned int i = 0; i < Wall_Cells_List.size(); i += 3)
	{

		//	 	cout<<Wall_Cells_List[i]<<"\t";
		Cell_Index = Wall_Cells_List[i];
		Face_No = Wall_Cells_List[i + 1];
		Ghost_Cell_Index = Wall_Cells_List[i + 2];
		//		cout<<Cell_Index<<"\t"<<Face_No<<"\t"<<Ghost_Cell_Index<<endl;
		Construct_Cell(Cell_Index, Face_No, Ghost_Cell_Index);
	}
	cout << "Number of Cells With Exit Boundary\t" << Exit_Cells_List.size() / 3 << endl;
	// 	Print(Cells_Inv_Volume);
	for (unsigned int i = 0; i < Exit_Cells_List.size(); i += 3)
	{
		// 		cout<<Exit_Cells_List[i]<<"\t";
		Cell_Index = Exit_Cells_List[i];
		Face_No = Exit_Cells_List[i + 1];
		Ghost_Cell_Index = Exit_Cells_List[i + 2];
		Construct_Cell(Cell_Index, Face_No, Ghost_Cell_Index);
	}
	// 	Print(Cells_Inv_Volume);
	cout << "Number of Cells with Symmetery Boundary\t" << Symmetry_Cells_List.size() / 3 << endl;
	for (unsigned int i = 0; i < Symmetry_Cells_List.size(); i += 3)
	{
		Cell_Index = Symmetry_Cells_List[i];
		Face_No = Symmetry_Cells_List[i + 1];
		Ghost_Cell_Index = Symmetry_Cells_List[i + 2];
		Construct_Cell(Cell_Index, Face_No, Ghost_Cell_Index);
	}
	// 	Print(Cells_Inv_Volume);
}

/*
Indicies in this function represent the indexes of neighbours of current cell. since all are hexahedrals each cell has 6 neighbours.
Indicies are named as follows
a0 - index of current cell,
								a4(i,j+1)
							|---------------|
							|				|
			a1	(i-1,j)		|	a0 (i,j)	|		a3(i+1,j)
							|				|
							|---------------|
								a2(i,j-1)
a1 - index of Left cell,  a2 - index of Bottom Cell,  a3 - index of Right cell,	 a4 - index of Top cell,
order in which neighbour indexes are defined in input file	((i,j),(i-1,j),(i+1,j),(i,j+1),(i,j-1))
*/

void Set_Indicies_of_Neighbour_Cells(const int &a0, const int &a1, const int &a2, const int &a3, const int &a4)
{
	/*
	//																			Cell index		face_no
	Cell_Neighbour_indicies.push_back(a0); //	current cell		(i,j,k)			0
	Cell_Neighbour_indicies.push_back(a1); // 	left 	cell		(i-1,j,k)		1			0
	Cell_Neighbour_indicies.push_back(a2); // 	bottom	cell		(i,j-1,k)		2			1
	Cell_Neighbour_indicies.push_back(a3); //	right 	cell		(i+1,j,k)		3			2
	Cell_Neighbour_indicies.push_back(a4); //	Top cell			(i,j+1,k)		4			3
	*/
	// cout<<a0<<"\t"<<a1<<"\t"<<a2<<"\t"<<a3<<"\t"<<a4<<endl;

	if (Cells[a0].Neighbours.empty())
		Cells[a0].Neighbours.resize(4, 0);
	Cells[a0].Neighbours[0] = a1; // Left Cell
	Cells[a0].Neighbours[1] = a2; // Bottom Cell
	Cells[a0].Neighbours[2] = a3; // Right Cell
	Cells[a0].Neighbours[3] = a4; // Top Cell
}

/*
 This Function takes 8 points and calculates all the features of a cell. Features include Face area, Face normal, volume for a given cell
 */
void Construct_Cell(V_D &o, V_D &a, V_D &b, V_D &c)
{
	V_D Temp(3, 0.0);
	V_D Diagonal_Vector1(2, 0.0);
	V_D Diagonal_Vector2(2, 0.0);
	double Area = 0.0;
	// Evaluating Cell Center
	Temp[0] = 0.25 * (o[0] + a[0] + b[0] + c[0]);
	Temp[1] = 0.25 * (o[1] + a[1] + b[1] + c[1]);
	Temp[2] = 0.25 * (o[2] + a[2] + b[2] + c[2]);
	// Cells_Cell_Center.push_back(Temp);
	Construct_Face(o, a, b, c); // Left Face		1	(i+1,j,k)	1,0,0
	// Flow is assumed to be in x-y plane and hence the other area and normal components in x direction are made zero to make sure that there is no flux in that direction
	/*	cout<<Face_Area_Components.size()<<endl;
		cout<<Face_Normal_Components.size()<<endl;*/
	// Cell_Face_Areas.push_back(Face_Area_Components);
	// Cell_Face_Normals.push_back(Face_Normal_Components);
	/*	Print(Face_Normal_Components);
		cout<<"-----------------------------"<<endl;
		Print(Face_Area_Components);
		cout<<"-----------------------------"<<endl;
	*/
	// Face_Area_Components.clear();
	// Face_Normal_Components.clear();
	// Cell_Neighbour_indicies.clear();
	Diagonal_Vector1[0] = b[0] - o[0];
	Diagonal_Vector1[1] = b[1] - o[1];
	Diagonal_Vector2[0] = c[0] - a[0];
	Diagonal_Vector2[1] = c[1] - a[1];
	// Print(Diagonal_Vector1);
	// cout<<endl;
	// Print(Diagonal_Vector2);
	// cout<<endl;
	Evaluate_Cross_Product(Diagonal_Vector1, Diagonal_Vector2, Area);
	// Cells_Area.push_back(Area);
	// Cells_Inv_Area.push_back((1.0 / Area));
	//	cout<<Area<<"\t"<<1.0/Area<<endl;
	//	cout<<"Construction of Cell done"<<endl;
}
// Constructing the Cell with grid cell structure as the passing variable

void Construct_Cell(Cell &Grid_Cell)
{
	V_D Temp(3, 0.0);
	V_D Diagonal_Vector1(2, 0.0);
	V_D Diagonal_Vector2(2, 0.0);
	double Area = 0.0;

	// Created the cell center and assign it to the cell center vector
	Grid_Cell.Cell_Center.resize(Temp.size(), 0);

	Grid_Cell.Cell_Center = Temp;
	Compute_Centroid(Grid_Cell);
	// Constructing Face Area and Face Normal
	Construct_Face(Grid_Cell);
	// Constructing Cell Area and Inverse Area
	Diagonal_Vector1[0] = Grid_Cell.Cell_Vertices[6] - Grid_Cell.Cell_Vertices[0];
	Diagonal_Vector1[1] = Grid_Cell.Cell_Vertices[7] - Grid_Cell.Cell_Vertices[1];
	Diagonal_Vector2[0] = Grid_Cell.Cell_Vertices[9] - Grid_Cell.Cell_Vertices[3];
	Diagonal_Vector2[1] = Grid_Cell.Cell_Vertices[10] - Grid_Cell.Cell_Vertices[4];
	Evaluate_Cross_Product(Diagonal_Vector1, Diagonal_Vector2, Area);
	Grid_Cell.Area = Area;
	Grid_Cell.Inv_Area = 1.0 / Area;
	// cout<<"Construction of Cell done"<<endl;
}

void Construct_Cell()
{
	R_Mid_dot_A = 0.0;
	double Area = 0.0;
	V_D Temp(3, 0.0);
	// cout<<"In construct cell()\t"<<Cells_Volume.size()<<"\t"<<Cells_Inv_Volume.size()<<endl;

	// cout<<"In construct cell()\t"<<Cell_Face_Areas.size()<<"\t"<<Cell_Face_Normals.size()<<endl;
	for (int i = 0; i < No_Ghost_Cells; i++)
	{
		// Face_Area_Components.resize(4, 0.0);
		// Face_Normal_Components.resize(8, 0.0);
		// Cell_Face_Areas.push_back(Face_Area_Components);
		// Cell_Face_Normals.push_back(Face_Normal_Components);
		// Cells_Area.push_back(Area);
		// Cells_Inv_Area.push_back(Area);
		// Cells_Cell_Center.push_back(Temp);
	}
	// cout << Cell_Face_Areas.size() << "\t" << Cell_Face_Normals.size() << endl;
}

// This Function creates a Ghost cell
void Construct_Cell(const int &Current_Cell_No, const int &Face_No, const int &Ghost_Cell_No)
{
	cout << "In Function Creating Ghost cell \t" << Current_Cell_No << "\t" << Face_No << "\t" << Ghost_Cell_No << endl;
	// check if the memory for ghost cell is created
	if (Cells.size() <= Ghost_Cell_No)
	{
		cout << "Memory for Ghost Cell is not created \n";
		//	exit(0);
	}
	//	cout<<"Current Cell Number\t"<<Current_Cell_No<<"\tFace Number\t"<<Face_No<<endl;
	V_D o(3, 0.0), a(3, 0.0), b(3, 0.0), c(3, 0.0), MP1(3, 0.0), MP2(3, 0.0), MP3(3, 0.0), MP4(3, 0.0), MP5(3, 0.0), Diagonal_Vector1(2, 0.0), Diagonal_Vector2(2, 0.0);
	V_D Vertex_Mid_Point(3, 0.0);
	V_D Co_Area(4, 0.0);

	double Mid_Point_Distance = 0.0;
	// cout << "Size of Cells \t" << Cells.size() << endl;
	// cout << "Number of Verticies for cell number \t" << Current_Cell_No << "\t" << Cells[Current_Cell_No].Cell_Vertices.size() << endl;
	// cout << Cells[Current_Cell_No].Cell_Vertices.size() << endl;

	o[0] = Cells[Current_Cell_No].Cell_Vertices[0];
	o[1] = Cells[Current_Cell_No].Cell_Vertices[1];
	o[2] = Cells[Current_Cell_No].Cell_Vertices[2];
	a[0] = Cells[Current_Cell_No].Cell_Vertices[3];
	a[1] = Cells[Current_Cell_No].Cell_Vertices[4];
	a[2] = Cells[Current_Cell_No].Cell_Vertices[5];
	b[0] = Cells[Current_Cell_No].Cell_Vertices[6];
	b[1] = Cells[Current_Cell_No].Cell_Vertices[7];
	b[2] = Cells[Current_Cell_No].Cell_Vertices[8];
	c[0] = Cells[Current_Cell_No].Cell_Vertices[9];
	c[1] = Cells[Current_Cell_No].Cell_Vertices[10];
	c[2] = Cells[Current_Cell_No].Cell_Vertices[11];
	//	cout<<"Current Cell Number\t"<<Current_Cell_No<<"\tFace Number\t"<<Face_No<<endl;

	MP1 = Cells[Current_Cell_No].Cell_Center;

	switch (Face_No)
	{
	case 0:
		// if(Ghost_Cell_No >= No_Physical_Cells)
		//{
		Vertex_Mid_Point[0] = 0.5 * (o[0] + c[0]);
		Vertex_Mid_Point[1] = 0.5 * (o[1] + c[1]);
		Vertex_Mid_Point[2] = 0.5 * (o[2] + c[2]);
		Distance_Between_Points(MP1, Vertex_Mid_Point, Mid_Point_Distance);
		MP2[0] = MP1[0] - 2.0 * Mid_Point_Distance;
		MP2[1] = MP1[1];
		MP2[2] = MP1[2];
		if (Cells[Ghost_Cell_No].Cell_Center.empty())
		{
			cout << "Ghost Cell Center is empty\n";
			Cells[Ghost_Cell_No].Cell_Center.resize(MP2.size(), 0.0);
		}
		Cells[Ghost_Cell_No].Cell_Center = MP2;
		//}
		break;
	case 1:
		//		cout<<Neighbour_2<<"\t"<<No_Physical_Cells<<endl;
		// if(Ghost_Cell_No >= No_Physical_Cells)
		//{
		Vertex_Mid_Point[0] = 0.5 * (o[0] + a[0]);
		Vertex_Mid_Point[1] = 0.5 * (o[1] + a[1]);
		Vertex_Mid_Point[2] = 0.5 * (o[2] + a[2]);
		Distance_Between_Points(MP1, Vertex_Mid_Point, Mid_Point_Distance);
		MP3[0] = MP1[0];
		MP3[1] = MP1[1] - 2.0 * Mid_Point_Distance;
		MP3[2] = MP1[2];
		//							cout<<Neighbour_2<<"\t"<<No_Physical_Cells<<endl;
		//							Print(MP1);
		//							Print(MP3);
		if (Cells[Ghost_Cell_No].Cell_Center.empty())
		{
			// cout << "Ghost Cell Center is empty\n";
			Cells[Ghost_Cell_No].Cell_Center.resize(MP3.size(), 0.0);
		}
		Cells[Ghost_Cell_No].Cell_Center = MP3;
		//}
		break;
	case 2:
		//	if(Ghost_Cell_No >= No_Physical_Cells)
		//	{
		Vertex_Mid_Point[0] = 0.5 * (b[0] + a[0]);
		Vertex_Mid_Point[1] = 0.5 * (b[1] + a[1]);
		Vertex_Mid_Point[2] = 0.5 * (b[2] + a[2]);
		Distance_Between_Points(MP1, Vertex_Mid_Point, Mid_Point_Distance);
		MP4[0] = MP1[0] + 2.0 * Mid_Point_Distance;
		MP4[1] = MP1[1];
		MP4[2] = MP1[2];
		if (Cells[Ghost_Cell_No].Cell_Center.empty())
		{
			// cout << "Ghost Cell Center is empty\n";
			Cells[Ghost_Cell_No].Cell_Center.resize(MP4.size(), 0.0);
		}
		Cells[Ghost_Cell_No].Cell_Center = MP4;
		break;
	case 3:
		//	if(Ghost_Cell_No >= No_Physical_Cells)
		//	{
		Vertex_Mid_Point[0] = 0.5 * (b[0] + c[0]);
		Vertex_Mid_Point[1] = 0.5 * (b[1] + c[1]);
		Vertex_Mid_Point[2] = 0.5 * (b[2] + c[2]);
		Distance_Between_Points(MP1, Vertex_Mid_Point, Mid_Point_Distance);
		MP5[0] = MP1[0];
		MP5[1] = MP1[1] + 2.0 * Mid_Point_Distance;
		MP5[2] = MP1[2];
		if (Cells[Ghost_Cell_No].Cell_Center.empty())
		{
			// cout << "Ghost Cell Center is empty\n";
			Cells[Ghost_Cell_No].Cell_Center.resize(MP5.size(), 0.0);
		}
		Cells[Ghost_Cell_No].Cell_Center = MP5;
		//	}

		break;
	}
	Cells[Ghost_Cell_No].Area = Cells[Current_Cell_No].Area;
	Cells[Ghost_Cell_No].Inv_Area = Cells[Current_Cell_No].Inv_Area;
}

// This function constructs the cell based on the vertices being passed to the function.
void Construct_Cell(V_D &Vertices)
{
	if (Vertices.size() == 12)
	{
		V_D o(3, 0.0), a(3, 0.0), b(3, 0.0), c(3, 0.0);
		o[0] = Vertices[0];
		o[1] = Vertices[1];
		o[2] = Vertices[2];
		a[0] = Vertices[3];
		a[1] = Vertices[4];
		a[2] = Vertices[5];
		b[0] = Vertices[6];
		b[1] = Vertices[7];
		b[2] = Vertices[8];
		c[0] = Vertices[9];
		c[1] = Vertices[10];
		c[2] = Vertices[11];
		Construct_Cell(o, a, b, c);
	}
	else if (Vertices.size() == 9)
	{
		V_D o(3, 0.0), a(3, 0.0), b(3, 0.0);
		o[0] = Vertices[0];
		o[1] = Vertices[1];
		o[2] = Vertices[2];
		a[0] = Vertices[3];
		a[1] = Vertices[4];
		a[2] = Vertices[5];
		b[0] = Vertices[6];
		b[1] = Vertices[7];
		b[2] = Vertices[8];
		Construct_Cell(o, a, b);
	}
	else
	{
		cout << "Vertices are not 12 in number\n";
		exit(0);
	}
}

// This function constructs a triangular cell based on the vertices being passed to the function
void Construct_Cell(V_D &o, V_D &a, V_D &b)
{
	// Area of the cell is calculated as the area of the triangle formed by the three points
	V_D Temp(3, 0.0);
	double Area = 0.0;
	Temp[0] = 0.25 * (o[0] + a[0] + b[0]);
	Temp[1] = 0.25 * (o[1] + a[1] + b[1]);
	Temp[2] = 0.25 * (o[2] + a[2] + b[2]);
	// Cells_Cell_Center.push_back(Temp);
	// Cell_Face_Areas.push_back(Face_Area_Components);
	// Cell_Face_Normals.push_back(Face_Normal_Components);
	// Face_Area_Components.clear();
	// Face_Normal_Components.clear();
	// Cell_Neighbour_indicies.clear();
	// Cells_Area.push_back(Area);
	// Cells_Inv_Area.push_back((1.0 / Area));
}

//	This function constructs faces and stores face area vectors and face normals for gride cell structure for a 2D cell
void Construct_Face(Cell &Grid_Cell)
{
	//	cout<<"In construct face function"<<endl;

	double L = 0, nx = 0, ny = 0;
	V_D o(3, 0.0), a(3, 0.0), b(3, 0.0), c(3, 0.0);
	o[0] = Grid_Cell.Cell_Vertices[0];
	o[1] = Grid_Cell.Cell_Vertices[1];
	o[2] = Grid_Cell.Cell_Vertices[2];
	a[0] = Grid_Cell.Cell_Vertices[3];
	a[1] = Grid_Cell.Cell_Vertices[4];
	a[2] = Grid_Cell.Cell_Vertices[5];
	b[0] = Grid_Cell.Cell_Vertices[6];
	b[1] = Grid_Cell.Cell_Vertices[7];
	b[2] = Grid_Cell.Cell_Vertices[8];
	c[0] = Grid_Cell.Cell_Vertices[9];
	c[1] = Grid_Cell.Cell_Vertices[10];
	c[2] = Grid_Cell.Cell_Vertices[11];
	// cout<<o[0]<<"\t"<<o[1]<<"\t"<<a[0]<<"\t"<<a[1]<<"\t"<<b[0]<<"\t"<<b[1]<<"\t"<<c[0]<<"\t"<<c[1]<<endl;
	/*Print(o);
	Print(a);
	Print(b);
	Print(c);*/

	Compute_Centroid(Grid_Cell);
	// cout << "Centroid\t" << Grid_Cell.Cell_Center[0] << "\t" << Grid_Cell.Cell_Center[1] << endl;
	Construct_Face(c, o, L, nx, ny); // Face 0  vector (co)
	Grid_Cell.Face_Areas.push_back(L);
	Grid_Cell.Face_Normals.push_back(nx);
	Grid_Cell.Face_Normals.push_back(ny);

	Construct_Face(o, a, L, nx, ny); // Face 1  vector (oa)
	Grid_Cell.Face_Areas.push_back(L);
	Grid_Cell.Face_Normals.push_back(nx);
	Grid_Cell.Face_Normals.push_back(ny);

	Construct_Face(a, b, L, nx, ny); // Face 2  vector (ab)
	Grid_Cell.Face_Areas.push_back(L);
	Grid_Cell.Face_Normals.push_back(nx);
	Grid_Cell.Face_Normals.push_back(ny);

	Construct_Face(b, c, L, nx, ny); // Face 3  vector (bc)
	Grid_Cell.Face_Areas.push_back(L);
	Grid_Cell.Face_Normals.push_back(nx);
	Grid_Cell.Face_Normals.push_back(ny);
	// cout<<Grid_Cell.Face_Area_Components.size()<<endl;
	// cout<<Grid_Cell.Face_Normal_Components.size()<<endl;
}

// This function takes two verticies at at time and constructs the face and normals of the face and length of the face for a 2D cell
void Construct_Face(V_D &a, V_D &b, double &dL, double &nx, double &ny)
{
	//	cout<<"In construct face function"<<endl;
	dL = 0.0;
	nx = 0.0;
	ny = 0.0;
	double dx = 0.0, dy = 0.0;
	dy = b[1] - a[1];
	dx = b[0] - a[0]; // Face  formed by a and b
	// Length of Each Side of the face
	dL = sqrt(dx * dx + dy * dy);
	if (dL == 0.0)
	{
		cout << "Length of the face is zero\n";
		Print(a);
		Print(b);
		exit(0);
	}
	else
	{
		nx = dy / dL;  //----------------- nx = dy/dl
		ny = -dx / dL; // --------------- ny = -dx/dl
	}
}

// This function takes four verticies at at time and constructs the face and normals of the face and length of the face for a 2D cell
void Construct_Face(V_D &o, V_D &a, V_D &b, V_D &c)
{
	//	   cout<<"In construct face function"<<endl;

	double A1_x = 0.0, A1_y = 0.0, A2_x = 0.0, A2_y = 0.0, A3_x = 0.0, A3_y = 0.0, A4_x = 0.0, A4_y = 0.0;
	double A_Mag1 = 0.0, A_Mag2 = 0.0, A_Mag3 = 0.0, A_Mag4 = 0.0;

	A1_y = o[1] - c[1];
	A1_x = o[0] - c[0]; // Face 0  vector (co)
	A2_y = a[1] - o[1];
	A2_x = a[0] - o[0]; // Face 1  vector (oa)
	A3_y = b[1] - a[1];
	A3_x = b[0] - a[0]; // Face 2  vector (ab)
	A4_y = c[1] - b[1];
	A4_x = c[0] - b[0]; // Face 3  vector (bc)

	// Length of Each Side of the face
	A_Mag1 = sqrt(A1_x * A1_x + A1_y * A1_y);
	A_Mag2 = sqrt(A2_x * A2_x + A2_y * A2_y);
	A_Mag3 = sqrt(A3_x * A3_x + A3_y * A3_y);
	A_Mag4 = sqrt(A4_x * A4_x + A4_y * A4_y);

	/*

		// In 2D these Represents the Lenghts of the sides of the control volume
		Face_Area_Components.push_back(A_Mag1);
		Face_Area_Components.push_back(A_Mag2);
		Face_Area_Components.push_back(A_Mag3);
		Face_Area_Components.push_back(A_Mag4);

		Face_Normal_Components.push_back(A1_y / A_Mag1);  //----------------- nx = dy/dl
		Face_Normal_Components.push_back(-A1_x / A_Mag1); // --------------- ny = -dx/dl

		Face_Normal_Components.push_back(A2_y / A_Mag2);  //----------------- nx = dy/dl
		Face_Normal_Components.push_back(-A2_x / A_Mag2); // --------------- ny = -dx/dl

		Face_Normal_Components.push_back(A3_y / A_Mag3);  //----------------- nx = dy/dl
		Face_Normal_Components.push_back(-A3_x / A_Mag3); // --------------- ny = -dx/dl

		Face_Normal_Components.push_back(A4_y / A_Mag4);  //----------------- nx = dy/dl
		Face_Normal_Components.push_back(-A4_x / A_Mag4); // --------------- ny = -dx/dl

		//	   Print(Face_Normal_Components);
		*/
}

// This function evaluates the distance between the current cell and its neighbours
void Calculate_Cell_Center_Distances()
{
	V_D Cell_Center(3, 0.0), N1_Cell_Center(3, 0.0), N2_Cell_Center(3, 0.0), N3_Cell_Center(3, 0.0), N4_Cell_Center(3, 0.0);
	double D1 = 0.0, D2 = 0.0, D3 = 0.0, D4 = 0.0;
	int N1, N2, N3, N4;
	//	cout<<"Calculating Cell Center distances"<<endl;
	for (int Cell_Index = 0; Cell_Index < No_Physical_Cells; Cell_Index++)
	{
		// cout << " am here with cell index\t" << Cell_Index << endl;
		//  Checking if the Cell Center distances is allocated or not
		//  if not allocated then allocate w.r.t the number of neighbours
		if (Cells[Cell_Index].Cell_Center_Distances.empty())
		{
			Cells[Cell_Index].Cell_Center_Distances.resize(Cells[Cell_Index].Neighbours.size(), 0.0);
		}
		if (Cells[Cell_Index].Cell_Center_Vector.empty())
		{
			Cells[Cell_Index].Cell_Center_Vector.resize(3 * Cells[Cell_Index].Neighbours.size(), 0.0);
		}
		// Fetching Cell Center
		Cell_Center = Cells[Cell_Index].Cell_Center;
		//  Fetching Neighbour Indicies
		N1 = Cells[Cell_Index].Neighbours[0];
		N2 = Cells[Cell_Index].Neighbours[1];
		N3 = Cells[Cell_Index].Neighbours[2];
		// Fetching the Neigbhours Cell Centers
		N1_Cell_Center = Cells[N1].Cell_Center;
		N2_Cell_Center = Cells[N2].Cell_Center;
		N3_Cell_Center = Cells[N3].Cell_Center;

		Distance_Between_Points(Cell_Center, N1_Cell_Center, D1);
		Distance_Between_Points(Cell_Center, N2_Cell_Center, D2);
		Distance_Between_Points(Cell_Center, N3_Cell_Center, D3);

		Cells[Cell_Index].Cell_Center_Distances[0] = D1;
		Cells[Cell_Index].Cell_Center_Distances[1] = D2;
		Cells[Cell_Index].Cell_Center_Distances[2] = D3;

		Cells[Cell_Index].Cell_Center_Vector[0] = N1_Cell_Center[0] - Cell_Center[0];
		Cells[Cell_Index].Cell_Center_Vector[1] = N1_Cell_Center[1] - Cell_Center[1];
		Cells[Cell_Index].Cell_Center_Vector[2] = N1_Cell_Center[2] - Cell_Center[2];
		Cells[Cell_Index].Cell_Center_Vector[3] = N2_Cell_Center[0] - Cell_Center[0];
		Cells[Cell_Index].Cell_Center_Vector[4] = N2_Cell_Center[1] - Cell_Center[1];
		Cells[Cell_Index].Cell_Center_Vector[5] = N2_Cell_Center[2] - Cell_Center[2];
		Cells[Cell_Index].Cell_Center_Vector[6] = N3_Cell_Center[0] - Cell_Center[0];
		Cells[Cell_Index].Cell_Center_Vector[7] = N3_Cell_Center[1] - Cell_Center[1];
		Cells[Cell_Index].Cell_Center_Vector[8] = N3_Cell_Center[2] - Cell_Center[2];

		if (Cells[Cell_Index].Neighbours.size() == 4)
		{
			N4 = Cells[Cell_Index].Neighbours[3];
			N4_Cell_Center = Cells[N4].Cell_Center;
			Cells[Cell_Index].Cell_Center_Vector[9] = N4_Cell_Center[0] - Cell_Center[0];
			Cells[Cell_Index].Cell_Center_Vector[10] = N4_Cell_Center[1] - Cell_Center[1];
			Cells[Cell_Index].Cell_Center_Vector[11] = N4_Cell_Center[2] - Cell_Center[2];
			Distance_Between_Points(Cell_Center, N4_Cell_Center, D4);
			Cells[Cell_Index].Cell_Center_Distances[3] = D4;
		}
		else
		{
			cout << "Neighbours are not 4 or 3 in number\n";
			exit(0);
		}
		// cout << "Cell Center Distances\t" << D1 << "\t" << D2 << "\t" << D3 << "\t" << D4 << endl;
	}
	cout << "Calculating Distance between Cell Centers Done" << endl;
}

/*void Create_Ghost_Cells_CellCenters(int & Current_Cell_No)
{
}*/

// Function for Evaluting Cross Product of Vectors, this function uses a Global vector to return the value, each time the global vector is cleared to ensure that the elements are zero
void Evaluate_Cross_Product(V_D &A, V_D &B, double &Area)
{
	Area = 0.0;
	Area = 0.5 * (B[1] * A[0] - B[0] * A[1]);
	if (isnan(Area) or (Area == 0.0))
	{
		cout << "In Evaluating Cross Product..... Check Area\t" << Area << endl;
		exit(0);
	}
	// cout<<"Face Area ="<<Area<<endl;
}

// Function Evaluating Face Normals
void Evaluate_Unit_Vector(V_D &AV)
{
	double temp = 0.0;
	for (unsigned int i = 0; i < AV.size(); i++)
		temp += AV[i] * AV[i];
	for (unsigned int i = 0; i < AV.size(); i++)
		AV[i] /= sqrt(temp);
}

// This function evaluates the dot product of two vectors
void Dot_Product(V_D &V1, V_D &V2, double &dotp)
{
	dotp = (V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]);
}

// This function evaluates distance between two points
void Distance_Between_Points(V_D &P1, V_D &P2, double &Distance)
{
	// This function calculates distance between two points
	double size = 0.0, temp = 0.0;
	Distance = 0.0;
	// Print(P1);
	// Print(P2);
	if (P1.size() == P2.size())
	{
		size = P1.size();
		temp = 0.0;
		for (unsigned int i = 0; i < size; i++)
		{
			temp = P2[i] - P1[i];
			Distance += temp * temp;
		}
		Distance = sqrt(Distance);
		//	cout << "Distance between two points\t" << Distance << endl;
	}
	else
	{
		Print(P1);
		Print(P2);
		cout << "In Distace between two points function, Passing Vectors should have same size" << endl;
		exit(0);
	}
}

// This Function evalutates the dot product of two vectors of equal size
void Evaluate_Dot_Product(V_D &A, V_D &B)
{
	global_temp = 0.0;
	for (unsigned int i = 0; i < A.size(); i++)
		global_temp += A[i] * B[i];
}

// This function finds the cell with smallest diagonal vector
void Find_Minimum_Length()
{
	double Min_Len = 100.0, a = 0.0, b = 0.0, c = 0.0, mag = 0.0;
	Cell_Minimum_Length = 100.0;
	for (int index = 0; index < No_Physical_Cells; index++)
	{
		a = fabs(Cells[index].Diagonal_Vector[0]);
		b = fabs(Cells[index].Diagonal_Vector[1]);
		c = fabs(Cells[index].Diagonal_Vector[2]);

		mag = sqrt(a * a + b * b + c * c);
		if (mag < 1e-5)
		{
			cout << "Diagonal Vector is zero\n";
			Print(Cells[index].Cell_Vertices);
			exit(0);
		}

		if (Cell_Minimum_Length > mag)
			Cell_Minimum_Length = mag;
		// 		cout<<a<<"\t"<<b<<"\t"<<c<<"\t"<<Cell_Minimum_Length<<endl;
	}
}

// This Function evaluates the area integral of a all the cells

void Check_Cells()
{
	double N1 = 0.0, N2 = 0.0;
	cout << "Checking Summation of areas to be zero for a given cell" << endl;
	cout << "no of physical cells\t" << No_Physical_Cells << endl;
	cout << "Cells in structured form\t" << Cells.size() << endl;
	if (No_Physical_Cells == Cells.size())
	{
		cout << "No of physical cells and cells in structured form are same\n";
		for (int Cell_Index = 0; Cell_Index < Cells.size(); Cell_Index++)
		{
			N1 = 0.0, N2 = 0.0;
			for (int Face_No = 0; Face_No < 4; Face_No++)
			{
				N1 += Cells[Cell_Index].Face_Normals[Face_No * 2 + 0] * Cells[Cell_Index].Face_Areas[Face_No];
				N2 += Cells[Cell_Index].Face_Normals[Face_No * 2 + 1] * Cells[Cell_Index].Face_Areas[Face_No];
			}
			if ((N1 > 1e-5) || (N2 > 1e-5))
				cout << "Cell No\t" << Cell_Index << "\t" << N1 << "\t" << N2 << endl;
		}
		cout << "Checking summation of Areas done\n";
	}
	else
	{
		cout << "No of physical cells and cells in structured form are not same\n";
		cout << "Check Read grid function\n";
		exit(0);
	}
}

// Function to compute centroid from vertices
void Compute_Centroid(V_D &Vertices, V_D &Centroid)
{
	if (Centroid.empty()) // Fixes NULL issue
	{
		Centroid.resize(3, 0.0);
	}

	double x = 0.0, y = 0.0, z = 0.0;
	int num_points = Vertices.size() / 3; // Assuming 3D points

	for (size_t i = 0; i < Vertices.size(); i += 3)
	{
		x += Vertices[i];
		y += Vertices[i + 1];
		z += Vertices[i + 2];
	}

	if (num_points > 0) // Prevent division by zero
	{
		Centroid[0] = x / num_points;
		Centroid[1] = y / num_points;
		Centroid[2] = z / num_points;
	}
}

// sorting points in anti clockwise direction based on angle from centroid
//  Function to compute angle of a point w.r.t. the centroid
double angleFromCentroid(const V_D &Centroid, const V_D &Point)
{
	double dx = Point[0] - Centroid[0];
	double dy = Point[1] - Centroid[1];
	return atan2(dy, dx); // Computes angle in radians
}

// Computes the centroid of a quadrilateral
void Compute_Centroid(V_D &o, V_D &a, V_D &b, V_D &c, V_D &Centroid)
{
	Centroid[0] = 0.25 * (o[0] + a[0] + b[0] + c[0]);
	Centroid[1] = 0.25 * (o[1] + a[1] + b[1] + c[1]);
	Centroid[2] = 0.25 * (o[2] + a[2] + b[2] + c[2]);
}

// Computes the centroid of a Cell based on the points
void Compute_Centroid(Cell &Grid_Cell)
{
	// cout<<"In Compute Centroid\n";
	if (Grid_Cell.Cell_Center.empty())
	{
		Grid_Cell.Cell_Center.resize(3, 0.0);
	}
	else
	{
		V_D Point1(3, 0.0), Point2(3, 0.0), Point3(3, 0.0);
		Point1[0] = Grid_Cell.Cell_Vertices[0];
		Point1[1] = Grid_Cell.Cell_Vertices[1];
		Point1[2] = Grid_Cell.Cell_Vertices[2];
		Point2[0] = Grid_Cell.Cell_Vertices[3];
		Point2[1] = Grid_Cell.Cell_Vertices[4];
		Point2[2] = Grid_Cell.Cell_Vertices[5];
		Point3[0] = Grid_Cell.Cell_Vertices[6];
		Point3[1] = Grid_Cell.Cell_Vertices[7];
		Point3[2] = Grid_Cell.Cell_Vertices[8];
		if (Grid_Cell.Cell_Vertices.size() == 12)
		{
			V_D Point4(3, 0.0);
			Point4[0] = Grid_Cell.Cell_Vertices[9];
			Point4[1] = Grid_Cell.Cell_Vertices[10];
			Point4[2] = Grid_Cell.Cell_Vertices[11];
			Compute_Centroid(Point1, Point2, Point3, Point4, Grid_Cell.Cell_Center);
		}
		else if (Grid_Cell.Cell_Vertices.size() == 9)
			Compute_Centroid(Point1, Point2, Point3, Grid_Cell.Cell_Center);
		else
		{
			cout << "Points are not 12 or 9 in number\n";
			exit(0);
		}
	}
}

// Computes the centroid of a triangle
void Compute_Centroid(V_D &o, V_D &a, V_D &b, V_D &Centroid)
{
	Centroid[0] = 0.25 * (o[0] + a[0] + b[0]);
	Centroid[1] = 0.25 * (o[1] + a[1] + b[1]);
	Centroid[2] = 0.25 * (o[2] + a[2] + b[2]);
}

/**
 * @brief Sorts a set of points in anti-clockwise order around their centroid.
 *
 * @param Points A flat vector of coordinates representing 3D points. The vector should have a size
 *               that is a multiple of 3, where each triplet (x, y, z) represents a point.
 *
 * Assumptions:
 * - The input vector `Points` is structured as a flat array of 3D points.
 * - The function computes the centroid of the points and sorts them based on their angle
 *   relative to the centroid in the XY plane.
 * - The Z-coordinate is ignored during sorting.
 */
// To check if the points are sorted in anti-clockwise order
bool Are_Points_Sorted_AntiClockWise(const std::vector<V_D> &groupedPoints, const V_D &Centroid)
{
	for (size_t i = 1; i < groupedPoints.size(); ++i)
	{
		double anglePrev = atan2(groupedPoints[i - 1][1] - Centroid[1], groupedPoints[i - 1][0] - Centroid[0]);
		double angleCurr = atan2(groupedPoints[i][1] - Centroid[1], groupedPoints[i][0] - Centroid[0]);
		if (anglePrev > angleCurr) // If any point is out of order
		{
			return false;
		}
	}
	return true;
}

void Sort_Points_AntiClockWise(V_D &Points)
{
	// Check if the input is valid
	if (Points.empty() || Points.size() % 3 != 0)
	{
		std::cerr << "Error: Points vector must be non-empty and have a size that is a multiple of 3." << std::endl;
		return;
	}

	V_D Centroid(3, 0.0);
	Compute_Centroid(Points, Centroid);

	// Grouping Points into 3D points
	std::vector<V_D> groupedPoints;
	for (size_t i = 0; i < Points.size(); i += 3)
	{
		groupedPoints.push_back({Points[i], Points[i + 1], Points[i + 2]});
	}

	// Check if points are already sorted
	if (!Are_Points_Sorted_AntiClockWise(groupedPoints, Centroid))
	{
		// Sorting based on angles from the centroid in the XY plane
		std::sort(groupedPoints.begin(), groupedPoints.end(),
				  [&Centroid](const V_D &a, const V_D &b)
				  {
					  double angleA = atan2(a[1] - Centroid[1], a[0] - Centroid[0]);
					  double angleB = atan2(b[1] - Centroid[1], b[0] - Centroid[0]);
					  if (std::abs(angleA - angleB) < 1e-9) // Handle floating-point precision issues
					  {
						  // If angles are the same, sort by Z-coordinate
						  return a[2] < b[2];
					  }
					  return angleA < angleB;
				  });
	}
	else
	{
		std::cout << "Points are already sorted in anti-clockwise order." << std::endl;
	}

	// Flattening the sorted points back into the original Points vector
	Points.clear();
	for (const auto &point : groupedPoints)
	{
		Points.insert(Points.end(), point.begin(), point.end());
	}
}

void Sort_Points_AntiClockWise(V_D &Points, V_I &Indices)
{
	// Check if the input is valid
	if (Points.empty() || Points.size() % 3 != 0 || Indices.size() * 3 != Points.size())
	{
		std::cerr << "Error: Mismatch between Points and Indices." << std::endl;
		return;
	}

	V_D Centroid(3, 0.0);
	Compute_Centroid(Points, Centroid);

	// Pair each point with its index
	std::vector<std::pair<V_D, int>> pointIndexPairs;
	for (size_t i = 0; i < Points.size(); i += 3)
	{
		V_D point = {Points[i], Points[i + 1], Points[i + 2]};
		pointIndexPairs.emplace_back(point, Indices[i / 3]);
	}

	// Check if points are already sorted
	std::vector<V_D> rawPoints;
	for (const auto &pair : pointIndexPairs)
		rawPoints.push_back(pair.first);

	if (!Are_Points_Sorted_AntiClockWise(rawPoints, Centroid))
	{
		std::sort(pointIndexPairs.begin(), pointIndexPairs.end(),
				  [&Centroid](const std::pair<V_D, int> &a, const std::pair<V_D, int> &b)
				  {
					  double angleA = atan2(a.first[1] - Centroid[1], a.first[0] - Centroid[0]);
					  double angleB = atan2(b.first[1] - Centroid[1], b.first[0] - Centroid[0]);
					  if (std::abs(angleA - angleB) < 1e-9)
						  return a.first[2] < b.first[2];
					  return angleA < angleB;
				  });
	}
	else
	{
		// std::cout << "Points are already sorted in anti-clockwise order." << std::endl;
	}

	// Rebuild the Points and Indices vector
	Points.clear();
	Indices.clear();
	for (const auto &pair : pointIndexPairs)
	{
		Points.insert(Points.end(), pair.first.begin(), pair.first.end());
		Indices.push_back(pair.second);
	}
}

bool Is_On_Boundary(double &x, double &y, double &minX, double &maxX, double &minY, double &maxY)
{
	return std::abs(x - minX) < EPSILON || std::abs(x - maxX) < EPSILON ||
		   std::abs(y - minY) < EPSILON || std::abs(y - maxY) < EPSILON;
}