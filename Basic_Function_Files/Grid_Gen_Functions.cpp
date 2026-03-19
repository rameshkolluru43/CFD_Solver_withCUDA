#include "Geometry_Header.h"

/* ======================================================================
   Initialize Grid from four boundary point lists.

   List1 = south (j=0),    List2 = east (i=nx-1),
   List3 = north (j=ny-1), List4 = west (i=0)

   Requires:
     - List1.size() == List3.size()  (matched xi-direction counts)
     - List2.size() == List4.size()  (matched eta-direction counts)
     - Corner connectivity forms a closed loop

   The lists define the boundary curves in counterclockwise order:
     List1: (0,0) -> (nx-1,0)
     List2: (nx-1,0) -> (nx-1,ny-1)
     List3: (nx-1,ny-1) -> (0,ny-1)   [stored in reverse xi order]
     List4: (0,ny-1) -> (0,0)         [stored in reverse eta order]
   ====================================================================== */
void Grid::operator()(vector<Point> &List1, vector<Point> &List2,
					  vector<Point> &List3, vector<Point> &List4)
{
	int No_P_Curve1 = List1.size();
	int No_P_Curve2 = List2.size();
	int No_P_Curve3 = List3.size();
	int No_P_Curve4 = List4.size();

	if (No_P_Curve1 != No_P_Curve3 || No_P_Curve2 != No_P_Curve4)
	{
		cerr << "ERROR: Mismatched point counts on opposing boundaries.\n"
			 << "  South: " << No_P_Curve1 << "  North: " << No_P_Curve3 << "\n"
			 << "  East:  " << No_P_Curve2 << "  West:  " << No_P_Curve4 << endl;
		exit(1);
	}

	cout << "Boundary corner check:\n";
	cout << "  West end  -> South start: ";
	List4[No_P_Curve4 - 1].Print();
	cout << " vs ";
	List1[0].Print();
	cout << "\n  South end -> East start:  ";
	List1[No_P_Curve1 - 1].Print();
	cout << " vs ";
	List2[0].Print();
	cout << "\n  East end  -> North start: ";
	List2[No_P_Curve2 - 1].Print();
	cout << " vs ";
	List3[0].Print();
	cout << "\n  North end -> West start:  ";
	List3[No_P_Curve3 - 1].Print();
	cout << " vs ";
	List4[0].Print();
	cout << endl;

	bool closed = (List4[No_P_Curve4 - 1] == List1[0])
			   && (List1[No_P_Curve1 - 1] == List2[0])
			   && (List2[No_P_Curve2 - 1] == List3[0])
			   && (List3[No_P_Curve3 - 1] == List4[0]);

	if (!closed)
	{
		cerr << "ERROR: Boundary point lists do not form a closed loop.\n"
			 << "Check corner connectivity.\n";
		exit(1);
	}

	cout << "Boundaries form a closed loop. Proceeding with grid generation.\n";

	nx = No_P_Curve1;
	ny = No_P_Curve2;

	x.assign(nx, vector<double>(ny, 0.0));
	y.assign(nx, vector<double>(ny, 0.0));
	xtemp.assign(nx, vector<double>(ny, 0.0));
	ytemp.assign(nx, vector<double>(ny, 0.0));
	erx.assign(nx, vector<double>(ny, 0.0));
	ery.assign(nx, vector<double>(ny, 0.0));
	P_source.assign(nx, vector<double>(ny, 0.0));
	Q_source.assign(nx, vector<double>(ny, 0.0));

	// South boundary (j=0): List1 in forward order
	for (int i = 0; i < nx; i++)
	{
		x[i][0] = List1[i].Get_x();
		y[i][0] = List1[i].Get_y();
	}

	// East boundary (i=nx-1): List2 in forward order
	for (int j = 0; j < ny; j++)
	{
		x[nx - 1][j] = List2[j].Get_x();
		y[nx - 1][j] = List2[j].Get_y();
	}

	// North boundary (j=ny-1): List3 in reverse xi order
	for (int i = 0; i < No_P_Curve3; i++)
	{
		x[nx - 1 - i][ny - 1] = List3[i].Get_x();
		y[nx - 1 - i][ny - 1] = List3[i].Get_y();
	}

	// West boundary (i=0): List4 in reverse eta order
	for (int j = 0; j < No_P_Curve4; j++)
	{
		x[0][ny - 1 - j] = List4[j].Get_x();
		y[0][ny - 1 - j] = List4[j].Get_y();
	}

	cout << "Grid boundaries initialized: " << nx << " x " << ny << " = "
		 << nx * ny << " points\n";
}
