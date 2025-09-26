#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"

// Function to compute the flux Jacobian matrix Ac
vector<V_D> Compute_Flux_Jacobian(int &Cell_No, vector<V_D> &Ac, int &Face_No)
{
	// Precomputed values
	//            cout<<"in Flux Jacobian function\t with Cell number "<<Cell_No<<"\tand face number \t"<<Face_No<<endl;
	Vector_Reset(Ac);

	//	    Print(Ac);
	double u = Primitive_Cells[Cell_No][1];
	double v = Primitive_Cells[Cell_No][2];
	double E = Primitive_Cells[Cell_No][6];
	double Vt = 0.0;

	//	    cout<<u<<"\t"<<v<<"\t"<<E<<"\t"<<gamma<<endl;

	int index = Face_No * 2;
	//      Obtaining the Normals of the face and face length
	/*nx = Cell_Face_Normals[Cell_No][index+0];//----------------- nx = dy/dl
	ny = Cell_Face_Normals[Cell_No][index+1];//----------------- ny = -dx/dl*/
	nx = Cells[Cell_No].Face_Normals[index + 0];
	ny = Cells[Cell_No].Face_Normals[index + 1];

	//	    cout<<"Normals \t"<<nx<<"\t"<<ny<<endl;

	double phi = 0.5 * (gamma - 1) * (u * u + v * v);
	double V = nx * u + ny * v;
	double a1 = gamma * E - phi;
	double a2 = gamma - 1;
	double a3 = gamma - 2;

	// Fill the matrix Ac
	Ac[0][0] = -Vt;
	Ac[0][1] = nx;
	Ac[0][2] = ny;
	Ac[0][3] = 0;

	Ac[1][0] = nx * phi - u * V;
	Ac[1][1] = V - Vt - a3 * nx * u;
	Ac[1][2] = ny * u - a2 * nx * v;
	Ac[1][3] = a2 * nx;

	Ac[2][0] = ny * phi - v * V;
	Ac[2][1] = nx * v - a2 * ny * u;
	Ac[2][2] = V - Vt - a3 * ny * v;
	Ac[2][3] = a2 * ny;

	Ac[3][0] = V * (phi - a1);
	Ac[3][1] = nx * a1 - a2 * u * V;
	Ac[3][2] = ny * a1 - a2 * v * V;
	Ac[3][3] = gamma * V - Vt;

	//	    cout<<"THe value of of flux Jacobian matrix is"<<endl;

	//	    Print(Ac);
	return Ac;
}

// Function to compute the flux Jacobian matrix Ac
vector<V_D> ComputeGhostCell_Flux_Jacobian(int &Ghost_Cell_No, int &Cell_No, vector<V_D> &Ac, int &Face_No)
{
	// Precomputed values
	//            cout<<"in Flux Jacobian function\t with Cell number "<<Cell_No<<"\tand face number \t"<<Face_No<<endl;
	Vector_Reset(Ac);

	//	    Print(Ac);
	double u = Primitive_Cells[Ghost_Cell_No][1];
	double v = Primitive_Cells[Ghost_Cell_No][2];
	double E = Primitive_Cells[Ghost_Cell_No][6];
	double Vt = 0.0;

	//	    cout<<u<<"\t"<<v<<"\t"<<E<<"\t"<<gamma<<endl;

	int index = Face_No * 2;
	//      Obtaining the Normals of the face and face length
	/* nx = Cell_Face_Normals[Cell_No][index+0];//----------------- nx = dy/dl
	 ny = Cell_Face_Normals[Cell_No][index+1];//----------------- ny = -dx/dl*/
	nx = Cells[Cell_No].Face_Normals[index + 0];
	ny = Cells[Cell_No].Face_Normals[index + 1];

	//	    cout<<"Normals \t"<<nx<<"\t"<<ny<<endl;

	double phi = 0.5 * (gamma - 1) * (u * u + v * v);
	double V = nx * u + ny * v;
	double a1 = gamma * E - phi;
	double a2 = gamma - 1;
	double a3 = gamma - 2;

	// Fill the matrix Ac
	Ac[0][0] = -Vt;
	Ac[0][1] = nx;
	Ac[0][2] = ny;
	Ac[0][3] = 0;

	Ac[1][0] = nx * phi - u * V;
	Ac[1][1] = V - Vt - a3 * nx * u;
	Ac[1][2] = ny * u - a2 * nx * v;
	Ac[1][3] = a2 * nx;

	Ac[2][0] = ny * phi - v * V;
	Ac[2][1] = nx * v - a2 * ny * u;
	Ac[2][2] = V - Vt - a3 * ny * v;
	Ac[2][3] = a2 * ny;

	Ac[3][0] = V * (phi - a1);
	Ac[3][1] = nx * a1 - a2 * u * V;
	Ac[3][2] = ny * a1 - a2 * v * V;
	Ac[3][3] = gamma * V - Vt;

	//	    cout<<"THe value of of flux Jacobian matrix is"<<endl;

	//	    Print(Ac);
	return Ac;
}

void CheckMatrixForErrors(vector<V_D> &A)
{
	bool found_nan = false;

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			double value = A[i][j];
			if (std::isnan(value))
			{
				cout << "NaN found at position (" << i << ", " << j << ")" << endl;
				found_nan = true;
				exit(0);
			}
			else if (value == INFINITY || value == -INFINITY)
			{
				cout << "Infinity found at position (" << i << ", " << j << ")" << endl;
				found_nan = true;
				exit(0);
			}
		}
	}
}

/*// Function to compute the Jacobian matrix for the Euler equations in 2D
vector<V_D> compute_Jacobian(double u, double v, double E, double nx, double ny, double gamma) {
	vector<V_D> Ac(4, vector<double>(4, 0.0));

	double phi = 0.5 * (gamma - 1) * (u * u + v * v);
	double V = nx * u + ny * v;
	double a1 = gamma * E - phi;
	double a2 = gamma - 1;
	double a3 = gamma - 2;

	// Fill the Jacobian matrix
	Ac[0][1] = nx;
	Ac[0][2] = ny;
	Ac[1][0] = nx * phi - u * V;
	Ac[1][1] = V - a3 * nx * u;
	Ac[1][2] = ny * u - a2 * nx * v;
	Ac[2][0] = ny * phi - v * V;
	Ac[2][1] = nx * v - a2 * ny * u;
	Ac[2][2] = V - a3 * ny * v;
	Ac[3][0] = V * (phi - a1);
	Ac[3][1] = nx * a1 - a2 * u * V;
	Ac[3][2] = ny * a1 - a2 * v * V;

	return Ac;
}
*/