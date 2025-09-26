#include "Geometry_Header.h"

Grid::Grid()
{
	ermax = 1e-12;
}

void Grid::operator()(Circle &C)
{
	int i, j, No_Points = C.get_Nop(), Total_Points;
	vector<Point> Point_List;
	Point_List = C.Get_Point_list();
	cout << "Number of Points on Circle\t" << Point_List.size() << endl;
	nx = ny = (No_Points / 4);
	cout << nx << "\t" << ny << endl;
	Total_Points = (nx + 1) * (ny + 1);
	// 	cout<<"=====================\n";
	x = vector<vector<double>>(nx + 1, vector<double>(ny + 1, 0.0));
	y = vector<vector<double>>(nx + 1, vector<double>(ny + 1, 0.0));
	xtemp = vector<vector<double>>(nx + 1, vector<double>(ny + 1, 0.0));
	ytemp = vector<vector<double>>(nx + 1, vector<double>(ny + 1, 0.0));
	erx = vector<vector<double>>(nx + 1, vector<double>(ny + 1, 0.0));
	ery = vector<vector<double>>(nx + 1, vector<double>(ny + 1, 0.0));

	// 	x = vector<double>(Total_Points,0.0);
	// 	y = vector<double>(Total_Points,0.0);
	// 	xtemp = vector<double>(Total_Points,0.0);
	// 	ytemp = vector<double>(Total_Points,0.0);

	for (int i = 0; i < nx; i++)
	{
		j = 0;
		// 		cout<<i<<"\t"<<j<<"\t"<<i+j<<endl;
		x[i][j] = Point_List[i].Get_x();
		y[i][j] = Point_List[i].Get_y();
		// 		cout<<x[i][j]<<"\t"<<y[i][j]<<endl;
	}
	// cout<<"=====================\n";
	for (j = 0; j < (ny + 1); j++)
	{
		i = nx;
		// 		cout<<i<<"\t"<<j<<"\t"<<i+j<<endl;
		x[i][j] = Point_List[i + j].Get_x();
		y[i][j] = Point_List[i + j].Get_y();
		// 		cout<<x[i][j]<<"\t"<<y[i][j]<<endl;
	}
	// 		cout<<"=====================\n";
	for (i = 0; i < (nx + 1); i++)
	{
		j = ny;
		// 		cout<<i<<"\t"<<j<<"\t"<<No_Points -( i+j)<<endl;
		x[i][j] = Point_List[No_Points - (i + j)].Get_x();
		y[i][j] = Point_List[No_Points - (i + j)].Get_y();
		// 		cout<<x[i][j]<<"\t"<<y[i][j]<<endl;
	}
	// 		cout<<"=====================\n";
	for (j = ny; j > 0; j--)
	{
		i = 0;
		// 		cout<<i<<"\t"<<j<<"\t"<<No_Points -( i+j)<<endl;
		x[i][j] = Point_List[No_Points - (i + j)].Get_x();
		y[i][j] = Point_List[No_Points - (i + j)].Get_y();
		// 		cout<<x[i][j]<<"\t"<<y[i][j]<<endl;
	}
	// 		cout<<"=====================\n";
}

void Grid::Read_data(string ipfile)
{
}

/*This Function generates grid using Transfinite interpolation for 2D */
void Grid::TFI()
{
	cout << "Generating Grid using Algebriac Transfinite Interpolation\n";
	cout << "Displaying the control points of the curves \n\n";
	cout << nx << "\t" << ny << endl;
	cout << "Control point 1 \t" << x[0][0] << "\t" << y[0][0] << "\n";
	cout << "Control point 2  \t" << x[nx - 1][0] << "\t" << y[nx - 1][0] << "\n";
	cout << "Control point 3 \t" << x[nx - 1][ny - 1] << "\t" << y[nx - 1][ny - 1] << "\n";
	cout << "Control point 4 \t" << x[0][ny - 1] << "\t" << y[0][ny - 1] << "\n";
	eta = 1.0 / (ny - 1);
	zi = 1.0 / (nx - 1);
	cout << eta << "\t" << zi << endl;

	for (int j = 1; j < ny - 1; j++)
	{
		t = j * eta;
		for (int i = 1; i < nx - 1; i++)
		{
			s = zi * i;
			x[i][j] = t * x[i][ny - 1] + (1.0 - t) * x[i][0] + s * x[nx - 1][j] + (1.0 - s) * x[0][j] - (1 - s) * (t)*x[0][ny - 1] - (1 - s) * (1 - t) * x[0][0] - (s) * (t)*x[nx - 1][ny - 1] - (s) * (1 - t) * x[nx - 1][0];

			y[i][j] = t * y[i][ny - 1] + (1.0 - t) * y[i][0] + s * y[nx - 1][j] + (1.0 - s) * y[0][j] - (1 - s) * (t)*y[0][ny - 1] - (1 - s) * (1 - t) * y[0][0] - (s) * (t)*y[nx - 1][ny - 1] - (s) * (1.0 - t) * y[nx - 1][0];
		}
	}
	xtemp = x;
	ytemp = y;
	cout << "Finished TFI, Matrix Ready to perform Elliptic Grid Generation\n";
}

void Grid::Elliptic()
{
	//	int timer = clock();
	// cout<<inv_eta<<"\t"<<inv_zi<<endl;
	for (int j = 1; j < ny - 1; j++) // loop increment in eta direction
	{
		for (int i = 1; i < nx - 1; i++) // loop  increment in zi direction
		{
			inv_eta = 1.0 / (eta * j);
			inv_zi = 1.0 / (zi * i);

			xzi = 0.5 * inv_zi * (x[i + 1][j] - xtemp[i - 1][j]);
			xeta = 0.5 * inv_eta * (x[i][j + 1] - xtemp[i][j - 1]);
			yzi = 0.5 * inv_zi * (y[i + 1][j] - ytemp[i - 1][j]);
			yeta = 0.5 * inv_eta * (y[i][j + 1] - ytemp[i][j - 1]);

			alpha = (xeta * xeta + yeta * yeta); //(g22)i,j alpha
			gamma1 = (xzi * xzi + yzi * yzi);	 //(g11)i,j gamma
			beta = (xzi * xeta + yzi * yeta);	 //(g12)i,j beta

			B = alpha * inv_zi * inv_zi;
			D = gamma1 * inv_eta * inv_eta;
			C = -0.5 * beta * inv_eta * inv_zi;
			A = B + D;

			xtemp[i][j] = (0.5 / A) * (B * (x[i + 1][j] + xtemp[i - 1][j]) + C * (x[i + 1][j + 1] - xtemp[i + 1][j - 1] - x[i - 1][j + 1] + xtemp[i - 1][j - 1]) + D * (x[i][j + 1] + xtemp[i][j - 1]));
			ytemp[i][j] = (0.5 / A) * (B * (y[i + 1][j] + ytemp[i - 1][j]) + C * (y[i + 1][j + 1] - ytemp[i + 1][j - 1] - y[i - 1][j + 1] + ytemp[i - 1][j - 1]) + D * (y[i][j + 1] + ytemp[i][j - 1]));

		} // end of zi direction loop
	} // end of eta direction loop
	//				timer = clock();
	//				cout<<timer<<endl;
}

void Grid::Periodic_Boundary_Condition()
{
	// 	cout<<"Applying periodic boundary condition\n";
	int i = nx - 1;
	for (int j = 1; j < ny - 1; j++)
	{
		inv_eta = 1.0 / (eta * j);
		inv_zi = 1.0 / (zi * i);

		xzi = 0.5 * inv_zi * (x[1][j] - x[i - 1][j]);
		xeta = 0.5 * inv_eta * (x[i][j + 1] - x[i][j - 1]);
		yzi = 0.5 * inv_zi * (y[1][j] - y[i - 1][j]);
		yeta = 0.5 * inv_eta * (y[i][j + 1] - y[i][j - 1]);
		alpha = (xeta * xeta + yeta * yeta); //(g22)i,j alpha
		gamma1 = (xzi * xzi + yzi * yzi);	 //(g11)i,j gamma
		beta = (xzi * xeta + yzi * yeta);	 //(g12)i,j beta

		B = alpha * inv_zi * inv_zi;
		D = gamma1 * inv_eta * inv_eta;
		C = -0.5 * beta * inv_eta * inv_zi;
		A = B + D;

		xtemp[i][j] = (0.5 / A) * (B * (x[1][j] + x[i - 1][j]) + C * (x[1][j + 1] - x[1][j - 1] - x[i - 1][j + 1] + x[i - 1][j - 1]) + D * (x[i][j + 1] + x[i][j - 1]));

		ytemp[i][j] = (0.5 / A) * (B * (y[1][j] + y[i - 1][j]) + C * (y[1][j + 1] - y[1][j + 1] - y[i - 1][j + 1] + y[i - 1][j - 1]) + D * (y[i][j + 1] + y[i][j - 1]));

		x[0][j] = xtemp[i][j];
		y[0][j] = ytemp[i][j];
	}
}

void Grid::Estimate_Error_and_Update()
{
	ertot = 0.0; // cout<<"ertot="<<ertot<<"\n";
	for (int j = 1; j < ny - 1; j++)
	{
		for (int i = 1; i < nx - 1; i++)
		{
			erx[i][j] = ((xtemp[i][j] - x[i][j]) / x[i][j]) * ((xtemp[i][j] - x[i][j]) / x[i][j]);
			ery[i][j] = ((ytemp[i][j] - y[i][j]) / y[i][j]) * ((ytemp[i][j] - y[i][j]) / y[i][j]);
			ertot += erx[i][j] + ery[i][j];

			x[i][j] = xtemp[i][j];
			y[i][j] = ytemp[i][j];
			xtemp[i][j] = 0.0;
			ytemp[i][j] = 0.0;
			erx[i][j] = 0.0;
			ery[i][j] = 0.0;
		}
	}
	ertot = sqrt(ertot);
}

void Grid::Generate_Grid(bool &With_Periodic_Domain)
{
	noit = 0;

	TFI();
	do
	{
		int timer = clock();
		Elliptic();
		if (With_Periodic_Domain)
			Periodic_Boundary_Condition();
		Estimate_Error_and_Update();
		noit = noit + 1;
		if (noit % 1000 == 0)
		{
			timer = clock();
			cout << noit << "\ttotal error=\t" << ertot << "\t" << timer / CLOCKS_PER_SEC << endl;
		}
	} while ((noit < 100000)); // end of iterative loops

	cout << "Elliptic Grid Generated\t......Generating Grid List\n";
	Generate_Grid_List();
}

void Grid::Generate_Grid(bool &With_Periodic_Domain, bool &EnableElliptic, int &iterations)
{
	noit = 0;
	TFI();
	cout << "Enalbe Elliptic solver\t" << EnableElliptic << "\t with number of iterations\t" << iterations << endl;
	if (EnableElliptic)
	{
		do
		{
			int timer = clock();
			Elliptic();
			if (With_Periodic_Domain)
				Periodic_Boundary_Condition();
			Estimate_Error_and_Update();
			noit = noit + 1;
			if (noit % 1000 == 0)
			{
				timer = clock();
				cout << noit << "\ttotal error=\t" << ertot << "\t" << timer / CLOCKS_PER_SEC << endl;
			}
		} while ((noit < iterations)); // end of iterative loops

		cout << "Elliptic Grid Generated\t......Generating Grid List\n";
	}
	Generate_Grid_List();
}

void Grid::Generate_Grid_List()
{
	Point Temp_Point;
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			Temp_Point(x[i][j], y[i][j], 0.0);
			Grid_List.push_back(Temp_Point);
		}
	}
	cout << "Planar Grid, List generated.....\n Leaving Grid Generator to Basic Program File" << endl;
}

vector<Point> &Grid::Get_Grid_List()
{
	return Grid_List;
}

void Grid::Stack_Grid(double &length, int &nop, vector<Point> &Stacked_Grid_List)
{
	cout << "Stacking grid \n";
	cout << "Grid plane size\t" << Grid_List.size() << endl;
	cout << "Length over which Grid to be stacked\t" << length << endl;
	cout << "Nop in stacking direction\t" << nop << endl;
	vector<Point> temp_list;
	double delx = length / (nop - 1);
	cout << "spacing in Z direction \t" << delx << endl;
	Point p;
	int k = 0;
	for (int i = 0; i < nop; i++)
	{
		for (unsigned int j = 0; j < Grid_List.size(); j++)
		{
			p(Grid_List[j].Get_x(), Grid_List[j].Get_y(), i * delx + Grid_List[j].Get_z());
			Stacked_Grid_List.push_back(p);
			// 			cout<<p.Get_x()<<"\t"<<p.Get_y()<<"\t"<<p.Get_z()<<endl;
		}
		k++;
		//		cout<<k<<"\t"<<"i\t"<<i<<endl;
	}
	cout << "Stacking Elliptic Core Grid Done\t Size of Grid list\t" << Stacked_Grid_List.size() << endl;
}
