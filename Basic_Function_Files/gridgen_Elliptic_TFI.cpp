#include "Geometry_Header.h"
#include <algorithm>
#include <numeric>

Grid::Grid()
{
	ermax = 1e-12;
	ertot = 0.0;
	noit = 0;
}

void Grid::operator()(Circle &C)
{
	int i, j, No_Points = C.get_Nop();
	vector<Point> Point_List;
	Point_List = C.Get_Point_list();
	cout << "Number of Points on Circle\t" << Point_List.size() << endl;
	nx = ny = (No_Points / 4);
	cout << nx << "\t" << ny << endl;

	x.assign(nx + 1, vector<double>(ny + 1, 0.0));
	y.assign(nx + 1, vector<double>(ny + 1, 0.0));
	xtemp.assign(nx + 1, vector<double>(ny + 1, 0.0));
	ytemp.assign(nx + 1, vector<double>(ny + 1, 0.0));
	erx.assign(nx + 1, vector<double>(ny + 1, 0.0));
	ery.assign(nx + 1, vector<double>(ny + 1, 0.0));
	P_source.assign(nx + 1, vector<double>(ny + 1, 0.0));
	Q_source.assign(nx + 1, vector<double>(ny + 1, 0.0));

	for (int i = 0; i < nx; i++)
	{
		j = 0;
		x[i][j] = Point_List[i].Get_x();
		y[i][j] = Point_List[i].Get_y();
	}
	for (j = 0; j < (ny + 1); j++)
	{
		i = nx;
		x[i][j] = Point_List[i + j].Get_x();
		y[i][j] = Point_List[i + j].Get_y();
	}
	for (i = 0; i < (nx + 1); i++)
	{
		j = ny;
		x[i][j] = Point_List[No_Points - (i + j)].Get_x();
		y[i][j] = Point_List[No_Points - (i + j)].Get_y();
	}
	for (j = ny; j > 0; j--)
	{
		i = 0;
		x[i][j] = Point_List[No_Points - (i + j)].Get_x();
		y[i][j] = Point_List[No_Points - (i + j)].Get_y();
	}
}

void Grid::Read_data(string ipfile)
{
}

/* ======================================================================
   Transfinite Interpolation (TFI) for algebraic grid generation.
   Generates interior points from four boundary curves using bilinear
   blending.
   ====================================================================== */
void Grid::TFI()
{
	cout << "Generating Grid using Algebraic Transfinite Interpolation\n";
	cout << "Grid dimensions: " << nx << " x " << ny << endl;
	cout << "Corner (0,0):       (" << x[0][0]           << ", " << y[0][0]           << ")\n";
	cout << "Corner (nx-1,0):    (" << x[nx - 1][0]      << ", " << y[nx - 1][0]      << ")\n";
	cout << "Corner (nx-1,ny-1): (" << x[nx - 1][ny - 1] << ", " << y[nx - 1][ny - 1] << ")\n";
	cout << "Corner (0,ny-1):    (" << x[0][ny - 1]      << ", " << y[0][ny - 1]      << ")\n";

	for (int j = 1; j < ny - 1; j++)
	{
		double eta_j = static_cast<double>(j) / (ny - 1);
		for (int i = 1; i < nx - 1; i++)
		{
			double xi_i = static_cast<double>(i) / (nx - 1);

			double x_proj_eta = (1.0 - eta_j) * x[i][0]      + eta_j * x[i][ny - 1];
			double x_proj_xi  = (1.0 - xi_i)  * x[0][j]      + xi_i  * x[nx - 1][j];
			double x_corners  = (1.0 - xi_i) * (1.0 - eta_j) * x[0][0]
							  + xi_i * (1.0 - eta_j) * x[nx - 1][0]
							  + (1.0 - xi_i) * eta_j * x[0][ny - 1]
							  + xi_i * eta_j * x[nx - 1][ny - 1];

			double y_proj_eta = (1.0 - eta_j) * y[i][0]      + eta_j * y[i][ny - 1];
			double y_proj_xi  = (1.0 - xi_i)  * y[0][j]      + xi_i  * y[nx - 1][j];
			double y_corners  = (1.0 - xi_i) * (1.0 - eta_j) * y[0][0]
							  + xi_i * (1.0 - eta_j) * y[nx - 1][0]
							  + (1.0 - xi_i) * eta_j * y[0][ny - 1]
							  + xi_i * eta_j * y[nx - 1][ny - 1];

			x[i][j] = x_proj_eta + x_proj_xi - x_corners;
			y[i][j] = y_proj_eta + y_proj_xi - y_corners;
		}
	}

	xtemp = x;
	ytemp = y;
	cout << "TFI complete. Ready for elliptic smoothing.\n";
}

/* ======================================================================
   Original Laplace elliptic smoother (P=Q=0).
   Kept for backward compatibility.
   ====================================================================== */
void Grid::Elliptic()
{
	double dxi = 1.0 / (nx - 1);
	double deta = 1.0 / (ny - 1);
	double dxi2 = dxi * dxi;
	double deta2 = deta * deta;
	double dxi_deta_4 = 4.0 * dxi * deta;

	for (int j = 1; j < ny - 1; j++)
	{
		for (int i = 1; i < nx - 1; i++)
		{
			double x_xi  = (x[i + 1][j] - xtemp[i - 1][j]) / (2.0 * dxi);
			double x_eta = (x[i][j + 1] - xtemp[i][j - 1]) / (2.0 * deta);
			double y_xi  = (y[i + 1][j] - ytemp[i - 1][j]) / (2.0 * dxi);
			double y_eta = (y[i][j + 1] - ytemp[i][j - 1]) / (2.0 * deta);

			double alp = x_eta * x_eta + y_eta * y_eta;
			double gam = x_xi  * x_xi  + y_xi  * y_xi;
			double bet = x_xi  * x_eta + y_xi  * y_eta;

			double coeff_B = alp / dxi2;
			double coeff_D = gam / deta2;
			double coeff_C = -bet / dxi_deta_4;
			double coeff_A = coeff_B + coeff_D;

			double cross_x = coeff_C * (x[i + 1][j + 1] - xtemp[i + 1][j - 1]
									   - x[i - 1][j + 1] + xtemp[i - 1][j - 1]);
			double cross_y = coeff_C * (y[i + 1][j + 1] - ytemp[i + 1][j - 1]
									   - y[i - 1][j + 1] + ytemp[i - 1][j - 1]);

			xtemp[i][j] = (1.0 / (2.0 * coeff_A)) *
				(coeff_B * (x[i + 1][j] + xtemp[i - 1][j])
				+ cross_x
				+ coeff_D * (x[i][j + 1] + xtemp[i][j - 1]));

			ytemp[i][j] = (1.0 / (2.0 * coeff_A)) *
				(coeff_B * (y[i + 1][j] + ytemp[i - 1][j])
				+ cross_y
				+ coeff_D * (y[i][j + 1] + ytemp[i][j - 1]));
		}
	}
}

/* ======================================================================
   Poisson elliptic grid generator with source terms P(xi,eta), Q(xi,eta).

   Solves the transformed Poisson equations:
     alpha * x_xixi - 2*beta * x_xieta + gamma * x_etaeta = -J^2*(P*x_xi + Q*x_eta)
     alpha * y_xixi - 2*beta * y_xieta + gamma * y_etaeta = -J^2*(P*y_xi + Q*y_eta)

   where:
     alpha = x_eta^2 + y_eta^2
     beta  = x_xi*x_eta + y_xi*y_eta
     gamma = x_xi^2 + y_xi^2
     J     = x_xi*y_eta - x_eta*y_xi  (Jacobian)

   Source terms P, Q control grid clustering and orthogonality:
   - Thomas-Middlecoff sources: enforce boundary point distribution
   - Boundary attraction sources: force orthogonality at walls
   ====================================================================== */
void Grid::Elliptic_Poisson()
{
	double dxi = 1.0 / (nx - 1);
	double deta = 1.0 / (ny - 1);
	double dxi2 = dxi * dxi;
	double deta2 = deta * deta;

	for (int j = 1; j < ny - 1; j++)
	{
		for (int i = 1; i < nx - 1; i++)
		{
			double x_xi  = (x[i + 1][j] - xtemp[i - 1][j]) / (2.0 * dxi);
			double x_eta = (x[i][j + 1] - xtemp[i][j - 1]) / (2.0 * deta);
			double y_xi  = (y[i + 1][j] - ytemp[i - 1][j]) / (2.0 * dxi);
			double y_eta = (y[i][j + 1] - ytemp[i][j - 1]) / (2.0 * deta);

			double alp = x_eta * x_eta + y_eta * y_eta;
			double gam = x_xi  * x_xi  + y_xi  * y_xi;
			double bet = x_xi  * x_eta + y_xi  * y_eta;

			double J = x_xi * y_eta - x_eta * y_xi;
			double J2 = J * J;

			double Pval = P_source[i][j];
			double Qval = Q_source[i][j];

			double rhs_x = -J2 * (Pval * x_xi + Qval * x_eta);
			double rhs_y = -J2 * (Pval * y_xi + Qval * y_eta);

			double coeff_B = alp / dxi2;
			double coeff_D = gam / deta2;
			double coeff_C = -bet / (4.0 * dxi * deta);
			double coeff_A = coeff_B + coeff_D;

			double cross_x = coeff_C * (x[i + 1][j + 1] - xtemp[i + 1][j - 1]
									   - x[i - 1][j + 1] + xtemp[i - 1][j - 1]);
			double cross_y = coeff_C * (y[i + 1][j + 1] - ytemp[i + 1][j - 1]
									   - y[i - 1][j + 1] + ytemp[i - 1][j - 1]);

			xtemp[i][j] = (1.0 / (2.0 * coeff_A)) *
				(coeff_B * (x[i + 1][j] + xtemp[i - 1][j])
				+ cross_x
				+ coeff_D * (x[i][j + 1] + xtemp[i][j - 1])
				+ rhs_x);

			ytemp[i][j] = (1.0 / (2.0 * coeff_A)) *
				(coeff_B * (y[i + 1][j] + ytemp[i - 1][j])
				+ cross_y
				+ coeff_D * (y[i][j + 1] + ytemp[i][j - 1])
				+ rhs_y);
		}
	}
}

/* ======================================================================
   Thomas-Middlecoff source term computation.
   Computes P(xi,eta) and Q(xi,eta) from the boundary point distribution
   so the interior grid inherits the boundary clustering while decaying
   smoothly toward the interior.

   On each boundary, the source is computed from the second derivative
   of arc-length with respect to the computational coordinate.
   Interior values are interpolated with exponential decay.
   ====================================================================== */
void Grid::ComputeThomasMiddlecoffSources()
{
	P_source.assign(nx, vector<double>(ny, 0.0));
	Q_source.assign(nx, vector<double>(ny, 0.0));

	double dxi  = 1.0 / (nx - 1);
	double deta = 1.0 / (ny - 1);
	double decay = source_params.decay_rate;

	vector<double> P_south(nx, 0.0), P_north(nx, 0.0);
	vector<double> Q_west(ny, 0.0),  Q_east(ny, 0.0);

	// P on j=0 (south) boundary from boundary point spacing
	for (int i = 1; i < nx - 1; i++)
	{
		double x_xi  = (x[i + 1][0] - x[i - 1][0]) / (2.0 * dxi);
		double y_xi  = (y[i + 1][0] - y[i - 1][0]) / (2.0 * dxi);
		double x_xixi = (x[i + 1][0] - 2.0 * x[i][0] + x[i - 1][0]) / (dxi * dxi);
		double y_xixi = (y[i + 1][0] - 2.0 * y[i][0] + y[i - 1][0]) / (dxi * dxi);
		double denom = x_xi * x_xi + y_xi * y_xi;
		if (fabs(denom) > 1e-30)
			P_south[i] = -(x_xi * x_xixi + y_xi * y_xixi) / denom;
	}

	// P on j=ny-1 (north) boundary
	for (int i = 1; i < nx - 1; i++)
	{
		int jb = ny - 1;
		double x_xi  = (x[i + 1][jb] - x[i - 1][jb]) / (2.0 * dxi);
		double y_xi  = (y[i + 1][jb] - y[i - 1][jb]) / (2.0 * dxi);
		double x_xixi = (x[i + 1][jb] - 2.0 * x[i][jb] + x[i - 1][jb]) / (dxi * dxi);
		double y_xixi = (y[i + 1][jb] - 2.0 * y[i][jb] + y[i - 1][jb]) / (dxi * dxi);
		double denom = x_xi * x_xi + y_xi * y_xi;
		if (fabs(denom) > 1e-30)
			P_north[i] = -(x_xi * x_xixi + y_xi * y_xixi) / denom;
	}

	// Q on i=0 (west) boundary
	for (int j = 1; j < ny - 1; j++)
	{
		double x_eta  = (x[0][j + 1] - x[0][j - 1]) / (2.0 * deta);
		double y_eta  = (y[0][j + 1] - y[0][j - 1]) / (2.0 * deta);
		double x_etaeta = (x[0][j + 1] - 2.0 * x[0][j] + x[0][j - 1]) / (deta * deta);
		double y_etaeta = (y[0][j + 1] - 2.0 * y[0][j] + y[0][j - 1]) / (deta * deta);
		double denom = x_eta * x_eta + y_eta * y_eta;
		if (fabs(denom) > 1e-30)
			Q_west[j] = -(x_eta * x_etaeta + y_eta * y_etaeta) / denom;
	}

	// Q on i=nx-1 (east) boundary
	for (int j = 1; j < ny - 1; j++)
	{
		int ib = nx - 1;
		double x_eta  = (x[ib][j + 1] - x[ib][j - 1]) / (2.0 * deta);
		double y_eta  = (y[ib][j + 1] - y[ib][j - 1]) / (2.0 * deta);
		double x_etaeta = (x[ib][j + 1] - 2.0 * x[ib][j] + x[ib][j - 1]) / (deta * deta);
		double y_etaeta = (y[ib][j + 1] - 2.0 * y[ib][j] + y[ib][j - 1]) / (deta * deta);
		double denom = x_eta * x_eta + y_eta * y_eta;
		if (fabs(denom) > 1e-30)
			Q_east[j] = -(x_eta * x_etaeta + y_eta * y_etaeta) / denom;
	}

	// Interpolate into interior with exponential decay
	for (int j = 1; j < ny - 1; j++)
	{
		double eta_j = static_cast<double>(j) / (ny - 1);
		double decay_south = exp(-decay * eta_j * (ny - 1));
		double decay_north = exp(-decay * (1.0 - eta_j) * (ny - 1));

		for (int i = 1; i < nx - 1; i++)
		{
			double xi_i = static_cast<double>(i) / (nx - 1);
			double decay_west = exp(-decay * xi_i * (nx - 1));
			double decay_east = exp(-decay * (1.0 - xi_i) * (nx - 1));

			P_source[i][j] = P_south[i] * decay_south
							+ P_north[i] * decay_north;

			Q_source[i][j] = Q_west[j] * decay_west
							+ Q_east[j] * decay_east;
		}
	}
}

/* ======================================================================
   Boundary attraction source terms.
   Adds additional source terms that attract grid lines toward wall
   boundaries, enforcing near-orthogonality and clustering near walls.
   Uses the boundary layer parameters (first_cell_height, growth_rate)
   to determine attraction strength.
   ====================================================================== */
void Grid::ComputeBoundaryAttractionSources()
{
	double amp = source_params.attraction_strength;

	// South boundary (j=0): attract in eta direction -> Q source
	if (bl_south.enabled)
	{
		double wall_decay = -log(0.01) / bl_south.num_layers;
		for (int j = 1; j < ny - 1; j++)
		{
			double dist = static_cast<double>(j);
			double attraction = amp * exp(-wall_decay * dist);
			for (int i = 1; i < nx - 1; i++)
				Q_source[i][j] -= attraction;
		}
	}

	// North boundary (j=ny-1): attract in eta direction
	if (bl_north.enabled)
	{
		double wall_decay = -log(0.01) / bl_north.num_layers;
		for (int j = 1; j < ny - 1; j++)
		{
			double dist = static_cast<double>(ny - 1 - j);
			double attraction = amp * exp(-wall_decay * dist);
			for (int i = 1; i < nx - 1; i++)
				Q_source[i][j] += attraction;
		}
	}

	// West boundary (i=0): attract in xi direction -> P source
	if (bl_west.enabled)
	{
		double wall_decay = -log(0.01) / bl_west.num_layers;
		for (int i = 1; i < nx - 1; i++)
		{
			double dist = static_cast<double>(i);
			double attraction = amp * exp(-wall_decay * dist);
			for (int j = 1; j < ny - 1; j++)
				P_source[i][j] -= attraction;
		}
	}

	// East boundary (i=nx-1): attract in xi direction
	if (bl_east.enabled)
	{
		double wall_decay = -log(0.01) / bl_east.num_layers;
		for (int i = 1; i < nx - 1; i++)
		{
			double dist = static_cast<double>(nx - 1 - i);
			double attraction = amp * exp(-wall_decay * dist);
			for (int j = 1; j < ny - 1; j++)
				P_source[i][j] += attraction;
		}
	}
}

void Grid::SetSourceTermParams(const SourceTermParams &params)
{
	source_params = params;
}

void Grid::SetBoundaryLayerParams(int boundary_id, const BoundaryLayerParams &params)
{
	switch (boundary_id)
	{
		case 0: bl_south = params; break;
		case 1: bl_east  = params; break;
		case 2: bl_north = params; break;
		case 3: bl_west  = params; break;
		default:
			cerr << "Invalid boundary_id: " << boundary_id
				 << " (0=south, 1=east, 2=north, 3=west)\n";
	}
}

void Grid::Periodic_Boundary_Condition()
{
	double dxi  = 1.0 / (nx - 1);
	double deta = 1.0 / (ny - 1);

	int i = nx - 1;
	for (int j = 1; j < ny - 1; j++)
	{
		double x_xi  = (x[1][j] - x[i - 1][j]) / (2.0 * dxi);
		double x_eta = (x[i][j + 1] - x[i][j - 1]) / (2.0 * deta);
		double y_xi  = (y[1][j] - y[i - 1][j]) / (2.0 * dxi);
		double y_eta = (y[i][j + 1] - y[i][j - 1]) / (2.0 * deta);

		double alp = x_eta * x_eta + y_eta * y_eta;
		double gam = x_xi  * x_xi  + y_xi  * y_xi;
		double bet = x_xi  * x_eta + y_xi  * y_eta;

		double coeff_B = alp / (dxi * dxi);
		double coeff_D = gam / (deta * deta);
		double coeff_C = -bet / (4.0 * dxi * deta);
		double coeff_A = coeff_B + coeff_D;

		xtemp[i][j] = (1.0 / (2.0 * coeff_A)) *
			(coeff_B * (x[1][j] + x[i - 1][j])
			+ coeff_C * (x[1][j + 1] - x[1][j - 1] - x[i - 1][j + 1] + x[i - 1][j - 1])
			+ coeff_D * (x[i][j + 1] + x[i][j - 1]));

		ytemp[i][j] = (1.0 / (2.0 * coeff_A)) *
			(coeff_B * (y[1][j] + y[i - 1][j])
			+ coeff_C * (y[1][j + 1] - y[1][j - 1] - y[i - 1][j + 1] + y[i - 1][j - 1])
			+ coeff_D * (y[i][j + 1] + y[i][j - 1]));

		x[0][j] = xtemp[i][j];
		y[0][j] = ytemp[i][j];
	}
}

/* ======================================================================
   Error estimation and grid update (standard Gauss-Seidel).
   Uses L2 norm of absolute change to avoid division-by-zero when
   grid coordinates pass through zero.
   ====================================================================== */
void Grid::Estimate_Error_and_Update()
{
	ertot = 0.0;
	for (int j = 1; j < ny - 1; j++)
	{
		for (int i = 1; i < nx - 1; i++)
		{
			double dx = xtemp[i][j] - x[i][j];
			double dy = ytemp[i][j] - y[i][j];
			ertot += dx * dx + dy * dy;

			x[i][j] = xtemp[i][j];
			y[i][j] = ytemp[i][j];
		}
	}
	ertot = sqrt(ertot / ((nx - 2) * (ny - 2)));
}

/* ======================================================================
   Error estimation with Successive Over-Relaxation (SOR).
   omega in (1,2) for over-relaxation; omega=1 is standard Gauss-Seidel.
   Typical optimal values are 1.4-1.8 for elliptic grid generation.
   ====================================================================== */
void Grid::Estimate_Error_and_Update_SOR(double omega)
{
	ertot = 0.0;
	for (int j = 1; j < ny - 1; j++)
	{
		for (int i = 1; i < nx - 1; i++)
		{
			double dx = xtemp[i][j] - x[i][j];
			double dy = ytemp[i][j] - y[i][j];
			ertot += dx * dx + dy * dy;

			x[i][j] = x[i][j] + omega * dx;
			y[i][j] = y[i][j] + omega * dy;
			xtemp[i][j] = x[i][j];
			ytemp[i][j] = y[i][j];
		}
	}
	ertot = sqrt(ertot / ((nx - 2) * (ny - 2)));
}

/* ======================================================================
   Original grid generation (backward compatible).
   ====================================================================== */
void Grid::Generate_Grid(bool &With_Periodic_Domain)
{
	noit = 0;
	TFI();
	do
	{
		Elliptic();
		if (With_Periodic_Domain)
			Periodic_Boundary_Condition();
		Estimate_Error_and_Update();
		noit++;
		if (noit % 1000 == 0)
			cout << "Iteration " << noit << "\tL2 error = " << ertot << endl;
	} while (noit < 100000);

	cout << "Elliptic Grid Generated\t......Generating Grid List\n";
	Generate_Grid_List();
}

void Grid::Generate_Grid(bool &With_Periodic_Domain, bool &EnableElliptic, int &iterations)
{
	noit = 0;
	TFI();
	cout << "Elliptic solver: " << (EnableElliptic ? "ON" : "OFF")
		 << ", iterations = " << iterations << endl;

	if (EnableElliptic)
	{
		if (source_params.enabled)
		{
			if (source_params.use_thomas_middlecoff)
				ComputeThomasMiddlecoffSources();
			if (source_params.use_boundary_attraction)
				ComputeBoundaryAttractionSources();
		}

		do
		{
			if (source_params.enabled)
				Elliptic_Poisson();
			else
				Elliptic();

			if (With_Periodic_Domain)
				Periodic_Boundary_Condition();
			Estimate_Error_and_Update();
			noit++;

			if (noit % 1000 == 0)
				cout << "Iteration " << noit << "\tL2 error = " << ertot << endl;

		} while (noit < iterations);

		cout << "Elliptic Grid Generated. Final L2 error = " << ertot << endl;
	}
	Generate_Grid_List();
}

/* ======================================================================
   Enhanced grid generation with SOR and convergence tolerance.
   Uses Poisson equations with source terms when enabled.
   Recomputes Thomas-Middlecoff sources periodically as interior moves.

   Parameters:
     periodic        - enable periodic BC on xi-max/xi-min
     enableElliptic  - enable elliptic smoothing (false = TFI only)
     iterations      - maximum number of iterations
     omega           - SOR relaxation factor (1.0 = Gauss-Seidel, 1.4-1.8 typical)
     convergence_tol - stop when L2 error drops below this
   ====================================================================== */
void Grid::Generate_Grid(bool &periodic, bool &enableElliptic, int &iterations,
						 double omega, double convergence_tol)
{
	noit = 0;
	TFI();

	cout << "Enhanced elliptic generator:\n"
		 << "  Elliptic: " << (enableElliptic ? "ON" : "OFF") << "\n"
		 << "  Max iterations: " << iterations << "\n"
		 << "  SOR omega: " << omega << "\n"
		 << "  Convergence tolerance: " << convergence_tol << "\n"
		 << "  Source terms: " << (source_params.enabled ? "ON" : "OFF") << "\n";

	if (!enableElliptic)
	{
		Generate_Grid_List();
		return;
	}

	if (source_params.enabled)
	{
		if (source_params.use_thomas_middlecoff)
			ComputeThomasMiddlecoffSources();
		if (source_params.use_boundary_attraction)
			ComputeBoundaryAttractionSources();
	}

	int source_recompute_interval = max(100, iterations / 20);

	do
	{
		if (source_params.enabled && source_params.use_thomas_middlecoff
			&& noit > 0 && noit % source_recompute_interval == 0)
		{
			ComputeThomasMiddlecoffSources();
			if (source_params.use_boundary_attraction)
				ComputeBoundaryAttractionSources();
		}

		if (source_params.enabled)
			Elliptic_Poisson();
		else
			Elliptic();

		if (periodic)
			Periodic_Boundary_Condition();

		Estimate_Error_and_Update_SOR(omega);
		noit++;

		if (noit % 500 == 0)
			cout << "Iteration " << noit << "\tL2 error = " << scientific << ertot << endl;

	} while (noit < iterations && ertot > convergence_tol);

	cout << fixed;
	cout << "Elliptic generation complete after " << noit << " iterations.\n"
		 << "  Final L2 error: " << scientific << ertot << endl;

	double min_jac, max_ar, max_skew, min_orth;
	ComputeGridQuality(min_jac, max_ar, max_skew, min_orth);
	cout << fixed << setprecision(4);
	cout << "Grid quality metrics:\n"
		 << "  Min Jacobian: " << min_jac << "\n"
		 << "  Max aspect ratio: " << max_ar << "\n"
		 << "  Max skewness: " << max_skew << "\n"
		 << "  Min orthogonality (deg): " << min_orth << "\n";

	Generate_Grid_List();
}

/* ======================================================================
   Grid quality metrics computation.
   Evaluates Jacobian, aspect ratio, skewness, and orthogonality
   for every interior cell in the structured grid.
   ====================================================================== */
void Grid::ComputeGridQuality(double &min_jac, double &max_ar,
							  double &max_skew, double &min_orth) const
{
	min_jac  = 1e30;
	max_ar   = 0.0;
	max_skew = 0.0;
	min_orth = 180.0;

	double dxi  = 1.0 / (nx - 1);
	double deta = 1.0 / (ny - 1);

	for (int j = 1; j < ny - 1; j++)
	{
		for (int i = 1; i < nx - 1; i++)
		{
			double x_xi  = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxi);
			double x_eta = (x[i][j + 1] - x[i][j - 1]) / (2.0 * deta);
			double y_xi  = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxi);
			double y_eta = (y[i][j + 1] - y[i][j - 1]) / (2.0 * deta);

			double J = x_xi * y_eta - x_eta * y_xi;
			if (J < min_jac) min_jac = J;

			double len_xi  = sqrt(x_xi * x_xi + y_xi * y_xi);
			double len_eta = sqrt(x_eta * x_eta + y_eta * y_eta);
			double ar = (len_xi > len_eta) ? len_xi / max(len_eta, 1e-30)
											: len_eta / max(len_xi, 1e-30);
			if (ar > max_ar) max_ar = ar;

			double dot = x_xi * x_eta + y_xi * y_eta;
			double cos_angle = dot / max(len_xi * len_eta, 1e-30);
			cos_angle = max(-1.0, min(1.0, cos_angle));
			double angle_deg = acos(fabs(cos_angle)) * 180.0 / M_PI;
			if (angle_deg < min_orth) min_orth = angle_deg;

			double skew = fabs(90.0 - angle_deg) / 90.0;
			if (skew > max_skew) max_skew = skew;
		}
	}
}

void Grid::Generate_Grid_List()
{
	Grid_List.clear();
	Grid_List.reserve(nx * ny);
	Point Temp_Point;
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			Temp_Point(x[i][j], y[i][j], 0.0);
			Grid_List.push_back(Temp_Point);
		}
	}
	cout << "Grid list generated: " << Grid_List.size() << " points\n";
}

vector<Point> &Grid::Get_Grid_List()
{
	return Grid_List;
}

/* ======================================================================
   Uniform Z-stacking (backward compatible).
   ====================================================================== */
void Grid::Stack_Grid(double &length, int &nop, vector<Point> &Stacked_Grid_List)
{
	cout << "Stacking grid:\n"
		 << "  Plane size: " << Grid_List.size() << "\n"
		 << "  Stack length: " << length << "\n"
		 << "  Points in Z: " << nop << endl;

	double delz = length / (nop - 1);
	Stacked_Grid_List.reserve(Grid_List.size() * nop);

	Point p;
	for (int i = 0; i < nop; i++)
	{
		for (unsigned int j = 0; j < Grid_List.size(); j++)
		{
			p(Grid_List[j].Get_x(), Grid_List[j].Get_y(),
			  i * delz + Grid_List[j].Get_z());
			Stacked_Grid_List.push_back(p);
		}
	}
	cout << "Stacking complete. Total points: " << Stacked_Grid_List.size() << endl;
}

/* ======================================================================
   Geometric-ratio Z-stacking for boundary layer resolution.
   First layer thickness is length*(1-r)/(1-r^nop) where r=growth_ratio.
   ====================================================================== */
void Grid::Stack_Grid(double &length, int &nop, vector<Point> &Stacked_Grid_List,
					  double growth_ratio)
{
	cout << "Stacking grid with geometric stretching:\n"
		 << "  Plane size: " << Grid_List.size() << "\n"
		 << "  Stack length: " << length << "\n"
		 << "  Points in Z: " << nop << "\n"
		 << "  Growth ratio: " << growth_ratio << endl;

	vector<double> z_coords(nop, 0.0);

	if (fabs(growth_ratio - 1.0) < 1e-10)
	{
		double delz = length / (nop - 1);
		for (int i = 0; i < nop; i++)
			z_coords[i] = i * delz;
	}
	else
	{
		double first_dz = length * (1.0 - growth_ratio)
						/ (1.0 - pow(growth_ratio, nop - 1));
		z_coords[0] = 0.0;
		double dz = first_dz;
		for (int i = 1; i < nop; i++)
		{
			z_coords[i] = z_coords[i - 1] + dz;
			dz *= growth_ratio;
		}
		double scale = length / z_coords[nop - 1];
		for (int i = 1; i < nop; i++)
			z_coords[i] *= scale;
	}

	cout << "  First layer dz: " << z_coords[1] - z_coords[0] << "\n"
		 << "  Last layer dz:  " << z_coords[nop - 1] - z_coords[nop - 2] << endl;

	Stacked_Grid_List.reserve(Grid_List.size() * nop);
	Point p;
	for (int i = 0; i < nop; i++)
	{
		for (unsigned int j = 0; j < Grid_List.size(); j++)
		{
			p(Grid_List[j].Get_x(), Grid_List[j].Get_y(),
			  z_coords[i] + Grid_List[j].Get_z());
			Stacked_Grid_List.push_back(p);
		}
	}
	cout << "Stacking complete. Total points: " << Stacked_Grid_List.size() << endl;
}
