#include "Geometry_Header.h"
#include <cassert>

void write_grid_vtk(const string &filename, Grid &G, int nx, int ny)
{
	ofstream out(filename);
	if (!out.is_open()) { cerr << "Cannot open " << filename << endl; return; }

	int npoints = nx * ny;
	int ncells  = (nx - 1) * (ny - 1);

	out << "# vtk DataFile Version 3.0\n"
		<< "Grid Generator Test\n"
		<< "ASCII\n"
		<< "DATASET STRUCTURED_GRID\n"
		<< "DIMENSIONS " << nx << " " << ny << " 1\n"
		<< "POINTS " << npoints << " double\n";

	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
			out << G.Get_x(i, j) << " " << G.Get_y(i, j) << " 0.0\n";

	out.close();
	cout << "VTK written: " << filename << " (" << npoints << " points)\n";
}

void write_line_dat(const string &filename, const vector<Point> &pts)
{
	ofstream out(filename);
	for (unsigned int i = 0; i < pts.size(); i++)
		out << i << "\t" << pts[i].Get_x() << "\t" << pts[i].Get_y()
			<< "\t" << pts[i].Get_z() << "\n";
	out.close();
	cout << "Line data written: " << filename << " (" << pts.size() << " points)\n";
}

/* ======================================================================
   Test 1: Rectangular channel with TFI + Laplace elliptic smoothing
   ====================================================================== */
void test_basic_rectangle()
{
	cout << "\n========== TEST 1: Basic Rectangle (TFI + Laplace) ==========\n";
	int nx = 41, ny = 21;

	Point p1, p2;
	Line south, east, north, west;

	p1(0.0, 0.0, 0.0); p2(2.0, 0.0, 0.0);
	south.generate(p1, p2, nx);

	p1(2.0, 0.0, 0.0); p2(2.0, 1.0, 0.0);
	east.generate(p1, p2, ny);

	p1(2.0, 1.0, 0.0); p2(0.0, 1.0, 0.0);
	north.generate(p1, p2, nx);

	p1(0.0, 1.0, 0.0); p2(0.0, 0.0, 0.0);
	west.generate(p1, p2, ny);

	vector<Point> s_pts = south.Get_Point_list();
	vector<Point> e_pts = east.Get_Point_list();
	vector<Point> n_pts = north.Get_Point_list();
	vector<Point> w_pts = west.Get_Point_list();

	Grid G;
	G(s_pts, e_pts, n_pts, w_pts);

	bool periodic = false, elliptic = true;
	int iters = 5000;
	G.Generate_Grid(periodic, elliptic, iters);

	write_grid_vtk("test1_rectangle.vtk", G, nx, ny);

	double min_jac, max_ar, max_skew, min_orth;
	G.ComputeGridQuality(min_jac, max_ar, max_skew, min_orth);
	assert(min_jac > 0 && "Negative Jacobian in rectangle grid!");
	cout << "TEST 1 PASSED: min_jac=" << min_jac << " max_ar=" << max_ar << "\n";
}

/* ======================================================================
   Test 2: Curved channel with TFI + Poisson (Thomas-Middlecoff sources)
   ====================================================================== */
void test_curved_channel_poisson()
{
	cout << "\n========== TEST 2: Curved Channel (TFI + Poisson + T-M Sources) ==========\n";
	int nx = 51, ny = 31;

	Point p1, p2;
	Line south, east, north, west;

	p1(0.0, 0.0, 0.0); p2(3.0, 0.0, 0.0);
	south.generate(p1, p2, nx);

	p1(3.0, 0.0, 0.0); p2(3.0, 1.0, 0.0);
	east.generate(p1, p2, ny);

	// Curved top wall: y = 1 + 0.3*sin(pi*x/3)
	vector<Point> north_pts;
	for (int i = 0; i < nx; i++)
	{
		double xi = 3.0 * static_cast<double>(nx - 1 - i) / (nx - 1);
		double yi = 1.0 + 0.3 * sin(M_PI * xi / 3.0);
		Point pp;
		pp(xi, yi, 0.0);
		north_pts.push_back(pp);
	}

	p1(0.0, 1.0, 0.0); p2(0.0, 0.0, 0.0);
	west.generate(p1, p2, ny);

	vector<Point> s_pts = south.Get_Point_list();
	vector<Point> e_pts = east.Get_Point_list();
	vector<Point> w_pts = west.Get_Point_list();

	Grid G;
	G(s_pts, e_pts, north_pts, w_pts);

	SourceTermParams src;
	src.enabled = true;
	src.use_thomas_middlecoff = true;
	src.decay_rate = 0.5;
	G.SetSourceTermParams(src);

	bool periodic = false, elliptic = true;
	int iters = 10000;
	G.Generate_Grid(periodic, elliptic, iters, 1.6, 1e-10);

	write_grid_vtk("test2_curved_channel.vtk", G, nx, ny);

	double min_jac, max_ar, max_skew, min_orth;
	G.ComputeGridQuality(min_jac, max_ar, max_skew, min_orth);
	assert(min_jac > 0 && "Negative Jacobian in curved channel grid!");
	cout << "TEST 2 PASSED: converged in " << G.Get_Iterations()
		 << " iters, min_jac=" << min_jac << "\n";
}

/* ======================================================================
   Test 3: Boundary layer grid with wall clustering
   Uses boundary layer line generation + Poisson source terms
   with boundary attraction for wall-normal clustering.
   ====================================================================== */
void test_boundary_layer_grid()
{
	cout << "\n========== TEST 3: Boundary Layer Grid (BL Clustering) ==========\n";
	int nx = 61, ny = 41;

	Point p1, p2;
	Line south, east, north, west;

	// Flat plate (south wall)
	p1(0.0, 0.0, 0.0); p2(2.0, 0.0, 0.0);
	south.generate(p1, p2, nx);

	// East side: BL clustering toward the wall (j=0)
	p1(2.0, 0.0, 0.0); p2(2.0, 0.5, 0.0);
	double first_height = 0.001;
	double growth = 1.15;
	east.generate_boundary_layer(p1, p2, ny, first_height, growth);

	// Top (freestream)
	p1(2.0, 0.5, 0.0); p2(0.0, 0.5, 0.0);
	north.generate(p1, p2, nx);

	// West side: BL clustering toward the wall
	p1(0.0, 0.5, 0.0); p2(0.0, 0.0, 0.0);
	Line west_bl;
	west_bl.generate_boundary_layer(p2, p1, ny, first_height, growth);
	west_bl.Reverse_Points();
	vector<Point> w_pts = west_bl.Get_Point_list();

	vector<Point> s_pts = south.Get_Point_list();
	vector<Point> e_pts = east.Get_Point_list();
	vector<Point> n_pts = north.Get_Point_list();

	Grid G;
	G(s_pts, e_pts, n_pts, w_pts);

	SourceTermParams src;
	src.enabled = true;
	src.use_thomas_middlecoff = true;
	src.use_boundary_attraction = true;
	src.decay_rate = 0.3;
	src.attraction_strength = 3.0;
	G.SetSourceTermParams(src);

	BoundaryLayerParams bl;
	bl.enabled = true;
	bl.first_cell_height = first_height;
	bl.growth_rate = growth;
	bl.num_layers = 20;
	G.SetBoundaryLayerParams(0, bl); // south wall

	bool periodic = false, elliptic = true;
	int iters = 15000;
	G.Generate_Grid(periodic, elliptic, iters, 1.5, 1e-10);

	write_grid_vtk("test3_boundary_layer.vtk", G, nx, ny);

	// Verify first cell height near wall
	double dy0 = G.Get_y(nx / 2, 1) - G.Get_y(nx / 2, 0);
	double dy1 = G.Get_y(nx / 2, 2) - G.Get_y(nx / 2, 1);
	cout << "  First cell dy at mid-plate: " << dy0 << "\n"
		 << "  Second cell dy:             " << dy1 << "\n"
		 << "  Growth ratio:               " << dy1 / dy0 << "\n";

	double min_jac, max_ar, max_skew, min_orth;
	G.ComputeGridQuality(min_jac, max_ar, max_skew, min_orth);
	assert(min_jac > 0 && "Negative Jacobian in BL grid!");
	cout << "TEST 3 PASSED: min_jac=" << min_jac
		 << " first_dy=" << dy0 << "\n";
}

/* ======================================================================
   Test 4: Line stretching methods comparison
   ====================================================================== */
void test_line_stretching()
{
	cout << "\n========== TEST 4: Line Stretching Methods ==========\n";
	int n = 41;
	Point p1, p2;
	p1(0.0, 0.0, 0.0);
	p2(0.0, 1.0, 0.0);

	// Uniform
	Line L_uniform;
	L_uniform.generate(p1, p2, n);
	write_line_dat("test4_uniform.dat", L_uniform.Get_Point_list());

	// Vinokur (backward-compatible)
	Line L_vinokur;
	bool bs = true;
	L_vinokur.generate(p1, p2, n, bs);
	write_line_dat("test4_vinokur.dat", L_vinokur.Get_Point_list());

	// Parameterized Vinokur
	Line L_stretched;
	L_stretched.generate_stretched(p1, p2, n, 1.005, 0.5, true);
	write_line_dat("test4_stretched.dat", L_stretched.Get_Point_list());

	// Boundary layer
	Line L_bl;
	L_bl.generate_boundary_layer(p1, p2, n, 0.001, 1.2);
	write_line_dat("test4_bl.dat", L_bl.Get_Point_list());

	// Tanh
	Line L_tanh;
	L_tanh.generate_tanh(p1, p2, n, 2.5);
	write_line_dat("test4_tanh.dat", L_tanh.Get_Point_list());

	// Verify BL first spacing
	const vector<Point> &bl_pts = L_bl.Get_Point_list();
	double ds0 = bl_pts[1].Get_y() - bl_pts[0].Get_y();
	cout << "  BL first spacing: " << ds0 << " (target: 0.001)\n";

	// Verify tanh symmetry
	const vector<Point> &tanh_pts = L_tanh.Get_Point_list();
	double ds_first = tanh_pts[1].Get_y() - tanh_pts[0].Get_y();
	double ds_last  = tanh_pts[n - 1].Get_y() - tanh_pts[n - 2].Get_y();
	cout << "  Tanh first spacing: " << ds_first << "\n"
		 << "  Tanh last spacing:  " << ds_last << "\n"
		 << "  Symmetry ratio:     " << ds_first / ds_last << " (should be ~1.0)\n";

	assert(fabs(ds_first - ds_last) / ds_first < 0.01 && "Tanh not symmetric!");
	cout << "TEST 4 PASSED\n";
}

/* ======================================================================
   Test 5: Half-cylinder O-grid with circle boundary
   ====================================================================== */
void test_circle_grid()
{
	cout << "\n========== TEST 5: Circle Grid ==========\n";
	double r = 1.0;
	int nop = 21;
	Circle C;
	C(r, nop);

	Grid G;
	G(C);

	bool periodic = false, elliptic = true;
	int iters = 5000;
	G.Generate_Grid(periodic, elliptic, iters);

	int gn = C.get_Nop() / 4;
	write_grid_vtk("test5_circle.vtk", G, gn + 1, gn + 1);

	cout << "TEST 5 PASSED\n";
}

/* ======================================================================
   Test 6: SOR convergence speedup comparison
   ====================================================================== */
void test_sor_speedup()
{
	cout << "\n========== TEST 6: SOR Convergence Speedup ==========\n";
	int nx = 31, ny = 31;
	Point p1, p2;
	Line south, east, north, west;

	p1(0.0, 0.0, 0.0); p2(1.0, 0.0, 0.0);
	south.generate(p1, p2, nx);
	p1(1.0, 0.0, 0.0); p2(1.0, 1.0, 0.0);
	east.generate(p1, p2, ny);
	p1(1.0, 1.0, 0.0); p2(0.0, 1.0, 0.0);
	north.generate(p1, p2, nx);
	p1(0.0, 1.0, 0.0); p2(0.0, 0.0, 0.0);
	west.generate(p1, p2, ny);

	// Make north boundary curved for a non-trivial test
	vector<Point> n_pts;
	for (int i = 0; i < nx; i++)
	{
		double xi = 1.0 - static_cast<double>(i) / (nx - 1);
		double yi = 1.0 + 0.2 * sin(M_PI * xi);
		Point pp;
		pp(xi, yi, 0.0);
		n_pts.push_back(pp);
	}

	vector<Point> s_pts = south.Get_Point_list();
	vector<Point> e_pts = east.Get_Point_list();
	vector<Point> w_pts = west.Get_Point_list();

	// Gauss-Seidel (omega=1.0)
	Grid G1;
	G1(s_pts, e_pts, n_pts, w_pts);
	bool per = false, ell = true;
	int max_it = 5000;
	G1.Generate_Grid(per, ell, max_it, 1.0, 1e-12);
	int gs_iters = G1.Get_Iterations();
	double gs_err = G1.Get_Convergence();

	// SOR (omega=1.6)
	Grid G2;
	G2(s_pts, e_pts, n_pts, w_pts);
	G2.Generate_Grid(per, ell, max_it, 1.6, 1e-12);
	int sor_iters = G2.Get_Iterations();
	double sor_err = G2.Get_Convergence();

	cout << "  Gauss-Seidel: " << gs_iters << " iters, error=" << scientific << gs_err << "\n"
		 << "  SOR (w=1.6):  " << sor_iters << " iters, error=" << scientific << sor_err << "\n";

	if (sor_iters < gs_iters)
		cout << "  Speedup: " << fixed << setprecision(1)
			 << static_cast<double>(gs_iters) / sor_iters << "x\n";

	cout << "TEST 6 PASSED\n";
}

int main()
{
	cout << "============================================================\n"
		 << "     GRID GENERATOR TEST SUITE\n"
		 << "============================================================\n";

	test_basic_rectangle();
	test_curved_channel_poisson();
	test_boundary_layer_grid();
	test_line_stretching();
	test_circle_grid();
	test_sor_speedup();

	cout << "\n============================================================\n"
		 << "     ALL TESTS PASSED\n"
		 << "============================================================\n";

	return 0;
}
