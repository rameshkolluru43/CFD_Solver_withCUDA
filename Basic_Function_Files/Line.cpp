#include "Geometry_Header.h"

Line::Line()
{
	Start_Point(0.0, 0.0, 0.0);
	End_Point(0.0, 0.0, 0.0);
	ds_x = 0.0;
	ds_y = 0.0;
	ds_z = 0.0;
	Length = 0.0;
	no_of_points = 0;
}

const Point &Line::get_Point(int i) const
{
	return point_list[i];
}

/* ======================================================================
   Uniform point distribution along a line.
   ====================================================================== */
void Line::generate(Point &sp, Point &ep, int &n)
{
	Point Temp_point;
	no_of_points = n;
	point_list.clear();
	point_list.reserve(n);

	ds_x = ep.Get_x() - sp.Get_x();
	ds_y = ep.Get_y() - sp.Get_y();
	ds_z = ep.Get_z() - sp.Get_z();
	Length = sqrt(ds_x * ds_x + ds_y * ds_y + ds_z * ds_z);

	delx = ds_x / (n - 1);
	dely = ds_y / (n - 1);
	delz = ds_z / (n - 1);

	for (int i = 0; i < n; i++)
	{
		Temp_point(sp.Get_x() + i * delx,
				   sp.Get_y() + i * dely,
				   sp.Get_z() + i * delz);
		point_list.push_back(Temp_point);
	}
	cout << "Created Line (uniform): " << point_list.size() << " points, L=" << Length << endl;
}

/* ======================================================================
   Vinokur stretching (backward compatible, hardcoded beta=1.01).
   ====================================================================== */
void Line::generate(Point &sp, Point &ep, int &n, bool &both_sides)
{
	double beta_v = 1.01, alpha_v, a, b, y_val, eta_v;
	no_of_points = n;
	point_list.clear();
	point_list.reserve(n);

	ds_x = ep.Get_x() - sp.Get_x();
	ds_y = ep.Get_y() - sp.Get_y();
	ds_z = ep.Get_z() - sp.Get_z();
	Length = sqrt(ds_x * ds_x + ds_y * ds_y + ds_z * ds_z);

	del_s = Length / (n - 1);
	vector<double> s_dist(n, 0.0);

	for (int j = 0; j < n; j++)
	{
		eta_v = static_cast<double>(j) / (n - 1);

		if (both_sides)
		{
			alpha_v = 0.5;
			a = (eta_v - alpha_v) / (1.0 - alpha_v);
			b = (beta_v + 1.0) / (beta_v - 1.0);
			y_val = Length * (((2 * alpha_v + beta_v) * pow(b, a) + 2 * alpha_v - beta_v)
					/ ((2 * alpha_v + 1.0) * (pow(b, a) + 1.0)));
		}
		else
		{
			alpha_v = 0.01;
			a = 1.0 - eta_v;
			b = (beta_v + 1.0) / (beta_v - 1.0);
			y_val = Length * ((beta_v + 1.0) - (beta_v - 1.0) * pow(b, a))
					/ (pow(b, a) + 1.0);
		}
		s_dist[j] = y_val;
	}

	distribute_points(sp, s_dist);
	cout << "Created Line (Vinokur): " << point_list.size() << " points, L=" << Length << endl;
}

/* ======================================================================
   Parameterized Vinokur stretching.
   beta_param: stretching intensity (>1.0, closer to 1 = more stretching)
   alpha_param: clustering location (0.5 = both ends, 0.0 = start only)
   both_sides: true = cluster at both ends, false = cluster at start
   ====================================================================== */
void Line::generate_stretched(Point &sp, Point &ep, int &n,
							  double beta_param, double alpha_param,
							  bool both_sides)
{
	Point Temp_point;
	no_of_points = n;
	point_list.clear();
	point_list.reserve(n);

	ds_x = ep.Get_x() - sp.Get_x();
	ds_y = ep.Get_y() - sp.Get_y();
	ds_z = ep.Get_z() - sp.Get_z();
	Length = sqrt(ds_x * ds_x + ds_y * ds_y + ds_z * ds_z);

	del_s = Length / (n - 1);

	vector<double> s_dist(n, 0.0);
	for (int j = 0; j < n; j++)
	{
		double eta_v = static_cast<double>(j) / (n - 1);
		double y_val;

		if (both_sides)
		{
			double a = (eta_v - alpha_param) / (1.0 - alpha_param);
			double b = (beta_param + 1.0) / (beta_param - 1.0);
			y_val = Length * (((2.0 * alpha_param + beta_param) * pow(b, a)
					+ 2.0 * alpha_param - beta_param)
					/ ((2.0 * alpha_param + 1.0) * (pow(b, a) + 1.0)));
		}
		else
		{
			double a = 1.0 - eta_v;
			double b = (beta_param + 1.0) / (beta_param - 1.0);
			y_val = Length * ((beta_param + 1.0) - (beta_param - 1.0) * pow(b, a))
					/ (pow(b, a) + 1.0);
		}
		s_dist[j] = y_val;
	}

	distribute_points(sp, s_dist);

	cout << "Created Line (stretched): " << point_list.size() << " points\n"
		 << "  beta=" << beta_param << " alpha=" << alpha_param
		 << " both_sides=" << both_sides << "\n"
		 << "  First spacing: " << s_dist[1] << "\n"
		 << "  Last spacing:  " << s_dist[n - 1] - s_dist[n - 2] << endl;
}

/* ======================================================================
   Boundary layer distribution using geometric growth.
   first_cell_height: physical height of the first cell
   growth_rate: ratio of successive cell heights (typically 1.1 - 1.3)

   The distribution grows geometrically from the start point up to the
   point where the total length is reached. If the geometric series
   doesn't reach the full length, the remaining points are uniformly
   distributed.
   ====================================================================== */
void Line::generate_boundary_layer(Point &sp, Point &ep, int &n,
								   double first_cell_height,
								   double growth_rate)
{
	Point Temp_point;
	no_of_points = n;
	point_list.clear();
	point_list.reserve(n);

	ds_x = ep.Get_x() - sp.Get_x();
	ds_y = ep.Get_y() - sp.Get_y();
	ds_z = ep.Get_z() - sp.Get_z();
	Length = sqrt(ds_x * ds_x + ds_y * ds_y + ds_z * ds_z);

	vector<double> s_dist(n, 0.0);
	s_dist[0] = 0.0;

	if (fabs(growth_rate - 1.0) < 1e-10)
	{
		double ds = Length / (n - 1);
		for (int i = 1; i < n; i++)
			s_dist[i] = i * ds;
	}
	else
	{
		double dh = first_cell_height;
		double total_geometric = first_cell_height * (pow(growth_rate, n - 1) - 1.0)
								/ (growth_rate - 1.0);

		if (total_geometric <= Length)
		{
			double scale = Length / total_geometric;
			dh = first_cell_height * scale;
			s_dist[1] = dh;
			for (int i = 2; i < n; i++)
			{
				dh *= growth_rate;
				s_dist[i] = s_dist[i - 1] + dh;
			}
		}
		else
		{
			dh = first_cell_height;
			for (int i = 1; i < n; i++)
			{
				s_dist[i] = s_dist[i - 1] + dh;
				if (s_dist[i] >= Length)
				{
					double uniform_ds = (Length - s_dist[i - 1]) / (n - i);
					for (int k = i; k < n; k++)
						s_dist[k] = s_dist[i - 1] + (k - i + 1) * uniform_ds;
					break;
				}
				dh *= growth_rate;
			}
			double scale = Length / s_dist[n - 1];
			for (int i = 1; i < n; i++)
				s_dist[i] *= scale;
		}
	}

	distribute_points(sp, s_dist);

	cout << "Created Line (boundary layer): " << point_list.size() << " points\n"
		 << "  First cell height: " << first_cell_height << "\n"
		 << "  Growth rate: " << growth_rate << "\n"
		 << "  Actual first dh: " << s_dist[1] << "\n"
		 << "  Actual last dh:  " << s_dist[n - 1] - s_dist[n - 2] << endl;
}

/* ======================================================================
   Hyperbolic tangent stretching.
   stretching_factor: controls clustering intensity.
     delta > 0: cluster at both ends
     Larger delta = more uniform; smaller delta = more clustering.
   Based on: s(eta) = 1 + tanh(delta*(eta - 0.5)) / tanh(delta/2)
   ====================================================================== */
void Line::generate_tanh(Point &sp, Point &ep, int &n,
						 double stretching_factor)
{
	Point Temp_point;
	no_of_points = n;
	point_list.clear();
	point_list.reserve(n);

	ds_x = ep.Get_x() - sp.Get_x();
	ds_y = ep.Get_y() - sp.Get_y();
	ds_z = ep.Get_z() - sp.Get_z();
	Length = sqrt(ds_x * ds_x + ds_y * ds_y + ds_z * ds_z);

	vector<double> s_dist(n, 0.0);
	double delta = stretching_factor;
	double tanh_half = tanh(delta * 0.5);

	for (int j = 0; j < n; j++)
	{
		double eta = static_cast<double>(j) / (n - 1);
		double s_norm = 0.5 * (1.0 + tanh(delta * (eta - 0.5)) / tanh_half);
		s_dist[j] = s_norm * Length;
	}

	distribute_points(sp, s_dist);

	cout << "Created Line (tanh): " << point_list.size() << " points\n"
		 << "  Stretching factor: " << delta << "\n"
		 << "  First spacing: " << s_dist[1] << "\n"
		 << "  Mid spacing:   " << s_dist[n / 2 + 1] - s_dist[n / 2] << "\n"
		 << "  Last spacing:  " << s_dist[n - 1] - s_dist[n - 2] << endl;
}

/* ======================================================================
   Helper: distribute points along the line direction given an
   arc-length distribution vector s_dist[0..n-1] where s_dist[0]=0
   and s_dist[n-1]=Length.
   ====================================================================== */
void Line::distribute_points(Point &sp, const vector<double> &s_distribution)
{
	Point Temp_point;
	int n = s_distribution.size();
	double ux = ds_x / Length;
	double uy = ds_y / Length;
	double uz = ds_z / Length;

	for (int j = 0; j < n; j++)
	{
		double s = s_distribution[j];
		Temp_point(sp.Get_x() + s * ux,
				   sp.Get_y() + s * uy,
				   sp.Get_z() + s * uz);
		point_list.push_back(Temp_point);
	}
}

void Line::Reverse_Points()
{
	int size = point_list.size();
	vector<Point> Reverselist;
	Reverselist.reserve(size);
	for (int i = size - 1; i >= 0; i--)
		Reverselist.push_back(point_list[i]);
	point_list = Reverselist;
}

void Line::Clear()
{
	point_list.clear();
	no_of_points = 0;
	Length = 0.0;
	ds_x = ds_y = ds_z = 0.0;
}

void Line::Merge(Line &Line1)
{
	cout << "Merging two lines\n";
	vector<Point> Temp_Point_List;
	Temp_Point_List = Line1.Get_Point_list();
	cout << point_list.size() << "\t" << Temp_Point_List.size() << endl;

	if (no_of_points > 0 && point_list[no_of_points - 1] == Temp_Point_List[0])
	{
		for (unsigned int i = 1; i < Temp_Point_List.size(); i++)
			point_list.push_back(Temp_Point_List[i]);
	}
	else
	{
		for (unsigned int i = 0; i < Temp_Point_List.size(); i++)
			point_list.push_back(Temp_Point_List[i]);
	}
	no_of_points = point_list.size();
}

void Line::Merge(vector<Point> &List1)
{
	cout << "Merging list of points to current line\n";
	cout << point_list.size() << "\t" << List1.size() << endl;

	if (no_of_points > 0 && point_list[no_of_points - 1] == List1[0])
	{
		for (unsigned int i = 1; i < List1.size(); i++)
			point_list.push_back(List1[i]);
	}
	else
	{
		for (unsigned int i = 0; i < List1.size(); i++)
			point_list.push_back(List1[i]);
	}
	no_of_points = point_list.size();
}

void Line::write_output()
{
	Point temp_Point;
	ofstream myfileout("line.txt", ios::app | ios::binary);
	if (myfileout.is_open())
	{
		for (unsigned int i = 0; i < point_list.size(); i++)
		{
			temp_Point = get_Point(i);
			myfileout << temp_Point.Get_x() << "\t"
					  << temp_Point.Get_y() << "\t"
					  << temp_Point.Get_z() << "\t"
					  << 0.0 << endl;
		}
	}
	else
	{
		cerr << "Unable to open output file line.txt\n";
	}
	myfileout.close();
}

const vector<Point> &Line::Get_Point_list() const
{
	return point_list;
}

void Line::Print()
{
	for (unsigned int i = 0; i < point_list.size(); i++)
		point_list[i].Print();
	cout << "---------------------------------\n";
}

double Line::Size()
{
	return point_list.size();
}
