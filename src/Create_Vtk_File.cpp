#include "definitions.h"
#include "Globals.h"
#include "Grid.h"
#include "IO_Write.h"
#include "MPI_Utils.h"

void Read_Write_Grid(const string &Grid_File, const string &Out_Put_File)
{
	if (!CFD_MPI_Is_Root())
		return;

	double x, y, z;
	long int a, b, c, d, e;
	string s1, s2;
	int nop, No_of_Cells;
	ifstream Ip_Grid_File(Grid_File.c_str(), ios::in);
	ofstream Op_Sol_File(Out_Put_File.c_str(), ios::out);
	//		cout<<Grid_File<<endl;
	//		cout<<Out_Put_File<<endl;
	if (Ip_Grid_File.is_open())
	{
		getline(Ip_Grid_File, s1);
		Op_Sol_File << s1 << endl;
		getline(Ip_Grid_File, s1);
		Op_Sol_File << s1 << endl;
		getline(Ip_Grid_File, s1);
		Op_Sol_File << s1 << endl;
		getline(Ip_Grid_File, s1);
		Op_Sol_File << s1 << endl;
		Ip_Grid_File >> s1 >> nop >> s2;
		Op_Sol_File << s1 << "\t" << nop << "\t" << s2 << endl;
		for (int i = 0; i < nop; i++)
		{
			Ip_Grid_File >> x >> y >> z;
			Op_Sol_File << x << "\t" << y << "\t" << z << endl;
		}
		Ip_Grid_File >> s1 >> No_of_Cells >> nop;
		Op_Sol_File << s1 << "\t" << No_of_Cells << "\t" << nop << endl;
		for (int i = 0; i < No_of_Cells; i++)
		{
			Ip_Grid_File >> a >> b >> c >> d >> e;
			Op_Sol_File << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << endl;
		}
		Ip_Grid_File >> s1 >> nop;
		Op_Sol_File << s1 << "\t" << nop << endl;
		int p;
		for (int i = 0; i < nop; i++)
		{
			Ip_Grid_File >> p;
			Op_Sol_File << p << endl;
		}
		//			cout<<"Sucessfully Created File with Grid Data"<<endl;
	}
	else
	{
		cout << "Could Not Open Input Grid File......Please Check the File Name\n";
	}
}

void Append_Solution(const string &Sol_File, const string &Update_Solution)
{
	if (!CFD_MPI_Is_Root())
		return;

	double P, T, Rho, u, v, M, dt, Po;
	int nx_c, ny_c, nz_c, Cells_In_Plane, No_of_Cells, iterations, Cindex;
	string text1;
	V_D Pressure, Temperature, Density, U_Velocity, V_Velocity, W_Velocity, Mach_No, Pt;
	V_I Cell_IDs;
	V_I Cartesian_Cells;
	ifstream solution_ipfile(Sol_File.c_str(), ios::in);
	ifstream grid_ipfile(Update_Solution.c_str(), ios::in);
	ofstream Solution_Update(Update_Solution.c_str(), ios::out | ios::app);
	int grid_cell_count = -1;

	// Read cell count from the existing VTK grid header (line starts with "CELLS").
	if (grid_ipfile.is_open())
	{
		string tok;
		while (grid_ipfile >> tok)
		{
			if (tok == "CELLS")
			{
				grid_ipfile >> grid_cell_count;
				break;
			}
		}
		grid_ipfile.close();
	}

	//	cout<<Sol_File<<endl;
	//	cout<<Update_Solution<<endl;

	if (solution_ipfile.is_open())
	{
		//		cout<<"Solution File Opened\t"<<Sol_File<<endl;
		if (Grid_Type == 0)
		{
			solution_ipfile >> nx_c >> ny_c >> nz_c;
		}
		else if (Grid_Type == 1)
		{
			solution_ipfile >> nx_1 >> ny_1 >> nx_2 >> ny_2;
		}
		else
			solution_ipfile >> nx_c >> ny_c >> nz_c;

		solution_ipfile >> No_of_Cells;
		//		cout<<nx_c<<"\t"<<ny_c<<"\t"<<nz_c<<endl;
		//		cout<<No_of_Cells<<endl;
		solution_ipfile >> iterations;
		//		cout<<"Solution after number of iterations\t"<<iterations<<endl;
		for (int i = 0; i < No_of_Cells; i++)
		{
			solution_ipfile >> Cindex >> dt >> Rho >> P >> T >> u >> v >> M >> Po;
			//		    cout<<dt<<"\t"<<Cindex<<endl;
			Cell_IDs.push_back(Cindex);
			Pressure.push_back(P);
			Density.push_back(Rho);
			Temperature.push_back(T);
			U_Velocity.push_back(u);
			V_Velocity.push_back(v);
			Mach_No.push_back(M);
			Pt.push_back(Po);
		}

		if (Enable_AMR)
		{
			auto midpoint3 = [](const V_D &a, const V_D &b) -> V_D
			{
				V_D m(3, 0.0);
				m[0] = 0.5 * (a[0] + b[0]);
				m[1] = 0.5 * (a[1] + b[1]);
				m[2] = 0.5 * (a[2] + b[2]);
				return m;
			};

			auto centroid4 = [](const V_D &a, const V_D &b, const V_D &c, const V_D &d) -> V_D
			{
				V_D q(3, 0.0);
				q[0] = 0.25 * (a[0] + b[0] + c[0] + d[0]);
				q[1] = 0.25 * (a[1] + b[1] + c[1] + d[1]);
				q[2] = 0.25 * (a[2] + b[2] + c[2] + d[2]);
				return q;
			};

			function<bool(int, V_D &)> resolve_quad_vertices;
			resolve_quad_vertices = [&](int cid, V_D &verts) -> bool
			{
				if (cid < 0 || cid >= static_cast<int>(Cells.size()))
					return false;
				const Cell &c = Cells[cid];
				if (c.Cell_Vertices.size() == 12 && c.AMR_Parent < 0)
				{
					verts = c.Cell_Vertices;
					return true;
				}

				if (c.AMR_Parent >= 0 && c.AMR_Parent < static_cast<int>(Cells.size()))
				{
					V_D pv;
					if (!resolve_quad_vertices(c.AMR_Parent, pv) || pv.size() != 12)
						return false;

					V_D v0 = {pv[0], pv[1], pv[2]};
					V_D v1 = {pv[3], pv[4], pv[5]};
					V_D v2 = {pv[6], pv[7], pv[8]};
					V_D v3 = {pv[9], pv[10], pv[11]};
					V_D m01 = midpoint3(v0, v1);
					V_D m12 = midpoint3(v1, v2);
					V_D m23 = midpoint3(v2, v3);
					V_D m30 = midpoint3(v3, v0);
					V_D pc = centroid4(v0, v1, v2, v3);

					V_D quads[4];
					// 0: SW, 1: SE, 2: NE, 3: NW
					quads[0] = {v0[0], v0[1], v0[2], m01[0], m01[1], m01[2], pc[0], pc[1], pc[2], m30[0], m30[1], m30[2]};
					quads[1] = {m01[0], m01[1], m01[2], v1[0], v1[1], v1[2], m12[0], m12[1], m12[2], pc[0], pc[1], pc[2]};
					quads[2] = {pc[0], pc[1], pc[2], m12[0], m12[1], m12[2], v2[0], v2[1], v2[2], m23[0], m23[1], m23[2]};
					quads[3] = {m30[0], m30[1], m30[2], pc[0], pc[1], pc[2], m23[0], m23[1], m23[2], v3[0], v3[1], v3[2]};

					double best_d2 = 1e300;
					int best_q = -1;
					const double cx = (c.Cell_Center.size() > 0) ? c.Cell_Center[0] : 0.0;
					const double cy = (c.Cell_Center.size() > 1) ? c.Cell_Center[1] : 0.0;
					for (int q = 0; q < 4; q++)
					{
						V_D qc = centroid4(
							V_D{quads[q][0], quads[q][1], quads[q][2]},
							V_D{quads[q][3], quads[q][4], quads[q][5]},
							V_D{quads[q][6], quads[q][7], quads[q][8]},
							V_D{quads[q][9], quads[q][10], quads[q][11]});
						const double dx = cx - qc[0];
						const double dy = cy - qc[1];
						const double d2 = dx * dx + dy * dy;
						if (d2 < best_d2)
						{
							best_d2 = d2;
							best_q = q;
						}
					}
					if (best_q < 0)
						return false;
					verts = quads[best_q];
					return true;
				}

				if (c.Cell_Vertices.size() == 12)
				{
					verts = c.Cell_Vertices;
					return true;
				}
				return false;
			};

			unordered_map<int, int> solIndexByCellId;
			solIndexByCellId.reserve(Cell_IDs.size() * 2 + 1);
			for (int i = 0; i < static_cast<int>(Cell_IDs.size()); i++)
				solIndexByCellId[Cell_IDs[i]] = i;

			V_I leafCells;
			Build_Leaf_Cell_List(leafCells);

			vector<V_D> leafVerts;
			vector<int> leafIds;
			leafVerts.reserve(leafCells.size());
			leafIds.reserve(leafCells.size());
			for (const int cid : leafCells)
			{
				if (cid < 0 || cid >= No_Physical_Cells)
					continue;
				V_D verts;
				if (!resolve_quad_vertices(cid, verts))
					continue;
				if (verts.size() != 12)
					continue;
				leafIds.push_back(cid);
				leafVerts.push_back(verts);
			}

			ofstream vtk_out(Update_Solution.c_str(), ios::out);
			if (!vtk_out.is_open())
			{
				cout << "Could not Open Final Data file for Updating Solution\n";
				return;
			}

			const int nCellsOut = static_cast<int>(leafVerts.size());
			const int nPtsOut = 4 * nCellsOut;
			vtk_out << "# vtk DataFile Version 3.0" << endl;
			vtk_out << "AMR leaf solution output" << endl;
			vtk_out << "ASCII" << endl;
			vtk_out << "DATASET UNSTRUCTURED_GRID" << endl;
			vtk_out << "POINTS\t" << nPtsOut << "\tfloat" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				const V_D &vtx = leafVerts[i];
				vtk_out << vtx[0] << "\t" << vtx[1] << "\t" << vtx[2] << endl;
				vtk_out << vtx[3] << "\t" << vtx[4] << "\t" << vtx[5] << endl;
				vtk_out << vtx[6] << "\t" << vtx[7] << "\t" << vtx[8] << endl;
				vtk_out << vtx[9] << "\t" << vtx[10] << "\t" << vtx[11] << endl;
			}
			vtk_out << "CELLS\t" << nCellsOut << "\t" << (5 * nCellsOut) << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				const int p = 4 * i;
				vtk_out << 4 << "\t" << p << "\t" << p + 1 << "\t" << p + 2 << "\t" << p + 3 << endl;
			}
			vtk_out << "CELL_TYPES\t" << nCellsOut << endl;
			for (int i = 0; i < nCellsOut; i++)
				vtk_out << 9 << endl; // VTK_QUAD

			vtk_out << "CELL_DATA\t" << nCellsOut << endl;
			vtk_out << "SCALARS\t Pressure\t double" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(Pressure.size())) ? Pressure[sidx] : 0.0) << endl;
			}

			vtk_out << "SCALARS\t Temperature\t double " << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(Temperature.size())) ? Temperature[sidx] : 0.0) << endl;
			}

			vtk_out << "SCALARS\t Density\t double " << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(Density.size())) ? Density[sidx] : 0.0) << endl;
			}

			vtk_out << "SCALARS\t U_Velocity\t double " << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(U_Velocity.size())) ? U_Velocity[sidx] : 0.0) << endl;
			}

			vtk_out << "SCALARS\t V_Velocity\t double" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(V_Velocity.size())) ? V_Velocity[sidx] : 0.0) << endl;
			}

			vtk_out << "SCALARS\t Total_Pressure\t double" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(Pt.size())) ? Pt[sidx] : 0.0) << endl;
			}

			vtk_out << "SCALARS\t Entropy\t double" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				const double pval = (sidx >= 0 && sidx < static_cast<int>(Pressure.size())) ? Pressure[sidx] : 1.0;
				const double rhoval = (sidx >= 0 && sidx < static_cast<int>(Density.size())) ? Density[sidx] : 1.0;
				vtk_out << log(pval / pow(max(rhoval, 1e-14), 1.4)) << endl;
			}

			vtk_out << "SCALARS\t Mach_Number\t double " << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				vtk_out << ((sidx >= 0 && sidx < static_cast<int>(Mach_No.size())) ? Mach_No[sidx] : 0.0) << endl;
			}

			vtk_out << "Vectors\t Velocity\t double " << endl;
			for (int i = 0; i < nCellsOut; i++)
			{
				auto it = solIndexByCellId.find(leafIds[i]);
				const int sidx = (it != solIndexByCellId.end()) ? it->second : leafIds[i];
				const double uo = (sidx >= 0 && sidx < static_cast<int>(U_Velocity.size())) ? U_Velocity[sidx] : 0.0;
				const double vo = (sidx >= 0 && sidx < static_cast<int>(V_Velocity.size())) ? V_Velocity[sidx] : 0.0;
				vtk_out << uo << "\t" << vo << "\t" << 0.0 << endl;
			}

			vtk_out << "SCALARS\t AMR_Level\t int" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
				vtk_out << ((leafIds[i] >= 0 && leafIds[i] < static_cast<int>(Cells.size())) ? Cells[leafIds[i]].AMR_Level : 0) << endl;

			vtk_out << "SCALARS\t AMR_Parent\t int" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
				vtk_out << ((leafIds[i] >= 0 && leafIds[i] < static_cast<int>(Cells.size())) ? Cells[leafIds[i]].AMR_Parent : -1) << endl;

			vtk_out << "SCALARS\t AMR_IsLeaf\t int" << endl;
			vtk_out << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < nCellsOut; i++)
				vtk_out << 1 << endl;
			vtk_out.close();
			return;
		}

		Cells_In_Plane = (nx_c - 1) * (ny_c - 1);
		int No_of_Cells_Out = No_of_Cells;
		if (grid_cell_count > 0 && grid_cell_count < No_of_Cells_Out)
			No_of_Cells_Out = grid_cell_count;
		//      cout<<"Cells in Plane\t"<<Cells_In_Plane<<endl;
		if (Solution_Update.is_open())
		{
			//			cout<<"updating solution file"<<endl;
			Solution_Update << "CELL_DATA\t" << No_of_Cells_Out << endl;
			Solution_Update << "SCALARS\t Pressure\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << Pressure[i] << endl;

			Solution_Update << "SCALARS\t Temperature\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << Temperature[i] << endl;

			Solution_Update << "SCALARS\t Density\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << Density[i] << endl;

			Solution_Update << "SCALARS\t U_Velocity\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << U_Velocity[i] << endl;

			Solution_Update << "SCALARS\t V_Velocity\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << V_Velocity[i] << endl;

			Solution_Update << "SCALARS\t Total_Pressure\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << Pt[i] << endl;

			Solution_Update << "SCALARS\t Entropy\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << log(Pressure[i] / pow(Density[i], 1.4)) << endl;

			Solution_Update << "SCALARS\t Mach_Number\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << Mach_No[i] << endl;

			Solution_Update << "Vectors\t Velocity\t double " << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
				Solution_Update << U_Velocity[i] << "\t" << V_Velocity[i] << "\t" << 0.0 << endl;

			// AMR metadata
			Solution_Update << "SCALARS\t AMR_Level\t int" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
			{
				const int cid = (i < static_cast<int>(Cell_IDs.size())) ? Cell_IDs[i] : i;
				int lvl = (cid >= 0 && cid < static_cast<int>(Cells.size())) ? Cells[cid].AMR_Level : 0;
				Solution_Update << lvl << endl;
			}

			Solution_Update << "SCALARS\t AMR_Parent\t int" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
			{
				const int cid = (i < static_cast<int>(Cell_IDs.size())) ? Cell_IDs[i] : i;
				int p = (cid >= 0 && cid < static_cast<int>(Cells.size())) ? Cells[cid].AMR_Parent : -1;
				Solution_Update << p << endl;
			}

			Solution_Update << "SCALARS\t AMR_IsLeaf\t int" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (int i = 0; i < No_of_Cells_Out; i++)
			{
				const int cid = (i < static_cast<int>(Cell_IDs.size())) ? Cell_IDs[i] : i;
				int leaf = (cid >= 0 && cid < static_cast<int>(Cells.size()) && Cells[cid].AMR_IsLeaf) ? 1 : 0;
				Solution_Update << leaf << endl;
			}
		}
		else
		{
			cout << "Could not Open Final Data file for Updating Solution\n";
		}
	}
	else
	{
		cout << "Could not Open Solution File\n";
	}
}
