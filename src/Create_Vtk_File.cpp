#include "definitions.h"
#include "Globals.h"
#include "IO_Write.h"

void Read_Write_Grid(const string &Grid_File, const string &Out_Put_File)
{
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
	double P, T, Rho, u, v, M, dt, Po;
	int nx_c, ny_c, nz_c, Cells_In_Plane, No_of_Cells, iterations, Cindex;
	string text1;
	V_D Pressure, Temperature, Density, U_Velocity, V_Velocity, W_Velocity, Mach_No, Pt;
	V_I Cartesian_Cells;
	ifstream solution_ipfile(Sol_File.c_str(), ios::in);
	ofstream Solution_Update(Update_Solution.c_str(), ios::out | ios::app);

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
			Pressure.push_back(P);
			Density.push_back(Rho);
			Temperature.push_back(T);
			U_Velocity.push_back(u);
			V_Velocity.push_back(v);
			Mach_No.push_back(M);
			Pt.push_back(Po);
		}
		Cells_In_Plane = (nx_c - 1) * (ny_c - 1);
		//      cout<<"Cells in Plane\t"<<Cells_In_Plane<<endl;
		if (Solution_Update.is_open())
		{
			//			cout<<"updating solution file"<<endl;
			Solution_Update << "CELL_DATA\t" << No_of_Cells << endl;
			Solution_Update << "SCALARS\t Pressure\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < Pressure.size(); i++)
				Solution_Update << Pressure[i] << endl;

			Solution_Update << "SCALARS\t Temperature\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < Temperature.size(); i++)
				Solution_Update << Temperature[i] << endl;

			Solution_Update << "SCALARS\t Density\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < Density.size(); i++)
				Solution_Update << Density[i] << endl;

			Solution_Update << "SCALARS\t U_Velocity\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < U_Velocity.size(); i++)
				Solution_Update << U_Velocity[i] << endl;

			Solution_Update << "SCALARS\t V_Velocity\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < V_Velocity.size(); i++)
				Solution_Update << V_Velocity[i] << endl;

			Solution_Update << "SCALARS\t Total_Pressure\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < Pt.size(); i++)
				Solution_Update << Pt[i] << endl;

			Solution_Update << "SCALARS\t Entropy\t double" << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < Pressure.size(); i++)
				Solution_Update << log(Pressure[i] / pow(Density[i], 1.4)) << endl;

			Solution_Update << "SCALARS\t Mach_Number\t double " << endl;
			Solution_Update << "LOOKUP_TABLE default" << endl;
			for (unsigned int i = 0; i < Mach_No.size(); i++)
				Solution_Update << Mach_No[i] << endl;

			Solution_Update << "Vectors\t Velocity\t double " << endl;
			for (int i = 0; i < No_of_Cells; i++)
				Solution_Update << U_Velocity[i] << "\t" << V_Velocity[i] << "\t" << 0.0 << endl;
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
