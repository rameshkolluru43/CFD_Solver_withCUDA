#include<iostream>
#include<string>
#include<vector>
#include<fstream>

using namespace std;

int main()
{
	string gridfilename,solutionfilename,final_vtkfile,final_structured_vtk;
	string s1,s2;
	double x,y,z,v1,v2,v3,P,T,rho;
	long int nop,a,b,c,d,e,f,g,h,j,timestep,No_of_Cells;
	int nx,ny,nz;
 	vector<double> pressure,temperature,density,vel1,vel2,vel3;
	gridfilename_core ="../Grid_Files/Pipe_Clustered2_ld10_61_21_41.vtk";
	solutionfilename="../NS_2O_Results/Pipe_AUSM_0.5__RK4_61_21_41_111325.txt";
	gridfilename_polar ="../Grid_Files/Pipe_Clustered2_ld10_61_21_41.vtk";
	solutionfilename_polar="../NS_2O_Results/Pipe_AUSM_0.5__RK4_61_21_41_111325.txt";
	final_vtkfile_core="../NS_2O_Results/Pipe_AUSM_61_21_41_Solution.vtk";
	final_vtkfile_polar="../NS_2O_Results/Pipe_AUSM_61_21_41_Solution.vtk";
// 	final_structured_vtk="/home2/ramesh/Desktop/Ph.D-Work/Research_Code/NS_Results/Pipe_Trial1_fd_struct_5000_RK4_31_31_71_111325.vtk";
	gridfilename ="../Euler_Cube_Results/channel_11_11_11.vtk";
	ifstream gridfile(gridfilename.c_str(),ios::in);
/*	solutionfilename="../Euler_Cube_Results/Cube_AUSM_11_11_11.txt";
	final_vtkfile="../Euler_Cube_Results/Cube_AUSM_11_11_11_Solution.vtk";
	ifstream solfile(solutionfilename.c_str(),ios::in);*/
	ofstream Pipefile(final_vtkfile.c_str(),ios::out|ios::app);
// 	ofstream Pipefile_struct(final_structured_vtk,ios::out|ios::app);
	if(Pipefile.is_open())
	{
		cout<<"Sucessfully opened Grid file\t"<<gridfilename<<endl;
		if(gridfile.is_open())
		{
			getline(gridfile,s1);
			Pipefile<<s1<<endl;
			getline(gridfile,s1);
			Pipefile<<s1<<endl;
			getline(gridfile,s1);
			Pipefile<<s1<<endl;
			getline(gridfile,s1);
			Pipefile<<s1<<endl;
			gridfile>>s1>>nop>>s2;
			Pipefile<<s1<<"\t"<<nop<<"\t"<<s2<<endl;
			for(int i=0;i<nop;i++)
			{
				gridfile>>x>>y>>z;
				Pipefile<<x<<"\t"<<y<<"\t"<<z<<endl;
			}
			gridfile>>s1>>No_of_Cells>>nop;
			Pipefile<<s1<<"\t"<<No_of_Cells<<"\t"<<nop<<endl;
			for(int i=0;i<No_of_Cells;i++)
			{
				gridfile>>a>>b>>c>>d>>e>>f>>g>>h>>j;
// 				Pipefile<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<e<<"\t"<<f<<"\t"<<g<<"\t"<<h<<"\t"<<j<<endl;
				Pipefile<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<e<<"\t"<<f<<"\t"<<g<<"\t"<<h<<"\t"<<j<<endl;
			}
			gridfile>>s1>>nop;
			Pipefile<<s1<<"\t"<<nop<<endl;
			int p;
			for(int i=0;i<nop;i++)
			{
				gridfile>>p;
				Pipefile<<p<<endl;
			}
		}
		else
		{
			cout<<"Could not open Grid datafile, Please check the file name\n";
		}
		gridfile.close();
		cout<<"Closed grid file \t"<<gridfilename<<endl;
		if(solfile.is_open())
		{
			cout<<"Solution file opened\t"<<solutionfilename<<endl;
			solfile>>s1>>s2>>timestep;
			cout<<s1<<" "<<s2<<"\tSolution at time step\t"<<timestep<<endl;
			solfile>>nx>>ny>>nz;
			cout<<nx<<"\t"<<ny<<"\t"<<nz<<endl;
			cout<<No_of_Cells<<endl;
			for(int i=0;i<No_of_Cells;i++)
			{
				solfile>>a>>P>>T>>rho>>v1>>v2>>v3;
				cout<<P<<"\t"<<T<<"\t"<<v1<<"\t"<<v2<<"\t"<<v3<<endl;
				pressure.push_back(P);
				temperature.push_back(T);
				density.push_back(rho);
				vel1.push_back(v1);
				vel2.push_back(v2);
				vel3.push_back(v3);
			}
			solfile.close();
			cout<<"size of presure data\t"<<pressure.size()<<endl;
		}
		else
		{
			cout<<"Could not open solution file, Please check the file name\n";
		}
		Pipefile<<"CELL_DATA\t"<<pressure.size()<<endl;
		Pipefile<<"SCALARS\t Pressure\t double \t 1"<<endl;
		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<pressure.size();i++)
		{
			Pipefile<<pressure[i]<<endl;
		}
		Pipefile<<"SCALARS\t temperature\t double\t1"<<endl;
		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<temperature.size();i++)
		{
			Pipefile<<temperature[i]<<endl;
		}
		Pipefile<<"SCALARS\t Density\tdouble\t1"<<endl;
		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<density.size();i++)
		{
			Pipefile<<density[i]<<endl;
		}
		Pipefile<<"SCALARS\t u_velocity\tdouble\t1"<<endl;
		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<vel1.size();i++)
		{
			Pipefile<<vel1[i]<<endl;
		}
		Pipefile<<"SCALARS\t v_velocity\t double\t1"<<endl;
		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<vel2.size();i++)
		{
			Pipefile<<vel2[i]<<endl;
		}
		Pipefile<<"SCALARS\tw_velocity\t double\t1"<<endl;
		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<vel3.size();i++)
		{
			Pipefile<<vel3[i]<<endl;
		}
		Pipefile<<"VECTORS\t Velocity\t double\t1"<<endl;
// 		Pipefile<<"LOOKUP_TABLE default"<<endl;
		for(int i=0;i<vel1.size();i++)
		{
			Pipefile<<vel2[i]<<"\t"<<vel3[i]<<"\t"<<vel1[i]<<endl;
		}
		Pipefile.close();
	}
	else
	{
		cout<<"Unable to open file to write the Pipe solution\n";
	}
cout<<"Generated solution vtk file with grid\n";
return 0;
}