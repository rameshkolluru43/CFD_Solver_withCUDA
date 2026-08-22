#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<cmath>
#include "gnuplot_i.hpp"

#define sleep_length  1
using namespace std;

int main()
{
	string gridfilename,solutionfilename,final_vtkfile,final_structured_vtk;
	string s1,s2;
	double x,y,z,v1,v2,v3,P,T,rho;
	long int nop,a,b,c,d,e,f,g,h,j,timestep,No_of_Cells;
	int nx,ny,nz;
 	vector<double> pressure,temparature,density,vel1,vel2,vel3;
//	solutionfilename="../NS_2O_Results/Pipe_AUSM_0.5__RK4_11_21_11.txt";
// 	solutionfilename="../NS_4O_RE_Results/Pipe_2DOFF_RK4_31_71_51_111325.txt";
// 	solutionfilename = "../Euler_Cube_Results/Cube_AUSM_0.2_11_11_11.txt";
	solutionfilename = "/home/ramesh/NS_Code_Plain/Solution_Files/Data_Primitve_variables_Channel10h_mm_11_11_11.txt";
	Gnuplot g1= Gnuplot("linespoints");
	Gnuplot g2= Gnuplot("linespoints");
	Gnuplot g3= Gnuplot("linespoints");
	ifstream solfile(solutionfilename.c_str(),ios::in);
	if(solfile.is_open())
	{
		cout<<"Solution file opened\t"<<solutionfilename<<endl;
// 		solfile>>s1>>s2>>timestep;
// 		cout<<s1<<" "<<s2<<"\tSolution at time step\t"<<timestep<<endl;
		solfile>>nx>>ny>>nz;
		cout<<nx<<"\t"<<ny<<"\t"<<nz<<endl;
		No_of_Cells = (nx-1)*(ny-1)*(nz-1);
		for(int i=0;i<No_of_Cells;i++)
		{
				solfile>>a>>P>>T>>rho>>v1>>v2>>v3;
// 				cout<<P<<"\t"<<T<<"\t"<<v1<<"\t"<<v2<<"\t"<<v3<<endl;
				pressure.push_back(P);
				temparature.push_back(T);
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
	vector<double> Average_Velocity,x_data;
	double vel_mag,temp=0.0;
	int index=0;
	for(int k=0;k<nz-1;k++)
	{	
		for(int j=0;j<(ny-1);j++)
		{
			vel_mag=0.0;
			x_data.push_back(0.5*0.019*j*(1.0/(ny-2)));
			for(int i=0;i<nx-1;i++)
			{
				index = i+j*(nx-1)+k*(ny-1)*(nx-1);
 				vel_mag += sqrt(vel1.at(index)*vel1.at(index)+vel2.at(index)*vel2.at(index)+vel3.at(index)*vel3.at(index));
			}
 			vel_mag /= (nx-1);
			Average_Velocity.push_back(vel_mag);
		}
		g1.plot_xy(Average_Velocity,x_data,"Velocity_Profile");
		sleep(sleep_length);
		Average_Velocity.clear();
		x_data.clear();
	}
// 	g1.reset_plot();
	vector<double> data,data1,data2;
	int i=0;
	for(int k=0;k<nz-1;k++)
	{
		for(int j=0;j<ny-1;j++)
		{
			index = i+j*(nx-1)+k*(ny-1)*(nx-1);
			data.push_back(density[index]);
			data1.push_back(pressure[index]);
			data2.push_back(temparature[index]);
		}
		g1.plot_x(data,"Density");
		g2.plot_x(data1,"Pressure");
		g3.plot_x(data2,"Temparature");
		sleep(sleep_length);
// 		data.clear();
// 		data1.clear();
// 		data2.clear();
// 		g1.reset_plot();
// 		g2.reset_plot();
// 		g3.reset_plot();
	}
return 0;
}