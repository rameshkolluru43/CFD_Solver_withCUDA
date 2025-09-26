#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"

extern "C" void Ellinga_Carbuncle()
{
    Directory_Name();
    File_Name();
	
	switch(Grid_Size)
	{
		case 1:
			Grid_File = "../Grid_Files/EllingCase/EllingCase_101_41.txt";
			Grid_Vtk_File = "../Grid_Files/EllingCase/EllingCase_101_41.vtk";
	        Error_File +="_101_41.txt";
	        Initial_Solution_File +="_101_41.txt";
	        Solution_File +="_101_41.txt";  
			Final_Solution_File +="_101_41.vtk"; 
			break;
			case 2:
			Grid_File = "../Grid_Files/EllingCase/EllingCase_201_81.txt";
			Grid_Vtk_File = "../Grid_Files/EllingCase/EllingCase_201_81.vtk";
	        Error_File +="_201_81.txt";
	        Initial_Solution_File +="_201_81.txt";
	        Solution_File +="_201_81.txt";  
			Final_Solution_File +="_201_81.vtk"; 
			
			break;
			case 3:
			Grid_File = "../Grid_Files/EllingCase/EllingCase_1001_401.txt";
			Grid_Vtk_File = "../Grid_Files/EllingCase/EllingCase_1001_401.vtk";
	        Error_File +="_1001_401.txt";
	        Initial_Solution_File +="_1001_401.txt";
	        Solution_File +="_1001_401.txt";  
			Final_Solution_File +="_1001_401.vtk"; 
			
			break;
	}
	
        
/* Reads Input grid file and does the preprocessing required for grid 
 * Calculates the Cell Normals, Cell Areas
 * Checks for Grid 
 */
	Form_Cells(Grid_File);
	cout<<"Grid_Type used \t"<<Grid_Type<<endl;
//	Write_Cell_Info("../Grid_Files/Forward_Step_Grid_Files/Forward_Step_241_81.txt");
        V_D V(2,0.0); 
        double Rho,P,M,u,C;
        cout<<" Initial Data for Solution\t"<<endl;
	cout<<"Initialize from a file or from Zero, Enter 1 to read data from file:\t";
	
	Is_Time_Dependent = true;
	Terminating_Time=100.0;
	cin>>Initialize_Type;
	int index=0;
	if(Initialize_Type==1)
		Initialize(Initial_Solution_File);	
	else
        {
            Initialize(Test_Case);
			cout<<"Initalizing Data with inlet conditions \t"<<endl;
			for(int j=0;j<ny_c;j++)
			{
            	for(int i=0;i<nx_c;i++)
            	{
					
					index = i+j*(nx_c-1);
//					cout<<i<<"\t"<<j<<"\t"<<index<<endl;
					if(i <=0.5*nx_c)
					{
						M = 3.0;
				    	Rho = 1.0;
						u   = 1.0;
						P_L   = 1.0/(gamma*M*M);
						P=P_L;
						C = sqrt(gamma*P_L/Rho);
						if(j == 0.5*(ny_c-1))
						{
//							 cout<<j<<"\t"<<i<<"\t"<<index<<endl;
							V_1 = 0.0;
							V_2 = 0.0;    
						}
						else
						{
							V_1 = M*C;
							V_2 = 0.0;    
						}				
					}
					else 
					{
						P_R   = P_L*(2.0*gamma*M*M-(gamma-1.0))/(gamma+1.0);
						P = P_R;
						Rho = ((gamma+1.0)/(gamma-1.0)*(P_R/P_L)+1.0)/((gamma+1.0)/(gamma-1.0)+(P_R/P_L));
						u   =sqrt(gamma*(((2.0+(gamma-1.0)*M*M)*P_R)/((2.0*gamma*M*M+(1.0-gamma))*Rho)));
						C = sqrt(gamma*P_R/Rho);
						V_1 = u;
						V_2 = 0.0;    
						
					}	
					V[0] = V_1;
					V[1] = V_2;
					Calculate_Computational_Variables(P,V,Rho,2);
					for(unsigned int jj=0;jj<Global_U.size();jj++)
						U_Cells[index][jj] = Global_U[jj];
					Calculate_Primitive_Variables(index,U_Cells[index]);
					for(unsigned int jj=0;jj<Global_Primitive.size();jj++)
                    	Primitive_Cells[index][jj] = Global_Primitive[jj];
            	}
        	}
		}
 
			Identify_Wall_Boundary_Faces(Grid_Type);   
			cout<<"Writing Initial Solution to file "<<Initial_Solution_File<<endl;
			Write_Solution(Initial_Solution_File,1);
//			Read_Write_Grid(Grid_Vtk_File,Final_Solution_File);
//			Append_Solution(Initial_Solution_File,Final_Solution_File);
			
			cout<<"Intialized Solution with inlet conditions, Identified Boundaries......... Ready to solve"<<endl;
//			exit(0);
}
