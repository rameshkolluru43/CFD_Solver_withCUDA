#include "../Basic_Function_Files/Geometry_Header.h"

vector<Point> Line_List,Grid_list,grid_cube_list;
//Inlet_Boundary_List -- list of Inlet boundary points

vector<int> Cell_Point_Data,Cell_Neighbour_Data,Inlet_Boundary_List,Wall_Boundary_List,Exit_Boundary_List,Symmetry_Boundary_List;
vector<int> Left_Boundary_Ghost_Cells, Right_Boundary_Ghost_Cells, Bottom_Boundary_Ghost_Cells, Top_Boundary_Ghost_Cells, Physical_Cells;


void Append_Boundary_Information(const string & Inputfile)
{
	ofstream File_Out(Inputfile.c_str(),ios::app);
	if(File_Out.is_open())
	{
		File_Out<<"Exit_Boundary_Cells\t"<<1<<"\t"<<Exit_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Exit_Boundary_List.size();i+=3)
			File_Out<<Exit_Boundary_List[i+0]<<"\t"<<Exit_Boundary_List[i+1]<<"\t"<<Exit_Boundary_List[i+2]<<endl;
		File_Out<<"Inlet_Boundary_Cells\t"<<0<<"\t"<<Inlet_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Inlet_Boundary_List.size();i+=3)
			File_Out<<Inlet_Boundary_List[i+0]<<"\t"<<Inlet_Boundary_List[i+1]<<"\t"<<Inlet_Boundary_List[i+2]<<endl;
		File_Out<<"Wall_Boundary_Cells\t"<<2<<"\t"<<Wall_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Wall_Boundary_List.size();i+=3)
			File_Out<<Wall_Boundary_List[i+0]<<"\t"<<Wall_Boundary_List[i+1]<<"\t"<<Wall_Boundary_List[i+2]<<endl;
		File_Out<<"Symmetry_Boundary_Cells\t"<<3<<"\t"<<Symmetry_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Symmetry_Boundary_List.size();i+=3)
			File_Out<<Symmetry_Boundary_List[i+0]<<"\t"<<Symmetry_Boundary_List[i+1]<<"\t"<<Symmetry_Boundary_List[i+2]<<endl;
	}
	else
	{
		cout<<"Could not open file to append Boundary points list\n";
	}
}

void Identify_Cells_and_Assign_Neighbours(int & Nx, int & Ny)
{
    // Number of points, physical cells, ghost cells, and total cells
    int Number_Of_Points = Nx * Ny;
    int Number_Physical_Cells = (Nx - 1) * (Ny - 1);
    int Number_Ghost_Cells = 2 * (Nx - 1) + 2 * (Ny - 1);
    int Total_Cells = Number_Physical_Cells + Number_Ghost_Cells;
    int a;
    cout << "Number of Points\t" << Number_Of_Points << endl;
    cout << "Number of Physical Cells\t" << Number_Physical_Cells << endl;
    cout << "Number of Ghost Cells\t" << Number_Ghost_Cells << endl;


    // Resize the vectors
    Bottom_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Left_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Top_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Right_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Physical_Cells.resize(Number_Physical_Cells, 0);
    Cell_Neighbour_Data.resize(Number_Physical_Cells*5,0);

    // Bottom boundary ghost cells
    for (int i = 0; i < (Nx - 1); i++)
    {
        Bottom_Boundary_Ghost_Cells[i] = i;
    }

    // Left boundary ghost cells
    for (int i = 0; i < (Ny - 1); i++)
    {
        Left_Boundary_Ghost_Cells[i] = (Nx - 1) + i * (Nx + 1);
    }

    // Physical cells
    for (int j = 0; j < (Ny - 1); j++)
    {
        for (int i = 0; i < (Nx - 1); i++)
        {
            Physical_Cells[i + j * (Nx - 1)] = i + j * (Nx - 1 + 2) + Nx;
        }
    }

    // Right boundary ghost cells
    for (int i = 0; i < (Ny - 1); i++)
    {
        Right_Boundary_Ghost_Cells[i] = 2 * (Nx - 1) + 1 + i * (Nx + 1);
    }

    // Top boundary ghost cells
    for (int i = 0; i < (Nx - 1); i++)
    {
        Top_Boundary_Ghost_Cells[i] = Number_Physical_Cells + (Nx - 1) + 2 * (Ny - 1) + i;
    }

    int index0, index1, index2, index3, index4;
    for (int j = 0; j < (Ny - 1); j++)
    {
        for (int i = 0; i < (Nx - 1); i++)
        {
            index0 = Physical_Cells[i + j * (Nx - 1)]; // Current cell (i,j,k)
            index1 = (i == 0) ? Left_Boundary_Ghost_Cells[j] : Physical_Cells[(i - 1) + j * (Nx - 1)]; // Left (i-1,j,k)
            index2 = (j == 0) ? Bottom_Boundary_Ghost_Cells[i] : Physical_Cells[i + (j - 1) * (Nx - 1)]; // Bottom (i,j-1,k)
            index3 = (i == (Nx - 2)) ? Right_Boundary_Ghost_Cells[j] : Physical_Cells[(i + 1) + j * (Nx - 1)]; //Right (i+1,j,k)
	    index4 = (j == (Ny - 2)) ? Top_Boundary_Ghost_Cells[i] : Physical_Cells[i + (j + 1) * (Nx - 1)]; // Top (i,j+1,k)
//	    cout<<"Am here"<<endl;
//	    cout <<i + j * (Nx - 1)<<"\t" << index0 << "\t" << index1 << "\t" << index2 << "\t" << index3 << "\t" << index4 << endl;
	  
	    a = (i + j * (Nx - 1))*5;
	  
//	    cout<<i + j * (Nx - 1)<<"\tValue of a \t"<<a<<endl;
	    
	    Cell_Neighbour_Data[a+0]=index0;	//Current Cell	     (i,j,k)
	    Cell_Neighbour_Data[a+1]=index1;	//Left Cell	     (i-1,j,k)
	    Cell_Neighbour_Data[a+2]=index2;	//Bottom Cell	     (i,j-1,k)
	    Cell_Neighbour_Data[a+3]=index3;	//Right Cell	     (i+1,j,k)
	    Cell_Neighbour_Data[a+4]=index4;	//TopCell	     (i,j+1,k)	
	    if(i == 0)
	    {
		    Exit_Boundary_List.push_back(index0);
		    Exit_Boundary_List.push_back(0);
		    Exit_Boundary_List.push_back(Left_Boundary_Ghost_Cells[j]);
	    }
	    if(i == (Nx - 2))
	    {
		    Exit_Boundary_List.push_back(index0);
		    Exit_Boundary_List.push_back(2);
		    Exit_Boundary_List.push_back(Right_Boundary_Ghost_Cells[j]);
	    }
	    if(j == 0)
	    {
		    Wall_Boundary_List.push_back(index0);
		    Wall_Boundary_List.push_back(1);
		    Wall_Boundary_List.push_back(Bottom_Boundary_Ghost_Cells[i]);
	    }  
	    if(j == (Ny - 2))
	    {
		    Inlet_Boundary_List.push_back(index0);
		    Inlet_Boundary_List.push_back(3);
		    Inlet_Boundary_List.push_back(Top_Boundary_Ghost_Cells[i]);
	    }	    
        }
    }
/*    cout<<"In The Identifying function"<<endl;
    for(int i=0;i<Number_Physical_Cells;i++)
    {
	    int k = i*5;
	    cout<<k<<"\t"<<Cell_Neighbour_Data[k+0]<<"\t"<<Cell_Neighbour_Data[k+1]<<"\t"<<
				Cell_Neighbour_Data[k+2]<<"\t"<<Cell_Neighbour_Data[k+3]<<"\t"<<
				Cell_Neighbour_Data[k+4]<<"\n";			
    	
    }*/
}

void Identify_Cells_and_Assign_Neighbours_Custom(int & Nx, int & Ny)
{
    // Number of cells
    int Number_Physical_Cells = (Nx - 1) * (Ny - 1);
    int Number_Ghost_Cells = 2 * (Nx - 1) + 2 * (Ny - 1);
    int Total_Cells = Number_Physical_Cells + Number_Ghost_Cells;

    // Data structures
    Bottom_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Left_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Top_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Right_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Physical_Cells.resize(Number_Physical_Cells, 0);
    Cell_Neighbour_Data.resize(Number_Physical_Cells * 5, 0);

    int currentIndex = 0;

    // Assign Bottom Ghost Cells
    for (int i = 0; i < (Nx - 1); i++)
    {
        Bottom_Boundary_Ghost_Cells[i] = currentIndex++;
    }

    // Assign Rows: Left Ghost Cell → Physical Cells → Right Ghost Cell
    for (int j = 0; j < (Ny - 1); j++)
    {
        Left_Boundary_Ghost_Cells[j] = currentIndex++; // Left ghost cell for the row

        for (int i = 0; i < (Nx - 1); i++)
        {
            Physical_Cells[i + j * (Nx - 1)] = currentIndex++; // Physical cells in the row
        }

        Right_Boundary_Ghost_Cells[j] = currentIndex++; // Right ghost cell for the row
    }

    // Assign Top Ghost Cells
    for (int i = 0; i < (Nx - 1); i++)
    {
        Top_Boundary_Ghost_Cells[i] = currentIndex++;
    }

    // Assign neighbors and boundary conditions
    int index0, index1, index2, index3, index4;
    for (int j = 0; j < (Ny - 1); j++)
    {
        for (int i = 0; i < (Nx - 1); i++)
        {
            index0 = Physical_Cells[i + j * (Nx - 1)]; // Current cell
            index1 = (i == 0) ? Left_Boundary_Ghost_Cells[j] : Physical_Cells[(i - 1) + j * (Nx - 1)]; // Left
            index2 = (j == 0) ? Bottom_Boundary_Ghost_Cells[i] : Physical_Cells[i + (j - 1) * (Nx - 1)]; // Bottom
            index3 = (i == (Nx - 2)) ? Right_Boundary_Ghost_Cells[j] : Physical_Cells[(i + 1) + j * (Nx - 1)]; // Right
            index4 = (j == (Ny - 2)) ? Top_Boundary_Ghost_Cells[i] : Physical_Cells[i + (j + 1) * (Nx - 1)]; // Top

            int a = (i + j * (Nx - 1)) * 5;

            // Assign neighbors
            Cell_Neighbour_Data[a + 0] = index0;
            Cell_Neighbour_Data[a + 1] = index1;
            Cell_Neighbour_Data[a + 2] = index2;
            Cell_Neighbour_Data[a + 3] = index3;
            Cell_Neighbour_Data[a + 4] = index4;

            // Boundary condition lists
            if (i == 0)
            {
                Inlet_Boundary_List.push_back(index0);
                Inlet_Boundary_List.push_back(0); // Face ID
                Inlet_Boundary_List.push_back(Left_Boundary_Ghost_Cells[j]);
            }
            if (i == (Nx - 2))
            {
                Exit_Boundary_List.push_back(index0);
                Exit_Boundary_List.push_back(2); // Face ID
                Exit_Boundary_List.push_back(Right_Boundary_Ghost_Cells[j]);
            }
	    if (j == 0)
	    {
	        if (i < 60 || i > 120)
	        {
	            Symmetry_Boundary_List.push_back(index0);
	            Symmetry_Boundary_List.push_back(1); // Face ID
	            Symmetry_Boundary_List.push_back(Bottom_Boundary_Ghost_Cells[i]);
	        }
	        else
	        {
	            Wall_Boundary_List.push_back(index0);
	            Wall_Boundary_List.push_back(1); // Face ID
	            Wall_Boundary_List.push_back(Bottom_Boundary_Ghost_Cells[i]);
	        }
	    }
            if (j == (Ny - 2))
            {
                Wall_Boundary_List.push_back(index0);
                Wall_Boundary_List.push_back(3); // Face ID
                Wall_Boundary_List.push_back(Top_Boundary_Ghost_Cells[i]);
            }
        }
    }
}
void Identify_Cells(int &nx,int &ny)
{
	
	int no_of_Cells,a,p=0,Points_In_Plane;
	no_of_Cells=(nx-1)*(ny-1);
	Points_In_Plane = nx*ny;
	cout<<no_of_Cells<<endl;
	Cell_Point_Data.resize(no_of_Cells*4,0);
	for(int j=0;j<(ny-1);j++)
	{
		for(int i=0;i<(nx-1);i++)
		{
			a=p*4;
			Cell_Point_Data[a+0]=i+j*nx;			// i,j,k
			Cell_Point_Data[a+1]=(i+1)+j*nx;		//i+1,j,k
			Cell_Point_Data[a+2]=(i+1)+(j+1)*nx;		//i+1,j+1,k
			Cell_Point_Data[a+3]=i+(j+1)*nx;		//i,j+1,k			
//			cout<<Cell_Point_Data[a+0]<<"\t"<<Cell_Point_Data[a+1]<<"\t"<<Cell_Point_Data[a+2];
//			cout<<"\t"<<Cell_Point_Data[a+3]<<endl ;
			p++;
		}
	}
	cout<<"Identifying Cells done \n";
	cout<<"Number of Cells\t"<<p<<endl;
}

void Identify_Cells(int &nx1,int &ny1,int & nx2,int & ny2)
{
	cout<<nx1<<"\t"<<ny1<<"\t"<<nx2<<"\t"<<ny2<<endl;
	int no_of_Cells,a,p=0,Points_In_Plane;
	no_of_Cells=(nx1-1)*(ny1-1) + (nx2-1)*(ny2-1);
	Points_In_Plane = nx1*ny1 + nx2*ny2;
	cout<<"Total Number of Points\t"<<Points_In_Plane<<endl;
	cout<<"Total Cells to be formed using Points in Plane \t"<<no_of_Cells<<endl;
	Cell_Point_Data.resize(no_of_Cells*4,0);
	for(int j=0;j<((ny1+ny2)-2);j++)
	{
		for(int i=0;i<(nx2-1);i++)
		{
			a=p*4;
			if(j < ny1-1) 
			{
				if(i < nx1-1)
				{	
//					cout<<"("<<i<<","<<j<<")\t";
					Cell_Point_Data[a+0]=i+j*nx1;			// i,j,k
					Cell_Point_Data[a+1]=(i+1)+j*nx1;		//i+1,j,k
					Cell_Point_Data[a+2]=(i+1)+(j+1)*nx1;	//i+1,j+1,k
					Cell_Point_Data[a+3]=i+(j+1)*nx1;		//i,j+1,k					
				}
				else{
					
					break;
				}
			}
			else
			{
//				cout<<"("<<i<<","<<j<<")\t";
				Cell_Point_Data[a+0] = i+(j-(ny1-1))*nx2 + nx1*ny1  ;			// i,j,k
				Cell_Point_Data[a+1] = (i+1)+(j-(ny1-1))*nx2 + nx1*ny1;		//i+1,j,k
				Cell_Point_Data[a+2] = (i+1)+((j+1)-(ny1-1))*nx2 + nx1*ny1;	//i+1,j+1,k
				Cell_Point_Data[a+3] = i+((j+1)-(ny1-1))*nx2 + nx1*ny1;		//i,j+1,k			
			}
//			cout<<Cell_Point_Data[a+0]<<"\t"<<Cell_Point_Data[a+1]<<"\t"<<Cell_Point_Data[a+2];
//			cout<<"\t"<<Cell_Point_Data[a+3]<<endl ;
			p++;
		}
//		cout<<"------------------------------------------------"<<endl;
	}
	cout<<"Identifying Cells done \n";
	cout<<"Number of Cells\t"<<p<<endl;
}

// This function identifies the neighbours for given cell
void Identify_Neighbours(int &nx,int &ny)
{
	
	int nbwgc,ntwgc,nlwgc,nrwgc,nwgc,nigc,nogc,ngc,a,index0,index1,index2,index3,index4;
	int Ghost_Cell_Top_Plane,Ghost_Cell_Bottom_Plane;
	int Ghost_Cells_Left_Plane,Ghost_Cells_Right_Plane,cells_in_Plane,no_of_Cells;
	Ghost_Cell_Top_Plane=0;Ghost_Cell_Bottom_Plane=0;
	Ghost_Cells_Left_Plane=0;Ghost_Cells_Right_Plane=0;
	cells_in_Plane=0;
	nwgc=0;nigc=0;nogc=0;ngc=0,a=0;
	nbwgc=0;nrwgc=0;ntwgc=0;nlwgc=0;
	no_of_Cells = (nx-1)*(ny-1);
	Cell_Neighbour_Data.resize(no_of_Cells*5,0);
	cells_in_Plane = (nx-1)*(ny-1);
	Ghost_Cell_Bottom_Plane=(nx-1);
	Ghost_Cell_Top_Plane=(nx-1);
	Ghost_Cells_Left_Plane=(ny-1);
	Ghost_Cells_Right_Plane=(ny-1);
// 	cout<<"nx\t"<<nx<<"\tny\t"<<ny<<"\tNumber of Cells \t"<<no_of_Cells<<endl;
	
	for(int j=0;j<(ny-1);j++)
	{
	  for(int i=0;i<(nx-1);i++)
	  {
// 		cout<<i<<"\t"<<j<<"\t"<<k<<endl;
		index0 = i+j*(nx-1);// current cell		(i,j,k)
		  if(i==(nx-2))
		    {
			  	index3 = no_of_Cells+Ghost_Cell_Bottom_Plane+nlwgc++;
			  	Exit_Boundary_List.push_back(index0);
			  	Exit_Boundary_List.push_back(2);
			 	Exit_Boundary_List.push_back(index3);
		    }
		else
		  index3 = (i+1)+j*(nx-1);	//Right cell	(i+1,j,k)
		// Left Wall Gohst Cells 
		if(i==0)
		{
		      index1 = 	no_of_Cells
			    +	Ghost_Cell_Bottom_Plane
			    +	Ghost_Cells_Right_Plane
			    +	Ghost_Cell_Top_Plane
			    +	nrwgc++;	//Number of Right Wall Ghost Cells
			Inlet_Boundary_List.push_back(index0);
			Inlet_Boundary_List.push_back(0);
			Inlet_Boundary_List.push_back(index1);
		}
		else
		  index1 = (i-1)+j*(nx-1);//left cell	(i-1,j,k)
		
		if(j==0)
		{
			index2 = no_of_Cells + nbwgc++;		//nbwgc - No of Bottom wall ghost cells
 			if(i<=60)
			{
				Symmetry_Boundary_List.push_back(index0);
				Symmetry_Boundary_List.push_back(1);
				Symmetry_Boundary_List.push_back(index2);				
			}
			if(i>60 & i<=120)
			{
				Wall_Boundary_List.push_back(index0);
				Wall_Boundary_List.push_back(1);
				Wall_Boundary_List.push_back(index2);
			}
			else
			{
				Symmetry_Boundary_List.push_back(index0);
				Symmetry_Boundary_List.push_back(1);
				Symmetry_Boundary_List.push_back(index2);
			}
			
		}
		else
			index2= i+(j-1)*(nx-1);//bottom cell (i,j-1,k)
		
		if(j==(ny-2))				// Top Face with Face Number 3
		{				
		  index4 =	no_of_Cells
			    +	Ghost_Cell_Bottom_Plane
			    +	Ghost_Cells_Left_Plane
			    +	ntwgc++;
 			if(i<=40)
			{
				Wall_Boundary_List.push_back(index0);
				Wall_Boundary_List.push_back(3);
				Wall_Boundary_List.push_back(index4);				
			}
/*			else if(i>=19 and i<= 100)
			{
				Inlet_Boundary_List.push_back(index0);
				Inlet_Boundary_List.push_back(3);
				Inlet_Boundary_List.push_back(index4);
			}*/
			else
			{
				Wall_Boundary_List.push_back(index0);
				Wall_Boundary_List.push_back(3);
				Wall_Boundary_List.push_back(index4);
			}
		}
		else
			index4 = i+(j+1)*(nx-1);//top cell	(i,j+1,k)
// 	cout<<index0<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\n";
				a = index0*5;
				Cell_Neighbour_Data[a+0]=index0;	//Current Cell		(i,j,k)
				Cell_Neighbour_Data[a+1]=index1;	//Left Cell		(i-1,j,k)
				Cell_Neighbour_Data[a+2]=index2;	//Bottom Cell		(i,j-1,k)
				Cell_Neighbour_Data[a+3]=index3;	//Right Cell		(i+1,j,k)
				Cell_Neighbour_Data[a+4]=index4;	//Top Cell		(i,j+1,k)
				
			}
		}
		cout<<"Identifying Neighbours done \n";
}

// This function identifies the neighbours for given cell
void Identify_Neighbours(int &nx1,int &ny1,int &nx2,int &ny2)
{
	
	int nbwgc,ntwgc,nlwgc,nrwgc,nwgc,nigc,nogc,ngc,a,index0,index1,index2,index3,index4;
	int Ghost_Cells_Top_Plane,Ghost_Cells_Bottom_Plane;
	int Ghost_Cells_Left_Plane,Ghost_Cells_Right_Plane,cells_in_Plane,no_of_Cells;
	Ghost_Cells_Top_Plane=0;Ghost_Cells_Bottom_Plane=0;
	Ghost_Cells_Left_Plane=0;Ghost_Cells_Right_Plane=0;
	cells_in_Plane=0;
	nwgc=0;nigc=0;nogc=0;ngc=0,a=0;
	nbwgc=0;nrwgc=0;ntwgc=0;nlwgc=0;
	no_of_Cells = (nx1-1)*(ny1-1) + (nx2-1)*(ny2-ny1);
	Cell_Neighbour_Data.resize(no_of_Cells*5,0);
	cells_in_Plane = (nx1-1)*(ny1-1) + (nx2-1)*(ny2-1);

	Ghost_Cells_Bottom_Plane=(nx2-1) + (ny1-1) ; // Total Ghost Cells in Bottom plane 
	Ghost_Cells_Top_Plane=(nx2-1);				// Total Ghost Cells in Top Plane
	Ghost_Cells_Left_Plane=(ny2-1);				// Total Ghost Cells in Left Plane
	Ghost_Cells_Right_Plane=(ny2-ny1);			// Total Ghost Cells in Right Plane

/*	cout<<nx1<<"\t"<<ny1<<"\t"<<nx2<<"\t"<<ny2<<endl;
	cout<<"Number of Physical Cells\t"<<no_of_Cells<<endl;
	cout<<"GC Bottom Plane\t"<<Ghost_Cells_Bottom_Plane<<endl;
	cout<<"GC Right Plane\t"<<Ghost_Cells_Right_Plane<<endl;
	cout<<"GC Top Plane \t"<<Ghost_Cells_Top_Plane<<endl;
	cout<<"GC Left Plane\t"<<Ghost_Cells_Left_Plane<<endl;
	
	 cout<<"i"<<"\t"<<"j"<<"\t"<<"index0"<<"\t"<<"index1"<<"\t"<<"index2"<<"\t"<<"index3"<<"\t"<<"index4"<<endl;
*/	
	for(int j=0;j<(ny2-1);j++)
	{
	  for(int i=0;i<(nx2-1);i++)
	  {
		  if(j<(ny1-1))
		  {
			  if(i<(nx1-1))
			  {
//				Current Cell index (i,j)				  
				  index0 = i + j*(nx1-1);
//*********************Check for Left Plane Cells indicated by index 1 *******************
				  if(i==0)
				  {
					  //Region 1 of block 1 for forward step case this is inlet boundary and for shock diffraction case this is wall boundary
					  index1 = no_of_Cells + Ghost_Cells_Bottom_Plane + Ghost_Cells_Right_Plane + Ghost_Cells_Top_Plane +nlwgc++; // Index for Left plane cells
					  Inlet_Boundary_List.push_back(index0);
					  Inlet_Boundary_List.push_back(0);
					  Inlet_Boundary_List.push_back(index1);
					  
//------------------- Condition for Bottom cells (i,j-1) Starts here --------------------------				    				 				  		    
					  if(j>0)
					  {
					  		index2 = i + (j-1)*(nx1-1); // Index for Bottom Plane cells
							index3 = (i+1) + j*(nx1-1); // Index for Right Plane cells
							index4 = i + (j+1)*(nx1-1);	// Index for Top Plane cells
					  }							
					  else
					  {
						  // Region 2 bottom side of block 1 which is wall for forward step case and inlet for shock diffraction case 
						 	index2 = no_of_Cells + nbwgc++; // Index for Bottom Plane cells	
							Wall_Boundary_List.push_back(index0);
							Wall_Boundary_List.push_back(1);
							Wall_Boundary_List.push_back(index2);				  	
							index3 = (i+1) + j*(nx1-1); // Index for Right Plane cells
							index4 = i + (j+1)*(nx1-1);	// Index for Top Plane cells							
					  }
//------------------- Condition for Bottom cells (i,j-1) Ends here --------------------------
				  }
				  else
				  {
					  index1 = (i-1) + j*(nx1-1);// Index for Left plane cells
					  if(j>0)
					  		index2 = i + (j-1)*(nx1-1); // Index for Bottom Plane cells
					  else
					  {
						  // Region 3 for block 1 case which is wall in Forward step case and inlet in shock diffraction case
						 	index2 = no_of_Cells + nbwgc++; // Index for Bottom Plane cells
							Wall_Boundary_List.push_back(index0);
							Wall_Boundary_List.push_back(1);
							Wall_Boundary_List.push_back(index2);					  	
					  }

					  index3 = (i+1) + j*(nx1-1); // Index for Right Plane cells
					  index4 = i + (j+1)*(nx1-1);// Index for Top Plane cells
				  }
//------------------- Condition for Right cells (i+1,j) Starts here --------------------------	
					  if(i==(nx1-2))
					  {
						  index3 = no_of_Cells + nbwgc++;
						  Wall_Boundary_List.push_back(index0);
						  Wall_Boundary_List.push_back(2);
						  Wall_Boundary_List.push_back(index3);
					  }
//------------------- Condition for Right cells (i+1,j) Ends here --------------------------	
					 			  
			  } // if condition ending for (i<(nx1-1))
			  else
				  break;
		  }
		  else
		  {
//				Current Cell index (i,j) in the second block after the step 			   
				  index0 = i + (j-(ny1-1))*(nx2-1) + (nx1-1)*(ny1-1);
//------------------- Condition for Left cells (i-1,j) Starts here --------------------------				    
				  if(i==0)
				  {
					  index1 = no_of_Cells + Ghost_Cells_Bottom_Plane + Ghost_Cells_Right_Plane + Ghost_Cells_Top_Plane +nlwgc++;// Index for Left plane cells
					  Inlet_Boundary_List.push_back(index0);
					  Inlet_Boundary_List.push_back(0);
					  Inlet_Boundary_List.push_back(index1);
				  }
				  else
				  {
					  index1 = (i-1) + (j-(ny1-1))*(nx2-1) + (nx1-1)*(ny1-1);
				  }		
//------------------- Condition for Left cells ends here -------------------------------------				    
/*
*/				  	  
//------------------- Condition for Bottom cells (i,j-1) Starts here --------------------------				    
				  if((j==(ny1-1)) and (i<=(nx1-2)))
				  {
					  index2 = i + (j-1)*(nx1-1);
				  }
				  else if((j==(ny1-1)) and (i<(nx2-1)))
				  {
					  index2 = no_of_Cells + nbwgc++;
					  Wall_Boundary_List.push_back(index0);
					  Wall_Boundary_List.push_back(1);
					  Wall_Boundary_List.push_back(index2);
				  }
				  else
					  index2 = i + ((j-1)-(ny1-1))*(nx2-1) + (nx1-1)*(ny1-1);
//------------------- Condition for Bottom cells Ends here --------------------------		
//------------------- Condition for Right cells (i+1,j) Starts here --------------------------	
				  if(i==(nx2-2))
				  {
					  // Region 5 exit boundary for forward step case and wall boundary for shock diffraction case 
					 index3 = no_of_Cells + Ghost_Cells_Bottom_Plane + nrwgc++;
					 Exit_Boundary_List.push_back(index0);
					 Exit_Boundary_List.push_back(2);
					 Exit_Boundary_List.push_back(index3);
				  }
				  else
				  {
					index3 = (i+1) + (j-(ny1-1))*(nx2-1) + (nx1-1)*(ny1-1);  	
				  }	    
					
//------------------- Condition for Right cells (i+1,j) Ends here --------------------------		
				  		    				  
//------------------- Condition for Top cells (i,j+1) Starts here --------------------------				    				 				  		    
				  if(j==(ny2-2))
				  {
					  //Region 6 Wall boundary for forward step case and exit for shock diffraction case 
					  index4 = no_of_Cells + Ghost_Cells_Bottom_Plane + Ghost_Cells_Right_Plane + ntwgc++;
					  Wall_Boundary_List.push_back(index0);
					  Wall_Boundary_List.push_back(3);
					  Wall_Boundary_List.push_back(index4);
				  }
				  else
				  {
					  index4 = i + ((j+1)-(ny1-1))*(nx2-1) + (nx1-1)*(ny1-1);
				  }
//------------------- Condition for Top cells (i,j+1) Ends here --------------------------				    				  			  
		  }
			  
//		  cout<<i<<"\t"<<j<<"\t"<<index0<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<endl;
			a = index0*5;
			Cell_Neighbour_Data[a+0]=index0;	//Current Cell		(i,j,k)
			Cell_Neighbour_Data[a+1]=index1;	//Left Cell			(i-1,j,k)
			Cell_Neighbour_Data[a+2]=index2;	//Bottom Cell		(i,j-1,k)
			Cell_Neighbour_Data[a+3]=index3;	//Right Cell		(i+1,j,k)
			Cell_Neighbour_Data[a+4]=index4;	//Top Cell			(i,j+1,k)
	  }
  }
		cout<<"Identifying Neighbours done \n";
}

void write_inputfile(int & nx,int & ny,vector<Point> & Point_List,string & opfilename)
{
	double size_of_pointlist;
	size_of_pointlist = Point_List.size();
	cout<<size_of_pointlist<<endl;
	cout<<Point_List.size()<<"\t"<<Cell_Neighbour_Data.size()<<"\t"<<Cell_Point_Data.size()<<endl;
	cout<<opfilename<<endl;
	ofstream outfile(opfilename.c_str(),ios::out);
	int no_of_Cells,k,j=0;
	no_of_Cells=(nx-1)*(ny-1);
	Point P;
	if(outfile.is_open())
	{
		outfile<<0<<endl;				//Grid Type
		outfile<<0<<endl;				//Conversion Type
		outfile<<nx<<"\t"<<ny<<endl;
		outfile<<no_of_Cells<<endl;
		for(int i=0;i<no_of_Cells;i++)
		{
			j=i*4;
			P=Point_List[Cell_Point_Data[j+0]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			P=Point_List[Cell_Point_Data[j+1]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			P=Point_List[Cell_Point_Data[j+2]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			P=Point_List[Cell_Point_Data[j+3]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			k=i*5;
// 			cout<<i<<"\t"<<j<<"\t"<<k<<endl;
			outfile<<Cell_Neighbour_Data[k+0]<<"\t"<<Cell_Neighbour_Data[k+1]<<"\t"<<
						Cell_Neighbour_Data[k+2]<<"\t"<<Cell_Neighbour_Data[k+3]<<"\t"<<
						Cell_Neighbour_Data[k+4]<<"\n";
		}
		outfile<<"Inlet_Boundary_Cells\t"<<0<<"\t"<<Inlet_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Inlet_Boundary_List.size();i+=3)
			outfile<<Inlet_Boundary_List[i+0]<<"\t"<<Inlet_Boundary_List[i+1]<<"\t"<<Inlet_Boundary_List[i+2]<<endl;
		outfile<<"Exit_Boundary_Cells\t"<<1<<"\t"<<Exit_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Exit_Boundary_List.size();i+=3)
			outfile<<Exit_Boundary_List[i+0]<<"\t"<<Exit_Boundary_List[i+1]<<"\t"<<Exit_Boundary_List[i+2]<<endl;
		outfile<<"Wall_Boundary_Cells\t"<<2<<"\t"<<Wall_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Wall_Boundary_List.size();i+=3)
			outfile<<Wall_Boundary_List[i+0]<<"\t"<<Wall_Boundary_List[i+1]<<"\t"<<Wall_Boundary_List[i+2]<<endl;
		outfile<<"Symmetry_Boundary_Cells\t"<<3<<"\t"<<Symmetry_Boundary_List.size()/3<<endl;
		cout<<"Symmetry_Boundary_Cells\t"<<3<<"\t"<<Symmetry_Boundary_List.size()/3<<endl;
		
		for(unsigned int i=0;i<Symmetry_Boundary_List.size();i+=3)
			outfile<<Symmetry_Boundary_List[i+0]<<"\t"<<Symmetry_Boundary_List[i+1]<<"\t"<<Symmetry_Boundary_List[i+2]<<endl;

	cout<<"Input file written for Cube\n";
	}
		else
		{
			cout<<"Could not open file for writing"<<endl;
		}
// 	Append_Boundary_Information(opfilename);
}


void write_inputfile(int & nx1,int & ny1,int & nx2 , int & ny2,vector<Point> & Point_List,string & opfilename)
{
	double size_of_pointlist;
	size_of_pointlist = Point_List.size();
	cout<<size_of_pointlist<<endl;
	cout<<Point_List.size()<<"\t"<<Cell_Neighbour_Data.size()<<"\t"<<Cell_Point_Data.size()<<endl;
	cout<<opfilename<<endl;
	ofstream outfile(opfilename.c_str(),ios::out);
	int no_of_Cells,k,j=0;
	no_of_Cells=(nx1-1)*(ny1-1) + (nx2-1)*(ny2-1);
	cout<<"Total Number of Cells\t"<<no_of_Cells<<endl;
	cout<<"Size of Points List\t"<<Point_List.size()<<endl;
	Point P;
	if(outfile.is_open())
	{
		outfile<<0<<endl;				//Grid Type
		outfile<<0<<endl;				//Conversion Type
		outfile<<nx2<<"\t"<<ny1+ny2-1<<endl;
		outfile<<no_of_Cells<<endl;
//		cout<<"While writing to the file"<<endl;
		for(int i=0;i<no_of_Cells;i++)
		{
			j=i*4;
			P=Point_List[Cell_Point_Data[j+0]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			P=Point_List[Cell_Point_Data[j+1]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			P=Point_List[Cell_Point_Data[j+2]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			P=Point_List[Cell_Point_Data[j+3]];
			outfile<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
			k=i*5;
// 			cout<<i<<"\t"<<j<<"\t"<<k<<endl;
			outfile<<Cell_Neighbour_Data[k+0]<<"\t"<<Cell_Neighbour_Data[k+1]<<"\t"<<
						Cell_Neighbour_Data[k+2]<<"\t"<<Cell_Neighbour_Data[k+3]<<"\t"<<
						Cell_Neighbour_Data[k+4]<<"\n";
//			cout<<Cell_Neighbour_Data[k+0]<<"\t"<<Cell_Neighbour_Data[k+1]<<"\t"<<
//						Cell_Neighbour_Data[k+2]<<"\t"<<Cell_Neighbour_Data[k+3]<<"\t"<<
//						Cell_Neighbour_Data[k+4]<<"\n";			
		}
		outfile<<"Inlet_Boundary_Cells\t"<<0<<"\t"<<Inlet_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Inlet_Boundary_List.size();i+=3)
			outfile<<Inlet_Boundary_List[i+0]<<"\t"<<Inlet_Boundary_List[i+1]<<"\t"<<Inlet_Boundary_List[i+2]<<endl;
		outfile<<"Exit_Boundary_Cells\t"<<1<<"\t"<<Exit_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Exit_Boundary_List.size();i+=3)
			outfile<<Exit_Boundary_List[i+0]<<"\t"<<Exit_Boundary_List[i+1]<<"\t"<<Exit_Boundary_List[i+2]<<endl;
		outfile<<"Wall_Boundary_Cells\t"<<2<<"\t"<<Wall_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Wall_Boundary_List.size();i+=3)
			outfile<<Wall_Boundary_List[i+0]<<"\t"<<Wall_Boundary_List[i+1]<<"\t"<<Wall_Boundary_List[i+2]<<endl;
		outfile<<"Symmetry_Boundary_Cells\t"<<3<<"\t"<<Symmetry_Boundary_List.size()/3<<endl;
		cout<<"Symmetry_Boundary_Cells\t"<<3<<"\t"<<Symmetry_Boundary_List.size()/3<<endl;
		for(unsigned int i=0;i<Symmetry_Boundary_List.size();i+=3)
			outfile<<Symmetry_Boundary_List[i+0]<<"\t"<<Symmetry_Boundary_List[i+1]<<"\t"<<Symmetry_Boundary_List[i+2]<<endl;

	cout<<"Input file written for Cube\n";
	}
		else
		{
			cout<<"Could not open file for writing"<<endl;
		}
// 	Append_Boundary_Information(opfilename);
}



void write_VTK(int &nx,int &ny,vector<Point> &Point_List,string & grid_opfile)
{
	ofstream out;
	double size_of_pointlist;
	size_of_pointlist = Point_List.size();
	cout<<grid_opfile<<endl;
	out.open(grid_opfile.c_str(),ios::out);
	int no_of_Cells;
	no_of_Cells=(nx-1)*(ny-1);
	Point P;
	if(out.is_open())
	{
		out<<"# vtk DataFile Version 2.0"<<endl;
		out<<"Grid Cube"<<endl;
		out<<"ASCII"<<endl;
		out<<"DATASET UNSTRUCTURED_GRID"<<endl;
		out<<"POINTS\t"<<size_of_pointlist<<"\t double"<<endl;
		for(int i=0;i<size_of_pointlist;i++)
		{
			P=Point_List[i];
			out<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
		}
		out<<"CELLS\t"<<no_of_Cells<<"\t"<<no_of_Cells*5<<endl;
		int j=0;
		for(int i=0;i<no_of_Cells;i++)
		{
			j=i*4;
			out<<4<<"\t"<<
			Cell_Point_Data[j+0]<<"\t"<<Cell_Point_Data[j+1]<<"\t"<<
			Cell_Point_Data[j+2]<<"\t"<<Cell_Point_Data[j+3]<<"\t"<<endl;
		}
		out<<"CELL_TYPES \t"<<no_of_Cells<<endl;
		for(int i=0;i<no_of_Cells;i++)
			out<<9<<endl;
		cout<<"Grid VTK file written for Cube\n";
 	}
	else
	{
		cout<<"Could not open file for writing"<<endl;
	}
}

void write_VTK(int &nx1,int &ny1,int &nx2,int &ny2,vector<Point> &Point_List,string & grid_opfile)
{
	ofstream out;
	double size_of_pointlist;
	size_of_pointlist = Point_List.size();
	cout<<grid_opfile<<endl;
	out.open(grid_opfile.c_str(),ios::out);
	int no_of_Cells;
	no_of_Cells=(nx1-1)*(ny1-1) + (nx2-1)*(ny2-1);
	cout<<"Total No of Cells\t"<<no_of_Cells<<endl;
	Point P;
	if(out.is_open())
	{
		out<<"# vtk DataFile Version 2.0"<<endl;
		out<<"Grid Cube"<<endl;
		out<<"ASCII"<<endl;
		out<<"DATASET UNSTRUCTURED_GRID"<<endl;
		out<<"POINTS\t"<<size_of_pointlist<<"\t double"<<endl;
		for(int i=0;i<size_of_pointlist;i++)
		{
			P=Point_List[i];
			out<<P.Get_x()<<"\t"<<P.Get_y()<<"\t"<<P.Get_z()<<endl;
		}
		out<<"CELLS\t"<<no_of_Cells<<"\t"<<no_of_Cells*5<<endl;
		int j=0;
		for(int i=0;i<no_of_Cells;i++)
		{
			j=i*4;
			out<<4<<"\t"<<
			Cell_Point_Data[j+0]<<"\t"<<Cell_Point_Data[j+1]<<"\t"<<
			Cell_Point_Data[j+2]<<"\t"<<Cell_Point_Data[j+3]<<"\t"<<endl;
		}
		out<<"CELL_TYPES \t"<<no_of_Cells<<endl;
		for(int i=0;i<no_of_Cells;i++)
			out<<9<<endl;
		cout<<"Grid VTK file \n";
 	}
	else
	{
		cout<<"Could not open file for writing"<<endl;
	}
}
