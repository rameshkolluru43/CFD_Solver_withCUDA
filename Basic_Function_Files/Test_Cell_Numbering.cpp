#include <iostream>
#include <vector>

using namespace std;

// Vectors to store ghost and physical cell indices
vector<int> Left_Boundary_Ghost_Cells, Right_Boundary_Ghost_Cells, Bottom_Boundary_Ghost_Cells, Top_Boundary_Ghost_Cells, Physical_Cells;

void Process_Cells_and_Assign_Neighbours(int &Nx, int &Ny)
{
    // Number of points, physical cells, ghost cells, and total cells
    int Number_Of_Points = Nx * Ny;
    int Number_Physical_Cells = (Nx - 1) * (Ny - 1);
    int Number_Ghost_Cells = 2 * (Nx - 1) + 2 * (Ny - 1);
    int Total_Cells = Number_Physical_Cells + Number_Ghost_Cells;

    cout << "Number of Points\t" << Number_Of_Points << endl;
    cout << "Number of Physical Cells\t" << Number_Physical_Cells << endl;
    cout << "Number of Ghost Cells\t" << Number_Ghost_Cells << endl;

    // Resize the vectors
    Bottom_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Left_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Top_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Right_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Physical_Cells.resize(Number_Physical_Cells, 0);

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

 /*   // Print the results
    cout << "\t";
    for (int i : Bottom_Boundary_Ghost_Cells)
    {
        cout << i << "\t";
    }
    cout << endl;

    for (int j = 0; j < (Ny - 1); j++)
    {
        cout << Left_Boundary_Ghost_Cells[j] << "\t";
        for (int i = 0; i < (Nx - 1); i++)
        {
            cout << Physical_Cells[i + j * (Nx - 1)] << "\t";
        }
        cout << Right_Boundary_Ghost_Cells[j] << endl;
    }

    cout << "\t";
    for (int i : Top_Boundary_Ghost_Cells)
    {
        cout << i << "\t";
    }
    cout << endl;

    // Assign and display neighbors
    cout << "Index and Neighbors:\n";*/
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

            cout << index0 << "\t" << index1 << "\t" << index2 << "\t" << index3 << "\t" << index4 << endl;
        }
    }
}

void Identify_Cells(int &Nx, int &Ny)
{
    // Number of points, physical cells, ghost cells, and total cells
    int Number_Of_Points = Nx * Ny;
    int Number_Physical_Cells = (Nx - 1) * (Ny - 1);
    int Number_Ghost_Cells = 2 * (Nx - 1) + 2 * (Ny - 1);
    int Total_Cells = Number_Physical_Cells + Number_Ghost_Cells;

    cout << "Number of Points\t" << Number_Of_Points << endl;
    cout << "Number of Physical Cells\t" << Number_Physical_Cells << endl;
    cout << "Number of Ghost Cells\t" << Number_Ghost_Cells << endl;

    // Resize the vectors
    Bottom_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Left_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Top_Boundary_Ghost_Cells.resize(Nx - 1, 0);
    Right_Boundary_Ghost_Cells.resize(Ny - 1, 0);
    Physical_Cells.resize(Number_Physical_Cells, 0);

    // Bottom boundary ghost cells
    for (int i = 0; i < (Nx - 1); i++)
    {
        Bottom_Boundary_Ghost_Cells[i] = i;
    }

    // Left boundary ghost cells
    for (int i = 0; i < (Ny - 1); i++)
    {
        Left_Boundary_Ghost_Cells[i] = (Nx - 1) + i *(Nx + 1);
    }

    // Physical cells
    for(int j = 0;j<(Ny - 1); j++)
    {
	    for (int i=0;i<(Nx-1);i++)
	    {
		    Physical_Cells[i + j*(Nx-1)] = i + j*(Nx-1 + 2) + (Nx) ;
	    }
    }

    // Right boundary ghost cells
    for (int i = 0; i < (Ny - 1); i++)
    {
        Right_Boundary_Ghost_Cells[i] = 2* (Nx - 1) + 1 + i *(Nx + 1) ;
    }

    // Top boundary ghost cells
    for (int i = 0; i < (Nx - 1); i++)
    {
        Top_Boundary_Ghost_Cells[i] = Number_Physical_Cells + (Nx - 1) + 2*(Ny -1) + i ;
    }

    // Print the results
//    cout << "Bottom Boundary Ghost Cells:\n";
    cout<<"\t";
    for (int i : Bottom_Boundary_Ghost_Cells)
    {
        cout << i << "\t ";
    }
/*    cout << "\nLeft Boundary Ghost Cells:\n";
    for (int i : Left_Boundary_Ghost_Cells)
    {
        cout << i << " ";
    }*/
//    cout << "\nPhysical Cells:\n";
    // Physical cells
    cout<<endl;
    for(int j = 0;j<(Ny - 1); j++)
    {
	    cout<<Left_Boundary_Ghost_Cells[j]<<"\t";
	    for (int i=0;i<(Nx-1);i++)
	    {
		    cout<<Physical_Cells[i + j*(Nx-1)]<<"\t";
	    }
	    cout<<Right_Boundary_Ghost_Cells[j]<<endl;
    }
    
/*    cout << "\nRight Boundary Ghost Cells:\n";
    for (int i : Right_Boundary_Ghost_Cells)
    {
        cout << i << " ";
    }*/
//    cout << "\nTop Boundary Ghost Cells:\n";
    cout<< "\t";
    for (int i : Top_Boundary_Ghost_Cells)
    {
        cout << i << "\t";
    }
    cout << endl;
}

void Assign_Neighbours(int & Nx, int & Ny)
{
	int index0,index1,index2,index3,index4;
	for(int j =0;j<(Ny-1);j++)
	{
		for(int i=0;i<(Nx-1);i++)
		{
			index0 = Physical_Cells[i + j*(Nx-1)];		// Current cell		(i,j,k)
			if(i==0)
			{
				index1 = Left_Boundary_Ghost_Cells[j];
			}
			else
			{
				index1 = Physical_Cells[(i-1) + j*(Nx-1)];		// Left Cell 		(i-1,j,k)
			}
			if(j==0)
			{
				index2 = Bottom_Boundary_Ghost_Cells[i];
			}
			else
			{
				index2 = Physical_Cells[i+ (j-1)*(Nx-1)];		// Bottom cell 		(i,j-1,k)
			}
			if(i==(Nx-2))
			{
				index3 = Right_Boundary_Ghost_Cells[j];
			}
			else
			{
				index3 = Physical_Cells[(i+1) + j*(Nx-1)];		// Right cell		(i+1,j,k)
			}
			if(j==(Ny-2))
			{
				index4 = Top_Boundary_Ghost_Cells[i];
			}
			else
			{
				index4 = Physical_Cells[i     + (j+1)*(Nx-1)];		// Top cell		(i,j+1,k)
			}
			cout<< index0<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<endl;
		}
	}
}

int main()
{
    int Nx = 65; // Number of points in x-direction
    int Ny = 33; // Number of points in y-direction
        Identify_Cells(Nx, Ny);
        Assign_Neighbours(Nx, Ny);
	Process_Cells_and_Assign_Neighbours(Nx,Ny);
    return 0;
}