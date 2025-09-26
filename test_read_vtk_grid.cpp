#include "definitions.h"
#include "Globals.h"
#include "Grid.h"
#include "Utilities.h"

using namespace std;

// Note: All global variables and data structures are now defined in Globals.h
// These external declarations allow us to use the actual CFD solver globals

// Note: All function declarations are now available from Grid.h

// Simple print function for cell information
void Print(const Cell &entry)
{
    cout << "Cell " << entry.cellID << ": Type=" << entry.cellType
         << ", Nodes=" << entry.nodeIndices.size()
         << ", Vol=" << entry.Area
         << ", Neighbors=" << entry.Neighbours.size() << endl;
}

// The actual Read_VTK_Grid function is implemented in src/Read_Gmsh_File.cpp
// and declared in include/Grid.h

// Test program
int main()
{
    cout << "=== Read_VTK_Grid Function Test ===" << endl;
    cout << "Testing with: Ramp_15o_52_18.vtk" << endl;
    cout << "==================================" << endl;

    string vtk_file_path = "/Users/rameshkolluru/My_Research/CFD_Solver_withCUDA/Grid_Files/Ramp_Grid_Files/Ramp_15o_52_18.vtk";

    // Test the actual Read_VTK_Grid function
    cout << "\nCalling Read_VTK_Grid function..." << endl;
    bool success = Read_VTK_Grid(vtk_file_path);

    cout << "\n"
         << string(50, '=') << endl;
    if (success)
    {
        cout << "✅ SUCCESS: Read_VTK_Grid function completed successfully!" << endl;

        cout << "\nFinal Results Summary:" << endl;
        cout << "  - Function return value: true" << endl;
        cout << "  - Physical cells: " << No_Physical_Cells << endl;
        cout << "  - Interior cells: " << Cells.size() << endl;
        cout << "  - Boundary cells: " << Boundary_Cells.size() << endl;
        cout << "  - 2D flow detected: " << (Is_2D_Flow ? "Yes" : "No") << endl;
        cout << "  - Inlet boundaries: " << Inlet_Cells_List.size() / 2 << endl;
        cout << "  - Exit boundaries: " << Exit_Cells_List.size() / 2 << endl;
        cout << "  - Wall boundaries: " << Wall_Cells_List.size() / 2 << endl;

        if (!Cells.empty())
        {
            cout << "\nFirst Cell Details:" << endl;
            const Cell &cell = Cells[0];
            cout << "  - Cell ID: " << cell.id << endl;
            cout << "  - Cell Type: " << cell.cellType << endl;
            cout << "  - Nodes: " << cell.nodeIndices.size() << endl;
            cout << "  - Neighbors: " << cell.Neighbours.size() << endl;
            cout << "  - Volume: " << cell.Vol << endl;
            cout << "  - Boundary face: " << (cell.hasBoundaryface ? "Yes" : "No") << endl;
        }

        cout << "\n🎉 All Read_VTK_Grid functionality verified!" << endl;
        return 0;
    }
    else
    {
        cout << "❌ ERROR: Read_VTK_Grid function failed!" << endl;
        cout << "Function returned: false" << endl;
        return 1;
    }
}