#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>
#include <vtkTetra.h>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

// Struct to store cell data (anticlockwise ordering)
struct CellInfo {
    int cellId;
    vector<int> points; // Anticlockwise O, A, B, C
    vector<int> neighbors; // Neighboring cell IDs
};

// Function to read VTK file and extract cell connectivity
void ReadVTKFile(const string &filename) {
    // Load VTK file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
    if (!mesh) {
        cerr << "Error: Could not read VTK file!" << endl;
        return;
    }

    cout << "Total number of cells: " << mesh->GetNumberOfCells() << endl;

    // Data structures to store cell connectivity
    map<int, CellInfo> cellData;

    // Iterate over each cell
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell *cell = mesh->GetCell(i);
        vtkIdList *pointIds = cell->GetPointIds();

        CellInfo cellInfo;
        cellInfo.cellId = i;

        // Store point IDs in anticlockwise order
        for (vtkIdType j = 0; j < pointIds->GetNumberOfIds(); j++) {
            cellInfo.points.push_back(pointIds->GetId(j));
        }

        // Ensure anticlockwise order for different cell types
        if (cell->GetCellType() == VTK_QUAD || cell->GetCellType() == VTK_TRIANGLE) {
            // Already in anticlockwise order
        } else if (cell->GetCellType() == VTK_HEXAHEDRON) {
            // For hexahedral cells, reorder face vertices
            vector<int> reordered = {cellInfo.points[0], cellInfo.points[1], cellInfo.points[3], cellInfo.points[2]};
            cellInfo.points = reordered;
        }

        // Find neighboring cells
        for (vtkIdType j = 0; j < mesh->GetNumberOfCells(); j++) {
            if (i == j) continue;
            vtkCell *neighbor = mesh->GetCell(j);
            vtkIdList *neighborIds = neighbor->GetPointIds();

            int sharedPoints = 0;
            for (vtkIdType a = 0; a < pointIds->GetNumberOfIds(); a++) {
                for (vtkIdType b = 0; b < neighborIds->GetNumberOfIds(); b++) {
                    if (pointIds->GetId(a) == neighborIds->GetId(b)) {
                        sharedPoints++;
                    }
                }
            }
            if (sharedPoints >= 2) { // A valid neighbor shares at least two points
                cellInfo.neighbors.push_back(j);
            }
        }

        cellData[i] = cellInfo;
    }

    // Print cell connectivity
    for (const auto &[id, info] : cellData) {
        cout << "Cell " << id << " (Type: " << mesh->GetCell(id)->GetCellType() << "): ";
        cout << "O=" << info.points[0] << ", A=" << info.points[1] << ", B=" << info.points[2] << ", C=" << info.points[3] << " | ";
        cout << "Neighbors: ";
        for (int n : info.neighbors) {
            cout << n << " ";
        }
        cout << endl;
    }
}

// Main function
int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <VTK filename>" << endl;
        return -1;
    }

    string vtkFile = argv[1];
    ReadVTKFile(vtkFile);
    return 0;
}
