#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <iostream>
#include <set>
#include <vector>
#include <map>

using namespace std;

// Function to find neighbors of each cell and print vertices
void FindNeighborsAndPrint(vtkSmartPointer<vtkUnstructuredGrid> mesh) {
    int numCells = mesh->GetNumberOfCells();
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

    // Store cell connectivity and neighbor information
    map<int, vector<int>> cellNeighbors;
    map<int, vector<array<double, 3>>> cellVertices; // Store cell vertex coordinates

    // Loop through all cells and store point connectivity
    map<int, set<int>> cellPoints;
    for (int cellId = 0; cellId < numCells; ++cellId) {
        vtkSmartPointer<vtkCell> cell = mesh->GetCell(cellId);
        set<int> pointSet;
        vector<array<double, 3>> vertexCoordinates;

        for (int i = 0; i < cell->GetNumberOfPoints(); ++i) {
            int pointId = cell->GetPointId(i);
            pointSet.insert(pointId);

            double coords[3];
            points->GetPoint(pointId, coords);
            vertexCoordinates.push_back({coords[0], coords[1], coords[2]});
        }

        cellPoints[cellId] = pointSet;
        cellVertices[cellId] = vertexCoordinates;
    }

    // Find neighbors by checking shared edges
    for (int cellId = 0; cellId < numCells; ++cellId) {
        set<int> currentCellPoints = cellPoints[cellId];

        for (int otherCellId = 0; otherCellId < numCells; ++otherCellId) {
            if (cellId == otherCellId) continue; // Skip self

            set<int> otherCellPoints = cellPoints[otherCellId];

            // Find shared points (common edge)
            set<int> sharedPoints;
            for (int pt : currentCellPoints) {
                if (otherCellPoints.find(pt) != otherCellPoints.end()) {
                    sharedPoints.insert(pt);
                }
            }

            // If exactly 2 points are shared, add as a neighbor (edge-based neighbor detection)
            if (sharedPoints.size() == 2) {
                cellNeighbors[cellId].push_back(otherCellId);
            }
        }
    }

    // Print results
    for (const auto& pair : cellNeighbors) {
        int cellId = pair.first;
        cout << "Cell " << cellId << " Vertices: ";

        // Print cell vertices
        for (const auto& vertex : cellVertices[cellId]) {
            cout << "(" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ") ";
        }
        cout << endl;

        // Print neighbors
        cout << "   Neighbors: ";
        for (int neighbor : pair.second) {
            cout << neighbor << " ";
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <vtk_file>" << endl;
        return -1;
    }

    // Read the VTK file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();

    // Find and print neighbors along with vertex coordinates
    FindNeighborsAndPrint(mesh);

    return 0;
}