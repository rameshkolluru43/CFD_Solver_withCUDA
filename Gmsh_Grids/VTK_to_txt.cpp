#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <map>
#include <set>
#include <vector>
#include <iostream>

// Function to sort edge (pair of nodes) for consistency
std::pair<vtkIdType, vtkIdType> sortEdge(vtkIdType id1, vtkIdType id2) {
    if (id1 > id2) std::swap(id1, id2);
    return {id1, id2};
}

int main(int argc, char *argv[]) {
    // Check for input file
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename.vtk>" << std::endl;
        return EXIT_FAILURE;
    }

    // Step 1: Read the .vtk file
    std::string filename = argv[1];
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // Step 2: Get the unstructured grid
    auto unstructuredGrid = reader->GetOutput();
    if (!unstructuredGrid) {
        std::cerr << "Error: Unable to read the unstructured grid from the file." << std::endl;
        return EXIT_FAILURE;
    }

    // Step 3: Access points and cells
    auto points = unstructuredGrid->GetPoints();
    auto cells = unstructuredGrid->GetCells();
    vtkIdType numCells = unstructuredGrid->GetNumberOfCells();
    vtkIdType numPoints = points->GetNumberOfPoints();

    std::cout << "Number of cells: " << numCells << std::endl;
    std::cout << "Number of points: " << numPoints << std::endl;

    // Step 4: Map edges to cells for neighbor detection
    std::map<std::pair<vtkIdType, vtkIdType>, std::vector<vtkIdType>> edgeToCells;
    for (vtkIdType cellId = 0; cellId < numCells; ++cellId) {
        auto cell = unstructuredGrid->GetCell(cellId);
        vtkIdType numPointsInCell = cell->GetNumberOfPoints();
        auto pointIds = cell->GetPointIds();

        // Create edges for the cell
        for (vtkIdType i = 0; i < numPointsInCell; ++i) {
            vtkIdType id1 = pointIds->GetId(i);
            vtkIdType id2 = pointIds->GetId((i + 1) % numPointsInCell);  // Wrap around for closed loop
            edgeToCells[sortEdge(id1, id2)].push_back(cellId);
        }
    }

    // Step 5: Identify neighbors for each cell
    std::vector<std::set<vtkIdType>> cellNeighbors(numCells);
    for (const auto& [edge, connectedCells] : edgeToCells) {
        if (connectedCells.size() > 1) {  // Shared edge indicates neighboring cells
            for (vtkIdType cell : connectedCells) {
                for (vtkIdType neighbor : connectedCells) {
                    if (cell != neighbor) {
                        cellNeighbors[cell].insert(neighbor);
                    }
                }
            }
        }
    }

    // Output neighbors for each cell
    std::cout << "\nCell Neighbors:\n";
    for (vtkIdType cellId = 0; cellId < numCells; ++cellId) {
        std::cout << "Cell " << cellId << " neighbors: ";
        for (vtkIdType neighbor : cellNeighbors[cellId]) {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;
    }

    // Step 6: Identify boundary edges
    std::cout << "\nBoundary Edges:\n";
    for (const auto& [edge, connectedCells] : edgeToCells) {
        if (connectedCells.size() == 1) {  // Boundary edge has only one connected cell
            std::cout << "Edge (" << edge.first << ", " << edge.second << ") belongs to Cell " << connectedCells[0] << std::endl;
        }
    }

    // Step 7: Access point and cell data for boundary conditions
    auto pointData = unstructuredGrid->GetPointData();
    auto cellData = unstructuredGrid->GetCellData();

    if (pointData) {
        std::cout << "\nPoint Data Arrays:\n";
        for (vtkIdType i = 0; i < pointData->GetNumberOfArrays(); ++i) {
            std::cout << "Array " << i << ": " << pointData->GetArrayName(i) << std::endl;
        }
    }

    if (cellData) {
        std::cout << "\nCell Data Arrays:\n";
        for (vtkIdType i = 0; i < cellData->GetNumberOfArrays(); ++i) {
            std::cout << "Array " << i << ": " << cellData->GetArrayName(i) << std::endl;
        }
    }

    return EXIT_SUCCESS;
}