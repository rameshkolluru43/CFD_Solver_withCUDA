/*#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <fstream>

// Define the Cell structure
struct Cell {
    vtkIdType id;                  // Cell ID
    std::vector<vtkIdType> points; // Points in the cell
    std::vector<vtkIdType> neighbors; // Neighboring cells
};

// Identify cells and neighbors
void identifyCellsAndNeighbors(vtkSmartPointer<vtkUnstructuredGrid> grid,
                                std::vector<Cell>& cells,
                                std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors) {
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); i++) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(i, cellPoints);

        Cell cell;
        cell.id = i;
        for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); j++) {
            cell.points.push_back(cellPoints->GetId(j));
        }

        // Find neighbors
        for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); j++) {
            vtkIdType pointId = cellPoints->GetId(j);
            vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
            grid->GetPointCells(pointId, pointCells);

            for (vtkIdType k = 0; k < pointCells->GetNumberOfIds(); k++) {
                vtkIdType neighborCellId = pointCells->GetId(k);
                if (neighborCellId != i &&
                    std::find(cell.neighbors.begin(), cell.neighbors.end(), neighborCellId) == cell.neighbors.end()) {
                    cell.neighbors.push_back(neighborCellId);
                }
            }
        }

        cells.push_back(cell);
        neighbors[i] = cell.neighbors;
    }
}

// Identify boundary cells
void identifyBoundaries(vtkSmartPointer<vtkUnstructuredGrid> grid,
                        std::unordered_set<vtkIdType>& boundaryCells) {
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); i++) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(i, cellPoints);

        int boundaryFlag = 0;
        for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); j++) {
            vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
            grid->GetPointCells(cellPoints->GetId(j), pointCells);

            if (pointCells->GetNumberOfIds() == 1) {
                boundaryFlag = 1; // This point belongs to the boundary
                break;
            }
        }
        if (boundaryFlag) {
            boundaryCells.insert(i);
        }
    }
}

// Write output to file
void writeOutput(const std::vector<Cell>& cells,
                 const std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors,
                 const std::unordered_set<vtkIdType>& boundaryCells,
                 const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Write cells and neighbors
    for (const auto& cell : cells) {
        file << "Cell " << cell.id << ": Vertices = ";
        for (const auto& point : cell.points) {
            file << point << " ";
        }
        file << "\n    Neighbors = ";
        for (const auto& neighbor : neighbors.at(cell.id)) {
            file << neighbor << " ";
        }
        file << "\n";
    }

    // Write boundary information
    file << "\nBoundary Information:\n";
    for (const auto& boundaryCell : boundaryCells) {
        file << "Boundary Cell " << boundaryCell << ": Faces and Ghost Cells\n";
        for (int face = 0; face < 4; ++face) {
            file << "    Face " << face << ": Ghost Cell ID = " << boundaryCell + 1000 + face << "\n";
        }
    }

    file.close();
    std::cout << "Output written to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.vtk>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string fileName = argv[1];

    // Read the .vtk file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();

    // Identify cells, neighbors, and boundaries
    std::vector<Cell> cells;
    std::unordered_map<vtkIdType, std::vector<vtkIdType>> neighbors;
    std::unordered_set<vtkIdType> boundaryCells;

    identifyCellsAndNeighbors(grid, cells, neighbors);
    identifyBoundaries(grid, boundaryCells);

    // Write the output
    writeOutput(cells, neighbors, boundaryCells, "output_cells_and_boundaries.txt");

    return EXIT_SUCCESS;
}*/