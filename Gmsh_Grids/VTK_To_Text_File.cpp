#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.vtk> <output.txt>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    // Read the VTK file
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(inputFileName.c_str());
    reader->Update();

    auto unstructuredGrid = reader->GetOutput();
    if (!unstructuredGrid) {
        std::cerr << "Error: Unable to read the unstructured grid from the file." << std::endl;
        return EXIT_FAILURE;
    }

    // Open the output file
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open the output file." << std::endl;
        return EXIT_FAILURE;
    }

    // Write number of points and cells
    auto points = unstructuredGrid->GetPoints();
    vtkIdType numPoints = points->GetNumberOfPoints();
    vtkIdType numCells = unstructuredGrid->GetNumberOfCells();

    outFile << "Number of Points: " << numPoints << std::endl;
    outFile << "Number of Cells: " << numCells << std::endl;

    // Write points
    outFile << "\nPoints:\n";
    for (vtkIdType i = 0; i < numPoints; ++i) {
        double point[3];
        points->GetPoint(i, point);
        outFile << i << " " << point[0] << " " << point[1] << " " << point[2] << "\n";
    }

    // Write cell connectivity
    outFile << "\nCells (Connectivity):\n";
    for (vtkIdType i = 0; i < numCells; ++i) {
        auto cell = unstructuredGrid->GetCell(i);
        auto pointIds = cell->GetPointIds();
        outFile << i << " ";
        for (vtkIdType j = 0; j < pointIds->GetNumberOfIds(); ++j) {
            outFile << pointIds->GetId(j) << " ";
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Output written to " << outputFileName << std::endl;

    return EXIT_SUCCESS;
}