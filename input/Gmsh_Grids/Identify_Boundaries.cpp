#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkFeatureEdges.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkGeometryFilter.h>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>
#include <string>

vtkSmartPointer<vtkUnstructuredGrid> readVTKFile(const std::string& filename) {
    std::ifstream fileCheck(filename);
    if (!fileCheck.good()) {
        std::cerr << "Error: VTK file not found: " << filename << std::endl;
        return nullptr;
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid;

    // Check file extension
    if (filename.substr(filename.find_last_of(".") + 1) == "vtu") {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        grid = reader->GetOutput();
    } else if (filename.substr(filename.find_last_of(".") + 1) == "vtk") {
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        grid = reader->GetOutput();
    } else {
        std::cerr << "Error: Unsupported file format. Use .vtu or .vtk" << std::endl;
        return nullptr;
    }

    if (!grid || grid->GetNumberOfCells() == 0) {
        std::cerr << "Error: Failed to load VTK file or grid is empty." << std::endl;
        return nullptr;
    }

    std::cout << "Loaded grid with " << grid->GetNumberOfCells() << " cells and " 
              << grid->GetNumberOfPoints() << " points." << std::endl;

    return grid;
}

void classifyBoundaries(vtkSmartPointer<vtkUnstructuredGrid> grid,
                        std::unordered_map<std::string, std::vector<std::tuple<double, double, double>>>& boundaryMapping) {
    // Convert UnstructuredGrid to PolyData using vtkGeometryFilter
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(grid);
    geometryFilter->Update();
    vtkSmartPointer<vtkPolyData> polyData = geometryFilter->GetOutput();

    // Extract boundary edges
    vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->SetInputData(polyData);
    featureEdges->BoundaryEdgesOn();
    featureEdges->FeatureEdgesOff();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->Update();

    // Get the boundary points
    vtkSmartPointer<vtkPolyData> boundaryData = featureEdges->GetOutput();
    vtkSmartPointer<vtkPoints> boundaryPoints = boundaryData->GetPoints();

    if (!boundaryPoints || boundaryPoints->GetNumberOfPoints() == 0) {
        std::cerr << "Error: No boundary points found in the grid." << std::endl;
        return;
    }

    // Classify boundary points
    std::vector<std::tuple<double, double, double>> walls, inlets, outlets;
    
    for (vtkIdType i = 0; i < boundaryPoints->GetNumberOfPoints(); i++) {
        double p[3];
        boundaryPoints->GetPoint(i, p);

        if (std::abs(p[0]) < 1e-6) { // Wall at x = 0.0
            walls.emplace_back(p[0], p[1], p[2]);
        } else if (std::abs(p[1]) < 1e-6) { // Inlet at y = 0.0
            inlets.emplace_back(p[0], p[1], p[2]);
        } else if (std::abs(p[1] - 1.0) < 1e-6) { // Outlet at y = 1.0
            outlets.emplace_back(p[0], p[1], p[2]);
        }
    }

    // Store classified boundaries
    boundaryMapping["wall"] = walls;
    boundaryMapping["inlet"] = inlets;
    boundaryMapping["outlet"] = outlets;

    std::cout << "Boundary classification completed:\n"
              << "  Walls: " << walls.size() << "\n"
              << "  Inlets: " << inlets.size() << "\n"
              << "  Outlets: " << outlets.size() << std::endl;
}

// Function to save classified boundary points as VTK file for visualization
void saveBoundaryToVTK(const std::string& filename, const std::vector<std::tuple<double, double, double>>& boundaryPoints) {
    vtkSmartPointer<vtkPolyData> boundaryData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for (const auto& p : boundaryPoints) {
        points->InsertNextPoint(std::get<0>(p), std::get<1>(p), std::get<2>(p));
    }

    boundaryData->SetPoints(points);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(boundaryData);
    writer->Write();
}

// Convert vtkUnstructuredGrid to vtkPolyData
vtkSmartPointer<vtkPolyData> ConvertToPolyData(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(unstructuredGrid);
    geometryFilter->Update();
    return geometryFilter->GetOutput();
}

// Example usage
//vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = ...; // Load your unstructured grid



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <VTK File>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string filename = argv[1];

    // Load the unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = readVTKFile(filename);
    if (!grid || grid->GetNumberOfCells() == 0) {
        std::cerr << "Error: Invalid or empty grid." << std::endl;
        return EXIT_FAILURE;
    }

    // Convert to vtkPolyData before using vtkFeatureEdges
    vtkSmartPointer<vtkPolyData> polyData = ConvertToPolyData(grid);
    if (!polyData || polyData->GetNumberOfPoints() == 0) {
        std::cerr << "Error: Failed to convert UnstructuredGrid to PolyData." << std::endl;
        return EXIT_FAILURE;
    }

    // Boundary classification
    std::unordered_map<std::string, std::vector<std::tuple<double, double, double>>> boundaryMapping;
    classifyBoundaries(grid, boundaryMapping);

    // Save classified boundaries to VTK files
    saveBoundaryToVTK("wall_boundary.vtp", boundaryMapping["wall"]);
    saveBoundaryToVTK("inlet_boundary.vtp", boundaryMapping["inlet"]);
    saveBoundaryToVTK("outlet_boundary.vtp", boundaryMapping["outlet"]);

    std::cout << "Boundary data saved as VTK files.\n";

    return EXIT_SUCCESS;
}