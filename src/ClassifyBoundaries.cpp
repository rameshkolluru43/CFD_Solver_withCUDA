#include "definitions.h"
#include "Globals.h"

// Function to classify and store boundary points
void classifyBoundaries(vtkSmartPointer<vtkUnstructuredGrid> grid,
                        std::unordered_map<std::string, std::vector<std::tuple<double, double, double>>> &boundaryMapping)
{

    // Step 1: Extract boundary edges (grid boundaries)
    vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->SetInputData(grid);
    featureEdges->BoundaryEdgesOn(); // Extract only boundary edges
    featureEdges->FeatureEdgesOff();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->Update();

    // Get the boundary points from extracted edges
    vtkSmartPointer<vtkPolyData> boundaryData = featureEdges->GetOutput();
    vtkSmartPointer<vtkPoints> boundaryPoints = boundaryData->GetPoints();

    if (!boundaryPoints)
    {
        std::cerr << "Error: No boundary points found in the grid." << std::endl;
        return;
    }

    // Step 2: Loop through boundary points and classify them
    std::vector<std::tuple<double, double, double>> walls, inlets, outlets;

    for (vtkIdType i = 0; i < boundaryPoints->GetNumberOfPoints(); i++)
    {
        double p[3];
        boundaryPoints->GetPoint(i, p);

        // Classify based on coordinate values
        if (std::abs(p[0]) < 1e-6)
        { // Wall at x = 0.0
            walls.emplace_back(p[0], p[1], p[2]);
        }
        else if (std::abs(p[1]) < 1e-6)
        { // Inlet at y = 0.0
            inlets.emplace_back(p[0], p[1], p[2]);
        }
        else if (std::abs(p[1] - 1.0) < 1e-6)
        { // Outlet at y = 1.0
            outlets.emplace_back(p[0], p[1], p[2]);
        }
    }

    // Step 3: Store classified boundaries in the hash map
    boundaryMapping["wall"] = walls;
    boundaryMapping["inlet"] = inlets;
    boundaryMapping["outlet"] = outlets;

    std::cout << "Boundary classification completed:\n"
              << "  Walls: " << walls.size() << "\n"
              << "  Inlets: " << inlets.size() << "\n"
              << "  Outlets: " << outlets.size() << std::endl;
}

// Function to save classified boundary points as VTK file for visualization
void saveBoundaryToVTK(const std::string &filename, const std::vector<std::tuple<double, double, double>> &boundaryPoints)
{
    vtkSmartPointer<vtkPolyData> boundaryData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for (const auto &p : boundaryPoints)
    {
        points->InsertNextPoint(std::get<0>(p), std::get<1>(p), std::get<2>(p));
    }

    boundaryData->SetPoints(points);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(boundaryData);
    writer->Write();
}