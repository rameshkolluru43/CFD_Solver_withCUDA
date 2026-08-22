
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <limits>
#include <cmath>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkPoints.h>

using namespace std;

// Structure to store a sorted edge representation
struct Edge {
    int p1, p2;
    Edge(int a, int b) {
        if (a < b) { p1 = a; p2 = b; }
        else { p1 = b; p2 = a; }
    }
    bool operator<(const Edge& other) const {
        return (p1 < other.p1) || (p1 == other.p1 && p2 < other.p2);
    }
};

// Function to identify boundary edges and classify them
void identifyBoundaryCells(
    const vector<vector<int>>& cells, 
    const vector<pair<double, double>>& points,
    const string& outputFile,
    const set<int>& cylinderEdges
) {
    map<Edge, int> edgeCount;
    set<int> boundaryCells;
    map<int, string> boundaryClassification;
    map<string, vector<int>> boundaryCellNumbers; // Store classified cells
    map<int, vector<int>> cellPoints; // Store points of each cell

    // Bounding box values
    double x_min = numeric_limits<double>::max();
    double x_max = numeric_limits<double>::lowest();
    double y_min = numeric_limits<double>::max();
    double y_max = numeric_limits<double>::lowest();

    // Compute bounding box
    for (const auto& p : points) {
        x_min = min(x_min, p.first);
        x_max = max(x_max, p.first);
        y_min = min(y_min, p.second);
        y_max = max(y_max, p.second);
    }

    // Count occurrences of each edge
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto& cell = cells[i];
        int num_pts = cell.size();
        for (int j = 0; j < num_pts; ++j) {
            Edge edge(cell[j], cell[(j + 1) % num_pts]);
            edgeCount[edge]++;
        }
    }

    // Identify boundary edges
    set<Edge> boundaryEdges;
    for (const auto& pair : edgeCount) {
        if (pair.second == 1) { // Edge appears only once (boundary)
            boundaryEdges.insert(pair.first);
        }
    }

    // Identify boundary cells
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto& cell = cells[i];
        int num_pts = cell.size();
        set<Edge> cellEdges;
        vector<double> x_coords, y_coords;
        bool isCylinderBoundary = false;

        for (int j = 0; j < num_pts; ++j) {
            Edge edge(cell[j], cell[(j + 1) % num_pts]);
            if (boundaryEdges.find(edge) != boundaryEdges.end()) {
                cellEdges.insert(edge);
            }
            if (cylinderEdges.find(cell[j]) != cylinderEdges.end()) {
                isCylinderBoundary = true;
            }
            x_coords.push_back(points[cell[j]].first);
            y_coords.push_back(points[cell[j]].second);
        }

        if (cellEdges.size() == 2) { // Cells with exactly 2 boundary edges
            boundaryCells.insert(i);
            string boundaryType;

            if (isCylinderBoundary) {
                boundaryType = "CylinderWall";
            } else if (all_of(x_coords.begin(), x_coords.end(), [&](double x) { return fabs(x - x_min) < 1e-3; })) {
                boundaryType = "Left (Inlet)";
            } else if (all_of(x_coords.begin(), x_coords.end(), [&](double x) { return fabs(x - x_max) < 1e-3; })) {
                boundaryType = "Right (Outlet)";
            } else if (all_of(y_coords.begin(), y_coords.end(), [&](double y) { return fabs(y - y_min) < 1e-3; })) {
                boundaryType = "BottomWall";
            } else if (all_of(y_coords.begin(), y_coords.end(), [&](double y) { return fabs(y - y_max) < 1e-3; })) {
                boundaryType = "TopWall";
            } else {
                double left_dist = *min_element(x_coords.begin(), x_coords.end()) - x_min;
                double right_dist = x_max - *max_element(x_coords.begin(), x_coords.end());
                double bottom_dist = *min_element(y_coords.begin(), y_coords.end()) - y_min;
                double top_dist = y_max - *max_element(y_coords.begin(), y_coords.end());

                double min_dist = min({fabs(left_dist), fabs(right_dist), fabs(bottom_dist), fabs(top_dist)});

                if (min_dist == fabs(left_dist)) boundaryType = "Left (Inlet)";
                else if (min_dist == fabs(right_dist)) boundaryType = "Right (Outlet)";
                else if (min_dist == fabs(bottom_dist)) boundaryType = "BottomWall";
                else if (min_dist == fabs(top_dist)) boundaryType = "TopWall";
            }

            boundaryClassification[i] = boundaryType;
            boundaryCellNumbers[boundaryType].push_back(i);
            cellPoints[i] = cell; // Store cell's points
        }
    }

    // Append boundary information to the output file
    ofstream outFile(outputFile, ios::app);
    if (!outFile) {
        cerr << "Error opening output file." << endl;
        return;
    }

    outFile << "
Boundary Cell Classification:
";
    for (const auto& entry : boundaryClassification) {
        outFile << "Cell " << entry.first << ": " << entry.second << "
";
    }

    // Write boundary cell counts
    outFile << "
Boundary Cell Counts:
";
    for (const auto& entry : boundaryCellNumbers) {
        outFile << entry.first << ": " << entry.second.size() << " cells
";
    }

    // Write cell points
    outFile << "
Boundary Cell Points:
";
    for (const auto& entry : boundaryCellNumbers) {
        outFile << entry.first << ":
";
        for (int cellId : entry.second) {
            outFile << "Cell " << cellId << " Points: ";
            for (int pt : cellPoints[cellId]) {
                outFile << pt << " ";
            }
            outFile << "
";
        }
    }

    outFile.close();
}

int main() {
    string geoFile = "cylinder_mesh.geo";
    string vtkFile = "cylinder_mesh.vtk";
    string outputFile = "output.txt";

    vector<vector<int>> cells;
    vector<pair<double, double>> points;
    
    // Read cylinder boundary from .geo file
    set<int> cylinderEdges = readGeoFile(geoFile);

    // Read mesh from .vtk file
    readVTKFile(vtkFile, cells, points);

    // Identify boundary cells and append results to output
    identifyBoundaryCells(cells, points, outputFile, cylinderEdges);

    cout << "Boundary classification completed. Check " << outputFile << " for results." << endl;
    return 0;
}
