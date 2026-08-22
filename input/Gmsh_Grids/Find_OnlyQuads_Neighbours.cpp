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
#include <fstream>
#include <iomanip> // For precise floating-point formatting
#include <limits>

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
    const string& outputFile
) {
    map<Edge, int> edgeCount;
    set<int> boundaryCells;
    
    // Bounding box values
    double x_min = numeric_limits<double>::max();
    double x_max = numeric_limits<double>::lowest();
    double y_min = numeric_limits<double>::max();
    double y_max = numeric_limits<double>::lowest();
    
    cout<<x_min<<"\t"<<x_max<<"\t"<<y_min<<"\t"<<y_max<<endl;

    // Compute bounding box
    for (const auto& p : points) {
        x_min = min(x_min, p.first);
        x_max = max(x_max, p.first);
        y_min = min(y_min, p.second);
        y_max = max(y_max, p.second);
    }
    
    // Determine boundary type
    if (all_of(x_coords.begin(), x_coords.end(), [&](double x) { return abs(x - x_min) < 1e-3; })) {
        boundaryClassification[i] = "Left";
    } else if (all_of(x_coords.begin(), x_coords.end(), [&](double x) { return abs(x - x_max) < 1e-3; })) {
        boundaryClassification[i] = "Right";
    } else if (all_of(y_coords.begin(), y_coords.end(), [&](double y) { return abs(y - y_min) < 1e-3; })) {
        boundaryClassification[i] = "Bottom";
    } else if (all_of(y_coords.begin(), y_coords.end(), [&](double y) { return abs(y - y_max) < 1e-3; })) {
        boundaryClassification[i] = "Top";
    } else {
        // Instead of Unknown, try finding the closest boundary
        double left_dist = *min_element(x_coords.begin(), x_coords.end()) - x_min;
        double right_dist = x_max - *max_element(x_coords.begin(), x_coords.end());
        double bottom_dist = *min_element(y_coords.begin(), y_coords.end()) - y_min;
        double top_dist = y_max - *max_element(y_coords.begin(), y_coords.end());

        // Find the smallest distance
        double min_dist = min({left_dist, right_dist, bottom_dist, top_dist});
    
        if (min_dist == left_dist) boundaryClassification[i] = "Left";
        else if (min_dist == right_dist) boundaryClassification[i] = "Right";
        else if (min_dist == bottom_dist) boundaryClassification[i] = "Bottom";
        else if (min_dist == top_dist) boundaryClassification[i] = "Top";
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
        if (pair.second == 1) { // Appears only once
            boundaryEdges.insert(pair.first);
        }
    }

    // Identify boundary cells
    map<int, string> boundaryClassification;
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto& cell = cells[i];
        int num_pts = cell.size();
        set<Edge> cellEdges;
        vector<double> x_coords, y_coords;

        for (int j = 0; j < num_pts; ++j) {
            Edge edge(cell[j], cell[(j + 1) % num_pts]);
            if (boundaryEdges.find(edge) != boundaryEdges.end()) {
                cellEdges.insert(edge);
            }
            x_coords.push_back(points[cell[j]].first);
            y_coords.push_back(points[cell[j]].second);
        }

        if (cellEdges.size() == 2) { // Cells with exactly 2 boundary edges
            boundaryCells.insert(i);

            // Determine boundary type
            if (all_of(x_coords.begin(), x_coords.end(), [&](double x) { return x <= x_min + 1e-6; })) {
                boundaryClassification[i] = "Left";
            } else if (all_of(x_coords.begin(), x_coords.end(), [&](double x) { return x >= x_max - 1e-6; })) {
                boundaryClassification[i] = "Right";
            } else if (all_of(y_coords.begin(), y_coords.end(), [&](double y) { return y <= y_min + 1e-6; })) {
                boundaryClassification[i] = "Bottom";
            } else if (all_of(y_coords.begin(), y_coords.end(), [&](double y) { return y >= y_max - 1e-6; })) {
                boundaryClassification[i] = "Top";
            } else {
                boundaryClassification[i] = "Unknown";
            }
        }
    }

    // Append boundary information to the output file
    ofstream outFile(outputFile, ios::app);
    if (!outFile) {
        cerr << "Error opening output file." << endl;
        return;
    }
    outFile << "\nBoundary Cell Classification:\n";
    for (const auto& entry : boundaryClassification) {
        outFile << "Cell " << entry.first << ": " << entry.second << "\n";
    }
    outFile.close();
}


// Struct to store Quad cell information
struct QuadCellInfo {
    vector<array<double, 3>> vertices; // Stores (x, y, z) coordinates of cell vertices
    vector<int> neighbors; // Stores neighboring quad cell IDs
};

/*
// Function to find neighbors among quads
map<int, QuadCellInfo> FindQuadNeighbors(vtkSmartPointer<vtkUnstructuredGrid> mesh) {
    int numCells = mesh->GetNumberOfCells();
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

    map<int, QuadCellInfo> quadCells;
    vector<int> quadCellIds;

    // Extract only quad cells
    for (int cellId = 0; cellId < numCells; ++cellId) {
        vtkSmartPointer<vtkCell> cell = mesh->GetCell(cellId);
        if (cell->GetNumberOfPoints() == 4) { // Only quads
            set<int> pointSet;
            vector<array<double, 3>> vertexCoordinates;

            for (int i = 0; i < cell->GetNumberOfPoints(); ++i) {
                int pointId = cell->GetPointId(i);
                pointSet.insert(pointId);

                double coords[3];
                points->GetPoint(pointId, coords);
                vertexCoordinates.push_back({coords[0], coords[1], coords[2]});
            }

            quadCells[cellId] = {vertexCoordinates, {}};
            quadCellIds.push_back(cellId);
        }
    }

    // Find neighbors among quads only
    for (int i = 0; i < quadCellIds.size(); ++i) {
        int cellId = quadCellIds[i];
        set<int> currentCellPoints;
        for (const auto& v : quadCells[cellId].vertices) {
            currentCellPoints.insert(v[0] * 1000 + v[1] * 100 + v[2]);
        }

        for (int j = 0; j < quadCellIds.size(); ++j) {
            int otherCellId = quadCellIds[j];
            if (cellId == otherCellId) continue;

            set<int> otherCellPoints;
            for (const auto& v : quadCells[otherCellId].vertices) {
                otherCellPoints.insert(v[0] * 1000 + v[1] * 100 + v[2]);
            }

            // Find shared points (common edge)
            set<int> sharedPoints;
            for (int pt : currentCellPoints) {
                if (otherCellPoints.find(pt) != otherCellPoints.end()) {
                    sharedPoints.insert(pt);
                }
            }

            // If exactly 2 points are shared, add as a neighbor
            if (sharedPoints.size() == 2) {
                quadCells[cellId].neighbors.push_back(otherCellId);
            }
        }
    }

    return quadCells;
}
*/

map<int, QuadCellInfo> FindQuadNeighbors(vtkSmartPointer<vtkUnstructuredGrid> mesh,
                                         vector<vector<int>>& cells,
                                         vector<pair<double, double>>& points) {
    int numCells = mesh->GetNumberOfCells();
    map<int, QuadCellInfo> quadCells;
    vector<int> quadCellIds;

    // Extract points
    int numPoints = mesh->GetNumberOfPoints();
    points.resize(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        double p[3];
        mesh->GetPoint(i, p);
        points[i] = {p[0], p[1]};
    }

    // Extract only quad cells
    for (int cellId = 0; cellId < numCells; ++cellId) {
        vtkSmartPointer<vtkCell> cell = mesh->GetCell(cellId);
        if (cell->GetNumberOfPoints() == 4) { // Only quads
            vector<int> cellPoints;
            vector<array<double, 3>> vertexCoordinates;
            for (int i = 0; i < 4; ++i) {
                int pointId = cell->GetPointId(i);
                cellPoints.push_back(pointId);
                double p[3];
                mesh->GetPoint(pointId, p);
                vertexCoordinates.push_back({p[0], p[1], p[2]});
            }
            quadCells[cellId] = {vertexCoordinates, {}};
            quadCellIds.push_back(cellId);
            cells.push_back(cellPoints); // Store cell connectivity
        }
    }
    // Find neighbors among quads only
    for (int i = 0; i < quadCellIds.size(); ++i) {
        int cellId = quadCellIds[i];
        set<int> currentCellPoints;
        for (const auto& v : quadCells[cellId].vertices) {
            currentCellPoints.insert(v[0] * 1000 + v[1] * 100 + v[2]);
        }

        for (int j = 0; j < quadCellIds.size(); ++j) {
            int otherCellId = quadCellIds[j];
            if (cellId == otherCellId) continue;

            set<int> otherCellPoints;
            for (const auto& v : quadCells[otherCellId].vertices) {
                otherCellPoints.insert(v[0] * 1000 + v[1] * 100 + v[2]);
            }

            // Find shared points (common edge)
            set<int> sharedPoints;
            for (int pt : currentCellPoints) {
                if (otherCellPoints.find(pt) != otherCellPoints.end()) {
                    sharedPoints.insert(pt);
                }
            }

            // If exactly 2 points are shared, add as a neighbor
            if (sharedPoints.size() == 2) {
                quadCells[cellId].neighbors.push_back(otherCellId);
            }
        }
    }
    return quadCells;
}

// Function to write results to a text file in the requested format
void WriteToFile(const string& filename, const map<int, QuadCellInfo>& quadCells) {
    ofstream outFile(filename);
    if (!outFile) {
        cerr << "Error: Unable to open file for writing." << endl;
        return;
    }

    // Set floating-point precision
    outFile << fixed << setprecision(5);

    // Write each quad cell's vertices and its neighbors in the requested format
    for (const auto& pair : quadCells) {
        int cellId = pair.first;

        // Write 4 vertices (one per line)
        for (const auto& vertex : pair.second.vertices) {
            outFile << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << endl;
        }

        // Write cell ID followed by its neighbors on the same line
        outFile << cellId;
        for (int neighbor : pair.second.neighbors) {
            outFile << " " << neighbor;
        }
        outFile << endl;
    }

    outFile.close();
    cout << "Results written to " << filename << endl;
}

/*
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

    // Find quad neighbors
    map<int, QuadCellInfo> quadCells = FindQuadNeighbors(mesh);

    // Write results to a text file
    WriteToFile("output.txt", quadCells);

    return 0;
}*/

int main() {
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName("cylinder_mesh.vtk");
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();

    vector<vector<int>> cells;
    vector<pair<double, double>> points;

    // Extract quad neighbors and prepare cells/points
    map<int, QuadCellInfo> quadCells = FindQuadNeighbors(mesh, cells, points);

    // Write quad neighbors to output
    WriteToFile("output.txt", quadCells);

    cout<<cells.size()<<endl;
    cout<<points.size()<<endl;
    // Identify boundary cells and append results to output
    identifyBoundaryCells(cells, points, "output.txt");

    

    return 0;
}


