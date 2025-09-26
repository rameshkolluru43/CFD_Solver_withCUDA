#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace std;

typedef vector<int> Triangle;
typedef vector<int> Quad;

// Custom hash function for pair<int, int>
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        return hash<T1>()(p.first) ^ hash<T2>()(p.second);
    }
};

// Function to read VTK file
void readVTK(const string& filename, vector<vector<double>>& points, vector<Triangle>& triangles, vector<Quad>& quads) {
    ifstream file(filename);
    string line;
    bool readingPoints = false, readingCells = false;
    int numPoints = 0, numCells = 0;
    
    while (getline(file, line)) {
        if (line.find("POINTS") != string::npos) {
            istringstream ss(line);
            string keyword;
            ss >> keyword >> numPoints;
            readingPoints = true;
            continue;
        }
        if (line.find("CELLS") != string::npos) {
            istringstream ss(line);
            string keyword;
            ss >> keyword >> numCells;
            readingPoints = false;
            readingCells = true;
            continue;
        }
        if (readingPoints) {
            istringstream ss(line);
            vector<double> point(3);
            ss >> point[0] >> point[1] >> point[2];
            points.push_back(point);
        }
        if (readingCells) {
            istringstream ss(line);
            int n, v1, v2, v3, v4;
            ss >> n >> v1 >> v2 >> v3;
            if (n == 3) { // Triangle
                triangles.push_back({v1, v2, v3});
            } else if (n == 4) { // Quadrilateral
                ss >> v4;
                quads.push_back({v1, v2, v3, v4});
            }
        }
    }
}

// Function to convert all triangles to quadrilateral mesh
vector<Quad> convertAllToQuads(const vector<Triangle>& triangles) {
    unordered_map<pair<int, int>, vector<int>, hash_pair> edgeMap;
    vector<Quad> quads;
    vector<bool> used(triangles.size(), false);
    
    // Populate edge map with triangle indices
    for (size_t i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int n1 = triangles[i][j];
            int n2 = triangles[i][(j + 1) % 3];
            if (n1 > n2) swap(n1, n2);
            edgeMap[make_pair(n1, n2)].push_back(i);
        }
    }
    
    // Find adjacent triangles that share an edge and convert all to quads
    for (auto& entry : edgeMap) {
        vector<int>& triIndices = entry.second;
        if (triIndices.size() == 2) { // Shared by two triangles
            int t1 = triIndices[0];
            int t2 = triIndices[1];
            if (!used[t1] && !used[t2]) {
                vector<int> q = {triangles[t1][0], triangles[t1][1], triangles[t2][1], triangles[t2][2]};
                quads.push_back(q);
                used[t1] = used[t2] = true;
            }
        }
    }
    
    return quads;
}

// Function to write VTK file with all quads
void writeVTK(const string& filename, const vector<vector<double>>& points, const vector<Quad>& quads) {
    ofstream file(filename);
    file << "# vtk DataFile Version 2.0\nConverted Mesh with All Quads\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    
    file << "POINTS " << points.size() << " double\n";
    for (auto &p : points) {
        file << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    
    file << "CELLS " << quads.size() << " " << quads.size() * 5 << "\n";
    for (auto &q : quads) {
        file << "4 " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "\n";
    }
    
    file << "CELL_TYPES " << quads.size() << "\n";
    for (size_t i = 0; i < quads.size(); ++i) {
        file << "9\n"; // VTK ID for quadrilateral elements
    }
}

int main() {
    vector<vector<double>> points;
    vector<Triangle> triangles;
    vector<Quad> quads;
    
    // Read the original VTK file
    readVTK("cylinder_mesh.vtk", points, triangles, quads);
    
    // Convert all triangles to quads and append existing quads
    vector<Quad> convertedQuads = convertAllToQuads(triangles);
    quads.insert(quads.end(), convertedQuads.begin(), convertedQuads.end());
    
    // Write the new VTK file with only quads
    writeVTK("cylinder_quad_mesh.vtk", points, quads);
    
    cout << "Converted mesh saved as cylinder_quad_mesh.vtk with all quadrilateral elements." << endl;
    return 0;
}