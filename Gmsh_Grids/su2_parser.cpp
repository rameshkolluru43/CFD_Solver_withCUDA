#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

// Structure to hold point coordinates
struct Point {
    double x, y, z;
};

// Structure to hold cell information
struct Cell {
    int id;
    std::vector<int> vertices;    // Indices of vertices
    std::vector<int> neighbors;  // Neighboring cells
};

// Function to parse the SU2 file
void parseSU2File(const std::string& filename, 
                  std::vector<Point>& points,
                  std::vector<Cell>& cells,
                  std::unordered_map<int, std::vector<int>>& elementNeighbors) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    bool inPointsSection = false, inElementsSection = false;
    int numPoints = 0, numElements = 0;

    while (std::getline(file, line)) {
        if (line.find("NPOIN=") != std::string::npos) {
            numPoints = std::stoi(line.substr(line.find("=") + 1));
            inPointsSection = true;
        } else if (line.find("NELEM=") != std::string::npos) {
            numElements = std::stoi(line.substr(line.find("=") + 1));
            inPointsSection = false;
            inElementsSection = true;
        } else if (inPointsSection) {
            std::istringstream iss(line);
            Point point;
            iss >> point.x >> point.y >> point.z;
            points.push_back(point);
        } else if (inElementsSection) {
            std::istringstream iss(line);
            int elementType, numVertices;
            iss >> elementType;
            Cell cell;
            cell.id = cells.size();

            while (iss >> numVertices) {
                cell.vertices.push_back(numVertices);
            }
            cells.push_back(cell);
        }
    }
    file.close();

    // Ensure only cells sharing a full edge (not just one node) are neighbors
    for (size_t i = 0; i < cells.size(); ++i) {
        for (size_t j = i + 1; j < cells.size(); ++j) {
            std::unordered_set<int> cell1(cells[i].vertices.begin(), cells[i].vertices.end());
            std::unordered_set<int> cell2(cells[j].vertices.begin(), cells[j].vertices.end());

            // Count how many vertices are shared
            std::vector<int> intersection;
            std::set_intersection(cell1.begin(), cell1.end(), cell2.begin(), cell2.end(),
                                  std::back_inserter(intersection));

            // Two quads are neighbors if they share exactly 2 vertices (1 full edge)
            if (intersection.size() == 2) {
                cells[i].neighbors.push_back(j);
                cells[j].neighbors.push_back(i);
            }
        }
    }
}

// Function to write output in the specified format
void writeOutput(const std::vector<Point>& points,
                 const std::vector<Cell>& cells,
                 const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    for (const auto& cell : cells) {
        // Write vertex coordinates
        for (const auto& vertex : cell.vertices) {
            const auto& point = points[vertex];
            file << point.x << "\t" << point.y << "\t" << point.z << "\n";
        }

        // Write cell ID and neighbors
        file << cell.id;
        for (const auto& neighbor : cell.neighbors) {
            file << "\t" << neighbor;
        }
        file << "\n";
    }

    file.close();
    std::cout << "Output written to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.su2> <output_file.txt>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];

    std::vector<Point> points;
    std::vector<Cell> cells;
    std::unordered_map<int, std::vector<int>> elementNeighbors;

    // Parse SU2 file
    parseSU2File(inputFilename, points, cells, elementNeighbors);

    // Write output file
    writeOutput(points, cells, outputFilename);

    return EXIT_SUCCESS;
}