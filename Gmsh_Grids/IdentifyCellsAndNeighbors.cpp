#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkIdList.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <tuple>
#include <string>
#include <sstream>
#include <iostream>

// Struct to store boundary information
struct BoundaryInfo {
    std::unordered_set<vtkIdType> leftBoundary;
    std::unordered_set<vtkIdType> rightBoundary;
    std::unordered_set<vtkIdType> topBoundary;
    std::unordered_set<vtkIdType> bottomBoundary;
    std::unordered_set<vtkIdType> cylinderBoundary;
};


struct Cell {
    vtkIdType id;
    std::vector<vtkIdType> points;
    std::vector<vtkIdType> neighbors;
    std::vector<std::tuple<double, double, double>> vertices;
};




// Function to extract vertex coordinates for each cell
void extractVertices(vtkSmartPointer<vtkUnstructuredGrid> grid, 
                     std::unordered_map<vtkIdType, std::vector<std::tuple<double, double, double>>>& cellVertices) {
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(i, cellPoints);

        std::vector<std::tuple<double, double, double>> vertices;

        for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); ++j) {
            vtkIdType pointId = cellPoints->GetId(j);
            double p[3];
            grid->GetPoint(pointId, p); // Get (x, y, z) coordinates
            vertices.emplace_back(p[0], p[1], p[2]);
        }

        cellVertices[i] = vertices;
    }
}



// Function to identify cells with 2 vertices (Line Cells)
void identifyLineCells(vtkSmartPointer<vtkUnstructuredGrid> grid, 
                        std::unordered_set<vtkIdType>& lineCells, 
                        double& minX, double& maxX, 
                        double& minY, double& maxY) {
    minX = minY = std::numeric_limits<double>::max();
    maxX = maxY = std::numeric_limits<double>::lowest();

    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(i, cellPoints);

        if (cellPoints->GetNumberOfIds() == 2) {  // Line cell
            lineCells.insert(i);

            // Compute bounding box
            for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); ++j) {
                double p[3];
                grid->GetPoint(cellPoints->GetId(j), p);
                minX = std::min(minX, p[0]);
                maxX = std::max(maxX, p[0]);
                minY = std::min(minY, p[1]);
                maxY = std::max(maxY, p[1]);
            }
        }
    }
}

// Function to classify boundary cells (Left, Right, Top, Bottom, Cylinder)
void classifyBoundaries(vtkSmartPointer<vtkUnstructuredGrid> grid, 
                        std::unordered_set<vtkIdType>& lineCells, 
                        BoundaryInfo& boundaries, 
                        double minX, double maxX, double minY, double maxY) {
    for (vtkIdType cellId : lineCells) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(cellId, cellPoints);

        double x1, y1, x2, y2;
        double p1[3], p2[3];
        grid->GetPoint(cellPoints->GetId(0), p1);
        grid->GetPoint(cellPoints->GetId(1), p2);

        x1 = p1[0]; y1 = p1[1];
        x2 = p2[0]; y2 = p2[1];

        // Classify based on bounding box
        if (x1 == minX && x2 == minX) {
            boundaries.leftBoundary.insert(cellId);
        } 
        else if (x1 == maxX && x2 == maxX) {
            boundaries.rightBoundary.insert(cellId);
        } 
        else if (y1 == maxY && y2 == maxY) {
            boundaries.topBoundary.insert(cellId);
        } 
        else if (y1 == minY && y2 == minY) {
            boundaries.bottomBoundary.insert(cellId);
        } 
        else {
            boundaries.cylinderBoundary.insert(cellId);  // Remaining lines are cylinder boundary
        }
    }
}

void identifyCellsAndNeighbors(vtkSmartPointer<vtkUnstructuredGrid> grid, std::vector<Cell>& cells,
                                std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors) {
					
	std::cout<<"Number of Cells are \t"<<grid->GetNumberOfCells()<<std::endl;
    	for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i)
	{
        	vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
		grid->GetCellPoints(i, cellPoints);
	
		vtkIdType numPoints = cellPoints->GetNumberOfIds(); // Get number of points in this cell
//		std::cout << "Cell " << i << " has " << numPoints << " points.\t";
		
		Cell cell;
		cell.id = i;
		for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); ++j) 
		{
			vtkIdType pointId = cellPoints->GetId(j);	
			cell.points.push_back(pointId);
			double p[3];
			grid->GetPoint(pointId, p);
			//std::cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ") ";
		}

		// Identify neighbors correctly (only for shared edges)
	for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); ++j) 
	{
		vtkIdType pointId1 = cellPoints->GetId(j);
		vtkIdType pointId2 = cellPoints->GetId((j + 1) % cellPoints->GetNumberOfIds());

		vtkSmartPointer<vtkIdList> pointCells1 = vtkSmartPointer<vtkIdList>::New();
		grid->GetPointCells(pointId1, pointCells1);

		vtkSmartPointer<vtkIdList> pointCells2 = vtkSmartPointer<vtkIdList>::New();
		grid->GetPointCells(pointId2, pointCells2);

		// Find common neighboring cells that share both points
		for (vtkIdType k = 0; k < pointCells1->GetNumberOfIds(); ++k) 
		{
		   vtkIdType neighborCandidate = pointCells1->GetId(k);
		   if (neighborCandidate != i && std::find(cell.neighbors.begin(), cell.neighbors.end(), neighborCandidate) == cell.neighbors.end()) 
		   {
		       for (vtkIdType m = 0; m < pointCells2->GetNumberOfIds(); ++m) 
		       {
		          if (pointCells2->GetId(m) == neighborCandidate) 
			  {
		             cell.neighbors.push_back(neighborCandidate);
		                    break;
		          }
		       }
		   }
		 }
	 }		

        cells.push_back(cell);
        neighbors[i] = cell.neighbors;
    }
}


void identifyBoundaries(vtkSmartPointer<vtkUnstructuredGrid> grid, std::unordered_set<vtkIdType>& boundaryCells) {
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(i, cellPoints);

        bool isBoundary = false;
        for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); ++j) {
            vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
            grid->GetPointCells(cellPoints->GetId(j), pointCells);

            if (pointCells->GetNumberOfIds() == 1) {
                isBoundary = true;
                break;
            }
        }

        if (isBoundary) {
            boundaryCells.insert(i);
        }
    }
}

void writeOutput(const std::vector<Cell>& cells, const std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors,
                 const std::unordered_set<vtkIdType>& boundaryCells, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Write header
    file << "Cell Data" << std::endl;

    // Write cell and neighbor information
    for (const auto& cell : cells) {
        file << "Cell " << cell.id << ": Vertices = ";
        for (const auto& point : cell.points) {
            file << point << " ";
        }
        file << "| Neighbors = ";
        for (const auto& neighbor : cell.neighbors) {
            file << neighbor << " ";
        }
        file << std::endl;
    }

    // Write boundary information
    file << "\nBoundary Information:" << std::endl;
    for (const auto& boundaryCell : boundaryCells) {
        file << "Boundary Cell " << boundaryCell << ": \n";
        for (int face = 0; face < 4; ++face) {
            file << "    Face " << face << " | Ghost Cell ID = " << boundaryCell + 1000 + face << "\n";
        }
    }

    file.close();
    std::cout << "Output written to " << filename << std::endl;
}

void write_Output(const std::vector<Cell>& cells, const std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors,
                 const std::unordered_set<vtkIdType>& boundaryCells, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Write header
    file << "Cell Data (Vertices, Cell ID, Neighbors)" << std::endl;

    // Write cell and neighbor information
    for (const auto& cell : cells) {
        // Write cell vertices
        for (const auto& point : cell.points) {
            file << point << " ";
        }
        // Write cell ID
        file << "| Cell ID = " << cell.id;

        // Write neighbors
        file << " | Neighbors = ";
        for (const auto& neighbor : cell.neighbors) {
            file << neighbor << " ";
        }
        file << std::endl;
    }

    // Write boundary information
    file << "\nBoundary Information:" << std::endl;
    for (const auto& boundaryCell : boundaryCells) {
        file << "Boundary Cell " << boundaryCell << ": \n";
        for (int face = 0; face < 4; ++face) {
            file << "    Face " << face << " | Ghost Cell ID = " << boundaryCell + 1000 + face << "\n";
        }
    }

    file.close();
    std::cout << "Output written to " << filename << std::endl;
}


// Function to write output (including vertices)
void write_Output(const std::vector<Cell>& cells, 
                  const std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors,
                  const std::unordered_set<vtkIdType>& boundaryCells, 
                  const std::unordered_map<vtkIdType, std::vector<std::tuple<double,
		  double, double> > >& cellVertices, 
                  const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Write header
    file << "Cell Data (Vertices, Cell ID, Neighbors)" << std::endl;

    // Write cell, vertex, and neighbor information
    for (const auto& cell : cells) {
        file << "Cell ID = " << cell.id << " | Vertices: ";

        // Check if cell ID exists in vertex map
        if (cellVertices.find(cell.id) != cellVertices.end()) {
            for (const auto& vertex : cellVertices.at(cell.id)) {
                file << "(" << std::get<0>(vertex) << ", " 
                     << std::get<1>(vertex) << ", " 
                     << std::get<2>(vertex) << ") ";
            }
        } else {
            file << "No vertex data ";
        }

        // Write cell neighbors
        file << " | Neighbors = ";
        for (const auto& neighbor : cell.neighbors) {
            file << neighbor << " ";
        }

        file << std::endl;
    }

    // Write boundary information
    file << "\nBoundary Information:" << std::endl;
    for (const auto& boundaryCell : boundaryCells) {
        file << "Boundary Cell " << boundaryCell << ": \n";
        for (int face = 0; face < 6; ++face) {  // Adjusted for 3D (6 faces instead of 4)
            file << "    Face " << face << " | Ghost Cell ID = " << boundaryCell + 1000 + face << "\n";
        }
    }

    file.close();
    std::cout << "Output written to " << filename << std::endl;
}


// Function to write results to a file
void writeBoundaryOutput(const std::unordered_set<vtkIdType>& lineCells, 
                 const BoundaryInfo& boundaries, 
                 double minX, double maxX, double minY, double maxY, 
                 const std::string& filename) {
    std::ofstream file(filename,std::ios::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    file << "Bounding Box: \n";
    file << "MinX = " << minX << ", MaxX = " << maxX << "\n";
    file << "MinY = " << minY << ", MaxY = " << maxY << "\n\n";

    file << "Line Cells (2 Vertices): \n";
    for (const auto& cell : lineCells) {
        file << "Cell " << cell << "\n";
    }

// Print number of cells in each boundary before listing them
    file << "Left Boundary Cells [" << boundaries.leftBoundary.size() << "]:\n ";
    for (const auto& cell : boundaries.leftBoundary) file << cell << " ";
    file << "\nRight Boundary Cells [" << boundaries.rightBoundary.size() << "]:\n ";
    for (const auto& cell : boundaries.rightBoundary) file << cell << " ";
    file << "\nTop Boundary Cells [" << boundaries.topBoundary.size() << "]: ";
    for (const auto& cell : boundaries.topBoundary) file << cell << " ";
    file << "\nBottom Boundary Cells [" << boundaries.bottomBoundary.size() << "]:\n ";
    for (const auto& cell : boundaries.bottomBoundary) file << cell << " ";
    file << "\nCylinder Boundary Cells [" << boundaries.cylinderBoundary.size() << "]:\n ";
    for (const auto& cell : boundaries.cylinderBoundary) file << cell << " ";
    file << "\n";
    
    file.close();
    std::cout << "Appended Boundary information to " << filename << std::endl;
}


// Function to write results to a file
void writeOutput(const std::vector<Cell>& cells, 
                 const std::unordered_map<vtkIdType, std::vector<vtkIdType>>& neighbors,
                 const std::unordered_set<vtkIdType>& boundaryCells,
                 const std::unordered_set<vtkIdType>& lineCells,
                 const std::unordered_map<vtkIdType,std::vector<std::tuple<double, double, double>>> cellVertices, 
                 const BoundaryInfo& boundaries, 
                 double minX, double maxX, double minY, double maxY, 
                 const std::string& filename, 
                 const std::unordered_map<vtkIdType, vtkIdType>& lineCellToParentQuad) {

    std::ofstream file(filename, std::ios::out | std::ios::trunc);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Compute total number of boundary cells
    size_t totalBoundaryCells = boundaries.leftBoundary.size()
                               + boundaries.rightBoundary.size()
                               + boundaries.topBoundary.size()
                               + boundaries.bottomBoundary.size()
                               + boundaries.cylinderBoundary.size();

    // Compute number of quads, triangles, and lines
    size_t numQuads = 0, numTriangles = 0, numLines = 0, totalCells = 0;
    for (const auto& cell : cells) {
        if (cell.points.size() == 4) numQuads++;
        else if (cell.points.size() == 3) numTriangles++;
        else if (cell.points.size() == 2) numLines++;
    }
    totalCells = cells.size();

    file << "Total Number of Cells: " << totalCells << "\n";
    file << "Total Number of Boundary Cells: " << totalBoundaryCells << "\n";
    file << "Number of Lines: " << numLines << "\n";
    file << "Number of Triangles: " << numTriangles << "\n";
    file << "Number of Quads: " << numQuads << "\n\n";

    file << "Bounding Box: \n";
    file << "MinX = " << minX << ", MaxX = " << maxX << "\n";
    file << "MinY = " << minY << ", MaxY = " << maxY << "\n\n";

    file << "Cell Data:\n";
    for (const auto& cell : cells) {
        if (cellVertices.find(cell.id) != cellVertices.end()) {
            for (const auto& vertex : cellVertices.at(cell.id)) {
                file << "    (" << std::get<0>(vertex) << ", "
                     << std::get<1>(vertex) << ", "
                     << std::get<2>(vertex) << ")\n";
            }
        } else {
            file << "    No vertex data\n";
        }
        file << cell.id << "\t";
        for (const auto& neighbor : cell.neighbors) {
            file << "    " << neighbor << "\t";
        }
        file << "\n";
    }

    file << "\nBoundary Classification:\n";
    auto writeBoundary = [&](const std::unordered_set<vtkIdType>& boundarySet, const std::string& name, int faceNumber) {
        file << name << " (" << boundarySet.size() << "): \n";
        for (const auto& cell : boundarySet) {
            vtkIdType parentCellId = lineCellToParentQuad.count(cell) ? lineCellToParentQuad.at(cell) : cell;
//            file << "    Parent Cell ID: " << parentCellId << " | Face Number: " << faceNumber << " | Neighbor: " << (totalCells++) << "\n";
            file << "    Line Cell ID: " << cell << " | Parent Cell ID: " << parentCellId << " | Face Number: " << faceNumber << " | Neighbor: " << (totalCells++) << "\n";
	    
        }
    };

    writeBoundary(boundaries.leftBoundary, "Left Boundary Cells", 0);
    writeBoundary(boundaries.rightBoundary, "Right Boundary Cells", 2);
    writeBoundary(boundaries.topBoundary, "Top Boundary Cells", 4);
    writeBoundary(boundaries.bottomBoundary, "Bottom Boundary Cells", 1);
    writeBoundary(boundaries.cylinderBoundary, "Cylinder Boundary Cells", 2);

    file.close();
    std::cout << "Written mesh and boundary information to " << filename << std::endl;
}

// Function to map line cells to their corresponding parent quad cells
void mapLineCellsToParentQuads(vtkSmartPointer<vtkUnstructuredGrid> grid, 
			const std::unordered_set<vtkIdType>& lineCells,
			std::unordered_map<vtkIdType, vtkIdType> lineCellToParentQuad) 
{
    for (const auto& lineCell : lineCells) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(lineCell, cellPoints);

        for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); ++j) {
            vtkIdType pointId = cellPoints->GetId(j);
            vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
            grid->GetPointCells(pointId, pointCells);

            for (vtkIdType k = 0; k < pointCells->GetNumberOfIds(); ++k) {
                vtkIdType candidateParent = pointCells->GetId(k);
                vtkSmartPointer<vtkIdList> parentPoints = vtkSmartPointer<vtkIdList>::New();
                grid->GetCellPoints(candidateParent, parentPoints);

                if (parentPoints->GetNumberOfIds() == 4 && candidateParent != lineCell) {
                    lineCellToParentQuad[lineCell] = candidateParent;
                    break;
                }
            }
        }
    }
}


// Function to map line cells to their corresponding parent quad cells (based on shared face)
std::unordered_map<vtkIdType, vtkIdType> mapLineCellsToParentQuads(vtkSmartPointer<vtkUnstructuredGrid> grid, 
                                                                   const std::unordered_set<vtkIdType>& lineCells) {
    std::unordered_map<vtkIdType, vtkIdType> lineCellToParentQuad;

    for (const auto& lineCell : lineCells) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellPoints(lineCell, cellPoints);

        if (cellPoints->GetNumberOfIds() != 2) continue; // Ensure it's a line cell

        vtkIdType point1 = cellPoints->GetId(0);
        vtkIdType point2 = cellPoints->GetId(1);

        vtkSmartPointer<vtkIdList> point1Cells = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> point2Cells = vtkSmartPointer<vtkIdList>::New();
        grid->GetPointCells(point1, point1Cells);
        grid->GetPointCells(point2, point2Cells);

        // Find common cells containing both points
        for (vtkIdType i = 0; i < point1Cells->GetNumberOfIds(); ++i) {
            vtkIdType candidateParent = point1Cells->GetId(i);
            if (candidateParent == lineCell) continue;

            vtkSmartPointer<vtkIdList> parentPoints = vtkSmartPointer<vtkIdList>::New();
            grid->GetCellPoints(candidateParent, parentPoints);

            if (parentPoints->GetNumberOfIds() == 4) { // It's a quad cell
                // Check if the quad cell contains both points of the line cell
                bool hasPoint1 = false, hasPoint2 = false;
                for (vtkIdType j = 0; j < parentPoints->GetNumberOfIds(); ++j) {
                    if (parentPoints->GetId(j) == point1) hasPoint1 = true;
                    if (parentPoints->GetId(j) == point2) hasPoint2 = true;
                }
                if (hasPoint1 && hasPoint2) {
                    lineCellToParentQuad[lineCell] = candidateParent;
                    break; // Assign the first quad found
                }
            }
        }
    }
    return lineCellToParentQuad;
}


void Read_OutputMeshData(const std::string& filename, 
                         int& totalCells, int& numQuads, int& numTriangles, int& numLines,
                         std::vector<Cell>& cells, std::vector<int>& boundaryCells) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::cout << "Reading " << filename << std::endl;

    std::string line;
    std::string discard;  // Fix for undeclared variable

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string keyword;
        ss >> keyword;

        if (keyword == "Total") {
            ss >> discard >> discard >> discard >> totalCells;
        } else if (keyword == "Number") {
            std::string type;
            ss >> type;
            if (type == "Quads:") ss >> numQuads;
            else if (type == "Triangles:") ss >> numTriangles;
            else if (type == "Lines:") ss >> numLines;
        } else if (keyword == "Bounding") {
            std::getline(file, line); // Skip MinX, MaxX
            std::getline(file, line); // Skip MinY, MaxY
        } else if (keyword == "Cell") {
            Cell cell;
            int cellID;
            ss >> discard >> cellID;

            for (int j = 0; j < 4; ++j) {
                double x, y, z;
                file >> x >> y >> z;
//                cell.vertices.push_back(x);
//                cell.vertices.push_back(y);
//                cell.vertices.push_back(z);
	        cell.vertices.push_back(std::make_tuple(x, y, z)); // ✅ Fixed
		
            }

            std::getline(file, line);
            std::istringstream ss_neighbors(line);
            ss_neighbors >> discard;  // Fix for undeclared variable
            int neighbor;
            while (ss_neighbors >> neighbor) {
                cell.neighbors.push_back(neighbor);
            }
            cells.push_back(cell);
        } else if (keyword == "Boundary") {
            std::getline(file, line);
            std::istringstream ss_boundary(line);
            int boundaryCell;
            while (ss_boundary >> boundaryCell) {
                boundaryCells.push_back(boundaryCell);
            }
        }
    }

    file.close();
    std::cout << "Successfully read output mesh data." << std::endl;
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.vtk>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string fileName = argv[1];
    std::string readFileName = argv[2];

    // Read the .vtk file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();

    // Process the grid
    std::vector<Cell> cells;
    std::unordered_map<vtkIdType, std::vector<vtkIdType>> neighbors;
    std::unordered_set<vtkIdType> boundaryCells;
    std::unordered_set<vtkIdType> lineCells;
    std::unordered_map<vtkIdType,std::vector<std::tuple<double, double, double>>> cellVertices; 

    identifyCellsAndNeighbors(grid, cells, neighbors);
    identifyBoundaries(grid, boundaryCells);
    extractVertices(grid, cellVertices);

    double minX, maxX, minY, maxY;
    identifyLineCells(grid, lineCells, minX, maxX, minY, maxY);

    BoundaryInfo boundaries;
    classifyBoundaries(grid, lineCells, boundaries, minX, maxX, minY, maxY);

    std::unordered_map<vtkIdType, vtkIdType> lineCellToParentQuad;

//    mapLineCellsToParentQuads(grid,lineCells,lineCellToParentQuad);
    lineCellToParentQuad = mapLineCellsToParentQuads(grid, lineCells);
    
    writeOutput(cells, neighbors, boundaryCells, lineCells, cellVertices, boundaries, minX, maxX, minY, maxY, "output_mesh_data.txt", lineCellToParentQuad);
    
    cells.clear();
    boundaryCells.clear();
    std::cout<<cells.size()<<std::endl;
    std::cout<<boundaryCells.size()<<std::endl;
    // **Convert unordered_set to vector<int>** ✅
        std::vector<int> boundaryCellsVector(boundaryCells.begin(), boundaryCells.end());
    
    int totalCells=0,numQuads=0,numTriangles=0,numLines=0;
    
    Read_OutputMeshData(readFileName, totalCells, numQuads,numTriangles, numLines,cells,boundaryCellsVector);
    std::cout<<cells.size()<<std::endl;
    std::cout<<boundaryCellsVector.size()<<std::endl;
				 
    return EXIT_SUCCESS;
}