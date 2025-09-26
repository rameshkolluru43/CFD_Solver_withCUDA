#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"
#include "Grid.h"

// The function reads the grid file and stores the points, cells and boundary cells in the respective vectors
// The function is incomplete and not generic to read all kind of cells
// Function to read the Gmsh grid file
// The function reads the grid file and stores the points, cells and boundary cells in the respective vectors

void Read_GmshVTK_Grid(const string &GridFileName)
{
    ifstream file(GridFileName);
    if (!file)
    {
        cerr << "Could not open the file! Check the file path: " << GridFileName << endl;
        exit(1);
    }

    Is_2D_Flow = true;
    string line;
    vector<double> Points;
    int No_of_Points = 0, No_of_Cells = 0, NumEntries = 0;

    // Mapping VTK cell types to node counts
    std::map<int, int> cellTypeToNodeCount = {
        {1, 1}, {3, 2}, {5, 3}, {8, 4}, {9, 8}, {19, 8}, {11, 10}, {12, 6}};

    cout << "Reading " << GridFileName << " ..." << endl;

    while (std::getline(file, line))
    {
        stringstream ss(line);
        string keyword;
        ss >> keyword;

        if (keyword == "POINTS")
        {
            ss >> No_of_Points;
            cout << "Number of Points: " << No_of_Points << endl;

            Points.reserve(3 * No_of_Points);
            double x, y, z;
            for (int i = 0; i < No_of_Points; i++)
            {
                file >> x >> y >> z;
                Points.insert(Points.end(), {x, y, z});
            }
        }
        else if (keyword == "CELLS")
        {
            ss >> No_of_Cells >> NumEntries;
            cout << "Number of Cells: " << No_of_Cells << endl;

            for (int i = 0; i < No_of_Cells; i++)
            {
                int numNodes;
                file >> numNodes;

                vector<int> nodeIndices(numNodes);
                for (int &nodeIndex : nodeIndices)
                {
                    file >> nodeIndex;
                }

                if (numNodes == 2)
                {
                    Boundary_Cells.push_back(Cell(numNodes, i, nodeIndices));
                }
                else if (numNodes == 3 || numNodes == 4)
                {
                    Is_2D_Flow = true;
                    Cells.push_back(Cell(numNodes, i - Boundary_Cells.size(), nodeIndices));
                }
            }
        }
        else if (keyword == "CELL_TYPES")
        {
            cout << "Reading Cell Types" << endl;

            for (int i = 0; i < No_of_Cells; i++)
            {
                int cellType;
                file >> cellType;

                if (cellType == 3)
                {
                    Boundary_Cells[i].cellType = cellType;
                }
                else if (cellType == 5 || cellType == 8)
                {
                    Cells[i - Boundary_Cells.size()].cellType = cellType;
                }
            }
        }
    }

    file.close();

    cout << "Number of Cells of each type: " << endl;
    cout << "Line: " << Boundary_Cells.size() << endl;
    cout << "Triangle/Quadrilateral: " << Cells.size() << endl;

    if (Is_2D_Flow)
    {
        cout << "2D Grid is being read" << endl;
        No_Physical_Cells = Cells.size();
    }
    else
    {
        cout << "3D Grid is being read" << endl;
    }

    Identify_Cells(Points, Cells, Is_2D_Flow, No_Physical_Cells);
    cout << "Identifying Cells is done" << endl;
    Identify_Neighbours(Points, Cells, Boundary_Cells);
    cout << "Identifying Neighbours is done" << endl;
    // Identifying and classifying the boundaries
    double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
    BoundingBox(Points, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
    Classify_Domain_Boundaries(Cells, Xmin, Xmax, Ymin, Ymax);
    Create_Boundary_Cells_Lists(Cells, Inlet_Cells_List, Exit_Cells_List, Wall_Cells_List);
    // Sort_Neighbours(Cells);
    for (size_t i = 0; i < Cells.size(); i++)
    {
        /* code */
        Print(Cells[i]);
    }
}

void Identify_Cells(V_D &Points, vector<Cell> &Cells, bool Is_2D_Flow, int &No_Physical_Cells)
{
    if (Cells.size() == No_Physical_Cells)
        cout << "Number of Physical Cells are correct\n";
    else
    {
        cout << "Number of Physical Cells are not correct\t" << Cells.size() << "\t" << No_Physical_Cells << endl;
        exit(0);
    }
    cout << "Identifying Cells" << endl;
    for (size_t i = 0; i < Cells.size(); i++)
    {
        Cells[i].NoBoundaryFaces = 0;
        if (Is_2D_Flow)
        {
            if (Cells[i].nodeIndices.size() == 3 & Cells[i].numNodes == 3)
            {
                Cells[i].cellType = 5; // Triangle
                Cells[i].cellID = i;
            }
            else if (Cells[i].nodeIndices.size() == 4 & Cells[i].numNodes == 4)
            {
                Cells[i].cellType = 8; // Quadrilateral
                Cells[i].cellID = i;
            }
            else if (Cells[i].nodeIndices.size() == 2 & Cells[i].numNodes == 2)
            {
                Cells[i].cellType = 3; // Line
                Cells[i].cellID = i;
            }
            else
            {
                cout << "Cell Type is not correct\t" << Cells[i].nodeIndices.size() << "\t" << Cells[i].numNodes << endl;
                exit(0);
            }
            Cells[i].Dimension = 2;
            if (Cells[i].nodeIndices.size() == 3)
            {
                Cells[i].numFaces = 3;
                int node1 = Cells[i].nodeIndices[0];
                int node2 = Cells[i].nodeIndices[1];
                int node3 = Cells[i].nodeIndices[2];
                Cells[i].Cell_Vertices.push_back(Points[3 * node1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node1 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node1 + 2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node2 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node2 + 2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node3]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node3 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node3 + 2]);
                // Print(Cells[i].nodeIndices);
                Sort_Points_AntiClockWise(Cells[i].Cell_Vertices, Cells[i].nodeIndices);
                // Print(Cells[i].nodeIndices);
                Construct_Cell(Cells[i]); // Construct the cell
                // Store the vertices of the cell for Triangle only three vertices will be there
            }
            else if (Cells[i].nodeIndices.size() == 4)
            {

                int node1 = Cells[i].nodeIndices[0];
                int node2 = Cells[i].nodeIndices[1];
                int node3 = Cells[i].nodeIndices[2];
                int node4 = Cells[i].nodeIndices[3];
                Cells[i].Cell_Vertices.push_back(Points[3 * node1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node1 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node1 + 2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node2 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node2 + 2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node3]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node3 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node3 + 2]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node4]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node4 + 1]);
                Cells[i].Cell_Vertices.push_back(Points[3 * node4 + 2]);
                // Print(Cells[i].nodeIndices);
                Sort_Points_AntiClockWise(Cells[i].Cell_Vertices, Cells[i].nodeIndices);
                // Print(Cells[i].nodeIndices);
                Construct_Cell(Cells[i]); // Construct the cell
            }
        }
        else
        {
            // 3D Flow
        }
    }
}

void Identify_Neighbours(V_D &Points, vector<Cell> &Cells, vector<Cell> &BLineCells)
{
    std::map<std::pair<int, int>, int> faceToCell;
    std::map<std::pair<int, int>, std::set<int>> faceToCells;
    std::unordered_map<std::pair<int, int>, bool, hash_pair> boundaryFaceMap;

    int tracker = 0;

    cout << "Identifying Neighbours" << endl;
    cout << "Cells Size is \t" << Cells.size() << endl;
    int GhostCellIndex = Cells.size();

    // Maps to store the face to cell mapping
    mapFacesToCells(Cells, faceToCells);

    // Check if all the cells have neighbours equal to the number of nodes otherwise insert ghost cell
    for (size_t i = 0; i < Cells.size(); i++)
    {
        if (Cells[i].Neighbours.size() < Cells[i].nodeIndices.size())
        {
            Cells[i].NoBoundaryFaces = Cells[i].nodeIndices.size() - Cells[i].Neighbours.size();
            do
            {
                Cells[i].Neighbours.push_back(GhostCellIndex);
                Cells[i].hasBoundaryface = true;
                GhostCellIndex++;
            } while (Cells[i].Neighbours.size() < Cells[i].nodeIndices.size());
        }
    }
}

void Identify_ParentCell(vector<Cell> &BCells, vector<Cell> &triQuadCells)
{
    // Identify the parent cell of the boundary cells
    // The parent cell of the boundary cell is the cell which has the boundary cell as its face
    // create a hash map for the boundary cell to its parent cell
    std::unordered_map<int, int> boundaryCellToParentCell;
    // get the faces of the boundary cells and compare it with the faces of the triQuadCells
    for (int i = 0; i < BCells.size(); i++)
    {
        vector<pair<int, int>> faces;
        // Extract faces (edges) for the current cell
        if (BCells[i].nodeIndices.size() == 2) // Line cell
        {
            faces = {
                {min(BCells[i].nodeIndices[0], BCells[i].nodeIndices[1]), max(BCells[i].nodeIndices[0], BCells[i].nodeIndices[1])}};
        }
        // Compare the faces of the boundary cell with the faces of the triQuadCells
        vector<pair<int, int>> triQuadFaces;
        for (int j = 0; j < triQuadCells.size(); j++)
        {
            if (triQuadCells[j].nodeIndices.size() == 3) // Triangle
            {
                triQuadFaces = {
                    {min(triQuadCells[j].nodeIndices[0], triQuadCells[j].nodeIndices[1]), max(triQuadCells[j].nodeIndices[0], triQuadCells[j].nodeIndices[1])},
                    {min(triQuadCells[j].nodeIndices[1], triQuadCells[j].nodeIndices[2]), max(triQuadCells[j].nodeIndices[1], triQuadCells[j].nodeIndices[2])},
                    {min(triQuadCells[j].nodeIndices[2], triQuadCells[j].nodeIndices[0]), max(triQuadCells[j].nodeIndices[2], triQuadCells[j].nodeIndices[0])}};
            }
            else if (triQuadCells[j].nodeIndices.size() == 4) // Quadrilateral
            {
                triQuadFaces = {
                    {min(triQuadCells[j].nodeIndices[0], triQuadCells[j].nodeIndices[1]), max(triQuadCells[j].nodeIndices[0], triQuadCells[j].nodeIndices[1])},
                    {min(triQuadCells[j].nodeIndices[1], triQuadCells[j].nodeIndices[2]), max(triQuadCells[j].nodeIndices[1], triQuadCells[j].nodeIndices[2])},
                    {min(triQuadCells[j].nodeIndices[2], triQuadCells[j].nodeIndices[3]), max(triQuadCells[j].nodeIndices[2], triQuadCells[j].nodeIndices[3])},
                    {min(triQuadCells[j].nodeIndices[3], triQuadCells[j].nodeIndices[0]), max(triQuadCells[j].nodeIndices[3], triQuadCells[j].nodeIndices[0])}};
            }
            // Compare the faces of the boundary cell with the faces of the triQuadCells
            // can this be done using hash map
            for (size_t k = 0; k < faces.size(); k++)
            {
                for (size_t l = 0; l < triQuadFaces.size(); l++)
                {
                    if (faces[k] == triQuadFaces[l])
                    {
                        boundaryCellToParentCell[i] = j;
                    }
                }
            }
        }
    }
    // Print the boundary cell and its parent cell
    //    for (auto it = boundaryCellToParentCell.begin(); it != boundaryCellToParentCell.end(); ++it)
    //    {
    //        cout << "Boundary Cell " << it->first << " Parent Cell " << it->second << endl;
    //    }
}
// Function to generate face-to-cell mapping and detect neighbors
void mapFacesToCells(vector<Cell> &cells, map<pair<int, int>, set<int>> &faceToCells)
{
    bool restartNeeded = false;
    set<pair<int, int>> nonManifoldFaces; // Store faces that are shared by more than two cells

    for (size_t cellID = 0; cellID < cells.size(); ++cellID)
    {
        const auto &cellNodes = cells[cellID].nodeIndices;
        // cout << "Cell ID: " << cellID << "\t";
        //  For each cell go through its nodes and create a faces
        for (size_t i = 0; i < cellNodes.size(); ++i)
        {
            int node1 = cellNodes[i];
            int node2 = cellNodes[(i + 1) % cellNodes.size()];

            // Ensure consistent ordering of face
            if (node1 > node2)
                swap(node1, node2);
            pair<int, int> face = {node1, node2};

            // Print face nodes
            // cout << "Face Nodes: " << node1 << "\t" << node2 << "\t";
            // search face in the list of faces
            auto it = faceToCells.find(face);
            if (it != faceToCells.end())
            {
                // The face is already stored, meaning it's shared
                set<int> &sharedCells = it->second;
                sharedCells.insert(cellID);

                if (sharedCells.size() == 2)
                {
                    // The face is now shared between exactly two cells (valid)
                    auto it = sharedCells.begin();
                    int firstCell = *it++;
                    int secondCell = *it;

                    // ✅ Use push_back() instead of overwriting Neighbours[i]
                    if (find(cells[firstCell].Neighbours.begin(), cells[firstCell].Neighbours.end(), secondCell) == cells[firstCell].Neighbours.end())
                    {
                        cells[firstCell].Neighbours.push_back(secondCell);
                    }
                    if (find(cells[secondCell].Neighbours.begin(), cells[secondCell].Neighbours.end(), firstCell) == cells[secondCell].Neighbours.end())
                    {
                        cells[secondCell].Neighbours.push_back(firstCell);
                    }
                    // cout << "First Cell: " << firstCell << "\tSecond Cell: " << secondCell << endl;
                }
                else if (sharedCells.size() > 2)
                {
                    // 🚨 Face is shared by more than two cells (Non-Manifold case)
                    cout << "Face is shared by more than two cells: ";
                    for (int c : sharedCells)
                        cout << c << " ";
                    cout << endl;
                    nonManifoldFaces.insert(face);
                }
            }
            else
            {
                // If the face is new, add it to the map
                faceToCells[face] = {static_cast<int>(cellID)};
                // cout << "New Face: " << cellID << endl;
            }
        }
    }

    // 🚨 Log non-manifold faces instead of causing infinite loop
    if (!nonManifoldFaces.empty())
    {
        cout << "Warning: Non-manifold faces detected! These faces are shared by more than two cells:\n";
        for (auto &face : nonManifoldFaces)
        {
            cout << "Face (" << face.first << ", " << face.second << ") shared by multiple cells.\n";
        }
    }
}
// Function to print a Cell
void printCellEntry(const Cell &entry)
{
    cout << "Cell Type: " << entry.cellType << ", Cell ID: " << entry.cellID << "  ";

    // Print node indices
    cout << "Node Indices: ";
    for (int node : entry.nodeIndices)
    {
        cout << node << "\t ";
    }
    // cout << endl;
    cout << "Gmsh Indices are : ";
    for (int node : entry.nodeIndices)
    {
        cout << node + 1 << "\t";
    }
    // Print neighbors
    cout << "Neighbours: \t" << entry.Neighbours.size() << "\t";
    if (entry.Neighbours.empty())
    {
        cout << "None";
    }
    else
    {
        for (int neighbor : entry.Neighbours)
        {
            cout << neighbor << "\t";
        }
    }
    cout << "Face Normals are : ";
    for (int i = 0; i < entry.Face_Normals.size(); i++)
    {
        cout << entry.Face_Normals[i] << "\t";
    }
    cout << endl;
}

void BoundingBox(V_D &Points, double &Xmin, double &Xmax, double &Ymin, double &Ymax, double &Zmin, double &Zmax)
{
    // obtain the bounding box of the domain by looping over all the points
    for (size_t i = 0; i < Points.size(); i += 3)
    {
        if (i == 0)
        {
            Xmin = Points[i];
            Xmax = Points[i];
            Ymin = Points[i + 1];
            Ymax = Points[i + 1];
            Zmin = Points[i + 2];
            Zmax = Points[i + 2];
        }
        else
        {
            if (Points[i] < Xmin)
                Xmin = Points[i];
            if (Points[i] > Xmax)
                Xmax = Points[i];
            if (Points[i + 1] < Ymin)
                Ymin = Points[i + 1];
            if (Points[i + 1] > Ymax)
                Ymax = Points[i + 1];
            if (Points[i + 2] < Zmin)
                Zmin = Points[i + 2];
            if (Points[i + 2] > Zmax)
                Zmax = Points[i + 2];
        }
    }
    cout << "Xmin: " << Xmin << " Xmax: " << Xmax << " Ymin: " << Ymin << " Ymax: " << Ymax << " Zmin: " << Zmin << " Zmax: " << Zmax << endl;
}

void Read_GmshMESH_Grid(const string &GridFileName)
{
    // This function reads a native .msh file generated by Gmsh
    ifstream file(GridFileName);
    if (!file)
    {
        cerr << "Could not open the file! Check the file path: " << GridFileName << endl;
        exit(1);
    }
    Is_2D_Flow = true;
    string line;
    vector<double> Points;
    vector<Cell> Cells, Physical_Cells, Boundary_LineCells;
    int No_of_Points = 0, No_of_Cells = 0, NumEntries = 0;
}

std::string Get_Boundary_Type(V_D &p1, V_D &p2, double &minX, double &maxX, double &minY, double &maxY, bool isBoundaryFace)
{
    if (p1.empty() || p2.empty())
    {
        std::cerr << "Error: One of the points is empty." << std::endl;
        return "INVALID";
    }
    // Print(p1);
    // Print(p2);
    if (!isBoundaryFace)
        return "INTERNAL";
    if (std::abs(p1[0] - minX) < EPSILON && std::abs(p2[0] - minX) < EPSILON)
        return "LEFT";
    else if (std::abs(p1[0] - maxX) < EPSILON && std::abs(p2[0] - maxX) < EPSILON)
        return "RIGHT";
    else if (std::abs(p1[1] - minY) < EPSILON && std::abs(p2[1] - minY) < EPSILON)
        return "BOTTOM";
    else if (std::abs(p1[1] - maxY) < EPSILON && std::abs(p2[1] - maxY) < EPSILON)
        return "TOP";

    return "WALL"; // Internal wall or embedded obstacle
}

void Classify_Domain_Boundaries(std::vector<Cell> &Cells, double &minX, double &maxX, double &minY, double &maxY)
{
    V_D p1(3, 0.0), p2(3, 0.0);

    for (size_t i = 0; i < Cells.size(); i++)
    {
        if (!Cells[i].hasBoundaryface)
            continue;

        // cout << "Cell ID: " << Cells[i].cellID << "\t";
        // cout << "Has boundary face : " << Cells[i].hasBoundaryface << "\t";
        // cout << "No of Boundary Faces : " << Cells[i].NoBoundaryFaces << "\t";
        //  Reset all face flags
        Cells[i].Left_Face = false;
        Cells[i].Right_Face = false;
        Cells[i].Top_Face = false;
        Cells[i].Bottom_Face = false;
        Cells[i].Interior_Face = false;

        Neighbour_1 = Cells[i].Neighbours[0];
        Neighbour_2 = Cells[i].Neighbours[1];
        Neighbour_3 = Cells[i].Neighbours[2];
        Neighbour_4 = Cells[i].Neighbours[3];

        // Which Neighhour index is greater than the Number_Physical_Cells or Cells.size()
        // cout << Neighbour_1 << "\t" << Neighbour_2 << "\t" << Neighbour_3 << "\t" << Neighbour_4 << "\t";
        size_t numNodes = Cells[i].nodeIndices.size();
        size_t numEdges = numNodes;

        bool isDomainBoundary = false;

        for (size_t e = 0; e < numEdges; e++)
        {
            size_t idx1 = (e % numNodes) * 3;
            size_t idx2 = ((e + 1) % numNodes) * 3;

            // Extract p1
            p1[0] = Cells[i].Cell_Vertices[idx1 + 0];
            p1[1] = Cells[i].Cell_Vertices[idx1 + 1];
            p1[2] = Cells[i].Cell_Vertices[idx1 + 2];

            // Extract p2
            p2[0] = Cells[i].Cell_Vertices[idx2 + 0];
            p2[1] = Cells[i].Cell_Vertices[idx2 + 1];
            p2[2] = Cells[i].Cell_Vertices[idx2 + 2];

            // Get the boundary type of the face
            std::string btype = Get_Boundary_Type(p1, p2, minX, maxX, minY, maxY, Cells[i].hasBoundaryface);

            if (btype == "LEFT")
            {
                Cells[i].Left_Face = true;
                isDomainBoundary = true;
            }
            else if (btype == "RIGHT")
            {
                Cells[i].Right_Face = true;
                isDomainBoundary = true;
            }
            else if (btype == "TOP")
            {
                Cells[i].Top_Face = true;
                isDomainBoundary = true;
            }
            else if (btype == "BOTTOM")
            {
                Cells[i].Bottom_Face = true;
                isDomainBoundary = true;
            }
        }
        // If it's a boundary cell but no face lies on the domain boundary, it's an internal boundary (like an obstacle)
        if (!isDomainBoundary)
        {
            Cells[i].Interior_Face = true;
        }
        // cout << Cells[i].Left_Face << "\t" << Cells[i].Bottom_Face << "\t" << Cells[i].Right_Face << "\t" << Cells[i].Top_Face << endl;
    }
}

void Create_Boundary_Cells_Lists(vector<Cell> &Cells, vector<int> &Inlet_Cells_List, vector<int> &Exit_Cells_List, vector<int> &Wall_Cells_List)
{
    // Mapping Cells with faces Left - face no 0 bottom - face no 1
    // right - face no 2 top - face no 3
    // Ghost cells are numberd after No_Physical_Cells
    cout << "Creating Boundary Cells Lists" << endl;
    cout << "Size of Cells is : " << Cells.size() << endl;
    for (size_t i = 0; i < Cells.size(); i++)
    {
        /* cout << "Cell ID: " << Cells[i].cellID << "\t";
         cout << "Has boundary face : " << Cells[i].hasBoundaryface << "\t";
         cout << "No of Boundary Faces : " << Cells[i].NoBoundaryFaces << "\t";
         cout << Cells[i].Left_Face << endl;
         cout << Cells[i].Right_Face << endl;
         cout << Cells[i].Top_Face << endl;
         cout << Cells[i].Bottom_Face << endl;*/

        if (Cells[i].Left_Face)
        {
            Inlet_Cells_List.push_back(i);
            Inlet_Cells_List.push_back(0);
        }
        if (Cells[i].Right_Face)
        {
            Exit_Cells_List.push_back(i);
            Exit_Cells_List.push_back(2);
        }
        if (Cells[i].Top_Face)
        {
            Wall_Cells_List.push_back(i);
            Wall_Cells_List.push_back(3);
        }
        if (Cells[i].Bottom_Face)
        {
            Wall_Cells_List.push_back(i);
            Wall_Cells_List.push_back(1);
        }
        if (Cells[i].Interior_Face)
        {
            Wall_Cells_List.push_back(i);
            Wall_Cells_List.push_back(1);
        }
    }
    cout << "Inlet Boundary Cells: " << Inlet_Cells_List.size() << endl;
    cout << "Exit Boundary Cells: " << Exit_Cells_List.size() << endl;
    cout << "Wall Boundary Cells: " << Wall_Cells_List.size() << endl;
}

// Helper: Given a face's normal (nx, ny), returns a priority integer:
//  1 -> left, 2 -> bottom, 3 -> right, 4 -> top.
// The assignment is made by comparing the normal's angle with the target angles:
// left: π (or –π), bottom: -π/2, right: 0, top: π/2.
int compute_face_priority(double nx, double ny)
{
    // Compute the angle in radians; range is (-PI, PI]
    double angle = std::atan2(ny, nx);

    // For the left face, we want to consider angles near π and -π as equivalent.
    double diff_left = fabs(angle - M_PI);
    if (angle < 0)
        diff_left = std::min(diff_left, fabs(angle + M_PI));

    double diff_bottom = fabs(angle - (-M_PI / 2));
    double diff_right = fabs(angle - 0);
    double diff_top = fabs(angle - (M_PI / 2));

    // Start by assuming left has the smallest difference.
    int priority = 1;
    double minDiff = diff_left;
    if (diff_bottom < minDiff)
    {
        minDiff = diff_bottom;
        priority = 2;
    }
    if (diff_right < minDiff)
    {
        minDiff = diff_right;
        priority = 3;
    }
    if (diff_top < minDiff)
    { /* minDiff = diff_top; */
        priority = 4;
    }
    return priority;
}

// Sort_Neighbours reorders each cell's Neighbours vector so that the
// indices (or ghost marker -1) are arranged in the following order:
// left face, bottom face, right face, top face.
// It uses the cell's Face_Normals (assumed to be a vector of 2D vectors)
// to decide each face's orientation.
void Sort_Neighbours(vector<Cell> &Cells)
{
    for (size_t i = 0; i < Cells.size(); i++)
    {
        Cell &cell = Cells[i];
        // Assume that the number of faces is stored in cell.numFaces
        // (or, equivalently, cell.nodeIndices.size() for many cell types).
        int nFaces = cell.numFaces;

        // Build a vector of (faceIndex, priority) pairs.
        vector<pair<int, int>> facePriority;
        for (int f = 0; f < nFaces; f++)
        {
            // Get the face normal. We assume that Face_Normals[f]
            // is a vector<double> of at least 2 values: [nx, ny].
            double nx = cell.Face_Normals[f * 2 + 0];
            double ny = cell.Face_Normals[f * 2 + 1];
            int priority = compute_face_priority(nx, ny);
            facePriority.push_back(make_pair(f, priority));
        }

        // Sort faces based on the computed priority.
        std::sort(facePriority.begin(), facePriority.end(),
                  [](const pair<int, int> &a, const pair<int, int> &b)
                  {
                      return a.second < b.second;
                  });

        // Rebuild a sorted neighbour vector. The ghost cell value (-1)
        // stays attached to its face.
        vector<int> sortedNeighbours(nFaces, -1);
        for (size_t idx = 0; idx < facePriority.size(); idx++)
        {
            int origFaceIndex = facePriority[idx].first;
            sortedNeighbours[idx] = cell.Neighbours[origFaceIndex];
        }

        // Update the cell's neighbour vector.
        cell.Neighbours = sortedNeighbours;
    }
}
