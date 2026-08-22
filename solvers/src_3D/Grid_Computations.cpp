#include "definitions.h"
#include "Grid.h"
#include "Globals.h"
#include "Utilities.h"

/**
 * @file Grid_Computations.cpp
 * @brief 3D Grid computation functions for CFD solver
 *
 * This file contains functions for reading 3D grid files, constructing hexahedral cells,
 * calculating face areas and normals, cell volumes, and geometric properties for 3D meshes.
 */

vector<Cell> Cells, Boundary_Cells, Co_Volume_Cells;

V_I Wall_Cells_List, Inlet_Cells_List, Exit_Cells_List, Symmetry_Cells_List;
V_I Far_Field_Out_Flow_List, Far_Field_InFlow_List;
V_I Back_Cells_List, Front_Cells_List; // Added for 3D boundaries

V_D Vertices;
int Total_No_Cells, No_Cartesian_Cells, No_Polar_Cells, No_Physical_Cells;
int No_Ghost_Cells, Cells_in_Volume, nx_c, ny_c, nz_c, nx_p, ny_p, nz_p, Grid_Type;

double global_temp, R_Mid_dot_A, Cell_Minimum_Length;
bool Is_Viscous_Wall, Is_3D_Flow, Is_Inlet_SubSonic, Is_Exit_SubSonic, Enable_Far_Field, has_Symmetry_BC;
double CFL;

// 3D cell type definitions
int numNodes, nodeIndex, PointCellType, LineCellType, TriangleCellType, QuadrilateralCellType;
int HexahedronCellType, TetrahedronCellType, WedgeCellType, PyramidCellType, PrismCellType;

int get_NoPhysical_Cells()
{
    return Total_No_Cells;
}

/**
 * @brief Reads 3D input grid file and constructs hexahedral cells
 * @param ipfile The path to the input grid file
 *
 * Function for reading 3D grid file containing hexahedral cell information.
 * The grid file contains vertices of cells followed by connectivity information.
 * Creates 3D physical cells, evaluates face areas, normals, and cell volumes.
 */
void Read_Grid(const string &ipfile)
{
    // 3D points for hexahedral cell (8 vertices)
    V_D P1(3, 0.0), P2(3, 0.0), P3(3, 0.0), P4(3, 0.0);
    V_D P5(3, 0.0), P6(3, 0.0), P7(3, 0.0), P8(3, 0.0);

    Cell Grid_Cell = {};
    double x = 0.0, y = 0.0, z = 0.0;
    int cellID, neighbor1, neighbor2, neighbor3, neighbor4, neighbor5, neighbor6, Conversion_Type;

    // Resize vertices for hexahedral cell (8 vertices × 3 coordinates = 24)
    Vertices.resize(24, 0.0);

    cout << "Reading 3D grid file: " << ipfile.c_str() << endl;
    ifstream Grid_File(ipfile.c_str(), ios::in);

    if (Grid_File.is_open())
    {
        cout << "3D Grid file opened for reading" << endl;
        Grid_File >> Grid_Type;
        Grid_File >> Conversion_Type;

        if (Conversion_Type == 1)
            cout << "Grid is created in mm - Converting to meters" << endl;
        else
            cout << "Grid is created in meters" << endl;

        switch (Grid_Type)
        {
        case GRID_CARTESIAN_3D: // 3D Cartesian Mesh
            cout << "Reading 3D Cartesian Grid" << endl;
            Grid_File >> nx_c >> ny_c >> nz_c; // 3D grid dimensions
            Grid_File >> No_Physical_Cells;

            cout << "Grid dimensions: nx=" << nx_c << ", ny=" << ny_c << ", nz=" << nz_c << endl;
            cout << "Physical cells: " << No_Physical_Cells << endl;

            Cells_in_Volume = (nx_c - 1) * (ny_c - 1) * (nz_c - 1);

            // Calculate ghost cells for 3D (6 boundary faces)
            No_Ghost_Cells = 2 * ((nx_c - 1) * (ny_c - 1)) + // front + back faces
                             2 * ((nx_c - 1) * (nz_c - 1)) + // top + bottom faces
                             2 * ((ny_c - 1) * (nz_c - 1));  // left + right faces

            cout << "Total ghost cells to be constructed: " << No_Ghost_Cells << endl;
            Total_No_Cells = No_Physical_Cells + No_Ghost_Cells;
            cout << "Total cells: " << Total_No_Cells << endl;
            break;

        case GRID_MULTI_BLOCK_3D: // 3D Multi-block mesh
            cout << "Reading 3D Multi-block Grid" << endl;
            Grid_File >> nx_1 >> ny_1 >> nz_1 >> nx_2 >> ny_2 >> nz_2;
            Grid_File >> No_Physical_Cells;

            cout << "Block 1: " << nx_1 << "x" << ny_1 << "x" << nz_1 << endl;
            cout << "Block 2: " << nx_2 << "x" << ny_2 << "x" << nz_2 << endl;

            Cells_in_Volume = (nx_1 - 1) * (ny_1 - 1) * (nz_1 - 1) +
                              (nx_2 - 1) * (ny_2 - 1) * (nz_2 - 1);

            // Calculate ghost cells for multi-block
            No_Ghost_Cells = 2 * ((nx_1 - 1) * (ny_1 - 1)) + 2 * ((nx_1 - 1) * (nz_1 - 1)) + 2 * ((ny_1 - 1) * (nz_1 - 1)) +
                             2 * ((nx_2 - 1) * (ny_2 - 1)) + 2 * ((nx_2 - 1) * (nz_2 - 1)) + 2 * ((ny_2 - 1) * (nz_2 - 1));

            Total_No_Cells = No_Physical_Cells + No_Ghost_Cells;
            break;

        default:
            cout << "Unknown 3D grid type: " << Grid_Type << endl;
            exit(1);
        }

        // Reading physical cell information (8 vertices per hexahedral cell)
        for (int i = 0; i < No_Physical_Cells; i++)
        {
            // Read 8 vertices of hexahedral cell
            Grid_File >> x >> y >> z;
            P1[0] = x;
            P1[1] = y;
            P1[2] = z;
            Grid_File >> x >> y >> z;
            P2[0] = x;
            P2[1] = y;
            P2[2] = z;
            Grid_File >> x >> y >> z;
            P3[0] = x;
            P3[1] = y;
            P3[2] = z;
            Grid_File >> x >> y >> z;
            P4[0] = x;
            P4[1] = y;
            P4[2] = z;
            Grid_File >> x >> y >> z;
            P5[0] = x;
            P5[1] = y;
            P5[2] = z;
            Grid_File >> x >> y >> z;
            P6[0] = x;
            P6[1] = y;
            P6[2] = z;
            Grid_File >> x >> y >> z;
            P7[0] = x;
            P7[1] = y;
            P7[2] = z;
            Grid_File >> x >> y >> z;
            P8[0] = x;
            P8[1] = y;
            P8[2] = z;

            // Apply unit conversion if needed
            if (Conversion_Type == 1)
            {
                Conversion_Factor(P1);
                Conversion_Factor(P2);
                Conversion_Factor(P3);
                Conversion_Factor(P4);
                Conversion_Factor(P5);
                Conversion_Factor(P6);
                Conversion_Factor(P7);
                Conversion_Factor(P8);
            }

            // Store vertices of current hexahedral cell (8 vertices × 3 coordinates)
            Vertices[0] = P1[0];
            Vertices[1] = P1[1];
            Vertices[2] = P1[2];
            Vertices[3] = P2[0];
            Vertices[4] = P2[1];
            Vertices[5] = P2[2];
            Vertices[6] = P3[0];
            Vertices[7] = P3[1];
            Vertices[8] = P3[2];
            Vertices[9] = P4[0];
            Vertices[10] = P4[1];
            Vertices[11] = P4[2];
            Vertices[12] = P5[0];
            Vertices[13] = P5[1];
            Vertices[14] = P5[2];
            Vertices[15] = P6[0];
            Vertices[16] = P6[1];
            Vertices[17] = P6[2];
            Vertices[18] = P7[0];
            Vertices[19] = P7[1];
            Vertices[20] = P7[2];
            Vertices[21] = P8[0];
            Vertices[22] = P8[1];
            Vertices[23] = P8[2];

            // Initialize grid cell
            Grid_Cell.cellID = i;
            Grid_Cell.Dimension = 3;
            Grid_Cell.numFaces = NUM_FACES_3D;  // 6 faces for hexahedron
            Grid_Cell.numNodes = NUM_NODES_HEX; // 8 nodes for hexahedron

            if (Grid_Cell.Cell_Vertices.empty())
                Grid_Cell.Cell_Vertices.resize(24, 0.0);
            Grid_Cell.Cell_Vertices = Vertices;

            if (Grid_Cell.Neighbours.empty())
                Grid_Cell.Neighbours.resize(NUM_FACES_3D, 0);

            // Read neighboring cell information (6 neighbors for 3D)
            Grid_File >> cellID >> neighbor1 >> neighbor2 >> neighbor3 >> neighbor4 >> neighbor5 >> neighbor6;
            Grid_Cell.Neighbours = {neighbor1, neighbor2, neighbor3,
                                    neighbor4, neighbor5, neighbor6};

            Grid_Cell.faceID.resize(NUM_FACES_3D, 0);
            Grid_Cell.nodeIndices.resize(NUM_NODES_HEX, 0);
            Grid_Cell.Volume = 0.0;
            Grid_Cell.Inv_Volume = 0.0;

            // Set node indices for hexahedral cell
            Grid_Cell.nodeIndices = {0, 1, 2, 3, 4, 5, 6, 7};

            // Ensure proper vertex ordering for 3D hexahedron
            Sort_Points_3D_Hexahedron(Grid_Cell.Cell_Vertices, Grid_Cell.nodeIndices);

            // Construct 3D cell geometry
            Construct_Cell_3D(Grid_Cell);

            Cells.push_back(Grid_Cell);

            // Clear cell data
            Grid_Cell.Cell_Vertices.clear();
            Grid_Cell.Neighbours.clear();
            Grid_Cell.faceID.clear();
            Grid_Cell.nodeIndices.clear();
            Grid_Cell.Face_Areas.clear();
            Grid_Cell.Face_Normals.clear();
            Grid_Cell.Face_Centers.clear();
        }

        cout << "Construction of 3D Physical Cells Completed" << endl;
        Check_Cells_3D();
        cout << "3D Cell checking completed" << endl;

        // Read boundary information for 3D
        Read_3D_Boundary_Information(Grid_File);

        Grid_File.close();
    }
    else
    {
        cout << "Error: Cannot open 3D grid file " << ipfile << endl;
        exit(1);
    }
}

/**
 * @brief Constructs 3D hexahedral cell geometry
 * @param Grid_Cell Reference to the cell being constructed
 *
 * Calculates cell center, face areas and normals, and cell volume for 3D hexahedral cells.
 */
void Construct_Cell_3D(Cell &Grid_Cell)
{
    V_D Temp(3, 0.0);
    double Volume = 0.0;

    // Initialize cell center
    Grid_Cell.Cell_Center.resize(3, 0);
    Grid_Cell.Cell_Center = Temp;

    // Compute 3D centroid
    Compute_Centroid_3D(Grid_Cell);

    // Construct face areas and normals for 6 faces
    Construct_Face_3D(Grid_Cell);

    // Calculate cell volume using divergence theorem
    Volume = Calculate_Hexahedron_Volume(Grid_Cell.Cell_Vertices);
    Grid_Cell.Volume = Volume;
    Grid_Cell.Inv_Volume = 1.0 / Volume;

    // Calculate face distances from cell center
    Calculate_Face_Distances_3D(Grid_Cell);
}

/**
 * @brief Constructs face geometry for 3D hexahedral cell
 * @param Grid_Cell Reference to the cell
 *
 * Calculates area, normal, and center for all 6 faces of hexahedral cell.
 */
void Construct_Face_3D(Cell &Grid_Cell)
{
    Grid_Cell.Face_Areas.resize(NUM_FACES_3D, 0.0);
    Grid_Cell.Face_Normals.resize(NUM_FACES_3D * 3, 0.0); // 3 components per face
    Grid_Cell.Face_Centers.resize(NUM_FACES_3D * 3, 0.0);

    // Face vertex indices for hexahedral cell (following standard convention)
    vector<vector<int>> face_vertices = {
        {0, 3, 7, 4}, // Face 0: Left face   (x-min)
        {1, 2, 6, 5}, // Face 1: Right face  (x-max)
        {0, 1, 5, 4}, // Face 2: Bottom face (y-min)
        {2, 3, 7, 6}, // Face 3: Top face    (y-max)
        {0, 1, 2, 3}, // Face 4: Back face   (z-min)
        {4, 5, 6, 7}  // Face 5: Front face  (z-max)
    };

    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Get vertices of current face
        V_D face_coords(12, 0.0); // 4 vertices × 3 coordinates
        for (int i = 0; i < 4; i++)
        {
            int vertex_idx = face_vertices[face][i];
            face_coords[3 * i] = Grid_Cell.Cell_Vertices[3 * vertex_idx];
            face_coords[3 * i + 1] = Grid_Cell.Cell_Vertices[3 * vertex_idx + 1];
            face_coords[3 * i + 2] = Grid_Cell.Cell_Vertices[3 * vertex_idx + 2];
        }

        // Calculate face area using cross product of diagonals
        V_D diagonal1(3), diagonal2(3), normal(3);
        diagonal1[0] = face_coords[6] - face_coords[0]; // P3 - P1
        diagonal1[1] = face_coords[7] - face_coords[1];
        diagonal1[2] = face_coords[8] - face_coords[2];

        diagonal2[0] = face_coords[9] - face_coords[3]; // P4 - P2
        diagonal2[1] = face_coords[10] - face_coords[4];
        diagonal2[2] = face_coords[11] - face_coords[5];

        // Cross product for normal
        CROSS_PRODUCT_3D(normal, diagonal1, diagonal2);
        double area = 0.5 * MAGNITUDE_3D(normal);

        // Normalize to get unit normal
        NORMALIZE_3D(normal);

        // Store face area and normal
        Grid_Cell.Face_Areas[face] = area;
        Grid_Cell.Face_Normals[3 * face] = normal[0];
        Grid_Cell.Face_Normals[3 * face + 1] = normal[1];
        Grid_Cell.Face_Normals[3 * face + 2] = normal[2];

        // Calculate face center
        double face_center[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < 4; i++)
        {
            face_center[0] += face_coords[3 * i];
            face_center[1] += face_coords[3 * i + 1];
            face_center[2] += face_coords[3 * i + 2];
        }
        Grid_Cell.Face_Centers[3 * face] = face_center[0] / 4.0;
        Grid_Cell.Face_Centers[3 * face + 1] = face_center[1] / 4.0;
        Grid_Cell.Face_Centers[3 * face + 2] = face_center[2] / 4.0;
    }
}

/**
 * @brief Computes centroid of 3D hexahedral cell
 * @param Grid_Cell Reference to the cell
 */
void Compute_Centroid_3D(Cell &Grid_Cell)
{
    Grid_Cell.Cell_Center[0] = 0.0;
    Grid_Cell.Cell_Center[1] = 0.0;
    Grid_Cell.Cell_Center[2] = 0.0;

    // Average all 8 vertices
    for (int i = 0; i < NUM_NODES_HEX; i++)
    {
        Grid_Cell.Cell_Center[0] += Grid_Cell.Cell_Vertices[3 * i];
        Grid_Cell.Cell_Center[1] += Grid_Cell.Cell_Vertices[3 * i + 1];
        Grid_Cell.Cell_Center[2] += Grid_Cell.Cell_Vertices[3 * i + 2];
    }

    Grid_Cell.Cell_Center[0] /= NUM_NODES_HEX;
    Grid_Cell.Cell_Center[1] /= NUM_NODES_HEX;
    Grid_Cell.Cell_Center[2] /= NUM_NODES_HEX;
}

/**
 * @brief Calculates volume of hexahedral cell
 * @param vertices Vector containing 8 vertices (24 coordinates)
 * @return Volume of the hexahedron
 */
double Calculate_Hexahedron_Volume(const V_D &vertices)
{
    // Decompose hexahedron into 6 tetrahedra and sum their volumes
    double total_volume = 0.0;

    // Define tetrahedron decomposition (vertex indices)
    vector<vector<int>> tetrahedra = {
        {0, 1, 2, 5}, {0, 2, 3, 7}, {0, 5, 7, 4}, {2, 5, 6, 7}, {0, 2, 5, 7}, {5, 6, 7, 2}};

    for (const auto &tet : tetrahedra)
    {
        // Get tetrahedron vertices
        V_D v0 = {vertices[3 * tet[0]], vertices[3 * tet[0] + 1], vertices[3 * tet[0] + 2]};
        V_D v1 = {vertices[3 * tet[1]], vertices[3 * tet[1] + 1], vertices[3 * tet[1] + 2]};
        V_D v2 = {vertices[3 * tet[2]], vertices[3 * tet[2] + 1], vertices[3 * tet[2] + 2]};
        V_D v3 = {vertices[3 * tet[3]], vertices[3 * tet[3] + 1], vertices[3 * tet[3] + 2]};

        // Calculate tetrahedron volume using scalar triple product
        V_D a = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        V_D b = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
        V_D c = {v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]};

        V_D cross_bc(3);
        CROSS_PRODUCT_3D(cross_bc, b, c);
        double tet_volume = fabs(DOT_PRODUCT_3D(a, cross_bc)) / 6.0;

        total_volume += tet_volume;
    }

    return total_volume;
}

/**
 * @brief Calculates distances from cell center to each face
 * @param Grid_Cell Reference to the cell
 */
void Calculate_Face_Distances_3D(Cell &Grid_Cell)
{
    Grid_Cell.Cell_Face_Distances.resize(NUM_FACES_3D, 0.0);

    for (int face = 0; face < NUM_FACES_3D; face++)
    {
        // Distance from cell center to face center
        double dx = Grid_Cell.Face_Centers[3 * face] - Grid_Cell.Cell_Center[0];
        double dy = Grid_Cell.Face_Centers[3 * face + 1] - Grid_Cell.Cell_Center[1];
        double dz = Grid_Cell.Face_Centers[3 * face + 2] - Grid_Cell.Cell_Center[2];

        Grid_Cell.Cell_Face_Distances[face] = sqrt(dx * dx + dy * dy + dz * dz);
    }
}

/**
 * @brief Reads 3D boundary information from grid file
 * @param Grid_File Reference to the input file stream
 */
void Read_3D_Boundary_Information(ifstream &Grid_File)
{
    cout << "Reading 3D Boundary Conditions Information" << endl;

    string boundary_name;
    int boundary_type, num_cells;

    // Read all 6 boundary types for 3D
    const vector<string> boundary_names = {"Left", "Right", "Bottom", "Top", "Back", "Front"};

    for (int i = 0; i < 6; i++)
    {
        Grid_File >> boundary_name >> boundary_type >> num_cells;
        cout << boundary_name << " boundary: type=" << boundary_type
             << ", cells=" << num_cells << endl;

        vector<int> boundary_cells;
        for (int j = 0; j < num_cells; j++)
        {
            int cell_id, face_id, ghost_id;
            Grid_File >> cell_id >> face_id >> ghost_id;
            boundary_cells.push_back(cell_id);
            boundary_cells.push_back(face_id);
            boundary_cells.push_back(ghost_id);
        }

        // Store boundary cells based on type
        switch (boundary_type)
        {
        case BC_INLET:
            Inlet_Cells_List.insert(Inlet_Cells_List.end(), boundary_cells.begin(), boundary_cells.end());
            break;
        case BC_OUTLET:
            Exit_Cells_List.insert(Exit_Cells_List.end(), boundary_cells.begin(), boundary_cells.end());
            break;
        case BC_WALL:
            Wall_Cells_List.insert(Wall_Cells_List.end(), boundary_cells.begin(), boundary_cells.end());
            break;
        case BC_SYMMETRY:
            Symmetry_Cells_List.insert(Symmetry_Cells_List.end(), boundary_cells.begin(), boundary_cells.end());
            break;
        case BC_FAR_FIELD:
            Far_Field_Out_Flow_List.insert(Far_Field_Out_Flow_List.end(), boundary_cells.begin(), boundary_cells.end());
            break;
        default:
            cout << "Unknown boundary type: " << boundary_type << endl;
        }
    }
}

/**
 * @brief Sorts vertices of hexahedral cell for proper orientation
 * @param vertices Vector of vertex coordinates
 * @param indices Vector of vertex indices
 */
void Sort_Points_3D_Hexahedron(V_D &vertices, V_I &indices)
{
    // For hexahedral cells, ensure vertices follow right-hand rule
    // Standard hexahedron vertex ordering:
    // Bottom face: 0-1-2-3 (counter-clockwise when viewed from top)
    // Top face: 4-5-6-7 (counter-clockwise when viewed from top)

    // This is a simplified version - in practice, you'd implement
    // proper geometric sorting based on connectivity

    // For now, assume input vertices are correctly ordered
    for (int i = 0; i < NUM_NODES_HEX; i++)
    {
        indices[i] = i;
    }
}

/**
 * @brief Checks validity of 3D cells
 */
void Check_Cells_3D()
{
    cout << "Checking 3D cell validity..." << endl;

    for (int i = 0; i < No_Physical_Cells; i++)
    {
        // Check positive volume
        if (Cells[i].Volume <= 0.0)
        {
            cout << "Error: Negative or zero volume in cell " << i
                 << ", Volume = " << Cells[i].Volume << endl;
        }

        // Check face areas
        for (int face = 0; face < NUM_FACES_3D; face++)
        {
            if (Cells[i].Face_Areas[face] <= 0.0)
            {
                cout << "Error: Negative or zero face area in cell " << i
                     << ", face " << face << endl;
            }
        }

        // Update minimum cell length for CFL condition
        double min_face_distance = *min_element(Cells[i].Cell_Face_Distances.begin(),
                                                Cells[i].Cell_Face_Distances.end());
        Cell_Minimum_Length = min(Cell_Minimum_Length, min_face_distance);
    }

    cout << "Minimum cell characteristic length: " << Cell_Minimum_Length << endl;
    cout << "3D Cell validation completed" << endl;
}

/**
 * @brief Unit conversion factor for coordinates
 * @param point Vector containing 3D coordinates
 */
void Conversion_Factor(V_D &point)
{
    const double mm_to_m = 0.001;
    point[0] *= mm_to_m;
    point[1] *= mm_to_m;
    point[2] *= mm_to_m;
}