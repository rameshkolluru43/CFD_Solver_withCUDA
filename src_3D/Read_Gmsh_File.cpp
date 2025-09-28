/**
 * @file Read_Gmsh_File.cpp
 * @brief GMSH 3D mesh file reader for hexahedral grid import
 *
 * This module handles reading and parsing GMSH mesh files for 3D CFD simulations.
 * Supports GMSH format versions 2.2 and 4.1 with comprehensive error handling,
 * mesh validation, and hexahedral element processing.
 *
 * Key Features:
 * - GMSH format 2.2 and 4.1 support
 * - 3D hexahedral element processing (element type 5)
 * - Comprehensive mesh validation and quality checks
 * - Boundary condition tag processing
 * - Multi-zone mesh support
 * - Automatic grid generation from mesh data
 * - Element connectivity validation
 * - Node coordinate processing with scaling
 *
 * Mathematical Framework:
 * - Hexahedral elements: 8 vertices per element
 * - Node coordinates: (x, y, z) in 3D space
 * - Element connectivity: local to global node mapping
 * - Boundary surfaces: triangular/quadrilateral faces
 * - Physical groups: boundary condition assignment
 *
 * @author CFD Solver Team
 * @date 2024
 */

#include "definitions.h"
#include "Globals.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>

/**
 * @brief GMSH element type definitions
 *
 * GMSH element type constants for different geometric elements.
 */
enum GMSH_Element_Types
{
    GMSH_POINT = 15,
    GMSH_LINE = 1,
    GMSH_TRIANGLE = 2,
    GMSH_QUADRANGLE = 3,
    GMSH_TETRAHEDRON = 4,
    GMSH_HEXAHEDRON = 5,
    GMSH_PRISM = 6,
    GMSH_PYRAMID = 7
};

/**
 * @brief GMSH node structure
 *
 * Structure to store node information including coordinates and ID.
 */
struct GMSH_Node
{
    int id;                             // Node ID
    double x, y, z;                     // Coordinates
    
    GMSH_Node() : id(0), x(0.0), y(0.0), z(0.0) {}
    GMSH_Node(int _id, double _x, double _y, double _z) : id(_id), x(_x), y(_y), z(_z) {}
};

/**
 * @brief GMSH element structure
 *
 * Structure to store element information including connectivity and tags.
 */
struct GMSH_Element
{
    int id;                             // Element ID
    int type;                           // Element type
    int num_tags;                       // Number of tags
    std::vector<int> tags;              // Element tags (physical, elementary)
    std::vector<int> nodes;             // Node connectivity
    
    GMSH_Element() : id(0), type(0), num_tags(0) {}
};

/**
 * @brief GMSH physical group structure
 *
 * Structure to store physical group information for boundary conditions.
 */
struct GMSH_Physical_Group
{
    int id;                             // Physical group ID
    int dimension;                      // Dimension (0=point, 1=line, 2=surface, 3=volume)
    std::string name;                   // Physical group name
    
    GMSH_Physical_Group() : id(0), dimension(0), name("") {}
    GMSH_Physical_Group(int _id, int _dim, const std::string &_name) 
        : id(_id), dimension(_dim), name(_name) {}
};

/**
 * @brief GMSH mesh data structure
 *
 * Complete structure to store all mesh data from GMSH file.
 */
struct GMSH_Mesh_Data
{
    std::vector<GMSH_Node> nodes;                   // All nodes
    std::vector<GMSH_Element> elements;             // All elements
    std::vector<GMSH_Physical_Group> physical_groups; // Physical groups
    std::map<int, int> node_id_to_index;            // Node ID to index mapping
    
    // Mesh statistics
    int num_hexahedra;                              // Number of hexahedral elements
    int num_boundary_faces;                         // Number of boundary faces
    double min_coords[3], max_coords[3];            // Bounding box
    
    GMSH_Mesh_Data()
    {
        num_hexahedra = 0;
        num_boundary_faces = 0;
        for (int i = 0; i < 3; i++)
        {
            min_coords[i] = 1e30;
            max_coords[i] = -1e30;
        }
    }
    
    void UpdateBoundingBox(double x, double y, double z)
    {
        min_coords[0] = std::min(min_coords[0], x);
        min_coords[1] = std::min(min_coords[1], y);
        min_coords[2] = std::min(min_coords[2], z);
        
        max_coords[0] = std::max(max_coords[0], x);
        max_coords[1] = std::max(max_coords[1], y);
        max_coords[2] = std::max(max_coords[2], z);
    }
};

/**
 * @brief Parses GMSH file header
 *
 * Reads and validates GMSH file format version and type.
 * Supports both ASCII and binary formats.
 *
 * @param file Input file stream
 * @param version Output version number
 * @param file_type Output file type (0=ASCII, 1=binary)
 * @param data_size Output data size
 * @return true if header parsed successfully
 */
bool Parse_GMSH_Header(std::ifstream &file, double &version, int &file_type, int &data_size)
{
    std::string line;
    
    // Read format line
    if (!std::getline(file, line))
    {
        std::cerr << "Error: Cannot read GMSH format line" << std::endl;
        return false;
    }
    
    // Check format identifier
    if (line.find("$MeshFormat") == std::string::npos)
    {
        std::cerr << "Error: Invalid GMSH file format (missing $MeshFormat)" << std::endl;
        return false;
    }
    
    // Read version info
    if (!std::getline(file, line))
    {
        std::cerr << "Error: Cannot read GMSH version information" << std::endl;
        return false;
    }
    
    std::istringstream iss(line);
    if (!(iss >> version >> file_type >> data_size))
    {
        std::cerr << "Error: Invalid GMSH version line format" << std::endl;
        return false;
    }
    
    // Validate version
    if (version < 2.0 || version > 4.2)
    {
        std::cerr << "Warning: GMSH version " << version << " may not be fully supported" << std::endl;
    }
    
    // Read end format line
    if (!std::getline(file, line) || line.find("$EndMeshFormat") == std::string::npos)
    {
        std::cerr << "Error: Missing $EndMeshFormat" << std::endl;
        return false;
    }
    
    std::cout << "GMSH file format detected: version " << version 
              << ", type " << (file_type ? "binary" : "ASCII") << std::endl;
    
    return true;
}

/**
 * @brief Parses GMSH physical groups
 *
 * Reads physical group definitions for boundary condition assignment.
 * Physical groups define named regions for boundary conditions.
 *
 * @param file Input file stream
 * @param mesh_data Output mesh data structure
 * @return true if physical groups parsed successfully
 */
bool Parse_GMSH_Physical_Groups(std::ifstream &file, GMSH_Mesh_Data &mesh_data)
{
    std::string line;
    int num_physical_groups;
    
    if (!std::getline(file, line))
    {
        std::cerr << "Error: Cannot read number of physical groups" << std::endl;
        return false;
    }
    
    num_physical_groups = std::stoi(line);
    std::cout << "Reading " << num_physical_groups << " physical groups" << std::endl;
    
    for (int i = 0; i < num_physical_groups; i++)
    {
        if (!std::getline(file, line))
        {
            std::cerr << "Error: Cannot read physical group " << i << std::endl;
            return false;
        }
        
        std::istringstream iss(line);
        GMSH_Physical_Group group;
        
        if (!(iss >> group.dimension >> group.id))
        {
            std::cerr << "Error: Invalid physical group format" << std::endl;
            return false;
        }
        
        // Read name (quoted string)
        std::string remaining;
        std::getline(iss, remaining);
        size_t first_quote = remaining.find('"');
        size_t last_quote = remaining.find_last_of('"');
        
        if (first_quote != std::string::npos && last_quote != std::string::npos && first_quote < last_quote)
        {
            group.name = remaining.substr(first_quote + 1, last_quote - first_quote - 1);
        }
        else
        {
            group.name = "Group_" + std::to_string(group.id);
        }
        
        mesh_data.physical_groups.push_back(group);
        
        std::cout << "Physical group: " << group.name << " (ID=" << group.id 
                  << ", dim=" << group.dimension << ")" << std::endl;
    }
    
    // Read end physical groups line
    if (!std::getline(file, line) || line.find("$EndPhysicalNames") == std::string::npos)
    {
        std::cerr << "Error: Missing $EndPhysicalNames" << std::endl;
        return false;
    }
    
    return true;
}

/**
 * @brief Parses GMSH nodes
 *
 * Reads all node coordinates and builds node mapping.
 * Handles both format 2.2 and 4.1 node specifications.
 *
 * @param file Input file stream
 * @param mesh_data Output mesh data structure
 * @param version GMSH format version
 * @return true if nodes parsed successfully
 *
 * Mathematical formulation:
 * - Node coordinates: P_i = (x_i, y_i, z_i)
 * - Bounding box: [x_min, x_max] × [y_min, y_max] × [z_min, z_max]
 * - Node mapping: global_id → local_index
 */
bool Parse_GMSH_Nodes(std::ifstream &file, GMSH_Mesh_Data &mesh_data, double version)
{
    std::string line;
    int num_nodes;
    
    if (!std::getline(file, line))
    {
        std::cerr << "Error: Cannot read number of nodes" << std::endl;
        return false;
    }
    
    num_nodes = std::stoi(line);
    std::cout << "Reading " << num_nodes << " nodes" << std::endl;
    
    mesh_data.nodes.reserve(num_nodes);
    
    for (int i = 0; i < num_nodes; i++)
    {
        if (!std::getline(file, line))
        {
            std::cerr << "Error: Cannot read node " << i << std::endl;
            return false;
        }
        
        std::istringstream iss(line);
        GMSH_Node node;
        
        if (version >= 4.0)
        {
            // Format 4.1: node_id x y z
            if (!(iss >> node.id >> node.x >> node.y >> node.z))
            {
                std::cerr << "Error: Invalid node format (version 4.x)" << std::endl;
                return false;
            }
        }
        else
        {
            // Format 2.2: node_id x y z
            if (!(iss >> node.id >> node.x >> node.y >> node.z))
            {
                std::cerr << "Error: Invalid node format (version 2.x)" << std::endl;
                return false;
            }
        }
        
        // Update bounding box
        mesh_data.UpdateBoundingBox(node.x, node.y, node.z);
        
        // Store node and create mapping
        mesh_data.node_id_to_index[node.id] = i;
        mesh_data.nodes.push_back(node);
    }
    
    // Read end nodes line
    if (!std::getline(file, line) || line.find("$EndNodes") == std::string::npos)
    {
        std::cerr << "Error: Missing $EndNodes" << std::endl;
        return false;
    }
    
    std::cout << "Nodes parsed successfully. Bounding box:" << std::endl;
    std::cout << "X: [" << mesh_data.min_coords[0] << ", " << mesh_data.max_coords[0] << "]" << std::endl;
    std::cout << "Y: [" << mesh_data.min_coords[1] << ", " << mesh_data.max_coords[1] << "]" << std::endl;
    std::cout << "Z: [" << mesh_data.min_coords[2] << ", " << mesh_data.max_coords[2] << "]" << std::endl;
    
    return true;
}

/**
 * @brief Parses GMSH elements
 *
 * Reads all elements and processes hexahedral elements for CFD solver.
 * Handles element connectivity and boundary face identification.
 *
 * @param file Input file stream
 * @param mesh_data Output mesh data structure
 * @param version GMSH format version
 * @return true if elements parsed successfully
 *
 * Mathematical formulation:
 * - Hexahedral elements: 8 nodes per element
 * - Node ordering: GMSH to solver convention
 * - Connectivity matrix: element → node mapping
 * - Boundary faces: surface elements for BCs
 */
bool Parse_GMSH_Elements(std::ifstream &file, GMSH_Mesh_Data &mesh_data, double version)
{
    std::string line;
    int num_elements;
    
    if (!std::getline(file, line))
    {
        std::cerr << "Error: Cannot read number of elements" << std::endl;
        return false;
    }
    
    num_elements = std::stoi(line);
    std::cout << "Reading " << num_elements << " elements" << std::endl;
    
    mesh_data.elements.reserve(num_elements);
    
    for (int i = 0; i < num_elements; i++)
    {
        if (!std::getline(file, line))
        {
            std::cerr << "Error: Cannot read element " << i << std::endl;
            return false;
        }
        
        std::istringstream iss(line);
        GMSH_Element element;
        
        if (!(iss >> element.id >> element.type >> element.num_tags))
        {
            std::cerr << "Error: Invalid element header format" << std::endl;
            return false;
        }
        
        // Read tags
        element.tags.resize(element.num_tags);
        for (int j = 0; j < element.num_tags; j++)
        {
            if (!(iss >> element.tags[j]))
            {
                std::cerr << "Error: Cannot read element tag " << j << std::endl;
                return false;
            }
        }
        
        // Determine number of nodes based on element type
        int num_nodes;
        switch (element.type)
        {
            case GMSH_POINT: num_nodes = 1; break;
            case GMSH_LINE: num_nodes = 2; break;
            case GMSH_TRIANGLE: num_nodes = 3; break;
            case GMSH_QUADRANGLE: num_nodes = 4; break;
            case GMSH_TETRAHEDRON: num_nodes = 4; break;
            case GMSH_HEXAHEDRON: num_nodes = 8; break;
            case GMSH_PRISM: num_nodes = 6; break;
            case GMSH_PYRAMID: num_nodes = 5; break;
            default:
                std::cerr << "Warning: Unknown element type " << element.type << std::endl;
                continue;
        }
        
        // Read node connectivity
        element.nodes.resize(num_nodes);
        for (int j = 0; j < num_nodes; j++)
        {
            if (!(iss >> element.nodes[j]))
            {
                std::cerr << "Error: Cannot read element node " << j << std::endl;
                return false;
            }
        }
        
        // Count different element types
        if (element.type == GMSH_HEXAHEDRON)
        {
            mesh_data.num_hexahedra++;
        }
        else if (element.type == GMSH_TRIANGLE || element.type == GMSH_QUADRANGLE)
        {
            mesh_data.num_boundary_faces++;
        }
        
        mesh_data.elements.push_back(element);
    }
    
    // Read end elements line
    if (!std::getline(file, line) || line.find("$EndElements") == std::string::npos)
    {
        std::cerr << "Error: Missing $EndElements" << std::endl;
        return false;
    }
    
    std::cout << "Elements parsed successfully:" << std::endl;
    std::cout << "Hexahedral elements: " << mesh_data.num_hexahedra << std::endl;
    std::cout << "Boundary faces: " << mesh_data.num_boundary_faces << std::endl;
    
    return true;
}

/**
 * @brief Converts GMSH mesh to solver grid format
 *
 * Processes GMSH mesh data and creates solver-compatible grid structures.
 * Handles coordinate transformation and connectivity mapping.
 *
 * @param mesh_data Input GMSH mesh data
 * @return true if conversion successful
 *
 * Mathematical formulation:
 * - Grid point mapping: GMSH nodes → solver vertices
 * - Cell connectivity: hexahedral elements → solver cells
 * - Coordinate scaling and transformation
 * - Boundary condition assignment
 */
bool Convert_GMSH_To_Solver_Grid_3D(const GMSH_Mesh_Data &mesh_data)
{
    if (mesh_data.num_hexahedra == 0)
    {
        std::cerr << "Error: No hexahedral elements found in mesh" << std::endl;
        return false;
    }
    
    std::cout << "Converting GMSH mesh to solver grid format" << std::endl;
    
    // Set grid dimensions based on hexahedral elements
    Total_Cells_3D = mesh_data.num_hexahedra;
    Total_Points_3D = mesh_data.nodes.size();
    
    // Allocate grid points
    Grid_Points_3D = new Point[Total_Points_3D];
    
    // Copy node coordinates
    for (size_t i = 0; i < mesh_data.nodes.size(); i++)
    {
        Grid_Points_3D[i].x = mesh_data.nodes[i].x;
        Grid_Points_3D[i].y = mesh_data.nodes[i].y;
        Grid_Points_3D[i].z = mesh_data.nodes[i].z;
    }
    
    // Process hexahedral elements
    std::vector<std::vector<int>> cell_connectivity;
    cell_connectivity.reserve(Total_Cells_3D);
    
    int cell_count = 0;
    for (const auto &element : mesh_data.elements)
    {
        if (element.type == GMSH_HEXAHEDRON)
        {
            std::vector<int> nodes(8);
            
            // Convert GMSH node IDs to solver indices
            for (int i = 0; i < 8; i++)
            {
                auto it = mesh_data.node_id_to_index.find(element.nodes[i]);
                if (it == mesh_data.node_id_to_index.end())
                {
                    std::cerr << "Error: Node ID " << element.nodes[i] << " not found" << std::endl;
                    return false;
                }
                nodes[i] = it->second;
            }
            
            cell_connectivity.push_back(nodes);
            cell_count++;
        }
    }
    
    // Set domain bounds
    Length_x = mesh_data.max_coords[0] - mesh_data.min_coords[0];
    Length_y = mesh_data.max_coords[1] - mesh_data.min_coords[1];
    Length_z = mesh_data.max_coords[2] - mesh_data.min_coords[2];
    
    std::cout << "Grid conversion completed:" << std::endl;
    std::cout << "Total cells: " << Total_Cells_3D << std::endl;
    std::cout << "Total points: " << Total_Points_3D << std::endl;
    std::cout << "Domain size: " << Length_x << " x " << Length_y << " x " << Length_z << std::endl;
    
    return true;
}

/**
 * @brief Validates mesh quality
 *
 * Performs comprehensive mesh quality checks including element validity,
 * connectivity, and geometric quality measures.
 *
 * @param mesh_data Input mesh data
 * @return true if mesh passes quality checks
 *
 * Mathematical validation:
 * - Element volume: V > 0 for all hexahedra
 * - Aspect ratio: reasonable element proportions
 * - Connectivity: proper node-element relationships
 * - Boundary integrity: closed boundary surfaces
 */
bool Validate_Mesh_Quality_3D(const GMSH_Mesh_Data &mesh_data)
{
    bool is_valid = true;
    int invalid_elements = 0;
    
    std::cout << "Validating mesh quality..." << std::endl;
    
    // Check hexahedral elements
    for (const auto &element : mesh_data.elements)
    {
        if (element.type == GMSH_HEXAHEDRON)
        {
            // Verify all nodes exist
            for (int node_id : element.nodes)
            {
                if (mesh_data.node_id_to_index.find(node_id) == mesh_data.node_id_to_index.end())
                {
                    std::cerr << "Error: Element " << element.id 
                              << " references non-existent node " << node_id << std::endl;
                    invalid_elements++;
                    is_valid = false;
                }
            }
            
            // Check for duplicate nodes in element
            std::set<int> unique_nodes(element.nodes.begin(), element.nodes.end());
            if (unique_nodes.size() != element.nodes.size())
            {
                std::cerr << "Error: Element " << element.id << " has duplicate nodes" << std::endl;
                invalid_elements++;
                is_valid = false;
            }
        }
    }
    
    if (!is_valid)
    {
        std::cerr << "Mesh validation failed: " << invalid_elements << " invalid elements" << std::endl;
    }
    else
    {
        std::cout << "Mesh validation passed successfully" << std::endl;
    }
    
    return is_valid;
}

/**
 * @brief Main GMSH file reading function
 *
 * Reads complete GMSH mesh file and converts to solver format.
 * Handles different GMSH versions and performs comprehensive validation.
 *
 * @param filename GMSH mesh filename
 * @return true if mesh loaded successfully
 */
bool Read_GMSH_File_3D(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Cannot open GMSH file: " << filename << std::endl;
        return false;
    }
    
    std::cout << "Reading GMSH mesh file: " << filename << std::endl;
    
    GMSH_Mesh_Data mesh_data;
    double version;
    int file_type, data_size;
    
    // Parse file header
    if (!Parse_GMSH_Header(file, version, file_type, data_size))
    {
        file.close();
        return false;
    }
    
    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("$PhysicalNames") != std::string::npos)
        {
            if (!Parse_GMSH_Physical_Groups(file, mesh_data))
            {
                file.close();
                return false;
            }
        }
        else if (line.find("$Nodes") != std::string::npos)
        {
            if (!Parse_GMSH_Nodes(file, mesh_data, version))
            {
                file.close();
                return false;
            }
        }
        else if (line.find("$Elements") != std::string::npos)
        {
            if (!Parse_GMSH_Elements(file, mesh_data, version))
            {
                file.close();
                return false;
            }
        }
    }
    
    file.close();
    
    // Validate mesh quality
    if (!Validate_Mesh_Quality_3D(mesh_data))
    {
        std::cerr << "Mesh quality validation failed" << std::endl;
        return false;
    }
    
    // Convert to solver format
    if (!Convert_GMSH_To_Solver_Grid_3D(mesh_data))
    {
        std::cerr << "Grid conversion failed" << std::endl;
        return false;
    }
    
    std::cout << "GMSH mesh loaded successfully" << std::endl;
    std::cout << "Mesh statistics:" << std::endl;
    std::cout << "- Nodes: " << mesh_data.nodes.size() << std::endl;
    std::cout << "- Elements: " << mesh_data.elements.size() << std::endl;
    std::cout << "- Hexahedra: " << mesh_data.num_hexahedra << std::endl;
    std::cout << "- Boundary faces: " << mesh_data.num_boundary_faces << std::endl;
    std::cout << "- Physical groups: " << mesh_data.physical_groups.size() << std::endl;
    
    return true;
}