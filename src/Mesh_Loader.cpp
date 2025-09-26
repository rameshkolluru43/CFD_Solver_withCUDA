#include "definitions.h"
#include "Globals.h"
#include "Grid.h"
#include <json/json.h>
#include <fstream>
#include <sstream>
#include <algorithm>

// Forward decl from this TU
static std::vector<double> readCsv1D(const std::string &path);
static std::string toLower(std::string s);
static bool hasExt(const std::string &path, const std::string &ext);

void Read_CSV_Mesh(const std::string &xnodesCsv, const std::string &ynodesCsv)
{
    std::vector<double> xnodes = readCsv1D(xnodesCsv);
    std::vector<double> ynodes = readCsv1D(ynodesCsv);
    if (xnodes.size() < 2 || ynodes.size() < 2)
    {
        std::cerr << "CSV mesh requires at least 2 xnodes and 2 ynodes\n";
        exit(1);
    }

    // Ensure unique sorted nodes (robustness)
    auto uniq_sort = [](std::vector<double> &v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end(), [](double a, double b) { return std::abs(a - b) < 1e-12; }), v.end());
    };
    uniq_sort(xnodes);
    uniq_sort(ynodes);

    const int Nx = static_cast<int>(xnodes.size());
    const int Ny = static_cast<int>(ynodes.size());
    Is_2D_Flow = true;

    // Build flat Points vector [x,y,z]*
    V_D Points;
    Points.reserve(3 * Nx * Ny);
    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            Points.push_back(xnodes[i]);
            Points.push_back(ynodes[j]);
            Points.push_back(0.0);
        }
    }

    // Create interior quad cells with anticlockwise node order
    Cells.clear();
    Boundary_Cells.clear();
    int id = 0;
    for (int j = 0; j < Ny - 1; ++j)
    {
        for (int i = 0; i < Nx - 1; ++i)
        {
            int n00 = j * Nx + i;
            int n10 = j * Nx + (i + 1);
            int n11 = (j + 1) * Nx + (i + 1);
            int n01 = (j + 1) * Nx + i;

            // Anticlockwise: (x_i,y_j)->(x_{i+1},y_j)->(x_{i+1},y_{j+1})->(x_i,y_{j+1})
            std::vector<int> nodeIdx = {n00, n10, n11, n01};
            Cell c(4, id, nodeIdx);
            c.Dimension = 2;

            // Populate Cell_Vertices
            for (int k = 0; k < 4; ++k)
            {
                int n = nodeIdx[k];
                c.Cell_Vertices.push_back(Points[3 * n + 0]);
                c.Cell_Vertices.push_back(Points[3 * n + 1]);
                c.Cell_Vertices.push_back(Points[3 * n + 2]);
            }

            // Enforce anticlockwise strictly using existing helper (also reorders indices if needed)
            Sort_Points_AntiClockWise(c.Cell_Vertices, c.nodeIndices);
            Construct_Cell(c);

            Cells.push_back(c);
            ++id;
        }
    }

    No_Physical_Cells = static_cast<int>(Cells.size());
    std::cout << "Constructed " << No_Physical_Cells << " cells from CSV mesh (" << Nx << " x " << Ny << " nodes)\n";

    // Identify neighbors and boundaries using existing utilities
    Identify_Neighbours(Points, Cells, Boundary_Cells);
    double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
    BoundingBox(Points, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
    Classify_Domain_Boundaries(Cells, Xmin, Xmax, Ymin, Ymax);
    Create_Boundary_Cells_Lists(Cells, Inlet_Cells_List, Exit_Cells_List, Wall_Cells_List);

    // Post checks
    Check_Cells();
}

void Load_Mesh(const std::string &configOrMeshPath)
{
    // If JSON config, parse it and delegate
    if (hasExt(configOrMeshPath, ".json"))
    {
        std::ifstream jf(configOrMeshPath);
        if (!jf.is_open())
        {
            std::cerr << "Could not open JSON config: " << configOrMeshPath << "\n";
            exit(1);
        }
        Json::CharReaderBuilder rb;
        Json::Value root;
        std::string errs;
        if (!Json::parseFromStream(rb, jf, &root, &errs))
        {
            std::cerr << "JSON parse error: " << errs << "\n";
            exit(1);
        }

        if (root.isMember("mesh"))
        {
            const auto &mesh = root["mesh"];
            if (mesh.isMember("xnodes") && mesh.isMember("ynodes"))
            {
                Read_CSV_Mesh(mesh["xnodes"].asString(), mesh["ynodes"].asString());
                return;
            }
            if (mesh.isMember("vtk"))
            {
                Read_GmshVTK_Grid(mesh["vtk"].asString());
                return;
            }
            if (mesh.isMember("txt"))
            {
                Read_Grid(mesh["txt"].asString());
                return;
            }
        }

        // Fallbacks if top-level provides direct path
        if (root.isMember("vtk"))
        {
            Read_GmshVTK_Grid(root["vtk"].asString());
            return;
        }
        if (root.isMember("grid"))
        {
            Read_Grid(root["grid"].asString());
            return;
        }

        std::cerr << "JSON did not contain recognized mesh keys. Expected mesh.xnodes/mesh.ynodes or mesh.vtk or mesh.txt\n";
        exit(1);
    }

    // Otherwise, use file extension to dispatch
    if (hasExt(configOrMeshPath, ".vtk"))
    {
        Read_GmshVTK_Grid(configOrMeshPath);
    }
    else if (hasExt(configOrMeshPath, ".txt"))
    {
        Read_Grid(configOrMeshPath);
    }
    else if (hasExt(configOrMeshPath, ".csv"))
    {
        // If a single CSV path is provided, assume naming convention *_x.csv and infer *_y.csv or vice versa
        std::string xpath = configOrMeshPath;
        std::string ypath = configOrMeshPath;
        auto pos = xpath.find_last_of("/\\");
        std::string dir = (pos == std::string::npos) ? std::string("") : xpath.substr(0, pos + 1);
        std::string file = (pos == std::string::npos) ? xpath : xpath.substr(pos + 1);
        if (file.find("x") != std::string::npos)
        {
            ypath = dir + std::string(file).replace(file.find("x"), 1, "y");
        }
        else if (file.find("y") != std::string::npos)
        {
            xpath = dir + std::string(file).replace(file.find("y"), 1, "x");
        }
        else
        {
            std::cerr << "Provide both xnodes and ynodes CSV via JSON config or use '*x.csv' naming.\n";
            exit(1);
        }
        Read_CSV_Mesh(xpath, ypath);
    }
    else
    {
        std::cerr << "Unsupported mesh/config path: " << configOrMeshPath << "\n";
        std::cerr << "Supported: .json (mesh.xnodes/mesh.ynodes or mesh.vtk or mesh.txt), .vtk, .txt, .csv (paired).\n";
        exit(1);
    }
}

static std::vector<double> readCsv1D(const std::string &path)
{
    std::ifstream f(path);
    if (!f.is_open())
    {
        std::cerr << "Failed to open CSV: " << path << "\n";
        exit(1);
    }
    std::vector<double> vals;
    std::string line;
    while (std::getline(f, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            if (tok.empty())
                continue;
            // allow whitespace separated values too
            std::stringstream ts(tok);
            double v;
            if ((ts >> v))
                vals.push_back(v);
        }
    }
    if (vals.empty())
    {
        std::cerr << "No numeric values found in CSV: " << path << "\n";
        exit(1);
    }
    return vals;
}

static std::string toLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

static bool hasExt(const std::string &path, const std::string &ext)
{
    std::string lp = toLower(path);
    return lp.size() >= ext.size() && lp.substr(lp.size() - ext.size()) == ext;
}
