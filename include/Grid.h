// File: Mesh.h
// Project: 2D Compressible Navier-Stokes Solver
// Created on: 2025-04-13

#ifndef MESH_H
#define MESH_H
#include "definitions.h"
#include "Globals.h"

/*-----------------Functions Required for Geometric Features------------------------------------------------*/
// Function for reading input grid file
void Read_Grid(const string &);
// Function for reading Gmsh Grid File
bool Read_VTK_Grid(const string &);
bool Read_GmshVTK_Grid(const string &);
void Read_GmshMESH_Grid(const string &);
// Function for reading a structured mesh from CSV node coordinate lists (x and y)
void Read_CSV_Mesh(const std::string &xnodesCsv, const std::string &ynodesCsv);
// Convenience loader that auto-detects by extension or JSON config (mesh.xnodes/mesh.ynodes)
bool Load_Mesh(const std::string &configOrMeshPath);
// Identifying Qaud and Traiangle cells from the grid points and
void Identify_Cells(V_D &, vector<Cell> &, bool, int &);
// Identifying the neighbours of the cells
void Identify_Neighbours(V_D &, vector<Cell> &, vector<Cell> &);
void Identify_ParentCell(vector<Cell> &, vector<Cell> &);
// Function for Constructing Cells from the data read from grid file
bool Form_Cells(const string &);
// Function to consturcting Cell
void Construct_Cell(V_D &, V_D &, V_D &, V_D &);
// Function to construct for Triangle Cells
void Construct_Cell(V_D &, V_D &, V_D &);
// Function to construct for using verticies
void Construct_Cell(V_D &);
void Construct_Cell(Cell &);
void Construct_Cell(V_D &, V_D &, V_D &, V_D &, V_D &, V_D &, V_D &, V_D &);
void Construct_Cell(const int &, const int &, const int &);
void Construct_Cell();
void TagRefinableCells(vector<Cell> &, double &);
void Compute_Gradient_Refinement_Indicator();
bool Apply_Adaptive_Refinement();

void Construct_Co_Volumes(int &);
// Function for Finding Cross Product of two vectors used in evaluating Face areas.
void Evaluate_Cross_Product(V_D &, V_D &, double &);
// Function for Dot Product of a Vector
void Evaluate_Dot_Product(V_D &, V_D &, double &);
// Function Evaluating Face Normals
void Evaluate_Unit_Vector(V_D &);
// Function for testing grid cells for /sum Ai = 0 and negative volumes
void Check_Cells();
// Function construct ghost cells
void Construct_Ghost_Cells();
// This Function Writes the Cells Physical Information to a file
void Write_Cell_Info(const string &);
void Compute_Centroid(Cell &);
void Compute_Centroid(V_D &, V_D &);
void Compute_Centroid(V_D &, V_D &, V_D &, V_D &);
void Compute_Centroid(V_D &, V_D &, V_D &, V_D &, V_D &);
// Function evaluates Face Area formed by pvoid Lax_Fedrichs()oints and also evaluates Face Normals,Face Center and R_Mid_dot_A for calculating Cell volume
void Construct_Face(V_D &, V_D &, V_D &, V_D &);
void Construct_Face(Cell &);
void Construct_Face(V_D &, V_D &, double &, double &, double &);
void Set_Indicies_of_Neighbour_Cells(const int &, const int &, const int &, const int &, const int &);
void mapFacesToCells(std::vector<Cell> &, std::map<std::pair<int, int>, std::set<int>> &);
void printFaceToCells(const std::map<std::pair<int, int>, std::pair<int, int>> &);
void Sort_Neighbours(vector<Cell> &);
int compute_face_priority(double, double);
int get_NoPhysical_Cells();
double Calculate_Vertex_Average(const V_D &, const V_D &);
void Conversion_Factor(V_D &);
// This Function calculates the distance between two cell centers
void Calculate_Cell_Center_Distances();
void Dot_Product(V_D &, V_D &, double &);
void Distance_Between_Points(V_D &, V_D &, double &);
void Sort_Points_AntiClockWise(V_D &);
void Sort_Points_AntiClockWise(V_D &, V_I &);
bool Are_Points_Sorted_AntiClockWise(const std::vector<V_D> &, const V_D &);
void Classify_Domain_Boundaries(std::vector<Cell> &, double &, double &, double &, double &);
bool Is_On_Boundary(double &, double &, double &, double &, double &, double &);
void BoundingBox(V_D &, double &, double &, double &, double &, double &, double &);
void Create_Boundary_Cells_Lists(vector<Cell> &Cells, vector<int> &, vector<int> &, vector<int> &);
string Get_Boundary_Type(V_D &p1, V_D &p2, double &, double &, double &, double &, bool);

#endif // MESH_H
       // --------------------------------------------------------------