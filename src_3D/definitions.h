#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <utility>
#include <functional>
#include <memory>
#include <cassert>
#include <chrono>
#include <numeric>

using namespace std;

// Mathematical constants
#define PI 3.14159265358979323846
#define GAMMA 1.4
#define R_GAS 287.0
#define PR 0.72
#define MU_REF 1.8e-5
#define T_REF 288.15
#define RHO_REF 1.225
#define P_REF 101325.0
#define C_P 1005.0
#define C_V 718.0

// Numerical constants
#define EPSILON 1e-12
#define LARGE_NUMBER 1e12
#define SMALL_NUMBER 1e-12
#define MAX_ITERATIONS 10000
#define CONVERGENCE_TOLERANCE 1e-8

// 3D Face definitions for hexahedral cells
#define NUM_FACES_3D 6
#define NUM_NODES_HEX 8
#define NUM_EDGES_HEX 12

// Face orientation constants for 3D hexahedral cells
#define Face_0 0 // Left face   (x- direction)
#define Face_1 1 // Right face  (x+ direction)
#define Face_2 2 // Bottom face (y- direction)
#define Face_3 3 // Top face    (y+ direction)
#define Face_4 4 // Back face   (z- direction)
#define Face_5 5 // Front face  (z+ direction)

// 3D Direction constants
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2
#define NUM_DIMENSIONS 3

// Boundary condition types for 3D
#define BC_INTERIOR 0
#define BC_WALL 1
#define BC_INLET 2
#define BC_OUTLET 3
#define BC_SYMMETRY 4
#define BC_FAR_FIELD 5
#define BC_PERIODIC 6

// Cell types for 3D meshes
#define CELL_TYPE_HEXAHEDRON 12
#define CELL_TYPE_TETRAHEDRON 10
#define CELL_TYPE_PRISM 13
#define CELL_TYPE_PYRAMID 14

// Conservative variable indices (3D Euler/Navier-Stokes)
#define RHO_INDEX 0
#define RHU_INDEX 1
#define RHV_INDEX 2
#define RHW_INDEX 3 // Additional momentum component for 3D
#define ENERGY_INDEX 4
#define NUM_CONSERVATIVE_VARS 5

// Primitive variable indices (3D)
#define PRIM_RHO 0
#define PRIM_U 1
#define PRIM_V 2
#define PRIM_W 3 // Additional velocity component for 3D
#define PRIM_P 4
#define NUM_PRIMITIVE_VARS 5

// Flux scheme types
#define FLUX_VAN_LEER 0
#define FLUX_ROE 1
#define FLUX_AUSM 2
#define FLUX_LLF 3
#define FLUX_HLLC 4

// Limiter types
#define LIMITER_NONE 0
#define LIMITER_MINMOD 1
#define LIMITER_VAN_LEER 2
#define LIMITER_SUPERBEE 3
#define LIMITER_VAN_ALBADA 4

// Time integration schemes
#define TIME_EULER 0
#define TIME_RK2 1
#define TIME_RK3 2
#define TIME_RK4 3

// Grid types for 3D
#define GRID_CARTESIAN_3D 0
#define GRID_UNSTRUCTURED_3D 1
#define GRID_CURVILINEAR_3D 2
#define GRID_MULTI_BLOCK_3D 3

// Solver types
#define SOLVER_EULER_3D 0
#define SOLVER_NAVIER_STOKES_3D 1
#define SOLVER_RANS_3D 2
#define SOLVER_LES_3D 3

// Turbulence models (for future extension)
#define TURB_NONE 0
#define TURB_SA 1
#define TURB_K_EPSILON 2
#define TURB_K_OMEGA 3
#define TURB_SST 4

// File format types
#define FORMAT_VTK 0
#define FORMAT_CGNS 1
#define FORMAT_PLOT3D 2
#define FORMAT_TECPLOT 3

// Memory alignment for CUDA optimization
#define MEMORY_ALIGNMENT 32
#define CUDA_BLOCK_SIZE 256
#define CUDA_MAX_THREADS_PER_BLOCK 1024

// Error codes
#define SUCCESS 0
#define ERROR_FILE_NOT_FOUND -1
#define ERROR_INVALID_INPUT -2
#define ERROR_MEMORY_ALLOCATION -3
#define ERROR_CONVERGENCE_FAILURE -4
#define ERROR_CUDA_FAILURE -5

// Utility macros for 3D operations
#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define MAX3(a, b, c) (max(max((a), (b)), (c)))
#define MIN3(a, b, c) (min(min((a), (b)), (c)))
#define CLAMP(x, min_val, max_val) (max(min_val, min(x, max_val)))

// Vector operations for 3D
#define DOT_PRODUCT_3D(a, b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])
#define MAGNITUDE_3D(v) (sqrt((v)[0] * (v)[0] + (v)[1] * (v)[1] + (v)[2] * (v)[2]))
#define NORMALIZE_3D(v)               \
    do                                \
    {                                 \
        double mag = MAGNITUDE_3D(v); \
        if (mag > EPSILON)            \
        {                             \
            (v)[0] /= mag;            \
            (v)[1] /= mag;            \
            (v)[2] /= mag;            \
        }                             \
    } while (0)

// Cross product macro for 3D
#define CROSS_PRODUCT_3D(result, a, b)                   \
    do                                                   \
    {                                                    \
        (result)[0] = (a)[1] * (b)[2] - (a)[2] * (b)[1]; \
        (result)[1] = (a)[2] * (b)[0] - (a)[0] * (b)[2]; \
        (result)[2] = (a)[0] * (b)[1] - (a)[1] * (b)[0]; \
    } while (0)

// Memory management macros
#define SAFE_DELETE(ptr)   \
    do                     \
    {                      \
        if (ptr)           \
        {                  \
            delete ptr;    \
            ptr = nullptr; \
        }                  \
    } while (0)
#define SAFE_DELETE_ARRAY(ptr) \
    do                         \
    {                          \
        if (ptr)               \
        {                      \
            delete[] ptr;      \
            ptr = nullptr;     \
        }                      \
    } while (0)

// Debug macros
#ifdef DEBUG
#define DBG_PRINT(x) cout << #x << " = " << x << endl
#define DBG_PRINT_3D(v) cout << #v << " = [" << (v)[0] << ", " << (v)[1] << ", " << (v)[2] << "]" << endl
#else
#define DBG_PRINT(x)
#define DBG_PRINT_3D(v)
#endif

// Performance timing macros
#define START_TIMER auto start_time = chrono::high_resolution_clock::now()
#define END_TIMER(name)                                                                     \
    do                                                                                      \
    {                                                                                       \
        auto end_time = chrono::high_resolution_clock::now();                               \
        auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time); \
        cout << name << " took: " << duration.count() << " microseconds" << endl;           \
    } while (0)

// Function declarations for 3D utility functions
double calculate_distance_3d(const vector<double> &p1, const vector<double> &p2);
double calculate_volume_hexahedron(const vector<vector<double>> &vertices);
double calculate_area_quadrilateral_3d(const vector<vector<double>> &vertices);
void calculate_face_normal_3d(const vector<vector<double>> &face_vertices, vector<double> &normal);
void calculate_cell_center_3d(const vector<vector<double>> &vertices, vector<double> &center);

// Template functions for vector operations
template <typename T>
inline T dot_product_3d(const vector<T> &a, const vector<T> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
inline T magnitude_3d(const vector<T> &v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

template <typename T>
inline void cross_product_3d(const vector<T> &a, const vector<T> &b, vector<T> &result)
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

template <typename T>
inline void normalize_3d(vector<T> &v)
{
    T mag = magnitude_3d(v);
    if (mag > static_cast<T>(EPSILON))
    {
        v[0] /= mag;
        v[1] /= mag;
        v[2] /= mag;
    }
}

#endif // DEFINITIONS_H