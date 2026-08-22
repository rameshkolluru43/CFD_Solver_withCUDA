#ifndef _Header_H
#define _Header_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <stdexcept>
#include <sstream>
#include <tuple>

#ifndef __CUDACC__
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#endif

#ifdef USE_VTK
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkIdList.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkFeatureEdges.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#endif
#ifndef __CUDACC__
#if __has_include(<boost/regex.hpp>)
  #include <boost/regex.hpp>
#elif __has_include(<regex>)
  #include <regex>
#endif
#endif

#if __has_include(<json/json.h>)
#include <json/json.h>
#elif __has_include(<jsoncpp/json.h>)
#include <jsoncpp/json.h>
#include <jsoncpp/reader.h>
#include <jsoncpp/value.h>
#include <jsoncpp/writer.h>
#else
#warning "jsoncpp headers not found; ensure JSON parsing components are available or provide stubs."
#endif

using namespace std;
// Filesystem disabled for CUDA prototype build to avoid nvcc host compiler issues.

#ifndef __CUDACC__
inline void CFD_EnsureInputFileExists(const std::string &path, const std::string &context)
{
  if (path.empty())
    throw std::runtime_error(context + ": empty input file path");

  std::error_code ec;
  if (!std::filesystem::exists(path, ec) || ec)
    throw std::runtime_error(context + ": input file not found: " + path);
  if (!std::filesystem::is_regular_file(path, ec) || ec)
    throw std::runtime_error(context + ": input path is not a regular file: " + path);
}

inline void CFD_EnsureParentDirectoryForOutput(const std::string &path, const std::string &context)
{
  if (path.empty())
    throw std::runtime_error(context + ": empty output file path");

  const std::filesystem::path outputPath(path);
  const std::filesystem::path parent = outputPath.parent_path();
  if (parent.empty())
    return;

  std::error_code ec;
  if (!std::filesystem::exists(parent, ec) || ec)
    throw std::runtime_error(context + ": output directory not found: " + parent.string());
  if (!std::filesystem::is_directory(parent, ec) || ec)
    throw std::runtime_error(context + ": output parent is not a directory: " + parent.string());
}

inline std::ifstream CFD_OpenInputFileOrThrow(const std::string &path, const std::string &context)
{
  CFD_EnsureInputFileExists(path, context);
  std::ifstream file(path, std::ios::in);
  if (!file.is_open())
    throw std::runtime_error(context + ": unable to open input file: " + path);
  return file;
}

inline std::ofstream CFD_OpenOutputFileOrThrow(const std::string &path, const std::string &context,
                                              std::ios_base::openmode mode = std::ios::out)
{
  CFD_EnsureParentDirectoryForOutput(path, context);
  std::ofstream file(path, mode);
  if (!file.is_open())
    throw std::runtime_error(context + ": unable to open output file: " + path);
  return file;
}

inline void CFD_RequireFinitePositive(const double value, const std::string &context)
{
  if (!std::isfinite(value) || value <= 0.0)
    throw std::runtime_error(context + ": expected finite positive value, got " + std::to_string(value));
}
#endif

// Global Constants to be used in The entire code
// Physical constants with _CONST suffix to avoid conflicts
const double RADIAN_CONST = M_PI / 180.0;
const double GAMMA_CONST = 1.4;
const double GAMMA_M_1_CONST = (GAMMA_CONST - 1.0);
const double GAMMA1_CONST = GAMMA_CONST / GAMMA_M_1_CONST;
const double R_GC_CONST = 287.0;
const double CP_CONST = GAMMA1_CONST * R_GC_CONST;
const double CV_CONST = R_GC_CONST / (GAMMA_CONST - 1.0);
const double GAMMA_R_CONST = GAMMA_CONST * R_GC_CONST; // gamma*R used to calculate speed of sound
const double EPSILON_CONST = 1e-6;

// Backward compatibility macros (to be gradually phased out)
#define radian RADIAN_CONST
#define gamma GAMMA_CONST
#define gamma_M_1 GAMMA_M_1_CONST
#define gamma1 GAMMA1_CONST
#define R_GC R_GC_CONST
#define cp CP_CONST
#define cv CV_CONST
#define gamma_R GAMMA_R_CONST
#define EPSILON EPSILON_CONST

/*---------------Terinary Functions to be used in Evaluating Limiters ----------------------------------------------*/

#define Phi(r) ((r >= 1) ? 1 : ((r < 1) & (r > 0)) ? r \
                                                   : 0)
#define Sign(a) ((a > 0.0) - (a < 0.0))

// Constants defined for Sutherland's law for evaluating viscosity
const double T_S_MU_CONST = 110.4; // Reference Temperature for Viscosity
const double T_S_K_CONST = 194.4;  // Reference Temperature for Thermal Conductivity
const double C_1_CONST = 1.458e-6; // This divided by sqrt of Temperature
const double T_REF_CONST = 288.15; // Reference Temperature

// Backward compatibility macros
#define T_S_Mu T_S_MU_CONST
#define T_S_K T_S_K_CONST
#define C_1 C_1_CONST
#define T_ref T_REF_CONST

//-------- Face Indicies to be used in the code---------
#define Face_0 0
#define Face_1 1
#define Face_2 2
#define Face_3 3

#ifdef __cplusplus
extern "C"
{
#endif

  string Get_ErrorFileName();
  string Get_Final_Solution_FileName();
  string Get_Initial_Solution_FileName();
  string Get_Grid_VtkFile();
  string Get_SolutionFile();

#ifdef __cplusplus
}
#endif

#endif // #ifndef _Header_H
