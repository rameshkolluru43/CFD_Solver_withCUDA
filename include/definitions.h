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
#include <algorithm>
#include <chrono>
#include <cfloat>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <set>

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
#include <boost/regex.hpp>

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

// Global Constants to be used in The entire code
#define radian M_PI / 180.0
#define gamma 1.4
#define gamma_M_1 (gamma - 1.0)
#define gamma1 gamma / gamma_M_1
#define R_GC 287.0
#define cp gamma1 *R_GC
#define cv R_GC / (gamma - 1.0)
#define gamma_R gamma *R_GC // This is gamm*R used to calculate the speed of sound
#define EPSILON 1e-6

/*---------------Terinary Functions to be used in Evaluating Limiters ----------------------------------------------*/

#define Phi(r) ((r >= 1) ? 1 : ((r < 1) & (r > 0)) ? r \
                                                   : 0)
#define Sign(a) ((a > 0.0) - (a < 0.0))

// Constants defined for Sutherland's law for evaluating viscosity
#define T_S_Mu 110.4 //	Reference Temperature for Viscosity
#define T_S_K 194.4  //	Reference Temperature for Thermal Conductivity
#define C_1 1.458e-6 //      This divided by sqrt of Temperature
#define T_ref 288.15

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
