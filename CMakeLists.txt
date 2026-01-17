cmake_minimum_required(VERSION 3.16)

project(CFD_solver LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Paths
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    /opt/homebrew/include
    /opt/homebrew/include/vtk-9.4
    /opt/homebrew/Cellar/jsoncpp/1.9.6/include
    /opt/homebrew/include/json
)

# Source Files (All .cpp from src/)
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Separate CUDA sources (avoid double-compiling .cu as part of generic list)
set(CUDA_SOURCES
    ${PROJECT_SOURCE_DIR}/src/Main_CUDA.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Time_Integration_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Flux_Calculations_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Gradient_Calculation_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Viscous_Flux_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Geometry_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Linear_Algebra_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Iterative_Solver_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Advanced_Optimizations_Cuda.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Grid_Cuda_Kernels.cu
    # New Priority 1 & 2 CUDA kernels (Added 2026-01-15)
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Roe_Flux_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/HLLC_LLF_Flux_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/MUSCL_WENO_Reconstruction_Cuda_Kernels.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Advanced_Flux_Schemes_Cuda_Wrappers.cu
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Reconstruction_Schemes_Cuda_Wrappers.cu
)

# Add Grid CUDA Host Wrappers (C++ source that needs CUDA compiler)
set(CUDA_HOST_SOURCES
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Grid_Cuda_Host_Wrappers.cpp
    ${PROJECT_SOURCE_DIR}/src/Globals_Stub_GPU.cpp
)

# Mark CUDA host sources to be compiled with nvcc
set_source_files_properties(
    ${CUDA_HOST_SOURCES}
    PROPERTIES LANGUAGE CUDA
)

# Link Directories
link_directories(
    /opt/homebrew/lib
)

# Executable
add_executable(CFD_solver ${SOURCES})
target_compile_definitions(CFD_solver PRIVATE USE_VTK=1)

# GPU executable (excludes regular main.cpp by replacing with Main_CUDA.cu if present)
add_executable(CFD_solver_gpu
    ${CUDA_SOURCES}
    ${CUDA_HOST_SOURCES}
)

set_target_properties(CFD_solver_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "60;70;75;80;86;90"
)

target_compile_definitions(CFD_solver_gpu PRIVATE USE_CUDA_DRIVER_MAIN=1)
target_include_directories(CFD_solver_gpu PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS
    /opt/homebrew/include
    /opt/homebrew/include/vtk-9.4
    /opt/homebrew/Cellar/jsoncpp/1.9.6/include
    /opt/homebrew/include/json
)

# Link Libraries
target_link_libraries(CFD_solver
    jsoncpp
    vtkCommonCore-9.4
    vtkCommonDataModel-9.4
    vtkCommonExecutionModel-9.4
    vtkFiltersCore-9.4
    vtkFiltersGeometry-9.4
    vtkIOXML-9.4
    vtkRenderingCore-9.4
    vtkIOCore-9.4
    vtkIOImage-9.4
    vtkInteractionStyle-9.4
    vtkRenderingOpenGL2-9.4
    vtksys-9.4
)

target_link_libraries(CFD_solver_gpu
    cuda cudart
    jsoncpp
    vtkCommonCore-9.4
    vtkCommonDataModel-9.4
    vtkCommonExecutionModel-9.4
    vtkFiltersCore-9.4
    vtkFiltersGeometry-9.4
    vtkIOXML-9.4
    vtkRenderingCore-9.4
    vtkIOCore-9.4
    vtkIOImage-9.4
    vtkInteractionStyle-9.4
    vtkRenderingOpenGL2-9.4
    vtksys-9.4
)