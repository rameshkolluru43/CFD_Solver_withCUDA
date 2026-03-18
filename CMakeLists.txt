cmake_minimum_required(VERSION 3.16)

project(CFD_solver LANGUAGES C CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect include paths for jsoncpp (Linux / macOS) -- optional for GPU-only builds
find_path(JSONCPP_INCLUDE_DIR json/json.h
    PATHS /usr/include/jsoncpp /usr/local/include/jsoncpp
          /usr/include /usr/local/include
          /opt/homebrew/include /opt/homebrew/Cellar/jsoncpp/1.9.6/include
)

# Include Paths
include_directories(${PROJECT_SOURCE_DIR}/include)
if(JSONCPP_INCLUDE_DIR)
    include_directories(${JSONCPP_INCLUDE_DIR})
    message(STATUS "Found jsoncpp headers at: ${JSONCPP_INCLUDE_DIR}")
else()
    message(STATUS "jsoncpp headers not found -- GPU build will use header guards to skip JSON includes")
endif()

# Source Files (All .cpp from src/)
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Exclude GPU-only files from CPU source list
list(FILTER SOURCES EXCLUDE REGEX "Globals_Stub_GPU\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "Main_CUDA\\.cu$")

# All CUDA kernel sources
file(GLOB CUDA_KERNEL_SOURCES "${PROJECT_SOURCE_DIR}/CUDA_KERNELS/*.cu")

set(CUDA_SOURCES
    ${PROJECT_SOURCE_DIR}/src/Main_CUDA.cu
    ${CUDA_KERNEL_SOURCES}
)

# Host wrapper .cpp files contain CUDA <<<>>> syntax, must be compiled by nvcc
set(CUDA_HOST_WRAPPERS
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Grid_Cuda_Host_Wrappers.cpp
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS/Matrix_Assembly_Cuda_Host_Wrappers.cpp
)
set_source_files_properties(${CUDA_HOST_WRAPPERS} PROPERTIES LANGUAGE CUDA)

# --- GPU executable: full solver with CUDA kernel support ---
# CPU sources for the GPU target (everything except Main.cpp and Globals_Stub_GPU.cpp)
file(GLOB GPU_CPU_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(FILTER GPU_CPU_SOURCES EXCLUDE REGEX "Main\\.cpp$")
list(FILTER GPU_CPU_SOURCES EXCLUDE REGEX "Globals_Stub_GPU\\.cpp$")

add_executable(CFD_solver_gpu
    ${CUDA_SOURCES}
    ${CUDA_HOST_WRAPPERS}
    ${GPU_CPU_SOURCES}
)

set_target_properties(CFD_solver_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "75"
)

target_compile_definitions(CFD_solver_gpu PRIVATE USE_CUDA=1)
target_include_directories(CFD_solver_gpu PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/CUDA_KERNELS
)
if(JSONCPP_INCLUDE_DIR)
    target_include_directories(CFD_solver_gpu PRIVATE ${JSONCPP_INCLUDE_DIR})
endif()

# --- Shared dependencies ---
find_package(VTK QUIET COMPONENTS
    CommonCore CommonDataModel CommonExecutionModel
    FiltersCore FiltersGeometry IOXML IOLegacy
    IOCore IOImage RenderingCore
)
find_package(Boost REQUIRED COMPONENTS regex)

if(VTK_FOUND)
    message(STATUS "Found VTK ${VTK_VERSION} at ${VTK_DIR}")
    target_compile_definitions(CFD_solver_gpu PRIVATE USE_VTK=1)
endif()

target_link_libraries(CFD_solver_gpu
    cudart
    ${VTK_LIBRARIES}
    jsoncpp
    Boost::regex
)

# --- CPU executable ---
if(VTK_FOUND)
    add_executable(CFD_solver ${SOURCES})
    target_compile_definitions(CFD_solver PRIVATE USE_VTK=1)
    target_link_libraries(CFD_solver ${VTK_LIBRARIES} jsoncpp Boost::regex)
    if(JSONCPP_INCLUDE_DIR)
        target_include_directories(CFD_solver PRIVATE ${JSONCPP_INCLUDE_DIR})
    endif()
else()
    message(STATUS "VTK not found -- skipping CPU target CFD_solver (GPU target will still build)")
endif()