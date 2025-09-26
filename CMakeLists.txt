cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Project information
project(CFD_Solver_withCUDA 
    VERSION 0.1.0
    DESCRIPTION "High-performance CFD solver leveraging NVIDIA CUDA technology"
    LANGUAGES CXX CUDA
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set CUDA standard
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific options
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# CUDA compiler options
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG")

# Find required packages
find_package(CUDA REQUIRED)
find_package(OpenMP)

# Optional packages
find_package(Doxygen)
find_package(GTest)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CUDA_INCLUDE_DIRS})

# Create library for core CFD functions
set(CFD_SOURCES
    src/cfd_solver.cpp
    src/cuda_kernels.cu
    src/utils.cpp
)

add_library(cfd_core ${CFD_SOURCES})

# Set CUDA properties for the library
set_target_properties(cfd_core PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

# Link libraries
target_link_libraries(cfd_core ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_CUSPARSE_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(cfd_core OpenMP::OpenMP_CXX)
endif()

# Main executable
add_executable(cfd_solver src/main.cpp)
target_link_libraries(cfd_solver cfd_core)

# Set CUDA properties for the executable
set_target_properties(cfd_solver PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

# Installation
install(TARGETS cfd_solver cfd_core
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/ DESTINATION include)

# Documentation with Doxygen
if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
    
    # Install documentation
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/docs/html
        DESTINATION share/doc/cfd_solver
        OPTIONAL
    )
endif()

# Testing
if(GTest_FOUND)
    enable_testing()
    
    # Test executable
    add_executable(unit_tests
        tests/test_main.cpp
        tests/test_cfd_solver.cpp
        tests/test_cuda_kernels.cpp
    )
    
    target_link_libraries(unit_tests cfd_core gtest gtest_main)
    
    # Register tests
    add_test(NAME unit_tests COMMAND unit_tests)
    
    # Custom test target
    add_custom_target(test_verbose
        COMMAND ${CMAKE_CTEST_COMMAND} --verbose
        DEPENDS unit_tests
    )
endif()

# Examples
add_subdirectory(examples)

# Packaging
set(CPACK_PACKAGE_NAME "cfd-solver-cuda")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VENDOR "Ramesh Kolluru")
set(CPACK_PACKAGE_CONTACT "rameshkolluru43@gmail.com")

# Source package
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/\\.git/"
    "/\\.gitignore"
    "/\\.github/"
    "\\.swp$"
    "\\.orig$"
    "/CMakeLists\\.txt\\.user$"
    "/private/"
)

include(CPack)

# Print configuration summary
message(STATUS "")
message(STATUS "CFD Solver with CUDA Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  CUDA compiler: ${CMAKE_CUDA_COMPILER}")
message(STATUS "  CUDA version: ${CUDA_VERSION}")
message(STATUS "  OpenMP support: ${OpenMP_CXX_FOUND}")
message(STATUS "  Doxygen found: ${DOXYGEN_FOUND}")
message(STATUS "  GTest found: ${GTest_FOUND}")
message(STATUS "")