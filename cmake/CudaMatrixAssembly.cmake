# CUDA Matrix Assembly Build Configuration
# Add this to your main CMakeLists.txt or include as a separate CMake file

# Enable CUDA language support
enable_language(CUDA)

# Find CUDA package
find_package(CUDA REQUIRED)

# Check for required CUDA version (minimum 10.0 for modern features)
if(CUDA_VERSION VERSION_LESS "10.0")
    message(WARNING "CUDA version ${CUDA_VERSION} is older than recommended (10.0+)")
endif()

# Set CUDA standard
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# CUDA compilation flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fPIC")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas -O3")

# Architecture-specific optimizations
# Auto-detect GPU architecture or set manually
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "60;70;75;80;86;89;90")  # Support common architectures
endif()

# Debug vs Release flags
set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -g -O0")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3 -DNDEBUG")

# Find required libraries
find_package(Thrust REQUIRED)
find_package(CUB REQUIRED)

# CUDA Matrix Assembly source files
set(CUDA_MATRIX_ASSEMBLY_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/CUDA_KERNELS/Matrix_Assembly_Cuda_Kernels.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/CUDA_KERNELS/Matrix_Assembly_Cuda_Host_Wrappers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CUDA_KERNELS/Matrix_Assembly_CUDA_Integration_Example.cpp
)

# Create CUDA Matrix Assembly library
add_library(cuda_matrix_assembly ${CUDA_MATRIX_ASSEMBLY_SOURCES})

# Set target properties
set_target_properties(cuda_matrix_assembly PROPERTIES
    CUDA_STANDARD 14
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
    POSITION_INDEPENDENT_CODE ON
)

# Include directories
target_include_directories(cuda_matrix_assembly PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/CUDA_KERNELS
    ${CUDA_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(cuda_matrix_assembly
    ${CUDA_LIBRARIES}
    ${CUDA_CUDART_LIBRARY}
    ${CUDA_curand_LIBRARY}
    ${CUDA_CUBLAS_LIBRARIES}
    ${CUDA_cusparse_LIBRARY}
)

# Add test executable
add_executable(test_matrix_assembly_cuda
    ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix_assembly_cuda.cpp
)

target_link_libraries(test_matrix_assembly_cuda
    cuda_matrix_assembly
    ${CUDA_LIBRARIES}
)

# Add integration to main CFD solver target
if(TARGET CFD_solver_gpu)
    target_link_libraries(CFD_solver_gpu cuda_matrix_assembly)
    target_compile_definitions(CFD_solver_gpu PRIVATE USE_CUDA_MATRIX_ASSEMBLY)
endif()

# Installation rules
install(TARGETS cuda_matrix_assembly
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/CUDA_KERNELS/Matrix_Assembly_Cuda_Kernels.h
    DESTINATION include/cuda_kernels
)

# Print configuration summary
message(STATUS "CUDA Matrix Assembly Configuration:")
message(STATUS "  CUDA Version: ${CUDA_VERSION}")
message(STATUS "  CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "  CUDA Compiler: ${CMAKE_CUDA_COMPILER}")
message(STATUS "  CUDA Flags: ${CMAKE_CUDA_FLAGS}")

# Optional: Add memory checking target
if(CUDA_FOUND)
    add_custom_target(cuda_memcheck
        COMMAND cuda-memcheck --tool=memcheck ./test_matrix_assembly_cuda
        DEPENDS test_matrix_assembly_cuda
        COMMENT "Running CUDA memory checker on matrix assembly tests"
    )
endif()

# Optional: Add profiling target
if(CUDA_FOUND)
    add_custom_target(cuda_profile
        COMMAND nvprof --print-gpu-trace ./test_matrix_assembly_cuda
        DEPENDS test_matrix_assembly_cuda
        COMMENT "Profiling CUDA matrix assembly kernels"
    )
endif()