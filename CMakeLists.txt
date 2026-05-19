cmake_minimum_required(VERSION 3.16)

project(CFD_solver LANGUAGES C CXX)

# Default to an optimized CPU build when using single-config generators (Ninja/Makefiles)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug|Release|RelWithDebInfo)" FORCE)
endif()

find_package(OpenMP)
set(CFD_OPENMP_IS_LIBOMP_BREW OFF)
set(CFD_LIBOMP_ROOT "" CACHE PATH "Optional: path to libomp (e.g. Homebrew: brew --prefix libomp)")
if(NOT OpenMP_CXX_FOUND AND APPLE)
	if(CFD_LIBOMP_ROOT STREQUAL "")
		execute_process(
			COMMAND brew --prefix libomp
			OUTPUT_VARIABLE _cfd_libomp_brew
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_QUIET
			RESULT_VARIABLE _cfd_brew_rc)
		if(_cfd_brew_rc EQUAL 0)
			set(CFD_LIBOMP_ROOT "${_cfd_libomp_brew}")
		endif()
	endif()
	if(CFD_LIBOMP_ROOT AND EXISTS "${CFD_LIBOMP_ROOT}/include/omp.h")
		set(OpenMP_CXX_FOUND TRUE)
		set(CFD_OPENMP_IS_LIBOMP_BREW ON)
		message(STATUS "OpenMP (libomp): ${CFD_LIBOMP_ROOT}")
	endif()
endif()
if(OpenMP_CXX_FOUND AND NOT CFD_OPENMP_IS_LIBOMP_BREW)
	message(STATUS "OpenMP enabled for CPU targets: ${OpenMP_CXX_FLAGS}")
elseif(NOT OpenMP_CXX_FOUND)
	message(STATUS "OpenMP not found — install libomp (e.g. brew install libomp) or use Linux GCC for parallel CPU")
endif()

# CPU-only machines (no nvcc): build CFD_solver without CUDA. GPU target requires CUDA toolkit.
find_program(CFD_NVCC_EXECUTABLE nvcc PATHS ENV CUDA_PATH ENV CUDA_HOME PATH_SUFFIXES bin)
if(CFD_NVCC_EXECUTABLE)
	set(CFD_HAVE_CUDA ON)
	set(CMAKE_CUDA_COMPILER ${CFD_NVCC_EXECUTABLE})
	enable_language(CUDA)
	message(STATUS "CUDA compiler found: ${CFD_NVCC_EXECUTABLE}")
else()
	set(CFD_HAVE_CUDA OFF)
	message(STATUS "CUDA not found — building CPU executable CFD_solver only")
endif()

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
# x86 AVX intrinsics (immintrin.h); omit on Apple Silicon / non-x86 builds
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
	list(FILTER SOURCES EXCLUDE REGEX "Grid_Computations_Optimized\\.cpp$")
	list(FILTER SOURCES EXCLUDE REGEX "Grid_Optimization_Integration_Guide\\.cpp$")
	list(FILTER SOURCES EXCLUDE REGEX "Grid_Performance_Benchmark\\.cpp$")
endif()

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

find_package(VTK QUIET COMPONENTS
    CommonCore CommonDataModel CommonExecutionModel
    FiltersCore FiltersGeometry IOXML IOLegacy
    IOCore IOImage RenderingCore
)
find_package(Boost REQUIRED COMPONENTS regex)

# --- GPU executable (optional): full solver with CUDA kernel support ---
if(CFD_HAVE_CUDA)
	file(GLOB GPU_CPU_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
	list(FILTER GPU_CPU_SOURCES EXCLUDE REGEX "Main\\.cpp$")
	list(FILTER GPU_CPU_SOURCES EXCLUDE REGEX "Globals_Stub_GPU\\.cpp$")

	add_executable(CFD_solver_gpu
		${CUDA_SOURCES}
		${CUDA_HOST_WRAPPERS}
		${GPU_CPU_SOURCES}
	)

	# sm_XX: override with -DCFD_CUDA_ARCH=86 (Ampere), 89 (Ada), etc.
	if(NOT DEFINED CFD_CUDA_ARCH)
		set(CFD_CUDA_ARCH "75" CACHE STRING "CUDA arch for -arch=sm_XX (e.g. 75, 86, 89)")
	endif()

	set_target_properties(CFD_solver_gpu PROPERTIES
		CUDA_SEPARABLE_COMPILATION ON
	)
	if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
		set_target_properties(CFD_solver_gpu PROPERTIES CUDA_ARCHITECTURES "${CFD_CUDA_ARCH}")
	else()
		message(STATUS "CMake ${CMAKE_VERSION}: CUDA_ARCHITECTURES ignored; using -arch=sm_${CFD_CUDA_ARCH}")
		target_compile_options(CFD_solver_gpu PRIVATE
			$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_${CFD_CUDA_ARCH}>)
	endif()

	target_compile_definitions(CFD_solver_gpu PRIVATE USE_CUDA=1)
	target_include_directories(CFD_solver_gpu PRIVATE
		${PROJECT_SOURCE_DIR}/include
		${PROJECT_SOURCE_DIR}/CUDA_KERNELS
	)
	if(JSONCPP_INCLUDE_DIR)
		target_include_directories(CFD_solver_gpu PRIVATE ${JSONCPP_INCLUDE_DIR})
	endif()

	if(VTK_FOUND)
		message(STATUS "Found VTK ${VTK_VERSION} at ${VTK_DIR}")
		target_compile_definitions(CFD_solver_gpu PRIVATE USE_VTK=1)
	endif()

	target_link_libraries(CFD_solver_gpu PRIVATE
		cudart
		${VTK_LIBRARIES}
		jsoncpp
		Boost::regex
	)
	if(OpenMP_CXX_FOUND)
		if(CFD_OPENMP_IS_LIBOMP_BREW)
			target_compile_options(CFD_solver_gpu PRIVATE -Xpreprocessor -fopenmp "-I${CFD_LIBOMP_ROOT}/include")
			target_link_libraries(CFD_solver_gpu PRIVATE "${CFD_LIBOMP_ROOT}/lib/libomp.dylib")
		else()
			target_link_libraries(CFD_solver_gpu PUBLIC OpenMP::OpenMP_CXX)
		endif()
	endif()
	target_compile_options(CFD_solver_gpu PRIVATE
		$<$<CONFIG:Release>:-O3 -DNDEBUG>
		$<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
	)
endif()

# --- CPU executable ---
if(VTK_FOUND)
    list(APPEND SOURCES "${PROJECT_SOURCE_DIR}/src/Globals_Config.cpp")
    add_executable(CFD_solver ${SOURCES})
    target_compile_definitions(CFD_solver PRIVATE USE_VTK=1)
    target_link_libraries(CFD_solver PRIVATE ${VTK_LIBRARIES} jsoncpp Boost::regex)
    target_compile_options(CFD_solver PRIVATE
    	$<$<CONFIG:Release>:-O3 -DNDEBUG>
    	$<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
    )
    if(OpenMP_CXX_FOUND)
    	if(CFD_OPENMP_IS_LIBOMP_BREW)
    		target_compile_options(CFD_solver PRIVATE -Xpreprocessor -fopenmp "-I${CFD_LIBOMP_ROOT}/include")
    		target_link_libraries(CFD_solver PRIVATE "${CFD_LIBOMP_ROOT}/lib/libomp.dylib")
    	else()
    		target_link_libraries(CFD_solver PUBLIC OpenMP::OpenMP_CXX)
    	endif()
    endif()
    if(JSONCPP_INCLUDE_DIR)
        target_include_directories(CFD_solver PRIVATE ${JSONCPP_INCLUDE_DIR})
    endif()
else()
    message(STATUS "VTK not found -- skipping CPU target CFD_solver (GPU target will still build)")
endif()
