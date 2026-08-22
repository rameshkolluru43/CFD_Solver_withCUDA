cmake_minimum_required(VERSION 3.16)

# Prefer -pthread over libpthread.a (nvlink cannot use static pthread during device link)
set(THREADS_PREFER_PTHREAD_FLAG ON)

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
set(CFD_CUDA_ARCHITECTURES "75" CACHE STRING "CUDA architecture list for CFD_solver_gpu, e.g. 50;61;75 or native")
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
          "$ENV{HOME}/deps/usr/include/jsoncpp"
          /opt/homebrew/include /opt/homebrew/Cellar/jsoncpp/1.9.6/include
)

# Include Paths
include_directories(${PROJECT_SOURCE_DIR}/solvers/include)
if(JSONCPP_INCLUDE_DIR)
    include_directories(${JSONCPP_INCLUDE_DIR})
    message(STATUS "Found jsoncpp headers at: ${JSONCPP_INCLUDE_DIR}")
else()
    message(STATUS "jsoncpp headers not found -- GPU build will use header guards to skip JSON includes")
endif()

# Source Files (All .cpp from src/)
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/solvers/src/*.cpp")

# Exclude GPU-only files from CPU source list
list(FILTER SOURCES EXCLUDE REGEX "Globals_Stub_GPU\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "Main_CUDA\\.cu$")
# x86 AVX intrinsics (immintrin.h); omit on Apple Silicon / non-x86 builds
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
	list(FILTER SOURCES EXCLUDE REGEX "Grid_Computations_Optimized\\.cpp$")
	list(FILTER SOURCES EXCLUDE REGEX "Grid_Optimization_Integration_Guide\\.cpp$")
	list(FILTER SOURCES EXCLUDE REGEX "Grid_Performance_Benchmark\\.cpp$")
endif()

# All CUDA kernel sources (2D path; exclude standalone 3D Euler kernels)
file(GLOB CUDA_KERNEL_SOURCES "${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS/*.cu")
list(FILTER CUDA_KERNEL_SOURCES EXCLUDE REGEX "Euler3D_Cuda\\.cu$")

set(CUDA_SOURCES
    ${PROJECT_SOURCE_DIR}/solvers/src/Main_CUDA.cu
    ${CUDA_KERNEL_SOURCES}
)

# Host wrapper .cpp files contain CUDA <<<>>> syntax, must be compiled by nvcc
set(CUDA_HOST_WRAPPERS
    ${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS/Grid_Cuda_Host_Wrappers.cpp
    ${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS/Matrix_Assembly_Cuda_Host_Wrappers.cpp
)
set_source_files_properties(${CUDA_HOST_WRAPPERS} PROPERTIES LANGUAGE CUDA)
set(CUDA_HOST_CPP_WRAPPERS
    ${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS/Inviscid_Flux_Cuda_MainLoop_Wrapper.cpp
)

find_package(VTK QUIET COMPONENTS
    CommonCore CommonDataModel CommonExecutionModel
    FiltersCore FiltersGeometry IOXML IOLegacy
    IOCore IOImage RenderingCore
)
set(BOOST_ROOT "$ENV{HOME}/deps/usr" CACHE PATH "Boost prefix")
find_package(Boost REQUIRED COMPONENTS regex)
# MPI: FindMPI fails on this machine (mpicxx points at missing -dev include).
# Use headers from a user prefix + the runtime SONAME if present.
set(MPI_CXX_FOUND FALSE)
find_path(CFD_MPI_INCLUDE_DIR mpi.h
	PATHS
		"$ENV{HOME}/deps/usr/lib/x86_64-linux-gnu/openmpi/include"
		"${PROJECT_SOURCE_DIR}/../deps/usr/lib/x86_64-linux-gnu/openmpi/include"
		/usr/lib/x86_64-linux-gnu/openmpi/include
		/usr/include
)
if(EXISTS "/usr/lib/x86_64-linux-gnu/libmpi.so.40")
	set(CFD_MPI_LIBRARY "/usr/lib/x86_64-linux-gnu/libmpi.so.40")
else()
	find_library(CFD_MPI_LIBRARY NAMES mpi mpi_cxx
		PATHS /usr/lib/x86_64-linux-gnu /usr/lib
		PATH_SUFFIXES openmpi/lib)
endif()
if(CFD_MPI_INCLUDE_DIR AND CFD_MPI_LIBRARY)
	add_library(CFD_MPI INTERFACE)
	target_include_directories(CFD_MPI INTERFACE "${CFD_MPI_INCLUDE_DIR}")
	target_compile_definitions(CFD_MPI INTERFACE OMPI_SKIP_MPICXX)
	target_link_libraries(CFD_MPI INTERFACE "${CFD_MPI_LIBRARY}")
	set(MPI_CXX_FOUND TRUE)
	message(STATUS "MPI enabled for CFD_solver_mpi: ${CFD_MPI_INCLUDE_DIR} ${CFD_MPI_LIBRARY}")
else()
	message(STATUS "MPI not found -- skipping optional CPU MPI target CFD_solver_mpi")
endif()

# --- GPU executable (optional): full solver with CUDA kernel support ---
if(CFD_HAVE_CUDA)
	file(GLOB GPU_CPU_SOURCES "${PROJECT_SOURCE_DIR}/solvers/src/*.cpp")
	list(FILTER GPU_CPU_SOURCES EXCLUDE REGEX "Main\\.cpp$")
	list(FILTER GPU_CPU_SOURCES EXCLUDE REGEX "Globals_Stub_GPU\\.cpp$")

	add_executable(CFD_solver_gpu
		${CUDA_SOURCES}
		${CUDA_HOST_WRAPPERS}
		${CUDA_HOST_CPP_WRAPPERS}
		${GPU_CPU_SOURCES}
	)

	# sm_XX: override with -DCFD_CUDA_ARCH=86 (Ampere), 89 (Ada), etc.
	if(NOT DEFINED CFD_CUDA_ARCH)
		set(CFD_CUDA_ARCH "75" CACHE STRING "CUDA arch for -arch=sm_XX (e.g. 75, 86, 89)")
	endif()

	# OFF: nvlink cannot link host static libs (libpthread.a); full clean rebuild required after changing this
	set_target_properties(CFD_solver_gpu PROPERTIES
		CUDA_SEPARABLE_COMPILATION OFF
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
		${PROJECT_SOURCE_DIR}/solvers/include
		${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS
		/usr/local/cuda/include
	)
	if(JSONCPP_INCLUDE_DIR)
		target_include_directories(CFD_solver_gpu PRIVATE ${JSONCPP_INCLUDE_DIR})
	endif()

	if(VTK_FOUND)
		message(STATUS "Found VTK ${VTK_VERSION} at ${VTK_DIR}")
		target_compile_definitions(CFD_solver_gpu PRIVATE USE_VTK=1)
	endif()

	# LINK_ONLY: host libs must not be passed to nvlink device link (libpthread.a fails)
	target_link_libraries(CFD_solver_gpu PRIVATE cudart)
	if(VTK_FOUND)
		target_link_libraries(CFD_solver_gpu PRIVATE $<LINK_ONLY:${VTK_LIBRARIES}>)
	endif()
	# Prefer absolute jsoncpp from deps tree (plain -ljsoncpp often misses non-system install).
	find_library(JSONCPP_LIBRARY NAMES jsoncpp
		PATHS /home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu
		NO_DEFAULT_PATH)
	if(NOT JSONCPP_LIBRARY)
		find_library(JSONCPP_LIBRARY NAMES jsoncpp)
	endif()
	if(JSONCPP_LIBRARY)
		target_link_libraries(CFD_solver_gpu PRIVATE $<LINK_ONLY:${JSONCPP_LIBRARY}>)
	else()
		target_link_libraries(CFD_solver_gpu PRIVATE $<LINK_ONLY:jsoncpp>)
	endif()
	target_link_libraries(CFD_solver_gpu PRIVATE $<LINK_ONLY:Boost::regex>)
	if(OpenMP_CXX_FOUND)
		if(CFD_OPENMP_IS_LIBOMP_BREW)
			target_compile_options(CFD_solver_gpu PRIVATE -Xpreprocessor -fopenmp "-I${CFD_LIBOMP_ROOT}/include")
			target_link_libraries(CFD_solver_gpu PRIVATE $<LINK_ONLY:"${CFD_LIBOMP_ROOT}/lib/libomp.dylib">)
		else()
			target_link_libraries(CFD_solver_gpu PRIVATE $<LINK_ONLY:OpenMP::OpenMP_CXX>)
		endif()
	endif()
	target_compile_options(CFD_solver_gpu PRIVATE
		$<$<CONFIG:Release>:-O3 -DNDEBUG>
		$<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
	)
endif()

# --- CPU executable ---
if(VTK_FOUND)
    set(CFD_CPU_SOURCES ${SOURCES})
    if(EXISTS "${PROJECT_SOURCE_DIR}/solvers/src/Globals_Config.cpp")
        list(APPEND CFD_CPU_SOURCES "${PROJECT_SOURCE_DIR}/solvers/src/Globals_Config.cpp")
    else()
        message(STATUS "src/Globals_Config.cpp not found — building CFD_solver without it")
    endif()
    add_executable(CFD_solver ${CFD_CPU_SOURCES})
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
    message(STATUS "VTK not found -- building CPU CFD_solver without VTK (TXT/legacy VTK I/O only)")
    add_executable(CFD_solver ${SOURCES})
    target_link_libraries(CFD_solver PRIVATE Boost::regex)
    if(JSONCPP_LIBRARY)
        target_link_libraries(CFD_solver PRIVATE ${JSONCPP_LIBRARY})
    else()
        target_link_libraries(CFD_solver PRIVATE jsoncpp)
    endif()
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
endif()

if(NOT JSONCPP_LIBRARY)
	find_library(JSONCPP_LIBRARY NAMES jsoncpp
		PATHS /home/ramesh.kolluru/deps/usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu)
endif()

if(MPI_CXX_FOUND)
	add_executable(CFD_solver_mpi ${SOURCES})
	target_compile_definitions(CFD_solver_mpi PRIVATE USE_MPI=1 OMPI_SKIP_MPICXX)
	if(VTK_FOUND)
		target_compile_definitions(CFD_solver_mpi PRIVATE USE_VTK=1)
		target_link_libraries(CFD_solver_mpi PRIVATE ${VTK_LIBRARIES})
	endif()
	target_link_libraries(CFD_solver_mpi PRIVATE Boost::regex CFD_MPI)
	if(JSONCPP_LIBRARY)
		target_link_libraries(CFD_solver_mpi PRIVATE ${JSONCPP_LIBRARY})
	else()
		target_link_libraries(CFD_solver_mpi PRIVATE jsoncpp)
	endif()
	target_compile_options(CFD_solver_mpi PRIVATE
		$<$<CONFIG:Release>:-O3 -DNDEBUG>
		$<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
	)
	if(OpenMP_CXX_FOUND)
		if(CFD_OPENMP_IS_LIBOMP_BREW)
			target_compile_options(CFD_solver_mpi PRIVATE -Xpreprocessor -fopenmp "-I${CFD_LIBOMP_ROOT}/include")
			target_link_libraries(CFD_solver_mpi PRIVATE "${CFD_LIBOMP_ROOT}/lib/libomp.dylib")
		else()
			target_link_libraries(CFD_solver_mpi PUBLIC OpenMP::OpenMP_CXX)
		endif()
	endif()
	if(JSONCPP_INCLUDE_DIR)
		target_include_directories(CFD_solver_mpi PRIVATE ${JSONCPP_INCLUDE_DIR})
	endif()
	message(STATUS "CPU MPI target: CFD_solver_mpi")
endif()

# --- 3D Euler/NS (clean src3d/ path; independent of incomplete src_3D/) ---
set(CFD_3D_HOST_SOURCES
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D.cpp
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D_Flux.cpp
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D_Viscous.cpp
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D_Turbulence.cpp
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D_Ramp.cpp
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D_SWBLI.cpp
	${PROJECT_SOURCE_DIR}/solvers/src3d/Euler3D_Cavity.cpp
)
add_executable(CFD_solver_3d
	${PROJECT_SOURCE_DIR}/solvers/src3d/Main_3D.cpp
	${CFD_3D_HOST_SOURCES}
)
target_include_directories(CFD_solver_3d PRIVATE ${PROJECT_SOURCE_DIR}/solvers/src3d)
target_compile_features(CFD_solver_3d PRIVATE cxx_std_17)
target_compile_options(CFD_solver_3d PRIVATE
	$<$<CONFIG:Release>:-O3 -DNDEBUG>
	$<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
)
if(OpenMP_CXX_FOUND)
	if(CFD_OPENMP_IS_LIBOMP_BREW)
		target_compile_options(CFD_solver_3d PRIVATE -Xpreprocessor -fopenmp "-I${CFD_LIBOMP_ROOT}/include")
		target_link_libraries(CFD_solver_3d PRIVATE "${CFD_LIBOMP_ROOT}/lib/libomp.dylib")
	else()
		target_link_libraries(CFD_solver_3d PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()
message(STATUS "3D CPU target: CFD_solver_3d")

if(CFD_HAVE_CUDA)
	add_executable(CFD_solver_3d_gpu
		${PROJECT_SOURCE_DIR}/solvers/src3d/Main_3D_CUDA.cu
		${CFD_3D_HOST_SOURCES}
		${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS/Euler3D_Cuda.cu
	)
	target_include_directories(CFD_solver_3d_gpu PRIVATE
		${PROJECT_SOURCE_DIR}/solvers/src3d
		${PROJECT_SOURCE_DIR}/solvers/CUDA_KERNELS
	)
	target_compile_features(CFD_solver_3d_gpu PRIVATE cxx_std_17)
	target_link_libraries(CFD_solver_3d_gpu PRIVATE cudart)
	set_target_properties(CFD_solver_3d_gpu PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
	if(NOT DEFINED CFD_CUDA_ARCH)
		set(CFD_CUDA_ARCH "75" CACHE STRING "CUDA arch for -arch=sm_XX (e.g. 75, 86, 89)")
	endif()
	if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
		set_target_properties(CFD_solver_3d_gpu PROPERTIES CUDA_ARCHITECTURES "${CFD_CUDA_ARCH}")
	else()
		target_compile_options(CFD_solver_3d_gpu PRIVATE
			$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_${CFD_CUDA_ARCH}>)
	endif()
	target_compile_options(CFD_solver_3d_gpu PRIVATE
		$<$<CONFIG:Release>:-O3 -DNDEBUG>
		$<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
	)
	if(OpenMP_CXX_FOUND)
		if(CFD_OPENMP_IS_LIBOMP_BREW)
			target_compile_options(CFD_solver_3d_gpu PRIVATE -Xpreprocessor -fopenmp "-I${CFD_LIBOMP_ROOT}/include")
			target_link_libraries(CFD_solver_3d_gpu PRIVATE $<LINK_ONLY:"${CFD_LIBOMP_ROOT}/lib/libomp.dylib">)
		else()
			target_link_libraries(CFD_solver_3d_gpu PRIVATE $<LINK_ONLY:OpenMP::OpenMP_CXX>)
		endif()
	endif()
	message(STATUS "3D CUDA target: CFD_solver_3d_gpu")
endif()
