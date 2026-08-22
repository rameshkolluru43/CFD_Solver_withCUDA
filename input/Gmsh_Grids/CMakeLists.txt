cmake_minimum_required(VERSION 3.10)
project(FindQuadNeighbors)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(FindQuadNeighbors IdentifyCellsAndNeighbors.cpp)
target_link_libraries(FindQuadNeighbors ${VTK_LIBRARIES})