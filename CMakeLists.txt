include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/TopoHexMeshLib_sheet_singularity/HMeshLib)
include_directories(${CMAKE_SOURCE_DIR}/TopoHexMeshLib_sheet_singularity/src)
set(CMAKE_GENERATOR "Visual Studio 17 2022")
set(CMAKE_GENERATOR_INSTANCE "F:/VS_2022")

cmake_minimum_required(VERSION 3.15)
project(RL_Mesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PYBIND11_FINDPYTHON ON)
set(pybind11_DIR "C:/Users/28551/AppData/Local/Programs/Python/Python311/Lib/site-packages/pybind11/share/cmake/pybind11")

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED CONFIG)

file(GLOB SOURCES
    main.cpp
    TopoHexMeshLib_sheet_singularity/HMeshLib/*.cpp)

add_executable(RL_Mesh_demo ${SOURCES})
target_link_libraries(RL_Mesh_demo PRIVATE 
    pybind11::embed
    Python::Python)