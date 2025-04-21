include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/TopoHexMeshLib_sheet_singularity/HMeshLib)
include_directories(${CMAKE_SOURCE_DIR}/TopoHexMeshLib_sheet_singularity/src)
set(CMAKE_GENERATOR "Visual Studio 17 2022")
set(CMAKE_GENERATOR_INSTANCE "F:/VS_2022")

cmake_minimum_required(VERSION 3.15)
project(RL_Mesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置调试配置
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDebugDLL")

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

# 为 Visual Studio 项目配置调试选项
if(MSVC)
    target_compile_options(RL_Mesh_demo PRIVATE /Zi /Od /MDd)
    set_target_properties(RL_Mesh_demo PROPERTIES 
        LINK_FLAGS_DEBUG "/DEBUG /INCREMENTAL"
        VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()