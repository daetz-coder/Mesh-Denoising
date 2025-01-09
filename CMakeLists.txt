# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

# 设置项目名称和版本
project(DenoiseMesh VERSION 1.0)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 查找 OpenMesh 库
find_package(OpenMesh REQUIRED COMPONENTS Core Tools)

if(NOT OpenMesh_FOUND)
    message(FATAL_ERROR "OpenMesh not found. Please install OpenMesh and try again.")
endif()

# 设置头文件目录
include_directories(${OpenMesh_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

# 添加源文件
add_executable(denoise_obj src/main.cpp src/denoise_obj.cpp)

# 链接 OpenMesh 库
target_link_libraries(denoise_obj OpenMeshCore OpenMeshTools)
