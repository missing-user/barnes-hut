cmake_minimum_required(VERSION 3.2)
project(barnes-hut)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(${SIM_HEADERS_DIR} ${PROJECT_SOURCE_DIR}/src/Simulation)
set(${QUAD_HEADERS_DIR} ${PROJECT_SOURCE_DIR}/src/QuadTree)

include_directories(${SIM_HEADERS_DIR})
include_directories(${QUAD_HEADERS_DIR})

add_subdirectory(src)
add_subdirectory(tests)