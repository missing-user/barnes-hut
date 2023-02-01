cmake_minimum_required(VERSION 3.3)
project(barnes-hut)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Ofast -fopenmp -g -pg") # -Ofast = -O3 -ffast-math -march=native -mtune=native
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native -fopenmp")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(${SIM_HEADERS_DIR} ${PROJECT_SOURCE_DIR}/src/Simulation)
set(${QUAD_HEADERS_DIR} ${PROJECT_SOURCE_DIR}/src/OctTree)

include_directories(${SIM_HEADERS_DIR})
include_directories(${QUAD_HEADERS_DIR})

add_subdirectory(src)
add_subdirectory(tests)