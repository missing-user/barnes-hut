cmake_minimum_required(VERSION 3.5)
project(barnes-hut)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native")
  if (OPENMP_FOUND)
    add_compile_definitions(_GLIBCXX_PARALLEL)
  endif()
endif()

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(${SIM_HEADERS_DIR} ${PROJECT_SOURCE_DIR}/src/Simulation)
set(${QUAD_HEADERS_DIR} ${PROJECT_SOURCE_DIR}/src/OctTree)

include_directories(${SIM_HEADERS_DIR})
include_directories(${QUAD_HEADERS_DIR})

add_subdirectory(src)
add_subdirectory(tests)