cmake_minimum_required(VERSION 3.13)
project(barnes-hut VERSION 1.0
                  DESCRIPTION "OpenMP Barnes-Hut Simulation with Out-of-core construction and Morton ordering of particles"
                  LANGUAGES CXX)

# Find packages and libraries
find_package(OpenMP)
find_package(Boost REQUIRED COMPONENTS program_options)
message(STATUS "Boost included at ${Boost_DIR}")
find_package(glm REQUIRED)
message(STATUS "GLM included at ${glm_DIR}")
find_package(xsimd REQUIRED)
message(STATUS "xsimd included at ${xsimd_DIR}")
find_package(libmorton REQUIRED)
message(STATUS "libmorton included at ${libmorton_DIR}")

set(CMAKE_EXPORT_COMPILE_COMMANDS true)
add_compile_definitions(_GLIBCXX_PARALLEL)

# Set Release by default
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native -Wall")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -qopt-report-phase=vec -qopt-report=3 -fast -Wall")
    add_link_options("-v")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} /openmp:experimental /Ox /fp:fast /arch:AVX2")
else()
  message(STATUS "Unknown CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
endif()

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Barnes-hut Library code
add_subdirectory(src)
# Example application using the library
add_subdirectory(simulation)
add_subdirectory(tests)