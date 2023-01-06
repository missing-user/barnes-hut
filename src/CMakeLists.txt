cmake_minimum_required(VERSION 3.2)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS "*.h" "*.cpp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
add_compile_options(-std=c++2a)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})
add_subdirectory(QuadTree)
add_subdirectory(Simulation)

set(glm_DIR "${PROJECT_SOURCE_DIR}/src/glm/cmake/glm") # if necessary
find_package(glm REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC glm::glm)
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC Simulation QuadTree)