cmake_minimum_required(VERSION 3.5)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS "*.cpp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
add_compile_options(-std=c++2a)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})
add_subdirectory(OctTree) #QuadTreeIndices
add_subdirectory(Simulation)

find_package(Boost CONFIG REQUIRED COMPONENTS program_options)
# find_package(Boost home\jimmy\Packages\boost_1_81_0)

find_package(glm CONFIG REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC glm::glm)
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
find_package(libmorton CONFIG REQUIRED)
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC libmorton)
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC Simulation OctTree)