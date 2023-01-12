cmake_minimum_required(VERSION 3.3)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS "*.cpp")

add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})
add_subdirectory(QuadTree) #QuadTreeIndices
add_subdirectory(Simulation)

find_package(Boost CONFIG REQUIRED COMPONENTS program_options timer )

find_package(glm CONFIG REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC glm::glm)
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC Simulation QuadTree)