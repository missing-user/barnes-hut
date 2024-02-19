cmake_minimum_required(VERSION 3.5)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS "*.cpp")

add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})
add_subdirectory(OctTree) #QuadTreeIndices
add_subdirectory(Simulation)

find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE Boost::boost Boost::program_options)
#target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
# find_package(Boost home\jimmy\Packages\boost_1_81_0)

find_package(glm CONFIG REQUIRED)
message(STATUS "GLM included at ${glm_DIR}")
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${GLM_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${GLM_INCLUDE_DIRS})

find_package(xsimd CONFIG REQUIRED)
message(STATUS "xsimd included at ${xsimd_DIR}")
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE xsimd)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE xsimd)

find_package(libmorton CONFIG REQUIRED)
message(STATUS "libmorton included at ${libmorton_DIR}")
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE libmorton::libmorton)

target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC Simulation OctTree)