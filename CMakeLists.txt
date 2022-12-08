cmake_minimum_required(VERSION 3.16)
project(barnes-hut)



set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES})
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)


set(glm_DIR "${PROJECT_SOURCE_DIR}/src/glm/cmake/glm") # if necessary
find_package(glm REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE glm::glm)



find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/tests/*.cpp)
add_executable(tests ${TEST_SRC_FILES})
target_include_directories(tests PRIVATE ${PROJECT_SOURCE_DIR}/src ${GTEST_INCLUDE_DIRS})
target_link_libraries(tests PRIVATE ${GTEST_LIBRARIES} Threads::Threads)
target_link_libraries(tests PRIVATE glm::glm)

enable_testing()
add_test(NAME tests COMMAND tests)
