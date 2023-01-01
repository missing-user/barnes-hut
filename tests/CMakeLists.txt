cmake_minimum_required(VERSION 3.2)
project(tests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

set(SOURCE_FILES main.cpp test.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE ../src ${GTEST_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GTEST_LIBRARIES} Threads::Threads)
target_link_libraries(tests PUBLIC Simulation QuadTree)

enable_testing()
add_test(NAME tests COMMAND tests)