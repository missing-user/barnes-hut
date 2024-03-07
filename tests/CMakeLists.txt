cmake_minimum_required(VERSION 3.13)
project(tests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

add_executable(tests main.cpp test.cpp)

target_include_directories(tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

find_package(GTest REQUIRED)
target_include_directories(tests PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(tests PRIVATE ${GTEST_LIBRARIES} Threads::Threads)

target_link_libraries(tests PUBLIC barneshut)

enable_testing()
add_test(NAME tests COMMAND tests)