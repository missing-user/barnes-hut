cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${CMAKE_PROJECT_NAME} main.cpp )

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE barneshut)

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ../include)

# Boost (for commandline program options)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})