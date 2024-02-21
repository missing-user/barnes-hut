cmake_minimum_required(VERSION 3.13)

add_executable(${CMAKE_PROJECT_NAME} main.cpp )
target_compile_features(${CMAKE_PROJECT_NAME} PRIVATE cxx_std_17)

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE barneshut)

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ../include)

# Boost (for commandline program options)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})