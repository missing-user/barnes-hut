cmake_minimum_required(VERSION 3.5)

include_directories(include)
add_subdirectory(src)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS  "src/*.cpp")

add_library(Simulation STATIC ${SOURCE_FILES})

target_include_directories(Simulation PUBLIC include)
target_link_directories(Simulation PRIVATE src)

target_link_libraries(Simulation PUBLIC OctTree)
