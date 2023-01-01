cmake_minimum_required(VERSION 3.2)

include_directories(include)
add_subdirectory(src)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS  "src/*.cpp")

add_library(QuadTree STATIC ${SOURCE_FILES})

target_include_directories(QuadTree PUBLIC include)
target_link_directories(QuadTree PRIVATE src)