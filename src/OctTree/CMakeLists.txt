cmake_minimum_required(VERSION 3.3)

include_directories(include)
add_subdirectory(src)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS  "src/*.cpp")

add_library(OctTree STATIC ${SOURCE_FILES})

target_include_directories(OctTree PUBLIC include)
target_link_directories(OctTree PRIVATE src)