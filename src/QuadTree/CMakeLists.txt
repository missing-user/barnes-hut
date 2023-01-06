cmake_minimum_required(VERSION 3.2)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS  "*.h" "*.cpp")

add_library(QuadTree SHARED STATIC ${SOURCE_FILES})
